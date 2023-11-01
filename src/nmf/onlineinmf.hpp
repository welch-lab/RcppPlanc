#pragma once

#ifdef _OPENMP
#include <omp.h>

#include <memory>
#endif
#include "bppnnls.hpp"
#include "inmf.hpp"

#define CHUNK_SIZE 1000

namespace planc {

// T1: Input data type, one of arma::mat, arma::sp_mat, H5Mat or H5SpMat
// T2: in-memory minibatch data type, matching T1,
// has to be arma::mat, arma::sp_mat, arma::mat or arma::sp_mat, respectively
template <typename T1, typename T2>
class ONLINEINMF : public INMF<T1> {
private:

    arma::mat giventGiven;
    std::vector<std::unique_ptr<arma::mat>> miniHi;    // each of size minibatchSizes[i] x k
    std::vector<std::unique_ptr<arma::mat>> Ai, Ai_old;    // each of size k x k
    std::vector<std::unique_ptr<arma::mat>> Bi, Bi_old;    // each of size m x k
    arma::uvec dataIdx;                            // {0...nDataset-1}
    arma::uvec dataIdxPrev;                        // The subset of dataIdx for data with existing factorization
    arma::uvec dataIdxNew;                         // The subset of dataIdx for data to be factorized
    arma::uvec nCellsNew;                          // With same length as dataIdxNew, the subset of ncol_E
    arma::uvec minibatchSizes, minibatchSizesOrig; // Totally `nDataset` elements
    arma::uvec epoch, epochPrev;                   // Number of epochs after/before the current iteration
    bool epochNext;                                // Used when need to remove old information
    std::vector<arma::uvec> samplingIdx;           // Totally `nDataset` vectors, each of size ncol_E[i];
    std::vector<arma::uvec> minibatchIdx;          // Totally `nDataset` vectors, each of size minibatchSizes[i] when regularly within an epoch;
    arma::uword iter;                                      // iter is the number of iterations performed after the current iteration
    arma::uword maxEpochs;                                 // The maximum number of epochs allowed to run
    std::vector<T2> E_mini;                 // contains the minibatches for each dataset, each of size m x minibatchSizes[i]

    // %%%%%%%%%%%%%%% Iteration helper functions %%%%%%%%%%%%%%%%%%%%%%%

    bool next() {
        // Update minibatchIdx and decide whether to stop the while loop
        // `minibatchIdx` collects the sample points to be used for an iteration
        this->iter++;
        arma::uword idx;
        arma::uword start, end;
        if (this->maxEpochs*this->nCellsNew[0] >= this->iter*this->minibatchSizes[this->dataIdxNew[0]]) {
            // This is not the very last iteration within maxEpoch
            for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
                idx = this->dataIdxNew[i];
                // epoch is uvec, automatically casted to the floor of the division
                this->epoch[idx] = this->iter * this->minibatchSizes[idx] / this->ncol_E[idx];
                start = ((this->iter - 1)*this->minibatchSizes[idx]) % this->ncol_E[idx];
                if (this->epoch[idx] == this->epochPrev[idx]) {
                    // This iteration fully stays within the same epoch
                    end = start + this->minibatchSizes[idx] - 1;
                    this->minibatchIdx[idx] = this->samplingIdx[idx].subvec(start, end);
                } else {
                    // This iteration is at the tail of the current epoch
                    // The epoch change indicator only looks at the first dataset
                    if (i == 0) this->epochNext = true;
                    this->epochPrev[idx] = this->epoch[idx];
                    end = this->ncol_E[idx] - 1;
                    arma::uvec oldTail = this->samplingIdx[idx].subvec(start, end);
                    this->permuteChunkIdx(idx);
                    if (oldTail.size() < this->minibatchSizes[idx]) {
                        arma::uvec newHead = this->samplingIdx[idx].subvec(0, this->minibatchSizes[idx] - oldTail.size() - 1);
                        this->minibatchIdx[idx] = arma::join_cols(oldTail, newHead);
                    } else {
                        this->minibatchIdx[idx] = oldTail;
                    }
                }
            }
        } else {
            // This is the very last iteration within maxEpoch
            for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
                idx = this->dataIdxNew[i];
                start = ((this->iter - 1)*this->minibatchSizes[idx]) % this->ncol_E[idx];
                end = this->ncol_E[idx] - 1;
                this->minibatchIdx[idx] = this->samplingIdx[idx].subvec(start, end);
            }
        }
        // Return bool value to indicate whether to stop the while loop
        if (this->minibatchSizes[this->dataIdxNew[0]] != this->minibatchIdx[this->dataIdxNew[0]].size()) {
            // This is the very last iteration within maxEpoch, if it is not as large as minibatchSize
            // then we just ignore it
            return false;
        }

        this->E_mini.clear(); // Need to see whether `clear` or `swap` is better
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            this->E_mini.push_back(T2());
        }

        for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
            idx = this->dataIdxNew[i];
            T1* Eptr = this->Ei[idx].get();
            this->createEmini(Eptr, idx);
        }
        // this->createEmini(); // Templated

        if (this->epoch[this->dataIdxNew[0]] < this->maxEpochs) return true;
        else return false;
    }

    void permuteChunkIdx(int i) {
        /*
        Overall, this function shuffles the indices of the original matrix by chunks
        For example, a dataset with 35 columns, chunkSize = 10:
        {0..9, 10..19, 20..29, 30..35}  becomes  {20..29, 0..9, 30..35, 10..19}
        The chunkSize is mainly considered for the case when we are accessing HDF5 data
        When it matches the chunking dimension of the HDF5 data,
        the access will be the most efficient.
        */
        // If T is H5Mat, need to get chunkSize to its `colChunkSize` attribute, the templated
        // function is defined at the end of this file
        arma::uword dataSize = this->ncol_E[i];
        arma::uword numChunks = dataSize / CHUNK_SIZE;
        if (numChunks * CHUNK_SIZE < dataSize) numChunks++;
        // Get a shuffled vector of numbers from 0 to numChunks-1
        arma::uvec shuffleOrder = arma::randperm<arma::uvec>(numChunks);
        this->samplingIdx[i] = arma::zeros<arma::uvec>(dataSize);
        unsigned int lastEnd = -1u;
        for (arma::uword j = 0; j < numChunks; ++j) {
            // origStart to origEnd marks a chunk of the original matrix
            arma::uword origStart = shuffleOrder[j] * CHUNK_SIZE;
            arma::uword origEnd = (shuffleOrder[j] + 1u) * CHUNK_SIZE - 1u;
            if (origEnd > dataSize - 1) origEnd = dataSize - 1u;
            arma::uvec origIdx = arma::linspace<arma::uvec>(origStart, origEnd, origEnd - origStart + 1);

            // permStart to permEnd marks the location of the chunk in the permuted matrix
            arma::uword permStart = lastEnd + 1u;
            arma::uword permEnd = lastEnd + origEnd - origStart + 1u;
            this->samplingIdx[i].subvec(permStart, permEnd) = origIdx;
            lastEnd = permEnd;
        }
    }

    void initMinibatch(int minibatchSize) {
        // Divide the user specified minibatchSize according to dataset sizes
        this->minibatchIdx.clear();
        this->minibatchIdx.resize(this->nDatasets);
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            this->minibatchIdx[i] = arma::zeros<arma::uvec>(1);
        }
        arma::uword idx;
        this->minibatchSizes.resize(this->nDatasets);
        for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
            idx = this->dataIdxNew[i];
            double ratio = (double)this->ncol_E[idx] / (double)arma::accu(this->nCellsNew);
            this->minibatchSizes[idx] = round(ratio * minibatchSize);
            try {
                if (this->minibatchSizes[idx] < 1) {
                    throw std::invalid_argument("Please set a larger `minibatchSize`.");
                } else if (this->minibatchSizes[idx] > this->ncol_E[idx]) {
                    throw std::invalid_argument("Please set a smaller `minibatchSize`.");
                }
            } catch(std::exception &ex) {
#ifdef USING_R
                std::string ex_str = ex.what();
                Rcpp::stop(ex_str);

#else
                throw ex;
#endif
            }
            this->minibatchIdx[idx] = arma::zeros<arma::uvec>(this->minibatchSizes[idx]);
        }
        this->minibatchSizesOrig = this->minibatchSizes;
    }

    void createEmini(T1* Eptr, arma::uword idx) {
        // By default do nothing, only create for H5Mat template
        this->E_mini[idx] = Eptr->cols(this->minibatchIdx[idx]);
    }

    // %%%%%%%%%%%%%%% Matrix Initiation %%%%%%%%%%%%%%%%%%%%%%%

    void initH() override {
        // Generate place holder of H matrices for all datasets (S1)
#ifdef _VERBOSE
        std::cout << "Initializing empty H matrices" << std::endl;
#endif
        std::unique_ptr<arma::mat> H;
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            H = std::make_unique<arma::mat>();
            *H = arma::zeros<arma::mat>(this->ncol_E[i], this->k);
            this->Hi.push_back(std::move(H));
        }
    }

    void initNewH() {
        // Generate place holder of H matrices only for new datasets (S2&3)
#ifdef _VERBOSE
        std::cout << "Initializing empty H matrices only for new datasets" << std::endl;
#endif
        std::unique_ptr<arma::mat> H;
        unsigned int idx;
        for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
            idx = this->dataIdxNew[i];
            H = std::make_unique<arma::mat>();
            arma::mat* Hptr = H.get();
            *Hptr = arma::zeros<arma::mat>(this->ncol_E[idx], this->k);
            this->Hi.push_back(std::move(H));
        }
    }

    void sampleV() {
#ifdef _VERBOSE
            std::cout << "Initializing V matrices by sampling from input" << std::endl;
#endif
        std::unique_ptr<arma::mat> V;
        std::unique_ptr<arma::mat> VT;
        arma::uword idx;
        for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
            idx = this->dataIdxNew[i];
            V = std::make_unique<arma::mat>();
            // *V taken from random sampled columns of Ei[i]
            arma::mat* Vptr = V.get();
            arma::uvec indices = arma::randperm(this->ncol_E[idx]).head(this->k);
            *Vptr = this->Ei[idx]->cols(indices);
            for (arma::uword j = 0; j < this->k; ++j) {
                Vptr->col(j) /= arma::norm(Vptr->col(j), 2);
            }
            this->Vi.push_back(std::move(V));
        }
    }

    void initW2() {
        // Initialization is different than regular iNMF, so "2"
#ifdef _VERBOSE
        std::cout << "Randomly initializing W matrix" << std::endl;
#endif
        // For scenario 1.
        // When scenario 2, call .initW(givenW, false) in wrapper
        this->W = std::make_unique<arma::mat>();
        arma::mat* Wptr = this->W.get();
        *Wptr = arma::randu<arma::mat>(this->m, this->k, arma::distr_param(0, 2));
        for (arma::uword i = 0; i < this->k; ++i) {
            Wptr->col(i) /= arma::norm(Wptr->col(i), 2);
        }
    }

    void initA() {
        std::unique_ptr<arma::mat> A, Aold;
        for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
            A = std::make_unique<arma::mat>();
            arma::mat* Aptr = A.get();
            Aold = std::make_unique<arma::mat>();
            arma::mat* Aoldptr = Aold.get();
            *Aptr = arma::zeros<arma::mat>(this->k, this->k);
            *Aoldptr = arma::zeros<arma::mat>(this->k, this->k);
            this->Ai.push_back(std::move(A));
            this->Ai_old.push_back(std::move(Aold));
        }
    }

    void initB() {
        std::unique_ptr<arma::mat> B, Bold;
        for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
            B = std::make_unique<arma::mat>();
            arma::mat *Bptr = B.get();
            Bold = std::make_unique<arma::mat>();
            arma::mat *Boldptr = Bold.get();
            *Bptr = arma::zeros<arma::mat>(this->m, this->k);
            *Boldptr = arma::zeros<arma::mat>(this->m, this->k);
            this->Bi.push_back(std::move(B));
            this->Bi_old.push_back(std::move(Bold));
        }
    }

    // %%%%%%%%%%%%%%% Main calculation %%%%%%%%%%%%%%%%%%%%%%%

    void solveHmini() {
        tic();
#ifdef _VERBOSE
        std::cout << "--Solving H of minibatches--  ";
#endif
        arma::mat* Wptr = this->W.get();
        arma::mat given(this->m, this->k);
        arma::uword idx;
        for (arma::uword i=0; i<this->dataIdxNew.size(); ++i) {
            idx = this->dataIdxNew[i];
            arma::mat* Vptr = this->Vi[idx].get();
            arma::mat* Hminiptr = this->miniHi[idx].get(); // `i` but not `idx`, nothing init for prev data
            T2 Emini = this->E_mini[idx];
            given = *Wptr + *Vptr;
            giventGiven = given.t() * given;
            giventGiven += (*Vptr).t() * (*Vptr) * this->lambda;
            arma::mat giventInput = given.t() * Emini;
            BPPNNLS<arma::mat, arma::vec> subProbH(giventGiven, giventInput, true);
            subProbH.solveNNLS();
            *Hminiptr = subProbH.getSolutionMatrix().t();
            giventInput.clear();
        }
        giventGiven.clear();
#ifdef _VERBOSE
        std::cout << toc() << " sec" << std::endl;
#endif
    }

    double scaleParam(arma::uword i) {
        // Helper used in .updateAandB()
        if (arma::any(this->dataIdxPrev == i)) return 0;
        else {
            if (this->iter == 1) return 0;
            if (this->iter == 2) return 1 / (double)this->minibatchSizes[i];
            else return (double)(this->iter - 2) / (double)(this->iter - 1);
        }
    }

    void updateAandB() {
        tic();
#ifdef _VERBOSE
        std::cout << "--Updating A and B--  ";
#endif
        arma::uword idx;
        for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
            idx = this->dataIdxNew[i];
            arma::mat* Aptr = this->Ai[idx].get();
            arma::mat* Aoldptr = this->Ai_old[idx].get();
            arma::mat* Bptr = this->Bi[idx].get();
            arma::mat* Boldptr = this->Bi_old[idx].get();
            arma::mat* Hminiptr = this->miniHi[idx].get(); // `i` but not `idx`, nothing init for prev data
            // T1* Eptr = this->Ei[idx].get();
            T2 Emini = this->E_mini[idx];
            if (this->epoch[this->dataIdxNew[0]] > 0 && this->epochNext) {
                // Remove information older than 2 epochs
                *Aptr -= *Aoldptr;
                *Aoldptr = this->scaleParam(idx) * *Aptr;
                *Bptr -= *Boldptr;
                *Boldptr = this->scaleParam(idx) * *Bptr;
            } else {
                *Aoldptr *= this->scaleParam(idx);
                *Boldptr *= this->scaleParam(idx);
            }

            // HiHit
            *Aptr *= this->scaleParam(idx);
            *Aptr += (*Hminiptr).t() * *Hminiptr / this->minibatchSizes[idx];
            for (arma::uword j = 0; j < this->k; ++j) {
                if ((*Aptr)(j, j) == 0) (*Aptr)(j, j) = 1e-15;
            }
            // XiHit
            *Bptr *= this->scaleParam(idx);
            *Bptr += Emini * *Hminiptr / this->minibatchSizes[idx];
        }
#ifdef _VERBOSE
        std::cout << toc() << " sec" << std::endl;
#endif
    }

    void updateW() {
        tic();
#ifdef _VERBOSE
        std::cout << "--Updating W--  ";
#endif
        arma::mat* Wptr = this->W.get();
        for (arma::uword j = 0; j < this->k; j++) {
            arma::vec numerator = arma::zeros<arma::vec>(this->m);
            double denominator = 0;
            for (arma::uword i = 0; i < this->nDatasets; i++) {
                arma::mat* Aptr = this->Ai[i].get();
                arma::mat* Bptr = this->Bi[i].get();
                arma::mat* Vptr = this->Vi[i].get();
                numerator += Bptr->col(j);
                numerator -= (*Wptr + *Vptr) * Aptr->col(j);
                denominator += (*Aptr)(j, j);
            }
            Wptr->col(j) += numerator / denominator;
            for (arma::uword i = 0; i < this->m; i++) {
                if ((*Wptr)(i, j) < 0) (*Wptr)(i, j) = 1e-16;
            }
        }
#ifdef _VERBOSE
        std::cout << toc() << " sec" << std::endl;
#endif
    }

    void updateV() {
        tic();
#ifdef _VERBOSE
        std::cout << "--Updating V--  ";
#endif
        arma::uword idx;
        arma::mat* Wptr = this->W.get();
        for (arma::uword j = 0; j < this->k; j++) {
            for (arma::uword i = 0; i < this->dataIdxNew.size(); i++) {
                idx = this->dataIdxNew[i];
                arma::mat* Aptr = this->Ai[idx].get();
                arma::mat* Bptr = this->Bi[idx].get();
                arma::mat* Vptr = this->Vi[idx].get();
                Vptr->col(j) += (Bptr->col(j) - (*Wptr + *Vptr*(1 + this->lambda)) * Aptr->col(j)) / ((1 + this->lambda) * (*Aptr)(j, j));
                for (arma::uword k = 0; k < this->m; k++) {
                    if ((*Vptr)(k, j) < 0) (*Vptr)(k, j) = 1e-16;
                }
            }
        }
#ifdef _VERBOSE
        std::cout << toc() << " sec" << std::endl;
#endif
    }

    void solveH() {
        // Solve H for all datasets (S1)
        tic();
#ifdef _VERBOSE
        std::cout << "--Solving H--  ";
#endif
        arma::mat* Wptr = this->W.get();
        arma::mat given(this->m, this->k);
        // arma::mat B;
        for (arma::uword i=0; i<this->nDatasets; ++i) {
            arma::mat* Vptr = this->Vi[i].get();
            arma::mat* Hptr = this->Hi[i].get();
            T1* Eptr = this->Ei[i].get();
            given = *Wptr + *Vptr;
            giventGiven = given.t() * given;
            giventGiven += (*Vptr).t() * (*Vptr) * this->lambda;
            unsigned int dataSize = this->ncol_E[i];
            unsigned int numChunks = dataSize / this->INMF_CHUNK_SIZE;
            if (numChunks * this->INMF_CHUNK_SIZE < dataSize) numChunks++;
#pragma omp parallel for schedule(dynamic) default(none) shared(dataSize, Eptr, numChunks, Hptr, given)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * this->INMF_CHUNK_SIZE;
                unsigned int spanEnd = (j + 1) * this->INMF_CHUNK_SIZE - 1;
                if (spanEnd > dataSize - 1) spanEnd = dataSize - 1;
                arma::mat giventInput = given.t() * (*Eptr).cols(spanStart, spanEnd);
                BPPNNLS<arma::mat, arma::vec> subProbH(giventGiven, giventInput, true);
                subProbH.solveNNLS();
                (*Hptr).rows(spanStart, spanEnd) = subProbH.getSolutionMatrix().t();
                giventInput.clear();
            }
        }
        giventGiven.clear();
#ifdef _VERBOSE
        std::cout << toc() << " sec" << std::endl;
#endif
    }

    // %%%%%%%%%%%%%%% Main loop functions %%%%%%%%%%%%%%%%%%%%%%%

    void solveHALS(arma::uword minibatchSize = 5000, arma::uword inputmaxEpochs = 5,
                   arma::uword maxHALSIter = 1, bool verbose = true) {
        // Main loop of online updating algorithm (S1&2)
        // Universal initialization
        this->maxEpochs = inputmaxEpochs;
        this->iter = 0;
        this->sampleV();
        this->initH();
        this->initA();
        this->initB();
        this->checkK();
        this->initMinibatch(minibatchSize);
        // Initialize miniHi
        std::unique_ptr<arma::mat> miniH;
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            miniH = std::unique_ptr<arma::mat>(new arma::mat(this->minibatchSizes[i], this->k));
            this->miniHi.push_back(std::move(miniH));
        }
        // Setup the progress bar
        int totalIters = arma::accu(this->nCellsNew) * maxEpochs / minibatchSize;
        Progress p(totalIters, verbose);
        // Initial shuffling

        for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
            int idx = dataIdxNew[i];
            this->permuteChunkIdx(idx);
        }
        // Start the main loop
        auto start = std::chrono::high_resolution_clock::now();
        while (this->next()) {
            if (this->iter % 5 == 0) Rcpp::checkUserInterrupt();
            // The `next()` function does:
            // 1. update the minibatch idx to be used for current iteration
            // 2. decide whether to stop the while loop
            // 3. prepare the in-memory minibatch data
            // Step 1: Solve H for the minibatches
            this->solveHmini();
            // Step 2: Update A and B
            this->updateAandB();
            // Step 3: Solve V and W with HALS
            for (arma::uword i = 0; i < maxHALSIter; ++i) {
                this->updateW();
                this->updateV();
            }
            // Reset epoch change indicator
            this->epochNext = false;

            if (!p.is_aborted()) p.increment();
            else break;
        }
        this->solveH();
        this->objective_err = this->computeObjectiveError();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        if (verbose) {
            Rcpp::Rcerr << "Total iterations:  " << totalIters << std::endl;
            Rcpp::Rcerr << "Total time:        " << duration.count() << " sec" << std::endl;
            Rcpp::Rcerr << "Objective error:   " << this->objective_err << std::endl;
        }
    }

    void projectNewData() {
        // Main loop of online updating algorithm (S3)
        tic();
#ifdef _VERBOSE
        std::cout << "--Solving H with only W--  ";
#endif
        this->initNewH();
        arma::mat* Wptr = this->W.get();
        arma::mat given(this->m, this->k);
        unsigned int idx;
        for (unsigned int i = 0; i < this->dataIdxNew.size(); ++i) {
            idx = this->dataIdxNew[i];
            arma::mat* Hptr = this->Hi[i].get(); // `i` but not `idx`, nothing init for prev data
            T1* Eptr = this->Ei[idx].get();
            giventGiven = (*Wptr).t() * (*Wptr);
            unsigned int dataSize = this->ncol_E[idx];
            unsigned int numChunks = dataSize / this->INMF_CHUNK_SIZE;
            if (numChunks * this->INMF_CHUNK_SIZE < dataSize) numChunks++;
#pragma omp parallel for schedule(dynamic) default(none) shared(dataSize, Wptr, Eptr, numChunks, Hptr)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * this->INMF_CHUNK_SIZE;
                unsigned int spanEnd = (j + 1) * this->INMF_CHUNK_SIZE - 1;
                if (spanEnd > dataSize - 1) spanEnd = dataSize - 1;
                arma::mat giventInput = (*Wptr).t() * (*Eptr).cols(spanStart, spanEnd);
                BPPNNLS<arma::mat, arma::vec> subProbH(giventGiven, giventInput, true);
                subProbH.solveNNLS();
                (*Hptr).rows(spanStart, spanEnd) = subProbH.getSolutionMatrix().t();
                giventInput.clear();
            }
            // iNMF is basically (W + V) * HT = E, now we solved W * HT = E, so V = 0
            // std::unique_ptr<arma::mat> V;
            // V = std::unique_ptr<arma::mat>(new arma::mat);
            // *V = arma::zeros<arma::mat>(this->m, this->k);
            // this->Vi.push_back(std::move(V));
        }
        giventGiven.clear();
#ifdef _VERBOSE
        std::cout << toc() << " sec" << std::endl;
#endif
    }

public:
    ONLINEINMF(std::vector<std::unique_ptr<T1>>& Ei, arma::uword k, double lambda) : INMF<T1>(Ei, k, lambda, false) {
        this->dataIdx = arma::linspace<arma::uvec>(0, this->nDatasets - 1, this->nDatasets);
        this->minibatchSizes = arma::zeros<arma::uvec>(this->nDatasets);
        this->epoch = arma::zeros<arma::uvec>(this->nDatasets);
        this->epochPrev = arma::zeros<arma::uvec>(this->nDatasets);
        this->epochNext = false;
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            this->samplingIdx.push_back(arma::zeros<arma::uvec>(this->ncol_E[i]));
        }
    }

    // %%%%%%%%%%%%%%% Public initializers %%%%%%%%%%%%%%%%%%%%%%%
    void initA(std::vector<arma::mat>& Ainit) {
        // Set A matrices for existing datasets (S2)
#ifdef _VERBOSE
        std::cout << "Taking initialized A matrices" << std::endl;
#endif
        std::unique_ptr<arma::mat> A;
        std::unique_ptr<arma::mat> Aold;
        try {
            if (Ainit.size() != this->nDatasets) {
                std::string msg = "Must provide " +
                                  std::string(std::to_string(this->nDatasets)) +
                                  " A matrices";
                throw std::invalid_argument(msg);

            }
        }
        catch(std::exception &ex) {
#ifdef USING_R
        std::string ex_str = ex.what();
        Rcpp::stop(ex_str);

#else
        throw ex;
#endif
    }

        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            try {
                if (Ainit[i].n_rows != this->k || Ainit[i].n_cols != this->k) {
                    std::string msg = "Given As must all be of size " +
                                      std::string(std::to_string(this->k)) + " x " +
                                      std::string(std::to_string(this->k));
                    throw std::invalid_argument(msg);
                }
            } catch(std::exception &ex) {
#ifdef USING_R
                std::string ex_str = ex.what();
                Rcpp::stop(ex_str);

#else
                throw ex;
#endif
            }
            A = std::make_unique<arma::mat>();
            arma::mat* Aptr = A.get();
            Aold = std::make_unique<arma::mat>();
            *Aptr = Ainit[i];
            this->Ai.push_back(std::move(A));
            this->Ai_old.push_back(std::move(Aold));
        }
    }

    void initB(std::vector<arma::mat>& Binit) {
        // Set B matrices for existing datasets (S2)
#ifdef _VERBOSE
    std::cout << "Taking initialized B matrices" << std::endl;
#endif
        try {
            if (Binit.size() != this->nDatasets) {
                std::string msg = "Must provide " +
                                  std::string(std::to_string(this->nDatasets)) +
                                  " B matrices";
                throw std::invalid_argument(msg);
            }
        } catch(std::exception &ex) {
#ifdef USING_R
            std::string ex_str = ex.what();
            Rcpp::stop(ex_str);

#else
            throw ex;
#endif
        }
        std::unique_ptr<arma::mat> B;
        std::unique_ptr<arma::mat> Bold;
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            try {
                if (Binit[i].n_rows != this->m || Binit[i].n_cols != this->k) {
                    std::string msg = "Given Bs must all be of size " +
                                      std::string(std::to_string(this->m)) + " x " +
                                      std::string(std::to_string(this->k));
                    throw std::invalid_argument(msg);
                }
            } catch(std::exception &ex) {
#ifdef USING_R
                std::string ex_str = ex.what();
                Rcpp::stop(ex_str);

#else
                throw ex;
#endif
            }
            B = std::make_unique<arma::mat>();
            arma::mat* Bptr = B.get();
            Bold = std::make_unique<arma::mat>();
            *Bptr = Binit[i];
            this->Bi.push_back(std::move(B));
            this->Bi_old.push_back(std::move(Bold));
        }
    }
    // For setting W and V for existing datasets in S2, use
    // .initW(W, false) and .initV(Vs, false) inherited from INMF class
    // (`false` for no transposition)

    // Scenario 1: Online iNMF on all data as new factorization
    void runOnlineINMF(arma::uword minibatchSize = 5000, arma::uword inputmaxEpochs = 5,
                       arma::uword maxHALSIter = 1, bool verbose = true) {
        if (verbose) {
            Rcpp::Rcerr << "Starting online iNMF scenario 1, factorize all datasets" << std::endl;
        }
        this->dataIdxNew = this->dataIdx;
        this->nCellsNew = this->ncol_E;
        this->initW2();
        this->solveHALS(minibatchSize, inputmaxEpochs, maxHALSIter, verbose);
    }

    // Scenario 2, project == false (default): Online iNMF on new data, factorized upon existing factorization
    // Scenario 3, project == true:  Project new datasets without updating existing factorization
    void runOnlineINMF(std::vector<std::unique_ptr<T1>>& E_new, bool project = false,
                       arma::uword minibatchSize = 5000, arma::uword inputmaxEpochs = 5,
                       arma::uword maxHALSIter = 1, bool verbose = true) {
        // Move new Es into Ei, and manage dataIdxNew, Prev, and nCellsNew
        this->dataIdxPrev = this->dataIdx;
        this->dataIdxNew = arma::linspace<arma::uvec>(this->nDatasets,
                                                      this->nDatasets + E_new.size() - 1,
                                                      E_new.size());
        this->nCellsNew = arma::zeros<arma::uvec>(E_new.size());
        for (arma::uword i = 0; i < E_new.size(); i++) {
            T1* EnewPtr = E_new[i].get();
            this->nCellsNew[i] = EnewPtr->n_cols;
            this->ncol_E.push_back(EnewPtr->n_cols);
            this->samplingIdx.push_back(arma::zeros<arma::uvec>(EnewPtr->n_cols));
            this->Ei.push_back(std::move(E_new[i]));
            this->nDatasets++;
        }
        this->epoch = arma::zeros<arma::uvec>(this->nDatasets);
        this->epochPrev = arma::zeros<arma::uvec>(this->nDatasets);
        this->dataIdx = arma::join_cols(this->dataIdxPrev, this->dataIdxNew);
        assert(this->dataIdx.size() == this->nDatasets);
        assert(arma::all(this->dataIdx == arma::linspace<arma::uvec>(0,
                                                                     this->nDatasets - 1,
                                                                     this->nDatasets)));
        if (!project) {
            if (verbose) {
                Rcpp::Rcerr << "Starting online iNMF scenario 2, " <<
                "update factorization with new datasets" << std::endl;
            }
            this->solveHALS(minibatchSize, inputmaxEpochs, maxHALSIter, verbose);
            this->objective_err = this->computeObjectiveError();
        } else {
            if (verbose) {
                Rcpp::Rcerr << "Starting online iNMF scenario 3, " <<
                "project new datasets without updating existing factorization" << std::endl;
            }
            this->projectNewData();
            // No factorization update, so no objective error
        }

    }

    // %%%%%%%%%%%%%%% Results Getters %%%%%%%%%%%%%%%%%%%%%%%
    arma::mat getAi(arma::uword i) {
        return *(this->Ai[i].get());
    }

    std::vector<std::unique_ptr<arma::mat>> getAllA() {
        return this->Ai;
    }

    arma::mat getBi(arma::uword i) {
        return *(this->Bi[i].get());
    }

    std::vector<std::unique_ptr<arma::mat>> getAllB() {
        return this->Bi;
    }
    // For getting H, W and V, use .getHi(), .getW(), .getVi() inherited from INMF class

}; // class ONLINEINMF

template<>
void ONLINEINMF<H5Mat, arma::mat>::permuteChunkIdx(int i) {
    // If T is H5Mat, need to get chunkSize to its `colChunkSize` attribute
    unsigned int colChunkSize = this->Ei[i]->colChunkSize;
    arma::uword dataSize = this->ncol_E[i];
    arma::uword numChunks = dataSize / colChunkSize;
    if (numChunks * colChunkSize < dataSize) numChunks++;
    // Get a shuffled vector of numbers from 0 to numChunks-1
    arma::uvec shuffleOrder = arma::randperm<arma::uvec>(numChunks);
    this->samplingIdx[i] = arma::zeros<arma::uvec>(dataSize);
    unsigned int lastEnd = -1u;
    for (arma::uword j = 0; j < numChunks; ++j) {
        // origStart to origEnd marks a chunk of the original matrix
        arma::uword origStart = shuffleOrder[j] * colChunkSize;
        arma::uword origEnd = (shuffleOrder[j] + 1u) * colChunkSize - 1u;
        if (origEnd > dataSize - 1u) origEnd = dataSize - 1;
        arma::uvec origIdx = arma::linspace<arma::uvec>(origStart, origEnd, origEnd - origStart + 1u);

        // permStart to permEnd marks the location of the chunk in the permuted matrix
        arma::uword permStart = lastEnd + 1u;
        arma::uword permEnd = lastEnd + origEnd - origStart + 1u;
        this->samplingIdx[i].subvec(permStart, permEnd) = origIdx;
        lastEnd = permEnd;
    }
}
} // namespace planc
