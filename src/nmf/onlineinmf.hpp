#pragma once

#ifdef _OPENMP
#include <omp.h>
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
    int iter;                                      // iter is the number of iterations performed after the current iteration
    int maxEpochs;                                 // The maximum number of epochs allowed to run
    std::vector<T2> E_mini;                 // contains the minibatches for each dataset, each of size m x minibatchSizes[i]

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
        if (this->minibatchSizes[0] != this->minibatchIdx[0].size()) {
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
        int lastEnd = -1;
        for (arma::uword j = 0; j < numChunks; ++j) {
            // origStart to origEnd marks a chunk of the original matrix
            arma::uword origStart = shuffleOrder[j] * CHUNK_SIZE;
            arma::uword origEnd = (shuffleOrder[j] + 1) * CHUNK_SIZE - 1;
            if (origEnd > dataSize - 1) origEnd = dataSize - 1;
            arma::uvec origIdx = arma::linspace<arma::uvec>(origStart, origEnd, origEnd - origStart + 1);

            // permStart to permEnd marks the location of the chunk in the permuted matrix
            arma::uword permStart = lastEnd + 1;
            arma::uword permEnd = lastEnd + origEnd - origStart + 1;
            this->samplingIdx[i].subvec(permStart, permEnd) = origIdx;
            lastEnd = permEnd;
        }
    }

    void initMinibatch(arma::uword minibatchSize) {
        // Divide the user specified minibatchSize according to dataset sizes
        this->minibatchIdx.clear();
        this->minibatchIdx.resize(this->nDatasets);
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            this->minibatchIdx[i] = arma::zeros<arma::uvec>(1);
        }
        arma::uword idx;
        for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
            idx = this->dataIdxNew[i];
            double ratio = (double)this->ncol_E[idx] / (double)this->nSum;
            this->minibatchSizes[idx] = round(ratio * minibatchSize);
            if (this->minibatchSizes[idx] < 1) {
                throw std::invalid_argument("Please set a larger `minibatchSize`.");
            } else if (this->minibatchSizes[idx] > this->ncol_E[idx]) {
                throw std::invalid_argument("Please set a smaller `minibatchSize`.");
            }
            this->minibatchIdx[idx] = arma::zeros<arma::uvec>(this->minibatchSizes[idx]);
        }
        this->minibatchSizesOrig = this->minibatchSizes;
    }

    void createEmini(T1* Eptr, arma::uword idx) {
        // By default do nothing, only create for H5Mat template
        this->E_mini[idx] = Eptr->cols(this->minibatchIdx[idx]);
    }

    double scaleParam(arma::uword i) {
        if (arma::any(this->dataIdxPrev == i)) return 0;
        else {
            if (this->iter == 1) return 0;
            if (this->iter == 2) return 1 / (double)this->minibatchSizes[i];
            else return (double)(this->iter - 2) / (double)(this->iter - 1);
        }
    }

    void solveHmini() {
        tic();
#ifdef _VERBOSE
        std::cout << "--Solving H of minibatches--  ";
#endif
        arma::mat* Wptr = this->W.get();
        arma::mat given(this->m, this->k);
        arma::uword idx;
        for (int i=0; i<this->dataIdxNew.size(); ++i) {
            idx = this->dataIdxNew[i];
            arma::mat* Vptr = this->Vi[idx].get();
            arma::mat* Hminiptr = this->miniHi[idx].get();
            // T1* Eptr = this->Ei[idx].get();
            T2 Emini = this->E_mini[idx];
            given = *Wptr + *Vptr;
            giventGiven = given.t() * given;
            giventGiven += (*Vptr).t() * (*Vptr) * this->lambda;
            // arma::mat giventInput = given.t() * Eptr->cols(this->minibatchIdx[idx]);
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
            arma::mat* Hminiptr = this->miniHi[idx].get();
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
            // *Bptr += Eptr->cols(this->minibatchIdx[idx]) * *Hminiptr / this->minibatchSizes[idx];
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
                for (arma::uword i = 0; i < this->m; i++) {
                    if ((*Vptr)(i, j) < 0) (*Vptr)(i, j) = 1e-16;
                }
            }
        }
#ifdef _VERBOSE
        std::cout << toc() << " sec" << std::endl;
#endif
    }

    void solveH() {
        tic();
#ifdef _VERBOSE
        std::cout << "--Solving H--  ";
#endif
        arma::mat* Wptr = this->W.get();
        arma::mat given(this->m, this->k);
        // arma::mat B;
        for (int i=0; i<this->nDatasets; ++i) {
            arma::mat* Vptr = this->Vi[i].get();
            arma::mat* Hptr = this->Hi[i].get();
            T1* Eptr = this->Ei[i].get();
            given = *Wptr + *Vptr;
            giventGiven = given.t() * given;
            giventGiven += (*Vptr).t() * (*Vptr) * this->lambda;
            unsigned int dataSize = this->ncol_E[i];
            unsigned int numChunks = dataSize / INMF_CHUNK_SIZE;
            if (numChunks * INMF_CHUNK_SIZE < dataSize) numChunks++;
#pragma omp parallel for schedule(auto)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * INMF_CHUNK_SIZE;
                unsigned int spanEnd = (j + 1) * INMF_CHUNK_SIZE - 1;
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

    void initH() {
#ifdef _VERBOSE
        std::cout << "Initializing empty H matrices" << std::endl;
#endif
        std::unique_ptr<arma::mat> H;
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            H = std::unique_ptr<arma::mat>(new arma::mat);
            *H = arma::zeros<arma::mat>(this->ncol_E[i], this->k);
            this->Hi.push_back(std::move(H));
        }
    }

    void initV() {
#ifdef _VERBOSE
            std::cout << "Initializing V matrices by sampling from input" << std::endl;
#endif
        std::unique_ptr<arma::mat> V;
        std::unique_ptr<arma::mat> VT;
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            V = std::unique_ptr<arma::mat>(new arma::mat);
            VT = std::unique_ptr<arma::mat>(new arma::mat);
            // *V taken from random sampled columns of Ei[i]
            arma::uvec indices = arma::randperm(this->ncol_E[i]).head(this->k);
            *V = this->Ei[i]->cols(indices);
            for (arma::uword j = 0; j < this->k; ++j) {
                V->col(j) /= arma::norm(V->col(j), 2);
            }
            *VT = (*V).t();
            this->Vi.push_back(std::move(V));
            this->ViT.push_back(std::move(VT));
        }
    }

    void initW() {
#ifdef _VERBOSE
        std::cout << "Randomly initializing W matrix" << std::endl;
#endif
        this->W = std::unique_ptr<arma::mat>(new arma::mat);
        this->WT = std::unique_ptr<arma::mat>(new arma::mat);
        *this->W = arma::randu<arma::mat>(this->m, this->k, arma::distr_param(0, 2));
        for (arma::uword i = 0; i < this->k; ++i) {
            this->W->col(i) /= arma::norm(this->W->col(i), 2);
        }
        *this->WT = (*this->W).t();
    }

    void initA(std::vector<arma::mat>& Ainit) {
        // TODO
    }

    void initA() {
        std::unique_ptr<arma::mat> A, Aold;
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            A = std::unique_ptr<arma::mat>(new arma::mat);
            Aold = std::unique_ptr<arma::mat>(new arma::mat);
            *A = arma::zeros<arma::mat>(this->k, this->k);
            *Aold = arma::zeros<arma::mat>(this->k, this->k);
            this->Ai.push_back(std::move(A));
            this->Ai_old.push_back(std::move(Aold));
        }
    }

    void initB(std::vector<arma::mat>& Binit) {
        // TODO
    }

    void initB() {
        std::unique_ptr<arma::mat> B, Bold;
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            B = std::unique_ptr<arma::mat>(new arma::mat);
            Bold = std::unique_ptr<arma::mat>(new arma::mat);
            *B = arma::zeros<arma::mat>(this->m, this->k);
            *Bold = arma::zeros<arma::mat>(this->m, this->k);
            this->Bi.push_back(std::move(B));
            this->Bi_old.push_back(std::move(Bold));
        }
    }

    void runScenario1(arma::uword minibatchSize = 5000, arma::uword maxEpochs = 5, arma::uword maxHALSIter = 1) {
        // Leave dataIdxPrev as empty
        std::cout << "Starting online iNMF scenario 1" << std::endl;
        this->maxEpochs = maxEpochs;
        this->iter = 0;
        this->dataIdxNew = this->dataIdx;
        this->nCellsNew = this->ncol_E;
        this->initMinibatch(minibatchSize);
        this->initW();
        this->initV();
        this->initH();
        this->initA();
        this->initB();
        // Initialize miniHi
        std::unique_ptr<arma::mat> miniH;
        for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
            int idx = this->dataIdxNew[i];
            miniH = std::unique_ptr<arma::mat>(new arma::mat(this->minibatchSizes[idx], this->k));
            this->miniHi.push_back(std::move(miniH));
        }
        // Setup the progress bar
        int totalIters = arma::accu(this->nCellsNew) * maxEpochs / minibatchSize;
        Progress p(totalIters, true);
        // Initial shuffling
        for (arma::uword i = 0; i < this->dataIdxNew.size(); ++i) {
            this->permuteChunkIdx(this->dataIdxNew[i]);
        }
        // Start the main loop
        auto start = std::chrono::high_resolution_clock::now();
        while (this->next()) {
            // The `next()` function will also update the minibatch idx to be
            // used for current iteration, and decide whether to stop the while loop
            // And also prepares the in-memory minibatch data

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
            p.increment();
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Average time per iteration: " << duration.count() / totalIters / 1e6 << " sec" << std::endl;

        this->solveH();
    }

    void scenario2() {

    }

    void scenario3() {

    }

    arma::mat getAi(int i) {
        return *(this->Ai[i].get());
    }

    std::vector<std::unique_ptr<arma::mat>> getAllA() {
        return this->Ai;
    }

    arma::mat getBi(int i) {
        return *(this->Bi[i].get());
    }

    std::vector<std::unique_ptr<arma::mat>> getAllB() {
        return this->Bi;
    }
}; // class ONLINEINMF

template<>
void ONLINEINMF<H5Mat, arma::mat>::permuteChunkIdx(int i) {
    // If T is H5Mat, need to get chunkSize to its `colChunkSize` attribute
    int colChunkSize = this->Ei[i]->colChunkSize;
    arma::uword dataSize = this->ncol_E[i];
    arma::uword numChunks = dataSize / colChunkSize;
    if (numChunks * colChunkSize < dataSize) numChunks++;
    // Get a shuffled vector of numbers from 0 to numChunks-1
    arma::uvec shuffleOrder = arma::randperm<arma::uvec>(numChunks);
    this->samplingIdx[i] = arma::zeros<arma::uvec>(dataSize);
    int lastEnd = -1;
    for (arma::uword j = 0; j < numChunks; ++j) {
        // origStart to origEnd marks a chunk of the original matrix
        arma::uword origStart = shuffleOrder[j] * colChunkSize;
        arma::uword origEnd = (shuffleOrder[j] + 1) * colChunkSize - 1;
        if (origEnd > dataSize - 1) origEnd = dataSize - 1;
        arma::uvec origIdx = arma::linspace<arma::uvec>(origStart, origEnd, origEnd - origStart + 1);

        // permStart to permEnd marks the location of the chunk in the permuted matrix
        arma::uword permStart = lastEnd + 1;
        arma::uword permEnd = lastEnd + origEnd - origStart + 1;
        this->samplingIdx[i].subvec(permStart, permEnd) = origIdx;
        lastEnd = permEnd;
    }
}
} // namespace planc
