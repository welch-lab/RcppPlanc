#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif
#include "bppnnls.hpp"
#include "inmf.hpp"

namespace planc {

template <typename T>
class UINMF : public INMF<T> {
private:
    arma::mat giventGiven;
    std::vector<std::unique_ptr<T>> ulist;
    std::vector<std::unique_ptr<T>> ulistT;
    std::vector<std::unique_ptr<arma::mat>> Ui;
    arma::uvec u;
    arma::vec lambda_i;
    bool uinmf_cleared = false;

    void sampleUandV() {
        // U and V must be sampled from the same random set of cells, thus put together
#ifdef _VERBOSE
            std::cout << "Initializing U and V matrices by sampling from input" << std::endl;
#endif
        std::unique_ptr<arma::mat> U;
        std::unique_ptr<arma::mat> V;
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            U = std::unique_ptr<arma::mat>(new arma::mat(this->u[i], this->k, arma::fill::zeros));
            V = std::unique_ptr<arma::mat>(new arma::mat(this->m, this->k, arma::fill::zeros));
            arma::uvec indices = arma::randperm(this->ncol_E[i]).head(this->k);
            *U = this->ulist[i]->cols(indices);
            *V = this->Ei[i]->cols(indices);
            this->Ui.push_back(std::move(U));
            this->Vi.push_back(std::move(V));
        }
    }

    void initW2() {
        // Initialization is different than regular iNMF, so "2"
        // NOTE that this is also different from the initW2 in onlineINMF
#ifdef _VERBOSE
        std::cout << "Randomly initializing W matrix" << std::endl;
#endif
        this->W = std::unique_ptr<arma::mat>(new arma::mat);
        *this->W = arma::randu<arma::mat>(this->m, this->k, arma::distr_param(0, 2));
    }

    void initH() {
#ifdef _VERBOSE
        std::cout << "Randomly initializing H matrices" << std::endl;
#endif
        std::unique_ptr<arma::mat> H;
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            H = std::unique_ptr<arma::mat>(new arma::mat);
            *H = arma::randu<arma::mat>(this->ncol_E[i], this->k, arma::distr_param(0, 2));
            this->Hi.push_back(std::move(H));
        }
    }

    double computeObjectiveError() {
        double obj = 0;
        arma::mat* Wptr = this->W.get();
        arma::mat L(this->m, this->k);
        for (arma::uword i = 0; i < this->nDatasets; ++i) {
            T* Eptr = this->Ei[i].get();
            T* usptr = this->ulist[i].get();
            arma::mat* Hptr = this->Hi[i].get();
            arma::mat* Vptr = this->Vi[i].get();
            arma::mat* Uptr = this->Ui[i].get();
            double sqnormE = arma::norm<T>(*Eptr, "fro");
            sqnormE *= sqnormE;
            double sqnormUs = arma::norm<T>(*usptr, "fro");
            sqnormUs *= sqnormUs;
            L = *Wptr + *Vptr;
            arma::mat LtL = L.t() * L;
            arma::mat HtH = Hptr->t() * *Hptr;
            arma::mat VtV = Vptr->t() * *Vptr;
            arma::mat UtU = Uptr->t() * *Uptr;
            arma::mat EtL = Eptr->t() * L;
            arma::mat UstU = usptr->t() * *Uptr;
            double TrLtLHtH = arma::trace(LtL * HtH);
            double TrUtUHtH = arma::trace(UtU * HtH);
            double TrHtEtL = arma::trace(Hptr->t() * EtL);
            double TrHtUstU = arma::trace(Hptr->t() * UstU);
            double TrVtVHtH = arma::trace(VtV * HtH);
            obj += sqnormE + sqnormUs + TrLtLHtH +
                (1+this->lambda_i[i]) * TrUtUHtH - 2 * TrHtEtL - 2 * TrHtUstU +
                this->lambda_i[i] * TrVtVHtH;
        }
        return obj;
    }

    void solveH() {
#ifdef _VERBOSE
        std::cout << "--Solving UINMF H--  ";
        std::chrono::system_clock::time_point iter_start_time = std::chrono::system_clock::now();
#endif
        arma::mat* Wptr = this->W.get();
        for (arma::uword i=0; i<this->nDatasets; ++i) {
            arma::mat given(this->m, this->k);
            arma::mat* Vptr = this->Vi[i].get();
            arma::mat* Hptr = this->Hi[i].get();
            arma::mat* Uptr = this->Ui[i].get();
            T* Eptr = this->Ei[i].get();
            T* usptr = this->ulist[i].get();
            arma::mat WV = *Wptr + *Vptr;
            giventGiven = WV.t() * WV;
            giventGiven += (*Vptr).t() * (*Vptr) * this->lambda_i[i];
            if (this->u[i] > 0) giventGiven += (*Uptr).t() * (*Uptr) * (1 + this->lambda_i[i]);

            unsigned int dataSize = this->ncol_E[i];
            unsigned int numChunks = dataSize / this->INMF_CHUNK_SIZE;
            if (numChunks * this->INMF_CHUNK_SIZE < dataSize) numChunks++;
#pragma omp parallel for schedule(dynamic)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * this->INMF_CHUNK_SIZE;
                unsigned int spanEnd = (j + 1) * this->INMF_CHUNK_SIZE - 1;
                if (spanEnd > dataSize - 1) spanEnd = dataSize - 1;
                arma::mat giventInput = WV.t() * (*Eptr).cols(spanStart, spanEnd);
                if (this->u[i] > 0) giventInput += (*Uptr).t() * usptr->cols(spanStart, spanEnd);
                BPPNNLS<arma::mat, arma::vec> subProbH(giventGiven, giventInput, true);
                subProbH.solveNNLS();
                (*Hptr).rows(spanStart, spanEnd) = subProbH.getSolutionMatrix().t();
                giventInput.clear();
            }
        }
        giventGiven.clear();
#ifdef _VERBOSE
        std::chrono::system_clock::time_point iter_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = iter_end_time - iter_start_time;
        std::cout << "Solving H took " << elapsed_seconds.count() << " sec" << std::endl;
#endif
    }


    void solveV() {
#ifdef _VERBOSE
        std::cout << "--Solving UINMF V--  ";
        std::chrono::system_clock::time_point iter_start_time = std::chrono::system_clock::now();
#endif
        arma::mat* Wptr = this->W.get();
        arma::mat giventInput(this->k, this->INMF_CHUNK_SIZE);;
        for (arma::uword i=0; i<this->nDatasets; ++i) {
            arma::mat* Hptr = this->Hi[i].get();\
            giventGiven = (*Hptr).t() * (*Hptr);
            giventGiven *= 1 + this->lambda_i[i];
            arma::mat* Vptr = this->Vi[i].get();
            T* ETptr = this->EiT[i].get();
            unsigned int numChunks = this->m / this->INMF_CHUNK_SIZE;
            if (numChunks * this->INMF_CHUNK_SIZE < this->m) numChunks++;
#pragma omp parallel for schedule(dynamic)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * this->INMF_CHUNK_SIZE;
                unsigned int spanEnd = (j + 1) * this->INMF_CHUNK_SIZE - 1;
                if (spanEnd > this->m - 1) spanEnd = this->m - 1;
                arma::mat giventInput;
                giventInput = (*Hptr).t() * ETptr->cols(spanStart, spanEnd);
                giventInput -= (*Hptr).t() * *Hptr * Wptr->rows(spanStart, spanEnd).t();
                BPPNNLS<arma::mat, arma::vec> subProbV(giventGiven, giventInput, true);
                subProbV.solveNNLS();
                (*Vptr).rows(spanStart, spanEnd) = subProbV.getSolutionMatrix().t();
            }
        }
        giventGiven.clear();
        giventInput.clear();
#ifdef _VERBOSE
        std::chrono::system_clock::time_point iter_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = iter_end_time - iter_start_time;
        std::cout << "Solving V took " << elapsed_seconds.count() << " sec" << std::endl;
#endif
    }

    void solveU() {
#ifdef _VERBOSE
        std::cout << "--Solving UINMF U--  ";
        std::chrono::system_clock::time_point iter_start_time = std::chrono::system_clock::now();
#endif
        arma::mat giventInput(this->k, this->INMF_CHUNK_SIZE);;
        for (arma::uword i=0; i<this->nDatasets; ++i) {
            if (this->u[i] == 0) continue; // skip if no U
            arma::mat* Hptr = this->Hi[i].get();
            arma::mat* Uptr = this->Ui[i].get();
            giventGiven = (*Hptr).t() * (*Hptr);
            giventGiven *= 1 + this->lambda_i[i];
            // T* USptr = this->ulist[i].get();
            T* USTptr = this->ulistT[i].get();
            unsigned int numChunks = this->u[i] / this->INMF_CHUNK_SIZE;
            if (numChunks * this->INMF_CHUNK_SIZE < this->u[i]) numChunks++;
#pragma omp parallel for schedule(dynamic)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * this->INMF_CHUNK_SIZE;
                unsigned int spanEnd = (j + 1) * this->INMF_CHUNK_SIZE - 1;
                if (spanEnd > this->u[i] - 1) spanEnd = this->u[i] - 1;
                arma::mat giventInput;
                giventInput = (*Hptr).t() * USTptr->cols(spanStart, spanEnd);
                BPPNNLS<arma::mat, arma::vec> subProbU(giventGiven, giventInput, true);
                subProbU.solveNNLS();
                (*Uptr).rows(spanStart, spanEnd) = subProbU.getSolutionMatrix().t();
            }
        }
        giventGiven.clear();
        giventInput.clear();
#ifdef _VERBOSE
        std::chrono::system_clock::time_point iter_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = iter_end_time - iter_start_time;
        std::cout << "Solving U took " << elapsed_seconds.count() << " sec" << std::endl;
#endif
    }

    void solveW() {
#ifdef _VERBOSE
        std::cout << "--Solving UINMF W--  ";
        std::chrono::system_clock::time_point iter_start_time = std::chrono::system_clock::now();
#endif
        arma::mat* Wptr = this->W.get();
        giventGiven = arma::zeros<arma::mat>(this->k, this->k); ///
        arma::mat giventInput;
        for (unsigned int j = 0; j < this->nDatasets; ++j) {
            arma::mat* Hptr = this->Hi[j].get();
            giventGiven += (*Hptr).t() * (*Hptr);
        }

        unsigned int numChunks = this->m / this->INMF_CHUNK_SIZE;
        if (numChunks * this->INMF_CHUNK_SIZE < this->m) numChunks++;
        for (unsigned int i = 0; i < numChunks; ++i) {
            unsigned int spanStart = i * this->INMF_CHUNK_SIZE;
            unsigned int spanEnd = (i + 1) * this->INMF_CHUNK_SIZE - 1;
            if (spanEnd > this->m - 1) spanEnd = this->m - 1;
            giventInput = arma::zeros<arma::mat>(this->k, spanEnd - spanStart + 1); ///
            #pragma omp parallel for ordered schedule(dynamic)
            for (unsigned int j = 0; j < this->nDatasets; ++j) {
                T* ETptr = this->EiT[j].get();
                arma::mat* Hptr = this->Hi[j].get();
                // arma::mat* VTptr = this->ViT[j].get();
                arma::mat* Vptr = this->Vi[j].get();
                arma::mat gtIpos = (*Hptr).t() * ETptr->cols(spanStart, spanEnd);
                arma::mat gtIneg = (*Hptr).t() * *Hptr * Vptr->rows(spanStart, spanEnd).t();
                #pragma omp ordered
                {
                giventInput += gtIpos;
                giventInput -= gtIneg;
                }
            }
            BPPNNLS<arma::mat, arma::vec> subProbW(giventGiven, giventInput, true); ///
            subProbW.solveNNLS();
            (*Wptr).rows(spanStart, spanEnd) = subProbW.getSolutionMatrix().t();
            // (*WTptr).cols(spanStart, spanEnd) = subProbW.getSolutionMatrix();
        }
        giventGiven.clear(); ///
        giventInput.clear(); ///
#ifdef _VERBOSE
        std::chrono::system_clock::time_point iter_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = iter_end_time - iter_start_time;
        std::cout << "Solving W took " << elapsed_seconds.count() << " sec" << std::endl;
#endif
    }

public:
    UINMF(std::vector<std::unique_ptr<T>>& Ei,
          std::vector<std::unique_ptr<T>>& ulist,
          arma::uword k, arma::vec lambda) : INMF<T>(Ei, k, 0, true) {
        this->ulist = std::move(ulist);
        this->lambda_i = lambda;
        u = arma::zeros<arma::uvec>(this->nDatasets);
        for (arma::uword i=0; i<this->nDatasets; ++i) {
            u[i] = this->ulist[i]->n_rows;
            T* us = this->ulist[i].get();
            T ut = us->t();
            std::unique_ptr<T> USTptr = std::make_unique<T>(ut);
            this->ulistT.push_back(std::move(USTptr));
        }
    }

    void optimizeUANLS(int niter=30, bool verbose=true) {
        if (verbose) {
            std::cerr << "UINMF started, niter=" << niter << std::endl;
        }
        auto start = std::chrono::high_resolution_clock::now();
        this->sampleUandV();
        this->initW2();
        this->initH();
        Progress p(niter, verbose);
        for (int iter=0; iter<niter; iter++) {
            Rcpp::checkUserInterrupt();
            this->solveH();
            this->solveV();
            this->solveU();
            this->solveW();
            if ( ! p.is_aborted() ) p.increment();
            else break;
        }
        this->objective_err = this->computeObjectiveError();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        if (verbose) {
            std::cerr << "Total time:      " << duration.count() << " sec" << std::endl;
            std::cerr << "Objective error: " << this->objective_err << std::endl;
        }
    }

    arma::mat getUi(int i) {
        return *(this->Ui[i].get());
    }

    std::vector<std::unique_ptr<arma::mat>> getAllU() {
        return this->Ui;
    }

    ~UINMF() {
        if (!this->uinmf_cleared) {
            for (unsigned int i = 0; i < this->Ui.size(); ++i) {
                this->Ui[i].reset();
                this->ulist[i].reset();
            }
            this->cleared = true;
        }
     }
}; // class UINMF

template<>
double UINMF<H5Mat>::computeObjectiveError() {
    double obj = 0;
    arma::mat* Wptr = this->W.get();
    arma::mat L(this->m, this->k);
    for (arma::uword i = 0; i < this->nDatasets; ++i) {
        H5Mat* Eptr = this->Ei[i].get();
        H5Mat* usptr = this->ulist[i].get();
        arma::mat* Hptr = this->Hi[i].get();
        arma::mat* Vptr = this->Vi[i].get();
        arma::mat* Uptr = this->Ui[i].get();
        double sqnormE = Eptr->normF();
        sqnormE *= sqnormE;
        double sqnormUs = usptr->normF();
        sqnormUs *= sqnormUs;
        L = *Wptr + *Vptr;
        arma::mat LtL = L.t() * L;
        arma::mat HtH = Hptr->t() * *Hptr;
        arma::mat VtV = Vptr->t() * *Vptr;
        arma::mat UtU = Uptr->t() * *Uptr;
        arma::mat EtL(this->ncol_E[i], this->k);
        arma::uword numChunks = Eptr->n_cols / Eptr->colChunkSize;
        if (numChunks * Eptr->colChunkSize < Eptr->n_cols) numChunks++;
        for (arma::uword j = 0; j < numChunks; ++j) {
            arma::uword spanStart = j * Eptr->colChunkSize;
            arma::uword spanEnd = (j + 1) * Eptr->colChunkSize - 1;
            if (spanEnd > Eptr->n_cols - 1) spanEnd = Eptr->n_cols - 1;
            EtL.rows(spanStart, spanEnd) = Eptr->cols(spanStart, spanEnd).t() * L;
        }

        arma::mat UstU(this->ncol_E[i], this->k);
        numChunks = usptr->n_cols / usptr->colChunkSize;
        if (numChunks * usptr->colChunkSize < usptr->n_cols) numChunks++;
        for (arma::uword j = 0; j < numChunks; ++j) {
            arma::uword spanStart = j * usptr->colChunkSize;
            arma::uword spanEnd = (j + 1) * usptr->colChunkSize - 1;
            if (spanEnd > usptr->n_cols - 1) spanEnd = usptr->n_cols - 1;
            UstU.rows(spanStart, spanEnd) = usptr->cols(spanStart, spanEnd).t() * *Uptr;
        }

        double TrLtLHtH = arma::trace(LtL * HtH);
        double TrUtUHtH = arma::trace(UtU * HtH);
        double TrHtEtL = arma::trace(Hptr->t() * EtL);
        double TrHtUstU = arma::trace(Hptr->t() * UstU);
        double TrVtVHtH = arma::trace(VtV * HtH);
        obj += sqnormE + sqnormUs + TrLtLHtH +
                (1+this->lambda_i[i]) * TrUtUHtH - 2 * TrHtEtL - 2 * TrHtUstU +
                this->lambda_i[i] * TrVtVHtH;
    }
    return obj;
}

template<>
double UINMF<H5SpMat>::computeObjectiveError() {
    double obj = 0;
    arma::mat* Wptr = this->W.get();
    arma::mat L(this->m, this->k);
    for (arma::uword i = 0; i < this->nDatasets; ++i) {
        H5SpMat* Eptr = this->Ei[i].get();
        H5SpMat* usptr = this->ulist[i].get();
        arma::mat* Hptr = this->Hi[i].get();
        arma::mat* Vptr = this->Vi[i].get();
        arma::mat* Uptr = this->Ui[i].get();
        double sqnormE = Eptr->normF();
        sqnormE *= sqnormE;
        double sqnormUs = usptr->normF();
        sqnormUs *= sqnormUs;
        L = *Wptr + *Vptr;
        arma::mat LtL = L.t() * L;
        arma::mat HtH = Hptr->t() * *Hptr;
        arma::mat VtV = Vptr->t() * *Vptr;
        arma::mat UtU = Uptr->t() * *Uptr;
        arma::mat EtL(this->ncol_E[i], this->k);
        arma::mat UstU(this->ncol_E[i], this->k);
        arma::uword numChunks = Eptr->n_cols / this->INMF_CHUNK_SIZE;
        if (numChunks * this->INMF_CHUNK_SIZE < Eptr->n_cols) numChunks++;
        for (arma::uword j = 0; j < numChunks; ++j) {
            arma::uword spanStart = j * this->INMF_CHUNK_SIZE;
            arma::uword spanEnd = (j + 1) * this->INMF_CHUNK_SIZE - 1;
            if (spanEnd > Eptr->n_cols - 1) spanEnd = Eptr->n_cols - 1;
            EtL.rows(spanStart, spanEnd) = Eptr->cols(spanStart, spanEnd).t() * L;
            UstU.rows(spanStart, spanEnd) = usptr->cols(spanStart, spanEnd).t() * *Uptr;
        }
        double TrLtLHtH = arma::trace(LtL * HtH);
        double TrUtUHtH = arma::trace(UtU * HtH);
        double TrHtEtL = arma::trace(Hptr->t() * EtL);
        double TrHtUstU = arma::trace(Hptr->t() * UstU);
        double TrVtVHtH = arma::trace(VtV * HtH);
        obj += sqnormE + sqnormUs + TrLtLHtH +
                (1+this->lambda_i[i]) * TrUtUHtH - 2 * TrHtEtL - 2 * TrHtUstU +
                this->lambda_i[i] * TrVtVHtH;
    }
    return obj;
}

} // namespace planc
