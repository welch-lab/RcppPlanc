#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif
#include "bppnnls.hpp"
#include "inmf.hpp"

namespace planc {

template <typename T>
class BPPINMF : public INMF<T> {
private:
    arma::mat giventGiven;

    void solveH() {
        tic();
#ifdef _VERBOSE
        std::cout << "--Solving H--  ";
#endif
        arma::mat* Wptr = this->W.get();
        arma::mat given(this->m, this->k);
        for (arma::uword i=0; i<this->nDatasets; ++i) {
            arma::mat* Vptr = this->Vi[i].get();
            arma::mat* Hptr = this->Hi[i].get();
            T* Eptr = this->Ei[i].get();
            given = *Wptr + *Vptr;
            giventGiven = given.t() * given;
            giventGiven += (*Vptr).t() * (*Vptr) * this->lambda;
            // giventInput = given.t() * (*Eptr);
            unsigned int dataSize = this->ncol_E[i];
            unsigned int numChunks = dataSize / this->INMF_CHUNK_SIZE;
            if (numChunks * this->INMF_CHUNK_SIZE < dataSize) numChunks++;
#pragma omp parallel for schedule(dynamic)
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

    void solveV() {
        tic();
#ifdef _VERBOSE
        std::cout << "--Solving V--  ";
#endif
        arma::mat* WTptr = this->WT.get();
        arma::mat giventInput(this->k, this->INMF_CHUNK_SIZE);;
        for (arma::uword i=0; i<this->nDatasets; ++i) {
            arma::mat* Hptr = this->Hi[i].get();\
            giventGiven = (*Hptr).t() * (*Hptr);
            giventGiven *= 1 + this->lambda;
            arma::mat* Vptr = this->Vi[i].get();
            arma::mat* VTptr = this->ViT[i].get();
            T* ETptr = this->EiT[i].get();
            unsigned int numChunks = this->m / this->INMF_CHUNK_SIZE;
            if (numChunks * this->INMF_CHUNK_SIZE < this->m) numChunks++;
#pragma omp parallel for schedule(dynamic)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * this->INMF_CHUNK_SIZE;
                unsigned int spanEnd = (j + 1) * this->INMF_CHUNK_SIZE - 1;
                if (spanEnd > this->m - 1) spanEnd = this->m - 1;
                arma::mat giventInput;
                giventInput = (*Hptr).t() * (*ETptr).cols(spanStart, spanEnd);
                giventInput -= (*Hptr).t() * *Hptr * (*WTptr).cols(spanStart, spanEnd);
                BPPNNLS<arma::mat, arma::vec> subProbV(giventGiven, giventInput, true);
                subProbV.solveNNLS();
                (*Vptr).rows(spanStart, spanEnd) = subProbV.getSolutionMatrix().t();
                (*VTptr).cols(spanStart, spanEnd) = subProbV.getSolutionMatrix();
            }
        }
        giventGiven.clear();
        giventInput.clear();
#ifdef _VERBOSE
        std::cout << toc() << " sec" << std::endl;
#endif
    }

    void solveW() {
        tic();
#ifdef _VERBOSE
        std::cout << "--Solving W--  ";
#endif
        arma::mat* Wptr = this->W.get();
        arma::mat* WTptr = this->WT.get();
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
            #pragma omp parallel for schedule(dynamic)
            for (unsigned int j = 0; j < this->nDatasets; ++j) {
                T* ETptr = this->EiT[j].get();
                arma::mat* Hptr = this->Hi[j].get();
                arma::mat* VTptr = this->ViT[j].get();
                giventInput += (*Hptr).t() * (*ETptr).cols(spanStart, spanEnd);
                giventInput -= (*Hptr).t() * *Hptr * (*VTptr).cols(spanStart, spanEnd);
            }
            BPPNNLS<arma::mat, arma::vec> subProbW(giventGiven, giventInput, true); ///
            subProbW.solveNNLS();
            (*Wptr).rows(spanStart, spanEnd) = subProbW.getSolutionMatrix().t();
            (*WTptr).cols(spanStart, spanEnd) = subProbW.getSolutionMatrix();
        }
        giventGiven.clear(); ///
        giventInput.clear(); ///
        *WTptr = (*Wptr).t();
#ifdef _VERBOSE
        std::cout << toc() << " sec" << std::endl;
#endif
    }

public:
    BPPINMF(std::vector<std::unique_ptr<T>>& Ei, arma::uword k, double lambda) : INMF<T>(Ei, k, lambda) {

    }

    void optimizeALS(unsigned int niter, bool verbose = true) {
        // execute private functions here
        if (verbose) {
            std::cerr << "INMF started, niter=" << niter << std::endl;
        }
        this->objective_err = this->computeObjectiveError();
        auto start = std::chrono::high_resolution_clock::now();
        unsigned int iter = 0;
        Progress p(niter, verbose);
        while (iter < niter ) {
            Rcpp::checkUserInterrupt();
            solveH();
            solveV();
            solveW();
            iter++;
            if ( ! p.is_aborted() ) p.increment();
            else break;
        }
        this->objective_err = this->computeObjectiveError();
        auto end = std::chrono::high_resolution_clock::now();
        if (verbose) {
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
            std::cerr << "Total time:      " << duration.count() << " sec" << std::endl;
            std::cerr << "Objective error: " << this->objective_err << std::endl;
        }
    }
};

}
