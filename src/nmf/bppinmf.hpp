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
        for (int i=0; i<this->nDatasets; ++i) {
            arma::mat* Vptr = this->Vi[i].get();
            arma::mat* Hptr = this->Hi[i].get();
            T* Eptr = this->Ei[i].get();
            given = *Wptr + *Vptr;
            giventGiven = given.t() * given;
            giventGiven += (*Vptr).t() * (*Vptr) * this->lambda;
            // giventInput = given.t() * (*Eptr);
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

    void solveV() {
        tic();
#ifdef _VERBOSE
        std::cout << "--Solving V--  ";
#endif
        arma::mat* WTptr = this->WT.get();
        arma::mat giventInput(this->k, INMF_CHUNK_SIZE);;
        for (int i=0; i<this->nDatasets; ++i) {
            arma::mat* Hptr = this->Hi[i].get();\
            giventGiven = (*Hptr).t() * (*Hptr);
            giventGiven *= 1 + this->lambda;
            arma::mat* Vptr = this->Vi[i].get();
            arma::mat* VTptr = this->ViT[i].get();
            T* ETptr = this->EiT[i].get();
            unsigned int numChunks = this->m / INMF_CHUNK_SIZE;
            if (numChunks * INMF_CHUNK_SIZE < this->m) numChunks++;
#pragma omp parallel for schedule(auto)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * INMF_CHUNK_SIZE;
                unsigned int spanEnd = (j + 1) * INMF_CHUNK_SIZE - 1;
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

        unsigned int numChunks = this->m / INMF_CHUNK_SIZE;
        if (numChunks * INMF_CHUNK_SIZE < this->m) numChunks++;
        for (unsigned int i = 0; i < numChunks; ++i) {
            unsigned int spanStart = i * INMF_CHUNK_SIZE;
            unsigned int spanEnd = (i + 1) * INMF_CHUNK_SIZE - 1;
            if (spanEnd > this->m - 1) spanEnd = this->m - 1;
            giventInput = arma::zeros<arma::mat>(this->k, spanEnd - spanStart + 1); ///
            #pragma omp parallel for schedule(auto)
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

    void optimizeALS(int niter, bool verbose = true) {
        // execute private functions here
        if (verbose) {
            std::cerr << "BPPINMF optimizeALS started, niter=" << niter << std::endl;
        }
        unsigned int iter = 0;
        double time_total=0;//, obj, delta=100, ;
        Progress p(niter, verbose);
        while (iter < niter ) {
        // while (delta > thresh && iter < niter ) {
            Rcpp::checkUserInterrupt();
            // tic();
#ifdef _VERBOSE
            std::cout << "========Staring iteration "
            << iter+1 << "========" << std::endl;
#endif
            solveH();
            solveV();
            solveW();
            // obj = this->computeObjectiveError();
            // delta = abs(this->objective_err - obj) / ((this->objective_err + obj) / 2);
            iter++;
            // this->objective_err = obj;
            // double time_iter = toc();
#ifdef _VERBOSE
            std::cout << "Objective:  " << obj << std::endl
                    << "Delta:      " << delta << std::endl
                    << "Total time: " << time_iter << " sec" << std::endl;
#endif
            // time_total += time_iter;

            if ( ! p.is_aborted() ) p.increment();
            else break;
        }
        this->objective_err = this->computeObjectiveError();
        if (verbose) {
            // std::cout << "Finished after " << iter << " iterations in " << time_total << " seconds." << std::endl
            std::cerr << "Final objective: " << this->objective_err << std::endl;
                // << "Final delta:     " << delta << std::endl;
        }
    }
};

}
