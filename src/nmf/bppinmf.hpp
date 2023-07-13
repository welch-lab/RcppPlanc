#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif
#include "bppnnls.hpp"
#include "inmf.hpp"

#define ONE_THREAD_MATRIX_SIZE 2000

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
        // arma::mat B;
        for (int i=0; i<this->nDatasets; ++i) {
            arma::mat* Vptr = this->Vi[i].get();
            arma::mat* Hptr = this->Hi[i].get();
            T* Eptr = this->Ei[i].get();
            given = *Wptr + *Vptr;
            giventGiven = given.t() * given;
            giventGiven += (*Vptr).t() * (*Vptr) * this->lambda;
            // giventInput = given.t() * (*Eptr);
            unsigned int dataSize = this->ncol_E[i];
            unsigned int numChunks = dataSize / ONE_THREAD_MATRIX_SIZE;
            if (numChunks * ONE_THREAD_MATRIX_SIZE < dataSize) numChunks++;
#pragma omp parallel for schedule(auto)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * ONE_THREAD_MATRIX_SIZE;
                unsigned int spanEnd = (j + 1) * ONE_THREAD_MATRIX_SIZE - 1;
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
        arma::mat giventInput(this->k, ONE_THREAD_MATRIX_SIZE);;
        for (int i=0; i<this->nDatasets; ++i) {
            arma::mat* Hptr = this->Hi[i].get();\
            giventGiven = (*Hptr).t() * (*Hptr);
            giventGiven *= 1 + this->lambda;
            arma::mat* Vptr = this->Vi[i].get();
            arma::mat* VTptr = this->ViT[i].get();
            T* ETptr = this->EiT[i].get();
            unsigned int numChunks = this->m / ONE_THREAD_MATRIX_SIZE;
            if (numChunks * ONE_THREAD_MATRIX_SIZE < this->m) numChunks++;
#pragma omp parallel for schedule(auto)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * ONE_THREAD_MATRIX_SIZE;
                unsigned int spanEnd = (j + 1) * ONE_THREAD_MATRIX_SIZE - 1;
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

        unsigned int numChunks = this->m / ONE_THREAD_MATRIX_SIZE;
        if (numChunks * ONE_THREAD_MATRIX_SIZE < this->m) numChunks++;
        for (unsigned int i = 0; i < numChunks; ++i) {
            unsigned int spanStart = i * ONE_THREAD_MATRIX_SIZE;
            unsigned int spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
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
    // BPPINMF(std::vector<std::unique_ptr<T>>& Ei, arma::uword k, double lambda,
    //         std::vector<std::unique_ptr<arma::mat>>& Hinit,
    //         std::vector<std::unique_ptr<arma::mat>>& Vinit,
    //         arma::mat& Winit) : INMF<T>(Ei, k, lambda, Hinit, Vinit, Winit) {

    // }

    void optimizeALS(int maxIter, const double thresh, bool verbose) {
        // execute private functions here
        if (verbose) {
            std::cerr << "BPPINMF optimizeALS started, maxIter="
                << maxIter << ", thresh=" << thresh << std::endl;
        }
#ifdef _VERBOSE
        std::cout << "BPPINMF optimizeALS started, maxIter="
            << maxIter << ", thresh=" << thresh << std::endl;
#endif
        unsigned int iter = 0;
        double delta=100, obj;
        Progress p(maxIter, verbose);
        while (delta > thresh && iter < maxIter ) {
            Rcpp::checkUserInterrupt();
            if ( ! p.is_aborted() ) {
#ifdef _VERBOSE
                tic();
                std::cout << "========Staring iteration "
                << iter+1 << "========" << std::endl;
#endif
                solveH();
                solveV();
                solveW();
                obj = this->computeObjectiveError();
                delta = abs(this->objective_err - obj) / ((this->objective_err + obj) / 2);
                iter++;
                this->objective_err = obj;
#ifdef _VERBOSE
                std::cout << "Objective:  " << obj << std::endl
                        << "Delta:      " << delta << std::endl
                        << "Total time: " << toc() << " sec" << std::endl;
#endif
                p.increment();
            } else {
                break;
            }
        }
        if (verbose) {
            std::cerr << "Finished after " << iter << " iterations." << std::endl
                << "Final objective: " << this->objective_err << std::endl
                << "Final delta:     " << delta << std::endl;
        }
    }

    arma::mat getHi(int i) {
        return *(this->Hi[i].get());
    }

    std::vector<std::unique_ptr<arma::mat>> getAllH() {
        return this->Hi;
    }

    arma::mat getVi(int i) {
        return *(this->Vi[i].get());
    }

    std::vector<std::unique_ptr<arma::mat>> getAllV() {
        return this->Vi;
    }

    arma::mat getW() {
        return *(this->W.get());
    }
};

}
