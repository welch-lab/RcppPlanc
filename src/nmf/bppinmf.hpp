#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif
#include "bppnnls.hpp"
#include "inmf.hpp"

#define ONE_THREAD_MATRIX_SIZE 2000

namespace planc {

template <class T>
class BPPINMF : INMF<T> {
private:
    void solveH() {
        // implement
        for (int i=0; i<this->nDatasets; i++) {
            this->updateC_solveH(i);
            T* Bptr = this->B_solveH[i].get();
            arma::mat* Hptr = this->Hi[i].get();

            unsigned int numChunks = this->ncol_E[i] / ONE_THREAD_MATRIX_SIZE;
            if (numChunks * ONE_THREAD_MATRIX_SIZE < this->ncol_E[i]) numChunks++;
#pragma omp parallel for schedule(auto)
            for (unsigned int i = 0; i < numChunks; i++) {
                unsigned int spanStart = i * ONE_THREAD_MATRIX_SIZE;
                unsigned int spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
                if (spanEnd > this->ncol_E[i] - 1) {
                    spanEnd = this->ncol_E[i] - 1;
                }
                BPPNNLS<arma::mat, arma::vec> subProbH(this->C_solveH,
                                                       (arma::mat)(*Bptr).cols(spanStart, spanEnd));
                subProbH.solveNNLS();
                (*Hptr).rows(spanStart, spanEnd) = subProbH.getSolutionMatrix().t();
            }
            this->updateC_solveW(i);
        }
    }

    void solveV() {
        // implement
        for (int i=0; i<this->nDatasets; i++) {
            this->updateC_solveV(i);
            unsigned int C_V_rows = 2 * this->E_ncols[i];
            T* Bptr = this->B_solveV->rows(0, C_V_rows).get();
            arma::mat* Vptr = this->Vi[i].get();

            unsigned int numChunks = this->m / ONE_THREAD_MATRIX_SIZE;
            if (numChunks * ONE_THREAD_MATRIX_SIZE < this->m) numChunks++;
#pragma omp parallel for schedule(auto)
            for (unsigned int i = 0; i < numChunks; i++) {
                unsigned int spanStart = i * ONE_THREAD_MATRIX_SIZE;
                unsigned int spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
                if (spanEnd > this->m - 1) {
                    spanEnd = this->m - 1;
                }
                BPPNNLS<arma::mat, arma::vec> subProbV(this->C_solveV->rows(0, C_V_rows),
                                                       (arma::mat)(*Bptr).cols(spanStart, spanEnd));
                subProbV.solveNNLS();
                (*Vptr).rows(spanStart, spanEnd) = subProbV.getSolutionMatrix().t();
            }
            this->updateB_solveW(i);
        }
    }

    void solveW() {
        // implement
        arma::mat* C_solveWptr = this->C_solveW.get();
        T* B_solveWptr = this->B_solveW.get();
        arma::mat* Wptr = this->W.get();
        unsigned int numChunks = this->m / ONE_THREAD_MATRIX_SIZE;
        if (numChunks * ONE_THREAD_MATRIX_SIZE < this->m) numChunks++;
#pragma omp parallel for schedule(auto)
        for (unsigned int i = 0; i < numChunks; i++) {
            unsigned int spanStart = i * ONE_THREAD_MATRIX_SIZE;
            unsigned int spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
            if (spanEnd > this->m - 1) {
                spanEnd = this->m - 1;
            }
            BPPNNLS<arma::mat, arma::vec> subProbW(C_solveWptr,
                                                   (arma::mat)B_solveWptr.cols(spanStart, spanEnd));
            subProbW.solveNNLS();
            (*Wptr).rows(spanStart, spanEnd) = subProbW.getSolutionMatrix().t();
        }
    }

public:
    BPPINMF(std::vector<std::unique_ptr<T>> &Ei, arma::uword k, double lambda) : INMF<T>(Ei, k, lambda) {

    }

    void optimizeALS(int maxIter, const double thresh) {
        // execute private functions here
        unsigned int iter = 0;
        double delta, obj;
        while (iter < maxIter && delta > thresh) {
            solveH();
            solveV();
            solveW();
            obj = this->objective();
            delta = arma::abs<double>(this->objective_err - obj) / ((this->objective_err + obj) / 2);
            iter++;
            this->objective_err = obj;
        }
    }

    arma::mat getHi(int i) {
        return *(this->Hi[i].get());
    }

    arma::mat getHi() {
        return this->Hi;
    }

    arma::mat getVi(int i) {
        return *(this->Vi[i].get());
    }

    arma::mat getVi() {
        return this->Vi;
    }

    arma::mat getW() {
        return *(this->W.get());
    }
};

}
