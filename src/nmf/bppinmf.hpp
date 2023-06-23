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
        auto time0 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapse;
        std::cout << "--Solving H--  ";
        for (int i=0; i<this->nDatasets; ++i) {
            this->updateC_solveH(i);
            arma::mat* Hptr = this->Hi[i].get();
            T* Eptr = this->Ei[i].get();
            unsigned int dataSize = this->ncol_E[i];
            unsigned int numChunks = dataSize / ONE_THREAD_MATRIX_SIZE;
            if (numChunks * ONE_THREAD_MATRIX_SIZE < dataSize) numChunks++;
#pragma omp parallel for schedule(auto)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * ONE_THREAD_MATRIX_SIZE;
                unsigned int spanEnd = (j + 1) * ONE_THREAD_MATRIX_SIZE - 1;
                if (spanEnd > dataSize - 1) spanEnd = dataSize - 1;
                arma::mat B_solveH_chunk = arma::zeros<arma::mat>(2 * this->m, spanEnd - spanStart + 1);
                B_solveH_chunk.rows(0, this->m - 1) = Eptr->cols(spanStart, spanEnd);
                // BPPNNLS<arma::mat, arma::vec> subProbH(this->C_solveH,
                //                                        (arma::mat)(*Bptr).cols(spanStart, spanEnd));
                BPPNNLS<arma::mat, arma::vec> subProbH(this->C_solveH, B_solveH_chunk);
                subProbH.solveNNLS();
                (*Hptr).rows(spanStart, spanEnd) = subProbH.getSolutionMatrix().t();
            }
            this->updateC_solveW(i);
        }
        auto time1 = std::chrono::high_resolution_clock::now();
        elapse = time1 - time0;
        std::cout << elapse.count() << " sec" << std::endl;
    }

    void solveV() {
        auto time0 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapse;
        std::cout << "--Solving V--  ";
        arma::mat* WTptr = this->WT.get();
        for (int i=0; i<this->nDatasets; ++i) {
            this->updateC_solveV(i);
            // this->updateB_solveV(i);
            unsigned int C_V_rows = 2 * this->ncol_E[i];
            // TODO: still creating new matrices here
            // arma::mat B = this->B_solveV.rows(0, C_V_rows - 1);
            arma::mat C = this->C_solveV.rows(0, C_V_rows - 1);
            arma::mat* Vptr = this->Vi[i].get();
            arma::mat* VTptr = this->ViT[i].get();
            arma::mat* Hptr = this->Hi[i].get();
            T* ETptr = this->EiT[i].get();
            unsigned int numChunks = this->m / ONE_THREAD_MATRIX_SIZE;
            if (numChunks * ONE_THREAD_MATRIX_SIZE < this->m) numChunks++;
#pragma omp parallel for schedule(auto)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * ONE_THREAD_MATRIX_SIZE;
                unsigned int spanEnd = (j + 1) * ONE_THREAD_MATRIX_SIZE - 1;
                if (spanEnd > this->m - 1) spanEnd = this->m - 1;
                arma::mat B = arma::zeros<arma::mat>(C_V_rows, spanEnd - spanStart + 1);
                B.rows(0, this->ncol_E[i] - 1) = *Hptr * (*WTptr).cols(spanStart, spanEnd);
                B.rows(0, this->ncol_E[i] - 1) -= (*ETptr).cols(spanStart, spanEnd);
                B *= -1;
                BPPNNLS<arma::mat, arma::vec> subProbV(C, B);
                subProbV.solveNNLS();
                (*Vptr).rows(spanStart, spanEnd) = subProbV.getSolutionMatrix().t();
            }
            *VTptr = (*Vptr).t();
            // this->updateB_solveW(i);
        }
        auto time1 = std::chrono::high_resolution_clock::now();
        elapse = time1 - time0;
        std::cout << elapse.count() << " sec" << std::endl;
    }

    void solveW() {
        auto time0 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapse;
        std::cout << "--Solving W--  ";
        arma::mat* Wptr = this->W.get();
        arma::mat* WTptr = this->WT.get();
        unsigned int numChunks = this->m / ONE_THREAD_MATRIX_SIZE;
        if (numChunks * ONE_THREAD_MATRIX_SIZE < this->m) numChunks++;
#pragma omp parallel for schedule(auto)
        for (unsigned int i = 0; i < numChunks; ++i) {
            unsigned int spanStart = i * ONE_THREAD_MATRIX_SIZE;
            unsigned int spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
            if (spanEnd > this->m - 1) spanEnd = this->m - 1;
            // See inmf.hpp protected members for dimensionality details
            arma::mat B = arma::zeros<arma::mat>(this->nSum, spanEnd - spanStart + 1);
            arma::uword rowStart = 0;
            arma::uword rowEnd = 0;
            for (unsigned int j = 0; j < this->nDatasets; ++j) {
                T* ETptr = this->EiT[j].get();
                arma::mat* Hptr = this->Hi[j].get();
                arma::mat* VTptr = this->ViT[j].get();
                rowEnd = rowStart + this->ncol_E[j] - 1;
                B.rows(rowStart, rowEnd) = *Hptr * (*VTptr).cols(spanStart, spanEnd);
                B.rows(rowStart, rowEnd) -= (*ETptr).cols(spanStart, spanEnd);
                rowStart = rowEnd + 1;
            }
            B *= -1;
            BPPNNLS<arma::mat, arma::vec> subProbW(this->C_solveW, B);
            subProbW.solveNNLS();
            (*Wptr).rows(spanStart, spanEnd) = subProbW.getSolutionMatrix().t();
        }
        *WTptr = (*Wptr).t();
        auto time1 = std::chrono::high_resolution_clock::now();
        elapse = time1 - time0;
        std::cout << elapse.count() << " sec" << std::endl;
    }

public:
    BPPINMF(std::vector<std::unique_ptr<T>>& Ei, arma::uword k, double lambda) : INMF<T>(Ei, k, lambda) {

    }

    void optimizeALS(int maxIter, const double thresh) {
        // execute private functions here
        std::cout << "BPPINMF optimizeALS started, maxIter=" << maxIter << ", thresh=" << thresh << std::endl;
        unsigned int iter = 0;
        double delta=100, obj;
        std::chrono::duration<double> elapse;

        while (delta > thresh && iter < maxIter ) {
            auto time0 = std::chrono::high_resolution_clock::now();
            std::cout << "========Staring iteration " << iter+1 << "========" << std::endl;
            solveH();
            solveV();
            solveW();
            obj = this->objective();
            std::cout << "Iteration: " << iter << " Objective: " << obj << std::endl;
            delta = abs(this->objective_err - obj) / ((this->objective_err + obj) / 2);
            std::cout << "Delta: " << delta << std::endl;
            iter++;
            this->objective_err = obj;
            auto time1 = std::chrono::high_resolution_clock::now();
            elapse = time1 - time0;
            std::cout << "Iteration " << iter << " took " << elapse.count() << " seconds" << std::endl;
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
