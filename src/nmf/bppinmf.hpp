#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif
#include "bppnnls.hpp"
#include "inmf.hpp"

#define ONE_THREAD_MATRIX_SIZE 2000

namespace planc {

template <class T>
class BPPINMF : public INMF<T> {
private:
    arma::mat giventGiven;

    void solveH() {
        auto time0 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapse;
        std::cout << "--Solving H--  ";
        arma::mat giventInput;
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
            // giventGiven = this->C_solveH.t() * this->C_solveH;
            // giventInput = this->C_solveH.t() * (*Bptr);
            giventInput = given.t() * (*Eptr);
            unsigned int dataSize = this->ncol_E[i];
            unsigned int numChunks = dataSize / ONE_THREAD_MATRIX_SIZE;
            if (numChunks * ONE_THREAD_MATRIX_SIZE < dataSize) numChunks++;
#pragma omp parallel for schedule(auto)
            for (unsigned int j = 0; j < numChunks; ++j) {
                unsigned int spanStart = j * ONE_THREAD_MATRIX_SIZE;
                unsigned int spanEnd = (j + 1) * ONE_THREAD_MATRIX_SIZE - 1;
                if (spanEnd > dataSize - 1) spanEnd = dataSize - 1;
                // B = arma::zeros<arma::mat>(2 * this->m, spanEnd - spanStart + 1);
                // B.rows(0, this->m - 1) = Eptr->cols(spanStart, spanEnd);
                // BPPNNLS<arma::mat, arma::vec> subProbH(this->C_solveH, (arma::mat)Bptr->cols(spanStart, spanEnd));
                BPPNNLS<arma::mat, arma::vec> subProbH(giventGiven, (arma::mat)giventInput.cols(spanStart, spanEnd), true);
                subProbH.solveNNLS();
                (*Hptr).rows(spanStart, spanEnd) = subProbH.getSolutionMatrix().t();
            }
        }
        giventGiven.clear();
        giventInput.clear();
        auto time1 = std::chrono::high_resolution_clock::now();
        elapse = time1 - time0;
        std::cout << elapse.count() << " sec" << std::endl;
    }

    void solveV() {
        auto time0 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapse;
        std::cout << "--Solving V--  ";
        arma::mat* WTptr = this->WT.get();
        arma::mat B;
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
                B = *Hptr * (*WTptr).cols(spanStart, spanEnd);
                B -= (*ETptr).cols(spanStart, spanEnd);
                B *= -1;
                giventInput = (*Hptr).t() * B;
                // BPPNNLS<arma::mat, arma::vec> subProbV(C, B);
                BPPNNLS<arma::mat, arma::vec> subProbV(giventGiven, giventInput, true);
                subProbV.solveNNLS();
                (*Vptr).rows(spanStart, spanEnd) = subProbV.getSolutionMatrix().t();
                (*VTptr).cols(spanStart, spanEnd) = subProbV.getSolutionMatrix();
            }
        }
        B.clear();
        giventGiven.clear();
        giventInput.clear();
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
        arma::mat B;//(this->nSum, ONE_THREAD_MATRIX_SIZE);
        arma::mat giventInput; //(this->k, ONE_THREAD_MATRIX_SIZE); ///
        // this->applyReg(this->regW(), &giventGiven); ///
        giventGiven = arma::zeros<arma::mat>(this->k, this->k); ///
        for (unsigned int j = 0; j < this->nDatasets; ++j) {
            arma::mat* Hptr = this->Hi[j].get();
            giventGiven += (*Hptr).t() * (*Hptr);
        }

        unsigned int numChunks = this->m / ONE_THREAD_MATRIX_SIZE;
        if (numChunks * ONE_THREAD_MATRIX_SIZE < this->m) numChunks++;
#pragma omp parallel for schedule(auto)
        for (unsigned int i = 0; i < numChunks; ++i) {
            unsigned int spanStart = i * ONE_THREAD_MATRIX_SIZE;
            unsigned int spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
            if (spanEnd > this->m - 1) spanEnd = this->m - 1;
            // See inmf.hpp protected members for dimensionality details
            // B = arma::zeros<arma::mat>(this->nSum, spanEnd - spanStart + 1);
            if (i == numChunks - 1) {
                // B.reshape(this->nSum, spanEnd - spanStart + 1);
                giventInput.reshape(this->k, spanEnd - spanStart + 1); ///
            }
            giventInput = arma::zeros<arma::mat>(this->k, spanEnd - spanStart + 1); ///
            for (unsigned int j = 0; j < this->nDatasets; ++j) {
                T* ETptr = this->EiT[j].get();
                arma::mat* Hptr = this->Hi[j].get();
                arma::mat* VTptr = this->ViT[j].get();
                B = *Hptr * (*VTptr).cols(spanStart, spanEnd);
                B -= (*ETptr).cols(spanStart, spanEnd);
                B *= -1;
                giventInput += (*Hptr).t() * B;
            }
            // BPPNNLS<arma::mat, arma::vec> subProbW(this->C_solveW, B);
            BPPNNLS<arma::mat, arma::vec> subProbW(giventGiven, giventInput, true); ///
            subProbW.solveNNLS();
            (*Wptr).rows(spanStart, spanEnd) = subProbW.getSolutionMatrix().t();
            (*WTptr).cols(spanStart, spanEnd) = subProbW.getSolutionMatrix();
        }
        B.clear();
        giventGiven.clear(); ///
        giventInput.clear(); ///
        *WTptr = (*Wptr).t();
        auto time1 = std::chrono::high_resolution_clock::now();
        elapse = time1 - time0;
        std::cout << elapse.count() << " sec" << std::endl;
    }

public:
    BPPINMF(std::vector<std::unique_ptr<T>>& Ei, arma::uword k, double lambda) : INMF<T>(Ei, k, lambda) {

    }
    // BPPINMF(std::vector<std::unique_ptr<T>>& Ei, arma::uword k, double lambda,
    //         std::vector<std::unique_ptr<arma::mat>>& Hinit,
    //         std::vector<std::unique_ptr<arma::mat>>& Vinit,
    //         arma::mat& Winit) : INMF<T>(Ei, k, lambda, Hinit, Vinit, Winit) {

    // }

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
            obj = this->computeObjectiveError();
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
