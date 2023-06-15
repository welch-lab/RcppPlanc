#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <memory>
#include "utils.hpp"

namespace planc {

    template <class T>
    class INMF {
    protected:
        arma::uword m, k, nDatasets, nMax, nSum;
        std::vector<arma::uword> ncol_E;             // vector of n_i
        std::vector<std::unique_ptr<T>> Ei;          // each of size mxn_i
        std::vector<std::unique_ptr<arma::mat>> Hi;  // each of size n_ixk
        std::vector<std::unique_ptr<arma::mat>> Vi;  // each of size mxk
        std::unique_ptr<arma::mat> W;                // mxk
        double lambda, sqrtLambda, objective_err;
        bool cleared;
        // std::vector<arma::mat> C_solveH;//(2*m, k);
        // T B_solveH;//(2*m, n_i);
        // std::vector<arma::mat> C_solveV;//(2*n_i, k);
        // std::vector<std::unique_ptr<T>> B_solveV; //(2*n_i, m);
        // arma::mat C_solveW;
        // T B_solveW;
        arma::mat C_solveH; //(2*m, k);
        std::vector<std::unique_ptr<T>> B_solveH;    //(2*m, n_i);
        arma::mat C_solveV; //(2*n_max, k);
        T B_solveV;         //(2*n_max, m);
        arma::mat C_solveW; //(nSum, k)
        T B_solveW;         //(nSum, m)

        // arma::mat HV;//(this->n, this->m); /// H(nxk) * Vt(kxm);
// #ifdef CMAKE_BUILD_SPARSE
//         arma::sp_mat B_solveW_i; //(this->n, this->m);
// #else
//         arma::mat B_solveW_i;   //(this->n, this->m); /// At(nxm) - HV(nxm);
// #endif

        void updateC_solveH(int i) {
            // Execute in constructor and after each update of W
            arma::mat* Wptr = this->W.get();
            arma::mat* Vptr = this->Vi[i].get();
            this->C_solveH.rows(0, this->m - 1) = *Wptr + *Vptr;
            this->C_solveH.rows(this->m, 2 * this->m - 1) = sqrtLambda * *Vptr;
        }

        void updateC_solveV(int i) {
            arma::mat* Hptr = this->Hi[i].get();
            this->C_solveV.rows(0, this->ncol_E[i] - 1) = *Hptr;
            this->C_solveV.rows(this->ncol_E[i], 2 * this->ncol_E[i] - 1) = sqrtLambda * *Hptr;
        }

        void updateC_solveW(int i) {
            // Only update slices of C_solveW that correspond to Hi[i]
            arma::mat* Hptr = this->Hi[i].get();
            unsigned int start = 0, end;
            // accumulate the number of columns of Ei up to i
            for (unsigned int j = 0; j < i; ++j) {
                start += this->ncol_E[j];
            }
            end = start + this->ncol_E[i] - 1;
            this->C_solveW.rows(start, end) = *Hptr;
        }

        void updateB_solveH(T* E) {
            std::unique_ptr<T> B_H;
            B_H = std::unique_ptr<T>(new T(2 * this->m, E->n_cols));
            // Copy the values from E to the top half of B_H
            B_H->rows(0, this->m - 1) = *E;
            B_solveH.push_back(std::move(B_H));
        }

        void updateB_solveV(int i) {
            T* Eptr = this->Ei[i].get();
            arma::mat* Wptr = W.get();
            arma::mat* Hptr = this->Hi[i].get();
            unsigned int dataSize = this->ncol_E[i];
            this->B_solveV.rows(0, dataSize - 1) = Eptr->t() - *Hptr * Wptr->t();
            this->B_solveV.rows(dataSize, 2 * dataSize - 1) = arma::zeros<T>(dataSize, this->m);
        }

        void updateB_solveW(int i) {
            T* Eptr = this->Ei[i].get();
            arma::mat* Hptr = this->Hi[i].get();
            arma::mat* Vptr = this->Vi[i].get();
            unsigned int start = 0, end;
            // accumulate the number of columns of Ei up to i
            for (unsigned int j = 0; j < i; ++j) {
                start += this->ncol_E[j];
            }
            end = start + this->ncol_E[i] - 1;
            this->B_solveW.rows(start, end) = Eptr->t() - *Hptr * Vptr->t();
        }

        double objective() {
            // Calculate the objective function
            double obj = 0, n1, n2;
            arma::mat* Wptr = W.get();
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                T* Eptr = this->Ei[i].get();
                arma::mat* Hptr = this->Hi[i].get();
                arma::mat* Vptr = this->Vi[i].get();
                arma::mat diff = Eptr->t() - *Hptr * (Vptr->t() + Wptr->t());
                n1 = arma::norm<arma::mat>(diff, "fro");
                n1 *= n1;
                obj += n1;
                n2 = arma::norm<arma::mat>(*Hptr * Vptr->t(), "fro");
                n2 *= n2;
                obj += this->lambda * n2;
            }
            return obj;
        }

        void initHWV() {
            // Initialize H, W, V
            std::unique_ptr<arma::mat> H;
            std::unique_ptr<arma::mat> V;
            // std::unique_ptr<arma::mat> W;
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                H = std::unique_ptr<arma::mat>(new arma::mat);
                *H = arma::randu<arma::mat>(this->ncol_E[i], this->k);
                this->Hi.push_back(std::move(H));

                V = std::unique_ptr<arma::mat>(new arma::mat);
                *V = arma::randu<arma::mat>(this->m, this->k);
                this->Vi.push_back(std::move(V));
            }
            this->W = std::unique_ptr<arma::mat>(new arma::mat);
            *this->W = arma::randu<arma::mat>(this->m, this->k);
            this->objective_err = this->objective();
            std::cout << "Initial objective: " << this->objective_err << std::endl;
        }
    public:
        INMF(std::vector<std::unique_ptr<T>>& Ei,
            arma::uword k, double lambda) {
            this->Ei = std::move(Ei);
            this->k = k;
            this->m = this->Ei[0].get()->n_rows;
            this->nMax = 0;
            this->nSum = 0;
            this->nDatasets = 0;
            std::cout << "k=" << k << "; m=" << m << std::endl;
            for (unsigned int i = 0; i < this->Ei.size(); ++i)
            {
                T* E = this->Ei[i].get();
                this->ncol_E.push_back(E->n_cols);
                if (E->n_cols > this->nMax) {
                    this->nMax = E->n_cols;
                }
                this->nSum += E->n_cols;
                this->nDatasets++;
                updateB_solveH(E);
            };
            std::cout << "nMax=" << this->nMax << "; nSum=" << this->nSum << std::endl;
            std::cout << "nDatasets=" << this->nDatasets << std::endl;
            this->initHWV();
            this->lambda = lambda;
            this->sqrtLambda = sqrt(lambda); //TODO
            //TODO implement common tasks i.e. norm, reg, etc
            this->C_solveH = arma::zeros(2 * this->m, this->k);
            this->C_solveV = arma::zeros(2 * this->nMax, this->k);
            this->C_solveW = arma::zeros(this->nSum, this->k);
            this->B_solveV = arma::zeros<T>(2 * this->nMax, this->m);
            this->B_solveW = arma::zeros<T>(this->nSum, this->m);
        };
        ~INMF() { clear(); }
        void clear() {
            if (!this->cleared) {
                this->C_solveH.clear();
                this->C_solveV.clear();
                this->C_solveW.clear();
                this->B_solveW.clear();
                this->B_solveV.clear();
                for (unsigned int i = 0; i < B_solveH.size(); ++i) {
                    B_solveH[i].reset();
                }
                for (unsigned int i = 0; i < Ei.size(); ++i) {
                    Ei[i].reset();
                }
                for (unsigned int i = 0; i < Hi.size(); ++i) {
                    Hi[i].reset();
                }
                for (unsigned int i = 0; i < Vi.size(); ++i) {
                    Vi[i].reset();
                }
                this->W.reset();
            }
            this->cleared = true;
        }
    };

}
