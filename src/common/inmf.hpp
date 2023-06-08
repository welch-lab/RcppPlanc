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
            arma::mat* Wptr = W.get();
            arma::mat* Vptr = Vi[i].get();
            arma::mat* C_solveHptr = &C_solveH;
            C_solveHptr->rows(0, this->m - 1) = *Wptr + *Vptr;
            C_solveHptr->rows(this->m, 2 * this->m - 1) = sqrtLambda * *Vptr;
            // this->C_solveH[i] = arma::join_cols(*Wptr + *Vptr, sqrtLambda * *Vptr);

        }

        void updateC_solveV(int i) {
            arma::mat* Hptr = Hi[i].get();
            arma::mat* C_solveVptr = &C_solveV;
            C_solveVptr->rows(0, this->ncol_E[i] - 1) = *Hptr;
            C_solveVptr->rows(this->ncol_E[i], 2 * this->ncol_E[i] - 1) = sqrtLambda * *Hptr;
            // this->C_solveV[i] = arma::join_cols(*Hptr, sqrtLambda * *Hptr);
        }

        void updateC_solveW(int i) {
            // Only update slices of C_solveW that correspond to Hi[i]
            arma::mat* C_solveWptr = &C_solveW;
            arma::mat* Hptr = Hi[i].get();
            unsigned int start = 0, end;
            // accumulate the number of columns of Ei up to i
            for (unsigned int j = 0; j < i; ++j) {
                start += this->ncol_E[j];
            }
            end = start + this->ncol_E[i] - 1;
            C_solveWptr->rows(start, end) = *Hptr;
        }

        void updateB_solveH(unsigned int i, T* E) { // call in constructor
            // TODO: Whether to take `i` as argument or directly take Ei[i]
            T B_H(2 * this->m, this->ncol_E[i]);
            // Copy the values from E to the top half of B_H
            B_H.rows(0, this->m - 1) = *E;
            B_solveH.push_back(std::make_unique<T>(B_H));
        }

        void updateB_solveV(int i) {
            T* Eptr = Ei[i].get();
            arma::mat* Wptr = W.get();
            arma::mat* Hptr = Hi[i].get();
            T* B_solveVptr = &B_solveV;
            B_solveVptr->rows(0, this->ncol_E[i] - 1) = Eptr->t() - *Hptr * Wptr->t();
            B_solveVptr->rows(this->ncol_E[i], 2 * this->ncol_E[i] - 1) = arma::zeros<T>(this->ncol_E[i], this->m);
            // this->B_solveV = arma::zeros(2 * this->ncol_E[i], this->m);
            // this->B_solveV.rows(0, this->ncol_E[i] - 1) = Eptr->t() - *Hptr * Wptr->t();
        }

        void updateB_solveW(int i) {
            T* Eptr = Ei[i].get();
            arma::mat* Hptr = Hi[i].get();
            arma::mat* Vptr = Vi[i].get();
            T* B_solveWptr = &B_solveW;
            unsigned int start = 0, end;
            // accumulate the number of columns of Ei up to i
            for (unsigned int j = 0; j < i; ++j) {
                start += this->ncol_E[j];
            }
            end = start + this->ncol_E[i] - 1;
            B_solveWptr->rows(start, end) = Eptr->t() - *Hptr * Vptr->t();
        }

        double objective() {
            // Calculate the objective function
            double obj = 0, n1, n2;
            arma::mat* Wptr = W.get();
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                T* Eptr = Ei[i].get();
                arma::mat* Hptr = Hi[i].get();
                arma::mat* Vptr = Vi[i].get();
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
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                this->Hi.push_back(std::make_unique<arma::mat>(this->ncol_E[i], this->k));
                this->Vi.push_back(std::make_unique<arma::mat>(this->m, this->k));
            }
            this->W = arma::randu<arma::mat>(this->m, this->k);
            this->objective_err = this->objective();
        }
    public:
        INMF(std::vector<std::unique_ptr<T>>& Ei,
            arma::uword k, double lambda) {
            this->k = k;
            this->m = Ei[0].get()->n_rows;
            this->nMax = 0;
            this->nSum = 0;
            for (unsigned int i = 0; i < Ei.size(); ++i)
            {
                T* E = Ei[i].get();
                this->ncol_E.push_back(E->n_cols);
                if (E->n_cols > this->nMax) {
                    this->nMax = E->n_cols;
                }
                this->nSum += E->n_cols;
                this->nDatasets++;
                updateB_solveH(i, E);
                // push random initialized H, V, W according to the size of Ei

            };
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
