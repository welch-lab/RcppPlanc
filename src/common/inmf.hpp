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
        std::vector<std::unique_ptr<T>> EiT;          // each of size n_ixm
        std::vector<std::unique_ptr<arma::mat>> Hi;  // each of size n_ixk
        std::vector<std::unique_ptr<arma::mat>> Vi;  // each of size mxk
        std::vector<std::unique_ptr<arma::mat>> ViT; // each of size kxm
        std::unique_ptr<arma::mat> W;                // mxk
        std::unique_ptr<arma::mat> WT;                // kxm
        double lambda, sqrtLambda, objective_err;
        double m_symm_reg;              /// Symmetric Regularization parameter
        arma::fvec m_regW;
        arma::fvec m_regH;
        bool cleared;
        // std::vector<arma::mat> C_solveH;//(2*m, k);
        // arma::mat B_solveH;//(2*m, n_i);
        // std::vector<arma::mat> C_solveV;//(2*n_i, k);
        // std::vector<std::unique_ptr<T>> B_solveV; //(2*n_i, m);
        // arma::mat C_solveW;
        // T B_solveW;
        // arma::mat C_solveH; //(2*m, k);
        // std::vector<std::unique_ptr<T>> B_solveH;    //(2*m, n_i);
        // arma::mat B_solveV;         //(2*n_max, m);
        // arma::mat C_solveW; //(nSum, k)
        // arma::mat B_solveW;         //(nSum, m)
        // arma::mat B_solveH_chunk; //(2*m, CHUNK_SIZE);
// #ifdef CMAKE_BUILD_SPARSE
//         arma::sp_mat B_solveW_i; //(this->n, this->m);
// #else
//         arma::mat B_solveW_i;   //(this->n, this->m); /// At(nxm) - HV(nxm);
// #endif
        void regW(const arma::fvec &iregW) { this->m_regW = iregW; }
        /// Sets the regularization on right low rank H
        void regH(const arma::fvec &iregH) { this->m_regH = iregH; }
        /// Returns the L2 and L1 regularization parameters of W as a vector
        arma::fvec regW() { return this->m_regW; }
        /// Returns the L2 and L1 regularization parameters of W as a vector
        arma::fvec regH() { return this->m_regH; }

        void applyReg(const arma::fvec &reg, arma::mat *AtA) {
            // Frobenius norm regularization
            if (reg(0) > 0) {
            arma::mat identity = arma::eye<arma::mat>(this->k, this->k);
            float lambda_l2 = reg(0);
            (*AtA) = (*AtA) + 2 * lambda_l2 * identity;
            }

            // L1 - norm regularization
            if (reg(1) > 0) {
            arma::mat onematrix = arma::ones<arma::mat>(this->k, this->k);
            float lambda_l1 = reg(1);
            (*AtA) = (*AtA) + 2 * lambda_l1 * onematrix;
            }
        }

        void applySymmetricReg(double sym_reg, arma::mat *lhs, arma::mat *fac, arma::mat *rhs) {
            if (sym_reg > 0) {
            arma::mat identity = arma::eye<arma::mat>(this->k, this->k);
            (*lhs) = (*lhs) + sym_reg * identity;
            (*rhs) = (*rhs) + sym_reg * (*fac);
            }
        }
        void symm_reg(const double &i_symm_reg) { this->m_symm_reg = i_symm_reg; }
        /// Returns the Symmetric regularization parameter
        double symm_reg() { return this->m_symm_reg; }

        // void updateC_solveH(int i) {
        //     // Execute in constructor and after each update of W
        //     arma::mat* Wptr = this->W.get();
        //     arma::mat* Vptr = this->Vi[i].get();
        //     this->C_solveH.rows(0, this->m - 1) = *Wptr + *Vptr;
        //     this->C_solveH.rows(this->m, 2 * this->m - 1) = sqrtLambda * *Vptr;
        // }

        // void updateC_solveW(int i) {
        //     // Only update slices of C_solveW that correspond to Hi[i]
        //     arma::mat* Hptr = this->Hi[i].get();
        //     unsigned int start = 0, end;
        //     // accumulate the number of columns of Ei up to i
        //     for (unsigned int j = 0; j < i; ++j) {
        //         start += this->ncol_E[j];
        //     }
        //     end = start + this->ncol_E[i] - 1;
        //     this->C_solveW.rows(start, end) = *Hptr;
        // }

        // void updateB_solveH(T* E) {
        //     std::unique_ptr<T> B_H;
        //     B_H = std::unique_ptr<T>(new T(2 * this->m, E->n_cols));
        //     // Copy the values from E to the top half of B_H
        //     B_H->rows(0, this->m - 1) = *E;
        //     B_solveH.push_back(std::move(B_H));
        // }

        double computeObjectiveError() {
            auto time0 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapse;
            // Calculate the objective function
            double obj = 0, norm;
            arma::mat* Wptr = this->W.get();
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                arma::uword dataSize = this->ncol_E[i];
                T* Eptr = this->Ei[i].get();
                arma::mat* Hptr = this->Hi[i].get();
                arma::mat* Vptr = this->Vi[i].get();
                unsigned int numChunks = dataSize / ONE_THREAD_MATRIX_SIZE;
                if (numChunks * ONE_THREAD_MATRIX_SIZE < dataSize) numChunks++;
#pragma omp parallel for schedule(auto)
                for (unsigned int j = 0; j < numChunks; ++j) {
                    arma::mat errMat(this->m, ONE_THREAD_MATRIX_SIZE);

                    unsigned int spanStart = j * ONE_THREAD_MATRIX_SIZE;
                    unsigned int spanEnd = (j + 1) * ONE_THREAD_MATRIX_SIZE - 1;
                    if (spanEnd > dataSize - 1) spanEnd = dataSize - 1;
                    if (j == numChunks - 1) errMat.resize(this->m, spanEnd - spanStart + 1);
                    errMat = (*Vptr + *Wptr) * arma::mat(Hptr->t()).cols(spanStart, spanEnd);
                    errMat -= Eptr->cols(spanStart, spanEnd);
                    norm = arma::norm<arma::mat>(errMat, "fro");
                    norm *= norm;
                    obj += norm;

                    errMat = *Vptr * arma::mat(Hptr->t()).cols(spanStart, spanEnd);
                    norm = arma::norm<arma::mat>(errMat, "fro");
                    norm *= norm;
                    obj += this->lambda * norm;
                }
            }
            auto time1 = std::chrono::high_resolution_clock::now();
            elapse = time1 - time0;
            std::cout << "Objective calculation time: " << elapse.count() << " sec" << std::endl;
            return obj;
        }

        // void initHWV() {
        //     // Initialize H, W, V
        //     std::unique_ptr<arma::mat> H;
        //     std::unique_ptr<arma::mat> V;
        //     std::unique_ptr<arma::mat> VT;
        //     // std::unique_ptr<arma::mat> W;
        //     for (arma::uword i = 0; i < this->nDatasets; ++i) {
        //         H = std::unique_ptr<arma::mat>(new arma::mat);
        //         *H = arma::randu<arma::mat>(this->ncol_E[i], this->k);
        //         this->Hi.push_back(std::move(H));

        //         V = std::unique_ptr<arma::mat>(new arma::mat);
        //         VT = std::unique_ptr<arma::mat>(new arma::mat);
        //         *V = arma::randu<arma::mat>(this->m, this->k);
        //         *VT = (*V).t();
        //         this->Vi.push_back(std::move(V));
        //         this->ViT.push_back(std::move(VT));
        //     }
        //     this->W = std::unique_ptr<arma::mat>(new arma::mat);
        //     this->WT = std::unique_ptr<arma::mat>(new arma::mat);
        //     *this->W = arma::randu<arma::mat>(this->m, this->k);
        //     *this->WT = (*this->W).t();
        //     this->objective_err = this->computeObjectiveError();
        //     std::cout << "Initial objective: " << this->objective_err << std::endl;
        // }

        void constructObject(std::vector<std::unique_ptr<T>>& Ei, arma::uword k, double lambda) {
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
                T ET = E->t();
                std::unique_ptr<T> ETptr = std::make_unique<T>(ET);
                this->EiT.push_back(std::move(ETptr));
                this->ncol_E.push_back(E->n_cols);
                if (E->n_cols > this->nMax) {
                    this->nMax = E->n_cols;
                }
                this->nSum += E->n_cols;
                this->nDatasets++;
                // updateB_solveH(E);
            };
            std::cout << "nMax=" << this->nMax << "; nSum=" << this->nSum << std::endl;
            std::cout << "nDatasets=" << this->nDatasets << std::endl;
            // this->initHWV();
            this->lambda = lambda;
            this->sqrtLambda = sqrt(lambda); //TODO
            this->m_regW = arma::zeros<arma::fvec>(2);
            this->m_regH = arma::zeros<arma::fvec>(2);
            //TODO implement common tasks i.e. norm, reg, etc
        }
    public:
        INMF(std::vector<std::unique_ptr<T>>& Ei, arma::uword k, double lambda) {
            this->constructObject(Ei, k, lambda);
            // this->initHWV();
        }
        // INMF(std::vector<std::unique_ptr<T>>& Ei, arma::uword k, double lambda,
        //      std::vector<std::unique_ptr<arma::mat>>& Hinit,
        //      std::vector<std::unique_ptr<arma::mat>>& Vinit,
        //      arma::mat& Winit) {
        //     this->constructObject(Ei, k, lambda);
        //     this->W = std::make_unique<arma::mat>(Winit);

        //     for (arma::uword i = 0; i < this->nDatasets; ++i) {
        //         this->Hi.push_back(std::move(Hinit[i]));
        //         arma::mat* Vptr = Vinit[i].get();
        //         std::unique_ptr<arma::mat> VTptr = std::unique_ptr<arma::mat>(new arma::mat);
        //         *VTptr = Vptr->t();
        //         this->Vi.push_back(std::move(Vinit[i]));
        //         this->ViT.push_back(std::move(VTptr));
        //     }
        //     std::unique_ptr<arma::mat> WTptr = std::unique_ptr<arma::mat>(new arma::mat);
        //     *WTptr = this->W->t();
        //     this->WT = std::move(WTptr);
        // }
        void initH(std::vector<arma::mat>& Hinit) {
            std::cout << "Taking initialized H matrices" << std::endl;
            std::unique_ptr<arma::mat> H;
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                H = std::unique_ptr<arma::mat>(new arma::mat);
                *H = Hinit[i];
                this->Hi.push_back(std::move(H));
            }
        }

        void initH() {
            std::cout << "Randomly initializing H matrices" << std::endl;
            std::unique_ptr<arma::mat> H;
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                H = std::unique_ptr<arma::mat>(new arma::mat);
                *H = arma::randu<arma::mat>(this->ncol_E[i], this->k);
                this->Hi.push_back(std::move(H));
            }
        }

        void initV(std::vector<arma::mat>& Vinit) {
            std::cout << "Taking initialized V matrices" << std::endl;
            std::unique_ptr<arma::mat> V;
            std::unique_ptr<arma::mat> VT;
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                V = std::unique_ptr<arma::mat>(new arma::mat);
                VT = std::unique_ptr<arma::mat>(new arma::mat);
                *V = Vinit[i];
                *VT = (*V).t();
                this->Vi.push_back(std::move(V));
                this->ViT.push_back(std::move(VT));
            }
        }

        void initV() {
            std::cout << "Randomly initializing V matrices" << std::endl;
            std::unique_ptr<arma::mat> V;
            std::unique_ptr<arma::mat> VT;
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                V = std::unique_ptr<arma::mat>(new arma::mat);
                VT = std::unique_ptr<arma::mat>(new arma::mat);
                *V = arma::randu<arma::mat>(this->m, this->k);
                *VT = (*V).t();
                this->Vi.push_back(std::move(V));
                this->ViT.push_back(std::move(VT));
            }
        }

        void initW(arma::mat Winit) {
            std::cout << "Taking initialized W matrix" << std::endl;
            this->W = std::unique_ptr<arma::mat>(new arma::mat);
            this->WT = std::unique_ptr<arma::mat>(new arma::mat);
            *this->W = Winit;
            *this->WT = (*this->W).t();
        }

        void initW() {
            std::cout << "Randomly initializing W matrix" << std::endl;
            this->W = std::unique_ptr<arma::mat>(new arma::mat);
            this->WT = std::unique_ptr<arma::mat>(new arma::mat);
            *this->W = arma::randu<arma::mat>(this->m, this->k);
            *this->WT = (*this->W).t();
        }

        ~INMF() { clear(); }
        void clear() {
            if (!this->cleared) {
                // this->C_solveH.clear();
                // this->C_solveW.clear();
                // this->B_solveW.clear();
                // this->B_solveV.clear();
                // this->B_solveH_chunk.clear();
                // for (unsigned int i = 0; i < B_solveH.size(); ++i) {
                //     B_solveH[i].reset();
                // }
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
