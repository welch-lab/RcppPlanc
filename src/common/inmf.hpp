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
        bool cleared;

        double computeObjectiveError() {
            tic();
#ifdef _VERBOSE
            std::cout << "--calc  obj--  ";
#endif
            double obj = 0;
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
                    unsigned int spanStart = j * ONE_THREAD_MATRIX_SIZE;
                    unsigned int spanEnd = (j + 1) * ONE_THREAD_MATRIX_SIZE - 1;
                    double norm;
                    if (spanEnd > dataSize - 1) spanEnd = dataSize - 1;
                    arma::mat errMat(this->m, spanEnd - spanStart + 1);
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
#ifdef _VERBOSE
            std::cout << toc() << " sec" << std::endl;
#endif
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
#ifdef _VERBOSE
            std::cout << "k=" << k << "; m=" << m << std::endl;
#endif
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
            };
#ifdef _VERBOSE
            std::cout << "nMax=" << this->nMax << "; nSum=" << this->nSum << std::endl;
            std::cout << "nDatasets=" << this->nDatasets << std::endl;
#endif
            // this->initHWV();
            this->lambda = lambda;
            this->sqrtLambda = sqrt(lambda); //TODO
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
#ifdef _VERBOSE
            std::cout << "Taking initialized H matrices" << std::endl;
#endif
            std::unique_ptr<arma::mat> H;
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                H = std::unique_ptr<arma::mat>(new arma::mat);
                *H = Hinit[i];
                this->Hi.push_back(std::move(H));
            }
        }

        void initH() {
#ifdef _VERBOSE
            std::cout << "Randomly initializing H matrices" << std::endl;
#endif
            std::unique_ptr<arma::mat> H;
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                H = std::unique_ptr<arma::mat>(new arma::mat);
                *H = arma::randu<arma::mat>(this->ncol_E[i], this->k);
                this->Hi.push_back(std::move(H));
            }
        }

        void initV(std::vector<arma::mat>& Vinit) {
#ifdef _VERBOSE
            std::cout << "Taking initialized V matrices" << std::endl;
#endif
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
#ifdef _VERBOSE
            std::cout << "Randomly initializing V matrices" << std::endl;
#endif
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
#ifdef _VERBOSE
            std::cout << "Taking initialized W matrix" << std::endl;
#endif
            this->W = std::unique_ptr<arma::mat>(new arma::mat);
            this->WT = std::unique_ptr<arma::mat>(new arma::mat);
            *this->W = Winit;
            *this->WT = (*this->W).t();
        }

        void initW() {
#ifdef _VERBOSE
            std::cout << "Randomly initializing W matrix" << std::endl;
#endif
            this->W = std::unique_ptr<arma::mat>(new arma::mat);
            this->WT = std::unique_ptr<arma::mat>(new arma::mat);
            *this->W = arma::randu<arma::mat>(this->m, this->k);
            *this->WT = (*this->W).t();
        }

        double objErr() {
            return this->objective_err;
        }

        ~INMF() { clear(); }
        void clear() {
            if (!this->cleared) {
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
