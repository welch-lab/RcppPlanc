#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <memory>
#include "utils.hpp"

#define INMF_CHUNK_SIZE 1000

namespace planc {

    template <typename T>
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
                unsigned int numChunks = dataSize / INMF_CHUNK_SIZE;
                if (numChunks * INMF_CHUNK_SIZE < dataSize) numChunks++;
#pragma omp parallel for schedule(auto)
                for (unsigned int j = 0; j < numChunks; ++j) {
                    unsigned int spanStart = j * INMF_CHUNK_SIZE;
                    unsigned int spanEnd = (j + 1) * INMF_CHUNK_SIZE - 1;
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

        void constructObject(std::vector<std::unique_ptr<T>>& Ei, arma::uword k, double lambda, bool makeTranspose) {
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
                if (makeTranspose) {
                    T ET = E->t();
                    std::unique_ptr<T> ETptr = std::make_unique<T>(ET);
                    this->EiT.push_back(std::move(ETptr));
                }
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

        void checkK() {
            for (unsigned int i = 0; i < this->nDatasets; ++i) {
                if (this->k != this->Hi[i]->n_cols) {
                    std::string msg = "Preset `k` (" + std::to_string(this->k) +
                                      ") does not match with H[" + std::to_string(i) + "]";
                    throw std::invalid_argument(msg);
                }
                if (this->k != this->Vi[i]->n_cols) {
                    std::string msg = "Preset `k` (" + std::to_string(this->k) +
                                      ") does not match with V[" + std::to_string(i) + "]";
                    throw std::invalid_argument(msg);
                }
            }
            if (this->k != this->W.get()->n_cols) {
                std::string msg = "Preset `k` (" + std::to_string(this->k) +
                                  ") does not match with W";
                throw std::invalid_argument(msg);
            }
        }
    public:
        INMF(std::vector<std::unique_ptr<T>>& Ei, arma::uword k, double lambda, bool makeTranspose = true) {
            this->constructObject(Ei, k, lambda, makeTranspose);
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
            if (this->Hi.size() != this->nDatasets) {
                std::string msg = "Must provide " +
                                  std::to_string(this->nDatasets) +
                                  " H matrices";
                throw std::invalid_argument(msg);
            }
            std::unique_ptr<arma::mat> H;
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                if (Hinit[i].n_cols != this->k || Hinit[i].n_rows != this->ncol_E[i]) {
                    std::string msg = "Each given H must be of size Ei[i].n_cols x " +
                                      std::to_string(this->k);
                    throw std::invalid_argument(msg);
                }
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

        void initV(std::vector<arma::mat>& Vinit, bool transpose = true) {
#ifdef _VERBOSE
            std::cout << "Taking initialized V matrices" << std::endl;
#endif
            if (Vinit.size() != this->nDatasets) {
                std::string msg = "Must provide " +
                                  std::to_string(this->nDatasets) +
                                  " V matrices";
                throw std::invalid_argument(msg);
            }
            std::unique_ptr<arma::mat> V;
            std::unique_ptr<arma::mat> VT;
            for (arma::uword i = 0; i < this->nDatasets; ++i) {
                if (Vinit[i].n_cols != this->k || Vinit[i].n_rows != this->m) {
                    std::string msg = "All given Vs must be of size " +
                                      std::to_string(this->m) + " x " +
                                      std::to_string(this->k);
                    throw std::invalid_argument(msg);
                }
                V = std::unique_ptr<arma::mat>(new arma::mat);
                *V = Vinit[i];
                this->Vi.push_back(std::move(V));
                if (transpose) {
                    VT = std::unique_ptr<arma::mat>(new arma::mat);
                    *VT = (*V).t();
                    this->ViT.push_back(std::move(VT));
                }
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

        void initW(arma::mat Winit, bool transpose = true) {
#ifdef _VERBOSE
            std::cout << "Taking initialized W matrix" << std::endl;
#endif
            if (Winit.n_cols != this->k || Winit.n_rows != this->m) {
                std::string msg = "Given W must be of size " +
                                  std::to_string(this->m) + " x " +
                                  std::to_string(this->k) + " but is " +
                                  std::to_string(Winit.n_rows) + " x " +
                                  std::to_string(Winit.n_cols);
                throw std::invalid_argument(msg);
            }
            this->W = std::unique_ptr<arma::mat>(new arma::mat);
            *this->W = Winit;
            if (transpose) {
                this->WT = std::unique_ptr<arma::mat>(new arma::mat);
                *this->WT = (*this->W).t();
            }
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

        ~INMF() { clear(); }
        void clear() {
            if (!this->cleared) {
                for (unsigned int i = 0; i < Ei.size(); ++i) {
                    Ei[i].reset();
                    // EiT[i].reset();
                }
                for (unsigned int i = 0; i < EiT.size(); ++i) {
                    EiT[i].reset();
                }
                for (unsigned int i = 0; i < Hi.size(); ++i) {
                    Hi[i].reset();
                }
                for (unsigned int i = 0; i < Vi.size(); ++i) {
                    Vi[i].reset();
                }
                for (unsigned int i = 0; i < ViT.size(); ++i) {
                    ViT[i].reset();
                }
                this->W.reset();
                if (this->WT.get() != nullptr) this->WT.reset();
            }
            this->cleared = true;
        }
    };

}
