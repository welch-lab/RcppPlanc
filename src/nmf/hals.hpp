#pragma once
/* Copyright 2016 Ramakrishnan Kannan */

#include "nmf.hpp"

namespace planc {

template <class T>
class HALSNMF : public NMF<T> {
 private:
  // Not happy with this design. However to avoid computing At again and again
  // making this as private variable.
  T At;
  arma::mat WtW;
  arma::mat HtH;
  arma::mat WtA;
  arma::mat AH;

  /*
   * Collected statistics are
   * iteration Htime Wtime totaltime normH normW densityH densityW relError
   */
  void allocateMatrices() {
    WtW = arma::zeros<arma::mat>(this->k, this->k);
    HtH = arma::zeros<arma::mat>(this->k, this->k);
    WtA = arma::zeros<arma::mat>(this->n, this->k);
    AH = arma::zeros<arma::mat>(this->m, this->k);
  }
  void freeMatrices() {
    this->At.clear();
    WtW.clear();
    HtH.clear();
    WtA.clear();
    AH.clear();
  }

 public:
  HALSNMF(const T &A, int lowrank, const int& ncores = 0) : NMF<T>(A, lowrank) {
    this->normalize_by_W();
    allocateMatrices();
    this->At = this->A.t();
  }
  HALSNMF(const T &A, const arma::mat &llf, const arma::mat &rlf, const int& ncores = 0) : NMF<T>(A, llf, rlf) {
    this->normalize_by_W();
    allocateMatrices();
    this->At = this->A.t();
  }
  void computeNMF() {
    unsigned int currentIteration = 0;
    #ifdef _VERBOSE
    INFO << "computed transpose At=" << PRINTMATINFO(this->At) << std::endl;
    #endif
    while (currentIteration < this->num_iterations()) {
      tic();
      // update H
      tic();
      WtA = this->W.t() * this->A;
      WtW = this->W.t() * this->W;
      this->applyReg(this->regH(), &this->WtW);
      #ifdef _VERBOSE
      INFO << "starting H Prereq for "
           << " took=" << toc() << PRINTMATINFO(WtW) << PRINTMATINFO(WtA)
           << std::endl;
      #endif
      // to avoid divide by zero error.
      tic();
      double normConst;
      arma::vec Hx;
      for (unsigned int x = 0; x < this->k; x++) {
        // H(i,:) = max(H(i,:) + WtA(i,:) - WtW_reg(i,:) * H,epsilon);
        Hx = this->H.col(x) + (((WtA.row(x)).t()) - (this->H * (WtW.col(x))));
        fixNumericalError<arma::vec>(&Hx, EPSILON_1EMINUS16, EPSILON_1EMINUS16);
        normConst = norm(Hx);
        if (normConst != 0) {
          this->H.col(x) = Hx;
        }
      }
      #ifdef _VERBOSE
      INFO << "Completed H (" << currentIteration << "/"
           << this->num_iterations() << ")"
           << " time =" << toc() << std::endl;
      #endif
      // update W;
      tic();
      AH = this->A * this->H;
      HtH = this->H.t() * this->H;
      this->applyReg(this->regW(), &this->HtH);
      #ifdef _VERBOSE
      INFO << "starting W Prereq for "
           << " took=" << toc() << PRINTMATINFO(HtH) << PRINTMATINFO(AH)
           << std::endl;
      #endif
      tic();
      arma::vec Wx;
      for (unsigned int x = 0; x < this->k; x++) {
        // arma::fvec Wx = W(:,x) + (AHt(:,x)-W*HHt(:,x))/HHtDiag(x);

        // W(:,i) = max(W(:,i) * HHt_reg(i,i) + AHt(:,i) - W * HHt_reg(:,i),
        //              epsilon);
        Wx = (this->W.col(x) * HtH(x, x)) +
             (((AH.col(x))) - (this->W * (HtH.col(x))));
        fixNumericalError<arma::vec>(&Wx, EPSILON_1EMINUS16, EPSILON_1EMINUS16);
        normConst = norm(Wx);
        if (normConst != 0) {
          Wx = Wx / normConst;
          this->W.col(x) = Wx;
        }
      }
      this->normalize_by_W();
#ifdef _VERBOSE
      INFO << "Completed W (" << currentIteration << "/"
           << this->num_iterations() << ")"
           << " time =" << toc() << std::endl;

      INFO << "Completed It (" << currentIteration << "/"
           << this->num_iterations() << ")"
           << " time =" << toc() << std::endl
#endif
      this->computeObjectiveError();
#ifdef _VERBOSE
      INFO << "Completed it = " << currentIteration
           << " HALSERR=" << sqrt(this->objective_err) / this->normA
           << std::endl;
#endif
      currentIteration++;
    }
  }
  ~HALSNMF() = default;
};

}  // namespace planc
