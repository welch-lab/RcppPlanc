#pragma once
/* Copyright 2016 Ramakrishnan Kannan */

#include "nmf.hpp"

namespace planc {

template <class T>
class MUNMF : public NMF<T> {
 private:
  // Not happy with this design. However to avoid computing At again and again
  // making this as private variable.
  T At;
  arma::mat WtW;
  arma::mat HtH;
  arma::mat AtW;
  arma::mat AH;

  /*
   * Collected statistics are
   * iteration Htime Wtime totaltime normH normW densityH densityW relError
   */
  void allocateMatrices() {
    WtW = arma::zeros<arma::mat>(this->k, this->k);
    HtH = arma::zeros<arma::mat>(this->k, this->k);
    AtW = arma::zeros<arma::mat>(this->n, this->k);
    AH = arma::zeros<arma::mat>(this->m, this->k);
  }
  void freeMatrices() {
    this->At.clear();
    WtW.clear();
    HtH.clear();
    AtW.clear();
    AH.clear();
  }

 public:
  MUNMF(const T &A, int lowrank) : NMF<T>(A, lowrank) {
    allocateMatrices();
    this->At = this->A.t();
  }
  MUNMF(const T &A, const arma::mat &llf, const arma::mat &rlf) : NMF<T>(A, llf, rlf) {
    allocateMatrices();
    this->At = this->A.t();
  }
  void computeNMF() {
    unsigned int currentIteration = 0;
    INFO << "computed transpose At=" << PRINTMATINFO(this->At) << std::endl;
    while (currentIteration < this->num_iterations()) {
      tic();
      // update H
      tic();
      AtW = this->At * this->W;
      WtW = this->W.t() * this->W;
      this->applyReg(this->regH(), &this->WtW);
      INFO << "starting H Prereq for "
           << " took=" << toc();
      INFO << PRINTMATINFO(WtW) << PRINTMATINFO(AtW) << std::endl;
      // to avoid divide by zero error.
      tic();
      // H = H.*AtW./(WtW_reg*H + epsilon);
      this->H = (this->H % AtW) / (this->H * WtW + EPSILON_1EMINUS16);
      INFO << "Completed H (" << currentIteration << "/"
           << this->num_iterations() << ")"
           << " time =" << toc() << std::endl;

      // update W;
      tic();
      AH = this->A * this->H;
      HtH = this->H.t() * this->H;
      this->applyReg(this->regW(), &this->HtH);
      INFO << "starting W Prereq for "
           << " took=" << toc() << PRINTMATINFO(HtH) << PRINTMATINFO(AH)
           << std::endl;
      tic();
      // W = W.*AH./(W*HtH_reg + epsilon);
      this->W = (this->W % AH) / ((this->W * HtH) + EPSILON_1EMINUS16);
      INFO << "Completed W (" << currentIteration << "/"
           << this->num_iterations() << ")"
           << " time =" << toc() << std::endl;
      INFO << "Completed It (" << currentIteration << "/"
           << this->num_iterations() << ")"
           << " time =" << toc() << std::endl;
      this->computeObjectiveError();
      INFO << "Completed it = " << currentIteration
           << " MUERR=" << sqrt(this->objective_err) / this->normA << std::endl;
      currentIteration++;
    }
    this->normalize_by_W();
  }
  ~MUNMF() { freeMatrices(); }
};

}  // namespace planc
