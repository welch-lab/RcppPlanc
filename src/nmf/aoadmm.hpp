#pragma once
/* Copyright 2016 Ramakrishnan Kannan */

#include "nmf.hpp"

namespace planc {

template <class T>
class AOADMMNMF : public NMF<T> {
 private:
  // Not happy with this design. However to avoid computing At again and again
  // making this as private variable.
  T At;
  arma::mat WtW;
  arma::mat HtH;
  arma::mat WtA;
  arma::mat AH;

  // Dual Variables
  arma::mat U;
  arma::mat V;

  // Auxiliary/Temporary Variables
  arma::mat Htaux;
  arma::mat tempHtaux;
  arma::mat H0;
  arma::mat Wtaux;
  arma::mat tempWtaux;
  arma::mat W0;
  arma::mat L;

  // Hyperparameters
  double alpha, beta, tolerance;
  int admm_iter;

  /*
   * Collected statistics are
   * iteration Htime Wtime totaltime normH normW densityH densityW relError
   */
  void allocateMatrices() {
    WtW = arma::zeros<arma::mat>(this->k, this->k);
    HtH = arma::zeros<arma::mat>(this->k, this->k);
    WtA = arma::zeros<arma::mat>(this->n, this->k);
    AH = arma::zeros<arma::mat>(this->m, this->k);

    // Dual Variables
    U.zeros(size(this->W));
    V.zeros(size(this->H));

    // Auxiliary/Temporary Variables
    Htaux.zeros(size(this->H.t()));
    H0.zeros(size(this->H));
    tempHtaux.zeros(size(this->H.t()));
    Wtaux.zeros(size(this->W.t()));
    W0.zeros(size(this->W));
    tempWtaux.zeros(size(this->W.t()));
    L.zeros(this->k, this->k);

    // Hyperparameters
    alpha = 0.0;
    beta = 0.0;
    tolerance = 0.01;
    admm_iter = 5;
  }
  void freeMatrices() {
    this->At.clear();
    WtW.clear();
    HtH.clear();
    WtA.clear();
    AH.clear();
  }

 public:
  AOADMMNMF(const T &A, int lowrank, const unsigned int ncores = 0) : NMF<T>(A, lowrank) {
    this->normalize_by_W();
    allocateMatrices();
  }
  AOADMMNMF(const T &A, const arma::mat &llf, const arma::mat &rlf, const unsigned int ncores = 0) : NMF<T>(A, llf, rlf) {
    this->normalize_by_W();
    allocateMatrices();
  }
  void computeNMF() {
    unsigned int currentIteration = 0;
    this->At = this->A.t();
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
      beta = trace(WtW) / this->k;
      beta = beta > 0 ? beta : 0.01;
      WtW.diag() += beta;
      #ifdef _VERBOSE
      INFO << "starting H Prereq for "
           << " took=" << toc() << PRINTMATINFO(WtW) << PRINTMATINFO(WtA)
           << std::endl;
      #endif
      // to avoid divide by zero error.
      tic();
      L = arma::chol(WtW, "lower");

      bool stop_iter = false;

      // Start ADMM loop from here
      for (int i = 0; i < admm_iter && !stop_iter; i++) {
        H0 = this->H;
        tempHtaux =
            arma::solve(arma::trimatl(L), WtA + (beta * (this->H.t() + V.t())));
        Htaux = arma::solve(arma::trimatu(L.t()), tempHtaux);

        this->H = Htaux.t();
        // Uncomment if numerical issues are seen
        // fixNumericalError<arma::mat>(&(this->H), EPSILON_1EMINUS16, 0.0);
        this->H = this->H - V;
        this->H.for_each(
            [](arma::mat::elem_type &val) { val = val > 0.0 ? val : 0.0; });
        V = V + this->H - Htaux.t();

        // Check stopping criteria
        double r = norm(this->H - Htaux.t(), "fro");
        double s = norm(this->H - H0, "fro");
        double normH = norm(this->H, "fro");
        double normV = norm(V, "fro");

        if (r < (tolerance * normH) && s < (tolerance * normV))
          stop_iter = true;
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
      alpha = trace(HtH) / this->k;
      alpha = alpha > 0 ? alpha : 0.01;
      HtH.diag() += alpha;
      #ifdef _VERBOSE
      INFO << "starting W Prereq for "
           << " took=" << toc() << PRINTMATINFO(HtH) << PRINTMATINFO(AH)
           << std::endl;
      #endif
      tic();
      L = arma::chol(HtH, "lower");

      stop_iter = false;

      // Start ADMM loop from here
      for (int i = 0; i < admm_iter && !stop_iter; i++) {
        W0 = this->W;
        tempWtaux = arma::solve(arma::trimatl(L),
                                AH.t() + alpha * (this->W.t() + U.t()));
        Wtaux = arma::solve(arma::trimatu(L.t()), tempWtaux);

        this->W = Wtaux.t();
        // Uncomment if numerical issues are seen
        // fixNumericalError<arma::mat>(&(this->W), EPSILON_1EMINUS16, 0.0);
        this->W = this->W - U;
        this->W.for_each(
            [](arma::mat::elem_type &val) { val = val > 0.0 ? val : 0.0; });

        U = U + this->W - Wtaux.t();

        // Check stopping criteria
        double r = norm(this->W - Wtaux.t(), "fro");
        double s = norm(this->W - W0, "fro");
        double normW = norm(this->W, "fro");
        double normU = norm(U, "fro");

        if (r < (tolerance * normW) && s < (tolerance * normU))
          stop_iter = true;
      }
      #ifdef _VERBOSE
      INFO << "Completed W (" << currentIteration << "/"
           << this->num_iterations() << ")"
           << " time =" << toc() << std::endl;

      INFO << "Completed It (" << currentIteration << "/"
           << this->num_iterations() << ")"
           << " time =" << toc() << std::endl;
      #endif
      this->computeObjectiveError();
      #ifdef _VERBOSE
      this->printObjective(currentIteration);
      #endif
      currentIteration++;
    }
  }
  ~AOADMMNMF() = default;
};

}  // namespace planc
