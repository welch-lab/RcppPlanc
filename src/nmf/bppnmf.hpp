#pragma once
/* Copyright 2016 Ramakrishnan Kannan */

#ifdef _OPENMP
#include <omp.h>
#endif
#include "nmf.hpp"
#include "bppnnls.hpp"
#include "hals.hpp"
#define ONE_THREAD_MATRIX_SIZE 2000

namespace planc {

template <typename T>
class BPPNMF : public NMF<T> {
 private:
  T At;
  arma::mat giventGiven;
  // designed as if W is given and H is found.
  // The transpose is the other problem.
  void updateOtherGivenOneMultipleRHS(const T &input, const arma::mat &given,
                                      char worh, arma::mat *othermat, arma::fvec reg) {
#if defined(_VERBOSE) || defined(COLLECTSTATS)
    double t2;
#endif
    unsigned int numChunks = input.n_cols / ONE_THREAD_MATRIX_SIZE;
    if (numChunks * ONE_THREAD_MATRIX_SIZE < input.n_cols) numChunks++;
#if defined(_VERBOSE) || defined(COLLECTSTATS)
    tic();
#endif
    arma::mat giventInput(this->k, input.n_cols);
    // This is WtW
    giventGiven = given.t() * given;
    this->applyReg(reg, &giventGiven);
    // This is WtA
    giventInput = given.t() * input;
    if (this->symm_reg() > 0) {
      arma::mat fac = given.t();
      this->applySymmetricReg(this->symm_reg(), &giventGiven, &fac,
                              &giventInput);
    }
#if defined(_VERBOSE) || defined(COLLECTSTATS)
    t2 = toc();
#endif
#ifdef _VERBOSE
    INFO << "starting " << worh << ". Prereq for " << worh << " took=" << t2
         << " NumChunks=" << numChunks << std::endl;
    INFO << "LHS::" << std::endl
         << giventGiven << std::endl
         << "RHS::" << std::endl
         << giventInput << std::endl;
#endif
#if defined(_VERBOSE) || defined(COLLECTSTATS)
    tic();
#endif
#pragma omp parallel for schedule(auto)
    for (unsigned int i = 0; i < numChunks; i++) {
      unsigned int spanStart = i * ONE_THREAD_MATRIX_SIZE;
      unsigned int spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
      if (spanEnd > input.n_cols - 1) {
        spanEnd = input.n_cols - 1;
      }

      BPPNNLS<arma::mat, arma::vec> subProblem(giventGiven,
                              (arma::mat)giventInput.cols(spanStart, spanEnd), true);
#ifdef _VERBOSE
      #pragma omp critical
      {
        INFO << "Scheduling " << worh << " start=" << spanStart
             << ", end=" << spanEnd
             // << ", tid=" << omp_get_thread_num()
             << std::endl
             << "LHS ::" << std::endl
             << giventGiven << std::endl
             << "RHS ::" << std::endl
             << (arma::mat)giventInput.cols(spanStart, spanEnd) << std::endl;
      }
#endif

      subProblem.solveNNLS();

#ifdef _VERBOSE
      INFO << "completed " << worh << " start=" << spanStart
           << ", end=" << spanEnd
           // << ", tid=" << omp_get_thread_num() << " cpu=" << sched_getcpu()
           << " time taken=" << t2 << std::endl;
#endif
      (*othermat).rows(spanStart, spanEnd) = subProblem.getSolutionMatrix().t();
    }
#if defined(_VERBOSE) || defined(COLLECTSTATS)
    double totalH2 = toc();
#endif
#ifdef _VERBOSE
    INFO << worh << " total time taken :" << totalH2 << std::endl;
#endif
    giventGiven.clear();
    giventInput.clear();
  }
void commonSolve() {
    unsigned int currentIteration = 0;
    while (currentIteration < this->num_iterations()) {
#ifdef COLLECTSTATS
      this->collectStats(currentIteration);
      this->stats(currentIteration + 1, 0) = currentIteration + 1;
#endif
#if defined(_VERBOSE) || defined(COLLECTSTATS)
      tic();
#endif
      updateOtherGivenOneMultipleRHS(this->At, this->H, 'W', &(this->W),
                                     this->regW());
#if defined(_VERBOSE) || defined(COLLECTSTATS)
      double totalW2 = toc();
      tic();
#endif
      updateOtherGivenOneMultipleRHS(this->A, this->W, 'H', &(this->H),
                                     this->regH());
#if defined(_VERBOSE) || defined(COLLECTSTATS)
      double totalH2 = toc();
#endif

#ifdef COLLECTSTATS
      // end of H and start of W are almost same.
      this->stats(currentIteration + 1, 1) = totalH2;
      this->stats(currentIteration + 1, 2) = totalW2;

      this->stats(currentIteration + 1, 3) = totalW2 + totalH2;
#endif
#ifdef _VERBOSE
      INFO << "Completed It (" << currentIteration << "/"
           << this->num_iterations() << ")"
           << " time =" << totalW2 + totalH2 << std::endl;
#endif
      this->computeObjectiveError();
#ifdef _VERBOSE
      this->printObjective(currentIteration);

#endif
      currentIteration++;
    }
    this->normalize_by_W();
#ifdef COLLECTSTATS
    this->collectStats(currentIteration);
    INFO << "NMF Statistics:" << std::endl << this->stats << std::endl;
#endif
  };
 public:
  BPPNMF(const T &A, int lowrank) : NMF<T>(A, lowrank) {
    giventGiven = arma::zeros<arma::mat>(lowrank, lowrank);
    this->At = A.t();
  }
  BPPNMF(const T &A, const arma::mat &llf, const arma::mat &rlf) : NMF<T>(A, llf, rlf) {
    this->At = A.t();
  }
  void computeNMFSingleRHS() {
    int currentIteration = 0;
    T At = this->A.t();
    this->computeObjectiveErr();
    while (currentIteration < this->num_iterations() &&
           this->objectiveErr > CONV_ERR) {
#ifdef COLLECTSTATS
      this->collectStats(currentIteration);
#endif
      // solve for H given W;
      arma::mat Wt = this->W.t();
      arma::mat WtW = Wt * this->W;
      this->applyReg(this->regH(), &this->WtW);
      arma::mat WtA = Wt * this->A;
      Wt.clear();
      {
#pragma omp parallel for schedule(auto)
        for (unsigned int i = 0; i < this->n; i++) {
          BPPNNLS<arma::mat, arma::vec> *subProblemforH =
              new BPPNNLS<arma::mat, arma::vec>(WtW, (arma::vec)WtA.col(i), true);
#ifdef _VERBOSE
          INFO << "Initialized subproblem and calling solveNNLS for "
               << "H(" << i << "/" << this->n << ")";
#endif
#if defined(_VERBOSE) || defined(COLLECTSTATS)
          tic();
          int numIter = subProblemforH->solveNNLS();
#else
          subProblemforH->solveNNLS();
#endif

#if defined(_VERBOSE) || defined(COLLECTSTATS)
          double t2 = toc();
#endif
#ifdef _VERBOSE
          INFO << subProblemforH->getSolutionVector();
#endif
          this->H.row(i) = subProblemforH->getSolutionVector().t();
#ifdef _VERBOSE
          INFO << "Comp H(" << i << "/" << this->n
               << ") of it=" << currentIteration << " time taken=" << t2
               << " num_iterations()=" << numIter << std::endl;
#endif
        }
      }
#ifdef _VERBOSE
      INFO << "H: at it = " << currentIteration << std::endl << this->H;
#endif
      {
        // clear previous allocations.
        WtW.clear();
        WtA.clear();
        arma::mat Ht = this->H.t();
        arma::mat HtH = Ht * this->H;
        this->applyReg(this->regW(), &this->HtH);
        arma::mat HtAt = Ht * At;
        Ht.clear();
// solve for W given H;
#pragma omp parallel for schedule(auto)
        for (unsigned int i = 0; i < this->m; i++) {
          BPPNNLS<arma::mat, arma::vec> *subProblemforW =
              new BPPNNLS<arma::mat, arma::vec>(HtH, (arma::vec)HtAt.col(i), true);
#ifdef _VERBOSE
          INFO << "Initialized subproblem and calling solveNNLS for "
               << "W(" << i << "/" << this->m << ")";
#endif
#if defined(_VERBOSE) || defined(COLLECTSTATS)
          tic();
          int numIter = subProblemforW->solveNNLS();
#else
          subProblemforW->solveNNLS();
#endif
      //     int numIter = subProblemforW->solveNNLS();
#if defined(_VERBOSE) || defined(COLLECTSTATS)
          double t2 = toc();
#endif
#ifdef _VERBOSE
          INFO << subProblemforW->getSolutionVector();
#endif

          this->W.row(i) = subProblemforW->getSolutionVector().t();
#ifdef _VERBOSE
          INFO << "Comp W(" << i << "/" << this->n
               << ") of it=" << currentIteration << " time taken=" << t2
               << " num_iterations()=" << numIter << std::endl;
#endif
        }
        HtH.clear();
        HtAt.clear();
      }
#ifdef _VERBOSE
      INFO << "W: at it = " << currentIteration << std::endl << this->W;
#endif
#ifdef COLLECTSTATS
      // INFO << "iteration = " << currentIteration << " currentObjectiveError="
      // << this->objective_err << std::endl;
#endif
      currentIteration++;
    }
  }


  void computeNMF() {
#ifdef COLLECTSTATS
    // this->objective_err;
#endif
#ifdef _VERBOSE
    INFO << PRINTMATINFO(this->At);
    INFO << "Starting BPP for num_iterations()=" << this->num_iterations()
         << std::endl;
#endif
this->commonSolve();
  };

  double getObjectiveError() { return this->objectiveErr; }

  /*
   * I dont like this function here. But this seems to be the
   * easy place for having it. This function really should have been
   * in BPPNNLS.hpp. It will take some time to refactor this.
   * Given, A and W, solve for H.
   */
  arma::mat solveScalableNNLS() {
    updateOtherGivenOneMultipleRHS(this->A, this->W, 'H', &(this->H));
    return this->H;
  }
  ~BPPNMF() { this->At.clear(); }
};  // class BPPNMF

template<>
void BPPNMF<arma::sp_mat>::computeNMF() {
//   unsigned int currentIteration = 0;
#ifdef COLLECTSTATS
  // this->objective_err;
#endif
  // tic();
  // this->At = this->A.t();  // do it once
  // INFO << "At time::" << toc() << std::endl;
  // run hals once to get proper initializations
  HALSNMF<arma::sp_mat> tempHals(this->A, this->W, this->H);
  tempHals.num_iterations(2);
  this->W = tempHals.getLeftLowRankFactor();
  this->H = tempHals.getRightLowRankFactor();
#ifdef _VERBOSE
  INFO << PRINTMATINFO(this->At);
  INFO << " nnz = " << this->At.n_nonzero << std::endl;
  INFO << "Starting BPP for num_iterations()=" << this->num_iterations()
       << std::endl;
#endif
  this->commonSolve();
}

}  // namespace planc
