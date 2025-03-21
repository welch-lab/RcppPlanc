#' Perform Integrative Non-negative Matrix Factorization
#' @description
#' Performs integrative non-negative matrix factorization (iNMF) (J.D. Welch,
#' 2019) to return factorized \eqn{H}, \eqn{W}, and \eqn{V} matrices. The
#' objective function is stated as
#'
#' \deqn{\arg\min_{H\ge0,W\ge0,V\ge0}\sum_{i}^{d}||E_i-(W+V_i)Hi||^2_F+
#' \lambda\sum_{i}^{d}||V_iH_i||_F^2}
#'
#' where \eqn{E_i} is the input non-negative matrix of the \eqn{i}'th dataset,
#' \eqn{d} is the total number of datasets. \eqn{E_i} is of size
#' \eqn{m \times n_i} for \eqn{m} features and \eqn{n_i} sample points,
#' \eqn{H_i} is of size \eqn{k \times n_i}, \eqn{V_i} is of size
#' \eqn{m \times k}, and \eqn{W} is of size \eqn{m \times k}.
#'
#' \code{inmf} optimizes the objective with ANLS strategy, while
#' \code{\link{onlineINMF}} optimizes the same objective with an online learning
#' strategy.
#' @param objectList list of input datasets. List elements should all be of the
#' same class. Viable classes include: matrix, dgCMatrix,
#' \link{H5Mat}, \link{H5SpMat}.
#' @param k Integer. Inner dimensionality to factorize the datasets into.
#' Default \code{20}.
#' @param lambda Regularization parameter. Larger values penalize
#' dataset-specific effects more strongly (i.e. alignment should increase as
#' \code{lambda} increases). Default \code{5}.
#' @param niter Integer. Total number of block coordinate descent iterations to
#' perform. Default \code{30}.
#' @param nCores The number of parallel tasks that will be spawned.
#' Default \code{2}
#' @param Hinit Initial values to use for \eqn{H} matrices. A list object where
#' each element is the initial \eqn{H} matrix of each dataset. Each should be
#' dense matrix of size \eqn{n_i \times k}. Default \code{NULL}.
#' @param Vinit Similar to \code{Hinit}, but each should be of size
#' \eqn{m \times k}.
#' @param Winit Initial values to use for \eqn{W} matrix. A matrix object of
#' size \eqn{m \times k}. Default \code{NULL}.
#' @param verbose Logical scalar. Whether to show information and progress.
#' Default \code{FALSE}.
#' @return A list of the following elements:
#' \itemize{
#'  \item{\code{H} - a list of result \eqn{H_i} matrices of size
#'  \eqn{n_i \times k}}
#'  \item{\code{V} - a list of result \eqn{V_i} matrices}
#'  \item{\code{W} - the result \eqn{W} matrix}
#'  \item{\code{objErr} - the final objective error value.}
#' }
#' @author Yichen Wang
#' @references Joshua D. Welch and et al., Single-Cell Multi-omic Integration
#' Compares and Contrasts Features of Brain Cell Identity, Cell, 2019
#' @examples
#' library(Matrix)
#' set.seed(1)
#' result <- inmf(list(ctrl.sparse, stim.sparse), k = 10, niter = 10, verbose = FALSE)
inmf <- function(
  objectList,
  k = 20,
  lambda = 5,
  niter = 30,
  nCores = 2,
  Hinit = NULL,
  Vinit = NULL,
  Winit = NULL,
  verbose = FALSE
) {
  mode <- .typeOfInput(objectList)
  res <- switch(
    mode,
    matrix = .bppinmf(objectList, k, nCores, lambda, niter, verbose,
                      Hinit, Vinit, Winit),
    dgCMatrix = .bppinmf(objectList, k, nCores, lambda, niter, verbose,
                         Hinit, Vinit, Winit),
    H5Mat = .bppinmf_h5dense(sapply(objectList, function(x) x$filename),
                             sapply(objectList, function(x) x$dataPath),
                             k, nCores, lambda, niter, verbose, Hinit, Vinit, Winit),
    H5SpMat = .bppinmf_h5sparse(filenames = sapply(objectList, function(x) x$filename),
                                valuePath = sapply(objectList, function(x) x$valuePath),
                                rowindPath = sapply(objectList, function(x) x$rowindPath),
                                colptrPath = sapply(objectList, function(x) x$colptrPath),
                                nrow = sapply(objectList, function(x) x$nrow),
                                ncol = sapply(objectList, function(x) x$ncol),
                                k = k, nCores = nCores, lambda = lambda, niter = niter,
                                verbose = verbose, Hinit = Hinit,
                                Vinit = Vinit, Winit = Winit)
  )
  names(res$H) <- names(res$V) <- names(objectList)
  return(res)
}

#' Perform Integrative Non-negative Matrix Factorization Using Online Learning
#' @description
#' Performs integrative non-negative matrix factorization (iNMF) (J.D. Welch,
#' 2019, C. Gao, 2021) using online learning approach to return factorized
#' \eqn{H}, \eqn{W}, and \eqn{V} matrices. The objective function is stated as
#'
#' \deqn{\arg\min_{H\ge0,W\ge0,V\ge0}\sum_{i}^{d}||E_i-(W+V_i)Hi||^2_F+
#' \lambda\sum_{i}^{d}||V_iH_i||_F^2}
#'
#' where \eqn{E_i} is the input non-negative matrix of the \eqn{i}'th dataset,
#' \eqn{d} is the total number of datasets. \eqn{E_i} is of size
#' \eqn{m \times n_i} for \eqn{m} features and \eqn{n_i} sample points,
#' \eqn{H_i} is of size \eqn{k \times n_i}, \eqn{V_i} is of size
#' \eqn{m \times k}, and \eqn{W} is of size \eqn{m \times k}.
#'
#' Different from \code{\link{inmf}} which optimizes the objective with ANLS
#' approach, \code{onlineINMF} optimizes the same objective with online learning
#' strategy, where it updates mini-batches of \eqn{H_i} solving the NNLS
#' problem, and updates \eqn{V_i} and \eqn{W} with HALS multiplicative method.
#'
#' This function allows online learning in 3 scenarios:
#'
#' \enumerate{
#'  \item Fully observed datasets;
#'  \item Iterative refinement using continually arriving datasets;
#'  \item Projection of new datasets without updating the existing factorization
#' }
#' @param objectList list of input datasets. List elements should all be of the
#' same class. Viable classes include: matrix, dgCMatrix,
#' \link{H5Mat}, \link{H5SpMat}.
#' @param newDatasets Same requirements as for new arriving datasets. Default
#' \code{NULL} for scenario 1, specify for scenario 2 or 3.
#' @param project Logical scalar, whether to run scenario 3. See description.
#' Default  \code{FALSE}.
#' @param k Integer. Inner dimensionality to factorize the datasets into.
#' Default \code{20}.
#' @param lambda Regularization parameter. Larger values penalize
#' dataset-specific effects more strongly (i.e. alignment should increase as
#' \code{lambda} increases). Default \code{5}.
#' @param maxEpoch The number of epochs to iterate through. Default \code{5}.
#' @param minibatchSize Total number of cells in each mini-batch. Default
#' \code{5000}.
#' @param maxHALSIter Maximum number of block coordinate descent (HALS
#' algorithm) iterations to perform for each update of \eqn{W} and \eqn{V}.
#' Default \code{1}. Changing this parameter is not recommended.
#' @param permuteChunkSize Number of cells in a chunk being shuffled before
#' subsetting to minibatches. Only appliable to in-memory data and for Scenario 
#' 1 and 2. Default \code{1000}.
#' @param nCores The number of parallel tasks that will be spawned.
#' Default \code{2}
#' @param Hinit,Vinit,Winit,Ainit,Binit Pass the previous factorization result
#' for datasets existing in \code{objectList}, in order to run scenario 2 or 3.
#' All should have \code{length(objectList)} matrices inside. See description
#' for dimensionality of \eqn{H_i}, \eqn{V_i} and \eqn{W_i}. \eqn{A_i} should
#' be of size \eqn{k \times k} and \eqn{B_i} should be of size \eqn{m \times k}
#' @param verbose Logical scalar. Whether to show information and progress.
#' Default \code{FALSE}.
#' @return A list of the following elements:
#' \itemize{
#'  \item{\code{H} - a list of result \eqn{H_i} matrices of size
#'  \eqn{n_i \times k}}
#'  \item{\code{V} - a list of result \eqn{V_i} matrices}
#'  \item{\code{W} - the result \eqn{W} matrix}
#'  \item{\code{A} - a list of result \eqn{A_i} matrices, \eqn{k \times k}}
#'  \item{\code{B} - a list of result \eqn{B_i} matrices, \eqn{m \times k}}
#'  \item{\code{objErr} - the final objective error value.}
#' }
#' @author Yichen Wang
#' @references Joshua D. Welch and et al., Single-Cell Multi-omic Integration
#' Compares and Contrasts Features of Brain Cell Identity, Cell, 2019
#'
#' Chao Gao and et al., Iterative single-cell multi-omic integration using
#' online learning, Nat Biotechnol., 2021
#' @examples
#' library(Matrix)
#'
#' # Scenario 1 with sparse matrices
#' set.seed(1)
#' res1 <- onlineINMF(list(ctrl.sparse, stim.sparse),
#'                    minibatchSize = 50, k = 10, verbose = FALSE)
#'
#' # Scenario 2 with H5 dense matrices
#' h5dense1 <- H5Mat(filename = system.file("extdata", "ctrl_dense.h5",
#'                              package = "RcppPlanc", mustWork = TRUE),
#'                                           dataPath = "scaleData")
#' h5dense2 <- H5Mat(filename = system.file("extdata", "stim_dense.h5",
#'                              package = "RcppPlanc", mustWork = TRUE),
#'                                           dataPath = "scaleData")
#' res2 <- onlineINMF(list(ctrl = h5dense1), minibatchSize = 50, k = 10, verbose = FALSE)
#' res3 <- onlineINMF(list(ctrl = h5dense1),
#'                    newDatasets = list(stim = h5dense2),
#'                    Hinit = res2$H, Vinit = res2$V, Winit = res2$W,
#'                    Ainit = res2$A, Binit = res2$B,
#'                    minibatchSize = 50, k = 10, verbose = FALSE)
#'
#' # Scenario 3 with H5 sparse matrices
#' h5sparse1 <- H5SpMat(filename = system.file("extdata", "ctrl_sparse.h5",
#'                                 package = "RcppPlanc", mustWork = TRUE),
#'                                 valuePath = "scaleDataSparse/data",
#'                                 rowindPath = "scaleDataSparse/indices",
#'                                 colptrPath = "scaleDataSparse/indptr",
#'                                 nrow = nrow(ctrl.sparse),
#'                                 ncol = ncol(ctrl.sparse))
#' h5sparse2 <- H5SpMat(filename = system.file("extdata", "stim_sparse.h5",
#'                                 package = "RcppPlanc", mustWork = TRUE),
#'                                 valuePath = "scaleDataSparse/data",
#'                                 rowindPath = "scaleDataSparse/indices",
#'                                 colptrPath = "scaleDataSparse/indptr",
#'                                 nrow = nrow(stim.sparse),
#'                                 ncol = nrow(stim.sparse))
#' res4 <- onlineINMF(list(ctrl = h5sparse1), minibatchSize = 50, k = 10, verbose = FALSE)
#' res5 <- onlineINMF(list(ctrl = h5sparse1),
#'                    newDatasets = list(stim = h5sparse2), project = TRUE,
#'                    Hinit = res4$H, Vinit = res4$V, Winit = res4$W,
#'                    Ainit = res4$A, Binit = res4$B,
#'                    minibatchSize = 50, k = 10, verbose = FALSE)
#'
onlineINMF <- function(
  objectList,
  newDatasets = NULL,
  project = FALSE,
  k = 20,
  lambda = 5,
  maxEpoch = 5,
  minibatchSize = 5000,
  maxHALSIter = 1,
  permuteChunkSize = 1000,
  nCores = 2,
  Hinit = NULL,
  Vinit = NULL,
  Winit = NULL,
  Ainit = NULL,
  Binit = NULL,
  verbose = FALSE
) {
  mode <- .typeOfInput(objectList)
  if (is.null(newDatasets)) {
    # Scenario 1
    res <- switch(
      mode,
      matrix = .onlineINMF(objectList, k, nCores, lambda, maxEpoch,
                           minibatchSize, maxHALSIter, permuteChunkSize, 
                           verbose),
      dgCMatrix = .onlineINMF(objectList, k, nCores, lambda, maxEpoch,
                              minibatchSize, maxHALSIter, permuteChunkSize,
                              verbose),
      H5Mat = .onlineINMF_h5dense(
        sapply(objectList, function(x) x$filename),
        sapply(objectList, function(x) x$dataPath),
        k, nCores, lambda, maxEpoch, minibatchSize,
        maxHALSIter, permuteChunkSize, verbose
      ),
      H5SpMat = .onlineINMF_h5sparse(
        sapply(objectList, function(x) x$filename),
        sapply(objectList, function(x) x$valuePath),
        sapply(objectList, function(x) x$rowindPath),
        sapply(objectList, function(x) x$colptrPath),
        sapply(objectList, function(x) x$nrow),
        sapply(objectList, function(x) x$ncol),
        k, nCores, lambda, maxEpoch, minibatchSize,
        maxHALSIter, permuteChunkSize, verbose
      )
    )
    names(res$H) <- names(res$V) <- names(res$A) <- names(res$B) <- names(objectList)
  } else {
    mode2 <- .typeOfInput(newDatasets)
    if (mode2 != mode) {
      stop("newDatasets should be of the same class as original datasets")
    }
    if (!isTRUE(project)) {
      # Scenario 2 result
      res <- switch(
        mode,
        matrix = .onlineINMF_withInitial(
          objectList, Hinit, Vinit, Winit, Ainit, Binit,
          newDatasets, k, nCores, lambda, maxEpoch,
          minibatchSize, maxHALSIter, permuteChunkSize, verbose
        ),
        dgCMatrix = .onlineINMF_withInitial(
          objectList, Hinit, Vinit, Winit, Ainit, Binit,
          newDatasets, k, nCores, lambda,
          maxEpoch, minibatchSize, maxHALSIter, permuteChunkSize, verbose
        ),
        H5Mat = .onlineINMF_h5dense_withInitial(
          sapply(objectList, function(x) x$filename),
          sapply(objectList, function(x) x$dataPath),
          sapply(newDatasets, function(x) x$filename),
          sapply(newDatasets, function(x) x$dataPath),
          Hinit, Vinit, Winit, Ainit, Binit, k, nCores, lambda, maxEpoch,
          minibatchSize, maxHALSIter, permuteChunkSize, verbose
        ),
        H5SpMat = .onlineINMF_h5sparse_withInitial(
          sapply(objectList, function(x) x$filename),
          sapply(objectList, function(x) x$valuePath),
          sapply(objectList, function(x) x$rowindPath),
          sapply(objectList, function(x) x$colptrPath),
          sapply(objectList, function(x) x$nrow),
          sapply(objectList, function(x) x$ncol),
          sapply(newDatasets, function(x) x$filename),
          sapply(newDatasets, function(x) x$valuePath),
          sapply(newDatasets, function(x) x$rowindPath),
          sapply(newDatasets, function(x) x$colptrPath),
          sapply(newDatasets, function(x) x$nrow),
          sapply(newDatasets, function(x) x$ncol),
          Hinit, Vinit, Winit, Ainit, Binit,
          k, nCores, lambda, maxEpoch, minibatchSize, maxHALSIter, 
          permuteChunkSize, verbose
        )
      )
      names(res$H) <- names(res$V) <- names(res$A) <- names(res$B) <-
        c(names(objectList), names(newDatasets))
    } else {
      res <- switch(
        mode,
        matrix = .onlineINMF_project(objectList, Winit, newDatasets, k, nCores, lambda),
        dgCMatrix = .onlineINMF_project(objectList, Winit, newDatasets, k, nCores, lambda),
        H5Mat = .onlineINMF_project_h5dense(
          sapply(objectList, function(x) x$filename),
          sapply(objectList, function(x) x$dataPath),
          sapply(newDatasets, function(x) x$filename),
          sapply(newDatasets, function(x) x$dataPath),
          Winit, k, nCores, lambda
        ),
        H5SpMat = .onlineINMF_project_h5sparse(
          sapply(objectList, function(x) x$filename),
          sapply(objectList, function(x) x$valuePath),
          sapply(objectList, function(x) x$rowindPath),
          sapply(objectList, function(x) x$colptrPath),
          sapply(objectList, function(x) x$nrow),
          sapply(objectList, function(x) x$ncol),
          sapply(newDatasets, function(x) x$filename),
          sapply(newDatasets, function(x) x$valuePath),
          sapply(newDatasets, function(x) x$rowindPath),
          sapply(newDatasets, function(x) x$colptrPath),
          sapply(newDatasets, function(x) x$nrow),
          sapply(newDatasets, function(x) x$ncol),
          Winit, k, nCores, lambda
        ),
        # Scenario 3 result
        names(res$H) <- names(newDatasets))
    }
  }
  return(res)
}

#' Perform Mosaic Integrative Non-negative Matrix Factorization with Unshared
#' Features
#' @description
#' Performs mosaic integrative non-negative matrix factorization (UINMF) (A.R.
#' Kriebel, 2022) to return factorized \eqn{H}, \eqn{W}, \eqn{V} and \eqn{U}
#' matrices. The objective function is stated as
#'
#' \deqn{\arg\min_{H\ge0,W\ge0,V\ge0,U\ge0}\sum_{i}^{d}
#' ||\begin{bmatrix}E_i \\ P_i \end{bmatrix} -
#' (\begin{bmatrix}W \\ 0 \end{bmatrix}+
#' \begin{bmatrix}V_i \\ U_i \end{bmatrix})Hi||^2_F+
#' \lambda_i\sum_{i}^{d}||\begin{bmatrix}V_i \\ U_i \end{bmatrix}H_i||_F^2}
#'
#' where \eqn{E_i} is the input non-negative matrix of the \eqn{i}'th dataset,
#' \eqn{P_i} is the input non-negative matrix for the unshared features,
#' \eqn{d} is the total number of datasets. \eqn{E_i} is of size
#' \eqn{m \times n_i} for \eqn{m} shared features and \eqn{n_i} sample points,
#' \eqn{P_i} is of size \eqn{u_i \times n_i} for \eqn{u_i} unshared feaetures,
#' \eqn{H_i} is of size \eqn{k \times n_i}, \eqn{V_i} is of size
#' \eqn{m \times k}, \eqn{W} is of size \eqn{m \times k} and \eqn{U_i} is of
#' size \eqn{u_i \times k}.
#'
#' Similar to \code{\link{inmf}}, \code{uinmf} also optimizes the objective with
#' ANLS algorithm.
#' @param objectList list of input datasets. List elements should all be of the
#' same class. Viable classes include: matrix, dgCMatrix,
#' \link{H5Mat}, \link{H5SpMat}.
#' @param unsharedList List of input unshared feature matrices, with the same
#' requirement as \code{objectList}.
#' @param k Integer. Inner dimensionality to factorize the datasets into.
#' Default \code{20}.
#' @param lambda Regularization parameter. Use one number for all datasets or a
#' vector to specify for each dataset. Larger values penalize dataset-specific
#' effects more strongly (i.e. alignment should increase as \code{lambda}
#' increases). Default \code{5}.
#' @param niter Integer. Total number of block coordinate descent iterations to
#' perform. Default \code{30}.
#' @param nCores The number of parallel tasks that will be spawned.
#' Default \code{2}.
#' @param verbose Logical scalar. Whether to show information and progress.
#' Default \code{FALSE}.
#' @return A list of the following elements:
#' \itemize{
#'  \item{\code{H} - a list of result \eqn{H_i} matrices of size
#'  \eqn{n_i \times k}}
#'  \item{\code{V} - a list of result \eqn{V_i} matrices}
#'  \item{\code{W} - the result \eqn{W} matrix}
#'  \item{\code{U} - a list of result \eqn{A_i} matrices}
#'  \item{\code{objErr} - the final objective error value.}
#' }
#' @author Yichen Wang
#' @references April R. Kriebel and Joshua D. Welch, UINMF performs mosaic
#' integration of single-cell multi-omic datasets using nonnegative matrix
#' factorization, Nat. Comm., 2022
#' @examples
#' # Fake matrices representing unshared features of the given datasets
#' # Real-life use should have features that are not presented in the
#' # intersection of features of all datasets involved.
#' ctrl.unshared <- ctrl.sparse[1:10,]
#' stim.unshared <- stim.sparse[11:30,]
#' set.seed(1)
#' result <- uinmf(list(ctrl.sparse, stim.sparse),
#'                 list(ctrl.unshared, stim.unshared), verbose = FALSE)
uinmf <- function(
  objectList,
  unsharedList,
  k = 20,
  lambda = 5,
  niter = 30,
  nCores = 2,
  verbose = FALSE
) {
  if (length(lambda) == 1) lambda <- rep(lambda, length(objectList))
  if (length(lambda) != length(objectList)) {
    stop("Must specify 1 lambda for all or each.")
  }
  mode <- .typeOfInput(objectList, null.rm = FALSE)
  unshareCleanUp <- .uinmf.matchDatasets(objectList, unsharedList)
  unsharedList <- unshareCleanUp[[1]]
  whichUnshared <- unshareCleanUp[[2]]
  res <- switch(
    mode,
    matrix = .uinmf_rcpp(objectList, unsharedList, whichUnshared,
                         k, nCores, lambda, niter, verbose),
    dgCMatrix = .uinmf_rcpp(objectList, unsharedList, whichUnshared,
                            k, nCores, lambda, niter, verbose),
    H5Mat = .uinmf_h5dense(sapply(objectList, function(x) x$filename),
                           sapply(objectList, function(x) x$dataPath),
                           sapply(unsharedList, function(x) x$filename),
                           sapply(unsharedList, function(x) x$dataPath),
                           whichUnshared, k, nCores, lambda, niter, verbose),
    H5SpMat = .uinmf_h5sparse(sapply(objectList, function(x) x$filename),
                              sapply(objectList, function(x) x$rowindPath),
                              sapply(objectList, function(x) x$colptrPath),
                              sapply(objectList, function(x) x$valuePath),
                              sapply(objectList, function(x) x$nrow),
                              sapply(objectList, function(x) x$ncol),
                              sapply(unsharedList, function(x) x$filename),
                              sapply(unsharedList, function(x) x$rowindPath),
                              sapply(unsharedList, function(x) x$colptrPath),
                              sapply(unsharedList, function(x) x$valuePath),
                              sapply(unsharedList, function(x) x$nrow),
                              sapply(unsharedList, function(x) x$ncol),
                              whichUnshared, k, nCores, lambda, niter, verbose)
  )
  names(res$H) <- names(res$V) <- names(objectList)
  names(res$U) <- names(unsharedList)
  return(res)
}

# Check and returns a processed `unsharedList` so it can be used with c++
# also returns a matching vector so we don't make empty dataset and have error
# If `objectList` has name a, b, c, d, and unsharedList has b, d, c,
# then matching vector has c(-1, 0, 2, 1). `-1` indicates no shared features,
# non-negs are 0-based index for which mat in `unsharedList` is the one for this
# object.
.uinmf.matchDatasets <- function(objectList, unsharedList) {
  if (is.null(names(objectList)) || is.null(names(unsharedList))) {
    if (length(unsharedList) != length(objectList)) {
      stop("Number of matrix in unshared feature list does not match, ",
           "please use named lists to indicate dataset matching.")
    }
    # No naming matching available, assume one-by-one matching
    names(objectList) <- names(unsharedList) <-
      as.character(seq_along(objectList))
  }
  # Remove unshared ones not existing in objectList First
  unsharedList <- unsharedList[names(unsharedList) %in% names(objectList)]
  # Remove empty unshared ones
  unsharedList <- unsharedList[sapply(unsharedList, function(x) {
    if (is.null(x)) return(FALSE)
    if (is.null(nrow(x))) return(FALSE)
    if (nrow(x) == 0) return(FALSE)
    return(TRUE)
  })]

  whichUnshared <- rep(-1, length(objectList))
  for (i in seq_along(objectList)) {
    d <- names(objectList)[i]
    if (d %in% names(unsharedList)) {
      # Guaranteed to cover all unshared data, b/c cleaned up above
      if (ncol(unsharedList[[d]]) != ncol(objectList[[d]])) {
        stop("Number of columns in each matrix from `unsharedList` ",
             "must match with the corresponding matrix from `object`")
      }
      whichUnshared[i] <- which(names(unsharedList) == d)[1] - 1
    }
  }

  mode <- .typeOfInput(objectList, null.rm = FALSE)
  mode2 <- .typeOfInput(unsharedList, null.rm = FALSE)
  if (mode != mode2) {
    stop("Data of unshared feature should be of the same class as ",
         "input data")
  }
  return(list(unsharedList, whichUnshared))
}
