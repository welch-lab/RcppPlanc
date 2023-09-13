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
#' \eqn{H_i} is of size \eqn{n_i \times k}, \eqn{V_i} is of size 
#' \eqn{m \times k}, and \eqn{W} is of size \eqn{m \times k}.
#' 
#' \code{inmf} optimizes the objective with ANLS strategy, while 
#' \code{\link{onlineINMF}} optimizes the same objective with an online learning
#' strategy. 
#' @param objectList list of input datasets. List elements should all be of the
#' same class. Viable classes include: matrix, \linkS4class{dgCMatrix}, 
#' \link{H5Mat}, \link{H5SpMat}. 
#' @param k Integer. Inner dimensionality to factorize the datasets into. 
#' Default \code{20}.
#' @param lambda Regularization parameter. Larger values penalize
#' dataset-specific effects more strongly (i.e. alignment should increase as
#' \code{lambda} increases). Default \code{5}.
#' @param niter Integer. Total number of block coordinate descent iterations to 
#' perform. Default \code{30}.
#' @param Hinit Initial values to use for \eqn{H} matrices. A list object where 
#' each element is the initial \eqn{H} matrix of each dataset. Each should be
#' dense matrix of size \eqn{n_i \times k}. Default \code{NULL}.
#' @param Vinit Similar to \code{Hinit}, but each should be of size 
#' \eqn{m \times k}.
#' @param Winit Initial values to use for \eqn{W} matrix. A matrix object of 
#' size \eqn{m \times k}. Default \code{NULL}.
#' @param verbose Logical scalar. Whether to show information and progress. 
#' Default \code{TRUE}.
#' @return A list of the following entries: \code{H} - a list of result 
#' \eqn{H_i} matrices, \code{V} - a list of result \eqn{V_i} matrices, 
#' \code{W} - the result \eqn{W} matrix, \code{objErr} - the final objective
#' error value.
#' @author Yichen Wang
#' @references Joshua D. Welch and et al., Single-Cell Multi-omic Integration 
#' Compares and Contrasts Features of Brain Cell Identity, Cell, 2019
#' @export
#' @examples
#' if (FALSE) {
#'     set.seed(1)
#'     res1 <- inmf(list(dense1, dense2))
#'     set.seed(1)
#'     res2 <- inmf(list(sparse1, sparse2))
#'     h5dense1 <- H5Mat("filepath1.h5", "dataPath")
#'     h5dense2 <- H5Mat("filepath2.h5", "dataPath")
#'     set.seed(1)
#'     res3 <- inmf(list(h5dense1, h5dense2))
#'     h5sparse1 <- H5SpMat("filepath1.h5", "group/data", "group/indices", 
#'                          "group/indptr", nrow = 4215, ncol = 6548)
#'     h5sparse2 <- H5SpMat("filepath2.h5", "group/data", "group/indices", 
#'                          "group/indptr", nrow = 4215, ncol = 7451)
#'     set.seed(1)
#'     res4 <- inmf(list(h5sparse1, h5sparse2))
#' }
inmf <- function(
    objectList,
    k = 20,
    lambda = 5,
    niter = 30,
    Hinit = NULL,
    Vinit = NULL,
    Winit = NULL,
    verbose = TRUE
) {
    mode <- .typeOfInput(objectList)
    res <- switch(
        mode,
        matrix = bppinmf(objectList, k, lambda, niter, verbose, 
                         Hinit, Vinit, Winit),
        dgCMatrix = bppinmf(objectList, k, lambda, niter, verbose, 
                            Hinit, Vinit, Winit),
        H5Mat = bppinmf_h5dense(sapply(objectList, function(x) x$filename),
                                sapply(objectList, function(x) x$dataPath),
                                k, lambda, niter, verbose, Hinit, Vinit, Winit),
        H5SpMat = bppinmf_h5sparse(sapply(objectList, function(x) x$filename),
                                   sapply(objectList, function(x) x$valuePath),
                                   sapply(objectList, function(x) x$rowindPath),
                                   sapply(objectList, function(x) x$colptrPath),
                                   sapply(objectList, function(x) x$nrow),
                                   sapply(objectList, function(x) x$ncol),
                                   k, lambda, niter, verbose, 
                                   Hinit, Vinit, Winit)
    )
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
#' \eqn{H_i} is of size \eqn{n_i \times k}, \eqn{V_i} is of size 
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
#' 
#' @param objectList list of input datasets. List elements should all be of the
#' same class. Viable classes include: matrix, \linkS4class{dgCMatrix}, 
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
#' @param maxEpochs The number of epochs to iterate through. Default \code{5}.
#' @param maxHALSIter Maximum number of block coordinate descent (HALS
#' algorithm) iterations to perform for each update of \eqn{W} and \eqn{V}.
#' Default \code{1}. Changing this parameter is not recommended.
#' @param miniBatchSize Total number of cells in each mini-batch. Default 
#' \code{5000}.
#' @param Vinit,Hinit,Ainit,Binit Pass the previous factorization result for
#' datasets existing in \code{objectList}, in order to run scenario 2 or 3. All 
#' should have \code{length(objectList)} matrices inside. See description for
#' dimensionality of \eqn{V_i} and \eqn{H_i}. \eqn{A_i} should be of size 
#' \eqn{k \times k} and \eqn{B_i} should be of size \eqn{m \times k}
#' @param verbose Logical scalar. Whether to show information and progress. 
#' Default \code{TRUE}.
#' @return A list of the following entries: \code{H} - a list of result 
#' \eqn{H_i} matrices, \code{V} - a list of result \eqn{V_i} matrices, 
#' \code{W} - the result \eqn{W} matrix, \code{A} - a list of result \eqn{A_i}
#' matrices, \code{B} - a list of result \eqn{B_i} matrices, \code{objErr} - the
#' final objective error value.
#' @author Yichen Wang
#' @references Joshua D. Welch and et al., Single-Cell Multi-omic Integration 
#' Compares and Contrasts Features of Brain Cell Identity, Cell, 2019
#' 
#' Chao Gao and et al., Iterative single-cell multi-omic integration using 
#' online learning, Nat Biotechnol., 2021
#' @export
#' @examples
#' if (FALSE) {
#'     # Scenario 1 with sparse matrices
#'     res1 <- onlineINMF(list(sparse1, sparse2))
#'     
#'     # Scenario 2 with H5 dense matrices
#'     h5dense1 <- H5Mat("filepath1.h5", "dataPath")
#'     h5dense2 <- H5Mat("filepath2.h5", "dataPath")
#'     res2 <- onlineINMF(list(h5dense1, h5dense2))
#'     h5dense3 <- H5Mat("filepath3.h5", "dataPath")
#'     res3 <- onlineINMF(list(h5dense1, h5dense2), 
#'                        newDatasets = list(h5dense3),
#'                        Vinit = res2$V, Winit = res2$W,
#'                        Ainit = res2$A, Binit = res2$B)
#'    
#'     # Scenario 3 with H5 sparse matrices
#'     h5sparse1 <- H5SpMat("filepath1.h5", "group/data", "group/indices", 
#'                          "group/indptr", nrow = 4215, ncol = 6548)
#'     h5sparse2 <- H5SpMat("filepath2.h5", "group/data", "group/indices", 
#'                          "group/indptr", nrow = 4215, ncol = 7451)
#'     res4 <- onlineINMF(list(h5sparse1, h5sparse2))
#'     h5sparse3 <- H5SpMat("filepath3.h5", "group/data", "group/indices", 
#'                          "group/indptr", nrow = 4215, ncol = 3456)
#'     res5 <- onlineINMF(list(h5sparse1, h5sparse2),
#'                        newDatasets = list(h5sparse3), project = TRUE,
#'                        Vinit = res4$V, Winit = res4$W,
#'                        Ainit = res4$A, Binit = res4$B)
#' }
onlineINMF <- function(
    objectList, 
    newDatasets = NULL,
    project = FALSE,
    k = 20, 
    lambda = 5, 
    maxEpoch = 5, 
    minibatchSize = 5000, 
    maxHALSIter = 1,
    Vinit = NULL,
    Winit = NULL,
    Ainit = NULL,
    Binit = NULL,
    verbose = TRUE
) {
    mode <- .typeOfInput(objectList)
    if (is.null(newDatasets)) {
        # Scenario 1
        res <- switch(
            mode,
            matrix = onlineINMF_S1(objectList, k, lambda, maxEpoch, 
                                   minibatchSize, maxHALSIter, verbose),
            dgCMatrix = onlineINMF_S1(objectList, k, lambda, maxEpoch, 
                                      minibatchSize, maxHALSIter, verbose),
            H5Mat = onlineINMF_S1_h5dense(sapply(objectList, 
                                                 function(x) x$filename),
                                          sapply(objectList, 
                                                 function(x) x$dataPath),
                                          k, lambda, maxEpoch, minibatchSize, 
                                          maxHALSIter, verbose),
            H5SpMat = onlineINMF_S1_h5sparse(sapply(objectList, 
                                                    function(x) x$filename),
                                             sapply(objectList, 
                                                    function(x) x$valuePath),
                                             sapply(objectList, 
                                                    function(x) x$rowindPath),
                                             sapply(objectList, 
                                                    function(x) x$colptrPath),
                                             sapply(objectList, 
                                                    function(x) x$nrow),
                                             sapply(objectList, 
                                                    function(x) x$ncol), 
                                             k, lambda, maxEpoch, minibatchSize,
                                             maxHALSIter, verbose)
        )
    } else {
        mode2 <- .typeOfInput(newDatasets)
        if (mode2 != mode) {
            stop("newDatasets should be of the same class as original datasets")
        }
        res <- switch(
            mode,
            matrix = onlineINMF_S23(objectList, Vinit, Winit, Ainit, Binit,
                                    newDatasets, k, lambda, project, maxEpoch,
                                    minibatchSize, maxHALSIter, verbose),
            dgCMatrix = onlineINMF_S23(objectList, Vinit, Winit, Ainit, Binit,
                                       newDatasets, k, lambda, project, 
                                       maxEpoch, minibatchSize, maxHALSIter, 
                                       verbose),
            H5Mat = onlineINMF_S23_h5dense(sapply(objectList, 
                                                  function(x) x$filename),
                                           sapply(objectList, 
                                                  function(x) x$dataPath),
                                           sapply(newDatasets, 
                                                  function(x) x$filename),
                                           sapply(newDatasets, 
                                                  function(x) x$dataPath),
                                           Vinit, Winit, Ainit, Binit,
                                           k, lambda, project, maxEpoch,
                                           minibatchSize, maxHALSIter, verbose),
            H5SpMat = onlineINMF_S23_h5sparse(sapply(objectList, 
                                                     function(x) x$filename),
                                              sapply(objectList, 
                                                     function(x) x$valuePath),
                                              sapply(objectList, 
                                                     function(x) x$rowindPath),
                                              sapply(objectList, 
                                                     function(x) x$colptrPath),
                                              sapply(objectList, 
                                                     function(x) x$nrow),
                                              sapply(objectList, 
                                                     function(x) x$ncol), 
                                              sapply(newDatasets, 
                                                     function(x) x$filename),
                                              sapply(newDatasets, 
                                                     function(x) x$valuePath),
                                              sapply(newDatasets, 
                                                     function(x) x$rowindPath),
                                              sapply(newDatasets, 
                                                     function(x) x$colptrPath),
                                              sapply(newDatasets, 
                                                     function(x) x$nrow),
                                              sapply(newDatasets, 
                                                     function(x) x$ncol), 
                                              Vinit, Winit, Ainit, Binit,
                                              k, lambda, project, maxEpoch,
                                              minibatchSize, maxHALSIter, 
                                              verbose)
        )
    }
    return(res)
}

#' Perform Mosaic Integrative Non-negative Matrix Factorization with Unshared 
#' Features
#' @description
#' Performs mosaic integrative non-negative matrix factorization (iNMF) (A.R. 
#' Kriebel, 2022) to return factorized \eqn{H}, \eqn{W}, \eqn{V} and \eqn{U} 
#' matrices. The objective function is stated as
#'
#' \deqn{\arg\min_{H\ge0,W\ge0,V\ge0,U\ge0}\sum_{i}^{d}
#' ||\binom{E_i}{P_i}-(\binom{W}{0}+\binom{V_i}{U_i})Hi||^2_F+
#' \lambda_i\sum_{i}^{d}||\binom{V_i}{U_i}H_i||_F^2}
#'
#' where \eqn{E_i} is the input non-negative matrix of the \eqn{i}'th dataset, 
#' \eqn{P_i} is the input non-negative matrix for the unshared features,
#' \eqn{d} is the total number of datasets. \eqn{E_i} is of size 
#' \eqn{m \times n_i} for \eqn{m} shared features and \eqn{n_i} sample points, 
#' \eqn{P_i} is of size \eqn{u_i \times n_i} for \eqn{u_i} unshared feaetures,
#' \eqn{H_i} is of size \eqn{n_i \times k}, \eqn{V_i} is of size 
#' \eqn{m \times k}, \eqn{W} is of size \eqn{m \times k} and \eqn{U_i} is of 
#' size \eqn{u_i \times k}.
#' 
#' Similar to \code{\link{inmf}}, \code{uinmf} also optimizes the objective with
#' ANLS algorithm.
#' @param objectList list of input datasets. List elements should all be of the
#' same class. Viable classes include: matrix, \linkS4class{dgCMatrix}, 
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
#' @param verbose Logical scalar. Whether to show information and progress. 
#' Default \code{TRUE}.
#' @return A list of the following entries: \code{H} - a list of result 
#' \eqn{H_i} matrices, \code{V} - a list of result \eqn{V_i} matrices, 
#' \code{W} - the result \eqn{W} matrix, \code{U} - the result \eqn{U_i} 
#' matrices, \code{objErr} - the final objective error value.
#' @author Yichen Wang
#' @references April R. Kriebel and Joshua D. Welch, UINMF performs mosaic 
#' integration of single-cell multi-omic datasets using nonnegative matrix 
#' factorization, Nat. Comm., 2022
#' @export
#' @examples
#' if (FALSE) {
#'     set.seed(1)
#'     res1 <- uinmf(list(dense1, dense2), 
#'                   list(dense.unshare1, dense.unshare2))
#'     set.seed(1)
#'     res2 <- uinmf(list(sparse1, sparse2), 
#'                   list(sparse.unshare1, sparse.unshare2))
#' }
uinmf <- function(
    objectList,
    unsharedList,
    k = 20,
    lambda = 5,
    niter = 30,
    verbose = TRUE
) {
    if (length(lambda) == 1) lambda <- rep(lambda, length(objectList))
    if (length(lambda) != length(objectList)) {
        stop("Must specify 1 lambda for all or each.")
    }
    mode <- .typeOfInput(objectList, null.rm = FALSE)
    unsharedList <- .uinmf.matchDatasets(objectList, unsharedList)
    res <- switch(
        mode,
        matrix = uinmf_rcpp(objectList, unsharedList, 
                            k, lambda, niter, verbose),
        dgCMatrix = uinmf_rcpp(objectList, unsharedList, 
                               k, lambda, niter, verbose),
        H5Mat = uinmf_h5dense(sapply(objectList, function(x) x$filename),
                              sapply(objectList, function(x) x$dataPath),
                              sapply(unsharedList, function(x) x$filename),
                              sapply(unsharedList, function(x) x$dataPath),
                              k, lambda, niter, verbose),
        H5SpMat = uinmf_h5sparse(sapply(objectList, function(x) x$filename),
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
                                 k, lambda, niter, verbose)
    )
    return(res)
}

.uinmf.matchDatasets <- function(objectList, unsharedList) {
    if (is.null(names(objectList)) || is.null(unsharedList)) {
        if (length(unsharedList) != length(objectList)) {
            stop("Number of matrix in unshared featire list does not match, ", 
                 "please use named lists to indicate dataset matching.")
        }
    } else {
        unsharedList <- lapply(names(objectList), function(n) unsharedList[[n]])
        names(unsharedList) <- names(objectList)
    }
    mode <- .typeOfInput(objectList, null.rm = FALSE)
    mode2 <- .typeOfInput(unsharedList, null.rm = TRUE)
    if (mode != mode2) {
        stop("Data of unshared feature should be of the same class as ", 
             "input data")
    }
    for (i in seq_along(unsharedList)) {
        if (is.null(unsharedList[[i]])) {
            # Assume no unshared feature for this dataset
            unsharedList[[i]] <- switch(
                mode,
                matrix = matrix(nrow = 0, ncol = ncol(objectList[[i]])),
                dgCMatrix = as(as(as(Matrix::Matrix(
                    nrow = 0, ncol = ncol(objectList[[i]])
                ), "dMatrix"), "generalMatrix"), "CsparseMatrix"),
                H5Mat = stop("Empty H5Mat not supported yet"),
                H5SpMat = stop("Empty H5SpMat not supported yet")
            )
        } else {
            if (ncol(unsharedList[[i]]) != ncol(objectList[[i]])) {
                stop("Number of columns in each matrix from `unsharedList` ",
                     "must match with the corresponding matrix from `object`")
            }
        }
    }
    return(unsharedList)
}
