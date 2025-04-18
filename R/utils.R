#' Argument list object for using a dense matrix stored in HDF5 file
#' @description
#' For running \code{\link{inmf}}, \code{\link{onlineINMF}} or
#' \code{\link{uinmf}} with dense matrix stored in HDF5 file, users will need
#' to construct an argument list for the filename of the HDF5 file as well as
#' the path in the file storing the matrix. \code{H5Mat} is provided as an
#' instructed constructor. Meanwhile, since the INMF functions require that
#' all datasets should be of the same type, \code{as.H5Mat} is provided for
#' writing in-memory data into a new HDF5 file on disk and returning the
#' constructed argument list.
#'
#' @param x For \code{as.H5Mat}, matrix of either dense or sparse type to be
#' written; for \code{print}, a \code{H5Mat} argument list object
#' @param filename Filename of the HDF5 file
#' @param dataPath Path in the HDF5 file that points to a 2D dense matrix.
#' Default \code{"data"} when using \code{as.H5Mat}.
#' @param overwrite Logical, whether to overwrite the file if already exists at
#' the given path. Default \code{FALSE}.
#' @param ... Passed down to hdf5r.Extra::h5Write
#' @rdname H5Mat
#' @return H5Mat object, indeed a list object.
#' @examples
#' if (require("withr")) {
#' H5MatEx <- function(){
#' withr::local_dir(withr::local_tempdir())
#' h <- H5Mat(system.file("extdata/ctrl_dense.h5", package = "RcppPlanc"),
#'            "data")
#' print(h)
#'
#' library(Matrix)
#' ctrl.dense <- as.matrix(ctrl.sparse)
#' h1 <- as.H5Mat(ctrl.dense, "ctrl_from_dense_to_dense.h5",
#'                dataPath = "data")
#' h1
#' h2 <- as.H5Mat(ctrl.sparse, "ctrl_from_sparse_to_dense.h5",
#'                dataPath = "data")
#' }
#' H5MatEx()
#' }
H5Mat <- function(
  filename,
  dataPath
) {
  if (!file.exists(filename)) stop("File not found: ", filename)
  argList <- list(filename = filename, dataPath = dataPath)
  class(argList) <- "H5Mat"
  return(argList)
}

#' @export
#' @rdname H5Mat
as.H5Mat <- function(x, filename, dataPath = "data", overwrite = FALSE, ...) {
  if (isFALSE(overwrite) && file.exists(filename)) {
    filename <- normalizePath(filename)
    stop("File already exists at the given path: ", filename)
  }
  UseMethod("as.H5Mat", x)
}

#' @export
#' @rdname H5Mat
#' @method as.H5Mat matrix
as.H5Mat.matrix <- function(x, filename, dataPath = "data", overwrite, ...) {
  res <- .rcpp_mat_to_h5mat(x, filename, dataPath, overwrite = overwrite, ...)
  H5Mat(res[1], res[2])
}

#' @export
#' @rdname H5Mat
#' @method as.H5Mat dgCMatrix
as.H5Mat.dgCMatrix <- function(x, filename, dataPath = "data",
                               overwrite = FALSE, ...) {
  res <- .rcpp_spmat_to_h5mat(x, filename, dataPath, overwrite = overwrite, ...)
  H5Mat(res[1], res[2])
}

#' @export
#' @rdname H5Mat
#' @method as.H5Mat default
as.H5Mat.default <- function(x, filename, dataPath = "data", ...) {
  if (methods::canCoerce(x, "matrix")) x <- as.matrix(x)
  else if (methods::canCoerce(x, "CsparseMatrix")) x <- methods::as(x, 'CsparseMatrix')
  else {
    stop("No method for coercing \"", class(x)[1], "\" to \"matrix\" or ",
         "\"dgCMatrix\"")
  }
  as.H5Mat(x, filename, dataPath, ...)
}

#' Prepare character information of a H5Mat object
#' @param x H5Mat argument list object
#' @param ... Not used.
#' @return A character scalar of the displayed message
#' @method format H5Mat
#' @export
#' @examples
#' h <- H5Mat(system.file("extdata/ctrl_dense.h5", package = "RcppPlanc"),
#'            "data")
#' format(h)
format.H5Mat <- function(x, ...) {
  msg <- paste0(
    "Argument list for constructing HDF5 dense matrix\n",
    "filename:  ", x$filename, "\n",
    "data path: ", x$dataPath, "\n"
  )
  return(msg)
}

#' Show information of a H5Mat object
#' @param x H5Mat argument list object
#' @param ... Not used.
#' @return NULL. Information displayed.
#' @method print H5Mat
#' @export
#' @examples
#' h <- H5Mat(system.file("extdata/ctrl_dense.h5", package = "RcppPlanc"),
#'            "data")
#' print(h)
print.H5Mat <- function(x, ...) {
  cat(format.H5Mat(x, ...))
}


#' Argument list object for using a sparse matrix stored in HDF5 file
#' @description
#' For running \code{\link{inmf}}, \code{\link{onlineINMF}} or
#' \code{\link{uinmf}} with sparse matrix stored in HDF5 file, users will need
#' to construct an argument list for the filename of the HDF5 file as well as
#' the paths in the file storing the arrays that construct the CSC (compressed
#' sparse column) matrix. \code{H5SpMat} is provided as an instructed
#' constructor. Meanwhile, since the INMF functions require that all datasets
#' should be of the same type, \code{as.H5SpMat} is provided for writing
#' in-memory data into a new HDF5 file on disk and returning the
#' constructed argument list.
#'
#' @param x For \code{as.H5SpMat}, matrix of either dense or sparse type to be
#' written; for \code{print}, a \code{H5SpMat} argument list object.
#' @param filename Filename of the HDF5 file
#' @param valuePath Path in the HDF5 file that points to a 1D array storing the
#' non-zero values of the sparse matrix. Default \code{"data"} when using
#' \code{as.H5SpMat}.
#' @param rowindPath Path in the HDF5 file that points to a 1D integer array
#' storing the row indices of non-zero values in each column of the sparse
#' matrix. Default \code{"indices"} when using \code{as.H5SpMat}.
#' @param colptrPath Path in the HDF5 file that points to a 1D integer array
#' storing the number of non-zero values in each column of the sparse matrix.
#' Default \code{"indptr"} when using \code{as.H5SpMat}.
#' @param dataPath For \code{as.H5SpMat} methods, the H5Group name for the
#' sparse matrix. Default \code{""}.
#' @param nrow,ncol Integer, the true dimensionality of the sparse matrix.
#' @param overwrite Logical, whether to overwrite the file if already exists at
#' the given path. Default \code{FALSE}.
#' @param ... not used
#' @rdname H5SpMat
#' @return H5SpMat object, indeed a list object.
#' @examples
#' if (require("withr")) {
#'   H5SpMatEx <- function() {
#'     withr::local_dir(withr::local_tempdir())
#'     h <- H5SpMat(system.file("extdata/ctrl_sparse.h5", package = "RcppPlanc"),
#'          "scaleDataSparse/data", "scaleDataSparse/indices", "scaleDataSparse/indptr", 173, 300)
#'     dim(h)
#'
#'     library(Matrix)
#'     ctrl.dense <- as.matrix(ctrl.sparse)
#'     h1 <- as.H5SpMat(ctrl.sparse, "ctrl_from_sparse_to_sparse.h5", "matrix")
#'     h1
#'     h2 <- as.H5SpMat(ctrl.dense, "ctrl_from_dense_to_sparse.h5", "matrix")
#'     h2
#'   }
#'   H5SpMatEx()
#' }
H5SpMat <- function(
  filename,
  valuePath,
  rowindPath,
  colptrPath,
  nrow,
  ncol
) {
  if (!file.exists(filename)) stop("File not found: ", filename)
  argList <- list(
    filename = filename,
    valuePath = valuePath,
    rowindPath = rowindPath,
    colptrPath = colptrPath,
    nrow = nrow,
    ncol = ncol
  )
  class(argList) <- "H5SpMat"
  return(argList)
}

#' @export
#' @rdname H5SpMat
as.H5SpMat <- function(x, filename, dataPath,
                       overwrite = FALSE) {
  if (isFALSE(overwrite) && file.exists(filename)) {
    filename <- normalizePath(filename)
    stop("File already exists at the given path: ", filename)
  }
  UseMethod("as.H5SpMat", x)
}

#' @export
#' @rdname H5SpMat
#' @method as.H5SpMat matrix
as.H5SpMat.matrix <- function(x, filename, dataPath = "",
                              overwrite = FALSE) {
  x <- methods::as(x, "CsparseMatrix")
  as.H5SpMat.dgCMatrix(x, filename = filename, dataPath = dataPath,
                       overwrite = overwrite)
}

#' @export
#' @rdname H5SpMat
#' @method as.H5SpMat dgCMatrix
as.H5SpMat.dgCMatrix <- function(x, filename, dataPath = "",
                                 overwrite = FALSE) {
  res <- .rcpp_spmat_to_h5spmat(x, filename,
                               dataPath, overwrite)
  H5SpMat(res[1], res[2], res[3], res[4], nrow(x), ncol(x))
}

#' @export
#' @rdname H5SpMat
#' @method as.H5SpMat default
as.H5SpMat.default <- function(x, filename, dataPath = "",
                               overwrite = FALSE, ...) {
  if (methods::canCoerce(x, "CsparseMatrix")) {
    x <- methods::as(x, "CsparseMatrix")
  } else if (methods::canCoerce(x, "matrix")) {
    x <- as.matrix(x)
    x <- methods::as(x, "CsparseMatrix")
  } else {
    stop("No method for coercing \"", class(x)[1], "\" to \"dgCMatrix\"")
  }
  as.H5SpMat(x, filename = filename, dataPath = dataPath,
             overwrite = overwrite, ...)
}

#' prepare character information of a H5SpMat object
#' @method format H5SpMat
#' @param x H5SpMat argument list object
#' @param ... Not used.
#' @return A character scalar of the displayed message
#' @export
#' @examples
#' h <- H5SpMat(system.file("extdata/ctrl_sparse.h5", package = "RcppPlanc"),
#'              "data", "indices", "indptr", 173, 300)
#' format(h)
format.H5SpMat <- function(x, ...) {
  data_path <- x$valuePath
  indices_path <- x$rowindPath
  indptr_path <- x$colptrPath
  if (dirname(data_path) == dirname(indices_path) &&
      dirname(data_path) == dirname(indptr_path)) {
    msg <- paste0(
    "Argument list for constructing HDF5 CSC sparse matrix\n",
    "filename:    ", x$filename, "\n",
    "data path:  ", dirname(x$valuePath), "\n",
    "dimension:   ", x$nrow, " x ", x$ncol, "\n"
    )
  } else {
    msg <- paste0(
    "Argument list for constructing HDF5 CSC sparse matrix\n",
    "filename:    ", x$filename, "\n",
    "value path:  ", x$valuePath, "\n",
    "rowind path: ", x$rowindPath, "\n",
    "colptr path: ", x$colptrPath, "\n",
    "dimension:   ", x$nrow, " x ", x$ncol, "\n"
    )
  }
  return(msg)
}

#' Show information of a H5SpMat object
#' @method print H5SpMat
#' @param x H5SpMat argument list object
#' @param ... Not used.
#' @return NULL. Information displayed.
#' @export
#' @examples
#' h <- H5SpMat(system.file("extdata/ctrl_sparse.h5", package = "RcppPlanc"),
#'              "data", "indices", "indptr", 173, 300)
#' print(h)
print.H5SpMat <- function(x, ...) {
  cat(format.H5SpMat(x, ...))
}

#' Retrieve the dimension of H5SpMat argument list
#' @method dim H5SpMat
#' @param x H5SpMat argument list object
#' @param value Numeric vector of two, for number of rows and number of columns.
#' @return Retriever returns a vector of two (nrow and ncol), setter sets the
#' value of that in the argument list.
#' @rdname dim-H5SpMat
#' @export
#' @examples
#' h <- H5SpMat(system.file("extdata/ctrl_sparse.h5", package = "RcppPlanc"),
#'              "data", "indices", "indptr", 173, 300)
#' dim(h)
#' nrow(h)
#' ncol(h)
#' dim(h) <- c(200, 200)
#' h
dim.H5SpMat <- function(x) { return(c(x$nrow, x$ncol)) }

#' @rdname dim-H5SpMat
#' @export
`dim<-.H5SpMat` <- function(x, value) {
  x$nrow <- value[1]
  x$ncol <- value[2]
  return(x)
}

.typeOfInput <- function(objectList, null.rm = TRUE) {
  classes <- sapply(objectList, function(x) class(x)[1])
  if (isTRUE(null.rm)) classes <- classes[classes != "NULL"]
  if (!all(classes == classes[1])) {
    stop("All datasets should be of the same class of input")
  }
  if (!classes[1] %in% c("matrix", "dgCMatrix", "H5Mat", "H5SpMat")) {
    stop("Datasets of class `", classes[1], "` is currently not supported.")
  }
  return(classes[1])
}

.makeUnsharedIdx <- function(EiNames, PiNames) {
  res <- rep(-1L, length(EiNames))
  for (i in seq_along(EiNames)) {
    res[i] <- match(EiNames[i], PiNames) - 1
  }
  return(res)
}
