#' Create argument list for constructing HDF5 based dense matrix
#' @param filename Filename of the HDF5 file
#' @param dataPath Path in the HDF5 file that points to a 2D dense matrix
#' @return H5Mat object, can be used like a list
H5Mat <- function(
    filename,
    dataPath
) {
    if (!file.exists(filename)) stop("File not found: ", filename)
    argList <- list(filename = filename, dataPath = dataPath)
    class(argList) <- "H5Mat"
    return(argList)
}

#' Show info of H5Mat object
#' @param x H5Mat object
#' @param ... not used
#' @method print H5Mat
print.H5Mat <- function(x, ...) {
    cat("Argument list for constructing HDF5 dense matrix\n",
        "filename:  ", x$filename, "\n",
        "data path: ", x$dataPath, sep = "")
}

#' Create argument list for constructing HDF5 based CSC sparse matrix
#' @param filename Filename of the HDF5 file
#' @param valuePath Path in the HDF5 file that points to the 1D value vector of
#' the sparse matrix.
#' @param rowindPath Path in the HDF5 file that points to the 1D rowind vector
#' of the sparse matrix.
#' @param colptrPath Path in the HDF5 file that points to the 1D colptr vector
#' of the sparse matrix.
#' @param ncol Integer, number of columns of the sparse matrix.
#' @param nrow Integer, number of rows of the sparse matrix.
#' @return H5SpMat object, can be used like a list
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

#' Show info of H5SpMat object
#' @param x H5SpMat object
#' @param ... not used
#' @method print H5SpMat
print.H5SpMat <- function(x, ...) {
  cat("Argument list for constructing HDF5 CSC sparse matrix\n",
      "filename:    ", x$filename, "\n",
      "value path:  ", x$valuePath, "\n",
      "rowind path: ", x$rowindPath, "\n",
      "colptr path: ", x$colptrPath, "\n",
      "dimension:   ", x$nrow, " x ", x$ncol, sep = "")
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
