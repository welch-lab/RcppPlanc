
rcpp_mat_to_h5mat <- function(x, filename, dataPath) {
  res <- list(2)
  res[1] <- filename
  res[2] <- dataPath
  # Set the dimensions of the dataset
  tDims <- list(2)
  tDims[1] <- x$ncol()
  tDims[2] <- x$nrow()
  hdf5r.Extra::h5Write(x, filename, dataPath, overwrite=TRUE, transpose=FALSE, block_size=1000L)
  return(res)
}
rcpp_spmat_to_h5mat <- function (data, rowind, colptr, nrow, ncol, filename, dataPath)
{
  res <- list(2)
  res[1] <- filename
  res[2] <- dataPath
  tDims <- list(2)
  tDims[1] <- ncol
  tDims[2] <- nrow
  x <- Matrix::sparseMatrix(i = rowind, p = colptr, x = data, dims = tDims, index1 = FALSE, repr = "C")
  x <- as.matrix(x)
  hdf5r.Extra::h5Write(x, filename, dataPath, overwrite = TRUE, transpose = FALSE)
  return(res)
}
rcpp_spmat_to_h5spmat <- function(value, rowind, colptr, nrow, ncol,
filename, dataPath) {
  res <- list(4)
  res[1] <- filename
  res[2] <- dataPath + "/x"
  res[3] <- dataPath + "/i"
  res[4] <- dataPath + "/p"
  tDims <- list(2)
  tDims[1] <- ncol
  tDims[2] <- nrow
  x <- Matrix::sparseMatrix(i = rowind, p = colptr, x = value, dims = tDims, index1 = FALSE, repr = "C")
  hdf5r.Extra::h5Write(x, filename, dataPath, overwrite = TRUE, transpose = TRUE)
  return(res)
}
