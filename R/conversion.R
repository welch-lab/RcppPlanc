.rcpp_mat_to_h5mat <- function(x, filename, dataPath, overwrite = TRUE, ...) {
  res <- c(filename, dataPath)
  # Set the dimensions of the dataset
  hdf5r.Extra::h5Write(x = x, file = filename, name = dataPath, overwrite = overwrite, transpose = FALSE, ...)
  return(res)
}

.rcpp_spmat_to_h5mat <- function(x, filename, dataPath, overwrite = TRUE, ...)
{
  res <- c(filename, dataPath)
  x <- as.matrix(x)
  hdf5r.Extra::h5Write(x = x, file = filename, name = dataPath, overwrite = overwrite, transpose = FALSE, ...)
  return(res)
}

.rcpp_spmat_to_h5spmat <- function(x, filename, dataPath, overwrite = TRUE) {
  res <- c(filename, paste0(dataPath, "/data"), paste0(dataPath, "/indices"), paste0(dataPath, "/indptr"))
  hdf5r.Extra::h5Write(x = x, file = filename, name = dataPath, overwrite = overwrite, transpose = FALSE)
  return(res)
}
