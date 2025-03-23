library(RcppPlanc)
library(testthat)
library(Matrix)
local_edition(3)
source("helper-nmf.R")
ctrl.h5ds <- H5Mat(filename = system.file("extdata", "ctrl_dense.h5",
                                          package = "RcppPlanc", mustWork = TRUE),
                   dataPath = "scaleData")
stim.h5ds <- H5Mat(filename = system.file("extdata", "stim_dense.h5",
                                          package = "RcppPlanc", mustWork = TRUE),
                   dataPath = "scaleData")
ctrl.h5sp <- H5SpMat(filename = system.file("extdata", "ctrl_sparse.h5",
                                            package = "RcppPlanc", mustWork = TRUE),
                     valuePath = "scaleDataSparse/data",
                     rowindPath = "scaleDataSparse/indices",
                     colptrPath = "scaleDataSparse/indptr",
                     nrow = nrow(ctrl.sparse),
                     ncol = ncol(ctrl.sparse))
stim.h5sp <- H5SpMat(filename = system.file("extdata", "stim_sparse.h5",
                                            package = "RcppPlanc", mustWork = TRUE),
                     valuePath = "scaleDataSparse/data",
                     rowindPath = "scaleDataSparse/indices",
                     colptrPath = "scaleDataSparse/indptr",
                     nrow = nrow(stim.sparse),
                     ncol = ncol(stim.sparse))
ctrl.dense <- as.matrix(ctrl.sparse)
stim.dense <- as.matrix(stim.sparse)
m <- nrow(ctrl.dense)
ni <- c(ncol(ctrl.dense), ncol(stim.dense))
k <- 20

test_that("inmf, w/o init", {
  skip_on_winbuilder()
  # Dense cases
  set.seed(1)
  res1 <- inmf(list(ctrl.dense, stim.dense), k = k, nCores = 2)
  expect_length(res1, 4)
  expect_length(res1$H, 2)
  expect_length(res1$V, 2)
  expect_true(all.equal(dim(res1$W), c(m, k)))
  expect_lte(res1$objErr, 3.8e4)

  # Error tests
  expect_error(inmf(list(ctrl.dense, stim.sparse)), 
               "All datasets should be of the same class")
  
  expect_error(inmf(list(ctrl.dense, ctrl.dense), k = 300),
               "k must be <= m")
  # Sparse cases
  set.seed(1)
  res2 <- inmf(list(ctrl.sparse, stim.sparse), k = k, nCores = 2)
  expect_true(all.equal(res1, res2))
  
  # dense h5 cases
  set.seed(1)
  res3 <- inmf(list(ctrl.h5ds, stim.h5ds), k = k)
  expect_true(all.equal(res1, res3))
  
  # sparse h5 cases
  set.seed(1)
  res4 <- inmf(objectList = list(ctrl.h5sp, stim.h5sp), k = k, nCores = 2)
  expect_true(all.equal(res1, res4, tolerance = 1e-6))
})

test_that("inmf, w/ init", {
  skip_on_winbuilder()
  set.seed(1)
  res0 <- inmf(list(ctrl.sparse, stim.sparse), k = k, nCores = 2)
  res1 <- inmf(list(ctrl.sparse, stim.sparse), k = k, Hinit = res0$H, Vinit = res0$V, Winit = res0$W, nCores = 2)
  res2 <- inmf(list(ctrl.dense, stim.dense), k = k, Hinit = res0$H, Vinit = res0$V, Winit = res0$W,, nCores = 2)
  expect_true(all.equal(res1, res2))
  
  expect_error(inmf(list(ctrl.sparse, stim.sparse), k = k, 
                    Hinit = lapply(res0$H, t), Vinit = res0$V, Winit = res0$W),
               "Each given H must be of size")
  expect_error(inmf(list(ctrl.sparse, stim.sparse), k = k, 
                    Hinit = res0$H, Vinit = lapply(res0$V, t), Winit = res0$W),
               "All given Vs must be of size")
  expect_error(inmf(list(ctrl.sparse, stim.sparse), k = k, 
                    Hinit = res0$H, Vinit = res0$V, Winit = t(res0$W)),
               "Given W must be of size")
})

test_that("onlineINMF, scenario 1", {
  skip_on_winbuilder()
  set.seed(1)
  res1 <- onlineINMF(list(ctrl.sparse, stim.sparse), k = k, minibatchSize = 50, nCores = 2)
  expect_length(res1, 6)
  expect_length(res1$H, 2)
  expect_length(res1$V, 2)
  expect_length(res1$A, 2)
  expect_length(res1$B, 2)
  expect_true(all.equal(dim(res1$W), c(m, k)))
  expect_lte(res1$objErr, 3.8e4)
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), k = k),
               "Please set a smaller")
  expect_error(onlineINMF(list(ctrl.sparse, stim.dense), k = k),
               "All datasets should be of the same class")
  # TODO could be more error cases
  set.seed(1)
  res2 <- onlineINMF(list(ctrl.dense, stim.dense), k = k, minibatchSize = 50, nCores = 2)
  expect_true(all.equal(res1, res2))
  # dense h5 cases
  set.seed(1)
  res3 <- onlineINMF(list(ctrl.h5ds, stim.h5ds), k = k, minibatchSize = 50, nCores = 2)
  expect_true(all.equal(res1, res3))
  # sparse h5 cases
  set.seed(1)
  res4 <- onlineINMF(list(ctrl.h5sp, stim.h5sp), k = k, minibatchSize = 50, nCores = 2)
  expect_true(all.equal(res1, res4, tolerance = 1e-6))
})
set.seed(233)
new.data <- ctrl.sparse
new.data@x <- new.data@x + runif(length(new.data@x))
new.data.dense <- as.matrix(new.data)
new.h5sp <- as.H5SpMat(new.data, "temp_new_sparse.h5", overwrite = TRUE)
new.h5m <- as.H5Mat(new.data.dense, "temp_new_dense.h5", overwrite = TRUE)
set.seed(42)
res0 <- onlineINMF(list(ctrl.sparse, stim.sparse), k = k, minibatchSize = 50)
# TODO create temp h5 file for new.data
test_that("onlineINMF, scenario 2", {
  skip_on_winbuilder()
  set.seed(1)
  res1 <- onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                     k = k, minibatchSize = 50, permuteChunkSize = 32,
                     Hinit = res0$H, Vinit = res0$V, Winit = res0$W,
                     Ainit = res0$A, Binit = res0$B, nCores = 2)
  set.seed(1)
  res2 <- onlineINMF(list(ctrl.dense, stim.dense), newDatasets = list(new.data.dense),
                     k = k, minibatchSize = 50, permuteChunkSize = 32,
                     Hinit = res0$H, Vinit = res0$V, Winit = res0$W,
                     Ainit = res0$A, Binit = res0$B, nCores = 2)
  expect_true(all.equal(res1, res2))
  expect_error(onlineINMF(list(ctrl.dense, stim.dense), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Hinit = res0$H, 
                          Vinit = res0$V, Winit = res0$W,
                          Ainit = res0$A, Binit = res0$B),
               "newDatasets should be of the same class as original datasets")
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Hinit = res0$H,
                          Vinit = lapply(res0$V, t), Winit = res0$W,
                          Ainit = res0$A, Binit = res0$B),
               "All given Vs must be of size 173 x 20")
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Hinit = res0$H,
                          Vinit = res0$V, Winit = t(res0$W),
                          Ainit = res0$A, Binit = res0$B),
               "Given W must be of size 173 x 20 but is 20 x 173")
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Hinit = res0$H,
                          Vinit = res0$V, Winit = res0$W,
                          Ainit = res0$A, Binit = lapply(res0$B, t)),
               "Given Bs must all be of size 173 x 20")
  # A's are square, and nearly symmetric (but not exactly), don't have a smart check here yet
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Hinit = res0$H, 
                          Winit = res0$W,
                          Ainit = res0$A, Binit = res0$B),
               "Must provide 2 V matrices")
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Hinit = res0$H, 
                          Vinit = res0$V, Winit = res0$W,
                          Binit = res0$B),
               "Must provide 2 A matrices")
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Hinit = res0$H, 
                          Vinit = res0$V, Winit = res0$W,
                          Ainit = res0$A),
               "Must provide 2 B matrices")
  set.seed(1)
  res3 <- onlineINMF(list(ctrl.h5sp, stim.h5sp), newDatasets = list(new.h5sp),
                     k = k, minibatchSize = 50, permuteChunkSize = 32,
                     Hinit = res0$H, Vinit = res0$V, Winit = res0$W,
                     Ainit = res0$A, Binit = res0$B, nCores = 2)
  expect_true(all.equal(res1, res3))
  set.seed(1)
  res4 <- onlineINMF(list(ctrl.h5ds, stim.h5ds), newDatasets = list(new.h5m),
                     k = k, minibatchSize = 50,  Hinit = res0$H, 
                     Vinit = res0$V, Winit = res0$W,
                     Ainit = res0$A, Binit = res0$B, nCores = 2)
  expect_true(all.equal(res1, res4))
})

test_that("onlineINMF, scenario 3", {
  skip_on_winbuilder()
  expect_no_error(onlineINMF(list(ctrl.sparse, stim.sparse),
                             newDatasets = list(new.data), project = TRUE,
                             k = k, minibatchSize = 50, Winit = res0$W, nCores = 2))
})

p1 <- ctrl.sparse[1:10,]
p2 <- stim.sparse[11:30,]
p1.dense <- as.matrix(p1)
p2.dense <- as.matrix(p2)
test_that("uinmf", {
 skip_on_winbuilder()
 set.seed(1)
 res1 <- uinmf(list(ctrl.sparse, stim.sparse), list(p1, p2), nCores = 2)
 expect_length(res1, 5)
 expect_length(res1$H, 2)
 expect_length(res1$V, 2)
 expect_length(res1$U, 2)
 expect_true(all.equal(dim(res1$W), c(173, 20)))
 expect_true(all.equal(dim(res1$U[[1]]), c(10, 20)))
 expect_true(all.equal(dim(res1$U[[2]]), c(20, 20)))
 expect_lte(res1$objErr, 4.5e4)
 set.seed(1)
 res2 <- uinmf(list(ctrl.dense, stim.dense), list(p1.dense, p2.dense), nCores = 2)
 expect_true(all.equal(res1, res2))
  expect_error(uinmf(list(ctrl.sparse, stim.sparse), list(p1, p2), lambda = 1:3),
               "Must specify 1 lambda for all or each.")
  expect_error(uinmf(list(ctrl.sparse, stim.sparse), list(p1)),
               "Number of matrix in unshared feature list does not match")
  expect_error(uinmf(list(ctrl.sparse, stim.sparse), list(p1.dense, p2.dense)),
               "Data of unshared feature should be of the same class as")
  expect_error(uinmf(list(ctrl.sparse, stim.sparse), list(p1, p2[,1:100])),
               "Number of columns in each matrix from")
  set.seed(1)
  res3 <- uinmf(list(a = ctrl.sparse, b = stim.sparse), list(a = p1), nCores = 2)
  set.seed(1)
  res4 <- uinmf(list(a = ctrl.dense, b = stim.dense), list(a = p1.dense), nCores = 2)
  expect_true(all.equal(res3, res4))
  expect_equal(length(res4$U), 1)
})

test_that("auxiliary", {
  skip_on_winbuilder()
  expect_error(H5Mat("thisFilenameShouldNotExist", "something"),
               "File not found")
  expect_error(as.H5Mat(ctrl.dense, "temp_new_dense.h5", overwrite = FALSE),
               "File already exists at the given path")

  # Test if `as.H5Mat.dgCMatrix` S3 method convert sparse to h5 dense correctly
  ctrl.h5ds_2 <- as.H5Mat(ctrl.sparse, "temp_ctrl_dense.h5", overwrite = TRUE, chunk_size = c(173, 300))
  set.seed(1)
  res1 <- onlineINMF(list(ctrl.sparse, stim.sparse), k = k, minibatchSize = 50, permuteChunkSize = 300)
  set.seed(1)
  res2 <- onlineINMF(list(ctrl.h5ds_2, stim.h5ds), k = k, minibatchSize = 50)
  expect_true(all.equal(res1, res2))

  expect_no_error(print.H5Mat(ctrl.h5ds))


  expect_error(H5SpMat("thisFilenameShouldNotExist", "something", "a", "b", 1, 2),
               "File not found")
  expect_error(as.H5SpMat(ctrl.dense, "temp_new_sparse.h5", overwrite = FALSE),
               "File already exists at the given path")

  # Test if `as.H5SpMat.matrix` S3 method convert dense to h5 sparse correctly
  ctrl.h5sp_2 <- as.H5SpMat(ctrl.dense, filename = "temp_ctrl_sparse.h5")
  set.seed(1)
  res1 <- onlineINMF(list(ctrl.sparse, stim.sparse), k = k, minibatchSize = 50)
  set.seed(1)
  res2 <- onlineINMF(list(ctrl.h5sp_2, stim.h5sp), k = k, minibatchSize = 50)
  expect_true(all.equal(res1, res2, tolerance = 1e-6))

  expect_no_error(print.H5SpMat(ctrl.h5sp))
})

unlink("temp_ctrl_dense.h5")
unlink("temp_ctrl_sparse.h5")
unlink("temp_new_sparse.h5")
unlink("temp_new_dense.h5")
