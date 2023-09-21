ctrl.dense <- as.matrix(ctrl.sparse)
stim.dense <- as.matrix(stim.sparse)
m <- nrow(ctrl.dense)
ni <- c(ncol(ctrl.dense), ncol(stim.dense))
k <- 20
set.seed(1)

## TODO: Test HDF5

test_that("inmf, w/o init", {
  # Dense cases
  set.seed(1)
  res1 <- inmf(list(ctrl.dense, stim.dense), k = k)
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
  res2 <- inmf(list(ctrl.sparse, stim.sparse), k = k)
  expect_true(all.equal(res1, res2))
})

test_that("inmf, w/ init", {
  set.seed(1)
  res0 <- inmf(list(ctrl.sparse, stim.sparse), k = k)
  res1 <- inmf(list(ctrl.sparse, stim.sparse), k = k, Hinit = res0$H, Vinit = res0$V, Winit = res0$W)
  res2 <- inmf(list(ctrl.dense, stim.dense), k = k, Hinit = res0$H, Vinit = res0$V, Winit = res0$W)
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
  set.seed(1)
  res1 <- onlineINMF(list(ctrl.sparse, stim.sparse), k = k, minibatchSize = 50)
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
  res2 <- onlineINMF(list(ctrl.dense, stim.dense), k = k, minibatchSize = 50)
  expect_true(all.equal(res1, res2))
})

set.seed(233)
new.data <- ctrl.sparse
new.data@x <- new.data@x + runif(length(new.data@x))
new.data.dense <- as.matrix(new.data)
test_that("onlineINMF, scenario 2", {
  res0 <- onlineINMF(list(ctrl.sparse, stim.sparse), k = k, minibatchSize = 50)
  set.seed(1)
  res1 <- onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                     k = k, minibatchSize = 50, Vinit = res0$V, Winit = res0$W, 
                     Ainit = res0$A, Binit = res0$B)
  set.seed(1)
  res2 <- onlineINMF(list(ctrl.dense, stim.dense), newDatasets = list(new.data.dense),
                     k = k, minibatchSize = 50, Vinit = res0$V, Winit = res0$W, 
                     Ainit = res0$A, Binit = res0$B)
  expect_true(all.equal(res1, res2))
  
  expect_error(onlineINMF(list(ctrl.dense, stim.dense), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Vinit = res0$V, Winit = res0$W, 
                          Ainit = res0$A, Binit = res0$B),
               "newDatasets should be of the same class as original datasets")
  
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Vinit = lapply(res0$V, t), Winit = res0$W, 
                          Ainit = res0$A, Binit = res0$B),
               "All given Vs must be of size 173 x 20")
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Vinit = res0$V, Winit = t(res0$W), 
                          Ainit = res0$A, Binit = res0$B),
               "Given W must be of size 173 x 20 but is 20 x 173")
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Vinit = res0$V, Winit = res0$W,
                          Ainit = res0$A, Binit = lapply(res0$B, t)),
               "Given Bs must all be of size 173 x 20")
  # A's are square, and nearly symmetric (but not exactly), don't have a smart check hear yet
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Winit = res0$W, 
                          Ainit = res0$A, Binit = res0$B),
               "Must provide 2 V matrices")
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Vinit = res0$V, Winit = res0$W, 
                          Binit = res0$B),
               "Must provide 2 A matrices")
  expect_error(onlineINMF(list(ctrl.sparse, stim.sparse), newDatasets = list(new.data),
                          k = k, minibatchSize = 50, Vinit = res0$V, Winit = res0$W, 
                          Ainit = res0$A),
               "Must provide 2 B matrices")
  
  # Scenario 3
  expect_no_error(onlineINMF(list(ctrl.sparse, stim.sparse), 
                             newDatasets = list(new.data), project = TRUE,
                             k = k, minibatchSize = 50, Winit = res0$W))
})


p1 <- ctrl.sparse[1:10,]
p2 <- stim.sparse[11:30,]
p1.dense <- as.matrix(p1)
p2.dense <- as.matrix(p2)
test_that("uinmf", {
  set.seed(1)
  res1 <- uinmf(list(ctrl.sparse, stim.sparse), list(p1, p2))
  expect_length(res1, 5)
  expect_length(res1$H, 2)
  expect_length(res1$V, 2)
  expect_length(res1$U, 2)
  expect_true(all.equal(dim(res1$W), c(173, 20)))
  expect_true(all.equal(dim(res1$U[[1]]), c(10, 20)))
  expect_true(all.equal(dim(res1$U[[2]]), c(20, 20)))
  expect_lte(res1$objErr, 4.5e4)
  set.seed(1)
  res2 <- uinmf(list(ctrl.dense, stim.dense), list(p1.dense, p2.dense))
  expect_true(all.equal(res1, res2))  
  
  expect_error(uinmf(list(ctrl.sparse, stim.sparse), list(p1, p2), lambda = 1:3),
               "Must specify 1 lambda for all or each.")
  expect_error(uinmf(list(ctrl.sparse, stim.sparse), list(p1)),
               "Number of matrix in unshared featire list does not match")
  expect_error(uinmf(list(ctrl.sparse, stim.sparse), list(p1.dense, p2.dense)),
               "Data of unshared feature should be of the same class as")
  expect_error(uinmf(list(ctrl.sparse, stim.sparse), list(p1, p2[,1:100])),
               "Number of columns in each matrix from")
  set.seed(1)
  res3 <- uinmf(list(a = ctrl.sparse, b = stim.sparse), list(a = p1))
  set.seed(1)
  res4 <- uinmf(list(a = ctrl.dense, b = stim.dense), list(a = p1.dense))
  expect_true(all.equal(res3, res4))
  expect_false(is.null(res3$U[[2]]))
  expect_equal(nrow(res4$U[[2]]), 0)
  expect_equal(ncol(res4$U[[2]]), 20)
})
