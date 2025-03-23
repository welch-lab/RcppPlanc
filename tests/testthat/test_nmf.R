source("helper-nmf.R")
set.seed(1)
m <- 100
n <- 300
k <- 10
mat <- matrix(runif(m*n, 0, 1), nrow = m, ncol = n)

test_that("dense, nmf, anlsbpp", {
    skip_on_winbuilder()
    # Working use
    res <- nmf(mat, k, niter = 100)
    expect_type(res, "list")
    expect_equal(nrow(res$W), m)
    expect_equal(ncol(res$W), k)
    expect_equal(nrow(res$H), n)
    expect_equal(ncol(res$H), k)
    expect_lte(res$objErr, 2025)
    expect_gte(res$objErr, 1)

    # Using init W and H
    res <- nmf(mat, k, niter = 100, Winit = res$W, Hinit = res$H, nCores = 2)
    expect_lte(res$objErr, 2025)
    expect_gte(res$objErr, 1)

    # Failing use
    expect_error({
      nmf(mat, k, algo = "hello")
    }, "Please choose `algo` from")


})

test_that("dense, nmf, admm", {
    skip_on_winbuilder()
    res <- nmf(mat, k, niter = 100, algo = "admm", nCores = 2)
    expect_lte(res$objErr, 2025)
    expect_gte(res$objErr, 1)

})

test_that("dense, nmf, hals", {
    skip_on_winbuilder()
    res <- nmf(mat, k, niter = 100, algo = "hals", nCores = 2)
    expect_lte(res$objErr, 2030)
    expect_gte(res$objErr, 1)
})

test_that("dense, nmf, mu", {
    skip_on_winbuilder()
    res <- nmf(mat, k, niter = 100, algo = "mu", nCores = 2)
    expect_lte(res$objErr, 2650)
    expect_gte(res$objErr, 1)
})

symmat <- t(mat) %*% mat
lambda <- 5

test_that("dense, symNMF, anlsbpp", {
    skip_on_winbuilder()
    res <- symNMF(symmat, k, niter = 100, lambda = lambda, algo = "anlsbpp", nCores = 2)
    expect_lte(res$objErr, 4e7)

    res <- symNMF(symmat, k, niter = 100, lambda = lambda, algo = "anlsbpp",
                  Hinit = res$H, nCores = 2)
    expect_lte(res$objErr, 4e7)

    expect_error({
      symNMF(mat, k, 100)
    }, "Input `x` is not square.")

    expect_error({
      symNMF(symmat, 1e4, 100)
    })

    expect_error({
      symNMF(symmat, k, 100, Hinit = t(res$H))
    }, "Hinit must be of size ")

    expect_error({
      symNMF(symmat, k, 100, algo = "hello")
    }, "Please choose `algo` from")
})

test_that("dense, symNMF, gnsym", {
  skip_on_winbuilder()
  res <- symNMF(symmat, k, niter = 100, algo = "gnsym", nCores = 2)
  expect_lte(res$objErr, 5.8e4)
})

library(Matrix)
# Sparsen the `mat`
sparsity <- .9
regenerate <- TRUE
while (regenerate) {
    zero.idx <- sample(length(mat), round(sparsity * length(mat)))
    mat.sp <- mat
    mat.sp[zero.idx] <- 0
    mat.sp <- as(mat.sp, "CsparseMatrix")
    # Make sure there is no col/row that has all zero
    if (sum(Matrix::colSums(mat.sp) == 0) == 0 &&
        sum(Matrix::rowSums(mat.sp) == 0) == 0) {
      regenerate <- FALSE
    }
}

test_that("sparse, nmf, anlsbpp", {
  skip_on_winbuilder()
  set.seed(1)
  res1 <- nmf(mat.sp, k, niter = 100, nCores = 2)
  expect_type(res1, "list")
  expect_equal(nrow(res1$W), m)
  expect_equal(ncol(res1$W), k)
  expect_equal(nrow(res1$H), n)
  expect_equal(ncol(res1$H), k)
  expect_lte(res1$objErr, 800)
  expect_gte(res1$objErr, 1)
  set.seed(1)
  res2 <- nmf(as.matrix(mat.sp), k, niter = 100, nCores = 2)
  expect_true(all.equal(res1, res2))
  # Using init W and H
  res <- nmf(mat.sp, k, niter = 100, Winit = res1$W, Hinit = res1$H, nCores = 2)
  # Expected max objective error
  expect_lte(res$objErr, 800)
  expect_gte(res$objErr, 1)
  # Expected min sparsity of W
  W.sparsity <- sum(res$W == 0) / length(res$W)
  cat("\nW sparsity:",W.sparsity,"\n")
  expect_gte(W.sparsity, .4)

  expect_error({
    nmf(mat, k, algo = "hello")
  }, "Please choose `algo` from")
})

test_that("sparse, nmf, admm", {
  skip_on_winbuilder()
  res <- nmf(mat.sp, k, niter = 100, algo = "admm", nCores = 2)
  expect_lte(res$objErr, 800)
  expect_gte(res$objErr, 1)
})

test_that("sparse, nmf, hals", {
  skip_on_winbuilder()
  res <- nmf(mat.sp, k, niter = 100, algo = "hals", nCores = 2)
  expect_lte(res$objErr, 800)
  expect_gte(res$objErr, 1)
})

test_that("sparse, nmf, mu", {
  skip_on_winbuilder()
  res <- nmf(mat.sp, k, niter = 100, algo = "mu", nCores = 2)
  expect_lte(res$objErr, 950)
  expect_gte(res$objErr, 1)
})

symmat.sp <- Matrix::t(mat.sp) %*% mat.sp

test_that("sparse, symNMF, anlsbpp", {
  skip_on_winbuilder()
  res <- symNMF(symmat.sp, k, niter = 100, lambda = lambda, algo = "anlsbpp", nCores = 2)
  expect_lte(res$objErr, 4e5)
})

test_that("sparse, symNMF, gnsym", {
  skip_on_winbuilder()
  res <- symNMF(symmat.sp, k, niter = 100, algo = "gnsym", nCores = 2)
  expect_lte(res$objErr, 1e4)
})
