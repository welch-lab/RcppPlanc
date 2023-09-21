set.seed(1)
m <- 100
n <- 300
k <- 10
mat <- matrix(runif(m*n, 0, 1), nrow = m, ncol = n)

test_that("dense, nmf, anlsbpp", {
    # Working use
    res <- nmf(mat, k, niter = 100)
    expect_type(res, "list")
    expect_equal(nrow(res$W), m)
    expect_equal(ncol(res$W), k)
    expect_equal(nrow(res$H), n)
    expect_equal(ncol(res$H), k)
    obj <- norm(mat - res$W %*% t(res$H), type = "F")^2
    cat(obj)
    expect_lte(obj, 2025)
    
    # Using init W and H
    res <- nmf(mat, k, niter = 100, Winit = res$W, Hinit = res$H)
    obj <- norm(mat - res$W %*% t(res$H), type = "F")^2
    cat(obj)
    expect_lte(obj, 2025)
    
    # Failing use
    expect_error({
      nmf(mat, k, algo = "hello")
    }, "Please choose `algo` from")
    
    
})

test_that("dense, nmf, admm", {
    res <- nmf(mat, k, niter = 100, algo = "admm")
    obj <- norm(mat - res$W %*% t(res$H), type = "F")^2
    cat(obj)
    expect_lte(obj, 2025)
})

test_that("dense, nmf, hals", {
    res <- nmf(mat, k, niter = 100, algo = "hals")
    obj <- norm(mat - res$W %*% t(res$H), type = "F")^2
    cat(obj)
    expect_lte(obj, 2030)
})

test_that("dense, nmf, mu", {
    res <- nmf(mat, k, niter = 100, algo = "mu")
    obj <- norm(mat - res$W %*% t(res$H), type = "F")^2
    cat(obj)
    expect_lte(obj, 2650)
})

symmat <- t(mat) %*% mat
lambda <- 5

test_that("dense, symNMF, anlsbpp", {
    res <- symNMF(symmat, k, niter = 100, lambda = lambda, algo = "anlsbpp")
    obj <- norm(symmat - res$W %*% t(res$H), type = "F")^2 + 
      lambda*norm(res$W - res$H, type = "F")^2
    cat(obj)
    expect_lte(obj, 4e7)
    
    res <- symNMF(symmat, k, niter = 100, lambda = lambda, algo = "anlsbpp", 
                  Hinit = res$H)
    obj <- norm(symmat - res$W %*% t(res$H), type = "F")^2 + 
      lambda*norm(res$W - res$H, type = "F")^2
    cat(obj)
    expect_lte(obj, 4e7)
    
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
  res <- symNMF(symmat, k, niter = 100, algo = "gnsym")
  obj <- norm(symmat - res$H %*% t(res$H), type = "F")^2
  cat(obj)
  expect_lte(obj, 5.8e4)
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
  set.seed(1)
  res1 <- nmf(mat.sp, k, niter = 100)
  expect_type(res1, "list")
  expect_equal(nrow(res1$W), m)
  expect_equal(ncol(res1$W), k)
  expect_equal(nrow(res1$H), n)
  expect_equal(ncol(res1$H), k)
  obj <- Matrix::norm(mat.sp - res1$W %*% t(res1$H), type = "F")^2
  cat(obj)
  expect_lte(obj, 800)
  set.seed(1)
  res2 <- nmf(as.matrix(mat.sp), k, niter = 100)
  expect_true(all.equal(res1, res2))
  # Using init W and H
  res <- nmf(mat.sp, k, niter = 100, Winit = res1$W, Hinit = res1$H)
  obj <- Matrix::norm(mat.sp - res$W %*% t(res$H), type = "F")^2
  cat(obj)
  # Expected max objective error
  expect_lte(obj, 800)
  # Expected min sparsity of W
  W.sparsity <- sum(res$W == 0) / length(res$W)
  cat("\nW sparsity:",W.sparsity,"\n")
  expect_gte(W.sparsity, .4)
  
  expect_error({
    nmf(mat, k, algo = "hello")
  }, "Please choose `algo` from")
})

test_that("sparse, nmf, admm", {
  res <- nmf(mat.sp, k, niter = 100, algo = "admm")
  obj <- Matrix::norm(mat.sp - res$W %*% t(res$H), type = "F")^2
  cat(obj)
  expect_lte(obj, 800)
})

test_that("sparse, nmf, hals", {
  res <- nmf(mat.sp, k, niter = 100, algo = "hals")
  obj <- Matrix::norm(mat.sp - res$W %*% t(res$H), type = "F")^2
  cat(obj)
  expect_lte(obj, 800)
})

test_that("sparse, nmf, mu", {
  res <- nmf(mat.sp, k, niter = 100, algo = "mu")
  obj <- Matrix::norm(mat.sp - res$W %*% t(res$H), type = "F")^2
  cat(obj)
  expect_lte(obj, 950)
})

symmat.sp <- Matrix::t(mat.sp) %*% mat.sp

test_that("sparse, symNMF, anlsbpp", {
  res <- symNMF(symmat.sp, k, niter = 100, lambda = lambda, algo = "anlsbpp")
  obj <- Matrix::norm(symmat.sp - res$W %*% t(res$H), type = "F")^2 + 
    lambda*Matrix::norm(res$W - res$H, type = "F")^2
  cat(obj)
  expect_lte(obj, 4e5)
})

test_that("sparse, symNMF, gnsym", {
  res <- symNMF(symmat.sp, k, niter = 100, algo = "gnsym")
  obj <- Matrix::norm(symmat.sp - res$H %*% t(res$H), type = "F")^2
  cat(obj)
  expect_lte(obj, 1e4)
})
