build_system <- function(n) {
  A <- matrix(0, n, n)
  b <- numeric(n)
  for (i in 1:n) {
    A[i,i] <- 11 * sqrt(i)
    for (j in 1:n) if (j != i) {
      sign_ij <- ifelse(((i + j) %% 2)==0, 1, -1)
      A[i,j] <- sign_ij * 1e-2 * sin(pi/2 / abs(i - j))
    }
    b[i] <- 8.5^(i/3)
  }
  list(A = A, b = b)
}

lu_doolittle <- function(A) {
  n <- nrow(A)
  L <- diag(1, n)
  U <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      s <- 0
      if (i > 1) s <- sum(L[i,1:(i-1)] * U[1:(i-1), j])
      U[i,j] <- A[i,j] - s
    }
    if (U[i,i] == 0) stop("Error")
    for (j in (i+1):n) {
      s <- 0
      if (i > 1) s <- sum(L[j,1:(i-1)] * U[1:(i-1), i])
      L[j,i] <- (A[j,i] - s) / U[i,i]
    }
  }
  Acomp <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (j < i) Acomp[i,j] <- L[i,j]
      else Acomp[i,j] <- U[i,j]
    }
  }
  list(L = L, U = U, compact = Acomp)
}

lu_partial_pivot <- function(A) {
  n <- nrow(A)
  P <- diag(1, n)
  L <- matrix(0, n, n)
  U <- A
  for (k in 1:(n-1)) {
    piv_row <- which.max(abs(U[k:n, k])) + k - 1
    if (abs(U[piv_row, k]) < .Machine$double.eps) stop("Error")
    if (piv_row != k) {
      U[c(k, piv_row), ] <- U[c(piv_row, k), ]
      P[c(k, piv_row), ] <- P[c(piv_row, k), ]
      if (k > 1) L[c(k, piv_row), 1:(k-1)] <- L[c(piv_row, k), 1:(k-1)]
    }
    for (i in (k+1):n) {
      L[i,k] <- U[i,k] / U[k,k]
      U[i, ] <- U[i, ] - L[i,k] * U[k, ]
    }
  }
  diag(L) <- 1
  list(L = L, U = U, P = P)
}

forward_subst <- function(L, b) {
  n <- length(b)
  y <- numeric(n)
  for (i in 1:n) {
    y[i] <- b[i] - sum(L[i,1:(i-1)] * y[1:(i-1)])
  }
  y
}

backward_subst <- function(U, y) {
  n <- length(y)
  x <- numeric(n)
  for (i in n:1) {
    x[i] <- (y[i] - sum(U[i,(i+1):n] * x[(i+1):n])) / U[i,i]
  }
  x
}

solve_lu <- function(A, b) {
  res <- lu_doolittle(A)
  y <- forward_subst(res$L, b)
  x <- backward_subst(res$U, y)
  list(x = x, L = res$L, U = res$U, compact = res$compact)
}

solve_lu_pivot <- function(A, b) {
  res <- lu_partial_pivot(A)
  Pb <- res$P %*% b
  y <- forward_subst(res$L, Pb)
  x <- backward_subst(res$U, y)
  list(x = x, L = res$L, U = res$U, P = res$P)
}

residual_norm <- function(A, x, b) {
  max(abs(A %*% x - b))
}

example_run <- function(n = 5, use_pivot = TRUE) {
  sys <- build_system(n)
  A <- sys$A; b <- sys$b
  if (use_pivot) {
    sol <- solve_lu_pivot(A, b)
  } else {
    sol <- solve_lu(A, b)
  }
  delta <- residual_norm(A, sol$x, b)
  list(n = n, A = A, b = b, x = sol$x, delta = delta, L = sol$L, U = sol$U, P = sol$P)
}

set.seed(1)
out <- example_run(n = 5, use_pivot = TRUE)
cat("n =", out$n, "\n")
cat("delta =", format(out$delta, scientific = TRUE), "\n")
cat("Решение x:\n")
print(out$x)
