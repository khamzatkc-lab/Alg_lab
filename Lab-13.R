make_uniform_grid <- function(a = 0, b = 1, N = 100) {
  x <- seq(a, b, length.out = N + 1)
  list(x = x, N = N, a = a, b = b)
}

make_nonuniform_grid <- function(a = 0, b = 1, N = 100, stretch = 2.0) {
  s <- seq(0, 1, length.out = N + 1)
  phi <- function(t) { (t^stretch) / (t^stretch + (1 - t)^stretch) }
  x <- a + (b - a) * phi(s)
  list(x = x, N = N, a = a, b = b)
}

assemble_variational_1d <- function(x, f_fun, alpha, beta) {
  n <- length(x)
  N <- n - 1
  h <- diff(x) 
  m <- N - 1
  A <- matrix(0, nrow = m, ncol = m)
  b <- numeric(m)
  
  for (j in 1:m) {
    hjm <- h[j]               
    hjp <- h[j + 1]   
    A[j, j] <- (1 / hjm) + (1 / hjp)
    if (j - 1 >= 1) {
      A[j, j - 1] <- A[j, j - 1] - 1 / hjm
    }
    if (j + 1 <= m) {
      A[j, j + 1] <- A[j, j + 1] - 1 / hjp
    }
  }
  
  f_vals <- f_fun(x)
  for (j in 1:m) {
    g <- j + 1
    hl <- h[g - 1]
    b[j] <- b[j] + hl * (f_vals[g - 1] + f_vals[g]) / 4
    hr <- h[g]
    b[j] <- b[j] + hr * (f_vals[g] + f_vals[g + 1]) / 4
  }
  
  b[1] <- b[1] + alpha / h[1]
  b[m] <- b[m] + beta / h[N]
  
  list(A = A, b = b)
}

add_delta_source <- function(x, rhs, xs, q) {
  n <- length(x)
  if (xs <= x[1]) {
    return(rhs)
  }
  if (xs >= x[n]) {
    return(rhs)
  }
  i <- max(which(x <= xs))
  if (i == n) i <- n - 1
  if (abs(xs - x[i]) < .Machine$double.eps) {
    j <- i - 1
    if (j >= 1 && j <= length(rhs)) {
      rhs[j] <- rhs[j] - q  
    }
    return(rhs)
  }
  if (abs(xs - x[i + 1]) < .Machine$double.eps) {
    j <- i
    if (j >= 1 && j <= length(rhs)) {
      rhs[j] <- rhs[j] - q
    }
    return(rhs)
  }
  
  hseg <- x[i + 1] - x[i]
  w_left <- (x[i + 1] - xs) / hseg
  w_right <- (xs - x[i]) / hseg
  j_left <- i - 1
  j_right <- i
  
  if (j_left >= 1 && j_left <= length(rhs)) rhs[j_left] <- rhs[j_left] - q * w_left
  if (j_right >= 1 && j_right <= length(rhs)) rhs[j_right] <- rhs[j_right] - q * w_right
  
  rhs
}

assemble_summation_identities <- function(x, f_fun, alpha, beta) {
  n <- length(x)
  N <- n - 1
  m <- N - 1
  h <- diff(x)
  A <- matrix(0, nrow = m, ncol = m)
  b <- numeric(m)
  
  for (j in 1:m) {
    hjm <- h[j]
    hjp <- h[j + 1]
    A[j, j] <- (1 / hjm) + (1 / hjp)
    if (j - 1 >= 1) A[j, j - 1] <- A[j, j - 1] - 1 / hjm
    if (j + 1 <= m) A[j, j + 1] <- A[j, j + 1] - 1 / hjp
  }
  
  f_vals <- f_fun(x)
  for (j in 1:m) {
    g <- j + 1
    hl <- h[g - 1]
    hr <- h[g]
    b[j] <- b[j] + hl * (f_vals[g - 1] + f_vals[g]) / 4
    b[j] <- b[j] + hr * (f_vals[g] + f_vals[g + 1]) / 4
  }
  
  b[1] <- b[1] + alpha / h[1]
  b[m] <- b[m] + beta / h[N]
  
  list(A = A, b = b)
}

solve_dirichlet_1d <- function(A, b, x, alpha, beta) {
  u_inner <- solve(A, b)
  u <- numeric(length(x))
  u[1] <- alpha
  u[length(x)] <- beta
  u[2:(length(x) - 1)] <- u_inner
  u
}

solve_bvp_variational <- function(grid, f_fun, alpha = 0, beta = 0, xs = NA, q = 0) {
  asm <- assemble_variational_1d(grid$x, f_fun, alpha, beta)
  if (!is.na(xs) && q != 0) {
    asm$b <- add_delta_source(grid$x, asm$b, xs, q)
  }
  u <- solve_dirichlet_1d(asm$A, asm$b, grid$x, alpha, beta)
  list(x = grid$x, u = u)
}

solve_bvp_summation <- function(grid, f_fun, alpha = 0, beta = 0, xs = NA, q = 0) {
  asm <- assemble_summation_identities(grid$x, f_fun, alpha, beta)
  if (!is.na(xs) && q != 0) {
    asm$b <- add_delta_source(grid$x, asm$b, xs, q)
  }
  u <- solve_dirichlet_1d(asm$A, asm$b, grid$x, alpha, beta)
  list(x = grid$x, u = u)
}

grid_uni <- make_uniform_grid(a = 0, b = 1, N = 100)
f1 <- function(x) rep(1, length(x))

res_var_1 <- solve_bvp_variational(grid_uni, f1, alpha = 0, beta = 0)
res_sum_1 <- solve_bvp_summation(grid_uni, f1, alpha = 0, beta = 0)

xs <- 0.5
q <- 1.0
res_var_delta <- solve_bvp_variational(grid_uni, f1, alpha = 0, beta = 0, xs = xs, q = q)
res_sum_delta <- solve_bvp_summation(grid_uni, f1, alpha = 0, beta = 0, xs = xs, q = q)

grid_nonuni <- make_nonuniform_grid(a = 0, b = 1, N = 120, stretch = 2.2)
f2 <- function(x) sin(pi * x)
res_var_nonuni <- solve_bvp_variational(grid_nonuni, f2, alpha = 0, beta = 0, xs = 0.53, q = 1.0)
res_sum_nonuni <- solve_bvp_summation(grid_nonuni, f2, alpha = 0, beta = 0, xs = 0.53, q = 1.0)

plot(res_var_1$x, res_var_1$u, type = "l", lwd = 2, col = "steelblue",
     main = "Решения: вариационный и сумматорные тождества",
     xlab = "x", ylab = "u(x)")
lines(res_sum_1$x, res_sum_1$u, lwd = 2, col = "firebrick")
legend("topright", legend = c("Variational", "Summation identities"),
       col = c("steelblue", "firebrick"), lwd = 2, bty = "n")

plot(res_var_delta$x, res_var_delta$u, type = "l", lwd = 2, col = "purple",
     main = "Дельта-источник в xs=0.5, q=1.0",
     xlab = "x", ylab = "u(x)")
lines(res_sum_delta$x, res_sum_delta$u, lwd = 2, col = "orange")
legend("topright", legend = c("Variational+delta", "Summation+delta"),
       col = c("purple", "orange"), lwd = 2, bty = "n")



