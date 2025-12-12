eps <- 1.0     
mu  <- 1.0     
lambda <- 1.0  
Tmax <- 1.0    
L <- 1.0       
dt <- 0.01     
dz <- 0.01     

time <- seq(0, Tmax, by = dt)
space <- seq(0, L, by = dz)
f <- function(t) sin(pi * t)   
p <- rep(1, length(space))

solve_forward <- function(p) {
  u <- outer(space, time, function(z,t) sin(pi*z)*sin(pi*t))
  return(u)
}

solve_adjoint <- function(p, u) {
  phi <- outer(space, time, function(z,t) cos(pi*z)*sin(pi*t))
  return(phi)
}

compute_gradient <- function(u, phi) {
  grad <- colSums(u * phi) * dt
  return(grad)
}

functional <- function(u, f) {
  J <- sum((u[1,] - f(time))^2) * dt
  return(J)
}

gradient_descent <- function(p0, alpha = 0.01, max_iter = 50) {
  p <- p0
  for (n in 1:max_iter) {
    u <- solve_forward(p)
    phi <- solve_adjoint(p, u)
    grad <- compute_gradient(u, phi)
    p <- p - alpha * grad
    J <- functional(u, f)
    cat("Iter:", n, "Functional:", J, "\n")
  }
  return(p)
}

p_solution <- gradient_descent(p)
