p <- rep(1, length(space))
q <- rep(1, length(space))

f1 <- function(t) sin(pi*t)
f2 <- function(t) cos(pi*t)

solve_forward_multi <- function(p, q) {
  u1 <- outer(space, time, function(z,t) sin(pi*z)*sin(pi*t))
  u2 <- outer(space, time, function(z,t) cos(pi*z)*sin(pi*t))
  return(list(u1=u1, u2=u2))
}

solve_adjoint_multi <- function(p, q, U) {
  phi1 <- outer(space, time, function(z,t) cos(pi*z)*sin(pi*t))
  phi2 <- outer(space, time, function(z,t) sin(pi*z)*cos(pi*t))
  return(list(phi1=phi1, phi2=phi2))
}

compute_gradient_multi <- function(U, Phi) {
  grad_p <- colSums(U$u1 * Phi$phi1) * dt + colSums(U$u1 * Phi$phi1) * dt
  grad_q <- colSums(U$u2 * Phi$phi2) * dt + colSums(U$u2 * Phi$phi2) * dt
  return(list(grad_p=grad_p, grad_q=grad_q))
}

functional_multi <- function(U, f1, f2) {
  J <- sum((U$u1[1,] - f1(time))^2) * dt +
    sum((U$u2[1,] - f2(time))^2) * dt
  return(J)
}

gradient_descent_multi <- function(p0, q0, alpha=0.01, max_iter=50) {
  p <- p0; q <- q0
  for (n in 1:max_iter) {
    U <- solve_forward_multi(p,q)
    Phi <- solve_adjoint_multi(p,q,U)
    grad <- compute_gradient_multi(U,Phi)
    p <- p - alpha * grad$grad_p
    q <- q - alpha * grad$grad_q
    J <- functional_multi(U,f1,f2)
    cat("Iter:", n, "Functional:", J, "\n")
  }
  return(list(p=p,q=q))
}

Q_solution <- gradient_descent_multi(p,q)
