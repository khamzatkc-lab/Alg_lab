c <- c(1.5, 1.0, 0.5)  
a <- -c                
b <- c                 

alpha <- c(0, 0, 0)   
beta <- c(2, 2, 2)    

N <- 100000

X <- matrix(runif(3 * N, min = a[1], max = b[1]), ncol = 3)
for (i in 2:3) {
  X[, i] <- runif(N, min = a[i], max = b[i])
}

inside <- rowSums((X / c)^2) <= 1

rho <- rowSums(abs(X[inside, ] - matrix(alpha, nrow = sum(inside), ncol = 3, byrow = TRUE))^matrix(beta, nrow = sum(inside), ncol = 3, byrow = TRUE))

W_volume <- prod(b - a)
V_volume <- W_volume * sum(inside) / N
Q <- mean(rho) * V_volume

cat("Приближённое значение интеграла Q =", Q, "\n")
