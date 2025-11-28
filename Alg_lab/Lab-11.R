a <- 1             
l <- 1             
T <- 1             
N <- 100           
M <- 100           
h <- l / N      
tau <- T / M       

if (a * tau / h > 1) {
  stop("Нарушено условие Куранта: a * tau / h должно быть ≤ 1")
}

x <- seq(0, l, length.out = N + 1)
t <- seq(0, T, length.out = M + 1)


mu1 <- function(x) sin(pi * x)      
mu2 <- function(x) 0                   
mu3 <- function(t) 0                
mu4 <- function(t) 0             
f <- function(x, t) 0                  

u <- matrix(0, nrow = N + 1, ncol = M + 1)

u[, 1] <- mu1(x)
for (i in 2:N) {
  u[i, 2] <- mu1(x[i]) + tau * mu2(x[i]) +
    (tau^2 / 2) * (a^2 * (u[i + 1, 1] - 2 * u[i, 1] + u[i - 1, 1]) / h^2 + f(x[i], t[1]))
}
u[1, 2] <- mu3(t[2])
u[N + 1, 2] <- mu4(t[2])

for (n in 2:(M - 1)) {
  for (i in 2:N) {
    u[i, n + 1] <- 2 * u[i, n] - u[i, n - 1] +
      (a^2 * tau^2 / h^2) * (u[i + 1, n] - 2 * u[i, n] + u[i - 1, n]) +
      tau^2 * f(x[i], t[n])
  }
  u[1, n + 1] <- mu3(t[n + 1])
  u[N + 1, n + 1] <- mu4(t[n + 1])
}

image(x, t, t(u), col = heat.colors(100), xlab = "x", ylab = "t", main = "Решение волнового уравнения")
