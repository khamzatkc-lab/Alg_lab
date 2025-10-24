f <- function(x) x^2 * sinh(x)
f_prime <- function(x) 2*x*sinh(x) + x^2*cosh(x)
f_double_prime <- function(x) 2*sinh(x) + 4*x*cosh(x) + x^2*sinh(x)

a <- 0
b <- 2
n <- 50
h <- (b - a) / n
x <- seq(a, b, length.out = n + 1)

y <- f(x)

f1_num <- numeric(length(x))
f1_num[2:n] <- (y[3:(n+1)] - y[1:(n-1)]) / (2*h)
f1_num[1] <- (y[2] - y[1]) / h
f1_num[n+1] <- (y[n+1] - y[n]) / h

f2_num <- numeric(length(x))
f2_num[2:n] <- (y[3:(n+1)] - 2*y[2:n] + y[1:(n-1)]) / (h^2)
f2_num[1] <- (y[3] - 2*y[2] + y[1]) / (h^2)
f2_num[n+1] <- (y[n+1] - 2*y[n] + y[n-1]) / (h^2)

f1_true <- f_prime(x)
f2_true <- f_double_prime(x)

error_f1 <- abs(f1_true - f1_num)
error_f2 <- abs(f2_true - f2_num)

par(mfrow = c(2, 2))

plot(x, f1_true, type = "l", col = "blue", lwd = 2, main = "Первая производная")
lines(x, f1_num, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Аналитическая", "Численная"), col = c("blue", "red"), lty = c(1,2))

plot(x, error_f1, type = "l", col = "darkgreen", lwd = 2, main = "Погрешность первой производной")

plot(x, f2_true, type = "l", col = "blue", lwd = 2, main = "Вторая производная")
lines(x, f2_num, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Аналитическая", "Численная"), col = c("blue", "red"), lty = c(1,2))

plot(x, error_f2, type = "l", col = "darkgreen", lwd = 2, main = "Погрешность второй производной")

###

f <- function(x) x^2 / 2
f_prime <- function(x) x
f_double_prime <- function(x) rep(1, length(x))

