a <- 0
b <- 1
n <- 30
h <- 1 / n
x <- seq(a, b, by = h)

q <- 3 / 5
set.seed(42)  
noise <- runif(length(x), -0.01 * x^q, 0.01 * x^q)
y <- x * (1 - x) + noise

phi <- cbind(1, x, x^2)

A <- solve(t(phi) %*% phi) %*% t(phi) %*% y
y_approx <- phi %*% A

error <- abs(y - y_approx)
epsilon_max <- max(error)
epsilon_mean <- mean(error)

plot(x, y, type = "p", col = "blue", pch = 16, main = "Аппроксимация методом наименьших квадратов",
     xlab = "x", ylab = "y", ylim = range(c(y, y_approx)))
lines(x, y_approx, col = "red", lwd = 2)
legend("topright", legend = c("Исходные данные", "Аппроксимация"), col = c("blue", "red"), pch = c(16, NA), lty = c(NA, 1))

cat("Максимальная погрешность:", epsilon_max, "\n")
cat("Средняя погрешность:", epsilon_mean, "\n")

###

y_test <- 1 + x^2
A_test <- solve(t(phi) %*% phi) %*% t(phi) %*% y_test
y_test_approx <- phi %*% A_test
error_test <- abs(y_test - y_test_approx)

cat("Погрешность тестовой функции (макс):", max(error_test), "\n")
cat("Погрешность тестовой функции (средн):", mean(error_test), "\n")

