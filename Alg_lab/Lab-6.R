E0 <- function(x) {
  x^3 * exp(-3 * x)
}

a <- 0
b <- 2

rectangle_method <- function(f, a, b, m) {
  h <- (b - a) / m
  F <- numeric(m)
  sum <- 0
  for (i in 1:m) {
    x_center <- a + (i - 0.5) * h
    F[i] <- f(x_center)
    sum <- sum + F[i]
  }
  return(h * sum)
}

gauss_method <- function(f, a, b, m) {
  gauss_data <- list(
    `3` = list(t = c(0.112701665, 0.5, 1 - 0.112701665),
               A = c(0.277777778, 0.444444444, 0.277777778)),
    `5` = list(t = c(0.046910077, 0.230765345, 0.5, 1 - 0.230765345, 1 - 0.046910077),
               A = c(0.118463443, 0.239314335, 0.284444444, 0.239314335, 0.118463443)),
    `7` = list(t = c(0.025446044, 0.129234407, 0.297077424, 0.5,
                     1 - 0.297077424, 1 - 0.129234407, 1 - 0.025446044),
               A = c(0.064742483, 0.139852696, 0.190915025, 0.208979592,
                     0.190915025, 0.139852696, 0.064742483)),
    `9` = list(t = c(0.015919880, 0.081934446, 0.193314284, 0.337873288, 0.5,
                     1 - 0.337873288, 1 - 0.193314284, 1 - 0.081934446, 1 - 0.015919880),
               A = c(0.040637194, 0.090324080, 0.130305348, 0.156173539, 0.165119678,
                     0.156173539, 0.130305348, 0.090324080, 0.040637194)),
    `11` = list(t = c(0.010885671, 0.056468700, 0.134923997, 0.240451935, 0.365228422, 0.5,
                      1 - 0.365228422, 1 - 0.240451935, 1 - 0.134923997, 1 - 0.056468700, 1 - 0.010885671),
                A = c(0.027834284, 0.062790185, 0.093145105, 0.116596882, 0.131402272, 0.136462543,
                      0.131402272, 0.116596882, 0.093145105, 0.062790185, 0.027834284))
  )
  
  T <- gauss_data[[as.character(m)]]$t
  A <- gauss_data[[as.character(m)]]$A
  
  sum <- 0
  for (i in 1:m) {
    x <- a + (b - a) * T[i]
    sum <- sum + A[i] * f(x)
  }
  return((b - a) * sum)
}

m_values <- seq(5, 100, by = 5)
rect_results <- sapply(m_values, function(m) rectangle_method(E0, a, b, m))

gauss_nodes <- c(3, 5, 7, 9, 11)
gauss_results <- sapply(gauss_nodes, function(m) gauss_method(E0, a, b, m))

true_value <- integrate(E0, a, b)$value

plot(m_values, rect_results, type = "l", col = "blue", lwd = 2,
     xlab = "Число узлов m", ylab = "Значение интеграла",
     main = "Метод прямоугольников и Гаусса")
abline(h = true_value, col = "black", lty = 2)
points(gauss_nodes, gauss_results, col = "red", pch = 19)
legend("bottomright", legend = c("Прямоугольники", "Гаусс", "Точное значение"),
       col = c("blue", "red", "black"), lty = c(1, NA, 2), pch = c(NA, 19, NA))
