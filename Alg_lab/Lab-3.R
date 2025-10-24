a <- 2
b <- 3
k <- 2
m <- 2
n <- 50

x <- seq(0, 1, length.out = n)
y <- (a + b * x^m)^(-1/k)

y_test <- x^2

function1 <- function(y) {
  n <- length(y)
  diff_table <- matrix(0, n, n)
  diff_table[,1] <- y
  for (j in 2:n) {
    for (i in 1:(n-j+1)) {
      diff_table[i,j] <- diff_table[i+1,j-1] - diff_table[i,j-1]
    }
  }
  return(diff_table)
}

diff_y <- function1(y)
diff_test <- function1(y_test)

interpolation <- function(x_nodes, diff_table, x_eval) {
  n <- length(x_nodes)
  h <- x_nodes[2] - x_nodes[1]
  s <- (x_eval - x_nodes[1]) / h
  result <- diff_table[1,1]
  mult_term <- 1
  for (i in 2:n) {
    mult_term <- mult_term * (s - (i - 2)) / (i - 1)
    result <- result + mult_term * diff_table[1,i]
  }
  return(result)
}

x_eval <- seq(0, 1, length.out = 200)
y_interp <- sapply(x_eval, function(xi) interpolation(x, diff_y, xi))
y_true <- (a + b * x_eval^m)^(-1/k)

errors <- abs(y_true - y_interp)
epsilon_max <- max(errors)
epsilon_mean <- mean(errors)

y_test_interp <- sapply(x_eval, function(xi) interpolation(x, diff_test, xi))
test_errors <- abs(x_eval^2 - y_test_interp)

library(ggplot2)
df <- data.frame(x = x_eval, y_true = y_true, y_interp = y_interp, error = errors)

ggplot(df, aes(x)) +
  geom_line(aes(y = y_true), color = "black", linewidth = 1.2, linetype = "solid") +
  geom_line(aes(y = y_interp), color = "green", linewidth = 1.2, linetype = "dashed") +
  geom_line(aes(y = error), color = "red", linewidth = 1, linetype = "solid") +
  labs(title = "Интерполяция Ньютона",
       subtitle = paste("εmax =", round(epsilon_max, 6), "εmean =", round(epsilon_mean, 6)),
       y = "y(x)", x = "x") +
  theme_minimal()
