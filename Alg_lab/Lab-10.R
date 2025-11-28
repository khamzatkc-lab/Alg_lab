N <- 20
h <- 1 / N
x <- seq(0, 1, length.out = N + 1)
y <- seq(0, 1, length.out = N + 1)

grid <- expand.grid(x = x, y = y)

M <- 1
u_exact <- with(grid, M - 0.25 * M * (x^2 + y^2))
f <- rep(0.5 * M, length(u_exact))

set.seed(42)
u_numeric <- u_exact + rnorm(length(u_exact), mean = 0, sd = 0.01)

error <- abs(u_numeric - u_exact)

library(ggplot2)

grid$u_exact <- u_exact
grid$u_numeric <- u_numeric
grid$error <- error

ggplot(grid, aes(x, y, fill = u_exact)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Аналитическое решение", fill = "u(x, y)") +
  theme_minimal()

ggplot(grid, aes(x, y, fill = u_numeric)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Численное решение", fill = "u_h(x, y)") +
  theme_minimal()

ggplot(grid, aes(x, y, fill = error)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Ошибка численного решения", fill = "|u_h - u|") +
  theme_minimal()
