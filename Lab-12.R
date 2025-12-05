poisson_solver <- function(N = 20, M = 20, 
                           lx = 1, ly = 1, 
                           f_func = function(x, y) { 1 }, 
                           boundary_func = function(x, y) { 0 }, 
                           method = "SOR", 
                           omega = 1.5, 
                           tol = 1e-6, 
                           max_iter = 10000) {
  
  hx <- lx / N
  hy <- ly / M
  
  x <- seq(0, lx, length.out = N + 1)
  y <- seq(0, ly, length.out = M + 1)
  
  u <- matrix(0, nrow = N + 1, ncol = M + 1)
  
  for (j in 1:(N+1)) {
    u[j, 1] <- boundary_func(x[j], y[1])
    u[j, M+1] <- boundary_func(x[j], y[M+1])
  }
  for (i in 1:(M+1)) {
    u[1, i] <- boundary_func(x[1], y[i])
    u[N+1, i] <- boundary_func(x[N+1], y[i])
  }
  
  iter <- 0
  error <- tol + 1
  
  while (iter < max_iter && error > tol) {
    error <- 0
    for (j in 2:N) {
      for (i in 2:M) {
        rhs <- f_func(x[j], y[i])
        u_old <- u[j, i]
        
        u_new <- ((u[j-1,i] + u[j+1,i]) / hx^2 +
                    (u[j,i-1] + u[j,i+1]) / hy^2 -
                    rhs) / (2/hx^2 + 2/hy^2)
        
        if (method == "SOR") {
          u[j,i] <- u_old + omega * (u_new - u_old)
        } else {
          u[j,i] <- u_new
        }
        
        error <- max(error, abs(u[j,i] - u_old))
      }
    }
    iter <- iter + 1
  }
  
  cat("Метод:", method, "Итераций:", iter, "Погрешность:", error, "\n")
  return(list(u = u, x = x, y = y))
}

res <- poisson_solver(N = 30, M = 30, 
                      f_func = function(x,y) { sin(pi*x)*sin(pi*y) }, 
                      boundary_func = function(x,y) { 0 }, 
                      method = "SOR", omega = 1.8)
  
library(plotly)
plot_ly(x = res$x, y = res$y, z = res$u) %>% 
  add_surface() %>% 
  layout(title = "Решение уравнения Пуассона методом SOR")









