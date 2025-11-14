slau_solver <- function(a, b, c, d) {
  n <- length(b)
  
  c_star <- numeric(n-1)
  d_star <- numeric(n)
  
  c_star[1] <- c[1] / b[1]
  d_star[1] <- d[1] / b[1]
  
  for (i in 2:(n-1)) {
    denom <- b[i] - a[i-1] * c_star[i-1]
    c_star[i] <- c[i] / denom
    d_star[i] <- (d[i] - a[i-1] * d_star[i-1]) / denom
  }
  
  denom <- b[n] - a[n-1] * c_star[n-1]
  d_star[n] <- (d[n] - a[n-1] * d_star[n-1]) / denom
  
  x <- numeric(n)
  x[n] <- d_star[n]
  
  for (i in (n-1):1) {
    x[i] <- d_star[i] - c_star[i] * x[i+1]
  }
  
  return(x)
}

a <- c(1, 1)          
b <- c(2, 3, 2)      
c <- c(1, 1)         
d <- c(5, 10, 6)    

solution <- slau_solver(a, b, c, d)
print(solution)


