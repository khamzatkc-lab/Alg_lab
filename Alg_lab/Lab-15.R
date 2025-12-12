l <- 1       
Tmax <- 1    
N1 <- 50      
N2 <- 100     

epsilon <- 1
mu <- 1
sigma <- 0.1
lambda <- 1

h <- 2*l / N1
tau <- Tmax / N2
z <- seq(-l, l, length.out = N1+1)
t <- seq(0, Tmax, length.out = N2+1)

y <- matrix(0, nrow = N1+1, ncol = N2+1)

F0 <- function(z, tt) {
  return(sin(pi*z)*sin(tt))
}

for (j in 2:N2) {
  for (i in 2:N1) {
    y[i, j+1] <- ( (epsilon*mu*(y[i,j] - 2*y[i,j-1] + y[i,j-2]) / tau^2) +
                     (mu*sigma*(y[i,j] - y[i,j-2]) / (2*tau)) +
                     ((y[i+1,j] - 2*y[i,j] + y[i-1,j]) / h^2) -
                     (lambda^2 * y[i,j]) +
                     F0(z[i], t[j]) )
  }
  y[1, j+1] <- 0
  y[N1+1, j+1] <- 0
}

image(t, z, y, main="Явная схема для уравнения геоэлектрики",
      xlab="t", ylab="z", col=heat.colors(50))
