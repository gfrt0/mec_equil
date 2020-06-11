# Day 2
# Port to R of Alfred Galichon's Python code 

# rm(list = ls())
if (!require("pacman")) install.packages("pacman")
library('pacman')

set.seed(777)
d = 8
nbx = 50
nby = 30

rg = .8
rs = .6
r = .7

n_x = rep(1, nbx)
m_y = rep(1, nby)
xi_x_k = matrix(runif(nbx * d), nrow = nbx, ncol = d)
zeta_y_k = matrix(runif(nby * d), nrow = nby, ncol = d)

gamma_x_y = alpha_x_y = matrix(0, nrow = nbx, ncol = nby)

c1 = expand.grid(xi_x_k[ , 7], zeta_y_k[ , 7])
c2 = expand.grid(xi_x_k[ , 8], zeta_y_k[ , 8])
alpha_x_y = - sqrt((c1[, 1] - c1[, 2])^2 + (c2[, 1] - c2[, 2])^2) / 50
alpha_x_y = matrix(alpha_x_y, nrow = nbx, ncol = nby)

d1 = do.call(cbind, lapply(1:6, FUN = function(x) expand.grid(xi_x_k[, x], zeta_y_k[, x])))
d2 = rowSums(do.call(cbind, lapply(c(1, 3, 5), function(x) (d1[, x] * d1[, x+1])^rg)))^(r / rg)
d3 = rowSums(do.call(cbind, lapply(c(7, 9, 11), function(x) (d1[, x] * d1[, x+1])^rs)))^(r / rs)
gamma_x_y = matrix((d2 + d3)^(1/r), nrow = nbx, ncol = nby)

Phi_x_y = alpha_x_y + gamma_x_y
K_x_y = temp * exp((Phi_x_y)/(2*temp))

# Choo Siow model as an optimisation problem

temp = 1

fxy_ChooSiow <- function(ab) {
  a <- ab[1:nbx]
  b <- ab[nbx+1:nby]

  criterion = sum(n_x * a) + sum(m_y * b)  + 
              2 * temp * sum(sum(exp((Phi_x_y - a - b)/(2*temp)))) + 
              temp * sum(exp(- a / temp)) +  temp * sum(exp(- b / temp))
  criterion
}

gradients <- function(ab) {
  a <- ab[1:nbx]
  b <- ab[nbx+1:nby]

  c(n_x - exp(- a / temp) - rowSums(exp((Phi_x_y - a - b)/(2*temp))),
    m_y - exp(- b / temp) - colSums(exp((Phi_x_y - a - b)/(2*temp))))
}

result <- stats::optim(par = rep(0, nbx + nby), fn = fxy_ChooSiow, method = "BFGS", gr = gradients)

a <- result$par[1:nbx]
b <- result$par[nbx+1:nby]

mu_x_0 = exp(- a / temp)
mu_x_0
mu_0_y = exp(- b / temp)
mu_0_y

# fundamental differences in the estimates by Alfred and mine.

mu_x_y = exp((Phi_x_y - a - b)/(2*temp))

A = sum(sum(alpha_x_y * mu_x_y))
B = sum(sum(gamma_x_y * mu_x_y))

# Choo Siow as an Equilibrium Problem with Gross Substitutes

# In this case we work with Gauss Seidel iteration in closed form

IPFP <- function(Phi_x_y, tol1 = 1e-12, tol2 = 1e-5) {

  K_x_y = exp((Phi_x_y)/(2*temp))
  
  A_x = rep(1, nbx)
  B_y = rep(1, nby)
  
  i = 1
  
  converge1 = converge2 = 10
  
  time <- microbenchmark::microbenchmark(
  while (converge1 > tol1 || converge2 > tol2) {
    A_x_ = sqrt(n_x + ((K_x_y %*% B_y)/2)^2) - (K_x_y %*% B_y)/2
    B_y_ = sqrt(m_y + ((t(K_x_y) %*% A_x)/2)^2) - (t(K_x_y) %*% A_x)/2
    converge1 = max(abs(c(A_x_ - A_x, B_y_ - B_y)))
    A_x = A_x_
    B_y = B_y_
    converge2 = max(c(abs(n_x - A_x^2 - (K_x_y %*% B_y) * A_x),
                      abs(m_y - B_y^2 - t(K_x_y) %*% A_x * B_y)))
    i = i + 1
  }, times = 1)
  
  cat("IPFP converged in", i, "iterations,", time$time/ 1e9, "seconds. \n")
  cat("Discrepancy for A_x:", converge1, "\n")
  cat("Discrepancy for B_y:", converge2, "\n")
  
  return(list('mu_x0' = as.vector(A_x^2), 'mu_0y' = as.vector(B_y^2), 'convertols' = c(converge1, converge2)))
}

IPFP(Phi_x_y)

# Linear Taxes

tau_k = seq(0, .5, .05)

PHI_k = ALPHA_k = GAMMA_k = rep(0, length(tau_k))

for (k in 1:length(tau_k)) {
  phitau_x_y <- alpha_x_y + (1 - tau_k[k]) * gamma_x_y
  AB <- IPFP(phitau_x_y)
  ax <- - temp * log(AB$mu_x0)
  by <- - temp * log(AB$mu_0y)
  mutau_x_y <- exp((phitau_x_y - ax - by) / 2*temp)
  ALPHA_k[k] <- sum(sum(mutau_x_y * alpha_x_y))
  GAMMA_k[k] <- sum(sum(mutau_x_y * gamma_x_y))
  PHI_k[k]   <- sum(sum(mutau_x_y * (gamma_x_y + alpha_x_y)))
}

plot(tau_k, ALPHA_k, ylab = "Alpha", xlab = "tau")
plot(tau_k, GAMMA_k, ylab = "Gamma", xlab = "tau")
plot(tau_k, PHI_k, ylab = "Phi", xlab = "tau")
