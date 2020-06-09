# Day 1
# Port to R of Alfred Galichon's Python code 

# rm(list = ls())
if (!require("pacman")) install.packages("pacman")
library('pacman')

#### 1. Generating Demand and Supply ####

# Pickup spots 

set.seed(777)
nbz = 3
z_df = as.data.frame(cbind(runif(nbz), runif(nbz)))
names(z_df) = c("h", "v")

# Drivers 

set.seed(778)
nbi = 200
i_df = as.data.frame(cbind(runif(nbi), runif(nbi), runif(nbi), 10 * runif(nbi)))
names(i_df) = c("h", "v", "tau", "lambda")

# Demand

set.seed(779)
nbj = 500
j_df = as.data.frame(cbind(runif(nbj), runif(nbj), 1/runif(nbj), 20 * runif(nbj), runif(nbj)))
names(j_df) = c("h", "v", "sigma", "epsilon", "eta")

# Generating demand

pacman::p_load('rdist')

avge_speed_drive = 25
avge_speed_walk = 4

T_iz = rdist::cdist(i_df[, c('h', 'v')], z_df[, c('h', 'v')]) / avge_speed_drive
T_jz = rdist::cdist(j_df[, c('h', 'v')], z_df[, c('h', 'v')]) / avge_speed_walk

# Aggregate supply

s_z <- function(p_z) {
  u_iz = t(matrix(p_z, nbz, nbi))  ^ (1 - i_df[, "tau"]) - i_df[, "lambda"] * T_iz
  u_i0 = rep(0, nbi)
  table(factor(apply(cbind(u_i0, u_iz), 1, which.max), levels = 2:(nbz+1), labels = 1:nbz))
}

# Aggregate demand

d_z <- function(p_z) {
  expos_j = 1 - 1/j_df[, 'sigma']
  c_jz = (t(matrix(p_z, nbz, nbj)) ^ expos_j + (j_df[, "epsilon"] * T_jz) ^ expos_j) ^ (1 / expos_j)
  c_j0 = j_df[, "eta"]
  table(factor(apply(cbind(c_j0, c_jz), 1, which.min), levels = 2:(nbz+1), labels = 1:nbz))
}

# Excess supply

e_z <- function(p_z) {
  s_z(p_z) - d_z(p_z)
}

#### 2. Smooth Max and a Trick to Improve Stability ####

# Smoothed Supply and Demand - very low stability of results.

ssmooth_z <- function(p_z, temperature = .01) {
  u_iz = t(matrix(p_z, nbz, nbi))  ^ (1 - i_df[, "tau"]) - i_df[, "lambda"] * T_iz
  u_i0 = rep(0, nbi)
  s_z = colSums(exp(cbind(u_i0, u_iz) / temperature) / rowSums(exp(cbind(u_i0, u_iz) / temperature)))
  names(s_z) = 0:nbz
  s_z[2:(nbz+1)]
}

dsmooth_z <- function(p_z, temperature = .01) {
  expos_j = 1 - 1/j_df[, 'sigma']
  c_jz = (t(matrix(p_z, nbz, nbj)) ^ expos_j + (j_df[, "epsilon"] * T_jz) ^ expos_j) ^ (1 / expos_j)
  c_j0 = j_df[, "eta"]
  d_z = colSums(exp(- cbind(c_j0, c_jz) / temperature) / rowSums(- exp(cbind(c_j0, c_jz) / temperature)))
  names(d_z) = 0:nbz
  d_z[2:(nbz+1)]
}

esmooth_z <- function(p_z, temp = .01) {
  ssmooth_z(p_z, temp) - dsmooth_z(p_z, temp)
}

# Smoothed Supply and Demand with a trick (overwirtes functions above as preferrable to use)

ssmooth_z <- function(p_z, temperature = .01) {
  u_iz = t(matrix(p_z, nbz, nbi))  ^ (1 - i_df[, "tau"]) - i_df[, "lambda"] * T_iz
  u_i0 = rep(0, nbi)
  s_z = colSums(exp((cbind(u_i0, u_iz) - apply(cbind(u_i0, u_iz), 1, max)) / temperature - 
        log(rowSums(exp((cbind(u_i0, u_iz) - apply(cbind(u_i0, u_iz), 1, max)) / temperature)))))
  names(s_z) = 0:nbz
  s_z[2:(nbz+1)]
}

dsmooth_z <- function(p_z, temperature = .01) {
  expos_j = 1 - 1/j_df[, 'sigma']
  c_jz = (t(matrix(p_z, nbz, nbj)) ^ expos_j + (j_df[, "epsilon"] * T_jz) ^ expos_j) ^ (1 / expos_j)
  c_j0 = j_df[, "eta"]
  d_z = colSums(exp((- cbind(c_j0, c_jz) + apply(cbind(c_j0, c_jz), 1, min)) / temperature - 
        log(rowSums(exp((- cbind(c_j0, c_jz) + apply(cbind(c_j0, c_jz), 1, min)) / temperature)))))
  names(d_z) = 0:nbz
  d_z[2:(nbz+1)]
}

esmooth_z <- function(p_z, temp = .01) {
  ssmooth_z(p_z, temp) - dsmooth_z(p_z, temp)
}

# Though there are concepts similar to python's classes in R, an R user would probably not use them commonly. I'll skip.

#### 3. Computation of Equilibrium Prices: Gauss-Seidel and Jacobi Algorithms ####

replaceentry <- function(vec, ind, x) {
  assign(deparse(substitute(vec)), replace(vec, ind, x), envir = .GlobalEnv)
}

pacman::p_load('pracma')

pmin = 0
pmax = 1e08
temperature = .001

cu <- function(vec, ind) {
  newprice <- pracma::brent(fun = function(x) esmooth_z(replaceentry(vec, ind, x), temperature)[ind], 
                            a = pmin, b = pmax)
  vectoreturn <- replaceentry(vec, ind, newprice$root)
  vectoreturn
}

GaussSeidel <- function(vec) {
  pp_z = vec
  for (zz in 1:nbz) {
    pp_z = cu(pp_z, zz)
  }
  pp_z
}

# Jacobi (allows for parallel computation)

pacman::p_load(foreach)
pacman::p_load(doParallel)

cores = floor(detectCores() / 2)
cl <- makeCluster(cores)
registerDoParallel(cl) # cluster closed at the very end.

Jacobi <- function(vec) {
  pp_z = vec
  pp_z <- foreach(z=1:nbz, .combine = c, .packages = "cmna", .export = ls(.GlobalEnv)) %dopar% {
    a = cu(pp_z, z)
    a[z]
  }
  pp_z
}

# Comparison of Gauss Seidel and Jacobi
p_z <- runif(nbz)

p_zGS <- p_zJ <- p_z
for (j in 1:10) {
  p_zGS <- GaussSeidel(p_zGS)
}
esmooth_z(p_zGS, temperature)

for (j in 1:10) {
  p_zJ <- Jacobi(p_zJ)
}
esmooth_z(p_zJ, temperature)

# Jacobi has converged much less than Gauss Seidel in 10 iterations.
# However, Jacobi can be parallelised for large nbz, whereas Gauss Seidel is constrained to be single-thread.

# Obtaining price vectors with the two methods

p_z <- runif(nbz)
p_zGS <- p_zJ <- p_z

valtol = 1e-5
steptol = 1e-12
criterione = criterionp = 10
iGS <- 1
a <- microbenchmark::microbenchmark(
while (criterione > valtol | criterionp > steptol) {
  p_zGS_ <- GaussSeidel(p_zGS)
  criterione <- max(abs(esmooth_z(p_zGS, temperature) - esmooth_z(p_zGS_, temperature)))
  criterionp <- max(abs(p_zGS - p_zGS_))
  p_zGS <- p_zGS_
  iGS <- iGS + 1
}, times = 1)

cat("Gauss Seidel Iterations: ", iGS, "\n")
cat("Gauss Seidel Excess Supply: ", esmooth_z(p_zGS, temperature), "\n")
cat("Gauss Seidel Runtime: ", a$time / 1e9, " seconds.")

criterione = criterionp = 10
iJ <- 1
b <- microbenchmark::microbenchmark(
  while (criterione > valtol | criterionp > steptol) {
  p_zJ_ <- Jacobi(p_zJ)
  criterione <- max(abs(esmooth_z(p_zJ, temperature) - esmooth_z(p_zJ_, temperature)))
  criterionp = max(abs(p_zJ - p_zJ_))
  p_zJ <- p_zJ_
  iJ <- iJ + 1 
}, times = 1)

cat("Jacobi Iterations: ", iJ, "\n")
cat("Jacobi Excess Supply: ", esmooth_z(p_zJ, temperature), "\n")
cat("Jacobi Runtime: ", b$time / 1e9, " seconds.")
cat("Jacobi Nodes Used: ", cores)

stopCluster(cl)

#### 99. Benchmarking Corner ####

# Brent's method and speed comparisons: a comparison of cmna::bisection, pracma::brent, pracma::bisect.

pacman::p_load("cmna")

cu <- function(vec, ind) {
  newprice <- cmna::bisection(f = function(x) esmooth_z(replaceentry(vec, ind, x), temperature)[ind], 
                              a = pmin, b = pmax, tol = 1e-6, m = 1e4)
  vectoreturn <- replaceentry(vec, ind, newprice)
  vectoreturn
}

cubrent <- function(vec, ind) {
  newprice <- pracma::brent(fun = function(x) esmooth_z(replaceentry(vec, ind, x), temperature)[ind], 
                              a = pmin, b = pmax)
  vectoreturn <- replaceentry(vec, ind, newprice$root)
  vectoreturn
}

cubisect <- function(vec, ind) {
  newprice <- pracma::bisect(fun = function(x) esmooth_z(replaceentry(vec, ind, x), temperature)[ind], 
                              a = pmin, b = pmax)
  vectoreturn <- replaceentry(vec, ind, newprice$root)
  vectoreturn
}

microbenchmark::microbenchmark(cu(p_z, 3), cubisect(p_z, 3), cubrent(p_z, 3), times = 100)

# My results: 
#    method       min       lq      mean     median     uq      max     neval
#  cu (cmna)   111.56017 118.93451 131.7361 122.7135 129.4310 388.7871   100
#  cubisect    386.94171 404.23572 426.8851 409.9341 431.5274 606.1484   100
#  cubrent     90.54658  95.73503  106.2200 100.7467 106.9277 269.5056   100
# 
# hence the choice to go with pracma::brent.
# cmna bisection is much faster than pracma's; unfortunately cmna does not provide the brent algorithm.
