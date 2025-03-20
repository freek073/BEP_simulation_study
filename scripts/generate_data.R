# This script will be used to generate the data for the simulation study 
# it will be stored in the 'data' folder for clearer overview

# specific seed which will allow for reproducibility
set.seed(45468)
library(MASS)

n <- c(20, 50, 200, 600, 1000)
i <- 1  
set.seed(45468)

p1g1 <- c(rep(sqrt(.6),10), rep(.4,4), rep(0,6))
p2g1 <- c(rep(.4,4), rep(0,6), rep(sqrt(.6),10))
P1 <- cbind(p1g1, p2g1)

p1g2 <- c(rep(sqrt(.6),10), rep(0,4), rep(.4,4), 0, 0)
p2g2 <- c(rep(0,4), rep(.4,4), 0, 0, rep(sqrt(.6),10))
P2 <- cbind(p1g2, p2g2)

PSI1 <- diag(1 - rowSums(P1^2))
SIGMA1 <- P1 %*% t(P1) + PSI1

PSI2 <- diag(1 - rowSums(P2^2))
SIGMA2 <- P2 %*% t(P2) + PSI2

# generate data
data_2factors <- mvrnorm(n = n[i], mu = rep(0,20), Sigma = SIGMA1, empirical = TRUE)

# factor complexity 5
p3g1 <- c(rep(0,5), rep(.4,5), rep(sqrt(.6),5), rep(0,5))
p4g1 <- c(rep(0,10), rep(.4,5), rep(sqrt(.6),5))
p5g1 <- c(rep(0,15), rep(.4,5))
P5 <- cbind(p1g1, p2g1, p3g1, p4g1, p5g1)

PSI5 <- diag(1 - rowSums(P5^2))
SIGMA5 <- P5 %*% t(P5) + PSI5

# generate data
data_5factors <- mvrnorm(n = n[i], mu = rep(0,20), Sigma = SIGMA5, empirical = TRUE)

# save the data set into the data folder
saveRDS(data_2factors, file = file.path("data", paste0("simulated_data_n", n[i], "_2factors.RDS")))
saveRDS(data_5factors, file = file.path("data", paste0("simulated_data_n", n[i], "_5factors.RDS")))

