# This script will be used to generate the data for the simulation study 
# it will be stored in the 'data' folder for clearer overview

# specific seed which will allow for reproducibility
library(MASS)
set.seed(45468)

# Define sample sizes to simulate
n_values <- c(20, 50, 200, 600, 1000)

# Define 2-factor model loading matrix
p1g1 <- c(rep(sqrt(.6),10), rep(.4,4), rep(0,6))
p2g1 <- c(rep(.4,4), rep(0,6), rep(sqrt(.6),10))
P1 <- cbind(p1g1, p2g1)
PSI1 <- diag(1 - rowSums(P1^2))
SIGMA1 <- P1 %*% t(P1) + PSI1

# Define 5-factor model loading matrix
p3g1 <- c(rep(0,5), rep(.4,5), rep(sqrt(.6),5), rep(0,5))
p4g1 <- c(rep(0,10), rep(.4,5), rep(sqrt(.6),5))
p5g1 <- c(rep(0,15), rep(.4,5))
P5 <- cbind(p1g1, p2g1, p3g1, p4g1, p5g1)

# Check for communality issue in P5
communalities <- rowSums(P5^2)
if (any(communalities >= 1)) {
  scale_factor <- sqrt(0.99 / max(communalities))
  P5 <- P5 * scale_factor
}

PSI5 <- diag(1 - rowSums(P5^2))
SIGMA5 <- P5 %*% t(P5) + PSI5

# Generate and save data for each sample size
for (n in n_values) {
  data_2factors <- mvrnorm(n = n, mu = rep(0, 20), Sigma = SIGMA1, empirical = TRUE)
  data_5factors <- mvrnorm(n = n, mu = rep(0, 20), Sigma = SIGMA5, empirical = TRUE)
  
  # Save datasets along with true loadings for evaluation
  saveRDS(list(data = data_2factors, true_loadings = P1), 
          file = file.path("data", paste0("simulated_data_n", n, "_2factors.RDS")))
  
  saveRDS(list(data = data_5factors, true_loadings = P5), 
          file = file.path("data", paste0("simulated_data_n", n, "_5factors.RDS")))
}
