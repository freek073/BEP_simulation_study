library(MASS)

# Define sample sizes and number of replications
n_values <- c(20, 50, 200, 600, 1000)
reps <- 100

# Seed for reproducibility

set.seed(45468) 

# 2 factor complexity 
p1g1 <- c(rep(sqrt(.6),10), rep(.4,4), rep(0,6))
p2g1 <- c(rep(.4,4), rep(0,6), rep(sqrt(.6),10))
P1 <- cbind(p1g1, p2g1)
PSI1 <- diag(1 - rowSums(P1^2))
SIGMA1 <- P1 %*% t(P1) + PSI1

# 5 factor complexity 
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


data_storage <- list()

# Data generation for the different sample sizes and complexities
# Multiple replications to capture sampling variations and estimate average performance

for (n in n_values) {
  for (r in 1:reps) {
    set.seed(100000 * n + r)  # reproducible per (n, rep)
    
    # Generate data
    data_2factor <- mvrnorm(n = n, mu = rep(0, 20), Sigma = SIGMA1, empirical = FALSE)
    data_5factor <- mvrnorm(n = n, mu = rep(0, 20), Sigma = SIGMA5, empirical = FALSE)
    
    # Store in list if needed for later use
    data_storage[[paste0("n", n, "_rep", r, "_2f")]] <- list(data = data_2factor, true_loadings = P1)
    data_storage[[paste0("n", n, "_rep", r, "_5f")]] <- list(data = data_5factor, true_loadings = P5)
  }
}

# Save data to be reproducible
saveRDS(data_storage, "simulated_datasets_all.RDS")
