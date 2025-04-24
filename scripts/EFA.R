# In this script EFA is performed 
# Perfornning EFA on multiple different sample sizes 

# Call the data set that corresponds with the correct sample size
sample_sizes <- c(20, 50, 200, 600, 1000)
for (n in sample_sizes) {
  file_path <- paste0("data/simulated_data_n", n, "_2factors.RDS")
  
  # Skips if a certain sample wasn't generated
  if (!file.exists(file_path)) next
  
  # Loading the datasets
  sim <- readRDS(file_path)
  data <- sim$data  
  
  # perform EFA without rotation
  # trycatch to not crash the script
  efa_none <- tryCatch({
    factanal(data, factors = 2, rotation = "none")
  }, error = function(e) {
    cat("EFA failed for n =", n, ":", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(efa_none)) {
    cat("\n EFA for n =", n, "\n")
    print(efa_result$loadings, digits = 2, cutoff = 0.3)
    cat("P-value:", efa_none$PVAL, "\n")
  }

  # Perform EFA with varimax rotation
  efa_varimax <- tryCatch({
    factanal(data, factors = 2, rotation = "varimax")
  }, error = function(e) {
    cat("EFA (varimax) failed for n =", n, ":", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(efa_varimax)) {
    cat("\n EFA for n =", n, "Rotation: varimax \n")
    print(efa_varimax$loadings, digits = 2, cutoff = 0.3)
    cat("P-value:", efa_varimax$PVAL, "\n")
  }
}

