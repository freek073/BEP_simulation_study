# Load necessary packages
library(lavaan)
library(regsem)

# Load simulated data, set sample sizes, factors, penalties, and reps
all_data <- readRDS("simulated_datasets_all.RDS")
sample_sizes <- c(20, 50, 200, 600, 1000)
factor_levels <- c(2, 5)
penalties <- c("ridge", "lasso")
reps <- 1:10

# Build model using lavaan from true loading matrix
build_lavaan_syntax <- function(loadings, prefix = "V") {
  lines <- character()
  for (f in 1:ncol(loadings)) {
    items <- which(loadings[, f] != 0)
    if (length(items) > 0) {
      lines <- c(lines, paste0("F", f, " =~ ", paste0(prefix, items, collapse = " + ")))
    }
  }
  paste(lines, collapse = "\n")
}

# Calculate congruence
calc_congruence <- function(est, true) {
  congruences <- sapply(1:ncol(true), function(i) {
    apply(est, 2, function(col) abs(cor(col, true[, i])))
  })
  mean(apply(congruences, 1, max))
}

# Calculate variance
calc_var_explained <- function(loadings) {
  sum(loadings^2) / nrow(loadings)
}

# Store output
results <- list()

# Main loop, iterate through the conditions
# Retrieve correct simlations
for (n in sample_sizes) {
  for (factors in factor_levels) {
    for (r in reps) {
      key <- paste0("n", n, "_rep", r, "_", ifelse(factors == 2, "2f", "5f"))
      if (!key %in% names(all_data)) next
      #prepare data to fit lavaan model
      sim <- all_data[[key]]
      data <- as.data.frame(sim$data)
      colnames(data) <- paste0("V", 1:ncol(data))
      true <- sim$true_loadings
      model <- build_lavaan_syntax(true)
      
      fit_lavaan <- tryCatch(cfa(model, data = data, std.lv = TRUE), error = function(e) NULL)
      if (is.null(fit_lavaan)) next
      if (!lavInspect(fit_lavaan, "converged")) next
      # Apply lasso and ridge
      for (penalty in penalties) {
        fit_pen <- tryCatch(
          regsem(fit_lavaan, type = penalty, lambda = 0.01,
                 pars_pen = "loadings", gradFun = "ram"),
          error = function(e) NULL
        )
        if (is.null(fit_pen)) next
        # Extract loadings
        est_params <- as.numeric(coef(fit_pen))
        lavaan_table <- parameterEstimates(fit_lavaan, standardized = TRUE)
        load_rows <- lavaan_table[lavaan_table$op == "=~", ]
        if (length(est_params) < nrow(load_rows)) next
        
        load_rows$est <- est_params[1:nrow(load_rows)]
        est_matrix <- matrix(0, nrow = nrow(true), ncol = ncol(true))
        rownames(est_matrix) <- paste0("V", 1:nrow(true))
        colnames(est_matrix) <- paste0("F", 1:ncol(true))
        
        for (i in 1:nrow(load_rows)) {
          row <- which(rownames(est_matrix) == load_rows$rhs[i])
          col <- which(colnames(est_matrix) == load_rows$lhs[i])
          if (length(row) > 0 && length(col) > 0) {
            est_matrix[row, col] <- load_rows$est[i]
          }
        }
        # Store results
        results[[length(results) + 1]] <- data.frame(
          Method = "regsem",
          Penalty = penalty,
          Factors = factors,
          SampleSize = n,
          Congruence = round(calc_congruence(est_matrix, true), 3),
          VarExplained = round(calc_var_explained(est_matrix), 3),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

# Final output
results_df <- do.call(rbind, results)
print(results_df)
saveRDS(results_df, "regsem_results_final_format.rds")




