library(psych)
install.packages("GPArotation")

# Load generated data
all_data <- readRDS("simulated_datasets_all.RDS")
sample_sizes <- c(20, 50, 200, 600, 1000)
reps <- 100

# Function for calculating the congruence
calc_congruence <- function(estimated, true_loadings) {
  congruences <- sapply(1:ncol(true_loadings), function(i) {
    apply(estimated, 2, function(est_col) abs(cor(est_col, true_loadings[, i], use = "pairwise.complete.obs")))
  })
  mean(apply(congruences, 1, max))
}

# Function to calculate RMSE of loadings
calc_rmse <- function(estimated, true_loadings) {
  sqrt(mean((as.matrix(estimated) - true_loadings)^2))
}

# Function to calculate SRMR (psych)
calc_srmr <- function(fa_result) {
  if(!is.null(fa_result$residual)) {
    residuals <- fa_result$residual
    sqrt(mean(residuals[lower.tri(residuals)]^2))
  } else NA
}

# -------------------------factanal method-------------------------------------

results_list_factanal_2f <- list()
results_list_factanal_5f <- list()
rotations <- c("varimax", "promax")

# 2 factor 
for (n in sample_sizes) {
  for (r in 1:reps) {
    for (rotation in rotations) {
      key_2f <- paste0("n", n, "_rep", r, "_2f")
      if (key_2f %in% names(all_data)) {
        sim <- all_data[[key_2f]]
        data <- sim$data
        true_loadings_2f <- sim$true_loadings
        
        if (any(is.na(data))) next
        
        result_template <- data.frame(
          sample_size = n,
          replication = r,
          method = "factanal",
          rotation = rotation,
          factors = 2,
          var_explained = NA,
          congruence = NA,
          rmse = NA,
          srmr = NA,  # Not used in factanal
          p_value = NA,
          converged = FALSE,
          stringsAsFactors = FALSE
        )
        
        fa_factanal <- NULL
        if (n > ncol(data) + 5) {
          fa_factanal <- tryCatch(
            factanal(data, factors = 2, rotation = rotation),
            error = function(e) NULL
          )
        }
        
        if (!is.null(fa_factanal)) {
          current_result <- result_template
          current_result$var_explained <- sum(fa_factanal$loadings^2) / ncol(data)
          current_result$congruence <- calc_congruence(as.matrix(fa_factanal$loadings), true_loadings_2f)
          current_result$rmse <- calc_rmse(fa_factanal$loadings, true_loadings_2f)
          current_result$p_value <- if (!is.null(fa_factanal$PVAL)) fa_factanal$PVAL else NA
          current_result$converged <- TRUE
          
          results_list_factanal_2f[[length(results_list_factanal_2f) + 1]] <- current_result
        }
      }
    }
  }
}

# 5 factor 
for (n in sample_sizes) {
  for (r in 1:reps) {
    for (rotation in rotations) {
      key_5f <- paste0("n", n, "_rep", r, "_5f")
      if (key_5f %in% names(all_data)) {
        sim <- all_data[[key_5f]]
        data <- sim$data
        true_loadings_5f <- sim$true_loadings
        
        if (any(is.na(data))) next
        
        result_template <- data.frame(
          sample_size = n,
          replication = r,
          method = "factanal",
          rotation = rotation,
          factors = 5,
          var_explained = NA,
          congruence = NA,
          rmse = NA,
          srmr = NA,  # Not used in factanal
          p_value = NA,
          converged = FALSE,
          stringsAsFactors = FALSE
        )
        
        fa_factanal <- NULL
        if (n > ncol(data) + 5) {
          fa_factanal <- tryCatch(
            factanal(data, factors = 5, rotation = rotation),
            error = function(e) NULL
          )
        }
        
        if (!is.null(fa_factanal)) {
          current_result <- result_template
          current_result$var_explained <- sum(fa_factanal$loadings^2) / ncol(data)
          current_result$congruence <- calc_congruence(as.matrix(fa_factanal$loadings), true_loadings_5f)
          current_result$rmse <- calc_rmse(fa_factanal$loadings, true_loadings_5f)
          current_result$p_value <- if (!is.null(fa_factanal$PVAL)) fa_factanal$PVAL else NA
          current_result$converged <- TRUE
          
          results_list_factanal_5f[[length(results_list_factanal_5f) + 1]] <- current_result
        }
      }
    }
  }
}

# generate table factanal 
efa_results_factanal <- rbind(
  do.call(rbind, results_list_factanal_2f),
  do.call(rbind, results_list_factanal_5f)
)
saveRDS(efa_results_factanal, "efa_results_factanal_with_rotations.RDS")
summary_stats_factanal <- aggregate(
  cbind(var_explained, congruence, rmse, p_value) ~ method + rotation + factors + sample_size,
  data = efa_results_factanal,
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
)
print(summary_stats_factanal)


# -------------------------pysch method-------------------------------------


results_list_psych_2f <- list()
results_list_psych_5f <- list()
rotations <- c("varimax", "promax")

# ------------------ 2-FACTOR BLOCK ------------------
for (n in sample_sizes) {
  for (r in 1:reps) {
    for (rotation in rotations) {
      key_2f <- paste0("n", n, "_rep", r, "_2f")
      if (key_2f %in% names(all_data)) {
        sim <- all_data[[key_2f]]
        data <- sim$data
        true_loadings_2f <- sim$true_loadings
        
        if (any(is.na(data))) next
        
        result_template <- data.frame(
          sample_size = n,
          replication = r,
          method = "psych",
          rotation = rotation,
          factors = 2,
          var_explained = NA,
          congruence = NA,
          rmse = NA,
          srmr = NA,
          p_value = NA,
          converged = FALSE,
          stringsAsFactors = FALSE
        )
        
        fa_psych <- tryCatch({
          R <- cor(data)
          if (any(eigen(R)$values <= 0)) stop("Non-positive definite matrix")
          fa_result <- fa(R, nfactors = 2, rotate = rotation, fm = "pa", n.obs = nrow(data))
          
          current_result <- result_template
          current_result$var_explained <- sum(fa_result$loadings^2) / ncol(data)
          current_result$congruence <- calc_congruence(as.matrix(fa_result$loadings), true_loadings_2f)
          current_result$rmse <- calc_rmse(fa_result$loadings, true_loadings_2f)
          current_result$srmr <- calc_srmr(fa_result)
          current_result$converged <- TRUE
          
          current_result
        }, error = function(e) NULL)
        
        if (!is.null(fa_psych)) {
          results_list_psych_2f[[length(results_list_psych_2f) + 1]] <- fa_psych
        }
      }
    }
  }
}

# ------------------ 5-FACTOR BLOCK ------------------
for (n in sample_sizes) {
  for (r in 1:reps) {
    for (rotation in rotations) {
      key_5f <- paste0("n", n, "_rep", r, "_5f")
      if (key_5f %in% names(all_data)) {
        sim <- all_data[[key_5f]]
        data <- sim$data
        true_loadings_5f <- sim$true_loadings
        
        if (any(is.na(data))) next
        
        result_template <- data.frame(
          sample_size = n,
          replication = r,
          method = "psych",
          rotation = rotation,
          factors = 5,
          var_explained = NA,
          congruence = NA,
          rmse = NA,
          srmr = NA,
          p_value = NA,
          converged = FALSE,
          stringsAsFactors = FALSE
        )
        
        fa_psych <- tryCatch({
          R <- cor(data)
          if (any(eigen(R)$values <= 0)) stop("Non-positive definite matrix")
          fa_result <- fa(R, nfactors = 5, rotate = rotation, fm = "pa", n.obs = nrow(data))
          
          current_result <- result_template
          current_result$var_explained <- sum(fa_result$loadings^2) / ncol(data)
          current_result$congruence <- calc_congruence(as.matrix(fa_result$loadings), true_loadings_5f)
          current_result$rmse <- calc_rmse(fa_result$loadings, true_loadings_5f)
          current_result$srmr <- calc_srmr(fa_result)
          current_result$converged <- TRUE
          
          current_result
        }, error = function(e) NULL)
        
        if (!is.null(fa_psych)) {
          results_list_psych_5f[[length(results_list_psych_5f) + 1]] <- fa_psych
        }
      }
    }
  }
}

# generate table for psych 
efa_results_psych <- rbind(
  do.call(rbind, results_list_psych_2f),
  do.call(rbind, results_list_psych_5f)
)
saveRDS(efa_results_psych, "efa_results_psych_with_rotations.RDS")
summary_stats_psych <- aggregate(
  cbind(var_explained, congruence, rmse, srmr) ~ method + rotation + factors + sample_size,
  data = efa_results_psych,
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
)
print(summary_stats_psych)






