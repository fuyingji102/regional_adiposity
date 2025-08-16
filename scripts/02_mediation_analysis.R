
# mediation_analysis.R (Improved Version)
# ---------------------------------------------------------
# X = regional adiposity (e.g. trunk fat)
# M = brain system-level BAG (mediator)
# Y = cognitive performance (task-specific)
# Covariates: demographic, lifestyle, brain volume, etc.
# Output: direct, indirect effect, CI, z-scores, FDR-p
# ---------------------------------------------------------

library(lavaan)
library(dplyr)
library(readr)

# Define task/covariate info
continuous_tasks <- c("ExecFunTR", "NumericMemory", "TrailMakingNum1", "FluidIntelligence",
                      "PairsMatching6", "Symbol", "VerbalMem", "NonverbalReason",
                      "ProcSpeedRT", "TrailMakingAlp2", "Vocab")
binary_tasks <- c("ProsMem")

continuous_covariates <- c("age", "education_years", "physical_activity",
                           "alcohol_frequency", "smoking_frequency", "MetS_Score", "BrainSegVol")
binary_covariates <- c("sex", "ethnicity", "employed", "handedness")
all_covariates <- c(continuous_covariates, binary_covariates)

# Dummy loop variables (to be defined in actual use)
task_list <- c(...)              # e.g. c("ExecFunTR", ...)
obesity_indicators <- c(...)     # e.g. c("trunk", "vat_mass")
brain_systems <- c(...)          # e.g. c("CER", "BST")

# FDR function
fdr_correct <- function(p_values) {
  p.adjust(p_values, method = "fdr")
}

# Main loop
for (task in task_list) {
  # Read cognition data
  cog_data <- read_csv(...) %>%
    rename(eid = 1, cognition_score = 2)

  if (task %in% continuous_tasks) {
    cog_data <- cog_data %>% mutate(cognition_score = scale(cognition_score))
  }

  for (indicator in obesity_indicators) {
    adiposity_data <- read_csv(...) %>%
      rename(eid = 1, obesity_metric = 2) %>%
      mutate(obesity_metric = scale(obesity_metric))

    covariates <- read_csv(...) %>%
      rename(eid = 1) %>%
      mutate(across(all_of(continuous_covariates), scale)) %>%
      select(eid, all_of(all_covariates))

    for (system in brain_systems) {
      mediator_data <- read_csv(...) %>%
        rename(eid = 1, mediator = "BAG_corrected") %>%
        mutate(mediator = scale(mediator))

      merged_data <- adiposity_data %>%
        inner_join(mediator_data, by = "eid") %>%
        inner_join(cog_data, by = "eid") %>%
        inner_join(covariates, by = "eid") %>%
        filter(!is.na(mediator), !is.na(obesity_metric))

      # Define mediation model
      model_string <- paste0(
        "mediator ~ obesity_metric + ", paste(all_covariates, collapse = " + "), "\n",
        "cognition_score ~ obesity_metric + mediator + ", paste(all_covariates, collapse = " + ")
      )

      # Fit model
      fit <- sem(model_string, data = merged_data, se = "bootstrap", bootstrap = 5000)
      params <- parameterEstimates(fit, standardized = TRUE)

      # Extract paths
      a <- filter(params, lhs == "mediator", rhs == "obesity_metric")
      b <- filter(params, lhs == "cognition_score", rhs == "mediator")
      c_prime <- filter(params, lhs == "cognition_score", rhs == "obesity_metric")

      # Compute indirect effect
      indirect <- a$est * b$est
      indirect_se <- sqrt((a$est^2 * b$se^2) + (b$est^2 * a$se^2))
      indirect_ci_lower <- indirect - 1.96 * indirect_se
      indirect_ci_upper <- indirect + 1.96 * indirect_se
      z_indirect <- indirect / indirect_se
      p_indirect <- 2 * (1 - pnorm(abs(z_indirect)))

      # Compute direct effect
      direct_se <- c_prime$se
      direct_ci_lower <- c_prime$est - 1.96 * direct_se
      direct_ci_upper <- c_prime$est + 1.96 * direct_se
      z_direct <- c_prime$est / direct_se
      p_direct <- 2 * (1 - pnorm(abs(z_direct)))

      # Save results
      results <- data.frame(
        Task = task,
        Indicator = indicator,
        System = system,
        Path_a = a$est,
        Path_b = b$est,
        Path_c_prime = c_prime$est,
        Indirect = indirect,
        Indirect_SE = indirect_se,
        Indirect_CI_Lower = indirect_ci_lower,
        Indirect_CI_Upper = indirect_ci_upper,
        Z_indirect = z_indirect,
        P_indirect = p_indirect,
        Direct_SE = direct_se,
        Direct_CI_Lower = direct_ci_lower,
        Direct_CI_Upper = direct_ci_upper,
        Z_direct = z_direct,
        P_direct = p_direct
      )

      # FDR
      results$FDR_P_indirect <- fdr_correct(results$P_indirect)
      results$FDR_P_direct <- fdr_correct(results$P_direct)

      # Output
      write_csv(results, paste0("mediation_", task, "_", indicator, "_", system, ".csv"))
    }
  }
}
