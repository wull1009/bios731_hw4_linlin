# source/p4_aggregate_results.R

library(dplyr)

files <- list.files("results/raw", pattern = "^res_n.*\\.csv$", full.names = TRUE)

if (length(files) == 0) {
  stop("No raw result files found in results/raw")
}

dat <- bind_rows(lapply(files, read.csv))

write.csv(dat, "results/summary/p4_all_results.csv", row.names = FALSE)

summary_bias_cov <- dat %>%
  group_by(n, method, parameter, true_mu) %>%
  summarise(
    nsim = n(),
    mean_bias = mean(bias),
    mcse_bias = sd(bias) / sqrt(n()),
    coverage = mean(cover),
    mcse_coverage = sqrt(coverage * (1 - coverage) / n()),
    .groups = "drop"
  )

write.csv(summary_bias_cov,
          "results/summary/p4_bias_coverage_summary.csv",
          row.names = FALSE)

summary_time <- dat %>%
  distinct(n, sim_id, method, time_sec) %>%
  group_by(n, method) %>%
  summarise(
    nsim = n(),
    mean_time = mean(time_sec),
    mcse_time = sd(time_sec) / sqrt(n()),
    .groups = "drop"
  )

write.csv(summary_time,
          "results/summary/p4_time_summary.csv",
          row.names = FALSE)
