# source/p4_make_plots.R

library(ggplot2)
library(dplyr)

bias_cov <- read.csv("results/summary/p4_bias_coverage_summary.csv")
time_sum <- read.csv("results/summary/p4_time_summary.csv")

bias_cov$n <- factor(bias_cov$n, levels = c(100, 1000, 10000))
time_sum$n <- factor(time_sum$n, levels = c(100, 1000, 10000))

p_bias <- ggplot(bias_cov, aes(x = n, y = mean_bias, color = method, group = method)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(
    aes(ymin = mean_bias - 1.96 * mcse_bias,
        ymax = mean_bias + 1.96 * mcse_bias),
    width = 0.15,
    position = position_dodge(width = 0.3)
  ) +
  facet_wrap(~ parameter, scales = "free_y") +
  labs(
    title = "Bias of mu-hat",
    x = "Sample size",
    y = "Monte Carlo mean bias"
  ) +
  theme_bw()

ggsave("results/figures/p4_bias_plot.png", p_bias, width = 10, height = 6, dpi = 300)

p_cov <- ggplot(bias_cov, aes(x = n, y = coverage, color = method, group = method)) +
  geom_hline(yintercept = 0.95, linetype = 2, color = "gray50") +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(
    aes(ymin = pmax(0, coverage - 1.96 * mcse_coverage),
        ymax = pmin(1, coverage + 1.96 * mcse_coverage)),
    width = 0.15,
    position = position_dodge(width = 0.3)
  ) +
  facet_wrap(~ parameter) +
  labs(
    title = "Coverage of mu-hat",
    x = "Sample size",
    y = "Coverage probability"
  ) +
  theme_bw()

ggsave("results/figures/p4_coverage_plot.png", p_cov, width = 10, height = 6, dpi = 300)

p_time <- ggplot(time_sum, aes(x = n, y = mean_time, color = method, group = method)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(
    aes(ymin = pmax(0, mean_time - 1.96 * mcse_time),
        ymax = mean_time + 1.96 * mcse_time),
    width = 0.15,
    position = position_dodge(width = 0.3)
  ) +
  scale_y_log10() +
  labs(
    title = "Computation time",
    x = "Sample size",
    y = "Mean elapsed time (seconds, log scale)"
  ) +
  theme_bw()

ggsave("results/figures/p4_time_plot.png", p_time, width = 8, height = 5, dpi = 300)
