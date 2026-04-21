library(tidyverse)

# =========================================================
# USER SETTINGS
# =========================================================

# Turn logistic prior for mu_alpha on or off
use_logistic_mu <- F

# If use_logistic_mu = TRUE:
#   mu_alpha ~ logistic(0, 1)
#
# If use_logistic_mu = FALSE:
#   mu_alpha ~ normal(mu_mean, mu_sd)
mu_mean <- 0
mu_sd   <- 1000

# sigma_alpha ~ half-normal(0, sigma_sd)
sigma_sd <- 100

# Simulation settings
n_hyper_draws   <- 5000
n_subj_per_draw <- 200

# =========================================================
# SIMULATION
# =========================================================

inv_logit <- plogis

if (use_logistic_mu) {
  mu_alpha <- rlogis(n_hyper_draws, location = 0, scale = 1)
  mu_label <- "mu_alpha ~ logistic(0, 1)"
} else {
  mu_alpha <- rnorm(n_hyper_draws, mean = mu_mean, sd = mu_sd)
  mu_label <- paste0("mu_alpha ~ normal(", mu_mean, ", ", mu_sd, ")")
}

sigma_alpha <- abs(rnorm(n_hyper_draws, mean = 0, sd = sigma_sd))
sigma_label <- paste0("sigma_alpha ~ half-normal(0, ", sigma_sd, ")")

sim_df <- map_dfr(seq_len(n_hyper_draws), function(d) {
  alpha_raw <- rnorm(n_subj_per_draw, 0, 1)
  alpha <- inv_logit(mu_alpha[d] + sigma_alpha[d] * alpha_raw)
  
  tibble(
    draw = d,
    mu_alpha = mu_alpha[d],
    sigma_alpha = sigma_alpha[d],
    alpha = alpha
  )
})

# =========================================================
# OUTPUT
# =========================================================

cat("\n-----------------------------------\n")
cat("Prior setup:\n")
cat(mu_label, "\n")
cat(sigma_label, "\n")
cat("-----------------------------------\n\n")

summary_tbl <- sim_df %>%
  summarise(
    mean_alpha = mean(alpha),
    sd_alpha   = sd(alpha),
    prop_lt_05 = mean(alpha < 0.05),
    prop_lt_10 = mean(alpha < 0.10),
    prop_gt_90 = mean(alpha > 0.90),
    prop_gt_95 = mean(alpha > 0.95)
  )

print(summary_tbl)

decile_tbl <- sim_df %>%
  mutate(decile = cut(alpha, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
  count(decile) %>%
  mutate(prop = n / sum(n))

print(decile_tbl)

ggplot(sim_df, aes(x = alpha)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60) +
  geom_hline(yintercept = 1, linetype = 2) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw() +
  labs(
    title = "Implied prior on subject-level alpha",
    subtitle = paste(mu_label, "|", sigma_label),
    x = "alpha",
    y = "Density"
  )