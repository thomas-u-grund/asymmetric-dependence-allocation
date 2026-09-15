# ------------------------------------------------------------------
# 22: Figure 2 -- the paper's central results figure. Predicted
#     probability that a relationship is the one that changes, as a
#     function of its trade-dependence asymmetry relative to an
#     otherwise-comparable alternative (same starting sign, same
#     structural position, same capability/mutual-dependence), from
#     the M4 conditional logit (results/primary_table_progressive.rds).
#
#     Because a conditional logit's odds depend only on DIFFERENCES in
#     covariates between alternatives in the same choice set, comparing
#     two relationships that are identical except for dependence
#     asymmetry reduces to a standard logistic curve in the
#     standardized asymmetry difference, with slope equal to the
#     dependence-asymmetry coefficient. The 95% band uses the delta
#     method on that same coefficient's clustered standard error.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(ggplot2) })


prog <- readRDS("results/primary_table_progressive.rds")
m4 <- prog$m4  # a coeftest matrix (triad-clustered), not a raw coxph fit
b <- m4["z_dep_asym", "Estimate"]
se_clustered <- m4["z_dep_asym", "Std. Error"]
cat("M4 dependence-asymmetry coefficient:", round(b, 4), " clustered SE:", round(se_clustered, 4), "\n")

# recover the mean/sd used to standardize dep_asymmetry, so the x-axis can be
# shown on the raw (percentile) scale rather than in z-units. Restricted to
# the 10th-90th percentile range -- the same comparison reported in the
# text -- rather than the full distribution's long tail, so the figure's
# visual impression matches the paper's own careful framing of this as a
# modest, not dramatic, effect.
h3b <- readRDS("results/bloc_fe_and_h3_clustering.rds")
long <- h3b$long_b
mu <- mean(long$dep_asymmetry); sdv <- sd(long$dep_asymmetry)
pctiles <- quantile(long$dep_asymmetry, seq(0.10, 0.90, by = 0.01))
cat("dep_asymmetry range used (10th-90th pctile):", round(range(pctiles), 3), "\n")

# Compare the FOCAL relationship at percentile p against a comparison
# relationship at the MIRROR percentile (100-p) -- e.g. focal at the 90th
# vs. comparison at the 10th, and vice versa. This is the exact comparison
# the text's "moving dependence asymmetry from the 10th to the 90th
# percentile shifts the predicted probability from 0.435 to 0.565"
# describes (verified: focal=10th/comparison=90th gives 0.435,
# focal=90th/comparison=10th gives 0.565); holding the comparison fixed at
# one end throughout, instead, would give a different (and unreported)
# pair of endpoints.
pctiles_mirror <- quantile(long$dep_asymmetry, rev(seq(0.10, 0.90, by = 0.01)))
z_focal <- (pctiles - mu) / sdv
z_comparison <- (pctiles_mirror - mu) / sdv
z_diff <- z_focal - z_comparison

logit <- b * z_diff
se_logit <- abs(z_diff) * se_clustered
p_hat <- plogis(logit)
p_lo <- plogis(logit - 1.96 * se_logit)
p_hi <- plogis(logit + 1.96 * se_logit)

df <- tibble(dep_asymmetry = as.numeric(pctiles), percentile = seq(10, 90, by = 1),
             p_hat = p_hat, p_lo = p_lo, p_hi = p_hi)

write.csv(df, "results/figure2_data.csv", row.names = FALSE)

marker_df <- df %>% filter(percentile %in% c(10, 90))

p <- ggplot(df, aes(x = percentile, y = p_hat)) +
  geom_ribbon(aes(ymin = p_lo, ymax = p_hi), fill = "grey70", alpha = 0.5) +
  geom_line(linewidth = 1) +
  geom_point(data = marker_df, size = 2.5) +
  geom_text(data = marker_df, aes(label = sprintf("%.3f", p_hat)), vjust = -1, size = 4) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  scale_y_continuous(limits = c(0.39, 0.61)) +
  labs(x = "Percentile of trade-dependence asymmetry (focal relationship)",
       y = "Predicted probability focal relationship is the one that changes",
       title = "Dependence asymmetry shifts which relationship changes") +
  theme_minimal(base_size = 13)

ggsave("results/fig2_main_effect.pdf", p, width = 7, height = 5)
ggsave("results/fig2_main_effect.png", p, width = 7, height = 5, dpi = 300)

cat("\nAt 10th/90th percentile of asymmetry (mirrored comparison):\n")
p10 <- df$p_hat[which.min(abs(df$percentile - 10))]
p90 <- df$p_hat[which.min(abs(df$percentile - 90))]
cat("  10th pctile: P =", round(p10, 3), "\n")
cat("  90th pctile: P =", round(p90, 3), "\n")

cat("\n22_figure2_main_effect.R complete. Saved results/fig2_main_effect.pdf, .png, and figure2_data.csv\n")
