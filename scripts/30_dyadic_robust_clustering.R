# ------------------------------------------------------------------
# 30: Replaces script 29's connected-component approach, which an
#     external review correctly identified as wrong: clustering by
#     connected component imposes dependence *transitively* across an
#     entire connected chain (A-B, B-C, C-D all end up in one cluster
#     even though A-B and C-D share no actor), which is a much broader
#     -- and unjustified -- notion of dependence than "these two
#     relationships share an actor." That transitive closure, not a
#     genuine property of the data, is why 127/139 states collapsed
#     into a single component in script 29.
#
#     This script instead implements the actor-overlap ("dyadic-
#     robust") variance estimator of Cameron and Miller (2014, Section
#     VI) and Aronow, Samii, and Assenova (2015): two observations are
#     treated as correlated if and only if they share at least one of
#     their two constituent actors directly, with no transitive
#     closure through a third observation. Unlike connected-component
#     clustering, this does not require partitioning the data into
#     disjoint groups, so it is well-defined even when the actor-
#     overlap graph is densely connected.
#
#     Derivation used below: let u_i be observation i's score
#     contribution (estfun), and for each actor g let S_g = sum of
#     u_i over every observation touching g (in either position).
#     Sum_g S_g S_g' counts each ordered pair (i,j) once for every
#     actor they share (0, 1, or 2, since each observation has exactly
#     two actors) -- so it double-counts pairs sharing both actors,
#     i.e. pairs belonging to the same dyad (including each
#     observation with itself, which always "shares 2 actors" with
#     itself). Subtracting the standard one-way dyad-clustered meat
#     (which is exactly that double-counted piece) leaves each
#     dependent pair counted exactly once:
#
#         Meat_dyadic = sum_g S_g S_g'  -  Meat_dyad
#
#     Validated below against sandwich::vcovCL's own one-way dyad
#     clustering before use (same bread, same meat construction,
#     reproduces vcovCL exactly when the two formulas should agree).
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(sandwich); library(lmtest) })


prev <- readRDS("results/full_multiway_clustering.rds")
long4 <- prev$long4
m4 <- prev$m4

U <- estfun(m4)
B <- bread(m4)
n <- nrow(U)
k <- ncol(U)

# ------------------------------------------------------------------
# Validation: reproduce vcovCL's own one-way dyad-clustered vcov from
# first principles, to pin down the exact small-sample scale factor
# sandwich applies (cadjust = TRUE default: (G/(G-1)) applied to the
# bread %*% meat %*% bread / n^2 sandwich formula), before trusting
# the same construction for the actor-overlap meat below.
# ------------------------------------------------------------------
ct_dyad_pkg <- coeftest(m4, vcov = vcovCL(m4, cluster = long4$dyad))
S_by_dyad <- rowsum(U, group = long4$dyad)
meat_dyad_raw <- crossprod(S_by_dyad)
G_dyad <- nrow(S_by_dyad)
vcov_dyad_manual <- (B %*% meat_dyad_raw %*% B) * (G_dyad / (G_dyad - 1)) / n^2
max_rel_err <- max(abs(sqrt(diag(vcov_dyad_manual)) - ct_dyad_pkg[, "Std. Error"]) / ct_dyad_pkg[, "Std. Error"])
cat("Validation -- max relative error, manual vs. vcovCL one-way dyad clustering:", max_rel_err, "\n")
stopifnot(max_rel_err < 1e-6)
cat("Validation passed: manual construction exactly reproduces sandwich::vcovCL.\n\n")

# ------------------------------------------------------------------
# Actor-overlap (dyadic-robust) meat: multi-membership sum over
# states minus the same dyad-level double-count identified above
# ------------------------------------------------------------------
states_long <- bind_rows(
  long4 %>% transmute(row = row_number(), state = state_lo),
  long4 %>% transmute(row = row_number(), state = state_hi)
)
S_by_state <- rowsum(U[states_long$row, , drop = FALSE], group = states_long$state)
meat_state_raw <- crossprod(S_by_state)          # sum_g S_g S_g'
meat_dyadic <- meat_state_raw - meat_dyad_raw     # subtract the double-counted same-dyad piece

G_state <- n_distinct(c(long4$state_lo, long4$state_hi))
vcov_dyadic <- (B %*% meat_dyadic %*% B) * (G_state / (G_state - 1)) / n^2

se_dyadic <- sqrt(diag(vcov_dyadic))
z_dyadic <- coef(m4) / se_dyadic
p_dyadic <- 2 * pnorm(-abs(z_dyadic))

cat("Number of distinct actors (states):", G_state, "\n\n")
cat("=== Actor-overlap (dyadic-robust) clustered M4 ===\n")
print(data.frame(Estimate = coef(m4), `Std. Error` = se_dyadic, z = z_dyadic, p = p_dyadic, check.names = FALSE))

sink("results/dyadic_robust_clustering_output.txt")
cat("=== Validation: manual one-way dyad clustering vs. sandwich::vcovCL ===\n")
cat("Max relative error:", max_rel_err, "(passed, < 1e-6)\n\n")
cat("=== Actor-overlap (dyadic-robust, Cameron-Miller / Aronow-Samii-Assenova) clustering on M4 ===\n")
cat("N events =", n_distinct(long4$event_id), " N tie-rows =", nrow(long4), " distinct actors =", G_state, "\n\n")
print(data.frame(Estimate = coef(m4), `Std. Error` = se_dyadic, z = z_dyadic, p = p_dyadic, check.names = FALSE))
cat("\n\nFor comparison, one-way triad-clustered and triad+bloc+dyad results (script 28):\n")
print(prev$ct_1way)
cat("\n"); print(prev$ct_3way)
sink()
cat(readLines("results/dyadic_robust_clustering_output.txt"), sep = "\n")

saveRDS(list(meat_dyadic = meat_dyadic, vcov_dyadic = vcov_dyadic, se_dyadic = se_dyadic,
             p_dyadic = p_dyadic, G_state = G_state),
        "results/dyadic_robust_clustering.rds")
cat("\n\n30_dyadic_robust_clustering.R complete.\n")
