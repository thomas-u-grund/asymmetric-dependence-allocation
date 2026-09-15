# ------------------------------------------------------------------
# 04: H1 (dissolution suppressed) and H2 (realignment suppressed) --
#     discrete-time cause-specific hazard models on the trade-valid
#     subsample, with spell-duration bins and era dummies,
#     cluster-robust SEs by triad.
#
#     Nested comparison per outcome:
#       (0) baseline structural controls only
#       (1) + symmetric trade dependence (mean across 3 ties) --
#           a baseline control, for a same-sample comparison
#       (2) + dep_asym_mean (mean dependence asymmetry across 3 ties)
#       (3) + dep_asym_max  (worst-tie dependence asymmetry) -- the
#           sharper "does this triad contain a trapped tie" test
#       (4) (3) + cap_asym_max, to ask whether the asymmetry effect
#           survives controlling for capability asymmetry (H4)
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(sandwich); library(lmtest)
})


py <- readRDS("data/person_year_risk_table.rds")

# Trade-valid subsample: requires ALL 3 ties to have a valid dep_asymmetry
# reading, since dep_asym_mean/max are only meaningful when computed over
# a complete set of 3 ties (a triad missing 2 of 3 ties' data would have
# its "mean"/"max" driven by a single, arbitrarily-observed tie).
py_trade <- py %>% filter(n_valid_ties == 3, year >= 1870) %>%
  mutate(
    z_trade_mean   = as.numeric(scale(trade_dependence_mean)),
    z_dep_asym_mean = as.numeric(scale(dep_asym_mean)),
    z_dep_asym_max  = as.numeric(scale(dep_asym_max)),
    z_cap_asym_max  = as.numeric(scale(cap_asym_max)),
    z_volume_mean   = as.numeric(scale(trade_volume_mean)),
    z_volume_max    = as.numeric(scale(trade_volume_max))
  )
cat("Trade-valid analysis sample:", nrow(py_trade), "person-years,",
    n_distinct(py_trade$triad_id), "triads,",
    n_distinct(py_trade$spell_id), "spells\n")
cat("Year range:", min(py_trade$year), "-", max(py_trade$year), "\n\n")

cat("Correlation, dep_asym_max vs. trade_volume_mean (worst-tie asymmetry vs. avg dyadic trade size):",
    round(cor(py_trade$dep_asym_max, py_trade$trade_volume_mean, use = "complete.obs"), 3), "\n")
cat("Correlation, dep_asym_max vs. trade_volume_max:",
    round(cor(py_trade$dep_asym_max, py_trade$trade_volume_max, use = "complete.obs"), 3), "\n\n")

fit_pair <- function(rhs, data) {
  f_dis  <- as.formula(paste("event_dissolve ~", rhs))
  f_flip <- as.formula(paste("event_flip ~", rhs))
  m_dis  <- glm(f_dis,  data = data, family = binomial())
  m_flip <- glm(f_flip, data = data, family = binomial())
  ct_dis  <- coeftest(m_dis,  vcov = vcovCL(m_dis,  cluster = data$triad_id))
  ct_flip <- coeftest(m_flip, vcov = vcovCL(m_flip, cluster = data$triad_id))
  list(m_dis = m_dis, m_flip = m_flip, ct_dis = ct_dis, ct_flip = ct_flip)
}

base_rhs   <- "z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin"
m0 <- fit_pair(base_rhs, py_trade)
m1 <- fit_pair(paste(base_rhs, "+ z_trade_mean"), py_trade)
m2 <- fit_pair(paste(base_rhs, "+ z_trade_mean + z_dep_asym_mean"), py_trade)
m3 <- fit_pair(paste(base_rhs, "+ z_trade_mean + z_dep_asym_max"), py_trade)
m4 <- fit_pair(paste(base_rhs, "+ z_trade_mean + z_dep_asym_max + z_cap_asym_max"), py_trade)
m5 <- fit_pair(paste(base_rhs, "+ z_trade_mean + z_dep_asym_max * z_cap_asym_max"), py_trade)

# --- does dep_asym_max survive controlling for the ABSOLUTE size of the
#     dyadic trade relationship, not just its share of each side's total
#     trade? This is the key check: is asymmetry genuine entanglement, or
#     a small-denominator artifact of trivial relationships? ---
m6 <- fit_pair(paste(base_rhs, "+ z_trade_mean + z_volume_mean + z_dep_asym_max"), py_trade)
m7 <- fit_pair(paste(base_rhs, "+ z_trade_mean + z_volume_max + z_dep_asym_max"), py_trade)

# --- H4 robustness: split sample by whether the triad contains at least
#     one tie where the trade-exposed side is also the capability-weaker
#     side (any_exposed_weaker) ---
py_weak <- py_trade %>% filter(any_exposed_weaker)
py_notweak <- py_trade %>% filter(!any_exposed_weaker)
cat("Split for H4: any_exposed_weaker TRUE:", nrow(py_weak),
    " | FALSE:", nrow(py_notweak), "\n\n")
m_weak    <- fit_pair(paste(base_rhs, "+ z_trade_mean + z_dep_asym_max"), py_weak)
m_notweak <- fit_pair(paste(base_rhs, "+ z_trade_mean + z_dep_asym_max"), py_notweak)

sink("results/hazard_models_output.txt")
cat("=== Trade-valid analysis sample ===\n")
cat("N =", nrow(py_trade), "person-years,", n_distinct(py_trade$triad_id), "triads,",
    n_distinct(py_trade$spell_id), "spells, years", min(py_trade$year), "-", max(py_trade$year), "\n\n")

report <- function(tag, m) {
  cat("\n\n============================================================\n")
  cat(tag, "\n============================================================\n")
  cat("--- DISSOLUTION ---\n"); print(m$ct_dis); cat("N =", nobs(m$m_dis), " AIC =", AIC(m$m_dis), "\n")
  cat("\n--- REALIGNMENT (flip) ---\n"); print(m$ct_flip); cat("N =", nobs(m$m_flip), " AIC =", AIC(m$m_flip), "\n")
}

report("MODEL 0: baseline structural controls only", m0)
report("MODEL 1: + symmetric trade dependence (mean, parent-project control)", m1)
report("MODEL 2 [H1/H2, mean asymmetry]: + dep_asym_mean", m2)
report("MODEL 3 [H1/H2, PRIMARY]: + dep_asym_max (worst-tie asymmetry)", m3)
report("MODEL 4 [H4]: + dep_asym_max + cap_asym_max", m4)
report("MODEL 5 [H4, interaction]: + dep_asym_max * cap_asym_max", m5)
report("MODEL 6 [volume check]: + trade_volume_mean + dep_asym_max", m6)
report("MODEL 7 [volume check]: + trade_volume_max (worst-tie) + dep_asym_max", m7)

cat("\n\n=== LRT: does adding dep_asym_max improve fit over trade-dependence-only model? ===\n")
cat("--- Dissolution ---\n"); print(anova(m1$m_dis, m3$m_dis, test = "Chisq"))
cat("--- Realignment ---\n"); print(anova(m1$m_flip, m3$m_flip, test = "Chisq"))

report("MODEL H4-split: any_exposed_weaker == TRUE (small/dependent state IS the weaker one)", m_weak)
report("MODEL H4-split: any_exposed_weaker == FALSE (asymmetric dependence without capability gap)", m_notweak)
sink()

cat(readLines("results/hazard_models_output.txt"), sep = "\n")

saveRDS(list(py_trade = py_trade, m0 = m0, m1 = m1, m2 = m2, m3 = m3, m4 = m4, m5 = m5,
             m6 = m6, m7 = m7, m_weak = m_weak, m_notweak = m_notweak),
        "results/hazard_models.rds")
cat("\n\n04_hazard_models.R complete. See results/hazard_models_output.txt\n")
