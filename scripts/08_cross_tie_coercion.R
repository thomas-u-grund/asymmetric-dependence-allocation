# ------------------------------------------------------------------
# 08: Cross-tie coercion test. H3 (script 05) showed the tie CARRYING
#     the trade-dependence asymmetry is disproportionately the one
#     that changes sign at realignment -- i.e. the exposed/dependent
#     relationship itself bends.
#
#     This tests a distinct, related mechanism: does having a highly
#     asymmetric tie ELSEWHERE in the triad make a DIFFERENT tie (not
#     itself asymmetric) more likely to be the one that changes? The
#     "US leans on a dependent client's economic exposure to force the
#     client to break with a third country" story -- the leverage sits
#     in one relationship, the pressure lands in another. Since every
#     triad has only 3 nodes, any two of its three edges automatically
#     share a vertex (there is no "unconnected" pair in a triangle), so
#     "another tie in the same triad" already is the adjacency the
#     coercion story needs -- no extra node-matching required.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest) })


# reuse the bloc-tagged single-tie-change sample already built in script 07
b <- readRDS("results/bloc_fe_and_h3_clustering.rds")
strict_b <- b$strict_b  # one row per event, all 3 ties trade-valid, primary_bloc attached

cat("Single-tie-change events (all 3 ties trade-valid, bloc-tagged):", nrow(strict_b), "\n\n")

# --- for each tie, the max dep_asymmetry of the OTHER TWO ties in the same triad ---
strict_b <- strict_b %>%
  rowwise() %>%
  mutate(
    other_dep_asym_12 = max(dep_asymmetry_13, dep_asymmetry_23, na.rm = TRUE),
    other_dep_asym_13 = max(dep_asymmetry_12, dep_asymmetry_23, na.rm = TRUE),
    other_dep_asym_23 = max(dep_asymmetry_12, dep_asymmetry_13, na.rm = TRUE)
  ) %>%
  ungroup()

long <- strict_b %>%
  select(event_id, triad_id, primary_bloc, ul12, ul13, ul23, emb12, emb13, emb23,
         dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23,
         other_dep_asym_12, other_dep_asym_13, other_dep_asym_23,
         changed12, changed13, changed23) %>%
  pivot_longer(
    cols = c(ul12, ul13, ul23, emb12, emb13, emb23,
             dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23,
             other_dep_asym_12, other_dep_asym_13, other_dep_asym_23,
             changed12, changed13, changed23),
    names_to = c(".value", "tie"),
    names_pattern = "(ul|emb|dep_asymmetry|other_dep_asym|changed)_?(12|13|23)"
  ) %>%
  rename(unbal_load = ul, embeddedness = emb) %>%
  mutate(z_unbal_load = as.numeric(scale(unbal_load)),
         z_dep_asym = as.numeric(scale(dep_asymmetry)),
         z_other_dep_asym = as.numeric(scale(other_dep_asym)),
         z_emb = as.numeric(scale(embeddedness)))

cat("Correlation, own dep_asymmetry vs. other_dep_asym_max (within tie rows):",
    round(cor(long$dep_asymmetry, long$other_dep_asym, use = "complete.obs"), 3), "\n\n")

# method="efron" for cluster-robust SEs (identical to "exact" point estimates
# here, since exactly one event per stratum -- see script 07's note)
clog_own_only  <- coxph(Surv(rep(1, nrow(long)), changed) ~ z_unbal_load + z_dep_asym + z_emb + strata(event_id),
                         data = long, method = "efron")
clog_cross     <- coxph(Surv(rep(1, nrow(long)), changed) ~ z_unbal_load + z_dep_asym + z_other_dep_asym + z_emb + strata(event_id),
                         data = long, method = "efron")
clog_cross_int <- coxph(Surv(rep(1, nrow(long)), changed) ~ z_unbal_load + z_dep_asym * z_other_dep_asym + z_emb + strata(event_id),
                         data = long, method = "efron")

ct_own_only <- coeftest(clog_own_only, vcov = vcovCL(clog_own_only, cluster = long$triad_id))
ct_cross_1way <- coeftest(clog_cross, vcov = vcovCL(clog_cross, cluster = long$triad_id))
cl2 <- long[, c("triad_id", "primary_bloc")]
ct_cross_2way <- coeftest(clog_cross, vcov = vcovCL(clog_cross, cluster = cl2, multi0 = TRUE))
ct_cross_int_1way <- coeftest(clog_cross_int, vcov = vcovCL(clog_cross_int, cluster = long$triad_id))
ct_cross_int_2way <- coeftest(clog_cross_int, vcov = vcovCL(clog_cross_int, cluster = cl2, multi0 = TRUE))

sink("results/cross_tie_coercion_output.txt")
cat("Single-tie-change events used:", nrow(strict_b), " | long-format rows:", nrow(long), "\n")
cat("Correlation, own dep_asymmetry vs. other_dep_asym_max:",
    round(cor(long$dep_asymmetry, long$other_dep_asym, use = "complete.obs"), 3), "\n\n")

cat("=== BASELINE: own dep_asymmetry only (reproduces H3, this subsample) ===\n")
print(ct_own_only)

cat("\n\n=== CROSS-TIE TEST: does asymmetry ELSEWHERE in the triad also predict this tie changing? ===\n")
cat("--- 1-way (triad) clustered SEs ---\n"); print(ct_cross_1way)
cat("\n--- 2-way (triad + bloc) clustered SEs ---\n"); print(ct_cross_2way)

cat("\n\n=== Interaction: does own asymmetry and elsewhere-asymmetry substitute/compound? ===\n")
cat("--- 1-way ---\n"); print(ct_cross_int_1way)
cat("\n--- 2-way ---\n"); print(ct_cross_int_2way)
sink()

cat(readLines("results/cross_tie_coercion_output.txt"), sep = "\n")

saveRDS(list(long = long, clog_own_only = clog_own_only, clog_cross = clog_cross, clog_cross_int = clog_cross_int,
             ct_own_only = ct_own_only, ct_cross_1way = ct_cross_1way, ct_cross_2way = ct_cross_2way,
             ct_cross_int_1way = ct_cross_int_1way, ct_cross_int_2way = ct_cross_int_2way),
        "results/cross_tie_coercion.rds")
cat("\n\n08_cross_tie_coercion.R complete.\n")
