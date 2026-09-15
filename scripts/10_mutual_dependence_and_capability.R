# ------------------------------------------------------------------
# 10: Two checks on the tie-choice model.
#
#     (1) The volume control (script 04/05/08) tested whether the
#     asymmetry effect was a small-denominator artifact, but absolute
#     dyadic trade volume is NOT the same thing as mutual dependence.
#     Construct mutual dependence directly from the same ingredients
#     already used for asymmetry:
#         d_ij = trade_share_low, d_ji = trade_share_high
#         mutual_dependence  = (d_ij + d_ji) / 2      (arithmetic mean)
#         mutual_dependence_min = min(d_ij, d_ji)     (weak-link version)
#         dep_asymmetry      = |d_ij - d_ji|           (already have this)
#     and test: holding overall mutual dependence constant, does
#     asymmetry still predict which tie changes?
#
#     (2) Does capability (CINC) asymmetry ALSO predict which tie
#     changes, the way trade-dependence asymmetry does (script 05/07),
#     or is its role confined to the triad-level hazard (raises
#     dissolution, lowers realignment, script 09)? If trade dependence
#     and capability asymmetry produce genuinely different adjustment
#     modes (one on WHICH tie changes, one on WHETHER the triad
#     resolves at all), that is a sharper theoretical claim than either
#     alone -- but it needs to be tested in the SAME model, not
#     asserted from two separate analyses.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest) })


b <- readRDS("results/bloc_fe_and_h3_clustering.rds")
strict_b <- b$strict_b  # one row per event: triad_id, t1, spell_id, primary_bloc,
                         # ul12/13/23, emb12/13/23, dep_asymmetry_12/13/23,
                         # cap_asymmetry_12/13/23, changed12/13/23

cat("Events (all 3 ties trade-valid, bloc-tagged):", nrow(strict_b), "\n")
stopifnot("spell_id" %in% names(strict_b))

# --- re-merge tie_trade for the raw shares, to build mutual dependence ---
tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high, trade_share_low, trade_share_high)

get_shares <- function(df, a, b_, suffix) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b_]]), hi = pmax(.data[[a]], .data[[b_]])) %>%
    left_join(tie_trade, by = c("t1" = "year", "lo" = "ccode_low", "hi" = "ccode_high")) %>%
    rename_with(~ paste0(., "_", suffix), c(trade_share_low, trade_share_high)) %>%
    select(-lo, -hi)
}

strict_b <- strict_b %>%
  get_shares("node1", "node2", "12") %>%
  get_shares("node1", "node3", "13") %>%
  get_shares("node2", "node3", "23") %>%
  mutate(
    mutual_dep_12 = (trade_share_low_12 + trade_share_high_12) / 2,
    mutual_dep_13 = (trade_share_low_13 + trade_share_high_13) / 2,
    mutual_dep_23 = (trade_share_low_23 + trade_share_high_23) / 2,
    mutual_dep_min_12 = pmin(trade_share_low_12, trade_share_high_12),
    mutual_dep_min_13 = pmin(trade_share_low_13, trade_share_high_13),
    mutual_dep_min_23 = pmin(trade_share_low_23, trade_share_high_23)
  )

long <- strict_b %>%
  select(event_id, triad_id, primary_bloc, ul12, ul13, ul23, emb12, emb13, emb23,
         dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23,
         cap_asymmetry_12, cap_asymmetry_13, cap_asymmetry_23,
         mutual_dep_12, mutual_dep_13, mutual_dep_23,
         mutual_dep_min_12, mutual_dep_min_13, mutual_dep_min_23,
         changed12, changed13, changed23) %>%
  pivot_longer(
    cols = c(ul12, ul13, ul23, emb12, emb13, emb23,
             dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23,
             cap_asymmetry_12, cap_asymmetry_13, cap_asymmetry_23,
             mutual_dep_12, mutual_dep_13, mutual_dep_23,
             mutual_dep_min_12, mutual_dep_min_13, mutual_dep_min_23,
             changed12, changed13, changed23),
    names_to = c(".value", "tie"),
    names_pattern = "(ul|emb|dep_asymmetry|cap_asymmetry|mutual_dep_min|mutual_dep|changed)_?(12|13|23)"
  ) %>%
  rename(unbal_load = ul, embeddedness = emb) %>%
  mutate(z_unbal_load = as.numeric(scale(unbal_load)),
         z_dep_asym = as.numeric(scale(dep_asymmetry)),
         z_cap_asym = as.numeric(scale(cap_asymmetry)),
         z_mutual = as.numeric(scale(mutual_dep)),
         z_mutual_min = as.numeric(scale(mutual_dep_min)),
         z_emb = as.numeric(scale(embeddedness)))

cat("Correlation, dep_asymmetry vs. mutual_dep (mean):", round(cor(long$dep_asymmetry, long$mutual_dep, use="complete.obs"),3), "\n")
cat("Correlation, dep_asymmetry vs. mutual_dep_min (weak-link):", round(cor(long$dep_asymmetry, long$mutual_dep_min, use="complete.obs"),3), "\n\n")

fit2way <- function(f, data) {
  # do.call, not a direct coxph(f, data=data) call, to avoid a known
  # coxph quirk where its internal match.call()-based re-evaluation of
  # the data argument breaks when coxph is invoked from inside a wrapper
  # function rather than directly at top level.
  m <- do.call(coxph, list(formula = f, data = data, method = "efron"))
  cl2 <- data[, c("triad_id", "primary_bloc")]
  list(m = m,
       ct_1way = coeftest(m, vcov = vcovCL(m, cluster = data$triad_id)),
       ct_2way = coeftest(m, vcov = vcovCL(m, cluster = cl2, multi0 = TRUE)))
}

surv <- Surv(rep(1, nrow(long)), long$changed)

# (1) mutual dependence (mean) as the "how important is this relationship"
#     control, replacing raw volume
r_mutual <- fit2way(surv ~ z_unbal_load + z_dep_asym + z_mutual + z_emb + strata(event_id), long)
# robustness: weak-link (min) version of mutual dependence
r_mutual_min <- fit2way(surv ~ z_unbal_load + z_dep_asym + z_mutual_min + z_emb + strata(event_id), long)

# (2) capability asymmetry in the SAME tie-choice model as trade asymmetry
r_cap_only <- fit2way(surv ~ z_unbal_load + z_cap_asym + z_emb + strata(event_id), long)
r_both     <- fit2way(surv ~ z_unbal_load + z_dep_asym + z_cap_asym + z_emb + strata(event_id), long)
r_full     <- fit2way(surv ~ z_unbal_load + z_dep_asym + z_cap_asym + z_mutual + z_emb + strata(event_id), long)
# the defensible "everything at once" spec, using the weak-link (min)
# mutual-dependence measure rather than the near-collinear mean version
r_full_min <- fit2way(surv ~ z_unbal_load + z_dep_asym + z_cap_asym + z_mutual_min + z_emb + strata(event_id), long)

sink("results/mutual_dependence_and_capability_output.txt")
cat("Events:", nrow(strict_b), " | long rows:", nrow(long), "\n")
cat("cor(dep_asymmetry, mutual_dep):", round(cor(long$dep_asymmetry, long$mutual_dep, use="complete.obs"),3), "\n")
cat("cor(dep_asymmetry, mutual_dep_min):", round(cor(long$dep_asymmetry, long$mutual_dep_min, use="complete.obs"),3), "\n\n")

cat("=== (1) dep_asymmetry holding MUTUAL DEPENDENCE (mean) constant ===\n")
cat("--- 1-way ---\n"); print(r_mutual$ct_1way)
cat("\n--- 2-way (triad+bloc) ---\n"); print(r_mutual$ct_2way)

cat("\n\n=== (1) robustness: mutual dependence as weak-link (min) ===\n")
cat("--- 1-way ---\n"); print(r_mutual_min$ct_1way)
cat("\n--- 2-way ---\n"); print(r_mutual_min$ct_2way)

cat("\n\n=== (2) capability asymmetry ALONE in the tie-choice model ===\n")
cat("--- 1-way ---\n"); print(r_cap_only$ct_1way)
cat("\n--- 2-way ---\n"); print(r_cap_only$ct_2way)

cat("\n\n=== (2) trade dep_asymmetry AND cap_asymmetry together ===\n")
cat("--- 1-way ---\n"); print(r_both$ct_1way)
cat("\n--- 2-way ---\n"); print(r_both$ct_2way)

cat("\n\n=== full model: dep_asymmetry + cap_asymmetry + mutual dependence (MEAN, near-collinear -- not informative, shown for completeness) ===\n")
cat("--- 1-way ---\n"); print(r_full$ct_1way)
cat("\n--- 2-way ---\n"); print(r_full$ct_2way)

cat("\n\n=== full model, DEFENSIBLE spec: dep_asymmetry + cap_asymmetry + mutual dependence (MIN, weak-link) ===\n")
cat("--- 1-way ---\n"); print(r_full_min$ct_1way)
cat("\n--- 2-way ---\n"); print(r_full_min$ct_2way)
sink()

cat(readLines("results/mutual_dependence_and_capability_output.txt"), sep = "\n")

saveRDS(list(long = long, r_mutual = r_mutual, r_mutual_min = r_mutual_min,
             r_cap_only = r_cap_only, r_both = r_both, r_full = r_full, r_full_min = r_full_min),
        "results/mutual_dependence_and_capability.rds")
cat("\n\n10_mutual_dependence_and_capability.R complete.\n")
