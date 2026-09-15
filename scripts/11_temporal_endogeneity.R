# ------------------------------------------------------------------
# 11: Response to external review, priority item 3 (temporal
#     endogeneity). Every result so far measures trade-dependence
#     asymmetry AT t1, the year immediately before the realignment
#     event itself -- i.e. contemporaneous with the resolution. That
#     leaves open reverse causality: if a relationship's trade is
#     already shifting BECAUSE political realignment is imminent (or
#     underway), "asymmetry predicts which tie changes" could partly
#     just be "the tie that's about to change already has unusual
#     trade numbers by the time we measure it."
#
#     Fix: measure dep_asymmetry (and the other tie-level covariates)
#     at the SPELL'S ONSET year -- before the imbalance that leads to
#     this realignment even began -- instead of at t1. Also test
#     stricter lags (onset year minus 1/3/5), so the measurement
#     predates not just the resolution but the imbalance itself.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest) })


b <- readRDS("results/bloc_fe_and_h3_clustering.rds")
strict_b <- b$strict_b  # triad_id, t1, spell_id, node1/2/3, ul12/13/23, emb12/13/23,
                         # dep_asymmetry_12/13/23 (AT t1), cap_asymmetry_*, changed*

stopifnot("spell_id" %in% names(strict_b))

# --- onset_year per spell, from the person-year panel ---
onset_lookup <- readRDS("data/person_year_risk_table.rds") %>%
  distinct(spell_id, onset_year)
# drop the pre-existing (AT t1) dep_asymmetry columns -- rebuilt fresh
# below at each timing specification, so the column names don't collide
strict_b <- strict_b %>%
  select(-dep_asymmetry_12, -dep_asymmetry_13, -dep_asymmetry_23) %>%
  left_join(onset_lookup, by = "spell_id")
cat("Events:", nrow(strict_b), " | missing onset_year:", sum(is.na(strict_b$onset_year)), "\n")
cat("Spell duration (t1 - onset_year + 1) summary:\n")
print(summary(strict_b$t1 - strict_b$onset_year + 1))

tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high, dep_asymmetry)

get_dep_at <- function(df, a, b_, year_col, suffix) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b_]]), hi = pmax(.data[[a]], .data[[b_]])) %>%
    left_join(tie_trade, by = setNames(c("year", "ccode_low", "ccode_high"), c(year_col, "lo", "hi"))) %>%
    rename(!!paste0("dep_asymmetry_", suffix) := dep_asymmetry) %>%
    select(-lo, -hi)
}

build_long <- function(df, label) {
  n_valid <- sum(!is.na(df$dep_asymmetry_12) & !is.na(df$dep_asymmetry_13) & !is.na(df$dep_asymmetry_23))
  cat(sprintf("%-30s events with all 3 ties valid: %d of %d (%.1f%%)\n",
              label, n_valid, nrow(df), 100 * n_valid / nrow(df)))
  df %>%
    filter(!is.na(dep_asymmetry_12), !is.na(dep_asymmetry_13), !is.na(dep_asymmetry_23)) %>%
    select(event_id, triad_id, primary_bloc, ul12, ul13, ul23, emb12, emb13, emb23,
           dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, changed12, changed13, changed23) %>%
    pivot_longer(
      cols = c(ul12, ul13, ul23, emb12, emb13, emb23,
               dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, changed12, changed13, changed23),
      names_to = c(".value", "tie"),
      names_pattern = "(ul|emb|dep_asymmetry|changed)_?(12|13|23)"
    ) %>%
    rename(unbal_load = ul, embeddedness = emb) %>%
    mutate(z_unbal_load = as.numeric(scale(unbal_load)),
           z_dep_asym = as.numeric(scale(dep_asymmetry)),
           z_emb = as.numeric(scale(embeddedness)))
}

fit2way <- function(data) {
  surv <- Surv(rep(1, nrow(data)), data$changed)
  m <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + z_emb + strata(event_id),
                           data = data, method = "efron"))
  cl2 <- data[, c("triad_id", "primary_bloc")]
  list(m = m,
       ct_1way = coeftest(m, vcov = vcovCL(m, cluster = data$triad_id)),
       ct_2way = coeftest(m, vcov = vcovCL(m, cluster = cl2, multi0 = TRUE)))
}

# --- (a) AT t1 (contemporaneous, reproduces the main H3 result on this subsample) ---
d_t1 <- strict_b %>%
  get_dep_at("node1", "node2", "t1", "12") %>%
  get_dep_at("node1", "node3", "t1", "13") %>%
  get_dep_at("node2", "node3", "t1", "23")
long_t1 <- build_long(d_t1, "at t1 (contemporaneous)")

# --- (b) AT spell onset (before the imbalance leading to this realignment began) ---
d_onset <- strict_b %>%
  mutate(lookup_year = onset_year) %>%
  get_dep_at("node1", "node2", "lookup_year", "12") %>%
  get_dep_at("node1", "node3", "lookup_year", "13") %>%
  get_dep_at("node2", "node3", "lookup_year", "23")
long_onset <- build_long(d_onset, "at spell onset")

# --- (c)-(e) lagged BEFORE onset: onset-1, onset-3, onset-5 ---
make_lag <- function(lag) {
  d <- strict_b %>%
    mutate(lookup_year = onset_year - lag) %>%
    get_dep_at("node1", "node2", "lookup_year", "12") %>%
    get_dep_at("node1", "node3", "lookup_year", "13") %>%
    get_dep_at("node2", "node3", "lookup_year", "23")
  build_long(d, sprintf("onset - %d years", lag))
}
long_lag1 <- make_lag(1)
long_lag3 <- make_lag(3)
long_lag5 <- make_lag(5)

r_t1     <- fit2way(long_t1)
r_onset  <- fit2way(long_onset)
r_lag1   <- fit2way(long_lag1)
r_lag3   <- fit2way(long_lag3)
r_lag5   <- fit2way(long_lag5)

sink("results/temporal_endogeneity_output.txt")
cat("=== dep_asymmetry measured AT t1 (contemporaneous with realignment) ===\n")
cat("N events:", n_distinct(long_t1$event_id), "\n")
cat("--- 1-way ---\n"); print(r_t1$ct_1way)
cat("\n--- 2-way ---\n"); print(r_t1$ct_2way)

cat("\n\n=== dep_asymmetry measured AT SPELL ONSET (before this realignment) ===\n")
cat("N events:", n_distinct(long_onset$event_id), "\n")
cat("--- 1-way ---\n"); print(r_onset$ct_1way)
cat("\n--- 2-way ---\n"); print(r_onset$ct_2way)

cat("\n\n=== dep_asymmetry measured at ONSET - 1 YEAR ===\n")
cat("N events:", n_distinct(long_lag1$event_id), "\n")
cat("--- 1-way ---\n"); print(r_lag1$ct_1way)
cat("\n--- 2-way ---\n"); print(r_lag1$ct_2way)

cat("\n\n=== dep_asymmetry measured at ONSET - 3 YEARS ===\n")
cat("N events:", n_distinct(long_lag3$event_id), "\n")
cat("--- 1-way ---\n"); print(r_lag3$ct_1way)
cat("\n--- 2-way ---\n"); print(r_lag3$ct_2way)

cat("\n\n=== dep_asymmetry measured at ONSET - 5 YEARS ===\n")
cat("N events:", n_distinct(long_lag5$event_id), "\n")
cat("--- 1-way ---\n"); print(r_lag5$ct_1way)
cat("\n--- 2-way ---\n"); print(r_lag5$ct_2way)

cat("\n\n=== SUMMARY: z_dep_asym coefficient across timing specifications ===\n")
summarize_one <- function(ct, label) {
  cat(sprintf("%-30s b=%8.4f  SE=%7.4f  p=%8.5f\n", label,
              ct["z_dep_asym","Estimate"], ct["z_dep_asym","Std. Error"], ct["z_dep_asym","Pr(>|z|)"]))
}
cat("--- 1-way clustered ---\n")
summarize_one(r_t1$ct_1way, "at t1 (contemporaneous)")
summarize_one(r_onset$ct_1way, "at spell onset")
summarize_one(r_lag1$ct_1way, "onset - 1yr")
summarize_one(r_lag3$ct_1way, "onset - 3yr")
summarize_one(r_lag5$ct_1way, "onset - 5yr")
cat("--- 2-way (triad+bloc) clustered ---\n")
summarize_one(r_t1$ct_2way, "at t1 (contemporaneous)")
summarize_one(r_onset$ct_2way, "at spell onset")
summarize_one(r_lag1$ct_2way, "onset - 1yr")
summarize_one(r_lag3$ct_2way, "onset - 3yr")
summarize_one(r_lag5$ct_2way, "onset - 5yr")
sink()

cat(readLines("results/temporal_endogeneity_output.txt"), sep = "\n")

saveRDS(list(r_t1 = r_t1, r_onset = r_onset, r_lag1 = r_lag1, r_lag3 = r_lag3, r_lag5 = r_lag5),
        "results/temporal_endogeneity.rds")
cat("\n\n11_temporal_endogeneity.R complete.\n")
