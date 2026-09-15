# ------------------------------------------------------------------
# 05: H3 (displacement) and H4 (capability). Identifies single-tie
#     realignment events and builds the unbalanced-load conditional
#     logit tie-choice model, with trade-dependence asymmetry and
#     capability asymmetry merged in per tie.
#
#     H3 claim: when the tie carrying the greatest unbalanced load is
#     ALSO highly asymmetric, it should be LESS likely (not equally
#     likely) to be the one that changes -- the trade-asymmetric tie
#     "absorbs" less of the balance-restoring pressure than its
#     structural load alone would predict, because the exposed side
#     can't afford to be the one that moves. Tested as an interaction:
#     unbal_load : dep_asymmetry should be NEGATIVE in a conditional
#     logit predicting which of a triad's 3 ties changed sign.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(purrr) })


spell_bounds <- readRDS("data/spell_bounds.rds")
tr <- readRDS("data/triads_all_years.rds") %>%
  select(year, triad_id, node1, node2, node3, sign12, sign13, sign23,
         tie12_emb, tie13_emb, tie23_emb, tie12_btw, tie13_btw, tie23_btw, balanced)

flipped <- spell_bounds %>% filter(event == "flipped")
cat("Realignment events:", nrow(flipped), "\n")

at_t1   <- tr %>% rename(sign12_t1 = sign12, sign13_t1 = sign13, sign23_t1 = sign23,
                          emb12 = tie12_emb, emb13 = tie13_emb, emb23 = tie23_emb)
at_t1p1 <- tr %>% select(year, triad_id, sign12, sign13, sign23) %>%
  rename(sign12_t2 = sign12, sign13_t2 = sign13, sign23_t2 = sign23)

events <- flipped %>%
  left_join(at_t1, by = c("triad_id", "t1" = "year")) %>%
  left_join(at_t1p1, by = c("triad_id", "lookup_year" = "year")) %>%
  mutate(
    changed12 = sign12_t1 != sign12_t2,
    changed13 = sign13_t1 != sign13_t2,
    changed23 = sign23_t1 != sign23_t2,
    n_changed = changed12 + changed13 + changed23
  )

single <- events %>% filter(n_changed == 1)
cat("Single-tie-change events:", nrow(single), "of", nrow(events),
    sprintf("(%.1f%%)\n", 100 * nrow(single) / nrow(events)))

# --- unbalanced load per (year, edge): for each edge, how many OTHER
# closed triads sharing it are themselves unbalanced that year ---
edge_long <- bind_rows(
  tr %>% transmute(year, triad_id, balanced, e_lo = pmin(node1, node2), e_hi = pmax(node1, node2)),
  tr %>% transmute(year, triad_id, balanced, e_lo = pmin(node1, node3), e_hi = pmax(node1, node3)),
  tr %>% transmute(year, triad_id, balanced, e_lo = pmin(node2, node3), e_hi = pmax(node2, node3))
)
edge_summary <- edge_long %>%
  group_by(year, e_lo, e_hi) %>%
  summarise(n_unbalanced = sum(!balanced), .groups = "drop")

get_unbal_load <- function(df, a, b, newcol) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(edge_summary, by = c("t1" = "year", "lo" = "e_lo", "hi" = "e_hi")) %>%
    rename(!!newcol := n_unbalanced) %>%
    mutate(!!newcol := !!sym(newcol) - 1) %>%
    select(-lo, -hi)
}

single <- single %>%
  get_unbal_load("node1", "node2", "ul12") %>%
  get_unbal_load("node1", "node3", "ul13") %>%
  get_unbal_load("node2", "node3", "ul23") %>%
  mutate(across(c(ul12, ul13, ul23), ~ replace_na(.x, 0)))

# --- trade-dependence asymmetry + capability asymmetry per tie, AT t1 ---
tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high, dep_asymmetry, cap_asymmetry, exposed_is_weaker, log_trade_volume)

get_trade_var <- function(df, a, b, suffix) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(tie_trade, by = c("t1" = "year", "lo" = "ccode_low", "hi" = "ccode_high")) %>%
    rename_with(~ paste0(., "_", suffix), c(dep_asymmetry, cap_asymmetry, exposed_is_weaker, log_trade_volume)) %>%
    select(-lo, -hi)
}

single <- single %>%
  get_trade_var("node1", "node2", "12") %>%
  get_trade_var("node1", "node3", "13") %>%
  get_trade_var("node2", "node3", "23")

single <- single %>%
  mutate(n_trade_valid = rowSums(!is.na(cbind(dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23))))
cat("\nSingle-tie-change events with all 3 ties having valid trade data:",
    sum(single$n_trade_valid == 3), "of", nrow(single),
    sprintf("(%.1f%%)\n", 100 * mean(single$n_trade_valid == 3)))

strict <- single %>% filter(n_trade_valid == 3)

# --- long format: one row per (event, tie) ---
long <- strict %>%
  mutate(event_id = row_number()) %>%
  select(event_id, ul12, ul13, ul23, emb12, emb13, emb23,
         dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23,
         cap_asymmetry_12, cap_asymmetry_13, cap_asymmetry_23,
         exposed_is_weaker_12, exposed_is_weaker_13, exposed_is_weaker_23,
         log_trade_volume_12, log_trade_volume_13, log_trade_volume_23,
         changed12, changed13, changed23) %>%
  pivot_longer(
    cols = -event_id,
    names_to = c(".value", "tie"),
    names_pattern = "(ul|emb|dep_asymmetry|cap_asymmetry|exposed_is_weaker|log_trade_volume|changed)_?(12|13|23)"
  ) %>%
  rename(unbal_load = ul, embeddedness = emb) %>%
  mutate(z_unbal_load = as.numeric(scale(unbal_load)),
         z_dep_asym = as.numeric(scale(dep_asymmetry)),
         z_cap_asym = as.numeric(scale(cap_asymmetry)),
         z_volume = as.numeric(scale(log_trade_volume)),
         z_emb = as.numeric(scale(embeddedness)))

cat("\nLong-format rows (3 per event):", nrow(long), " | events:", n_distinct(long$event_id), "\n")

# --- H3: does the changed tie have the max unbalanced load, on this
#     trade-valid subsample? ---
clog_baseline <- clogit(changed ~ z_unbal_load + strata(event_id), data = long)

# --- does dep_asymmetry alone predict a LOWER chance of being the changed tie? ---
clog_asym <- clogit(changed ~ z_unbal_load + z_dep_asym + strata(event_id), data = long)

# --- H3 PRIMARY: interaction. A negative unbal_load:dep_asym coefficient
#     means high asymmetry suppresses the normal "high load -> changes"
#     relationship (displacement) ---
clog_interact <- clogit(changed ~ z_unbal_load * z_dep_asym + z_emb + strata(event_id), data = long)

# --- H4: does capability asymmetry do independent work, or moderate H3? ---
clog_h4_add <- clogit(changed ~ z_unbal_load * z_dep_asym + z_emb + z_cap_asym + strata(event_id), data = long)
clog_h4_int <- clogit(changed ~ z_unbal_load * z_dep_asym * z_cap_asym + z_emb + strata(event_id), data = long)

# --- volume check: does the (reversed) positive dep_asymmetry effect survive
#     controlling for the ABSOLUTE size of the tie's trade relationship?
#     Competing explanation to rule out: maybe it's not asymmetry that makes
#     a tie more likely to be the one that changes, but simply that low-
#     volume ties (which also tend to be low-asymmetry, see correlation
#     check) are "cheaper" to flip regardless of asymmetry. ---
clog_volume <- clogit(changed ~ z_unbal_load + z_dep_asym + z_volume + z_emb + strata(event_id), data = long)

sink("results/h3_displacement_output.txt")
cat("Realignment events:", nrow(events), " | single-tie-change:", nrow(single),
    sprintf(" (%.1f%%)\n", 100*nrow(single)/nrow(events)))
cat("Single-tie-change events with all 3 ties trade-valid:", nrow(strict),
    sprintf(" (%.1f%% of single-tie-change events)\n", 100*nrow(strict)/nrow(single)))
cat("Long-format rows:", nrow(long), " events:", n_distinct(long$event_id), "\n\n")

cat("=== BASELINE: does the changed tie have the highest unbalanced load? ===\n")
print(summary(clog_baseline))

cat("\n\n=== does dep_asymmetry alone predict which tie changes? ===\n")
print(summary(clog_asym))

cat("\n\n=== H3 PRIMARY: unbal_load x dep_asymmetry interaction ===\n")
print(summary(clog_interact))

cat("\n\n=== H4: + capability asymmetry (additive) ===\n")
print(summary(clog_h4_add))

cat("\n\n=== H4: full 3-way interaction (unbal_load x dep_asym x cap_asym) ===\n")
print(summary(clog_h4_int))

cat("\n\n=== volume check: does dep_asym's effect survive controlling for absolute trade volume? ===\n")
cat("Correlation, dep_asymmetry vs. log_trade_volume (within long-format tie rows):",
    round(cor(long$dep_asymmetry, long$log_trade_volume, use = "complete.obs"), 3), "\n")
print(summary(clog_volume))
sink()

cat(readLines("results/h3_displacement_output.txt"), sep = "\n")

saveRDS(list(single = single, strict = strict, long = long,
             clog_baseline = clog_baseline, clog_asym = clog_asym,
             clog_interact = clog_interact, clog_h4_add = clog_h4_add, clog_h4_int = clog_h4_int,
             clog_volume = clog_volume),
        "results/h3_displacement.rds")
cat("\n\n05_h3_displacement.R complete. See results/h3_displacement_output.txt\n")
