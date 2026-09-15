# ------------------------------------------------------------------
# 03: Build the person-year competing-risks panel (spell construction,
#     cause-specific events, duration bins), with per-tie trade-
#     dependence ASYMMETRY and capability asymmetry merged in as
#     TIME-VARYING covariates (updated every year of the spell, not
#     frozen at onset).
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(tidyr)
})


tr_full <- readRDS("data/triads_all_years.rds")
tr <- tr_full %>%
  select(year, triad_id, node1, node2, node3, balanced, actor_load_mean,
         tie_btw_mean, tie_emb_mean, triad_log_cinc_mean) %>%
  arrange(triad_id, year)

LAST_YEAR <- max(tr$year)

# --- identify spells ---
tr <- tr %>%
  group_by(triad_id) %>%
  mutate(
    prev_year = lag(year), prev_balanced = lag(balanced),
    is_onset = !balanced & (is.na(prev_year) | prev_year != year - 1 | prev_balanced),
    spell_id_local = cumsum(is_onset)
  ) %>%
  ungroup() %>%
  mutate(spell_id = paste(triad_id, spell_id_local, sep = "_S"))

spells <- tr %>% filter(!balanced)
cat("Person-year candidate rows (all unbalanced triad-years):", nrow(spells), "\n")
cat("Unique spells:", n_distinct(spells$spell_id), "\n")

spell_bounds <- spells %>% group_by(spell_id, triad_id) %>% summarise(t1 = max(year), .groups = "drop")
onset <- spells %>% group_by(spell_id) %>% summarise(onset_year = min(year), .groups = "drop")

next_status <- tr %>% select(year, triad_id, balanced) %>% rename(next_year = year, balanced_next = balanced)

spell_bounds <- spell_bounds %>%
  mutate(lookup_year = t1 + 1) %>%
  left_join(next_status, by = c("triad_id" = "triad_id", "lookup_year" = "next_year")) %>%
  mutate(
    event = case_when(
      t1 >= LAST_YEAR                 ~ "censored_end_of_panel",
      is.na(balanced_next)            ~ "dissolved",
      balanced_next                   ~ "flipped",
      TRUE                            ~ "ERROR_still_unbalanced"
    )
  )
stopifnot(!any(spell_bounds$event == "ERROR_still_unbalanced"))
cat("\nSpell terminal events:\n"); print(table(spell_bounds$event))

person_years <- spells %>%
  left_join(spell_bounds %>% select(spell_id, t1, event), by = "spell_id") %>%
  left_join(onset, by = "spell_id") %>%
  mutate(
    is_terminal_row = year == t1,
    row_event = ifelse(is_terminal_row, event, "continues"),
    duration = year - onset_year + 1,
    dur_bin = cut(duration, breaks = c(0, 1, 2, 3, 5, 10, Inf),
                  labels = c("yr1", "yr2", "yr3", "yr4-5", "yr6-10", "yr11+"))
  ) %>%
  filter(row_event != "censored_end_of_panel")

cat("\nUsable person-year rows (transition to Y+1 observed):", nrow(person_years), "\n")
print(table(person_years$row_event))

# ------------------------------------------------------------------
# Merge time-varying trade-dependence asymmetry + capability asymmetry,
# per tie, at each (year, triad) row -- not frozen at spell onset.
# ------------------------------------------------------------------
tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high, trade_dependence, dep_exposed, dep_asymmetry,
         log_dep_ratio, cap_asymmetry, exposed_is_weaker, log_trade_volume)

get_tie_var <- function(df, a, b, cols, suffix) {
  df %>%
    mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(tie_trade, by = c("year" = "year", "lo" = "ccode_low", "hi" = "ccode_high")) %>%
    rename_with(~ paste0(., "_", suffix), all_of(cols)) %>%
    select(-lo, -hi)
}

trade_cols <- c("trade_dependence", "dep_exposed", "dep_asymmetry", "log_dep_ratio",
                 "cap_asymmetry", "exposed_is_weaker", "log_trade_volume")

person_years <- person_years %>%
  get_tie_var("node1", "node2", trade_cols, "12") %>%
  get_tie_var("node1", "node3", trade_cols, "13") %>%
  get_tie_var("node2", "node3", trade_cols, "23")

# --- triad-year summaries across the 3 ties ---
person_years <- person_years %>%
  rowwise() %>%
  mutate(
    trade_dependence_mean = mean(c(trade_dependence_12, trade_dependence_13, trade_dependence_23), na.rm = TRUE),
    dep_asym_mean  = mean(c(dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23), na.rm = TRUE),
    dep_asym_max   = max(c(dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23), na.rm = TRUE),
    cap_asym_mean  = mean(c(cap_asymmetry_12, cap_asymmetry_13, cap_asymmetry_23), na.rm = TRUE),
    cap_asym_max   = max(c(cap_asymmetry_12, cap_asymmetry_13, cap_asymmetry_23), na.rm = TRUE),
    # absolute size of the dyadic relationship, independent of trade_share --
    # controls for whether dep_asymmetry is genuine entanglement or a
    # small-denominator artifact
    trade_volume_mean = mean(c(log_trade_volume_12, log_trade_volume_13, log_trade_volume_23), na.rm = TRUE),
    trade_volume_max  = max(c(log_trade_volume_12, log_trade_volume_13, log_trade_volume_23), na.rm = TRUE),
    n_valid_ties   = sum(!is.na(c(dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23))),
    any_exposed_weaker = any(c(exposed_is_weaker_12, exposed_is_weaker_13, exposed_is_weaker_23), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(across(c(dep_asym_mean, dep_asym_max, cap_asym_mean, cap_asym_max, trade_volume_mean, trade_volume_max),
                ~ ifelse(is.infinite(.) | is.nan(.), NA, .)))

cat("\nTrade-asymmetry coverage in person-year panel:\n")
cat("  rows with >=1 valid tie:", sum(person_years$n_valid_ties >= 1),
    sprintf(" (%.1f%%)\n", 100*mean(person_years$n_valid_ties >= 1)))
cat("  rows with all 3 ties valid:", sum(person_years$n_valid_ties == 3),
    sprintf(" (%.1f%%)\n", 100*mean(person_years$n_valid_ties == 3)))

zscore <- function(x) as.numeric(scale(x))
person_years <- person_years %>% mutate(
  z_actor_load = zscore(actor_load_mean),
  z_tie_btw    = zscore(tie_btw_mean),
  z_tie_emb    = zscore(tie_emb_mean),
  z_cinc       = zscore(triad_log_cinc_mean),
  era = cut(year, breaks = c(1815, 1945, 1989, 2015), labels = c("pre-1945", "1945-1989", "post-1989")),
  event_dissolve = as.numeric(row_event == "dissolved"),
  event_flip     = as.numeric(row_event == "flipped")
)

write_csv(person_years, "data/person_year_risk_table.csv")
saveRDS(person_years, "data/person_year_risk_table.rds")
saveRDS(spell_bounds, "data/spell_bounds.rds")

cat("\n03_person_year_panel.R complete.\n")
