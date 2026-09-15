# ------------------------------------------------------------------
# 20: Two checks on the M4 specification.
#     (1) add a control for how long a candidate relationship has
#     already held its current sign to M4, since alliances and MIDs
#     have very different typical durations and duration could
#     correlate with dependence asymmetry independently of starting
#     sign itself;
#     (2) rerun the temporal-endogeneity (pre-spell) specifications
#     using the M4 covariate set throughout, not the M1 baseline.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest) })


# ==================================================================
# Build "years in current sign, ending at year Y" for every (dyad, year)
# ==================================================================
st <- read_csv("data/signed_ties.csv", show_col_types = FALSE) %>%
  filter(year >= 1816, year <= 2012) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  arrange(ccode_low, ccode_high, year)

st <- st %>% group_by(ccode_low, ccode_high) %>%
  mutate(prev_year = lag(year), prev_sign = lag(sign),
         is_run_start = is.na(prev_year) | prev_year != year - 1 | prev_sign != sign,
         run_id = cumsum(is_run_start)) %>%
  group_by(ccode_low, ccode_high, run_id) %>%
  mutate(duration_in_sign = row_number()) %>%
  ungroup() %>%
  select(year, ccode_low, ccode_high, duration_in_sign)

get_duration <- function(df, a, b_, year_col, newcol) {
  df %>% mutate(lo = as.character(pmin(as.numeric(.data[[a]]), as.numeric(.data[[b_]]))),
                hi = as.character(pmax(as.numeric(.data[[a]]), as.numeric(.data[[b_]])))) %>%
    left_join(st, by = setNames(c("year","ccode_low","ccode_high"), c(year_col,"lo","hi"))) %>%
    rename(!!newcol := duration_in_sign) %>% select(-lo,-hi)
}

# ==================================================================
# PART A: add duration to M4, at t1 (contemporaneous)
# ==================================================================
h3 <- readRDS("results/h3_displacement.rds")
strict <- h3$strict %>% mutate(event_id = row_number())
tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high, trade_share_low, trade_share_high)
get_shares <- function(df, a, b_, suffix) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b_]]), hi = pmax(.data[[a]], .data[[b_]])) %>%
    left_join(tie_trade, by = c("t1"="year","lo"="ccode_low","hi"="ccode_high")) %>%
    rename_with(~ paste0(., "_", suffix), c(trade_share_low, trade_share_high)) %>% select(-lo,-hi)
}
strict <- strict %>% get_shares("node1","node2","12") %>% get_shares("node1","node3","13") %>% get_shares("node2","node3","23") %>%
  mutate(mutual_min_12 = pmin(trade_share_low_12, trade_share_high_12),
         mutual_min_13 = pmin(trade_share_low_13, trade_share_high_13),
         mutual_min_23 = pmin(trade_share_low_23, trade_share_high_23)) %>%
  get_duration("node1","node2","t1","dur12") %>%
  get_duration("node1","node3","t1","dur13") %>%
  get_duration("node2","node3","t1","dur23")

cat("Duration-in-sign summary (t1, pooled across 3 ties):\n")
print(summary(c(strict$dur12, strict$dur13, strict$dur23)))

long <- strict %>%
  select(event_id, triad_id, sign12_t1, sign13_t1, sign23_t1,
         ul12, ul13, ul23, emb12, emb13, emb23, dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23,
         cap_asymmetry_12, cap_asymmetry_13, cap_asymmetry_23,
         mutual_min_12, mutual_min_13, mutual_min_23, dur12, dur13, dur23, changed12, changed13, changed23) %>%
  pivot_longer(cols = c(ul12,ul13,ul23,emb12,emb13,emb23,dep_asymmetry_12,dep_asymmetry_13,dep_asymmetry_23,
                        cap_asymmetry_12,cap_asymmetry_13,cap_asymmetry_23,
                        mutual_min_12,mutual_min_13,mutual_min_23,dur12,dur13,dur23,changed12,changed13,changed23),
               names_to = c(".value","tie"),
               names_pattern = "(ul|emb|dep_asymmetry|cap_asymmetry|mutual_min|dur|changed)_?(12|13|23)") %>%
  left_join(strict %>% select(event_id, sign12_t1, sign13_t1, sign23_t1) %>%
              pivot_longer(cols = c(sign12_t1,sign13_t1,sign23_t1), names_to = "tie2", values_to = "sign_t1") %>%
              mutate(tie2 = gsub("sign|_t1","",tie2)),
            by = c("event_id","tie" = "tie2")) %>%
  mutate(initial_negative = as.numeric(sign_t1 == -1),
         z_unbal_load = as.numeric(scale(ul)), z_dep_asym = as.numeric(scale(dep_asymmetry)),
         z_emb = as.numeric(scale(emb)), z_cap_asym = as.numeric(scale(cap_asymmetry)),
         z_mutual_min = as.numeric(scale(mutual_min)),
         z_log_dur = as.numeric(scale(log(pmax(dur, 1)))))

long_m4dur <- long %>% filter(!is.na(z_dep_asym), !is.na(z_cap_asym), !is.na(z_mutual_min), !is.na(z_log_dur))
cat("\nM4 + duration sample: N events =", n_distinct(long_m4dur$event_id), "\n")

surv <- Surv(rep(1, nrow(long_m4dur)), long_m4dur$changed)
m4 <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + initial_negative + z_mutual_min + z_cap_asym + z_emb + strata(event_id), data = long_m4dur, method = "efron"))
m4dur <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + initial_negative + z_mutual_min + z_cap_asym + z_log_dur + z_emb + strata(event_id), data = long_m4dur, method = "efron"))
ct_m4 <- coeftest(m4, vcov = vcovCL(m4, cluster = long_m4dur$triad_id))
ct_m4dur <- coeftest(m4dur, vcov = vcovCL(m4dur, cluster = long_m4dur$triad_id))

# two-way (triad+bloc) clustered SE for M4 (for Table 1 completeness)
h3b <- readRDS("results/bloc_fe_and_h3_clustering.rds")
bloc_lookup <- h3b$strict_b %>% mutate(event_id_b = row_number()) %>% select(event_id_b, primary_bloc, t1, triad_id)
# match by triad_id+t1 (unique per event in both objects, built from the same underlying `strict`)
strict_key <- strict %>% select(event_id, triad_id, t1)
bloc_key <- bloc_lookup %>% select(triad_id, t1, primary_bloc) %>% distinct()
event_bloc <- strict_key %>% left_join(bloc_key, by = c("triad_id","t1"))
long_m4_bloc <- long %>% filter(!is.na(z_dep_asym), !is.na(z_cap_asym), !is.na(z_mutual_min)) %>%
  left_join(event_bloc %>% select(event_id, primary_bloc), by = "event_id") %>%
  mutate(primary_bloc = ifelse(is.na(primary_bloc), "none", primary_bloc))
m4_forbloc <- do.call(coxph, list(formula = Surv(rep(1,nrow(long_m4_bloc)), long_m4_bloc$changed) ~
                                     z_unbal_load + z_dep_asym + initial_negative + z_mutual_min + z_cap_asym + z_emb + strata(event_id),
                                   data = long_m4_bloc, method = "efron"))
cl2 <- long_m4_bloc[, c("triad_id","primary_bloc")]
ct_m4_2way <- coeftest(m4_forbloc, vcov = vcovCL(m4_forbloc, cluster = cl2, multi0 = TRUE))

sink("results/duration_control_and_m4_consistency_output.txt")
cat("=== PART A: duration-in-current-sign control added to M4 ===\n")
cat("Duration (years, t1) summary:\n"); print(summary(c(strict$dur12, strict$dur13, strict$dur23)))
cat("\n--- M4 (no duration), matched sample ---\n"); print(ct_m4)
cat("\n--- M4 + duration control ---\n"); print(ct_m4dur)
cat("\n--- M4, two-way (triad+bloc) clustered SE ---\n"); print(ct_m4_2way)
sink()
cat(readLines("results/duration_control_and_m4_consistency_output.txt"), sep = "\n")

saveRDS(list(m4 = m4, m4dur = m4dur, ct_m4 = ct_m4, ct_m4dur = ct_m4dur, ct_m4_2way = ct_m4_2way),
        "results/duration_control_part_a.rds")
cat("\n\nPart A complete.\n")

# ==================================================================
# PART B: re-run the temporal-endogeneity (pre-spell) checks using the
# M4 covariate set (starting sign, weak-link dependence, capability),
# not just the M1 baseline used in script 11.
# ==================================================================
onset_lookup <- readRDS("data/person_year_risk_table.rds") %>% distinct(spell_id, onset_year)
# drop the pre-existing (AT t1) versions of these columns -- Part B re-merges
# them fresh at each lag timing, so the names would otherwise collide
strict_b <- strict %>%
  select(-dep_asymmetry_12, -dep_asymmetry_13, -dep_asymmetry_23,
         -cap_asymmetry_12, -cap_asymmetry_13, -cap_asymmetry_23,
         -mutual_min_12, -mutual_min_13, -mutual_min_23,
         -dur12, -dur13, -dur23) %>%
  left_join(onset_lookup, by = "spell_id")

tie_trade_full <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high, dep_asymmetry, cap_asymmetry, trade_share_low, trade_share_high)

get_m4vars_at <- function(df, a, b_, year_col, suffix) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b_]]), hi = pmax(.data[[a]], .data[[b_]])) %>%
    left_join(tie_trade_full, by = setNames(c("year","ccode_low","ccode_high"), c(year_col,"lo","hi"))) %>%
    mutate(mutual_min = pmin(trade_share_low, trade_share_high)) %>%
    rename_with(~ paste0(., "_", suffix), c(dep_asymmetry, cap_asymmetry, mutual_min)) %>%
    select(-lo, -hi, -trade_share_low, -trade_share_high)
}

bloc_lookup <- h3b$strict_b %>% select(triad_id, t1, primary_bloc) %>% distinct()
event_bloc <- strict %>% select(event_id, triad_id, t1) %>% left_join(bloc_lookup, by = c("triad_id","t1")) %>%
  mutate(primary_bloc = ifelse(is.na(primary_bloc), "none", primary_bloc))
strict_b <- strict_b %>% left_join(event_bloc %>% select(event_id, primary_bloc), by = "event_id")

build_m4_at <- function(lag_years, label) {
  d <- strict_b %>% mutate(lookup_year = onset_year - lag_years) %>%
    get_m4vars_at("node1","node2","lookup_year","12") %>%
    get_m4vars_at("node1","node3","lookup_year","13") %>%
    get_m4vars_at("node2","node3","lookup_year","23") %>%
    get_duration("node1","node2","lookup_year","dur12") %>%
    get_duration("node1","node3","lookup_year","dur13") %>%
    get_duration("node2","node3","lookup_year","dur23")
  d <- d %>% filter(!is.na(dep_asymmetry_12), !is.na(dep_asymmetry_13), !is.na(dep_asymmetry_23),
                     !is.na(cap_asymmetry_12), !is.na(cap_asymmetry_13), !is.na(cap_asymmetry_23),
                     !is.na(mutual_min_12), !is.na(mutual_min_13), !is.na(mutual_min_23))
  long_d <- d %>% select(event_id, triad_id, primary_bloc, sign12_t1, sign13_t1, sign23_t1, ul12,ul13,ul23,emb12,emb13,emb23,
                          dep_asymmetry_12,dep_asymmetry_13,dep_asymmetry_23,cap_asymmetry_12,cap_asymmetry_13,cap_asymmetry_23,
                          mutual_min_12,mutual_min_13,mutual_min_23,changed12,changed13,changed23) %>%
    pivot_longer(cols = c(ul12,ul13,ul23,emb12,emb13,emb23,dep_asymmetry_12,dep_asymmetry_13,dep_asymmetry_23,
                          cap_asymmetry_12,cap_asymmetry_13,cap_asymmetry_23,mutual_min_12,mutual_min_13,mutual_min_23,
                          changed12,changed13,changed23),
                 names_to = c(".value","tie"),
                 names_pattern = "(ul|emb|dep_asymmetry|cap_asymmetry|mutual_min|changed)_?(12|13|23)") %>%
    left_join(d %>% select(event_id, sign12_t1, sign13_t1, sign23_t1) %>%
                pivot_longer(cols = c(sign12_t1,sign13_t1,sign23_t1), names_to = "tie2", values_to = "sign_t1") %>%
                mutate(tie2 = gsub("sign|_t1","",tie2)),
              by = c("event_id","tie" = "tie2")) %>%
    mutate(initial_negative = as.numeric(sign_t1 == -1),
           z_unbal_load = as.numeric(scale(ul)), z_dep_asym = as.numeric(scale(dep_asymmetry)),
           z_emb = as.numeric(scale(emb)), z_cap_asym = as.numeric(scale(cap_asymmetry)),
           z_mutual_min = as.numeric(scale(mutual_min)))
  surv_d <- Surv(rep(1, nrow(long_d)), long_d$changed)
  m <- do.call(coxph, list(formula = surv_d ~ z_unbal_load + z_dep_asym + initial_negative + z_mutual_min + z_cap_asym + z_emb + strata(event_id),
                            data = long_d, method = "efron"))
  ct1 <- coeftest(m, vcov = vcovCL(m, cluster = long_d$triad_id))
  cl2 <- long_d[, c("triad_id","primary_bloc")]
  ct2 <- coeftest(m, vcov = vcovCL(m, cluster = cl2, multi0 = TRUE))
  cat(sprintf("%-20s N events=%5d  b=%8.4f  p(triad)=%8.5f  p(triad+bloc)=%8.5f\n", label, n_distinct(long_d$event_id),
              ct1["z_dep_asym","Estimate"], ct1["z_dep_asym","Pr(>|z|)"], ct2["z_dep_asym","Pr(>|z|)"]))
  invisible(list(ct1 = ct1, ct2 = ct2))
}

sink("results/m4_temporal_endogeneity_output.txt")
cat("=== M4 specification at each measurement timing, 1-way and 2-way (triad+bloc) clustered ===\n")
r0 <- build_m4_at(0, "Spell onset")
r1 <- build_m4_at(1, "Onset - 1yr")
r3 <- build_m4_at(3, "Onset - 3yr")
r5 <- build_m4_at(5, "Onset - 5yr")
sink()
cat(readLines("results/m4_temporal_endogeneity_output.txt"), sep = "\n")
saveRDS(list(r0=r0,r1=r1,r3=r3,r5=r5), "results/m4_temporal_endogeneity.rds")
cat("\n\n20_duration_control_and_m4_consistency.R complete.\n")
