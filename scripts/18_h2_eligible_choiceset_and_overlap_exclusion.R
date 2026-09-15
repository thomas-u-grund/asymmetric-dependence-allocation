# ------------------------------------------------------------------
# 18: Two choice-set corrections.
#
#     (1) H2's choice set has the same problem H1's originally did: in
#     a one-negative triad, the negative (MID) relationship cannot
#     possibly be "the alliance that terminated." Re-estimate H2 with
#     the choice set restricted to relationships that were actually
#     alliances (positive) at t1 -- for one-negative triads this is a
#     clean two-way choice between the two eligible alliances.
#
#     (2) Exclude every single-tie-change event where the CHANGED tie's
#     dyad-year, at t1, had an alliance and a MID coexisting (coded
#     negative under this paper's convention, but a neg_to_pos
#     transition there can just mean the MID ended while the alliance
#     continued throughout -- not two hostile states newly allying).
#     Re-estimate H1 on the events that remain.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest) })


# ==================================================================
# PART A: H2, eligible-alliance choice set only
# ==================================================================
spell_bounds <- readRDS("data/spell_bounds.rds")
tr <- readRDS("data/triads_all_years.rds") %>%
  select(year, triad_id, node1, node2, node3, sign12, sign13, sign23, tie12_emb, tie13_emb, tie23_emb)
dissolved <- spell_bounds %>% filter(event == "dissolved")
at_t1 <- tr %>% rename(sign12_t1 = sign12, sign13_t1 = sign13, sign23_t1 = sign23,
                        emb12 = tie12_emb, emb13 = tie13_emb, emb23 = tie23_emb)
diss_events <- dissolved %>% left_join(at_t1, by = c("triad_id", "t1" = "year"))

signed_ties <- read_csv("data/signed_ties.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high))
tie_present_at <- function(df, a, b_, year_col, newcol) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b_]]), hi = pmax(.data[[a]], .data[[b_]])) %>%
    left_join(signed_ties %>% select(year, ccode_low, ccode_high) %>% mutate(present = TRUE),
              by = setNames(c("year","ccode_low","ccode_high"), c(year_col,"lo","hi"))) %>%
    rename(!!newcol := present) %>% mutate(!!newcol := !is.na(!!sym(newcol))) %>% select(-lo,-hi)
}
diss_events <- diss_events %>% mutate(lookup_year = t1 + 1) %>%
  tie_present_at("node1","node2","lookup_year","present12") %>%
  tie_present_at("node1","node3","lookup_year","present13") %>%
  tie_present_at("node2","node3","lookup_year","present23") %>%
  mutate(gone12 = !present12, gone13 = !present13, gone23 = !present23, n_gone = gone12+gone13+gone23)
diss_single <- diss_events %>% filter(n_gone == 1) %>%
  mutate(gone_sign = case_when(gone12 ~ sign12_t1, gone13 ~ sign13_t1, gone23 ~ sign23_t1))

tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high, dep_asymmetry)
get_tv <- function(df, a, b_, suffix) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b_]]), hi = pmax(.data[[a]], .data[[b_]])) %>%
    left_join(tie_trade, by = c("t1"="year","lo"="ccode_low","hi"="ccode_high")) %>%
    rename(!!paste0("dep_asymmetry_", suffix) := dep_asymmetry) %>% select(-lo,-hi)
}
diss_single <- diss_single %>% get_tv("node1","node2","12") %>% get_tv("node1","node3","13") %>% get_tv("node2","node3","23")

# restrict to alliance-termination events (the disappearing tie was positive)
alliance_end <- diss_single %>% filter(gone_sign == 1)
cat("Alliance-termination events:", nrow(alliance_end), "\n")
cat("Sign composition at t1 (how many of the 3 candidate ties were alliances):\n")
alliance_end <- alliance_end %>% mutate(n_alliances_t1 = (sign12_t1==1)+(sign13_t1==1)+(sign23_t1==1))
print(table(alliance_end$n_alliances_t1))

# eligible choice set: only ties that were alliances (positive) at t1 can be
# "the alliance that terminated" -- drop the negative (MID) tie, if any
build_eligible_long <- function(df) {
  df %>% mutate(event_id = row_number()) %>%
    select(event_id, triad_id, sign12_t1, sign13_t1, sign23_t1, emb12, emb13, emb23,
           dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, gone12, gone13, gone23) %>%
    pivot_longer(cols = c(emb12, emb13, emb23, dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, gone12, gone13, gone23),
                 names_to = c(".value", "tie"), names_pattern = "(emb|dep_asymmetry|gone)_?(12|13|23)") %>%
    left_join(df %>% mutate(event_id = row_number()) %>% select(event_id, sign12_t1, sign13_t1, sign23_t1) %>%
                pivot_longer(cols = c(sign12_t1, sign13_t1, sign23_t1), names_to = "tie2", values_to = "sign_t1") %>%
                mutate(tie2 = gsub("sign|_t1", "", tie2)),
              by = c("event_id", "tie" = "tie2")) %>%
    filter(sign_t1 == 1) %>%  # ELIGIBLE CHOICE SET: only alliances at t1
    mutate(z_dep_asym = as.numeric(scale(dep_asymmetry)), z_emb = as.numeric(scale(emb)))
}

long_elig <- build_eligible_long(alliance_end %>% filter(!is.na(dep_asymmetry_12), !is.na(dep_asymmetry_13), !is.na(dep_asymmetry_23)))
cat("\nLong-format rows in eligible-choice-set model:", nrow(long_elig), " | events:", n_distinct(long_elig$event_id), "\n")
cat("Alternatives per event (should be 2 for one-negative triads, 3 for three-alliance triads):\n")
print(table(table(long_elig$event_id)))

surv_e <- Surv(rep(1, nrow(long_elig)), long_elig$gone)
m_elig <- do.call(coxph, list(formula = surv_e ~ z_dep_asym + z_emb + strata(event_id), data = long_elig, method = "efron"))
ct_elig <- coeftest(m_elig, vcov = vcovCL(m_elig, cluster = long_elig$triad_id))

# compare: original pooled-choice-set model (all 3 ties eligible) on the SAME alliance-termination events
long_pooled <- alliance_end %>% filter(!is.na(dep_asymmetry_12), !is.na(dep_asymmetry_13), !is.na(dep_asymmetry_23)) %>%
  mutate(event_id = row_number()) %>%
  select(event_id, triad_id, emb12, emb13, emb23, dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, gone12, gone13, gone23) %>%
  pivot_longer(cols = c(emb12, emb13, emb23, dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, gone12, gone13, gone23),
               names_to = c(".value", "tie"), names_pattern = "(emb|dep_asymmetry|gone)_?(12|13|23)") %>%
  mutate(z_dep_asym = as.numeric(scale(dep_asymmetry)), z_emb = as.numeric(scale(emb)))
surv_p <- Surv(rep(1, nrow(long_pooled)), long_pooled$gone)
m_pooled <- do.call(coxph, list(formula = surv_p ~ z_dep_asym + z_emb + strata(event_id), data = long_pooled, method = "efron"))
ct_pooled <- coeftest(m_pooled, vcov = vcovCL(m_pooled, cluster = long_pooled$triad_id))

# ==================================================================
# PART B: exclude alliance+MID overlap events, re-estimate H1
# ==================================================================
h3 <- readRDS("results/h3_displacement.rds")
strict <- h3$strict

alliances_only <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(pmin(ccode1, ccode2)), ccode_high = as.character(pmax(ccode1, ccode2))) %>%
  filter(defense == 1 | neutrality == 1 | nonaggression == 1 | entente == 1) %>%
  select(year, ccode_low, ccode_high) %>% distinct() %>% mutate(has_alliance = TRUE)
mids_only <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(pmin(statea, stateb)), ccode_high = as.character(pmax(statea, stateb))) %>%
  filter(hihost >= 3) %>% select(year, ccode_low, ccode_high) %>% distinct() %>% mutate(has_mid = TRUE)
overlap <- alliances_only %>% inner_join(mids_only, by = c("year","ccode_low","ccode_high")) %>%
  mutate(overlap_key = paste(year, ccode_low, ccode_high))

is_overlap_tie <- function(df, a, b_, year_col) {
  # node1/node2/node3 are CHARACTER ccodes -- pmin/pmax on them would compare
  # lexicographically (e.g. "9" > "10"), silently swapping low/high relative
  # to the numeric convention used everywhere else in this pipeline. Convert
  # to numeric before pmin/pmax, then back to character to match the
  # reference tables' key format.
  df %>% mutate(lo = as.character(pmin(as.numeric(.data[[a]]), as.numeric(.data[[b_]]))),
                hi = as.character(pmax(as.numeric(.data[[a]]), as.numeric(.data[[b_]])))) %>%
    mutate(key = paste(.data[[year_col]], lo, hi)) %>%
    mutate(is_overlap = key %in% overlap$overlap_key) %>% pull(is_overlap)
}
# compute per-tie overlap flags as ordinary columns FIRST, then select among
# them with case_when -- calling is_overlap_tie(cur_data(), ...) directly
# inside case_when's branches silently mis-evaluates cur_data()'s scope and
# was caught producing an implausible ~96% "overlap" rate; verified against
# a direct call that the per-tie rate for tie12 alone is a sane 37.5%
strict <- strict %>% mutate(
  overlap12 = is_overlap_tie(., "node1", "node2", "t1"),
  overlap13 = is_overlap_tie(., "node1", "node3", "t1"),
  overlap23 = is_overlap_tie(., "node2", "node3", "t1"),
  changed_tie_is_overlap = case_when(
    changed12 ~ overlap12,
    changed13 ~ overlap13,
    changed23 ~ overlap23
  ))
cat("\n\nEvents where the CHANGED tie was an alliance+MID overlap dyad at t1:",
    sum(strict$changed_tie_is_overlap, na.rm = TRUE), "of", nrow(strict), "\n")

strict_no_overlap <- strict %>% filter(!changed_tie_is_overlap | is.na(changed_tie_is_overlap))
long_no_overlap <- strict_no_overlap %>% mutate(event_id = row_number()) %>%
  select(event_id, triad_id, ul12, ul13, ul23, emb12, emb13, emb23,
         dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, changed12, changed13, changed23) %>%
  pivot_longer(cols = c(ul12,ul13,ul23,emb12,emb13,emb23,dep_asymmetry_12,dep_asymmetry_13,dep_asymmetry_23,changed12,changed13,changed23),
               names_to = c(".value","tie"), names_pattern = "(ul|emb|dep_asymmetry|changed)_?(12|13|23)") %>%
  rename(unbal_load = ul, embeddedness = emb) %>%
  mutate(z_unbal_load = as.numeric(scale(unbal_load)), z_dep_asym = as.numeric(scale(dep_asymmetry)), z_emb = as.numeric(scale(embeddedness)))
surv_no <- Surv(rep(1, nrow(long_no_overlap)), long_no_overlap$changed)
m_no_overlap <- do.call(coxph, list(formula = surv_no ~ z_unbal_load + z_dep_asym + z_emb + strata(event_id), data = long_no_overlap, method = "efron"))
ct_no_overlap <- coeftest(m_no_overlap, vcov = vcovCL(m_no_overlap, cluster = long_no_overlap$triad_id))

sink("results/h2_eligible_and_overlap_exclusion_output.txt")
cat("=== PART A: H2 with eligible-alliance choice set ===\n")
cat("Alliance-termination events:", nrow(alliance_end), "\n")
cat("N alliances at t1 (composition):\n"); print(table(alliance_end$n_alliances_t1))
cat("\n--- Pooled choice set (all 3 candidate ties, original spec) ---\n"); print(ct_pooled)
cat("\n--- ELIGIBLE choice set (alliances only) ---\n"); print(ct_elig)

cat("\n\n=== PART B: excluding alliance+MID overlap events, H1 ===\n")
cat("Overlap dyad-years in the full panel:", nrow(overlap), "\n")
cat("Single-tie-change events where the CHANGED tie was an overlap dyad at t1:",
    sum(strict$changed_tie_is_overlap, na.rm = TRUE), "of", nrow(strict), "\n")
cat("Events remaining after exclusion:", nrow(strict_no_overlap), "\n")
print(ct_no_overlap)
sink()

cat(readLines("results/h2_eligible_and_overlap_exclusion_output.txt"), sep = "\n")
cat("\n\n18_h2_eligible_choiceset_and_overlap_exclusion.R complete.\n")
