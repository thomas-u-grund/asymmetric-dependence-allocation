# ------------------------------------------------------------------
# 15: Response to second round of external review, priority items 1-2.
#     (1) Split single-tie realignment events by transition DIRECTION:
#     alliance-to-hostility (pos_to_neg) vs. hostility-to-alliance
#     (neg_to_pos). Does trade-dependence asymmetry predict tie choice
#     in BOTH directions, or only one? Pooling them, as the main draft
#     does, implicitly assumes the mechanism is direction-agnostic.
#
#     (2) Build an EDGE-LEVEL model of dissolution: among triads that
#     resolve via a tie disappearing entirely, which of the 3 candidate
#     ties is the one that disappears? This is the direct analogue of
#     the realignment tie-choice model, needed to actually test the
#     manuscript's claim that capability asymmetry governs "exit" --
#     that claim currently rests only on a TRIAD-level hazard (does the
#     triad dissolve at all), not on WHICH tie dissolves. Also split by
#     whether the disappearing tie was positive or negative beforehand.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest) })


# ==================================================================
# PART 1: transition-direction split of the realignment tie-choice model
# ==================================================================
h3 <- readRDS("results/h3_displacement.rds")
s <- h3$strict %>%
  mutate(
    trans12 = case_when(changed12 & sign12_t1 == 1 & sign12_t2 == -1 ~ "pos_to_neg",
                         changed12 & sign12_t1 == -1 & sign12_t2 == 1 ~ "neg_to_pos"),
    trans13 = case_when(changed13 & sign13_t1 == 1 & sign13_t2 == -1 ~ "pos_to_neg",
                         changed13 & sign13_t1 == -1 & sign13_t2 == 1 ~ "neg_to_pos"),
    trans23 = case_when(changed23 & sign23_t1 == 1 & sign23_t2 == -1 ~ "pos_to_neg",
                         changed23 & sign23_t1 == -1 & sign23_t2 == 1 ~ "neg_to_pos"),
    transition = coalesce(trans12, trans13, trans23)
  )
cat("Transition-type counts:\n"); print(table(s$transition))

build_long <- function(df) {
  df %>% mutate(event_id = row_number()) %>%
    select(event_id, triad_id, ul12, ul13, ul23, emb12, emb13, emb23,
           dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, changed12, changed13, changed23) %>%
    pivot_longer(cols = c(ul12, ul13, ul23, emb12, emb13, emb23, dep_asymmetry_12, dep_asymmetry_13,
                           dep_asymmetry_23, changed12, changed13, changed23),
                 names_to = c(".value", "tie"), names_pattern = "(ul|emb|dep_asymmetry|changed)_?(12|13|23)") %>%
    rename(unbal_load = ul, embeddedness = emb) %>%
    mutate(z_unbal_load = as.numeric(scale(unbal_load)), z_dep_asym = as.numeric(scale(dep_asymmetry)),
           z_emb = as.numeric(scale(embeddedness)))
}
fit_report <- function(data, label) {
  n_ev <- n_distinct(data$event_id)
  surv <- Surv(rep(1, nrow(data)), data$changed)
  m <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + z_emb + strata(event_id),
                            data = data, method = "efron"))
  ct <- coeftest(m, vcov = vcovCL(m, cluster = data$triad_id))
  cat(sprintf("%-30s N events=%5d  b=%8.4f  SE=%7.4f  p=%8.5f\n",
              label, n_ev, ct["z_dep_asym","Estimate"], ct["z_dep_asym","Std. Error"], ct["z_dep_asym","Pr(>|z|)"]))
  list(m = m, ct = ct, n_ev = n_ev)
}

long_pos_to_neg <- build_long(s %>% filter(transition == "pos_to_neg"))
long_neg_to_pos <- build_long(s %>% filter(transition == "neg_to_pos"))

cat("\n=== Tie-choice model, BY TRANSITION DIRECTION ===\n")
r_p2n <- fit_report(long_pos_to_neg, "alliance -> hostility")
r_n2p <- fit_report(long_neg_to_pos, "hostility -> alliance")

# ==================================================================
# PART 2: edge-level dissolution model -- which tie disappears?
# ==================================================================
spell_bounds <- readRDS("data/spell_bounds.rds")
tr <- readRDS("data/triads_all_years.rds") %>%
  select(year, triad_id, node1, node2, node3, sign12, sign13, sign23,
         tie12_emb, tie13_emb, tie23_emb)

dissolved <- spell_bounds %>% filter(event == "dissolved")
cat("\nDissolution events (triad no longer a closed triad at t1+1):", nrow(dissolved), "\n")

at_t1 <- tr %>% rename(sign12_t1 = sign12, sign13_t1 = sign13, sign23_t1 = sign23,
                        emb12 = tie12_emb, emb13 = tie13_emb, emb23 = tie23_emb)
diss_events <- dissolved %>% left_join(at_t1, by = c("triad_id", "t1" = "year"))

# which specific edge is gone at t1+1? check each of the 3 dyads' presence
# in the FULL dyad-year signed-ties table (not just among closed triads,
# since the triad itself is no longer closed -- the missing edge could be
# any of the 3, or (rarely) more than one)
signed_ties <- read_csv("data/signed_ties.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high))
tie_present_at <- function(df, a, b_, year_col, newcol) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b_]]), hi = pmax(.data[[a]], .data[[b_]])) %>%
    left_join(signed_ties %>% select(year, ccode_low, ccode_high) %>% mutate(present = TRUE),
              by = setNames(c("year","ccode_low","ccode_high"), c(year_col,"lo","hi"))) %>%
    rename(!!newcol := present) %>% mutate(!!newcol := !is.na(!!sym(newcol))) %>% select(-lo,-hi)
}
diss_events <- diss_events %>%
  mutate(lookup_year = t1 + 1) %>%
  tie_present_at("node1","node2","lookup_year","present12") %>%
  tie_present_at("node1","node3","lookup_year","present13") %>%
  tie_present_at("node2","node3","lookup_year","present23") %>%
  mutate(gone12 = !present12, gone13 = !present13, gone23 = !present23,
         n_gone = gone12 + gone13 + gone23)

cat("Distribution of number of ties gone at dissolution:\n"); print(table(diss_events$n_gone))
diss_single <- diss_events %>% filter(n_gone == 1)
cat("Single-tie dissolutions:", nrow(diss_single), "of", nrow(diss_events), "\n")

diss_single <- diss_single %>%
  mutate(gone_sign = case_when(gone12 ~ sign12_t1, gone13 ~ sign13_t1, gone23 ~ sign23_t1))
cat("Sign of the disappearing tie just before it disappeared:\n"); print(table(diss_single$gone_sign))

# merge trade asymmetry AND capability asymmetry per candidate tie, at t1
tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high, dep_asymmetry, cap_asymmetry)
get_tv <- function(df, a, b_, suffix) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b_]]), hi = pmax(.data[[a]], .data[[b_]])) %>%
    left_join(tie_trade, by = c("t1"="year","lo"="ccode_low","hi"="ccode_high")) %>%
    rename_with(~ paste0(., "_", suffix), c(dep_asymmetry, cap_asymmetry)) %>% select(-lo,-hi)
}
diss_single <- diss_single %>%
  get_tv("node1","node2","12") %>% get_tv("node1","node3","13") %>% get_tv("node2","node3","23")

diss_valid <- diss_single %>%
  filter(!is.na(cap_asymmetry_12), !is.na(cap_asymmetry_13), !is.na(cap_asymmetry_23))
cat("Single-tie dissolutions with capability data on all 3 ties:", nrow(diss_valid), "\n")
diss_valid_trade <- diss_single %>%
  filter(!is.na(dep_asymmetry_12), !is.na(dep_asymmetry_13), !is.na(dep_asymmetry_23))
cat("...and with trade data on all 3 ties:", nrow(diss_valid_trade), "\n")

build_long_diss <- function(df, extra_cols) {
  df %>% mutate(event_id = row_number()) %>%
    select(event_id, triad_id, all_of(paste0(rep(c("emb","gone"), each=3), c(12,13,23))),
           all_of(paste0(rep(extra_cols, each=3), "_", c(12,13,23)))) %>%
    pivot_longer(cols = -c(event_id, triad_id),
                 names_to = c(".value","tie"),
                 names_pattern = paste0("(emb|gone|", paste(extra_cols, collapse="|"), ")_?(12|13|23)")) %>%
    rename(embeddedness = emb, gone_tie = gone)
}

sink("results/transition_split_and_edge_dissolution_output.txt")
cat("=== PART 1: tie-choice model by transition direction ===\n")
cat("alliance -> hostility (pos_to_neg): N events =", r_p2n$n_ev, "\n"); print(r_p2n$ct)
cat("\nhostility -> alliance (neg_to_pos): N events =", r_n2p$n_ev, "\n"); print(r_n2p$ct)

cat("\n\n=== PART 2: edge-level dissolution model ===\n")
cat("Dissolution events (triad-level):", nrow(dissolved), "\n")
cat("Single-tie dissolutions:", nrow(diss_single), "of", nrow(diss_events),
    sprintf(" (%.1f%%)\n", 100*nrow(diss_single)/nrow(diss_events)))
cat("Sign of disappearing tie beforehand:\n"); print(table(diss_single$gone_sign))
cat("With capability data on all 3 ties:", nrow(diss_valid), "\n")
cat("With trade data on all 3 ties:", nrow(diss_valid_trade), "\n")

if (nrow(diss_valid) > 30) {
  long_diss_cap <- build_long_diss(diss_valid, c("cap_asymmetry"))
  m_diss_cap <- do.call(coxph, list(formula = Surv(rep(1,nrow(long_diss_cap)), long_diss_cap$gone_tie) ~
                                       as.numeric(scale(cap_asymmetry)) + as.numeric(scale(embeddedness)) + strata(event_id),
                                     data = long_diss_cap, method = "efron"))
  ct_diss_cap <- coeftest(m_diss_cap, vcov = vcovCL(m_diss_cap, cluster = long_diss_cap$triad_id))
  cat("\n--- Does capability asymmetry predict WHICH tie disappears? (all single-tie dissolutions) ---\n")
  print(ct_diss_cap)

  # split by sign of the disappearing tie
  for (sgn in c(1, -1)) {
    sub <- diss_valid %>% filter(gone_sign == sgn)
    lbl <- ifelse(sgn == 1, "alliance disappearing (+ -> 0)", "MID ending (- -> 0)")
    cat("\n---", lbl, ": N =", nrow(sub), "---\n")
    if (nrow(sub) > 20) {
      ld <- build_long_diss(sub, c("cap_asymmetry"))
      mm <- tryCatch(do.call(coxph, list(formula = Surv(rep(1,nrow(ld)), ld$gone_tie) ~
                                            as.numeric(scale(cap_asymmetry)) + as.numeric(scale(embeddedness)) + strata(event_id),
                                          data = ld, method = "efron")), error = function(e) NULL)
      if (!is.null(mm)) print(coeftest(mm, vcov = vcovCL(mm, cluster = ld$triad_id)))
    } else cat("(too few events for a stable estimate)\n")
  }
} else {
  cat("\nToo few single-tie dissolutions with capability data for a stable edge-level model.\n")
}

if (nrow(diss_valid_trade) > 30) {
  long_diss_trade <- build_long_diss(diss_valid_trade, c("dep_asymmetry"))
  m_diss_trade <- do.call(coxph, list(formula = Surv(rep(1,nrow(long_diss_trade)), long_diss_trade$gone_tie) ~
                                         as.numeric(scale(dep_asymmetry)) + as.numeric(scale(embeddedness)) + strata(event_id),
                                       data = long_diss_trade, method = "efron"))
  ct_diss_trade <- coeftest(m_diss_trade, vcov = vcovCL(m_diss_trade, cluster = long_diss_trade$triad_id))
  cat("\n--- Does TRADE asymmetry predict WHICH tie disappears? (pooled, alliance+MID together) ---\n")
  print(ct_diss_trade)

  # split by prior sign: alliance ending (a genuine relationship loss) vs. MID
  # ending (arguably de-escalation, not loss) -- these are not substantively
  # equivalent and the "protected from disappearing" interpretation may only
  # cleanly apply to one of them
  for (sgn in c(1, -1)) {
    sub <- diss_valid_trade %>% filter(gone_sign == sgn)
    lbl <- ifelse(sgn == 1, "alliance disappearing (+ -> 0)", "MID ending (- -> 0)")
    cat("\n---", lbl, ": N =", nrow(sub), "---\n")
    if (nrow(sub) > 20) {
      ld <- build_long_diss(sub, c("dep_asymmetry"))
      mm <- tryCatch(do.call(coxph, list(formula = Surv(rep(1,nrow(ld)), ld$gone_tie) ~
                                            as.numeric(scale(dep_asymmetry)) + as.numeric(scale(embeddedness)) + strata(event_id),
                                          data = ld, method = "efron")), error = function(e) NULL)
      if (!is.null(mm)) print(coeftest(mm, vcov = vcovCL(mm, cluster = ld$triad_id)))
    } else cat("(too few events for a stable estimate)\n")
  }
} else {
  cat("\nToo few single-tie dissolutions with trade data for a stable edge-level model.\n")
}
sink()

cat(readLines("results/transition_split_and_edge_dissolution_output.txt"), sep = "\n")
cat("\n\n15_transition_split_and_edge_dissolution.R complete.\n")
