# ------------------------------------------------------------------
# 28: Full multiway clustering on the M4 tie-choice model (Table 1),
#     combining triad, primary bloc, dyad, and both constituent states
#     simultaneously via the Cameron-Gelbach-Miller sandwich estimator
#     (sandwich::vcovCL with a multi-column cluster data frame).
#
#     Appendix C.2 already clusters by dyad and by each of the two
#     constituent states, separately and pairwise with triad -- this
#     script is the single combined estimate that subsumes all of
#     those as special cases, on the exact M4 specification and sample
#     (n=3,456).
#
#     NOTE: the "state" dimension used here (splitting each tie into a
#     lower-ccode and higher-ccode column) is flawed -- see script 29's
#     header. Only this script's triad+bloc+dyad results are reported
#     in the paper (Appendix C.2); the state-clustering rows are
#     superseded by script 30's dyadic-robust estimator.
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest)
})


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
         mutual_min_23 = pmin(trade_share_low_23, trade_share_high_23))

# dyad and state identifiers for each of the triad's three candidate ties
strict <- strict %>%
  mutate(
    dyad_12 = paste0(pmin(node1, node2), "_", pmax(node1, node2)),
    dyad_13 = paste0(pmin(node1, node3), "_", pmax(node1, node3)),
    dyad_23 = paste0(pmin(node2, node3), "_", pmax(node2, node3)),
    state_lo_12 = pmin(node1, node2), state_hi_12 = pmax(node1, node2),
    state_lo_13 = pmin(node1, node3), state_hi_13 = pmax(node1, node3),
    state_lo_23 = pmin(node2, node3), state_hi_23 = pmax(node2, node3)
  )

# bloc tagging, same construction as script 07 (primary multilateral
# alliance bloc a triad belongs to at t1, "none" if it belongs to none
# of the 19 blocs with >=8 members)
mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE)
big <- mem %>% group_by(version4id) %>% summarise(n_members = n_distinct(ccode), .groups = "drop") %>% filter(n_members >= 8)
ALL_IDS <- big$version4id
get_members_at <- function(id, yr) {
  sub <- mem %>% filter(version4id == id, mem_st_year <= yr, (mem_end_year >= yr | is.na(mem_end_year) | mem_end_year == 0))
  unique(sub$ccode)
}
tr <- readRDS("data/triads_all_years.rds") %>% select(year, triad_id, node1, node2, node3)
tag_triads <- function(triad_year_df) {
  years_present <- sort(unique(triad_year_df$year))
  lapply(years_present, function(y) {
    want_ids <- triad_year_df$triad_id[triad_year_df$year == y]
    ty <- tr %>% filter(year == y, triad_id %in% want_ids)
    if (nrow(ty) == 0) return(NULL)
    primary_bloc <- rep(NA_integer_, nrow(ty))
    for (id in ALL_IDS) {
      m <- get_members_at(id, y)
      if (length(m) == 0) next
      hit <- rowSums(cbind(ty$node1 %in% m, ty$node2 %in% m, ty$node3 %in% m)) == 3
      primary_bloc[is.na(primary_bloc) & hit] <- id
    }
    tibble(triad_id = ty$triad_id, year = y,
           primary_bloc = ifelse(is.na(primary_bloc), "none", as.character(primary_bloc)))
  }) %>% bind_rows()
}
tagB <- tag_triads(strict %>% transmute(triad_id, year = t1) %>% distinct())
strict <- strict %>% left_join(tagB, by = c("triad_id", "t1" = "year")) %>%
  mutate(primary_bloc = ifelse(is.na(primary_bloc), "none", primary_bloc))

# long format -- identical construction to script 19's M4, with the
# cluster identifiers carried along for each candidate tie
long <- strict %>%
  select(event_id, triad_id, primary_bloc, sign12_t1, sign13_t1, sign23_t1,
         ul12, ul13, ul23, emb12, emb13, emb23, dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23,
         cap_asymmetry_12, cap_asymmetry_13, cap_asymmetry_23,
         mutual_min_12, mutual_min_13, mutual_min_23, changed12, changed13, changed23,
         dyad_12, dyad_13, dyad_23, state_lo_12, state_hi_12, state_lo_13, state_hi_13, state_lo_23, state_hi_23) %>%
  pivot_longer(cols = c(ul12,ul13,ul23,emb12,emb13,emb23,dep_asymmetry_12,dep_asymmetry_13,dep_asymmetry_23,
                        cap_asymmetry_12,cap_asymmetry_13,cap_asymmetry_23,
                        mutual_min_12,mutual_min_13,mutual_min_23,changed12,changed13,changed23,
                        dyad_12,dyad_13,dyad_23,state_lo_12,state_hi_12,state_lo_13,state_hi_13,state_lo_23,state_hi_23),
               names_to = c(".value","tie"),
               names_pattern = "(ul|emb|dep_asymmetry|cap_asymmetry|mutual_min|changed|dyad|state_lo|state_hi)_?(12|13|23)") %>%
  left_join(strict %>% select(event_id, sign12_t1, sign13_t1, sign23_t1) %>%
              pivot_longer(cols = c(sign12_t1,sign13_t1,sign23_t1), names_to = "tie2", values_to = "sign_t1") %>%
              mutate(tie2 = gsub("sign|_t1","",tie2)),
            by = c("event_id","tie" = "tie2")) %>%
  mutate(initial_negative = as.numeric(sign_t1 == -1),
         z_unbal_load = as.numeric(scale(ul)), z_dep_asym = as.numeric(scale(dep_asymmetry)),
         z_emb = as.numeric(scale(emb)), z_cap_asym = as.numeric(scale(cap_asymmetry)),
         z_mutual_min = as.numeric(scale(mutual_min)))

# M4 sample: complete cases on mutual_min and cap_asym, exactly as script 19
long4 <- long %>% filter(!is.na(z_mutual_min), !is.na(z_cap_asym))
cat("M4 sample: N events =", n_distinct(long4$event_id), " N tie-rows =", nrow(long4), "\n")
cat("Distinct triads:", n_distinct(long4$triad_id),
    " blocs:", n_distinct(long4$primary_bloc),
    " dyads:", n_distinct(long4$dyad),
    " states (lo):", n_distinct(long4$state_lo),
    " states (hi):", n_distinct(long4$state_hi), "\n\n")

m4 <- coxph(Surv(rep(1, nrow(long4)), long4$changed) ~
              z_unbal_load + z_dep_asym + initial_negative + z_mutual_min + z_cap_asym + z_emb + strata(event_id),
            data = long4, method = "efron")

ct_1way   <- coeftest(m4, vcov = vcovCL(m4, cluster = long4$triad_id))
ct_2way   <- coeftest(m4, vcov = vcovCL(m4, cluster = long4[, c("triad_id","primary_bloc")], multi0 = TRUE))
ct_3way   <- coeftest(m4, vcov = vcovCL(m4, cluster = long4[, c("triad_id","primary_bloc","dyad")], multi0 = TRUE))
ct_4way   <- coeftest(m4, vcov = vcovCL(m4, cluster = long4[, c("triad_id","primary_bloc","dyad","state_lo")], multi0 = TRUE))
ct_5way   <- tryCatch(
  coeftest(m4, vcov = vcovCL(m4, cluster = long4[, c("triad_id","primary_bloc","dyad","state_lo","state_hi")], multi0 = TRUE)),
  error = function(e) { cat("5-way clustering FAILED:", conditionMessage(e), "\n"); NULL }
)

sink("results/full_multiway_clustering_output.txt")
cat("=== Full multiway clustering on M4 (Table 1), progressively adding cluster dimensions ===\n")
cat("N events =", n_distinct(long4$event_id), " N tie-rows =", nrow(long4), "\n\n")
cat("--- 1-way: triad ---\n"); print(ct_1way)
cat("\n--- 2-way: triad + bloc ---\n"); print(ct_2way)
cat("\n--- 3-way: triad + bloc + dyad ---\n"); print(ct_3way)
cat("\n--- 4-way: triad + bloc + dyad + state(lower ccode) ---\n"); print(ct_4way)
if (!is.null(ct_5way)) { cat("\n--- 5-way: triad + bloc + dyad + state(lower) + state(higher) ---\n"); print(ct_5way) }
sink()
cat(readLines("results/full_multiway_clustering_output.txt"), sep = "\n")

saveRDS(list(long4 = long4, m4 = m4, ct_1way = ct_1way, ct_2way = ct_2way,
             ct_3way = ct_3way, ct_4way = ct_4way, ct_5way = ct_5way),
        "results/full_multiway_clustering.rds")
cat("\n\n28_full_multiway_clustering.R complete.\n")
