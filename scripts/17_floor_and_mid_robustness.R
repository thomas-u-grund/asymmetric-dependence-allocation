# ------------------------------------------------------------------
# 17: Two measurement-sensitivity checks.
#
#     (A) The $50m minimum-national-trade floor used to compute trade
#     shares (script 01) is a fixed NOMINAL threshold applied across
#     1870-2012 -- if it is not inflation-adjusted, its effective
#     stringency changes enormously over 140 years. Test sensitivity of
#     the main tie-choice result to this floor: no floor, $10m, $50m
#     (main analysis), $200m.
#
#     (B) Quantify how often a dyad-year is coded negative because an
#     alliance and a MID coexist (rather than because two states are
#     hostile with no alliance at all) -- a neg_to_pos transition in
#     the former case can just mean the MID ended while a pre-existing
#     alliance continued, a different event from two hostile states
#     allying. Also test sensitivity to the MID hostility threshold
#     (hihost >= 3, main analysis, vs >= 4 and >= 5).
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(igraph); library(survival); library(sandwich); library(lmtest) })

raw <- "raw"

# ==================================================================
# PART A: trade-floor sensitivity (reuses the EXISTING triads/spells/
# events -- only the trade merge changes, so no network rebuild needed)
# ==================================================================
dyadic_trade <- read_csv(file.path(raw, "trade/COW_Trade_4.0/Dyadic_COW_4.0.csv"), show_col_types = FALSE) %>%
  mutate(across(c(flow1, flow2, smoothtotrade), ~ na_if(., -9))) %>%
  mutate(ccode_low = pmin(ccode1, ccode2), ccode_high = pmax(ccode1, ccode2)) %>%
  group_by(year, ccode_low, ccode_high) %>%
  summarise(total_dyadic_trade = sum(smoothtotrade, na.rm = TRUE), .groups = "drop")
national_trade_raw <- read_csv(file.path(raw, "trade/COW_Trade_4.0/National_COW_4.0.csv"), show_col_types = FALSE) %>%
  mutate(total_trade_raw = imports + exports) %>%
  select(ccode, year, total_trade_raw)

build_tie_trade <- function(floor) {
  national_trade <- national_trade_raw %>%
    mutate(total_trade = ifelse(total_trade_raw < floor, NA, total_trade_raw)) %>%
    select(ccode, year, total_trade)
  dyadic_trade %>%
    left_join(national_trade, by = c("ccode_low" = "ccode", "year")) %>% rename(total_trade_low = total_trade) %>%
    left_join(national_trade, by = c("ccode_high" = "ccode", "year")) %>% rename(total_trade_high = total_trade) %>%
    mutate(
      overflow = total_dyadic_trade > pmin(total_trade_low, total_trade_high, na.rm = TRUE),
      trade_share_low  = ifelse(overflow, NA, total_dyadic_trade / total_trade_low),
      trade_share_high = ifelse(overflow, NA, total_dyadic_trade / total_trade_high),
      dep_asymmetry = ifelse(!is.na(trade_share_low) & !is.na(trade_share_high),
                              abs(trade_share_low - trade_share_high), NA)
    ) %>%
    select(year, ccode_low, ccode_high, dep_asymmetry)
}

h3 <- readRDS("results/h3_displacement.rds")
# single already has ul12/13/23, emb12/13/23, and a dep_asymmetry_* built with the
# main analysis's $50m floor -- drop those so each floor variant below can be
# re-merged under its own column names without collision
strict_base <- h3$single %>%
  select(-dep_asymmetry_12, -dep_asymmetry_13, -dep_asymmetry_23,
         -cap_asymmetry_12, -cap_asymmetry_13, -cap_asymmetry_23,
         -exposed_is_weaker_12, -exposed_is_weaker_13, -exposed_is_weaker_23,
         -log_trade_volume_12, -log_trade_volume_13, -log_trade_volume_23,
         -n_trade_valid)

fit_floor <- function(floor, label) {
  tt <- build_tie_trade(floor) %>% mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high))
  get_dep <- function(df, a, b_, suffix) {
    df %>% mutate(lo = pmin(.data[[a]], .data[[b_]]), hi = pmax(.data[[a]], .data[[b_]])) %>%
      left_join(tt, by = c("t1" = "year", "lo" = "ccode_low", "hi" = "ccode_high")) %>%
      rename(!!paste0("dep_asymmetry_", suffix) := dep_asymmetry) %>% select(-lo, -hi)
  }
  d <- strict_base %>% get_dep("node1", "node2", "12") %>% get_dep("node1", "node3", "13") %>% get_dep("node2", "node3", "23") %>%
    filter(!is.na(dep_asymmetry_12), !is.na(dep_asymmetry_13), !is.na(dep_asymmetry_23))
  long <- d %>% mutate(event_id = row_number()) %>%
    select(event_id, triad_id, ul12, ul13, ul23, emb12, emb13, emb23,
           dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, changed12, changed13, changed23) %>%
    pivot_longer(cols = c(ul12, ul13, ul23, emb12, emb13, emb23, dep_asymmetry_12, dep_asymmetry_13,
                          dep_asymmetry_23, changed12, changed13, changed23),
                 names_to = c(".value", "tie"), names_pattern = "(ul|emb|dep_asymmetry|changed)_?(12|13|23)") %>%
    rename(unbal_load = ul, embeddedness = emb) %>%
    mutate(z_unbal_load = as.numeric(scale(unbal_load)), z_dep_asym = as.numeric(scale(dep_asymmetry)), z_emb = as.numeric(scale(embeddedness)))
  surv <- Surv(rep(1, nrow(long)), long$changed)
  m <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + z_emb + strata(event_id), data = long, method = "efron"))
  ct <- coeftest(m, vcov = vcovCL(m, cluster = long$triad_id))
  cat(sprintf("%-20s N events=%5d  b=%7.4f  SE=%6.4f  p=%8.5f\n", label, n_distinct(long$event_id),
              ct["z_dep_asym","Estimate"], ct["z_dep_asym","Std. Error"], ct["z_dep_asym","Pr(>|z|)"]))
  invisible(ct)
}

sink("results/floor_and_mid_robustness_output.txt")
cat("=== PART A: trade-floor sensitivity (national-trade minimum, in COW-reported units) ===\n")
fit_floor(0, "No floor")
fit_floor(10, "$10m floor")
fit_floor(50, "$50m floor (main)")
fit_floor(200, "$200m floor")

# ==================================================================
# PART B: alliance+MID overlap quantification, and MID threshold sensitivity
# ==================================================================
cat("\n\n=== PART B: how often does 'negative' mean alliance+MID coexisting, not hostility alone? ===\n")
alliances_raw <- read_csv(file.path(raw, "alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv"), show_col_types = FALSE) %>%
  mutate(ccode_low = pmin(ccode1, ccode2), ccode_high = pmax(ccode1, ccode2)) %>%
  filter(defense == 1 | neutrality == 1 | nonaggression == 1 | entente == 1) %>%
  select(year, ccode_low, ccode_high) %>% distinct() %>% mutate(has_alliance = TRUE)
mids_raw <- read_csv(file.path(raw, "mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv"), show_col_types = FALSE) %>%
  mutate(ccode_low = pmin(statea, stateb), ccode_high = pmax(statea, stateb))

for (thresh in c(3, 4, 5)) {
  m <- mids_raw %>% filter(hihost >= thresh) %>% select(year, ccode_low, ccode_high) %>% distinct() %>% mutate(has_mid = TRUE)
  overlap <- m %>% inner_join(alliances_raw, by = c("year", "ccode_low", "ccode_high"))
  cat(sprintf("hihost >= %d: %d negative dyad-years total, %d (%.1f%%) coexist with an alliance\n",
              thresh, nrow(m), nrow(overlap), 100 * nrow(overlap) / nrow(m)))
}

cat("\n=== MID threshold sensitivity: rebuilding the signed network and tie-choice test for hihost >= 4 ===\n")
build_and_test_mid <- function(mid_thresh) {
  alliances <- alliances_raw %>% mutate(sign = 1) %>% select(year, ccode_low, ccode_high, sign)
  mids <- mids_raw %>% filter(hihost >= mid_thresh) %>% mutate(sign = -1) %>% select(year, ccode_low, ccode_high, sign) %>% distinct()
  signed_ties <- bind_rows(alliances, mids) %>% group_by(year, ccode_low, ccode_high) %>%
    summarise(sign = ifelse(any(sign == -1), -1, 1), .groups = "drop") %>% filter(year >= 1816, year <= 2012)

  years <- sort(unique(signed_ties$year))
  process_year <- function(yr) {
    d <- signed_ties %>% filter(year == yr) %>% mutate(ccode1 = as.character(ccode_low), ccode2 = as.character(ccode_high)) %>%
      select(ccode1, ccode2, sign) %>% distinct(ccode1, ccode2, .keep_all = TRUE)
    if (nrow(d) < 3) return(NULL)
    g <- graph_from_data_frame(d[, c("ccode1","ccode2")], directed = FALSE); E(g)$sign <- d$sign
    tri <- igraph::triangles(g); if (length(tri) == 0) return(NULL)
    tri_mat <- matrix(V(g)$name[tri], ncol = 3, byrow = TRUE)
    el <- as_data_frame(g, what = "edges") %>% mutate(e1 = pmin(from, to), e2 = pmax(from, to))
    embeddedness <- sapply(seq_len(ecount(g)), function(i) { ends_i <- ends(g, i); length(intersect(neighbors(g, ends_i[1]), neighbors(g, ends_i[2]))) })
    edge_tbl <- tibble(e1 = el$e1, e2 = el$e2, sign = E(g)$sign, embeddedness = embeddedness)
    ids_sorted <- t(apply(tri_mat, 1, sort)); s1 <- ids_sorted[,1]; s2 <- ids_sorted[,2]; s3 <- ids_sorted[,3]
    idx12 <- match(paste(s1,s2), paste(edge_tbl$e1, edge_tbl$e2)); idx13 <- match(paste(s1,s3), paste(edge_tbl$e1, edge_tbl$e2)); idx23 <- match(paste(s2,s3), paste(edge_tbl$e1, edge_tbl$e2))
    sign12 <- edge_tbl$sign[idx12]; sign13 <- edge_tbl$sign[idx13]; sign23 <- edge_tbl$sign[idx23]
    balanced <- (sign12*sign13*sign23) > 0
    tibble(year = yr, node1 = s1, node2 = s2, node3 = s3, sign12, sign13, sign23, balanced = balanced,
           tie12_emb = edge_tbl$embeddedness[idx12], tie13_emb = edge_tbl$embeddedness[idx13], tie23_emb = edge_tbl$embeddedness[idx23],
           triad_id = paste(s1,s2,s3,sep="-"))
  }
  triads_all <- bind_rows(lapply(years, function(yr) tryCatch(process_year(yr), error=function(e) NULL)))
  tr <- triads_all %>% select(year, triad_id, node1, node2, node3, balanced, sign12, sign13, sign23, tie12_emb, tie13_emb, tie23_emb) %>% arrange(triad_id, year)
  tr <- tr %>% group_by(triad_id) %>% mutate(prev_year = lag(year), prev_balanced = lag(balanced),
    is_onset = !balanced & (is.na(prev_year) | prev_year != year-1 | prev_balanced), spell_id_local = cumsum(is_onset)) %>%
    ungroup() %>% mutate(spell_id = paste(triad_id, spell_id_local, sep="_S"))
  LAST_YEAR <- max(tr$year)
  spell_bounds <- tr %>% filter(!balanced) %>% group_by(spell_id, triad_id) %>% summarise(t1 = max(year), .groups="drop")
  next_status <- tr %>% select(year, triad_id, balanced) %>% rename(next_year=year, balanced_next=balanced)
  spell_bounds <- spell_bounds %>% mutate(lookup_year = t1+1) %>% left_join(next_status, by=c("triad_id","lookup_year"="next_year")) %>%
    mutate(event = case_when(t1>=LAST_YEAR ~ "censored_end_of_panel", is.na(balanced_next) ~ "dissolved", balanced_next ~ "flipped", TRUE ~ "ERROR"))
  flipped <- spell_bounds %>% filter(event=="flipped")
  at_t1 <- tr %>% rename(sign12_t1=sign12, sign13_t1=sign13, sign23_t1=sign23, emb12=tie12_emb, emb13=tie13_emb, emb23=tie23_emb)
  at_t1p1 <- tr %>% select(year, triad_id, sign12, sign13, sign23) %>% rename(sign12_t2=sign12, sign13_t2=sign13, sign23_t2=sign23)
  events <- flipped %>% left_join(at_t1, by=c("triad_id","t1"="year")) %>% left_join(at_t1p1, by=c("triad_id","lookup_year"="year")) %>%
    mutate(changed12 = sign12_t1!=sign12_t2, changed13 = sign13_t1!=sign13_t2, changed23 = sign23_t1!=sign23_t2, n_changed = changed12+changed13+changed23)
  single <- events %>% filter(n_changed==1)
  edge_long <- bind_rows(tr %>% transmute(year, triad_id, balanced, e_lo=pmin(node1,node2), e_hi=pmax(node1,node2)),
                          tr %>% transmute(year, triad_id, balanced, e_lo=pmin(node1,node3), e_hi=pmax(node1,node3)),
                          tr %>% transmute(year, triad_id, balanced, e_lo=pmin(node2,node3), e_hi=pmax(node2,node3)))
  edge_summary <- edge_long %>% group_by(year, e_lo, e_hi) %>% summarise(n_unbalanced=sum(!balanced), .groups="drop")
  get_unbal_load <- function(df, a, b_, newcol) {
    df %>% mutate(lo=pmin(.data[[a]],.data[[b_]]), hi=pmax(.data[[a]],.data[[b_]])) %>%
      left_join(edge_summary, by=c("t1"="year","lo"="e_lo","hi"="e_hi")) %>%
      rename(!!newcol := n_unbalanced) %>% mutate(!!newcol := !!sym(newcol) - 1) %>% select(-lo,-hi)
  }
  single <- single %>% get_unbal_load("node1","node2","ul12") %>% get_unbal_load("node1","node3","ul13") %>% get_unbal_load("node2","node3","ul23") %>%
    mutate(across(c(ul12,ul13,ul23), ~replace_na(.x, 0)))
  tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>% mutate(ccode_low=as.character(ccode_low), ccode_high=as.character(ccode_high)) %>% select(year,ccode_low,ccode_high,dep_asymmetry)
  get_trade_var <- function(df, a, b_, suffix) {
    df %>% mutate(lo=pmin(.data[[a]],.data[[b_]]), hi=pmax(.data[[a]],.data[[b_]])) %>%
      left_join(tie_trade, by=c("t1"="year","lo"="ccode_low","hi"="ccode_high")) %>%
      rename(!!paste0("dep_asymmetry_",suffix) := dep_asymmetry) %>% select(-lo,-hi)
  }
  single <- single %>% get_trade_var("node1","node2","12") %>% get_trade_var("node1","node3","13") %>% get_trade_var("node2","node3","23")
  strict <- single %>% filter(!is.na(dep_asymmetry_12), !is.na(dep_asymmetry_13), !is.na(dep_asymmetry_23))
  long <- strict %>% mutate(event_id=row_number()) %>%
    select(event_id, triad_id, ul12,ul13,ul23,emb12,emb13,emb23,dep_asymmetry_12,dep_asymmetry_13,dep_asymmetry_23,changed12,changed13,changed23) %>%
    pivot_longer(cols=c(ul12,ul13,ul23,emb12,emb13,emb23,dep_asymmetry_12,dep_asymmetry_13,dep_asymmetry_23,changed12,changed13,changed23),
                 names_to=c(".value","tie"), names_pattern="(ul|emb|dep_asymmetry|changed)_?(12|13|23)") %>%
    rename(unbal_load=ul, embeddedness=emb) %>%
    mutate(z_unbal_load=as.numeric(scale(unbal_load)), z_dep_asym=as.numeric(scale(dep_asymmetry)), z_emb=as.numeric(scale(embeddedness)))
  surv <- Surv(rep(1,nrow(long)), long$changed)
  m2 <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + z_emb + strata(event_id), data=long, method="efron"))
  ct <- coeftest(m2, vcov=vcovCL(m2, cluster=long$triad_id))
  cat(sprintf("hihost >= %d:  N events=%5d  b=%7.4f  SE=%6.4f  p=%8.5f\n", mid_thresh, n_distinct(long$event_id),
              ct["z_dep_asym","Estimate"], ct["z_dep_asym","Std. Error"], ct["z_dep_asym","Pr(>|z|)"]))
}
build_and_test_mid(4)
build_and_test_mid(5)
sink()

cat(readLines("results/floor_and_mid_robustness_output.txt"), sep = "\n")
cat("\n\n17_floor_and_mid_robustness.R complete.\n")
