# ------------------------------------------------------------------
# 14: Response to external review, priority item 5/10 (partial):
#     does the main result depend on coding ANY of the four ATOP
#     alliance types (defense, neutrality, non-aggression, entente) as
#     a positive tie, or does it survive the stricter, more common
#     definition of "ally" as a DEFENSE PACT specifically? Self-
#     contained rebuild: signed ties -> triads -> spells -> single-tie
#     realignment events -> tie-choice model, condensed into one script
#     since only the final coefficient, not the full battery, is
#     needed for this check.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(igraph); library(survival); library(sandwich); library(lmtest) })

raw <- "raw"

# --- 1. signed ties, DEFENSE PACTS ONLY (vs. all 4 ATOP types in the main analysis) ---
alliances <- read_csv(file.path(raw, "alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv"), show_col_types = FALSE) %>%
  mutate(ccode_low = pmin(ccode1, ccode2), ccode_high = pmax(ccode1, ccode2)) %>%
  filter(defense == 1) %>%
  mutate(sign = 1) %>% select(year, ccode_low, ccode_high, sign) %>% distinct()
mids <- read_csv(file.path(raw, "mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv"), show_col_types = FALSE) %>%
  mutate(ccode_low = pmin(statea, stateb), ccode_high = pmax(statea, stateb)) %>%
  filter(hihost >= 3) %>% mutate(sign = -1) %>% select(year, ccode_low, ccode_high, sign) %>% distinct()
signed_ties <- bind_rows(alliances, mids) %>%
  group_by(year, ccode_low, ccode_high) %>%
  summarise(sign = ifelse(any(sign == -1), -1, 1), .groups = "drop") %>%
  filter(year >= 1816, year <= 2012)
cat("Defense-only signed ties:", nrow(signed_ties), " (positive:", sum(signed_ties$sign==1), " negative:", sum(signed_ties$sign==-1), ")\n")

# --- 2. triads, balance, embeddedness, unbalanced load (condensed) ---
years <- sort(unique(signed_ties$year))
process_year <- function(yr) {
  d <- signed_ties %>% filter(year == yr) %>%
    mutate(ccode1 = as.character(ccode_low), ccode2 = as.character(ccode_high)) %>%
    select(ccode1, ccode2, sign) %>% distinct(ccode1, ccode2, .keep_all = TRUE)
  if (nrow(d) < 3) return(NULL)
  g <- graph_from_data_frame(d[, c("ccode1","ccode2")], directed = FALSE)
  E(g)$sign <- d$sign
  tri <- igraph::triangles(g)
  if (length(tri) == 0) return(NULL)
  tri_mat <- matrix(V(g)$name[tri], ncol = 3, byrow = TRUE)
  el <- as_data_frame(g, what = "edges") %>% mutate(e1 = pmin(from, to), e2 = pmax(from, to))
  embeddedness <- sapply(seq_len(ecount(g)), function(i) {
    ends_i <- ends(g, i); length(intersect(neighbors(g, ends_i[1]), neighbors(g, ends_i[2])))
  })
  edge_tbl <- tibble(e1 = el$e1, e2 = el$e2, sign = E(g)$sign, embeddedness = embeddedness)
  ids_sorted <- t(apply(tri_mat, 1, sort))
  s1 <- ids_sorted[,1]; s2 <- ids_sorted[,2]; s3 <- ids_sorted[,3]
  idx12 <- match(paste(s1,s2), paste(edge_tbl$e1, edge_tbl$e2))
  idx13 <- match(paste(s1,s3), paste(edge_tbl$e1, edge_tbl$e2))
  idx23 <- match(paste(s2,s3), paste(edge_tbl$e1, edge_tbl$e2))
  sign12 <- edge_tbl$sign[idx12]; sign13 <- edge_tbl$sign[idx13]; sign23 <- edge_tbl$sign[idx23]
  balanced <- (sign12*sign13*sign23) > 0
  edge_in_triangle <- table(c(idx12,idx13,idx23))
  tibble(year = yr, node1 = s1, node2 = s2, node3 = s3, sign12, sign13, sign23, balanced = balanced,
         tie12_emb = edge_tbl$embeddedness[idx12], tie13_emb = edge_tbl$embeddedness[idx13], tie23_emb = edge_tbl$embeddedness[idx23],
         tie_emb_mean = (edge_tbl$embeddedness[idx12]+edge_tbl$embeddedness[idx13]+edge_tbl$embeddedness[idx23])/3,
         triad_id = paste(s1,s2,s3,sep="-"))
}
results <- lapply(years, function(yr) tryCatch(process_year(yr), error=function(e) NULL))
triads_all <- bind_rows(results)
cat("Defense-only triads:", nrow(triads_all), " unbalanced:", sum(!triads_all$balanced), "\n")

# --- 3. spells + single-tie realignment events (condensed) ---
tr <- triads_all %>% select(year, triad_id, node1, node2, node3, balanced, sign12, sign13, sign23,
                             tie12_emb, tie13_emb, tie23_emb) %>% arrange(triad_id, year)
tr <- tr %>% group_by(triad_id) %>%
  mutate(prev_year = lag(year), prev_balanced = lag(balanced),
         is_onset = !balanced & (is.na(prev_year) | prev_year != year-1 | prev_balanced),
         spell_id_local = cumsum(is_onset)) %>% ungroup() %>%
  mutate(spell_id = paste(triad_id, spell_id_local, sep="_S"))
LAST_YEAR <- max(tr$year)
spell_bounds <- tr %>% filter(!balanced) %>% group_by(spell_id, triad_id) %>% summarise(t1 = max(year), .groups="drop")
next_status <- tr %>% select(year, triad_id, balanced) %>% rename(next_year=year, balanced_next=balanced)
spell_bounds <- spell_bounds %>% mutate(lookup_year = t1+1) %>%
  left_join(next_status, by=c("triad_id","lookup_year"="next_year")) %>%
  mutate(event = case_when(t1>=LAST_YEAR ~ "censored_end_of_panel", is.na(balanced_next) ~ "dissolved",
                            balanced_next ~ "flipped", TRUE ~ "ERROR"))
flipped <- spell_bounds %>% filter(event=="flipped")
cat("Defense-only realignment events:", nrow(flipped), "\n")

at_t1 <- tr %>% rename(sign12_t1=sign12, sign13_t1=sign13, sign23_t1=sign23, emb12=tie12_emb, emb13=tie13_emb, emb23=tie23_emb)
at_t1p1 <- tr %>% select(year, triad_id, sign12, sign13, sign23) %>% rename(sign12_t2=sign12, sign13_t2=sign13, sign23_t2=sign23)
events <- flipped %>% left_join(at_t1, by=c("triad_id","t1"="year")) %>% left_join(at_t1p1, by=c("triad_id","lookup_year"="year")) %>%
  mutate(changed12 = sign12_t1!=sign12_t2, changed13 = sign13_t1!=sign13_t2, changed23 = sign23_t1!=sign23_t2,
         n_changed = changed12+changed13+changed23)
single <- events %>% filter(n_changed==1)
cat("Defense-only single-tie-change events:", nrow(single), "\n")

# --- 4. unbalanced load ---
edge_long <- bind_rows(
  tr %>% transmute(year, triad_id, balanced, e_lo=pmin(node1,node2), e_hi=pmax(node1,node2)),
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

# --- 5. trade-dependence asymmetry per tie ---
tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high, dep_asymmetry)
get_trade_var <- function(df, a, b_, suffix) {
  df %>% mutate(lo=pmin(.data[[a]],.data[[b_]]), hi=pmax(.data[[a]],.data[[b_]])) %>%
    left_join(tie_trade, by=c("t1"="year","lo"="ccode_low","hi"="ccode_high")) %>%
    rename(!!paste0("dep_asymmetry_",suffix) := dep_asymmetry) %>% select(-lo,-hi)
}
single <- single %>% get_trade_var("node1","node2","12") %>% get_trade_var("node1","node3","13") %>% get_trade_var("node2","node3","23")
strict <- single %>% filter(!is.na(dep_asymmetry_12), !is.na(dep_asymmetry_13), !is.na(dep_asymmetry_23))
cat("Defense-only events, all 3 ties trade-valid:", nrow(strict), "\n")

long <- strict %>% mutate(event_id = row_number()) %>%
  select(event_id, triad_id, ul12, ul13, ul23, emb12, emb13, emb23,
         dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, changed12, changed13, changed23) %>%
  pivot_longer(cols = c(ul12,ul13,ul23,emb12,emb13,emb23,dep_asymmetry_12,dep_asymmetry_13,dep_asymmetry_23,changed12,changed13,changed23),
               names_to = c(".value","tie"), names_pattern = "(ul|emb|dep_asymmetry|changed)_?(12|13|23)") %>%
  rename(unbal_load = ul, embeddedness = emb) %>%
  mutate(z_unbal_load = as.numeric(scale(unbal_load)), z_dep_asym = as.numeric(scale(dep_asymmetry)),
         z_emb = as.numeric(scale(embeddedness)))

surv <- Surv(rep(1, nrow(long)), long$changed)
m <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + z_emb + strata(event_id), data = long, method = "efron"))
ct <- coeftest(m, vcov = vcovCL(m, cluster = long$triad_id))

sink("results/defense_only_network_output.txt")
cat("=== DEFENSE-PACTS-ONLY alliance definition (vs. all 4 ATOP types in main analysis) ===\n")
cat("Signed ties:", nrow(signed_ties), " | triads:", nrow(triads_all), " | unbalanced:", sum(!triads_all$balanced), "\n")
cat("Realignment events:", nrow(flipped), " | single-tie:", nrow(single), " | trade-valid:", nrow(strict), "\n\n")
print(ct)
sink()
cat(readLines("results/defense_only_network_output.txt"), sep = "\n")
cat("\n\n14_defense_only_network.R complete.\n")
