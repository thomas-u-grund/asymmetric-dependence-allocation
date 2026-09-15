# ------------------------------------------------------------------
# 06: Robustness battery for the primary H1/H2 hazard result
#     (dep_asym_max): contiguity/regime-type controls and two-way
#     triad x bloc clustering.
#
#     (1) contiguity + Polity5 regime-type controls added on top of
#         the primary spec (base + z_trade_mean + z_dep_asym_max)
#     (2) two-way clustering by triad_id AND primary multilateral-bloc
#         membership, to check the triad-clustered SEs aren't
#         understating uncertainty given how many triads share bloc
#         membership (NATO, Warsaw Pact, etc.)
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(readxl); library(sandwich); library(lmtest); library(tidyr)
})


py <- readRDS("data/person_year_risk_table.rds")
py_trade <- py %>% filter(n_valid_ties == 3, year >= 1870) %>%
  mutate(
    z_trade_mean   = as.numeric(scale(trade_dependence_mean)),
    z_dep_asym_max = as.numeric(scale(dep_asym_max))
  )
cat("Trade-valid analysis sample:", nrow(py_trade), "person-years\n\n")

zscore <- function(x) as.numeric(scale(x))
base_rhs <- "z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin + z_trade_mean + z_dep_asym_max"

# ------------------------------------------------------------------
# (1) Contiguity + regime type
# ------------------------------------------------------------------
# node1/node2/node3 are already present in py_trade (carried through from
# the person-year panel), no need to re-join them.
py2 <- py_trade

contig <- read_csv("raw/contiguity/DirectContiguity320/contdird.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(pmin(state1no, state2no)), ccode_high = as.character(pmax(state1no, state2no))) %>%
  distinct(year, ccode_low, ccode_high) %>% mutate(contiguous = 1L)
get_contig <- function(df, a, b, newcol) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(contig, by = c("year" = "year", "lo" = "ccode_low", "hi" = "ccode_high")) %>%
    rename(!!newcol := contiguous) %>% select(-lo, -hi)
}
py2 <- py2 %>%
  get_contig("node1", "node2", "c12") %>% get_contig("node1", "node3", "c13") %>% get_contig("node2", "node3", "c23") %>%
  mutate(across(c(c12, c13, c23), ~ replace_na(.x, 0))) %>%
  rowwise() %>% mutate(triad_contig_share = mean(c(c12, c13, c23))) %>% ungroup()

polity <- read_excel("raw/polity/p5v2018.xls") %>%
  transmute(ccode = as.character(ccode), year, polity2) %>%
  filter(!is.na(polity2)) %>% distinct(ccode, year, .keep_all = TRUE)
get_polity <- function(df, node, newcol) {
  df %>% left_join(polity, by = c("year" = "year", setNames("ccode", node))) %>% rename(!!newcol := polity2)
}
py2 <- py2 %>%
  get_polity("node1", "p1") %>% get_polity("node2", "p2") %>% get_polity("node3", "p3") %>%
  rowwise() %>%
  mutate(triad_min_polity2 = min(c(p1, p2, p3), na.rm = ifelse(all(is.na(c(p1,p2,p3))), TRUE, FALSE))) %>%
  ungroup() %>% mutate(triad_min_polity2 = ifelse(is.infinite(triad_min_polity2), NA, triad_min_polity2))

py_ir <- py2 %>% filter(!is.na(triad_min_polity2)) %>%
  mutate(z_contig = zscore(triad_contig_share), z_polity = zscore(triad_min_polity2))
cat("IR-controls robustness sample:", nrow(py_ir), "person-years,", n_distinct(py_ir$triad_id), "triads\n")

m_d_ir <- glm(as.formula(paste("event_dissolve ~", base_rhs, "+ z_contig + z_polity")), data = py_ir, family = binomial())
m_f_ir <- glm(as.formula(paste("event_flip ~",     base_rhs, "+ z_contig + z_polity")), data = py_ir, family = binomial())
ct_d_ir <- coeftest(m_d_ir, vcov = vcovCL(m_d_ir, cluster = py_ir$triad_id))
ct_f_ir <- coeftest(m_f_ir, vcov = vcovCL(m_f_ir, cluster = py_ir$triad_id))

# ------------------------------------------------------------------
# (2) Two-way clustering: triad_id x primary multilateral-bloc membership
# ------------------------------------------------------------------
mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE)
big <- mem %>% group_by(version4id) %>% summarise(n_members = n_distinct(ccode), .groups = "drop") %>% filter(n_members >= 8)
ALL_IDS <- big$version4id
cat("\nMultilateral alliances with >=8 members:", length(ALL_IDS), "\n")

get_members_at <- function(id, yr) {
  sub <- mem %>% filter(version4id == id, mem_st_year <= yr, (mem_end_year >= yr | is.na(mem_end_year) | mem_end_year == 0))
  unique(sub$ccode)
}

tr <- readRDS("data/triads_all_years.rds") %>% select(year, triad_id, node1, node2, node3)
years_present <- sort(unique(py_trade$year))
tag <- lapply(years_present, function(y) {
  ty <- tr %>% filter(year == y, triad_id %in% py_trade$triad_id[py_trade$year == y])
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

py_bloc <- py_trade %>% left_join(tag, by = c("triad_id", "year")) %>%
  mutate(primary_bloc = ifelse(is.na(primary_bloc), "none", primary_bloc))
cat("Triad-years in a big bloc:", sum(py_bloc$primary_bloc != "none"),
    sprintf("(%.1f%%)\n", 100 * mean(py_bloc$primary_bloc != "none")))

m_dis_primary <- glm(as.formula(paste("event_dissolve ~", base_rhs)), data = py_bloc, family = binomial())
m_flip_primary <- glm(as.formula(paste("event_flip ~", base_rhs)), data = py_bloc, family = binomial())
cl2 <- py_bloc[, c("triad_id", "primary_bloc")]

ct_dis_1way <- coeftest(m_dis_primary, vcov = vcovCL(m_dis_primary, cluster = py_bloc$triad_id))
ct_dis_2way <- coeftest(m_dis_primary, vcov = vcovCL(m_dis_primary, cluster = cl2, multi0 = TRUE))
ct_flip_1way <- coeftest(m_flip_primary, vcov = vcovCL(m_flip_primary, cluster = py_bloc$triad_id))
ct_flip_2way <- coeftest(m_flip_primary, vcov = vcovCL(m_flip_primary, cluster = cl2, multi0 = TRUE))

# ------------------------------------------------------------------
sink("results/robustness_output.txt")
cat("=== (1) DISSOLUTION, + contiguity + regime type (n=", nrow(py_ir), ") ===\n", sep = ""); print(ct_d_ir)
cat("\n=== (1) REALIGNMENT, + contiguity + regime type ===\n"); print(ct_f_ir)

cat("\n\n=== (2) Multilateral blocs (>=8 members):", length(ALL_IDS), "===\n")
cat("Triad-years in a big bloc:", sum(py_bloc$primary_bloc != "none"),
    sprintf("(%.1f%%)\n", 100 * mean(py_bloc$primary_bloc != "none")))

cat("\n=== (2) DISSOLUTION: z_dep_asym_max under 1-way vs 2-way (triad+bloc) clustering ===\n")
cat("1-way: b=", coef(m_dis_primary)["z_dep_asym_max"],
    " SE=", ct_dis_1way["z_dep_asym_max","Std. Error"], " p=", ct_dis_1way["z_dep_asym_max","Pr(>|z|)"], "\n")
cat("2-way: b=", coef(m_dis_primary)["z_dep_asym_max"],
    " SE=", ct_dis_2way["z_dep_asym_max","Std. Error"], " p=", ct_dis_2way["z_dep_asym_max","Pr(>|z|)"], "\n")

cat("\n=== (2) REALIGNMENT: z_dep_asym_max under 1-way vs 2-way (triad+bloc) clustering ===\n")
cat("1-way: b=", coef(m_flip_primary)["z_dep_asym_max"],
    " SE=", ct_flip_1way["z_dep_asym_max","Std. Error"], " p=", ct_flip_1way["z_dep_asym_max","Pr(>|z|)"], "\n")
cat("2-way: b=", coef(m_flip_primary)["z_dep_asym_max"],
    " SE=", ct_flip_2way["z_dep_asym_max","Std. Error"], " p=", ct_flip_2way["z_dep_asym_max","Pr(>|z|)"], "\n")

cat("\n\n=== Full 2-way clustered coefficient tables (for reference) ===\n")
cat("--- Dissolution ---\n"); print(ct_dis_2way)
cat("\n--- Realignment ---\n"); print(ct_flip_2way)
sink()

cat(readLines("results/robustness_output.txt"), sep = "\n")

saveRDS(list(m_d_ir = m_d_ir, m_f_ir = m_f_ir, ct_d_ir = ct_d_ir, ct_f_ir = ct_f_ir,
             m_dis_primary = m_dis_primary, m_flip_primary = m_flip_primary,
             ct_dis_1way = ct_dis_1way, ct_dis_2way = ct_dis_2way,
             ct_flip_1way = ct_flip_1way, ct_flip_2way = ct_flip_2way),
        "results/robustness.rds")
cat("\n\n06_robustness.R complete. See results/robustness_output.txt\n")
