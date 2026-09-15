# ------------------------------------------------------------------
# 07: Follow-up to script 06's two-way clustering result, which showed
#     the H2 realignment-suppression effect (dep_asym_max) survives
#     contiguity/regime controls but loses significance under two-way
#     (triad + bloc) clustered SEs -- unsurprising given 85% of the
#     trade-valid sample sits inside one of 19 big multilateral blocs.
#
#     Two-way clustering asks "is the SE understated given bloc-level
#     non-independence" -- it does not ask whether the effect is real
#     WITHIN blocs vs. only an artifact of BETWEEN-bloc differences
#     (e.g. some blocs are just more trade-asymmetric AND less prone
#     to realignment for unrelated reasons). This script asks that
#     more informative question directly, via bloc fixed effects, and
#     separately extends the two-way clustering check to H3 (the tie-
#     choice conditional logit), which hasn't been tested against it.
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest); library(survival)
})


mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE)
big <- mem %>% group_by(version4id) %>% summarise(n_members = n_distinct(ccode), .groups = "drop") %>% filter(n_members >= 8)
ALL_IDS <- big$version4id
get_members_at <- function(id, yr) {
  sub <- mem %>% filter(version4id == id, mem_st_year <= yr, (mem_end_year >= yr | is.na(mem_end_year) | mem_end_year == 0))
  unique(sub$ccode)
}
tr <- readRDS("data/triads_all_years.rds") %>% select(year, triad_id, node1, node2, node3)

tag_triads <- function(triad_year_df) {
  # triad_year_df: data frame with columns triad_id, year (the years to tag)
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

# ==================================================================
# PART A: H2, bloc fixed effects (does dep_asym_max hold WITHIN blocs?)
# ==================================================================
py <- readRDS("data/person_year_risk_table.rds")
py_trade <- py %>% filter(n_valid_ties == 3, year >= 1870) %>%
  mutate(
    z_trade_mean   = as.numeric(scale(trade_dependence_mean)),
    z_dep_asym_max = as.numeric(scale(dep_asym_max))
  )

tagA <- tag_triads(py_trade %>% distinct(triad_id, year))
py_bloc <- py_trade %>% left_join(tagA, by = c("triad_id", "year")) %>%
  mutate(primary_bloc = ifelse(is.na(primary_bloc), "none", primary_bloc))

bloc_counts <- py_bloc %>% count(primary_bloc, sort = TRUE)
cat("=== PART A: bloc fixed effects ===\n")
cat("Distinct bloc categories (incl. 'none'):", nrow(bloc_counts), "\n")
print(head(bloc_counts, 10))
# collapse blocs with too few person-years to identify their own FE into
# an "other_bloc" bucket, so a handful of tiny blocs don't produce
# unstable/undefined coefficients
small_blocs <- bloc_counts %>% filter(n < 20, primary_bloc != "none") %>% pull(primary_bloc)
py_bloc <- py_bloc %>% mutate(bloc_fe = ifelse(primary_bloc %in% small_blocs, "other_bloc", primary_bloc))
cat("Bloc categories after collapsing rare ones (<20 obs) into 'other_bloc':",
    n_distinct(py_bloc$bloc_fe), "\n\n")

base_rhs <- "z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin + z_trade_mean + z_dep_asym_max"

m_dis_nofe <- glm(as.formula(paste("event_dissolve ~", base_rhs)), data = py_bloc, family = binomial())
m_flip_nofe <- glm(as.formula(paste("event_flip ~", base_rhs)), data = py_bloc, family = binomial())
m_dis_fe <- glm(as.formula(paste("event_dissolve ~", base_rhs, "+ bloc_fe")), data = py_bloc, family = binomial())
m_flip_fe <- glm(as.formula(paste("event_flip ~", base_rhs, "+ bloc_fe")), data = py_bloc, family = binomial())

ct_dis_nofe <- coeftest(m_dis_nofe, vcov = vcovCL(m_dis_nofe, cluster = py_bloc$triad_id))
ct_flip_nofe <- coeftest(m_flip_nofe, vcov = vcovCL(m_flip_nofe, cluster = py_bloc$triad_id))
ct_dis_fe <- coeftest(m_dis_fe, vcov = vcovCL(m_dis_fe, cluster = py_bloc$triad_id))
ct_flip_fe <- coeftest(m_flip_fe, vcov = vcovCL(m_flip_fe, cluster = py_bloc$triad_id))

lrt_dis <- anova(m_dis_nofe, m_dis_fe, test = "Chisq")
lrt_flip <- anova(m_flip_nofe, m_flip_fe, test = "Chisq")

# ==================================================================
# PART B: H3, two-way (triad x bloc) clustering on the tie-choice model
# ==================================================================
h3 <- readRDS("results/h3_displacement.rds")
strict <- h3$strict  # one row per single-tie-change event, all 3 ties trade-valid

tagB <- tag_triads(strict %>% transmute(triad_id, year = t1) %>% distinct())
strict_b <- strict %>% mutate(event_id = row_number()) %>%
  left_join(tagB, by = c("triad_id", "t1" = "year")) %>%
  mutate(primary_bloc = ifelse(is.na(primary_bloc), "none", primary_bloc))

bloc_counts_events <- strict_b %>% count(primary_bloc, sort = TRUE)
cat("\n\n=== PART B: two-way clustering on H3 (event-level bloc tagging) ===\n")
cat("Realignment events by bloc:\n"); print(head(bloc_counts_events, 10))

# rebuild the long format WITH triad_id and primary_bloc carried along
long_b <- strict_b %>%
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

# clogit()'s default "exact" method has no score residuals, so sandwich
# cluster-robust SEs are unavailable for it. With exactly one event per
# stratum (one changed tie per event), "exact" and "efron" give IDENTICAL
# point estimates -- no tied event times exist to break differently -- so
# refitting with method="efron" for this diagnostic only is exact, not an
# approximation, and unlocks sandwich's cluster-robust/two-way machinery.
clog_primary <- coxph(Surv(rep(1, nrow(long_b)), changed) ~ z_unbal_load + z_dep_asym + z_emb + strata(event_id),
                       data = long_b, method = "efron")

cl2_h3 <- long_b[, c("triad_id", "primary_bloc")]
ct_h3_1way <- coeftest(clog_primary, vcov = vcovCL(clog_primary, cluster = long_b$triad_id))
ct_h3_2way <- tryCatch(
  coeftest(clog_primary, vcov = vcovCL(clog_primary, cluster = cl2_h3, multi0 = TRUE)),
  error = function(e) { cat("2-way clustering on clogit FAILED:", conditionMessage(e), "\n"); NULL }
)

# ------------------------------------------------------------------
sink("results/bloc_fe_and_h3_clustering_output.txt")
cat("=== PART A: H2 bloc fixed effects ===\n")
cat("Bloc categories (after collapsing rare ones):", n_distinct(py_bloc$bloc_fe), "\n\n")
cat("--- DISSOLUTION, no bloc FE ---\n"); print(ct_dis_nofe)
cat("\n--- DISSOLUTION, WITH bloc FE ---\n"); print(ct_dis_fe)
cat("\n--- LRT: bloc FE jointly significant (dissolution)? ---\n"); print(lrt_dis)
cat("\n\n--- REALIGNMENT, no bloc FE ---\n"); print(ct_flip_nofe)
cat("\n--- REALIGNMENT, WITH bloc FE ---\n"); print(ct_flip_fe)
cat("\n--- LRT: bloc FE jointly significant (realignment)? ---\n"); print(lrt_flip)

cat("\n\n=== PART B: H3 two-way (triad + bloc) clustering ===\n")
cat("Single-tie-change events used:", nrow(strict_b), "\n")
cat("\n--- 1-way (triad only) ---\n"); print(ct_h3_1way)
if (!is.null(ct_h3_2way)) { cat("\n--- 2-way (triad + bloc) ---\n"); print(ct_h3_2way) }
sink()

cat(readLines("results/bloc_fe_and_h3_clustering_output.txt"), sep = "\n")

saveRDS(list(py_bloc = py_bloc, m_dis_fe = m_dis_fe, m_flip_fe = m_flip_fe,
             ct_dis_fe = ct_dis_fe, ct_flip_fe = ct_flip_fe, lrt_dis = lrt_dis, lrt_flip = lrt_flip,
             strict_b = strict_b, long_b = long_b, clog_primary = clog_primary,
             ct_h3_1way = ct_h3_1way, ct_h3_2way = ct_h3_2way),
        "results/bloc_fe_and_h3_clustering.rds")
cat("\n\n07_bloc_fe_and_h3_clustering.R complete.\n")
