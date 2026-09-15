# ------------------------------------------------------------------
# 27: Response to external review point 3 -- does trade-dependence
#     asymmetry simply proxy for a difference in the two states' overall
#     economic size (one large trading economy, one small one) rather
#     than a distinct relational mechanism? Section 6.2 already shows
#     dependence asymmetry survives controlling for the WEAK-LINK level
#     of mutual dependence between the two specific states in a tie.
#     This is a different, harsher test: a state-level measure of each
#     side's overall economic size, independent of anything about their
#     bilateral relationship. Using each state's total national trade
#     (imports + exports, Correlates of War Trade v4.0 national file --
#     already the denominator used to build the dependence-share
#     measure, so no new external dataset is required), construct an
#     economic-size-asymmetry control: the absolute log difference in
#     the two states' own total trade. Add it to M4 and check whether
#     dependence asymmetry survives.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest) })


national <- read_csv("raw/trade/COW_Trade_4.0/National_COW_4.0.csv", show_col_types = FALSE) %>%
  mutate(ccode = as.character(ccode),
         total_trade = ifelse(imports < 0 | exports < 0, NA, imports + exports)) %>%
  filter(!is.na(total_trade)) %>%
  select(ccode, year, total_trade)
cat("National trade coverage:", n_distinct(national$ccode), "states,",
    min(national$year), "-", max(national$year), "\n")

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

long <- strict %>%
  select(event_id, triad_id, sign12_t1, sign13_t1, sign23_t1, t1, node1, node2, node3,
         ul12, ul13, ul23, emb12, emb13, emb23, dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23,
         cap_asymmetry_12, cap_asymmetry_13, cap_asymmetry_23,
         mutual_min_12, mutual_min_13, mutual_min_23, changed12, changed13, changed23) %>%
  pivot_longer(cols = c(ul12,ul13,ul23,emb12,emb13,emb23,dep_asymmetry_12,dep_asymmetry_13,dep_asymmetry_23,
                        cap_asymmetry_12,cap_asymmetry_13,cap_asymmetry_23,
                        mutual_min_12,mutual_min_13,mutual_min_23,changed12,changed13,changed23),
               names_to = c(".value","tie"),
               names_pattern = "(ul|emb|dep_asymmetry|cap_asymmetry|mutual_min|changed)_?(12|13|23)") %>%
  mutate(n1 = case_when(tie=="12"~node1, tie=="13"~node1, tie=="23"~node2),
         n2 = case_when(tie=="12"~node2, tie=="13"~node3, tie=="23"~node3)) %>%
  left_join(strict %>% select(event_id, sign12_t1, sign13_t1, sign23_t1) %>%
              pivot_longer(cols = c(sign12_t1,sign13_t1,sign23_t1), names_to = "tie2", values_to = "sign_t1") %>%
              mutate(tie2 = gsub("sign|_t1","",tie2)),
            by = c("event_id","tie" = "tie2"))

# merge each candidate relationship's two states' own total national trade at t1
long <- long %>%
  left_join(national, by = c("n1" = "ccode", "t1" = "year")) %>% rename(total_trade_1 = total_trade) %>%
  left_join(national, by = c("n2" = "ccode", "t1" = "year")) %>% rename(total_trade_2 = total_trade)

cat("Rows with both states' total trade available:", sum(!is.na(long$total_trade_1) & !is.na(long$total_trade_2)),
    "of", nrow(long), "\n")

long <- long %>% mutate(
  econ_size_asym = abs(log(total_trade_1 + 1) - log(total_trade_2 + 1)),
  initial_negative = as.numeric(sign_t1 == -1),
  z_unbal_load = as.numeric(scale(ul)), z_dep_asym = as.numeric(scale(dep_asymmetry)),
  z_emb = as.numeric(scale(emb)), z_cap_asym = as.numeric(scale(cap_asymmetry)),
  z_mutual_min = as.numeric(scale(mutual_min)),
  z_econ_size_asym = as.numeric(scale(econ_size_asym))
)

# compare the paper's actual full M4 specification (starting sign + weak-
# link mutual dependence + capability asymmetry) with and without economic-
# size asymmetry added on top, on the identical sample
long_c <- long %>% filter(!is.na(z_dep_asym), !is.na(z_cap_asym), !is.na(z_mutual_min), !is.na(z_econ_size_asym))
cat("\nFinal sample: N events =", n_distinct(long_c$event_id), "\n")
cat("Correlation(dep_asymmetry, econ_size_asym):", round(cor(long_c$dep_asymmetry, long_c$econ_size_asym, use="complete.obs"), 3), "\n")

surv <- Surv(rep(1, nrow(long_c)), long_c$changed)
m_base <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + initial_negative + z_mutual_min + z_cap_asym + z_emb + strata(event_id), data = long_c, method = "efron"))
m_size <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + initial_negative + z_mutual_min + z_cap_asym + z_emb + z_econ_size_asym + strata(event_id), data = long_c, method = "efron"))
ct_base <- coeftest(m_base, vcov = vcovCL(m_base, cluster = long_c$triad_id))
ct_size <- coeftest(m_size, vcov = vcovCL(m_size, cluster = long_c$triad_id))

# diagnostics: econ_size_asym correlates 0.70 with cap_asymmetry (larger
# trading economies tend to be more materially capable) -- check whether
# the drop in dep_asymmetry's coefficient is driven by genuine shared
# explanatory content or by three-way multicollinearity among correlated
# regressors, by (a) dropping the possibly-redundant capability control,
# and (b) two-way (triad+bloc) clustering the full specification
bloc_map <- readRDS("results/bloc_fe_and_h3_clustering.rds")$long_b %>% distinct(event_id, primary_bloc)
long_c2 <- long_c %>% left_join(bloc_map, by = "event_id")
cat("\nMissing bloc after join:", sum(is.na(long_c2$primary_bloc)), "of", nrow(long_c2), "\n")

m_size_nocap <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + initial_negative + z_mutual_min + z_emb + z_econ_size_asym + strata(event_id), data = long_c, method = "efron"))
ct_size_nocap <- coeftest(m_size_nocap, vcov = vcovCL(m_size_nocap, cluster = long_c$triad_id))

ct_size_2way <- coeftest(m_size, vcov = vcovCL(m_size,
  cluster = data.frame(triad_id = long_c2$triad_id, primary_bloc = long_c2$primary_bloc), multi0 = TRUE))

sink("results/economic_size_confound_output.txt")
cat("=== Without economic-size-asymmetry control (= M4) ===\n"); print(ct_base)
cat("\n=== With economic-size-asymmetry control added to M4 ===\n"); print(ct_size)
cat("\n=== Same, two-way (triad+bloc) clustered ===\n"); print(ct_size_2way)
cat("\n=== With economic-size-asymmetry, capability asymmetry DROPPED (collinearity check) ===\n"); print(ct_size_nocap)
sink()
cat(readLines("results/economic_size_confound_output.txt"), sep = "\n")

saveRDS(list(m_base=m_base, m_size=m_size, m_size_nocap=m_size_nocap, ct_base=ct_base, ct_size=ct_size,
             ct_size_2way=ct_size_2way, ct_size_nocap=ct_size_nocap, n_events=n_distinct(long_c$event_id)),
        "results/economic_size_confound.rds")
cat("\n\n27_economic_size_confound.R complete.\n")
