# ------------------------------------------------------------------
# 19: Builds Table 1's progressive model sequence (baseline -> +
#     starting sign -> + weak-link mutual dependence -> + capability
#     asymmetry), so the starting-sign-adjusted specification, not the
#     unadjusted baseline, is the primary reported model.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest) })


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
  select(event_id, triad_id, sign12_t1, sign13_t1, sign23_t1,
         ul12, ul13, ul23, emb12, emb13, emb23, dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23,
         cap_asymmetry_12, cap_asymmetry_13, cap_asymmetry_23,
         mutual_min_12, mutual_min_13, mutual_min_23, changed12, changed13, changed23) %>%
  pivot_longer(cols = c(ul12,ul13,ul23,emb12,emb13,emb23,dep_asymmetry_12,dep_asymmetry_13,dep_asymmetry_23,
                        cap_asymmetry_12,cap_asymmetry_13,cap_asymmetry_23,
                        mutual_min_12,mutual_min_13,mutual_min_23,changed12,changed13,changed23),
               names_to = c(".value","tie"),
               names_pattern = "(ul|emb|dep_asymmetry|cap_asymmetry|mutual_min|changed)_?(12|13|23)") %>%
  left_join(strict %>% select(event_id, sign12_t1, sign13_t1, sign23_t1) %>%
              pivot_longer(cols = c(sign12_t1,sign13_t1,sign23_t1), names_to = "tie2", values_to = "sign_t1") %>%
              mutate(tie2 = gsub("sign|_t1","",tie2)),
            by = c("event_id","tie" = "tie2")) %>%
  mutate(initial_negative = as.numeric(sign_t1 == -1),
         z_unbal_load = as.numeric(scale(ul)), z_dep_asym = as.numeric(scale(dep_asymmetry)),
         z_emb = as.numeric(scale(emb)), z_cap_asym = as.numeric(scale(cap_asymmetry)),
         z_mutual_min = as.numeric(scale(mutual_min)))

fit <- function(rhs, data, label) {
  surv <- Surv(rep(1, nrow(data)), data$changed)
  m <- do.call(coxph, list(formula = as.formula(paste("surv ~", rhs, "+ strata(event_id)")), data = data, method = "efron"))
  ct <- coeftest(m, vcov = vcovCL(m, cluster = data$triad_id))
  cat("\n===", label, "=== N events =", n_distinct(data$event_id), "\n"); print(ct)
  invisible(ct)
}

sink("results/primary_table_progressive_output.txt")
m1 <- fit("z_unbal_load + z_dep_asym + z_emb", long, "M1: baseline")
long2 <- long %>% filter(!is.na(z_mutual_min))
m2 <- fit("z_unbal_load + z_dep_asym + initial_negative + z_emb", long2, "M2: + starting sign")
m3 <- fit("z_unbal_load + z_dep_asym + initial_negative + z_mutual_min + z_emb", long2, "M3: + weak-link mutual dependence")
long4 <- long2 %>% filter(!is.na(z_cap_asym))
m4 <- fit("z_unbal_load + z_dep_asym + initial_negative + z_mutual_min + z_cap_asym + z_emb", long4, "M4: + capability asymmetry (full model)")
sink()
cat(readLines("results/primary_table_progressive_output.txt"), sep = "\n")

saveRDS(list(m1=m1,m2=m2,m3=m3,m4=m4), "results/primary_table_progressive.rds")
cat("\n\n19_primary_table_progressive_models.R complete.\n")
