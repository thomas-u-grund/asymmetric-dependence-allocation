# ------------------------------------------------------------------
# 25: Figure 4 -- H2, termination protection. Among alliances eligible
#     to terminate (restricted choice set, script 18 Part A), the more
#     asymmetrically dependent alliance is LESS likely to be the one
#     that dissolves. Uses the same delta-method, mirrored-percentile
#     construction as Figure 2 (script 22): two otherwise-comparable
#     alliances differing only in dependence asymmetry, focal at
#     percentile p vs. comparison at the mirror percentile 100-p.
#     (An earlier version of this figure also included a Panel A
#     reproducing Figure 2's H1 curve; dropped as a pure duplicate with
#     no added information -- Figure 2 already shows that result.)
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest); library(ggplot2) })


# ==================================================================
# Rebuild the Part-A eligible-choice-set H2 model (script 18), which
# is not saved as its own .rds -- only console text was saved there.
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

alliance_end <- diss_single %>% filter(gone_sign == 1) %>%
  mutate(n_alliances_t1 = (sign12_t1==1)+(sign13_t1==1)+(sign23_t1==1))

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
    filter(sign_t1 == 1) %>%
    mutate(z_dep_asym = as.numeric(scale(dep_asymmetry)), z_emb = as.numeric(scale(emb)))
}
long_elig <- build_eligible_long(alliance_end %>% filter(!is.na(dep_asymmetry_12), !is.na(dep_asymmetry_13), !is.na(dep_asymmetry_23)))

surv_e <- Surv(rep(1, nrow(long_elig)), long_elig$gone)
m_elig <- do.call(coxph, list(formula = surv_e ~ z_dep_asym + z_emb + strata(event_id), data = long_elig, method = "efron"))
ct_elig <- coeftest(m_elig, vcov = vcovCL(m_elig, cluster = long_elig$triad_id))
b_h2 <- ct_elig["z_dep_asym", "Estimate"]; se_h2 <- ct_elig["z_dep_asym", "Std. Error"]
cat("H2 eligible-choice-set dependence-asymmetry coefficient:", round(b_h2,4), " SE:", round(se_h2,4),
    " p:", round(ct_elig["z_dep_asym","Pr(>|z|)"],4), "\n")
cat("N events in eligible-choice-set model:", n_distinct(long_elig$event_id), "\n")

mu2 <- mean(long_elig$dep_asymmetry); sdv2 <- sd(long_elig$dep_asymmetry)
pctiles2 <- quantile(long_elig$dep_asymmetry, seq(0.10, 0.90, by = 0.01))
# mirrored-percentile comparison, matching Panel A's construction (script 22):
# focal at percentile p vs. comparison at percentile 100-p.
pctiles2_mirror <- quantile(long_elig$dep_asymmetry, rev(seq(0.10, 0.90, by = 0.01)))
z_focal2 <- (pctiles2 - mu2) / sdv2
z_comp2 <- (pctiles2_mirror - mu2) / sdv2
z_diff2 <- z_focal2 - z_comp2
logit2 <- b_h2 * z_diff2
se_logit2 <- abs(z_diff2) * se_h2
p_hat2 <- plogis(logit2); p_lo2 <- plogis(logit2 - 1.96*se_logit2); p_hi2 <- plogis(logit2 + 1.96*se_logit2)

df_h2 <- tibble(percentile = seq(10, 90, by = 1), p_hat = p_hat2, p_lo = p_lo2, p_hi = p_hi2)
write.csv(df_h2, "results/figure4_data.csv", row.names = FALSE)

marker_df <- df_h2 %>% filter(percentile %in% c(10, 90))

# Single panel only -- Panel A duplicated Figure 2 (H1) with no added
# information; this figure now shows H2 alone. Title kept non-causal
# ("associated with") to match the paper's own inferential language.
p <- ggplot(df_h2, aes(x = percentile, y = p_hat)) +
  geom_ribbon(aes(ymin = p_lo, ymax = p_hi), fill = "grey70", alpha = 0.5) +
  geom_line(linewidth = 1) +
  geom_point(data = marker_df, size = 2.5) +
  geom_text(data = marker_df, aes(label = sprintf("%.3f", p_hat)), vjust = -1, size = 3.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  labs(x = "Percentile of trade-dependence asymmetry (focal alliance)",
       y = "Predicted probability focal alliance is the one that terminates",
       title = "Dependence asymmetry is associated with\nprotection from alliance termination") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(size = 15))

ggsave("results/fig4_termination_protection.pdf", p, width = 7, height = 5)
ggsave("results/fig4_termination_protection.png", p, width = 7, height = 5, dpi = 300)

cat("\n25_figure4_reorient_vs_terminate.R complete. Saved results/fig4_termination_protection.pdf, .png, and figure4_data.csv\n")
