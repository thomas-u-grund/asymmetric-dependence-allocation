# ------------------------------------------------------------------
# 24: Figure 3 -- temporal-ordering coefficient plot. The M4
#     dependence-asymmetry coefficient measured contemporaneously (at
#     the event year, Table 1's preferred specification) and at
#     increasing lags before the triad becomes inconsistent (spell
#     onset, onset-1yr, onset-3yr, onset-5yr, Table 2), each shown with
#     95% CIs under both triad clustering and triad+bloc two-way
#     clustering. Point estimates are nearly flat across the lag
#     sequence while the two-way-clustered interval widens at longer
#     lags -- the paper's point that the effect is not an artifact of
#     contemporaneous measurement, without overstating precision at
#     the longest lag.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest); library(ggplot2) })


# --- rebuild the M4 (contemporaneous) sample exactly as in script 19, then
#     add primary_bloc (a per-event constant, confirmed) so a two-way
#     clustered SE can be computed for the contemporaneous point too. ---
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

long4 <- long %>% filter(!is.na(z_mutual_min), !is.na(z_cap_asym))

bloc_map <- readRDS("results/bloc_fe_and_h3_clustering.rds")$long_b %>%
  distinct(event_id, primary_bloc)
long4 <- long4 %>% left_join(bloc_map, by = "event_id")
cat("Contemporaneous M4 sample: N events =", n_distinct(long4$event_id),
    " (missing bloc:", sum(is.na(long4$primary_bloc)), " rows)\n")

surv <- Surv(rep(1, nrow(long4)), long4$changed)
m4_contemp <- do.call(coxph, list(
  formula = surv ~ z_unbal_load + z_dep_asym + initial_negative + z_mutual_min + z_cap_asym + z_emb + strata(event_id),
  data = long4, method = "efron"))
ct1_contemp <- coeftest(m4_contemp, vcov = vcovCL(m4_contemp, cluster = long4$triad_id))
ct2_contemp <- coeftest(m4_contemp, vcov = vcovCL(m4_contemp,
  cluster = data.frame(triad_id = long4$triad_id, primary_bloc = long4$primary_bloc), multi0 = TRUE))

cat("\nContemporaneous M4, recomputed here:\n")
cat("  triad-clustered:   "); print(ct1_contemp["z_dep_asym", ])
cat("  triad+bloc 2-way:  "); print(ct2_contemp["z_dep_asym", ])

# --- pull the lag results already fit in script 11 ---
lag <- readRDS("results/m4_temporal_endogeneity.rds")

rows <- list(
  list(label = "Event year\n(contemporaneous)", order = 1, ct1 = ct1_contemp, ct2 = ct2_contemp),
  list(label = "Spell onset",  order = 2, ct1 = lag$r0$ct1, ct2 = lag$r0$ct2),
  list(label = "Onset − 1 yr", order = 3, ct1 = lag$r1$ct1, ct2 = lag$r1$ct2),
  list(label = "Onset − 3 yr", order = 4, ct1 = lag$r3$ct1, ct2 = lag$r3$ct2),
  list(label = "Onset − 5 yr", order = 5, ct1 = lag$r5$ct1, ct2 = lag$r5$ct2)
)

df <- bind_rows(lapply(rows, function(r) {
  b <- r$ct1["z_dep_asym", "Estimate"]
  se1 <- r$ct1["z_dep_asym", "Std. Error"]
  se2 <- r$ct2["z_dep_asym", "Std. Error"]
  tibble(label = r$label, order = r$order, b = b,
         lo1 = b - 1.96 * se1, hi1 = b + 1.96 * se1,
         lo2 = b - 1.96 * se2, hi2 = b + 1.96 * se2,
         p1 = r$ct1["z_dep_asym", "Pr(>|z|)"], p2 = r$ct2["z_dep_asym", "Pr(>|z|)"])
}))
df$label <- factor(df$label, levels = df$label[order(df$order)])

write.csv(df, "results/figure3_data.csv", row.names = FALSE)
cat("\nFigure 3 data:\n"); print(df, digits = 3)

# two-way CI drawn wide (behind), one-way CI drawn narrow (in front), so both
# are visible without one occluding the other -- makes the widening-at-longer-
# lags point legible at a glance.
p <- ggplot(df, aes(x = label, y = b)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lo2, ymax = hi2), width = 0.12, linewidth = 0.7, color = "grey45") +
  geom_errorbar(aes(ymin = lo1, ymax = hi1), width = 0.12, linewidth = 1.1, color = "black") +
  geom_point(size = 3) +
  labs(x = NULL, y = "Dependence-asymmetry coefficient (M4, standardized)",
       title = "Dependence asymmetry measured before the triad\nbecomes inconsistent",
       subtitle = "Black bars: 95% CI, triad-clustered. Grey bars: 95% CI, triad+bloc two-way clustered.\nPoint estimates are stable across measurement timing; the two-way-clustered interval widens at longer lags.") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.text.x = element_text(size = 10))

ggsave("results/fig3_temporal_ordering.pdf", p, width = 8, height = 5.5)
ggsave("results/fig3_temporal_ordering.png", p, width = 8, height = 5.5, dpi = 300)

cat("\n24_figure3_temporal_ordering.R complete. Saved results/fig3_temporal_ordering.pdf, .png, and figure3_data.csv\n")
