# ------------------------------------------------------------------
# 01: Build signed dyad-year panel from Correlates of War data, plus
# asymmetric trade-dependence and capability-asymmetry measures.
#
# Adapted from ../2025 Social Balance/scripts/01_build_dyad_year.R.
# Same alliance/MID signed-tie construction; the trade block is changed
# to KEEP trade_share_low/trade_share_high separately instead of
# collapsing them into a single symmetric geometric-mean dependence
# score. That symmetric score (trade_dependence) is still produced,
# for continuity with the original project's control variable, but the
# asymmetry measures are the point of this project:
#   - dep_exposed  = max(share_low, share_high)   how exposed is the
#                    more-dependent side (Hirschman "vulnerability")
#   - dep_asymmetry = |share_low - share_high|     how lopsided the
#                    dependence is between the two sides
#   - dep_ratio    = max(share)/min(share)         same idea, scale-free
# ------------------------------------------------------------------
library(dplyr)
library(readr)
library(tidyr)

raw <- "raw"

# 1. Alliances -> positive ties (defense, neutrality, nonaggression, entente)
alliances <- read_csv(file.path(raw, "alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv"),
                       show_col_types = FALSE) %>%
  mutate(ccode_low = pmin(ccode1, ccode2), ccode_high = pmax(ccode1, ccode2)) %>%
  filter(defense == 1 | neutrality == 1 | nonaggression == 1 | entente == 1) %>%
  mutate(sign = 1) %>%
  select(year, ccode_low, ccode_high, sign) %>%
  distinct()

# 2. MIDs -> negative ties (hostility level >= 3: display of force or higher)
mids <- read_csv(file.path(raw, "mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv"),
                  show_col_types = FALSE) %>%
  mutate(ccode_low = pmin(statea, stateb), ccode_high = pmax(statea, stateb)) %>%
  filter(hihost >= 3) %>%
  mutate(sign = -1) %>%
  select(year, ccode_low, ccode_high, sign) %>%
  distinct()

# 3. Merge: if a dyad-year has both an alliance and a MID, code as negative
#    (conflict dominates cooperation in the same year -- conservative choice,
#    matches the parent project's convention, DECISIONS.md there)
signed_ties <- bind_rows(alliances, mids) %>%
  group_by(year, ccode_low, ccode_high) %>%
  summarise(sign = ifelse(any(sign == -1), -1, 1),
            multiplex = n_distinct(sign) > 1, .groups = "drop")

cat("Signed ties built:", nrow(signed_ties), "\n")
cat("  positive:", sum(signed_ties$sign == 1), " negative:", sum(signed_ties$sign == -1), "\n")
cat("  multiplex (both alliance+MID same year):", sum(signed_ties$multiplex), "\n")

# 4. State system membership -> restrict universe to actual states in that year
states <- read_csv(file.path(raw, "states/States2024/statelist2024.csv"), show_col_types = FALSE)
state_years <- states %>%
  rowwise() %>%
  mutate(yrs = list(seq(styear, endyear))) %>%
  ungroup() %>%
  select(ccode, yrs) %>%
  unnest(yrs) %>%
  rename(year = yrs) %>%
  distinct()

write_csv(signed_ties, "data/signed_ties.csv")
write_csv(state_years, "data/state_years.csv")

# ------------------------------------------------------------------
# 5. Trade dependence, kept ASYMMETRIC (1870-2014 coverage only, per
#    COW Trade 4.0 -- this bounds this project's usable sample to
#    roughly 1870-2012 once intersected with the alliance/MID window,
#    versus the parent project's full 1816-2012).
# ------------------------------------------------------------------
dyadic_trade <- read_csv(file.path(raw, "trade/COW_Trade_4.0/Dyadic_COW_4.0.csv"), show_col_types = FALSE) %>%
  mutate(across(c(flow1, flow2, smoothtotrade), ~ na_if(., -9))) %>%
  mutate(ccode_low = pmin(ccode1, ccode2), ccode_high = pmax(ccode1, ccode2)) %>%
  group_by(year, ccode_low, ccode_high) %>%
  summarise(total_dyadic_trade = sum(smoothtotrade, na.rm = TRUE), .groups = "drop")

national_trade <- read_csv(file.path(raw, "trade/COW_Trade_4.0/National_COW_4.0.csv"), show_col_types = FALSE) %>%
  mutate(total_trade = imports + exports,
         # same floor as the parent project: require >= $50m recorded
         # national trade before using it as a denominator, else shares
         # from tiny/misreported totals blow up above 1
         total_trade = ifelse(total_trade < 50, NA, total_trade)) %>%
  select(ccode, year, total_trade)

# 6. National capability (CINC), for the capability-asymmetry side of H4
nmc <- read_csv(file.path(raw, "nmc/NMCv7/abridged/NMC-70-abridged.csv"), show_col_types = FALSE) %>%
  select(ccode, year, cinc)

tie_trade <- dyadic_trade %>%
  left_join(national_trade, by = c("ccode_low" = "ccode", "year")) %>%
  rename(total_trade_low = total_trade) %>%
  left_join(national_trade, by = c("ccode_high" = "ccode", "year")) %>%
  rename(total_trade_high = total_trade) %>%
  left_join(nmc, by = c("ccode_low" = "ccode", "year")) %>%
  rename(cinc_low = cinc) %>%
  left_join(nmc, by = c("ccode_high" = "ccode", "year")) %>%
  rename(cinc_high = cinc) %>%
  mutate(
    overflow = total_dyadic_trade > pmin(total_trade_low, total_trade_high, na.rm = TRUE),
    trade_share_low  = ifelse(overflow, NA, total_dyadic_trade / total_trade_low),
    trade_share_high = ifelse(overflow, NA, total_dyadic_trade / total_trade_high),

    # symmetric measure, kept for continuity with the parent project's control
    trade_dependence = ifelse(!is.na(trade_share_low) & !is.na(trade_share_high),
                               sqrt(pmax(trade_share_low, 0) * pmax(trade_share_high, 0)), NA),
    log_trade_dependence = ifelse(!is.na(trade_dependence) & trade_dependence > 0,
                                   log(trade_dependence), NA),

    # asymmetry measures -- the point of this project
    dep_exposed    = ifelse(!is.na(trade_share_low) & !is.na(trade_share_high),
                             pmax(trade_share_low, trade_share_high), NA),
    dep_sheltered  = ifelse(!is.na(trade_share_low) & !is.na(trade_share_high),
                             pmin(trade_share_low, trade_share_high), NA),
    dep_asymmetry  = dep_exposed - dep_sheltered,
    dep_ratio      = ifelse(dep_sheltered > 0, dep_exposed / dep_sheltered, NA),
    log_dep_ratio  = ifelse(!is.na(dep_ratio) & dep_ratio > 0, log(dep_ratio), NA),
    # which side (low/high ccode) is the exposed one -- needed later to link
    # dependence-asymmetry to a SPECIFIC state, not just the dyad in the abstract
    exposed_is_low = trade_share_low >= trade_share_high,

    # capability asymmetry (Hirschman/dependency-theory companion to H4):
    # is the trade-exposed side also the materially weaker side?
    log_cinc_low   = log(pmax(cinc_low, 1e-6)),
    log_cinc_high  = log(pmax(cinc_high, 1e-6)),
    cap_asymmetry  = abs(log_cinc_low - log_cinc_high),
    # TRUE if the trade-exposed state is also the lower-capability state
    exposed_is_weaker = ifelse(exposed_is_low, cinc_low < cinc_high, cinc_high < cinc_low),

    # ABSOLUTE size of the dyadic trade relationship, independent of either
    # side's total trade. Needed to check whether dep_asymmetry is picking
    # up genuine economic entanglement or is partly a small-denominator
    # artifact (a tiny economy shows a high trade_share with almost any
    # partner regardless of the dollar value of the relationship).
    log_trade_volume = ifelse(!is.na(total_dyadic_trade) & total_dyadic_trade > 0,
                               log(total_dyadic_trade), NA)
  ) %>%
  select(year, ccode_low, ccode_high,
         trade_dependence, log_trade_dependence,
         trade_share_low, trade_share_high,
         dep_exposed, dep_sheltered, dep_asymmetry, dep_ratio, log_dep_ratio, exposed_is_low,
         cinc_low, cinc_high, cap_asymmetry, exposed_is_weaker,
         total_dyadic_trade, log_trade_volume,
         overflow)

n_overflow <- sum(tie_trade$overflow, na.rm = TRUE)
cat("Trade dyad-years:", nrow(tie_trade), " | overflow (excluded, share>1):", n_overflow,
    sprintf("(%.2f%%)\n", 100 * n_overflow / nrow(tie_trade)))

n_valid <- sum(!is.na(tie_trade$dep_asymmetry))
cat("Dyad-years with valid asymmetry measure:", n_valid,
    sprintf("(%.1f%% of all trade dyad-years)\n", 100 * n_valid / nrow(tie_trade)))

write_csv(tie_trade, "data/tie_trade.csv")

cat("01_build_dyad_year.R complete.\n")
