# ------------------------------------------------------------------
# 12: Response to external review, priority item 11. Translate the
#     conditional-logit coefficient into a substantive quantity: how
#     often is the most trade-asymmetric relationship in a realigning
#     triad the one that actually changes, against the 1-in-3 baseline
#     random allocation would produce?
# ------------------------------------------------------------------
suppressMessages({ library(dplyr) })

h3 <- readRDS("results/h3_displacement.rds")
strict <- h3$strict

s <- strict %>%
  mutate(
    changed_asym = case_when(changed12 ~ dep_asymmetry_12, changed13 ~ dep_asymmetry_13, changed23 ~ dep_asymmetry_23),
    max_asym = pmax(dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, na.rm = TRUE),
    is_max = changed_asym == max_asym
  )
nvals <- function(a, b, c) length(unique(round(c(a, b, c), 8)))
s$n_distinct <- mapply(nvals, s$dep_asymmetry_12, s$dep_asymmetry_13, s$dep_asymmetry_23)

bt_naive <- binom.test(sum(s$is_max, na.rm = TRUE), sum(!is.na(s$is_max)), p = 1/3)
strict3 <- s %>% filter(n_distinct == 3)
bt_strict <- binom.test(sum(strict3$is_max), nrow(strict3), p = 1/3)

sink("results/substantive_magnitude_output.txt")
cat("N events:", nrow(s), "\n\n")
cat("=== Naive: changed tie has MAX dep_asymmetry, all events ===\n")
cat(sprintf("%.1f%% (null = 33.3%%)\n", 100 * mean(s$is_max, na.rm = TRUE)))
print(bt_naive)

cat("\nDistinctness of dep_asymmetry across a triad's 3 ties:\n")
cat("all 3 distinct:", round(mean(s$n_distinct == 3) * 100, 1), "%\n")
cat("exactly 2 distinct:", round(mean(s$n_distinct == 2) * 100, 1), "%\n")
cat("all equal:", round(mean(s$n_distinct == 1) * 100, 1), "%\n")

cat("\n=== PRIMARY: restricted to events with all 3 values distinct, n=", nrow(strict3), " ===\n", sep = "")
cat(sprintf("%.1f%% (null = 33.3%%)\n", 100 * mean(strict3$is_max)))
print(bt_strict)
sink()

cat(readLines("results/substantive_magnitude_output.txt"), sep = "\n")

saveRDS(list(s = s, strict3 = strict3, bt_naive = bt_naive, bt_strict = bt_strict),
        "results/substantive_magnitude.rds")
cat("\n\n12_substantive_magnitude.R complete.\n")
