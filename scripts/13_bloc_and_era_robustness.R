# ------------------------------------------------------------------
# 13: Response to external review, priority item 4 (partial): leave-
#     one-bloc-out and pre/post-1945 robustness for the main tie-choice
#     result. Uses the bloc-tagged sample already built in script 07.
#
#     Bloc identification (via raw/alliances member lists, not assumed
#     from memory): the two largest blocs in the trade-valid sample are
#     NOT NATO/Warsaw Pact -- they are the Arab League (version4id 199,
#     1,228 of 3,456 events) and the Rio Pact/inter-American system
#     (version4id 210, 1,008 events). NATO (227) has 355 events; the
#     Warsaw Pact (243) has only 18. This corrects the "NATO, Warsaw
#     Pact" framing used loosely in the manuscript's data section and
#     should be fixed there.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest) })


b <- readRDS("results/bloc_fe_and_h3_clustering.rds")
long <- b$long_b  # event_id, triad_id, primary_bloc, z_unbal_load, z_dep_asym, z_emb, changed
strict_b <- b$strict_b  # for t1 (year), to build the era split

# attach t1 (year of realignment) to long via event_id -> strict_b mapping
year_lookup <- strict_b %>% mutate(event_id = row_number()) %>% select(event_id, t1)
long <- long %>% left_join(year_lookup, by = "event_id")

cat("Bloc sizes (events, i.e. 1/3 of long rows):\n")
print(strict_b %>% count(primary_bloc, sort = TRUE) %>% mutate(pct = round(100*n/sum(n),1)))

BLOC_NAMES <- c("199" = "Arab League", "210" = "Rio Pact / inter-American",
                "227" = "NATO", "243" = "Warsaw Pact")

fit_and_report <- function(data, label) {
  n_ev <- n_distinct(data$event_id)
  surv <- Surv(rep(1, nrow(data)), data$changed)
  m <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + z_emb + strata(event_id),
                           data = data, method = "efron"))
  ct <- coeftest(m, vcov = vcovCL(m, cluster = data$triad_id))
  cat(sprintf("%-40s N events=%5d  b=%7.4f  SE=%6.4f  p=%8.5f\n",
              label, n_ev, ct["z_dep_asym","Estimate"], ct["z_dep_asym","Std. Error"], ct["z_dep_asym","Pr(>|z|)"]))
  invisible(ct)
}

sink("results/bloc_and_era_robustness_output.txt")
cat("=== Leave-one-bloc-out (dependence-asymmetry coefficient, triad-clustered SEs) ===\n")
fit_and_report(long, "Full sample (all blocs)")
for (id in c("199", "210", "227", "243")) {
  sub <- long %>% filter(primary_bloc != id)
  fit_and_report(sub, sprintf("Excluding %s (id %s)", BLOC_NAMES[id], id))
}
sub_nonbloc <- long %>% filter(primary_bloc == "none")
fit_and_report(sub_nonbloc, "ONLY triads outside any big bloc")

cat("\n\n=== Era split ===\n")
fit_and_report(long %>% filter(t1 < 1945), "Pre-1945")
fit_and_report(long %>% filter(t1 >= 1945), "1945-present")
sink()

cat(readLines("results/bloc_and_era_robustness_output.txt"), sep = "\n")
cat("\n\n13_bloc_and_era_robustness.R complete.\n")
