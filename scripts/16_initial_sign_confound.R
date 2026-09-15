# ------------------------------------------------------------------
# 16: Response to external review -- the single most important
#     remaining concern. 3,437 of 3,456 single-tie-change realignment
#     events (99.4%) occur in triads with exactly one negative tie
#     (balance requires an odd number of negatives, so 1 or 3); only 19
#     have three. In a one-negative-two-positive triad, a hostility ->
#     alliance resolution MECHANICALLY identifies the pre-existing
#     negative tie as "the one that changed" -- there is no real 3-way
#     choice for that transition direction. If negative ties are simply
#     more trade-asymmetric than positive ties for unrelated reasons,
#     the main result could be an artifact of tie-type composition
#     rather than a genuine choice mechanism.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(tidyr); library(survival); library(sandwich); library(lmtest) })


h3 <- readRDS("results/h3_displacement.rds")
s <- h3$strict %>% mutate(n_neg_t1 = (sign12_t1 == -1) + (sign13_t1 == -1) + (sign23_t1 == -1))
cat("Triad composition at t1 (number of negative ties before the change):\n")
print(table(s$n_neg_t1))

s <- s %>% mutate(
  trans12 = case_when(changed12 & sign12_t1 == 1 & sign12_t2 == -1 ~ "pos_to_neg", changed12 & sign12_t1 == -1 & sign12_t2 == 1 ~ "neg_to_pos"),
  trans13 = case_when(changed13 & sign13_t1 == 1 & sign13_t2 == -1 ~ "pos_to_neg", changed13 & sign13_t1 == -1 & sign13_t2 == 1 ~ "neg_to_pos"),
  trans23 = case_when(changed23 & sign23_t1 == 1 & sign23_t2 == -1 ~ "pos_to_neg", changed23 & sign23_t1 == -1 & sign23_t2 == 1 ~ "neg_to_pos"),
  transition = coalesce(trans12, trans13, trans23))

# --- diagnostic: is the pre-existing negative tie simply more trade-asymmetric? ---
one_neg <- s %>% filter(n_neg_t1 == 1) %>% rowwise() %>% mutate(
  neg_asym = case_when(sign12_t1 == -1 ~ dep_asymmetry_12, sign13_t1 == -1 ~ dep_asymmetry_13, sign23_t1 == -1 ~ dep_asymmetry_23),
  pos_asym_mean = mean(c(dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23)[c(sign12_t1, sign13_t1, sign23_t1) == 1], na.rm = TRUE)
) %>% ungroup()
diag_test <- t.test(one_neg$neg_asym, one_neg$pos_asym_mean, paired = TRUE)

# --- TEST A: pos_to_neg events, choice restricted to the 2 originally-positive ties
#     (the pre-existing negative tie cannot logically be "the one that changed to
#     negative" in a pos_to_neg event, so dropping it gives a genuine 2-way choice) ---
p2n <- s %>% filter(transition == "pos_to_neg", n_neg_t1 == 1)
long_p2n <- p2n %>% mutate(event_id = row_number()) %>%
  select(event_id, triad_id, sign12_t1, sign13_t1, sign23_t1, ul12, ul13, ul23, emb12, emb13, emb23,
         dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, changed12, changed13, changed23) %>%
  pivot_longer(cols = c(ul12, ul13, ul23, emb12, emb13, emb23, dep_asymmetry_12, dep_asymmetry_13,
                        dep_asymmetry_23, changed12, changed13, changed23),
               names_to = c(".value", "tie"), names_pattern = "(ul|emb|dep_asymmetry|changed)_?(12|13|23)") %>%
  left_join(p2n %>% mutate(event_id = row_number()) %>%
              select(event_id, sign12_t1, sign13_t1, sign23_t1) %>%
              pivot_longer(cols = c(sign12_t1, sign13_t1, sign23_t1), names_to = "tie2", values_to = "sign_t1") %>%
              mutate(tie2 = gsub("sign|_t1", "", tie2)),
            by = c("event_id", "tie" = "tie2")) %>%
  filter(sign_t1 == 1) %>%
  mutate(z_unbal_load = as.numeric(scale(ul)), z_dep_asym = as.numeric(scale(dep_asymmetry)), z_emb = as.numeric(scale(emb)))

surv_a <- Surv(rep(1, nrow(long_p2n)), long_p2n$changed)
m_a <- do.call(coxph, list(formula = surv_a ~ z_unbal_load + z_dep_asym + z_emb + strata(event_id), data = long_p2n, method = "efron"))
ct_a <- coeftest(m_a, vcov = vcovCL(m_a, cluster = long_p2n$triad_id))

p2n_check <- p2n %>% rowwise() %>% mutate(
  pos_asyms = list(c(dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23)[which(c(sign12_t1, sign13_t1, sign23_t1) == 1)]),
  changed_flags = list(c(changed12, changed13, changed23)[which(c(sign12_t1, sign13_t1, sign23_t1) == 1)])
) %>% ungroup()
p2n_check$higher_asym_changed <- sapply(seq_len(nrow(p2n_check)), function(i) {
  a <- p2n_check$pos_asyms[[i]]; cc <- p2n_check$changed_flags[[i]]
  if (length(a) != 2 || any(is.na(a)) || a[1] == a[2]) return(NA)
  which.max(a) == which(cc)
})
bt_a_naive <- binom.test(sum(p2n_check$higher_asym_changed, na.rm = TRUE), sum(!is.na(p2n_check$higher_asym_changed)), p = 0.5)

# --- TEST B: three-negative triads (n=19) -- genuine 3-way choice, no sign confound at all ---
three_neg <- s %>% filter(n_neg_t1 == 3) %>% rowwise() %>% mutate(
  changed_asym = case_when(changed12 ~ dep_asymmetry_12, changed13 ~ dep_asymmetry_13, changed23 ~ dep_asymmetry_23),
  max_asym = max(dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23, na.rm = TRUE),
  is_max = changed_asym == max_asym
) %>% ungroup()
bt_b <- binom.test(sum(three_neg$is_max, na.rm = TRUE), sum(!is.na(three_neg$is_max)), p = 1/3)

# --- MAIN TEST: initial_negative added directly to the full 3,456-event model ---
h3b <- readRDS("results/bloc_fe_and_h3_clustering.rds")
long <- h3b$long_b
sign_long <- s %>% mutate(event_id = row_number()) %>%
  select(event_id, sign12_t1, sign13_t1, sign23_t1) %>%
  pivot_longer(cols = c(sign12_t1, sign13_t1, sign23_t1), names_to = "tie", values_to = "sign_t1") %>%
  mutate(tie = gsub("sign|_t1", "", tie))
long2 <- long %>% left_join(sign_long, by = c("event_id", "tie")) %>%
  mutate(initial_negative = as.numeric(sign_t1 == -1))

surv <- Surv(rep(1, nrow(long2)), long2$changed)
m1 <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + z_emb + initial_negative + strata(event_id), data = long2, method = "efron"))
ct1 <- coeftest(m1, vcov = vcovCL(m1, cluster = long2$triad_id))
m2 <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym * initial_negative + z_emb + strata(event_id), data = long2, method = "efron"))
ct2 <- coeftest(m2, vcov = vcovCL(m2, cluster = long2$triad_id))

# --- DYAD-LEVEL CLUSTERING: the same state-pair (dyad) recurs across many triads;
#     cluster by the specific pair of states forming each candidate tie, not just
#     the triad as a whole ---
strict2 <- s %>% mutate(event_id = row_number(),
                         dyad12 = paste(pmin(node1, node2), pmax(node1, node2)),
                         dyad13 = paste(pmin(node1, node3), pmax(node1, node3)),
                         dyad23 = paste(pmin(node2, node3), pmax(node2, node3)))
dyad_long <- strict2 %>% select(event_id, dyad12, dyad13, dyad23) %>%
  pivot_longer(cols = c(dyad12, dyad13, dyad23), names_to = "tie", values_to = "dyad_id") %>%
  mutate(tie = gsub("dyad", "", tie))
long3 <- long2 %>% left_join(dyad_long, by = c("event_id", "tie"))
cat("\nUnique dyad_ids among candidate ties:", n_distinct(long3$dyad_id), "of", nrow(long3), "long rows",
    "| max repeats of one dyad:", max(table(long3$dyad_id)), "\n")

m0 <- do.call(coxph, list(formula = surv ~ z_unbal_load + z_dep_asym + z_emb + strata(event_id), data = long3, method = "efron"))
ct_triad_only <- coeftest(m0, vcov = vcovCL(m0, cluster = long3$triad_id))
ct_dyad       <- coeftest(m0, vcov = vcovCL(m0, cluster = long3$dyad_id))
cl3 <- long3[, c("triad_id", "dyad_id")]
ct_triad_dyad_2way <- coeftest(m0, vcov = vcovCL(m0, cluster = cl3, multi0 = TRUE))

# --- STATE-LEVEL CLUSTERING: dyad clustering catches the same exact pair
#     recurring, but not two DIFFERENT dyads sharing one common state
#     (e.g. a triad's US-UK and US-France relationships). Cluster by each
#     tie's lower- and higher-ccode state separately, and two-way. ---
strict3 <- strict2 %>% mutate(
  s1_12 = pmin(node1, node2), s2_12 = pmax(node1, node2),
  s1_13 = pmin(node1, node3), s2_13 = pmax(node1, node3),
  s1_23 = pmin(node2, node3), s2_23 = pmax(node2, node3))
state_low_long <- strict3 %>% select(event_id, s1_12, s1_13, s1_23) %>%
  pivot_longer(cols = c(s1_12, s1_13, s1_23), names_to = "tie", values_to = "state_low") %>%
  mutate(tie = gsub("s1_", "", tie))
state_high_long <- strict3 %>% select(event_id, s2_12, s2_13, s2_23) %>%
  pivot_longer(cols = c(s2_12, s2_13, s2_23), names_to = "tie", values_to = "state_high") %>%
  mutate(tie = gsub("s2_", "", tie))
long4 <- long3 %>% left_join(state_low_long, by = c("event_id", "tie")) %>%
  left_join(state_high_long, by = c("event_id", "tie"))
cat("\nUnique state_low:", n_distinct(long4$state_low), " state_high:", n_distinct(long4$state_high), "\n")

ct_state_low  <- coeftest(m0, vcov = vcovCL(m0, cluster = long4$state_low))
ct_state_high <- coeftest(m0, vcov = vcovCL(m0, cluster = long4$state_high))
cl4 <- long4[, c("state_low", "state_high")]
ct_state_2way <- coeftest(m0, vcov = vcovCL(m0, cluster = cl4, multi0 = TRUE))

sink("results/initial_sign_confound_output.txt")
cat("=== Triad composition (n_neg_t1) ===\n"); print(table(s$n_neg_t1))
cat("\n=== Diagnostic: is the pre-existing negative tie more trade-asymmetric than the positive ties? ===\n")
cat("Mean asymmetry, negative tie:", round(mean(one_neg$neg_asym, na.rm=TRUE),4),
    " | positive ties (mean):", round(mean(one_neg$pos_asym_mean, na.rm=TRUE),4), "\n")
print(diag_test)

cat("\n\n=== TEST A: pos_to_neg events, choice restricted to the 2 originally-positive ties (n=", nrow(p2n), " events) ===\n", sep="")
cat("--- Conditional logit ---\n"); print(ct_a)
cat("\n--- Naive rank comparison (of the 2 positive ties, did the more-asymmetric one change?) ---\n")
cat(round(mean(p2n_check$higher_asym_changed, na.rm=TRUE)*100,1), "% (n=", sum(!is.na(p2n_check$higher_asym_changed)), ", null=50%)\n")
print(bt_a_naive)

cat("\n\n=== TEST B: three-negative triads (n=", nrow(three_neg), " events) ===\n", sep="")
cat(round(mean(three_neg$is_max, na.rm=TRUE)*100,1), "% (null=33.3%)\n")
print(bt_b)

cat("\n\n=== MAIN TEST: full sample (n=3,456), + initial_negative control ===\n")
print(ct1)
cat("\n=== Interaction: asymmetry x initial_negative ===\n")
print(ct2)

cat("\n\n=== DYAD-LEVEL CLUSTERING (main baseline model, no initial_negative) ===\n")
cat("Unique dyads:", n_distinct(long3$dyad_id), " | max repeats:", max(table(long3$dyad_id)), "\n")
cat("--- triad-clustered ---\n"); print(ct_triad_only)
cat("\n--- dyad-clustered ---\n"); print(ct_dyad)
cat("\n--- 2-way (triad + dyad) clustered ---\n"); print(ct_triad_dyad_2way)

cat("\n\n=== STATE-LEVEL CLUSTERING (catches different dyads sharing one state) ===\n")
cat("--- clustered by state_low ---\n"); print(ct_state_low)
cat("\n--- clustered by state_high ---\n"); print(ct_state_high)
cat("\n--- 2-way (state_low + state_high) clustered ---\n"); print(ct_state_2way)
sink()

cat(readLines("results/initial_sign_confound_output.txt"), sep = "\n")

saveRDS(list(ct1 = ct1, ct2 = ct2, ct_a = ct_a, bt_a_naive = bt_a_naive, bt_b = bt_b,
             ct_triad_only = ct_triad_only, ct_dyad = ct_dyad, ct_triad_dyad_2way = ct_triad_dyad_2way,
             ct_state_low = ct_state_low, ct_state_high = ct_state_high, ct_state_2way = ct_state_2way,
             diag_test = diag_test),
        "results/initial_sign_confound.rds")
cat("\n\n16_initial_sign_confound.R complete.\n")
