# ------------------------------------------------------------------
# 26: Mine the actual triad data for concrete historical examples of
#     the modal "realignment" event -- a militarized dispute subsiding
#     within a standing multilateral alliance framework -- to ground
#     Sections 3/5.3's theoretical language in real cases. Candidates
#     only; the cases used in the manuscript (Section 5.3) were
#     independently verified against historical sources afterward.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr) })


h3 <- readRDS("results/h3_displacement.rds")
strict <- h3$strict %>% mutate(event_id = row_number())
# dep_asymmetry_12/13/23 already present in h3$strict -- no re-merge needed

# restrict to the trade-valid sample (all three ties have trade data),
# hostility-to-alliance direction (changed tie was negative, becomes positive)
valid <- strict %>% filter(!is.na(dep_asymmetry_12), !is.na(dep_asymmetry_13), !is.na(dep_asymmetry_23))

changed_was_neg <- valid %>% mutate(
  changed_sign_t1 = case_when(changed12 ~ sign12_t1, changed13 ~ sign13_t1, changed23 ~ sign23_t1),
  changed_asym = case_when(changed12 ~ dep_asymmetry_12, changed13 ~ dep_asymmetry_13, changed23 ~ dep_asymmetry_23),
  max_asym = pmax(dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23),
  changed_is_max = changed_asym == max_asym,
  changed_n1 = case_when(changed12 ~ node1, changed13 ~ node1, changed23 ~ node2),
  changed_n2 = case_when(changed12 ~ node2, changed13 ~ node3, changed23 ~ node3)
) %>% filter(changed_sign_t1 == -1, event == "flipped")

cat("Hostility-to-alliance events (trade-valid):", nrow(changed_was_neg), "\n")
cat("Of these, changed tie has the max asymmetry:", sum(changed_was_neg$changed_is_max), "\n\n")

# candidates: changed tie is the max-asymmetry tie (matches H1 exactly),
# spread across different decades for variety
candidates <- changed_was_neg %>% filter(changed_is_max) %>%
  select(triad_id, t1, changed_n1, changed_n2, changed_asym, max_asym,
         node1, node2, node3, emb12, emb13, emb23) %>%
  arrange(desc(changed_asym))

states <- read_csv("raw/states/States2024/statelist2024.csv", show_col_types = FALSE) %>%
  mutate(ccode = as.character(ccode)) %>% distinct(ccode, .keep_all = TRUE) %>%
  select(ccode, statenme)

named <- candidates %>%
  left_join(states, by = c("changed_n1" = "ccode")) %>% rename(name1 = statenme) %>%
  left_join(states, by = c("changed_n2" = "ccode")) %>% rename(name2 = statenme) %>%
  left_join(states, by = c("node1" = "ccode")) %>% rename(nname1 = statenme) %>%
  left_join(states, by = c("node2" = "ccode")) %>% rename(nname2 = statenme) %>%
  left_join(states, by = c("node3" = "ccode")) %>% rename(nname3 = statenme)

# spread across decades for variety; print full triad membership so each
# candidate's context (the two OTHER, unchanged relationships) is visible
out <- named %>%
  mutate(decade = floor(t1/10)*10) %>%
  group_by(decade) %>% slice_max(changed_asym, n = 2, with_ties = FALSE) %>% ungroup() %>%
  select(t1, name1, name2, changed_asym, triad_third = nname3, nname1, nname2, nname3, triad_id) %>%
  arrange(t1)

print(out, n = 100)
write.csv(named %>% arrange(desc(changed_asym)) %>% head(60),
          "results/historical_candidates_top60.csv", row.names = FALSE)
write.csv(out, "results/historical_candidates_by_decade.csv", row.names = FALSE)
cat("\nSaved results/historical_candidates_top60.csv and _by_decade.csv\n")
