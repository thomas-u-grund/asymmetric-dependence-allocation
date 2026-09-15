# ------------------------------------------------------------------
# 29: Corrects a flaw in script 28's "state clustering": splitting
#     each tie into a "lower-ccode state" and "higher-ccode state"
#     column and clustering on those two columns separately is an
#     artifact of numeric ccode ordering, not a real grouping. The
#     same state (say, France) sits in the
#     lower-ccode column for some ties and the higher-ccode column for
#     others, so two observations that both involve France are only
#     recognized as dependent if France happens to occupy the SAME
#     column in both -- understating the true actor-overlap structure
#     script 28's Appendix C.2 table claimed to capture.
#
#     Fix: cluster on the state's IDENTITY, not its column position.
#     Two ties are correlated if they share ANY constituent state,
#     regardless of position. Since every tie is, by construction, an
#     edge between exactly two states, the graph formed by these edges
#     partitions the states into connected components -- and every
#     tie's two states necessarily fall in the same component (they
#     are joined by the tie itself). Assigning each observation to its
#     states' connected-component ID therefore gives a clean partition
#     (unlike the raw multi-membership structure), so it can be used
#     directly as one more clustering dimension in the same
#     Cameron-Gelbach-Miller multiway sandwich (vcovCL) machinery
#     already used throughout this project, alongside triad and bloc.
#     This also subsumes dyad-level clustering: two observations
#     sharing a dyad necessarily share both states, hence the same
#     component, so a separate dyad dimension is redundant once this
#     is included.
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest); library(igraph)
})


prev <- readRDS("results/full_multiway_clustering.rds")
long4 <- prev$long4
m4 <- prev$m4

cat("Reusing M4 sample from script 28: N events =", n_distinct(long4$event_id),
    " N tie-rows =", nrow(long4), "\n\n")

# ------------------------------------------------------------------
# Build the state-overlap graph and its connected components
# ------------------------------------------------------------------
edges <- long4 %>% transmute(a = as.character(state_lo), b = as.character(state_hi)) %>% distinct()
g <- graph_from_data_frame(edges, directed = FALSE)
comp <- components(g)
cat("States (graph nodes):", vcount(g), "\n")
cat("Connected components:", comp$no, "\n")
cat("Component size distribution (number of STATES per component):\n")
print(sort(table(comp$csize), decreasing = TRUE))

comp_id <- tibble(state = names(comp$membership), state_component = as.integer(comp$membership))

long4 <- long4 %>%
  mutate(state_lo = as.character(state_lo), state_hi = as.character(state_hi)) %>%
  left_join(comp_id, by = c("state_lo" = "state"))

# sanity check: every tie's two states must fall in the same component,
# since the tie itself is an edge joining them
check <- long4 %>% left_join(comp_id, by = c("state_hi" = "state"), suffix = c("_lo","_hi")) %>%
  summarise(mismatch = sum(state_component_lo != state_component_hi))
cat("\nSanity check -- ties whose two states land in different components (should be 0):", check$mismatch, "\n")

n_events_per_component <- long4 %>% distinct(event_id, state_component) %>% count(state_component, sort = TRUE)
cat("\nEvents per component (top 10):\n"); print(head(n_events_per_component, 10))
cat("Number of components with >=1 event:", n_distinct(long4$state_component), "\n")
cat("Share of events in the single largest component:",
    round(100 * max(n_events_per_component$n) / n_distinct(long4$event_id), 1), "%\n\n")

# ------------------------------------------------------------------
# Re-run clustering with the corrected, position-agnostic state
# grouping in place of script 28's lower/higher-ccode columns
# ------------------------------------------------------------------
ct_1way    <- coeftest(m4, vcov = vcovCL(m4, cluster = long4$triad_id))
ct_2way    <- coeftest(m4, vcov = vcovCL(m4, cluster = long4[, c("triad_id","primary_bloc")], multi0 = TRUE))
ct_3way_sc <- coeftest(m4, vcov = vcovCL(m4, cluster = long4[, c("triad_id","primary_bloc","state_component")], multi0 = TRUE))

sink("results/actor_overlap_clustering_output.txt")
cat("=== Corrected actor-overlap (connected-component) state clustering on M4 ===\n")
cat("N events =", n_distinct(long4$event_id), " N tie-rows =", nrow(long4), "\n")
cat("Distinct states:", vcount(g), " connected components:", comp$no, "\n")
cat("Sanity check (tie's two states in different components; should be 0):", check$mismatch, "\n")
cat("Share of events in the single largest component:",
    round(100 * max(n_events_per_component$n) / n_distinct(long4$event_id), 1), "%\n\n")
cat("--- 1-way: triad ---\n"); print(ct_1way)
cat("\n--- 2-way: triad + bloc ---\n"); print(ct_2way)
cat("\n--- 3-way: triad + bloc + state-connected-component (subsumes dyad and both state-position checks) ---\n")
print(ct_3way_sc)
sink()
cat(readLines("results/actor_overlap_clustering_output.txt"), sep = "\n")

saveRDS(list(long4 = long4, comp = comp, ct_1way = ct_1way, ct_2way = ct_2way, ct_3way_sc = ct_3way_sc),
        "results/actor_overlap_clustering.rds")
cat("\n\n29_actor_overlap_clustering.R complete.\n")
