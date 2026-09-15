# Scripts

Run in numeric order from the `replication_pack/` root
(`for f in scripts/*.R; do Rscript "$f"; done`). Each script only depends
on outputs from earlier-numbered scripts.

| # | Script | What it does | Maps to |
|---|---|---|---|
| 01 | `build_dyad_year.R` | Builds the signed dyad-year panel from Correlates of War alliance/MID data, plus directed trade-dependence shares and capability asymmetry | Section 5.1–5.2 |
| 02 | `triads_and_load.R` | Enumerates closed triads per year, classifies balance, computes unbalanced load | Section 5.1, unbalanced load |
| 03 | `person_year_panel.R` | Builds the person-year competing-risks panel (spell construction, duration bins) | Section 5.3 |
| 04 | `hazard_models.R` | Triad-level dissolution/realignment hazard models as a function of trade dependence | Section 6.5, 6.8 |
| 05 | `h3_displacement.R` | Identifies single-tie realignment events; builds the tie-choice conditional-logit sample | Table 1 construction, Section 6.1 |
| 06 | `robustness.R` | Contiguity/regime-type controls and two-way clustering on the triad-level hazard result | Section 6.8 background |
| 07 | `bloc_fe_and_h3_clustering.R` | Bloc fixed effects and two-way (triad+bloc) clustering | Appendix A; Section 6.7/Appendix C.2 groundwork |
| 08 | `cross_tie_coercion.R` | Tests whether asymmetry in one tie predicts change in a *different* tie | Section 6.6 (discriminating implication) |
| 09 | `cap_asym_bloc_fe.R` | Bloc-FE robustness for the capability-asymmetry hazard result | Section 7.2 background |
| 10 | `mutual_dependence_and_capability.R` | Constructs mutual dependence (mean and weak-link) alongside asymmetry | Section 6.2 |
| 11 | `temporal_endogeneity.R` | Re-measures covariates at spell onset and 1/3/5 years before onset | Section 6.4, Table 2, Figure 3 |
| 12 | `substantive_magnitude.R` | Rank-comparison translation of the M4 coefficient (52.1% figure) | Section 6.1 |
| 13 | `bloc_and_era_robustness.R` | Leave-one-bloc-out and pre/post-1945 robustness | Appendix C.1 |
| 14 | `defense_only_network.R` | Independently rebuilds the network restricting alliances to defense pacts only | Appendix C.1 |
| 15 | `transition_split_and_edge_dissolution.R` | Splits realignment by direction (hostility-to-alliance vs. reverse); builds the edge-level dissolution model | Section 6.3; Section 7.3–7.4 |
| 16 | `initial_sign_confound.R` | Tests whether starting sign explains the main result | Section 6.3 (decisive test) |
| 17 | `floor_and_mid_robustness.R` | Sensitivity to the $50m trade-reporting floor and the dispute-hostility threshold | Appendix C.3, C.4 |
| 18 | `h2_eligible_choiceset_and_overlap_exclusion.R` | Restricts H2's choice set to eligible alliances; excludes the alliance/dispute overlap sample | Section 7.3; Appendix C.4 |
| 19 | `primary_table_progressive_models.R` | Builds Table 1's progressive M1–M4 sequence | Table 1 |
| 20 | `duration_control_and_m4_consistency.R` | Adds a relationship-duration control; reruns temporal-endogeneity checks on M4 | Appendix C.5; Table 2 |
| 21 | `contiguity_check.R` | Adds dyadic contiguity to M4 | Appendix C.6 |
| 22 | `figure2_main_effect.R` | Builds Figure 2 (predicted probability, 10th–90th percentile) | Figure 2 |
| 23 | `figure1_triad_diagram.R` | Builds Figure 1 (conceptual triad diagram) | Figure 1 |
| 24 | `figure3_temporal_ordering.R` | Builds Figure 3 (coefficient stability across measurement timing) | Figure 3 |
| 25 | `figure4_reorient_vs_terminate.R` | Builds Figure 4 (H2 termination-protection predicted probability) | Figure 4 |
| 26 | `find_historical_examples.R` | Mines the event data for concrete historical illustrations of the modal realignment event | Section 5.3 (Harlan County, Trujillo, Saudi–Yemen cases) |
| 27 | `economic_size_confound.R` | Exploratory check of whether dependence asymmetry proxies for state-level economic size; not referenced in the manuscript (see script header for why) | Not cited in main text |
| 28 | `full_multiway_clustering.R` | First attempt at combining triad, bloc, dyad, and state clustering; the triad+bloc+dyad portion is valid and reported, the state-clustering portion (lower/higher-ccode split) is superseded by script 30 | Appendix C.2, Table (triad+bloc+dyad rows) |
| 29 | `actor_overlap_clustering.R` | Attempted fix via connected-component clustering; identified as over-correcting (transitive closure) on further review and **not used** in the paper — kept on disk for transparency | Appendix C.2 (discussed, not the reported estimator) |
| 30 | `dyadic_robust_clustering.R` | Correct actor-overlap (dyadic-robust) variance estimator, validated against `vcovCL`'s one-way dyad clustering; this is the estimator reported in the paper | Appendix C.2 (actor-overlap row), Section 6.7, Section 8 |
