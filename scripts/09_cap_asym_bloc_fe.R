# ------------------------------------------------------------------
# 09: Close the open item flagged in the outline (5.4): does the
#     capability (CINC) asymmetry effect on the hazard models survive
#     the same bloc fixed-effects stress test that overturned the
#     dep_asym_max triad-level "suppression" result in script 07?
#     Reuses the bloc-tagged person-year panel already built there.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(sandwich); library(lmtest) })


b <- readRDS("results/bloc_fe_and_h3_clustering.rds")
py_bloc <- b$py_bloc  # already has bloc_fe (collapsed rare blocs), z_trade_mean, z_dep_asym_max

# cap_asym_max is already present in py_bloc (carried through from
# person_year_risk_table via py_trade in script 07) -- no need to rejoin.
py_bloc <- py_bloc %>% mutate(z_cap_asym_max = as.numeric(scale(cap_asym_max)))

cat("Sample:", nrow(py_bloc), "person-years,", n_distinct(py_bloc$triad_id), "triads\n")
cat("Bloc FE categories:", n_distinct(py_bloc$bloc_fe), "\n\n")

base_rhs <- "z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin + z_trade_mean + z_dep_asym_max + z_cap_asym_max"

m_dis_nofe  <- glm(as.formula(paste("event_dissolve ~", base_rhs)), data = py_bloc, family = binomial())
m_flip_nofe <- glm(as.formula(paste("event_flip ~",     base_rhs)), data = py_bloc, family = binomial())
m_dis_fe    <- glm(as.formula(paste("event_dissolve ~", base_rhs, "+ bloc_fe")), data = py_bloc, family = binomial())
m_flip_fe   <- glm(as.formula(paste("event_flip ~",     base_rhs, "+ bloc_fe")), data = py_bloc, family = binomial())

ct_dis_nofe  <- coeftest(m_dis_nofe,  vcov = vcovCL(m_dis_nofe,  cluster = py_bloc$triad_id))
ct_flip_nofe <- coeftest(m_flip_nofe, vcov = vcovCL(m_flip_nofe, cluster = py_bloc$triad_id))
ct_dis_fe    <- coeftest(m_dis_fe,    vcov = vcovCL(m_dis_fe,    cluster = py_bloc$triad_id))
ct_flip_fe   <- coeftest(m_flip_fe,   vcov = vcovCL(m_flip_fe,   cluster = py_bloc$triad_id))

lrt_dis  <- anova(m_dis_nofe, m_dis_fe, test = "Chisq")
lrt_flip <- anova(m_flip_nofe, m_flip_fe, test = "Chisq")

# two-way (triad + bloc) clustering too, for the same reason it mattered before
cl2 <- py_bloc[, c("triad_id", "primary_bloc")]
ct_dis_2way  <- coeftest(m_dis_nofe,  vcov = vcovCL(m_dis_nofe,  cluster = cl2, multi0 = TRUE))
ct_flip_2way <- coeftest(m_flip_nofe, vcov = vcovCL(m_flip_nofe, cluster = cl2, multi0 = TRUE))

sink("results/cap_asym_bloc_fe_output.txt")
cat("=== DISSOLUTION, no bloc FE ===\n"); print(ct_dis_nofe)
cat("\n=== DISSOLUTION, WITH bloc FE ===\n"); print(ct_dis_fe)
cat("\n=== DISSOLUTION, two-way (triad+bloc) clustered, no FE ===\n"); print(ct_dis_2way)
cat("\n=== LRT: bloc FE jointly significant (dissolution)? ===\n"); print(lrt_dis)

cat("\n\n=== REALIGNMENT, no bloc FE ===\n"); print(ct_flip_nofe)
cat("\n=== REALIGNMENT, WITH bloc FE ===\n"); print(ct_flip_fe)
cat("\n=== REALIGNMENT, two-way (triad+bloc) clustered, no FE ===\n"); print(ct_flip_2way)
cat("\n=== LRT: bloc FE jointly significant (realignment)? ===\n"); print(lrt_flip)

cat("\n\n=== SUMMARY: z_cap_asym_max across specifications ===\n")
summarize_coef <- function(ct, label) {
  cat(sprintf("%-35s b=%8.4f  SE=%7.4f  p=%8.5f\n", label,
              ct["z_cap_asym_max","Estimate"], ct["z_cap_asym_max","Std. Error"], ct["z_cap_asym_max","Pr(>|z|)"]))
}
cat("--- Dissolution ---\n")
summarize_coef(ct_dis_nofe, "1-way clustered, no FE")
summarize_coef(ct_dis_2way, "2-way clustered, no FE")
summarize_coef(ct_dis_fe,   "1-way clustered, WITH bloc FE")
cat("--- Realignment ---\n")
summarize_coef(ct_flip_nofe, "1-way clustered, no FE")
summarize_coef(ct_flip_2way, "2-way clustered, no FE")
summarize_coef(ct_flip_fe,   "1-way clustered, WITH bloc FE")
sink()

cat(readLines("results/cap_asym_bloc_fe_output.txt"), sep = "\n")

saveRDS(list(py_bloc = py_bloc, m_dis_fe = m_dis_fe, m_flip_fe = m_flip_fe,
             ct_dis_nofe = ct_dis_nofe, ct_flip_nofe = ct_flip_nofe,
             ct_dis_fe = ct_dis_fe, ct_flip_fe = ct_flip_fe,
             ct_dis_2way = ct_dis_2way, ct_flip_2way = ct_flip_2way,
             lrt_dis = lrt_dis, lrt_flip = lrt_flip),
        "results/cap_asym_bloc_fe.rds")
cat("\n\n09_cap_asym_bloc_fe.R complete.\n")
