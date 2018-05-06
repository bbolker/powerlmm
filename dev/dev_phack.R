p <- study_parameters(n1 = 11,
                      n2 = 40, # per treatment
                      icc_pre_subject = 0.5,
                      cor_subject = -0.5,
                      var_ratio = 0.02,
                      effect_size = 0)

# Formulas --------------------------------------------------------------------
# OLS (t-test)
f_PT <- sim_formula("y ~ treatment",
                    test = "treatment",
                    data_transform = transform_to_posttest)

# ANCOVA
f_PT_pre <- sim_formula("y ~ treatment + pretest",
                        test = "treatment",
                        data_transform = transform_to_posttest)

# analyze as 2-level longitudinal
f_LMM <- sim_formula("y ~ time * treatment +
                     (1 + time | subject)")

# constrain treatment differences at pretest
f_LMM_c <- sim_formula("y ~ time + time:treatment +
                         (1 + time | subject)")

# combine formulas
f <- sim_formula_compare("posttest" = f_PT,
                         "ANCOVA" = f_PT_pre,
                         "LMM" = f_LMM,
                         "LMM_c" = f_LMM_c)

# Run sim --------------------------------------------------------------------
res <- simulate(p,
                formula = f,
                nsim = 1000,
                cores = 4,
                satterthwaite = TRUE,
                batch_progress = FALSE)

# need to specify what parameter estimates the treatment effect.
tests <-  list("posttest" = "treatment",
               "ANCOVA" = "treatment",
               "LMM" = "time:treatment",
               "LMM_c" = "time:treatment")

summary(res, para = tests)

x <- p_hack.plcp_sim(res, para = tests)
mean(x$pval < 0.05)
mean(x$estimate)
