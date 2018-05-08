library(powerlmm)

p <- study_parameters(n1 = 11,
                      n2 = 20, # per cluster
                      n3 = 3,
                      icc_pre_subject = 0.5,
                      icc_pre_cluster = 0.01,
                      cor_subject = -0.2,
                      icc_slope = 0.01,
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
                nsim = 3000,
                cores = 16,
                satterthwaite = TRUE,
                batch_progress = FALSE)

# need to specify what parameter estimates the treatment effect.
tests <-  list("posttest" = "treatment",
               "ANCOVA" = "treatment",
               "LMM" = "time:treatment",
               "LMM_c" = "time:treatment")

summary(res, para = tests)


print(summary(res, model_selection = "p-hack", para = tests), verbose = TRUE)

# should throw error
summary(res, model_selection = "p-hack")



###

# combine formulas
f2 <- sim_formula_compare("posttest" = f_PT,
                         "ANCOVA" = f_PT_pre,
                         "diff" = f_diff)



# Run sim --------------------------------------------------------------------
res2 <- simulate(p,
                formula = f2,
                nsim = 1000,
                cores = 16,
                satterthwaite = FALSE,
                batch_progress = FALSE)

# need to specify what parameter estimates the treatment effect.
summary(res, model_selection = "p-hack", para = "treatment")

#should  give error
summary(res, model_selection = "p-hack", para =  list("posttest" = "treatment",
                                                      "ANCOVA" = "treatment",
                                                      "diff" = "treatment"))


## TODO: LMMs are ignored if not Satterth dfs used
