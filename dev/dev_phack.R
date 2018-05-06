library(powerlmm)

p <- study_parameters(n1 = 11,
                      n2 = 60, # per treatment
                      icc_pre_subject = 0.5,
                      cor_subject = -0.4,
                      var_ratio = c(0.02),
                      effect_size = 0)

transform_to_pre <- function(data) {
    tmp <- data[data$time != min(data$time), ]
    tmp$pretest <- data[data$time == min(data$time), "y"]

    tmp
}


transform_to_posttest_diff <- function(data) {
    tmp <- data[data$time == max(data$time), ]
    tmp$pretest <- data[data$time == min(data$time), "y"]
    tmp$y <- tmp$y - tmp$pretest
    tmp
}




# Formulas --------------------------------------------------------------------
# OLS (t-test)
f_PT <- sim_formula("y ~ treatment",
                    test = "treatment",
                    data_transform = transform_to_posttest)

# diff score
f_diff <- sim_formula("y ~ treatment",
                    test = "treatment",
                    data_transform = transform_to_posttest_diff)


# ANCOVA
f_PT_pre <- sim_formula("y ~ treatment + pretest",
                        test = "treatment",
                        data_transform = transform_to_posttest)

# diff score
f_diff_pre <- sim_formula("y ~ treatment + pretest",
                        test = "treatment",
                        data_transform = transform_to_posttest_diff)

# analyze as 2-level longitudinal
f_LMM <- sim_formula("y ~ time * treatment +
                     (1 + time | subject)")

# constrain treatment differences at pretest
f_LMM_c <- sim_formula("y ~ time + time:treatment +
                         (1 + time | subject)")

f_LMM_pre <- sim_formula("y ~ pretest + time * treatment +
                         (1 + time | subject)",
                         test = c("time:treatment"),
                         data_transform = transform_to_pre)


f_LMM_pre2 <- sim_formula("y ~ pretest + pretest:time + time * treatment +
                         (1 + time | subject)",
                          test = c("time:treatment"),
                          data_transform = transform_to_pre)


# combine formulas
f <- sim_formula_compare("posttest" = f_PT,
                         "ANCOVA" = f_PT_pre,
                         "diff" = f_diff,
                         "diff_pre" = f_diff_pre,
                         "LMM" = f_LMM,
                         "LMM_c" = f_LMM_c,
                         "LMM_pre" = f_LMM_pre,
                         "LMM_pre2" = f_LMM_pre2)



# Run sim --------------------------------------------------------------------
res <- simulate(p,
                formula = f,
                nsim = 1000,
                cores = 16,
                satterthwaite = TRUE,
                batch_progress = FALSE)

# need to specify what parameter estimates the treatment effect.
tests <-  list("posttest" = "treatment",
               "ANCOVA" = "treatment",
               "diff" = "treatment",
               "diff_pre" = "treatment",
               "LMM" = "time:treatment",
               "LMM_c" = "time:treatment",
               "LMM_pre" = "time:treatment",
               "LMM_pre2" = "time:treatment")

summary(res, para = tests)

mean(p_hack.plcp_sim(res, para = tests)$pval < 0.05)


x <- lapply(res, function(x) { mean(p_hack.plcp_sim(x, para = tests)$pval < 0.05)})
x

mean(x$pval < 0.05)
mean(x$estimate)
