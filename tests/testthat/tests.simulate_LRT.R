## single
test_that("stepwise LRT", {
    models <- list(list("label" = "m0",
                        "ll" = -10,
                        "df" = 1),
                   list("label" = "m1",
                        "ll" = -8,
                        "df" = 2),
                   list("label" = "m2",
                        "ll" = -8,
                        "df" = 3),
                   list("label" = "m3",
                        "ll"= -6,
                        "df" = 4)
    )
    expect_equal(step_bw.plcp_sim(models), "m3") # should pick 'm3'
    expect_equal(step_fw.plcp_sim(models), "m1") # should pick 'm1'

    ## multiple
    models <- list(list("label" = "m0",
                        "ll" = c(-10, -10),
                        "df" = c(1,1)),
                   list("label" = "m1",
                        "ll" = c(-10, -8),
                        "df" = c(2,2)),
                   list("label" = "m2",
                        "ll" = c(-10, -10),
                        "df" = c(3,3)),
                   list("label" = "m3",
                        "ll"= c(-10, -10),
                        "df" = c(4,4))
    )
    expect_equal(step_bw.plcp_sim(models), c("m0", "m1")) # m0, m1
    expect_equal(step_fw.plcp_sim(models), c("m0", "m1")) # m0, m1

    #
    models <- list(list("label" = "m0",
                        "ll" = c(-10, -10),
                        "df" = c(1,1)),
                   list("label" = "m1",
                        "ll" = c(-10, -8),
                        "df" = c(2,2)),
                   list("label" = "m2",
                        "ll" = c(-6, -10),
                        "df" = c(3,3)),
                   list("label" = "m3",
                        "ll"= c(-10, -10),
                        "df" = c(4,4))
    )
    expect_equal(step_bw.plcp_sim(models), c("m2", "m1")) # m2, m1
    expect_equal(step_fw.plcp_sim(models), c("m0", "m1")) # m0, m1

    models <- list(list("label" = "m0",
                        "ll" = c(-10, -10),
                        "df" = c(1, 1)),
                   list("label" = "m1",
                        "ll" = c(-8, -8),
                        "df" = c(2, 2)),
                   list("label" = "m2",
                        "ll" = c(-6, -6),
                        "df" = c(3, 3)),
                   list("label" = "m3",
                        "ll"= c(-4, -5),
                        "df" = c(4, 4))
    )
    expect_equal(step_bw.plcp_sim(models), c("m3", "m2")) # m3, m2
    expect_equal(step_fw.plcp_sim(models), c("m3", "m2")) # m3, m2
})


test_that("sim LRT", {
    ## sim
    p <- study_parameters(n1 = 5,
                          n2 = 5,
                          icc_pre_subject = 0.5,
                          cor_subject = -0.5,
                          var_ratio = 0.03)

    f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
    f1 <- sim_formula("y ~ time * treatment + (1 + time || subject)")
    f2 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
    f <- sim_formula_compare("m0" = f0, "m1" = f1, "m2" = f2)


    res <- simulate(p, formula = f, nsim = 4, satterthwaite = FALSE, cores = 1, CI = FALSE)

    x <-  summary(res)
    expect_equal(names(x$summary), c("m0", "m1", "m2"))

    x <- summary(res, model_selection = "FW")
    expect_equal(names(x$summary), "model_selection")
    expect_equal(x$model_direction, "FW")
})

test_that("sim LRT multi", {
    ## sim
    p <- study_parameters(n1 = 5,
                          n2 = 5:6,
                          icc_pre_subject = 0.5,
                          cor_subject = -0.5,
                          var_ratio = 0.03)

    f0 <- sim_formula("y ~ time * treatment + (1 | subject)")
    f1 <- sim_formula("y ~ time * treatment + (1 + time || subject)")
    f2 <- sim_formula("y ~ time * treatment + (1 + time | subject)")
    f <- sim_formula_compare("m0" = f0, "m1" = f1, "m2" = f2)


    res <- simulate(p, formula = f, nsim = 4, satterthwaite = FALSE, cores = 1, CI = FALSE, batch_progress = FALSE)

    x <-  summary(res)
    expect_output(print(x), "^Model: 'All' | Type: fixed$")

    x <- summary(res, model_selection = "FW")
    expect_output(print(x), "^Model: 'model_selection' | Type: fixed$")
})

test_that("LRT calcs", {
    p <- study_parameters(n1 = 5,
                          n2 = 5,
                          T_end = 10,
                          icc_pre_subject = 0.5,
                          cor_subject = -0.5,
                          var_ratio = 0.05)

    d <- simulate_data(p)
    fit0 <- lme4::lmer(y ~ time * treatment + (1 | subject), data = d)
    fit1 <- lme4::lmer(y ~ time * treatment + (1 + time | subject), data = d)

    av <- anova(fit0, fit1, refit = FALSE)

    ll0 <- stats::logLik(fit0)
    df0 <- attr(ll0, "df")
    ll1 <- stats::logLik(fit1)
    df1 <- attr(ll1, "df")

    m0 <- list("label" = "m0",
               "ll" = as.numeric(ll0),
               "df" = df0)

    m1 <- list("label" = "m1",
               "ll" = as.numeric(ll1),
               "df" = df1)

    expect_equivalent(comp_LRT(m0, m1), av$`Pr(>Chisq)`[2])
})


test_that("p-hack", {

    p <- study_parameters(n1 = 3,
                          n2 = 10:11, # per treatment
                          icc_pre_subject = 0.5,
                          cor_subject = -0.4,
                          var_ratio = 0.02,
                          effect_size = 0)



    # Formulas --------------------------------------------------------------------
    # OLS (t-test)
    f_PT <- sim_formula("y ~ treatment",
                        test = "treatment",
                        data_transform = transform_to_posttest)


    # diff score
    f_PT_pre <- sim_formula("y ~ treatment + pretest",
                              test = "treatment",
                              data_transform = transform_to_posttest)

    # analyze as 2-level longitudinal
    f_LMM <- sim_formula("y ~ time * treatment +
                         (1 + time | subject)")



    # combine formulas
    f <- sim_formula_compare("posttest" = f_PT,
                             "ANCOVA" = f_PT_pre,
                             "LMM" = f_LMM)



    # Run sim --------------------------------------------------------------------
    res <- simulate(p,
                    formula = f,
                    nsim = 2,
                    cores = 1,
                    satterthwaite = FALSE,
                    batch_progress = FALSE)

    # need to specify what parameter estimates the treatment effect.
    tests <-  list("posttest" = "treatment",
                   "ANCOVA" = "treatment",
                   "LMM" = "time:treatment")


    # should throw error
    expect_error(summary(res[[1]], model_selection = "p-hack"), "'para' can't be NULL.")
    expect_error(summary(res, model_selection = "p-hack", para = "treatment1"), "No 'para': treatment1, in model: posttest")
    expect_error(summary(res, model_selection = "p-hack", para = list("posttest"= "treatment")), "When 'para' is a list it must contain a parameter name for each model.")
    expect_error(summary(res, model_selection = "p-hack", para = list("posttest"= "treatment",
                                                                      "ANCOVA"= "treatment",
                                                                      "LMM"= "treatment1")), "No 'para': treatment1, in model: LMM")


    x <- summary(res[[1]], model_selection = "p-hack", para = list("posttest"= "treatment",
                                                         "ANCOVA"= "treatment",
                                                         "LMM"= "time:treatment"))
    expect_output(print(x), "Fixed effects: 'p-hack'")
    expect_output(print(x), "model M_est")
    expect_output(print(x), "Results based on p-hacking")

    # print multi
    x <- summary(res, model_selection = "p-hack", para = list("posttest"= "treatment",
                                                                   "ANCOVA"= "treatment",
                                                                   "LMM"= "time:treatment"))
    expect_equal(nrow(x), 2)
    expect_output(print(x), "Model: 'model_selection' \\| Type: 'fixed' \\| Parameter\\(s\\): 'p-hack'")
    expect_output(print(x), "M_est theta M_se")
    expect_output(print(x), "nsim:  2")

})



