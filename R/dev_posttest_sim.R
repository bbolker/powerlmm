library(lmerTest)
library(powerlmm)
library(parallel)
library(dplyr)

p <- study_parameters(n1 = 11,
                      n2 = 10,
                      n3 = 5,
                      icc_pre_subject =  0.5,
                      cor_subject = -0.5,
                      cor_cluster = -0.5,
                      icc_pre_cluster = 0.1,
                      var_ratio = 0.02,
                      icc_slope = 0.1,
                      cohend = 0.8)


posttest_sim <- function(i, p) {
    d <- simulate_data(p)
    d_post <- d[d$time == 10, ]
    d_post$pre <- d[d$time == 0, "y"]
    fit <- lmer(y ~ treatment + (1 | cluster), data = d_post)
    x <- summary(fit)
    vv <- as.data.frame(VarCorr(fit))

    ind <- which(rownames(x$coefficients) %in% "treatment")
    x <- as.data.frame(x$coefficients)
    x <- x[ind, c(1,2, 3, 5)]
    x$sigma_cluster2 <- vv$vcov[1]
    x$sigma_error2 <- vv$vcov[2]
    x$icc <- x$sigma_cluster2/(x$sigma_cluster2 + x$sigma_error2)
    x
}

cl <- makeCluster(10)
clusterEvalQ(cl, {
    library(lmerTest)
    library(powerlmm)
    library(parallel)
    library(dplyr)
})
clusterExport(cl, c("posttest_sim", "p"))
res <- parLapply(cl, 1:2000, posttest_sim, p = p)
res <- do.call(rbind, res) %>%
    as.data.frame
stopCluster(cl)

res_sum <- res %>%
    summarise(est = mean(Estimate),
              se = mean(`Std. Error`),
              sd_est = sd(Estimate),
              df = median(df),
              sigma_cluster = mean(sigma_cluster2),
              cluster_prop_zero = mean(near(sigma_cluster2, 0)),
              sigma_error = mean(sigma_error2),
              icc = mean(icc),
              power = mean(`Pr(>|t|)` < 0.05))
res_sum

