
library(lmerTest)
library(dplyr)
library(parallel)
library(powerlmm)

p <- study_parameters(n1 = 11,
                       n2 = 10,
                       n3 = 4,
                       icc_pre_subject = 0.6,
                       icc_pre_cluster = 0.1,
                       cor_cluster = 0,
                       cor_subject = -0.4,
                       var_ratio = 0.03,
                       icc_slope = 0.05,
                       sigma_error = 30,
                       effect_size = 0)

sim_posttest <- function(i, p) {
    d <- simulate_data(p)

    d_post <- d[d$time == 10, ]
    d_post$pre <-  d[d$time == 0, "y"]

    d_post <- d_post %>%
        group_by(cluster) %>%
        mutate(pre_cluster = mean(pre))
    d_post$pre <- d_post$pre - d_post$pre_cluster

    fit_post <- lmer(y ~ treatment  + (1 | cluster), data = d_post)
    x <- summary(fit_post)

    vv <- as.data.frame(VarCorr(fit_post))

    ind <- which(rownames(x$coefficients) == "treatment")

    data.frame(est = x$coefficients[ind, 1],
               se = x$coefficients[ind, 2],
               df = x$coefficients[ind, 3],
               pval = x$coefficients[ind, 5],
               sigma2_v0 = vv$vcov[1],
               sigma2_error = vv$vcov[2])

}

cl <- makeCluster(15)
clusterEvalQ(cl, {
    library(powerlmm)
    library(lmerTest)
    library(dplyr)
})

res <- parLapply(cl, 1:2000, sim_posttest, p = p)

res <- do.call(rbind, res)
stopCluster(cl)

res %>% summarise(mean(est),
                  mean(se),
                  sd(est),
                  median(df),
                  mean(pval < 0.05),
                  mean(sigma2_v0),
                  sigma2_v0_zero = mean(near(sigma2_v0, 0)),
                  mean(sigma2_error))

get_posttest_power(p)
convert(p, adjust_pre = TRUE)
get_autocor(p)

## Work okay when ICC_pre_cluster = 0

