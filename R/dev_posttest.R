## DEV
#' convert to posttest
#' calc power for 2-level model

convert <- function(object, time = NULL) {
    paras <- powerlmm:::NA_to_zero(object)
    if(is.null(time)) time <- paras$T_end

    u0 <- paras$sigma_subject_intercept
    u1 <- paras$sigma_subject_slope
    u01 <- paras$cor_subject * u0 * u1
    v0 <- paras$sigma_cluster_intercept
    v1 <- paras$sigma_cluster_slope
    v01 <- v0 * v1 * paras$cor_cluster
    error <- paras$sigma_error

    # v0 at posttest
    v0 <- (v0^2 + 2*v01*time + v1^2*time^2)

    # error posttest
    error <- (u0^2 + 2*u01*time + u1^2*time^2) + error^2

    list("v0" = v0, "error"=error)
}

get_posttest_power <- function(p) {
    x <- convert(p)

    diff <- powerlmm:::get_slope_diff(p)
    se <- sqrt(2 * (x$error + p$n2 * x$v0) / ( p$n2 * p$n3))

    lambda <- diff/se

    alpha <- 0.05
    df <- get_n3(p)$total - 2
    power <- pt(qt(1-alpha/2, df = df), df = df, ncp = lambda, lower.tail = FALSE) +
        pt(qt(alpha/2, df = df), df = df, ncp = lambda)
    list("power" = power, "se" = se)
}

p <- study_parameters(n1 = 11,
                      n2 = 4,
                      n3 = 10,
                      T_end = 10,
                      icc_pre_subject =  0.8,
                      cor_subject = -0.3,
                      icc_pre_cluster = NA,
                      var_ratio = 0.03,
                      icc_slope = 0.05,
                      cohend = -0.8)

get_posttest_power(p)
get_power(p, df = "satterth")$power


