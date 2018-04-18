## DEV
#' convert to posttest
#' calc power for 2-level model

convert <- function(object, adjust_pre = TRUE, time = NULL) {
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


    if(adjust_pre) {
        n2 <- paras$n2
        rho_W <- get_autocor(object)[[1]]
        # This is incorrect,
        rho_B <- get_autocor(object)[[2]]
        error <- (1 - n2/(n2-1) * rho_W^2) * error
        v0 <- (1 - rho_B^2 + 1/(n2-1) * rho_W^2) * v0
    }

    list("v0" = v0, "error"=error)
}

get_posttest_power <- function(p, adjust_pre = TRUE) {
    x <- convert(p, adjust_pre = adjust_pre)

    diff <- powerlmm:::get_slope_diff(p)
    se <- sqrt(2 * (x$error + p$n2 * x$v0) / ( p$n2 * p$n3))

    lambda <- diff/se

    alpha <- 0.05
    df <- get_n3(p)$total - 2 - 1
    power <- pt(qt(1-alpha/2, df = df), df = df, ncp = lambda, lower.tail = FALSE) +
        pt(qt(alpha/2, df = df), df = df, ncp = lambda)
    list("power" = power, "se" = se, "df" = df, "v0" = x$v0, "error" = x$error)
}
get_autocor <- function(object) {
    paras <- NA_to_zero(object)
    u0 <- paras$sigma_subject_intercept
    u1 <- paras$sigma_subject_slope
    v0 <- paras$sigma_cluster_intercept
    v1 <- paras$sigma_cluster_slope
    v01 <- v0 * v1 * paras$cor_cluster
    error <- paras$sigma_error
    T_end <- paras$T_end

    u01 <- paras$cor_subject * u0 * u1

    time <- get_time_vector(paras)

    num <- (u0^2 + u01*(time[1] + T_end) + u1^2*time[1]*T_end)
    tot0 <- (u0^2 + 2*u01*time[1] + u1^2*time[1]^2 + error^2)
    tot1 <- (u0^2 + 2*u01*T_end + u1^2*T_end^2 + error^2)

    lvl2 <- num/sqrt(tot0 * tot1)

    num <- (v0^2 + v01*(time[1] + T_end) + v1^2*time[1]*T_end)
    tot0 <- (v0^2 + 2*v01*time[1] + v1^2*time[1]^2)
    tot1 <- (v0^2 + 2*v01*T_end + v1^2*T_end^2)

    lvl3 <- num/sqrt(tot0 * tot1)
    lvl3 <- ifelse(is.nan(lvl3), 0, lvl3)

    list(lvl2, lvl3)

}

