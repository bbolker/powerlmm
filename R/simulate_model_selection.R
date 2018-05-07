# Stepwise ----------------------------------------------------------------
# postprocessing functions
# ll and df are saved each simulation

LRT_test <- function(m0, m1, alpha) {
    pvals <- comp_LRT(m0, m1)
    x <- which(pvals < alpha)
    if(length(x) == 0) x <- NA

    x
}
comp_LRT <- function(m0, m1) {



    if(unique(m0$df) < unique(m1$df)) {
        dev0 <- -2 * m0$ll
        df0 <- m0$df

        dev1 <- -2 * m1$ll
        df1 <- m1$df
    } else {
        dev0 <- -2 * m1$ll
        df0 <- m1$df

        dev1 <- -2 * m0$ll
        df1 <- m0$df
    }


    pval <- 1 - stats::pchisq(dev0-dev1, df1-df0)

    pval
}

prepare_LRT_models <- function(object) {
    models <- lapply(seq_along(object$res), function(i) {
        m <- object$res[i]
        list("label" = names(m),
             "ll" = m[[1]]$logLik,
             "df" = m[[1]]$df)
    })

    models
}

# @return vector of picked models; character vector
step_fw.plcp_sim <- function(models, alpha = 0.1) {

    m0 <- models[[1]]
    m1 <- models[[2]]

    res0 <- rep(m0$label, length(m0$ll))
    winners <- LRT_test(m0, m1, alpha = alpha)
    res0[winners] <- m1$label

    if(length(models) > 2) {
        for(i in seq_along(models)[-(1:2)]) {
            new_winners <- LRT_test(models[[i-1]], models[[i]], alpha = alpha)
            # only want to keep models that were also selected in previous steps
            winners <- new_winners[new_winners %in% winners]

            res0[winners] <- models[[i]]$label
        }
    }
    res0
}
# @return vector of picked models; character vector
step_bw.plcp_sim <- function(models, alpha = 0.1) {

    K <- length(models)
    m0 <- models[[K-1]]
    m1 <- models[[K]]

    res0 <- rep(models[[1]]$label, length(m0$ll))
    # keep these models
    winners <- LRT_test(m0, m1, alpha = alpha)
    res0[winners] <- m1$label

    if(length(models) > 2) {
        for(i in rev(seq_along(models)[-c(K-1, K)])) {
            new_winners <- LRT_test(models[[i]], models[[i + 1]], alpha = alpha)
            # only want to keep testing models that were NOT selected in previous steps
            tmp <- new_winners[!new_winners %in% winners]
            res0[tmp] <- models[[i + 1]]$label

            # only continue elimination models that are non-sig
            winners <- c(tmp, winners)
            #
        }
    }
    res0
}

## p-hacking
p_hack.plcp_sim <- function(object, para) {

    tmp <- vector("list", length(object$res))
    names(tmp) <- names(object$res)
    for(i in seq_along(tmp)) {
        d <- object$res[[i]]$FE
        test <- para[[names(tmp)[i]]]
        tmp[[i]] <- d[d$parameter == test, ]
        tmp[[i]]$model <- names(tmp)[[i]]
    }

    tmp <- do.call(rbind, tmp)
    out <- vector("list", object$nsim)
    for(i in 1:object$nsim) {
        d <- tmp[tmp$sim == i, ]
        out[[i]] <- d[which.min(d$pval), ]

    }
    x <- do.call(rbind, out)
    x$parameter <- "p-hack"
    x

}
# Update object -----------------------------------------------------------


# helper that returns the relevant estimates from each sim
# from the model picked by LRT test
get_sim_para <- function(i, effect, object, mod) {
    x <- object$res[[mod]][[effect]]
    x <- x[x$sim == i, ]
    x$mod <- mod
    x
}

do_model_selection <- function(object, direction = "FW", alpha = 0.1, para = NULL) {
    models <- prepare_LRT_models(object)
    if(direction == "FW") {
        winners <- step_fw.plcp_sim(models, alpha = alpha)
    } else if(direction == "BW") {
        winners <- step_bw.plcp_sim(models, alpha = alpha)
    }

    nsim <- object$nsim

    if(direction %in% c("FW", "BW")) {   # skeleton
        x <- list("RE" = vector("list", nsim),
                  "FE" = vector("list", nsim),
                  "tot_n" = vector("numeric", nsim),
                  "convergence" = vector("numeric", nsim))

        for(i in seq_along(winners)) {
            mod <- winners[i]
            x$RE[[i]] <- get_sim_para(i, "RE", object, mod)
            x$FE[[i]] <- get_sim_para(i, "FE", object, mod)
            x$tot_n[i] <- object$res[[mod]][["tot_n"]][i]
            x$convergence[i] <- object$res[[mod]][["convergence"]][i]
        }
        x$RE <- do.call(rbind, x$RE)
        x$FE <- do.call(rbind, x$FE)
    } else if(direction == "p-hack") {
        x <- p_hack.plcp_sim(object, para = para)
        winners <- x$model
        x <- list("RE" = data.frame(grp = numeric(), var1 = numeric(), var2 = numeric(), vcov = numeric(), parameter = numeric(), sim = numeric()),
                  "FE" = x,
                  "tot_n" = object$res[[1]]$tot_n,
                  "convergence" = NA)

    }


    object$res <- list("model_selection" = x)
    object$model_selected <- table(winners,
                                   dnn = "Model selected (proportion)")/nsim
    object$model_direction <- direction
    object$LRT_alpha <- alpha
    object
}



