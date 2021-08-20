mlogit_ame <- function(thetas) {

    model <- using_mixed_multi #nolint
    model$coefficients <- thetas

    end <- length(using_mixed_multi$model) #nolint
    init_coef_names <- names(using_mixed_multi$model)[-c(1, end - 1, end)]

    first_var <- which(init_coef_names == "transaction_cost")
    last_var <- which(init_coef_names == "hh_size")
    coef_names <- init_coef_names[first_var:last_var]

    # I am copying the code provided by Julius Vainora on StackOverflow.
    num_cores <- 5
    mixl_ame <- coef_names %>%
      parallel::mclapply(
        X = .,
        FUN = function(variable) {
          out <- variable %>%
            stats::effects(
              object = model,
              covariate = .,
              type = "aa",
              data = mlogit_cpc_v1,
            )
          return(out)
        },
        mc.cores = num_cores,
        mc.cleanup = TRUE,
        mc.preschedule = TRUE
      )
    me_vector <- mixl_ame %>%
      lapply(
        X = .,
        FUN = colMeans
      ) %>%
      unlist() %>%
      as.numeric()

    return(me_vector)
}


# This is the code to compute the numerical gradient from the numDeriv package.
mixl_jacobian <- function(func,
                          x,
                          method, # "simple" or "richardson"
                          method_args = list()
                          ) {

  f <- func(x)
  n <- length(x) # number of variables.

  side <- rep(NA, n)

  if (tolower(method) == "simple") {
    args <- list(eps = 1e-4) # default
    args[names(method_args)] <- method_args

    side[is.na(side)] <- 1
    eps <- rep(args$eps, n) * side

    df <- matrix(NA, nrow = length(f), ncol = n)
    for (i in 1:n) {
      if (i == 1) cat("i iterations left: ", sep = "")
      cat(n - i + 1, " ")

      dx <- x
      dx[i] <- dx[i] + eps[i]
      df[, i] <- (func(dx) - f) / eps[i]

     }

    return(df)

  } else if (tolower(method) == "richardson") {

    args <- list(
      eps = 1e-4,
      d = 0.0001,
      zero.tol = sqrt(.Machine$double.eps / 7e-7),
      r = 4,
      v = 2,
      show.details = FALSE
    ) # default
    args[names(method_args)] <- method_args
    d <- args$d
    r <- args$r
    v <- args$v
    a <- array(NA, c(length(f), r, n))

    h <- abs(d * x) + args$eps * (abs(x) < args$zero.tol)
    pna <- (side == 1)  & !is.na(side) # double these on plus side
    mna <- (side == -1) & !is.na(side) # double these on minus side

    for (k in 1:r) { # successively reduce h
        if (k == 1) cat("h iterations left: ", sep = "")
        cat(r - k, " ", sep = "")
        ph <- mh <- h
        ph[pna] <- 2 * ph[pna]
        ph[mna] <- 0
        mh[mna] <- 2 * mh[mna]
        mh[pna] <- 0

        for (i in 1:n) {
            fd_plus <- func(x + ph * (i == seq(n)))
            fd_minus <- func(x - mh * (i == seq(n)))
            a[, k, i] <- (fd_plus - fd_minus) / (2 * h[i])

        }

        h <- h / v     # Reduced h by 1/v.

    }

   for (m in 1:(r - 1)) {
        if (m == 1) cat("m iterations left: ", sep = "")
        cat((r - 1) - m, " ", sep = "")

       a_1 <- a[, 2:(r + 1 - m), , drop = FALSE]
       a_2 <- a[, 1:(r - m), , drop = FALSE]
       a <- (a_1 * (4^m) - a_2) / (4^m - 1)

    }
    # drop second dim of a, which is now 1
    # (but not other dim's even if they are 1
    return(array(a, dim(a)[c(1, 3)]))

  } else {
    stop("This method is not supported.")
  }

}