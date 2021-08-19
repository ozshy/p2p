#===============================================================================
                                  # PREAMBLE #
#===============================================================================

# Written by: Brian Prescott

# Purpose: This script estimates the mixed logit reported in the paper
# How People Pay Each Other: Data, Theory, and Calibrations.

# Notes: In the script we test multiple forms of the mixed logit. However, we
# only report the one which had the minimum AIC value. Since the primary
# contribution of the article is not the empirical analysis, we leave all
# remaining specifications in the script for the interested reader. The code
# is styled so that a line does not exceed 80 characters and follows the 
# tidyverse style guide elsewhere.

# WARNING: The code to compute the numerical gradient takes approximately one
# day to run if estimated without using parallel computing. Else, it will take
# a few hours. If replicating from scratch, then take this into consideration.
# The same applies to running the imputations.

#===============================================================================
                                  # FILE PATHS #
#===============================================================================

# Loading in my packages
library(plyr) # Data manipulation
library(tidyverse) # Data manipulation
library(writexl) # exporting excel data
library(mlogit) # For mixed multinomial logits
library(xtable) # for latex tables
library(numDeriv) # for computation of the Jacobian of marginal effects
library(pracma) # alternative to numDeriv for computing the Jacobian 
library(doParallel) # for the parJacobian() function I will define
library(parallel) # for the detectCores() function I will use
library(mice) # for the multiple imputations
library(VIM) # for the imputation robustness checks

#===============================================================================
                            # PARAMETERS #
#===============================================================================
# INPUT YOUR WD HERE
setwd(paste0(
  "/home/bp/PrescottLibrary/economics/papers/",
  "academic-papers/p2p/RnR-r1-github-materials/"
))

# These variables will be needed if you choose not to run the script from 
# scratch.
data_using_date <- "2021-07-08"
models_using_date <- "2021-08-03"
grad_using_date <- "2021-08-05"
imputations_using_date <- "2021-08-08"
me_imputations_date <- "2021-08-12"

est_models <- FALSE
run_other_specifications <- FALSE
est_with_imputations <- FALSE

#===============================================================================
                            # FUNCTIONS #
#===============================================================================
mean_omit_NA <- function(y) {
  return(mean(x = y, na.rm = TRUE))
}

#===============================================================================
                        # DATA IMPORT AND EDITING #
#===============================================================================

# Final year of the DCPC data. We will not use 2020 data as we are interested
# in establishing a pre-COVID baseline. 
year <- 2019

# This is the vintage of the data we will be using.
init_cpc_data <- read_rds(
  file = str_c("p2p_df_", data_using_date, ".rds")
 ) %>%
  mutate(
    context = case_when(
      context_repay == 1 ~ "repay", 
      context_lend == 1 ~ "lend", 
      context_donation == 1 ~ "donate", 
      context_goodServ == 1 ~ "goods or services", 
      context_repay == 0 & context_lend == 0 & 
        context_donation == 0 & context_goodServ == 0 ~ "other"
    ),
    payee = case_when(
      payee_goodServ.nonbusi == 1 ~ "goods/services (non-business)",
      payee_goodServ.busi == 1 ~ "goods/services (business)",
      payee_friFam == 1 ~ "Friends or family",
      payee_profPeer == 1 ~ "Coworker, classmate, or fellow military",
      payee_other == 1 ~ "Other",
    ),
    context = factor(context), 
    payee = factor(payee)
  ) %>%
  dplyr::select(
    starts_with("cost."), starts_with("convenience."), starts_with("security."),
    payment_choice, transaction_cost, age, gender, work_employed, race_white, 
    race_black, hispanic, year2019, year2018, year2017, year2016, married, 
    yrs_of_schooling, household_income, hh_size, id, starts_with("context"), 
    starts_with("payee")
  ) %>%
  dplyr::select(-c(
    str_c("cost.", c("cash", "check", "elect")), 
    str_c("convenience.", c("cash", "check", "elect")), 
    str_c("security.", c("cash", "check", "elect")), 
  ))
  
cpc_data <- na.omit(init_cpc_data) # Imposing the missing at random assumption

#===============================================================================
                # TRANSFORMING DATA FOR LOGIT ESTIMATION #
#===============================================================================

data_cols <- colnames(cpc_data)
col_1 <- which(data_cols == "cost.Cash")
col_2 <- which(data_cols == "security.Electronic")
alt_specific_vars <- c(col_1:col_2)

# Setting up the data in the format required by mlogit()
mlogit_cpc_v1 <- cpc_data %>%
  mutate(choiceid = 1:nrow(cpc_data)) %>%
  arrange(id, choiceid) %>%
  data.frame() %>%
  dfidx::dfidx(
    data = .,
    idx = list(c("choiceid", "id")),
    choice = "payment_choice",
    varying = alt_specific_vars,
    opposite = "cost"
  )

# These are additional variable to test if it is desired. 
mlogit_cpc_v1$goodServ_convenience <- mlogit_cpc_v1$context_goodServ * 
  mlogit_cpc_v1$convenience
mlogit_cpc_v1$goodServ_security <- mlogit_cpc_v1$context_goodServ * 
  mlogit_cpc_v1$security
mlogit_cpc_v1$goodServ_cost <- mlogit_cpc_v1$context_goodServ * 
  mlogit_cpc_v1$cost

mlogit_cpc_v1$repay_convenience <- mlogit_cpc_v1$context_repay * 
  mlogit_cpc_v1$convenience
mlogit_cpc_v1$repay_security <- mlogit_cpc_v1$context_repay * 
  mlogit_cpc_v1$security
mlogit_cpc_v1$repay_cost <- mlogit_cpc_v1$context_repay * 
  mlogit_cpc_v1$cost

mlogit_cpc_v1$lend_convenience <- mlogit_cpc_v1$context_lend * 
  mlogit_cpc_v1$convenience
mlogit_cpc_v1$lend_security <- mlogit_cpc_v1$context_lend * 
  mlogit_cpc_v1$security
mlogit_cpc_v1$lend_cost <- mlogit_cpc_v1$context_lend * 
  mlogit_cpc_v1$cost

mlogit_cpc_v1$donate_convenience <- mlogit_cpc_v1$context_donation * 
  mlogit_cpc_v1$convenience
mlogit_cpc_v1$donate_security <- mlogit_cpc_v1$context_donation * 
  mlogit_cpc_v1$security
mlogit_cpc_v1$donate_cost <- mlogit_cpc_v1$context_donation * 
  mlogit_cpc_v1$cost

#===============================================================================
                      # MIXED DISCRETE CHOICE ESTIMATIONS #
#===============================================================================
set.seed(1995) # setting the seed for the random draws of the mlogit() function

# RnR Note: These specifications were changed to account for the payee side of a
# p2p transaction. We do this by including indicator variables which signal what
# the p2p transaction was for.


if (est_models == TRUE) {
    primary_logit_model <- formula(
    "payment_choice ~ convenience + security + cost |
        transaction_cost + context_goodServ + context_repay + context_lend +
        context_donation +
        payee_goodServ.busi + payee_goodServ.nonbusi + payee_friFam +
        payee_profPeer + age +
        gender + work_employed + race_white +
        race_black + hispanic + year2019 + year2018 + year2017 +
        year2016 + married + yrs_of_schooling + household_income + hh_size | 1"
    )
    # Estimation (1) [THIS IS THE SPECIFICATION WE USE]
    primary_model <- mlogit::mlogit(
    formula = primary_logit_model,
    rpar = c(
        "(Intercept):Check" = "n",
        "(Intercept):Electronic" = "n"
    ),
    reflevel = "Cash",
    data = mlogit_cpc_v1,
    halton = NA, # using halton sequences
    panel = TRUE, # accounting for individual heterogeneity
    R = 2000
    )
    summary(primary_model)
    AIC(primary_model)

    # ALTERNATIVE WAY TO MODEL THE P2P CONTEXT. WE CHOSE VERSION (1) 
    # ABOVE SINCE THE MODEL HAS A BETTER FIT.
    secondary_logit_model <- formula(
    "payment_choice ~ convenience + security + cost +
        goodServ_convenience + goodServ_security + goodServ_cost + 
        repay_convenience + repay_security + repay_cost + 
        lend_convenience + lend_security + lend_cost + 
        donate_convenience + donate_security + donate_cost |
        transaction_cost + age + gender + work_employed + race_white + 
        race_black + hispanic + year2019 + year2018 + year2017 +
        year2016 + married + yrs_of_schooling + household_income + hh_size + 
        context_goodServ + context_repay + context_lend + context_donation +
        payee_goodServ.busi + payee_goodServ.nonbusi + payee_friFam + 
        payee_profPeer | 1"
    )
    secondary_model <- mlogit::mlogit(
    formula = secondary_logit_model,
    rpar = c(
        "(Intercept):Check" = "n",
        "(Intercept):Electronic" = "n"
    ),
    reflevel = "Cash",
    data = mlogit_cpc_v1,
    halton = NA, # using halton sequences
    panel = TRUE, # accounting for individual heterogeneity 
    R = 2000
    ) 
    summary(secondary_model)
    AIC(secondary_model)

    # Saving the R&R models in a separate list
    model_list_RnR <- list(
    primary_model,
    secondary_model
    )
    names(model_list_RnR) <- c("primary_model", "secondary_model")
    save(
    list = c("model_list_RnR"),
    file = str_c("output/mixed_models_RnR_", Sys.Date(), ".Rdata")
    )

    # Only run this if you want to see the other specifications.
    # WARNING: The models are computationally intensive and will 
    # take time to run.
    if (run_other_specifications == TRUE) {
    # Estimation (2): This is the primary model but without the context 
    # specific information.
    mixed_multi <- mlogit(
        payment_choice ~ convenience + security + cost |
        transaction_cost + age + gender + work_employed + 
        race_white + race_black + hispanic + year2019 + 
        year2018 + year2017 + year2016 + married + 
        yrs_of_schooling + household_income + hh_size + 
        context_goodServ + context_repay + context_lend + context_donation | 
        1,
        rpar = c(
        "convenience" = "n",
        "security" = "n",
        "cost" = "n"
        ),
        reflevel = "Cash",
        data = mlogit_cpc_v1,
        halton = NA, # using halton sequences
        panel = TRUE, # accounting for individual heterogeneity
        R = 2000
    ) 
    print(summary(mixed_multi)) 
    print(AIC(mixed_multi))
    
    # Estimation (3)
    mixed_multi.wCorr <- update(
        mixed_multi, 
        correlation = TRUE
    )
    print(summary(mixed_multi.wCorr))
    print(AIC(mixed_multi.wCorr))
    
    # Estimation (4): This is the primary model which allows for correlation
    # across the random intercepts. 
    primary_model.wCorr <- update(
        primary_model, 
        correlation = TRUE
    )
    print(summary(primary_model.wCorr))
    print(AIC(primary_model.wCorr))
    
    # The predicted shares
    predShares <- fitted(
        mixed_multi, 
        outcome = FALSE
    ) %>%
        apply(
        X = ., 
        MARGIN = 2,
        FUN = function(y) {
            return(mean(x = y, na.rm = TRUE))
        }
        )
    print(predShares)

    # Selecting the model which minimizes the AIC
    model_list <- list(
      primary_model,  # The revised models
      secondary_model, # The original models
      mixed_multi,
      mixed_multi.wCorr,
      primary_model.wCorr,
      primary_model.wCorr.wCorr
    )

    save(
        list = c("model_list"),
        file = str_c("mixed_models_", Sys.Date(), ".Rdata")
    )

    # Choosing the mixed logit which minimizes the AIC
    model_choice <- which.min(unlist(purrr::map(.x = model_list, .f = AIC)))
    using_mixed_multi <- model_list[[model_choice]]
    summary(using_mixed_multi)
    coef(using_mixed_multi)

    }

} else {
  # These will be needed later on if you choose to run the imputations. 
  # Since they take up minimal RAM, I do not worry about brining them in.
    primary_logit_model <- formula(
    "payment_choice ~ convenience + security + cost |
        transaction_cost + context_goodServ + context_repay +
        context_lend + context_donation +
        payee_goodServ.busi + payee_goodServ.nonbusi +
        payee_friFam + payee_profPeer +
        age + gender + work_employed + race_white +
        race_black + hispanic + year2019 + year2018 + year2017 +
        year2016 + married + yrs_of_schooling + household_income + hh_size | 1"
    )
    # Loading the R&R models only back in. If you want the other models you
    # will have to run them from scratch.
    load(
      file = stringr::str_c("mixed_models_RnR_", models_using_date, ".Rdata")
    )

}

#===============================================================================
                    # AVERAGE MARGINAL EFFECT ESTIMATIONS #
#===============================================================================

# Only run this if you have already run the models above. 
using_mixed_multi <- model_list_RnR$primary_model
print(summary(using_mixed_multi))

# Coefficient names plus manually creating the coefficient names
end <- length(using_mixed_multi$model)
init_coef_names <- names(using_mixed_multi$model)[-c(1, end - 1, end)]

first_var <- which(init_coef_names == "transaction_cost")
last_var <- which(init_coef_names == "hh_size")
coef_names <- init_coef_names[first_var:last_var]

# The number of cores I use on my machine to run the computation.
# Note, this only works if you are running a Linux or macOS machine.
# If not, then you will need to re-format the below mclapply() to
# a standard lapply() function. I also had to set this manually b/c
# I was running into RAM problems. You may have to run the code a few
# times in case your cores do not communicate properly.
num_cores <- 5
{ #nolint
  start_time <- Sys.time()
  p2p_ame <- as.list(coef_names) %>%
    parallel::mclapply(
      X = .,
      FUN = function(x) {
        out <- x %>%
          stats::effects(
            object = using_mixed_multi,
            covariate = .,
            type = "aa",
            data = mlogit_cpc_v1,
          ) %>%
          apply(
            X = .,
            MARGIN = 2,
            FUN = mean_omit_NA
          )

        return(out)
      },
      mc.preschedule = TRUE,
      mc.cores = num_cores,
      mc.cleanup = TRUE
    ) %>%
    unlist() %>%
    as.numeric()
    print(Sys.time() - start_time)
}
print(p2p_ame)

# Putting the point estimates for the AME into a data.frame for the SE
# estimation
me_matrix <- p2p_ame %>%
  matrix(
    data = .,
    ncol = 3,
    byrow = TRUE
  ) %>%
  data.frame() %>%
  as_tibble() %>%
  dplyr::rename(
    Cash = X1,
    Check = X2,
    Electronic = X3
  ) %>%
  mutate(
    Covariate = coef_names
  ) %>%
  dplyr::select(Covariate, everything())

# Printing the point estimates
me_matrix %>%
  mutate_at(
    .vars = vars(-Covariate),
    .funs = ~ round(x = ., digits = 4)
  ) %>%
  data.frame() %>%
  print()

if (est_models == TRUE) {
    #==========================================================================
                    # MARGINAL EFFECT STANDARD ERROR ESTIMATIONS #
    #==========================================================================
    # This section will compute the SE and CI for the marginal effects of the
    # demographic and transaction variables in the discrete choice model.
    # This section takes about 3 hours to run in its entirety. 

    mlogit_ame <- function(thetas) {

        model <- using_mixed_multi #nolint 
        model$coefficients <- thetas

        end <- length(using_mixed_multi$model) #nolint
        init_coef_names <- names(using_mixed_multi$model)[-c(1, end - 1, end)] 

        first_var <- which(init_coef_names == "transaction_cost")
        last_var <- which(init_coef_names == "hh_size")
        coef_names <- init_coef_names[first_var:last_var]

        # I am copying the code provided by Julius Vainora on StackOverflow.
        num_cores <- 3
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
          lapply(X = ., FUN = colMeans) %>%
          unlist() %>%
          as.numeric()

        return(me_vector)
    }


    # This is the code to compute the  numerical gradient from the
    # numDeriv package. We did the calculation with the "simple" method.
    mixl_jacobian <- function(func,
                              x,
                              method, # "simple" or "richardson"
                              method_args = list()
                              ) {

      f <- func(x)
      n <- length(x)	 # number of variables.

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

    # This function is just an edited version of the numDeriv::jacobian()
    # It allows for parallel computation of the numerical derivative which
    # results in significant speed gains.
    {
    grad_start_time <- Sys.time()
    grad <- mixl_jacobian(
        func = mlogit_ame,
        x = using_mixed_multi$coefficients,
        method = "simple",
        method_args = list()
    )
    print(Sys.time() - grad_start_time)
    }

    save(
      list = c("grad"),
      file = str_c("ame_numGradient_", Sys.Date(), ".Rdata")
    )

} else {
  load(file = stringr::str_c("ame_numGradient_", grad_using_date, ".Rdata"))

}

#===============================================================================

# Matrix of standard errors which are computed using the delta method
multinomial_me_se <- (grad %*% vcov(using_mixed_multi) %*% t(grad)) %>%
  diag() %>%
  sqrt() %>%
  matrix(
    data = .,
    ncol = 3,
    byrow = TRUE
  )

me_matrix_se <- data.frame(
  Covariate = coef_names,
  Cash.SE = multinomial_me_se[, 1],
  Check.SE = multinomial_me_se[, 2],
  Electronic.SE = multinomial_me_se[, 3],
  stringsAsFactors = FALSE
)

marginal_effects_data <- inner_join(
  x = me_matrix,
  y = me_matrix_se,
  by = "Covariate"
)

# Creating a new marginal effects matrix for each banking status variable 
for (b in c("Cash", "Check", "Electronic")) {

  temp_data <- marginal_effects_data %>%
    dplyr::select(Covariate, starts_with(b)) %>%
    data.frame()

  temp_data <- temp_data %>%
    mutate(
      lower.99ci = temp_data[, 2] - (2.576 * temp_data[, 3]),
      upper.99ci = temp_data[, 2] + (2.576 * temp_data[, 3]),
      lower.95ci = temp_data[, 2] - (1.96 * temp_data[, 3]),
      upper.95ci = temp_data[, 2] + (1.96 * temp_data[, 3]),
      lower.90ci = temp_data[, 2] - (1.64 * temp_data[, 3]),
      upper.90ci = temp_data[, 2] + (1.64 * temp_data[, 3])
    )

  assign(str_c(b, "_me_data"), temp_data)

}


sig_gen <- function(data) {

  Variables <- data.frame(
    Covariate = data$Covariate,
    stringsAsFactors = FALSE
  )

  new_data <- data %>%
      dplyr::select(-Covariate) %>%
      round(digits = 4) %>%
      mutate(
        significance = ifelse(
          test = (lower.90ci < 0 & upper.90ci < 0) |
          (lower.90ci > 0 & upper.90ci > 0),
          yes = "*",
          no = ""
        ),
        significance = ifelse(
          test = (lower.95ci < 0 & upper.95ci < 0) |
            (lower.95ci > 0 & upper.95ci > 0),
          yes = "**",
          no = significance
        ),
        significance = ifelse(
          test = (lower.99ci < 0 & upper.99ci < 0) |
            (lower.99ci > 0 & upper.99ci > 0),
          yes = "***",
          no = significance
        )
      )

  output <- cbind(Variables, new_data)

  return(output)
}

marginal_effects <- list(
  Cash_me_data,
  Check_me_data,
  Electronic_me_data
) %>%
  lapply(
    X = .,
    FUN = sig_gen
  )
me_frame <- marginal_effects %>%
  bind_cols() %>%
  dplyr::select(-c(
    starts_with("lower"), starts_with("upper")
  )) %>%
  mutate(
    Covariate...1 = factor(
      Covariate...1,
      levels = c(
        "transaction_cost",
        str_c("context_", c("goodServ", "repay", "lend", "donation")),
        str_c("payee_", c("goodServ.busi", "goodServ.nonbusi",
                          "friFam", "profPeer")),
        "year2016", "year2017",
        "year2018", "year2019", "age",
        "yrs_of_schooling", "gender", "race_black",
        "race_white", "hispanic", "household_income",
        "hh_size", "married", "work_employed"
        )
      )
    ) %>%
  arrange(Covariate...1)
print(me_frame)
View(me_frame)

#===============================================================================
                # INVESTIGATING THE MISSING AT RANDOM ASSUMPTION #
#===============================================================================
# At 100 imputations, the below code took about 5 hours to run in its entirety.
# If you choose to replicate this section, take this into consideration.

if (est_models == TRUE) {
  number_imputations <- 100
  df_imputations <- mice(
    data = init_cpc_data %>%
      dplyr::select(-c(
        starts_with("context_"), starts_with("payee_"),
      )),
    method = c(rep("pmm", 26), rep("polyreg", 2)),
    seed = 1995,
    m = number_imputations
  )

  # This does not need to be computed in parallel
  imputed_dfs <- 1:number_imputations %>%
    as.list() %>%
    purrr::map(
      .x = .,
      .f = function(k) {
        output <- mice::complete(df_imputations, k) %>%
          mutate(
            context_goodServ = as.numeric(context == "goods or services"),
            context_lend = as.numeric(context == "lend"),
            context_repay = as.numeric(context == "repay"),
            context_donation = as.numeric(context == "donate"),
            payee_goodServ.busi = as.numeric(
              payee == "goods/services (business)"
            ),
            payee_goodServ.nonbusi = as.numeric(
              payee == "goods/services (non-business)"
            ),
            payee_other = as.numeric(payee == "Other"),
            payee_profPeer = as.numeric(
              payee == "Coworker, classmate, or fellow military"
            ),
            payee_friFam = as.numeric(payee == "Friends or family")
          ) %>%
          dplyr::select(-c(context, payee))

          return(output)
      }
    )

  # This code takes a very long time to run. You can turn this into a parallel
  # computation to try and speed things up, however, it was killing my RAM so
  # I abandoned such an approach on my own machine.
  imputed_mlogits <- seq_along(imputed_dfs) %>%
    lapply(
      X = .,
      FUN = function(i) {
        cat(i, " ")

        df <- imputed_dfs[[i]]
        data_cols <- colnames(df)
        imp_col_1 <- which(data_cols == "cost.Cash")
        imp_col_2 <- which(data_cols == "security.Electronic")
        alt_specific_vars <- c(imp_col_1:imp_col_2)

        imputed_mlogit_df <- df %>%
          mutate(choiceid = 1:nrow(df)) %>%
          arrange(id, choiceid) %>%
          data.frame() %>%
          dfidx::dfidx(
            data = .,
            idx = list(c("choiceid", "id")),
            choice = "payment_choice",
            varying = alt_specific_vars,
            opposite = "cost"
          )

        model_est <- mlogit::mlogit(
          formula = primary_logit_model,
          rpar = c(
            "(Intercept):Check" = "n",
            "(Intercept):Electronic" = "n"
          ),
          reflevel = "Cash",
          data = imputed_mlogit_df,
          halton = NA, # using halton sequences
          panel = TRUE, # accounting for individual heterogeneity
          R = 2000
        )

        return(model_est)
      }
    )
  save(
    list = c("imputed_dfs", "imputed_mlogits"),
    file = str_c("logit_imputations_", Sys.Date(), ".RData")
  )

  imputed_marginalEffects <- seq_along(imputed_mlogits) %>%
    purrr::map(
      .x = .,
      .f = function(i) {
        cat(i, " ")

        using_model <- imputed_mlogits[[i]]
        df <- imputed_dfs[[i]] 
        data_cols <- colnames(df)
        alt_specific_vars <- c(
          which(data_cols == "cost.Cash"):
          which(data_cols == "security.Electronic")
        )

        using_data <- df %>%
          mutate(choiceid = 1:nrow(df)) %>%
          arrange(id, choiceid) %>%
          data.frame() %>%
          dfidx::dfidx(
            data = .,
            idx = list(c("choiceid", "id")),
            choice = "payment_choice",
            varying = alt_specific_vars,
            opposite = "cost"
          )

        # Coefficient names plus manually creating the coefficient names
        end <- length(using_model$model)
        init_coef_names <- names(using_model$model)[-c(1, end - 1, end)]

        first_var <- which(init_coef_names == "transaction_cost")
        last_var <- which(init_coef_names == "hh_size")
        coef_names <- init_coef_names[first_var:last_var]

        # A quick way to compute the marginal effects for each variable.
        num_cores <- 3
        me <- parallel::mclapply(
          X = as.list(coef_names),
          FUN = function(x) {
            me_calc <- stats::effects(
                object = using_model,
                covariate = x,
                type = "aa",
                data = using_data
              ) %>%
              data.frame() %>%
              apply(
                X = .,
                MARGIN = 2,
                FUN = mean_omit_NA
              )

          return(me_calc)
          },
          mc.cores = 2, # a constraint imposed due to my machine
          mc.cleanup = TRUE,
          mc.preschedule = TRUE
        ) %>%
          bind_rows()

        return(me)
      }
    )
  imputed_marginalEffects <- imputed_marginalEffects %>%
    purrr::map(
      .x = ., 
      .f = function(df) {
        output <- df %>%
        data.frame()
        rownames(output) <- coef_names

      return(output)
      }
    )
  save(
    list = c("imputed_marginalEffects"),
    file = stringr::str_c("marginalEffects_imputations_", Sys.Date(), ".RData")
  )

} else {
  load(
    file = stringr::str_c(
      "logit_imputations_", imputations_using_date, ".RData"
    )
  )

  load(
    file = stringr::str_c(
      "marginalEffects_imputations_", me_imputations_date, ".RData"
    )
  )

}

using_p2p_me <- 1:3 %>%
  as.list() %>%
  purrr::map(
    .x = ., 
    .f = function(i) {
      df <- marginal_effects[[i]]
      method <- c("Cash", "Check", "Electronic")[i]

      output_df <- df %>%
        dplyr::select(Covariate, all_of(method), contains("95")) %>%
        dplyr::rename(covariate = Covariate)

      return(output_df)
    }
  ) %>%
  bind_rows()

me_figureDF <- imputed_marginalEffects %>%
  purrr::map(
    .x = ., 
    .f = function(df) {
      output <- df %>%
        rownames_to_column(
          .data = ., 
          var = "covariate"
        ) %>%
        mutate(which = "imputed")
    }
  ) %>%
  bind_rows() %>%
  bind_rows(
    x = .,
    y = using_p2p_me %>%
      mutate(which = "not imputed")
  )

me_mar_figure_df <- me_figureDF %>%
  dplyr::select(covariate, Cash, which, contains("95")) %>%
  dplyr::filter(!is.na(Cash)) %>%
  dplyr::rename(me = Cash) %>%
  mutate(method = "Cash") %>%
  bind_rows(
    x = ., 
    y = me_figureDF %>%
      dplyr::select(covariate, Check, which, contains("95")) %>%
      dplyr::filter(!is.na(Check)) %>%
      dplyr::rename(me = Check) %>%
      mutate(method = "Check")
  ) %>%
  bind_rows(
    x = .,
    y = me_figureDF %>%
      dplyr::select(covariate, Electronic, which, contains("95")) %>%
      dplyr::filter(!is.na(Electronic)) %>%
      dplyr::rename(me = Electronic) %>%
      mutate(method = "Electronic")
  ) %>%
  mutate(
    covariate = factor(
      x = covariate,
      levels = rev(c(
        "transaction_cost",
        str_c("context_", c("goodServ", "repay", "lend", "donation")),
        str_c("payee_", c("goodServ.busi", "goodServ.nonbusi",
                          "friFam", "profPeer")),
        "year2016", "year2017",
        "year2018", "year2019", "age",
        "yrs_of_schooling", "household_income",
        "work_employed", "hh_size", "married", "gender", "hispanic",
        "race_black", "race_white"
        )),
      labels = rev(c(
        "ln(Transaction value)", "Context:Goods or services",
        "Context:Donation", "Context:Lending", "Context:Repayment",
        "Payee:Goods or services (business)", 
        "Payee:Good or services (non-business)", 
        "Payee:Friends or family", "Payee:Coworker or classmate",
        "Year:2019", "Year:2018", "Year:2017", "Year:2016", "Age",
        "Education", "ln(Household income)", "Employed", "Household size",
        "Married", "Gender", "Latinx", "Race:Black", "Race:White"
      ))
    )
  )

me_mar_figure_df %>%
  ggplot(
    data = ., 
    aes(
      x = covariate, 
      color = which
    )
  ) + 
  geom_point(
    aes(
      y = me, 
      color = "color1", 
      shape = "shape1"
    ),
    size = 4,
    position = position_jitter(width = 0.10), 
    data = me_mar_figure_df %>%
      dplyr::filter(which == "imputed")
  ) + 
  geom_linerange(
    data = me_mar_figure_df %>%
      dplyr::filter(which == "not imputed"),
    aes(
      ymin = lower.95ci, 
      ymax = upper.95ci
    ),
    color = "midnightblue",
    show.legend = FALSE
  ) +
  geom_point(
    aes(
      y = me, 
      color = "color2", 
      shape = "shape2",
    ),
    size = 4,
    data = me_mar_figure_df %>%
      dplyr::filter(which == "not imputed")
  ) + 
  scale_y_continuous(name = "Average marginal effect") + 
  scale_x_discrete(name = NULL) + 
  scale_shape_manual(
    name = NULL, 
    values = c(
      "shape1" = 1, 
      "shape2" = 16
    ),
    labels = c(
      "shape1" = "Imputed estimates", 
      "shape2" = "Non-imputed estimate"
    )
  ) + 
  scale_color_manual(
    name = NULL, 
    values = c(
      "color1" = "skyblue1", 
      "color2" = "midnightblue"
    ),
    labels = c(
      "color1" = "Imputed estimates", 
      "color2" = "Non-imputed estimate"
    )
  ) +
  coord_flip() +
  facet_wrap(~ method) +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.background = element_rect(color = "black"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )