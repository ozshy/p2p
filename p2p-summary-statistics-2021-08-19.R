#===============================================================================
                            # PREAMBLE #
#===============================================================================

# Written by: Brian Prescott

# Purpose: This script constructs the summary statistics table reported in the
# paper How People Pay Each Other: Data, Theory, and Calibrations by Claire
# Greene, Brian Prescott, and Oz Shy.

# Notes: In the script we test multiple forms of the mixed logit. However, we
# only report the one which had the minimum AIC value. Since the primary
# contribution of the article is not the empirical analysis, we leave all
# remaining specifications in the script for the interested reader.

#===============================================================================
                            # PACKAGES #
#===============================================================================
library(tidyverse)
library(xtable)

#===============================================================================
                            # DATA IMPORT #
#===============================================================================
setwd(paste0(
  "/home/bp/PrescottLibrary/economics/papers/academic-papers/",
  "p2p/RnR-r1-github-materials/"
))

data_using_date <- "2021-07-08"

cpc_data <- read_rds(file = str_c("p2p_df_", data_using_date, ".rds"))

#===============================================================================
                      # SUMMARY STATISTICS GENERATION #
#===============================================================================
  
# summary statistics table generation
methods <- c("cash", "elect", "check")
list(mean, sd, median) %>%
  purrr::map(
    .x = .,
    .f = function(h) {
      cohorts <- c("Pre_50s", "50s", "60s", "70s", "80s", "90s_00s")
      cpc_data %>%
        dplyr::select(-c(
          id, pi, moblie_pay_app, age_cohort, payment_choice,
          ind_payee, contains("year"),
          str_c("convenience.", c("Cash", "Check")),
          str_c("security.", c("Cash", "Check")),
          str_c("cost.", c("Cash", "Check")),
          ends_with("Electronic")
        )) %>%
         mutate(
           household_income = exp(household_income) - 1,
           transaction_cost = exp(transaction_cost) - 1
          ) %>%
          summarise_all(
            .tbl = .,
            .funs = ~ h(x = ., na.rm = TRUE)
          ) %>%
          t() %>%
          data.frame() %>%
        rownames_to_column(
          .data = .,
          var = "Covariate"
        ) %>%
        mutate(
          Covariate = factor(
            x = Covariate,
            levels = c(
              "age", "gender", "hh_size", "race_white", "race_black",
              "work_employed", "ind_payee",
              str_c("cost.", methods), str_c("convenience.", methods),
              str_c("security.", methods), "transaction_cost",
              "household_income",
              "yrs_of_schooling", "hispanic", "married",
              str_c("context_", c("goodServ", "repay", "lend", "donation")),
              str_c("payee_", c("goodServ.busi", "goodServ.nonbusi",
                                "friFam", "profPeer", "other"))
            )
          )
        ) %>%
        arrange(Covariate) %>%
        dplyr::rename(h = ".")
    }
  ) %>%
  purrr::reduce(
    .x = .,
    .f = inner_join,
    by = "Covariate"
  ) %>%
  dplyr::rename(
    mean = h.x,
    sd = h.y,
    median = h
  ) %>%
  xtable::xtable() %>%
  print(include.rownames = FALSE)