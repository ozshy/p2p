#==========================================================================================
                                    # PREAMBLE #
#==========================================================================================

# Written by: Brian Prescott

# Purpose: This script constructs the summary statistics table reported in the paper
# How People Pay Each Other: Data, Theory, and Calibrations by
# Claire Greene, Brian Prescott, and Oz Shy.

# Notes: In the script we test multiple forms of the mixed logit. However, 
# we only report the one which had the minimum AIC value. Since the primary 
# contribution of the article is not the empirical analysis, we leave all remaining
# specifications in the script for the interested reader. 

#==========================================================================================
                                    # LIBRARIES #
#==========================================================================================

library(tidyverse)
library(xtable)

#==========================================================================================
                                    # DATA IMPORT #
#==========================================================================================

setwd("a:/papers/p2p paper/code/p2p-wp-code/")
cpc_data <- read_rds("p2p_df.rds")

#==========================================================================================
                            # SUMMARY STATISTICS GENERATION #
#==========================================================================================
  
# summary statistics table generation
methods <- c("cash", "elect", "check")

list(mean, sd, median) %>%
  purrr::map(.x = ., 
             .f = function(h) {
               
               cpc_data %>%
                 dplyr::select(-c(id, birthYear, pi, moblie_pay_app, 
                                  in_person, nests, payment_choice, uasid.f, 
                                  contains("year"), device, 
                                  str_c("convenience.", c("Cash", "Check")), 
                                  str_c("security.", c("Cash", "Check")), 
                                  str_c("cost.", c("Cash", "Check")), ends_with("Electronic"))) %>%
                  mutate(household_income = exp(household_income), 
                         transaction_cost = exp(transaction_cost)) %>%
                  summarise_all(.tbl = ., 
                                .funs = ~ h(x = ., na.rm = TRUE)) %>%
                  t() %>%
                  data.frame() %>%
                  rownames_to_column(.data = ., 
                                     var = "Covariate") %>%
                  mutate(Covariate = factor(Covariate, 
                                            levels = c(str_c("cost.", methods), 
                                                       str_c("convenience.", methods), 
                                                       str_c("security.", methods), 
                                                       "transaction_cost", "age", "yrs_of_schooling", 
                                                       "work_employed", "gender", "hispanic", 
                                                       "household_income", "hh_size", "married", 
                                                       "race_black", "race_white", 
                                                       str_c("age_cohort", 
                                                             c("Pre_50s", "50s", "60s", "70s", 
                                                               "80s", "90s_00s"))))) %>%
                  arrange(Covariate) %>%
                  dplyr::rename(h = ".")
               
             }) %>%
  purrr::reduce(.x = ., 
                .f = inner_join, 
                by = "Covariate") %>%
  dplyr::rename(mean = h.x, 
                sd = h.y, 
                median = h) %>%
  xtable::xtable() %>%
  print(include.rownames = FALSE)
            
