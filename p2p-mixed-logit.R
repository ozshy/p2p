#==========================================================================================
                                    # PREAMBLE #
#==========================================================================================

# Written by: Brian Prescott

# Purpose: This script estimates the mixed logit reported in the paper
# How People Pay Each Other: Data, Theory, and Calibrations.

# Notes: In the script we test multiple forms of the mixed logit. However, 
# we only report the one which had the minimum AIC value. Since the primary 
# contribution of the article is not the empirical analysis, we leave all remaining
# specifications in the script for the interested reader. 

# WARNING: The code to compute the numerical gradient takes approximately one day to run. 
# If replicating take this into consideration.

#==========================================================================================
                                    # FILE PATHS #
#==========================================================================================

# Setting my libraries
library(plyr) # Data manipulation
library(tidyverse) # Data manipulation
library(haven) # Importing Stata data
library(questionr) # Tabulations
library(writexl) # exporting excel data
library(mlogit) # For mixed multinomial logits
library(xtable) # for latex tables
library(numDeriv) # for computation of the gradient of marginal effects

#==========================================================================================
                                    # PARAMETERS #
#==========================================================================================

setwd("c:/users/brian/dropbox/research/p2p-wp-materials/") # INPUT YOUR WD HERE

run_other_specifications <- FALSE

# Setting the seed for any random numbers which need to be drawn. 
# I chose 10 since it was my favorite number growing up.
set.seed(10) 

#==========================================================================================
                                      # FUNCTIONS #
#==========================================================================================

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# A convenience function
not.nan <- function(x) { 
  
  return(!is.nan(x)) 
  
}

#==========================================================================================
                                # DATA IMPORT AND EDITING #
#==========================================================================================

year <- 2019
cpc_data <- read_rds("p2p_df.rds")

alt_specific_vars <- c(which(colnames(cpc_data) == "cost.Cash"): 
                       which(colnames(cpc_data) == "security.Electronic"))

mlogit_cpc_v1 <- cpc_data %>%
  na.omit() %>%
  mutate(choiceid = 1:nrow(cpc_data)) %>%
  arrange(id, choiceid) %>%
  data.frame() %>%
  dfidx::dfidx(data = .,
               idx = list(c("choiceid", "id")),
               choice = "payment_choice",
               varying = alt_specific_vars, 
               opposite = "cost")

#==========================================================================================
                            # MIXED DISCRETE CHOICE ESTIMATIONS #
#==========================================================================================

# Estimation (3) [THIS IS THE SPECIFICATION WE USE]
mixed_multi_noRC <- mlogit(payment_choice ~ convenience + security + cost |
                                            transaction_cost + age + gender + 
                                            work_employed + race_white + race_black + 
                                            hispanic + year2019 + year2018 + year2017 +
                                            year2016 + married + yrs_of_schooling + 
                                            household_income + hh_size | 1,
                           rpar = c("(Intercept):Check" = "n",
                                    "(Intercept):Electronic" = "n"),
                           reflevel = "Cash",
                           data = mlogit_cpc_v1,
                           halton = NA, # using halton sequences
                           panel = TRUE, # accounting for individual heterogeneity
                           R = 2000) 
print(summary(mixed_multi_noRC))
print(AIC(mixed_multi_noRC))

# Only if you want to see the other specifications
if (run_other_specifications == TRUE) {
  
  # Estimation (1)
  mixed_multi <- mlogit(payment_choice ~ convenience + security + cost |
                                         transaction_cost + age + gender +
                                         work_employed + race_white + race_black +
                                         hispanic + year2019 + year2018 + year2017 +
                                         year2016 + married + yrs_of_schooling +
                                         household_income + hh_size | 1,
                         rpar = c("convenience" = "n",
                                  "security" = "n",
                                  "cost" = "n"
                                  ),
                         reflevel = "Cash",
                         data = mlogit_cpc_v1,
                         halton = NA, # using halton sequences
                         panel = TRUE, # accounting for individual heterogeneity
                         R = 2000) 
  print(summary(mixed_multi)) 
  print(AIC(mixed_multi))
  
  # Estimation (2)
  mixed_multi.wCorr <- update(mixed_multi,
                              correlation = TRUE)
  print(summary(mixed_multi.wCorr)) 
  print(AIC(mixed_multi.wCorr))
  
  # Estimation (4)
  mixed_multi_noRC.wCorr <- update(mixed_multi_noRC, 
                                   correlation = TRUE)
  print(summary(mixed_multi_noRC.wCorr))
  print(AIC(mixed_multi_noRC.wCorr))
  
  # The predicted shares
  predShares <- fitted(mixed_multi, outcome = F) %>%
                  apply(., 2, 
                        function(y) {
                          
                          return(mean(x = y, na.rm = TRUE))
                          
                          })
  print(predShares)
  
  # Selecting the model which minimizes the AIC
  model_list <- list(mixed_multi,
                     mixed_multi.wCorr,
                     mixed_multi_noRC,
                     mixed_multi_noRC.wCorr)
  
  save(list = c("model_list"),
       file = str_c("mixed_models_", Sys.Date() ,".Rdata"))
  
  # Choosing the mixed logit which minimizes the AIC
  using_mixed_multi <- model_list[[which.min(unlist(purrr::map(model_list, AIC)))]]
  summary(using_mixed_multi)
  coef(using_mixed_multi)
  
  } else {
    
    using_mixed_multi <- mixed_multi_noRC
    summary(using_mixed_multi)
    coef(using_mixed_multi)
    
}

#==========================================================================================
                          # AVERAGE MARGINAL EFFECT ESTIMATIONS #
#==========================================================================================

# Coefficient names plus manually creating the coefficient names
end <- length(using_mixed_multi$model)
coef_names <- names(using_mixed_multi$model)[-c(1, end - 1, end)]

first_var <- which(coef_names == "transaction_cost")
last_var <- which(coef_names == "hh_size")
coef_names <- coef_names[first_var:last_var]

# A quick way to compute the marginal effects for each variable. 
p2p_ME <- purrr::map(.x = as.list(coef_names),
                     .f = function(x) {
                       
                       me_calc <- x %>%
                         stats::effects(using_mixed_multi,
                                        covariate = .,
                                        type = "aa",
                                        data = mlogit_cpc_v1) %>%
                         data.frame() %>%
                         summarise_all(.tbl = .,
                                       .funs = ~ mean(x = ., na.rm = TRUE))

                       return(me_calc)
                       
                       }) %>%
            bind_rows() 
rownames(p2p_ME) <- coef_names
print(p2p_ME)

p2p_ME %>%
  rownames_to_column(.data = ., 
                     var = "Covariate") %>%
  mutate(Covariate = factor(Covariate, 
                            levels = c("transaction_cost", "year2016", "year2017", 
                                       "year2018", "year2019", "age", "yrs_of_schooling", 
                                       "gender", "race_black", "race_white", "hispanic", 
                                       "household_income", "hh_size", "married", 
                                       "work_employed"))) %>%
  arrange(Covariate) %>%
  xtable(digits = 3) %>%
  print(include.rownames = FALSE)

#==========================================================================================
                          # MARGINAL EFFECT STANDARD ERROR ESTIMATIONS #
#==========================================================================================

# This section will compute the SE and CI for the marginal effects of the demograhic
# and transaction variables in the discrete choice model.

# Computing the matrix of marginal effects point estimates.
for (c in 1:length(coef_names)) {
  
  if (c == 1) me_matrix <- matrix(NA, 
                                  nrow = length(coef_names), 
                                  ncol = 4) %>% 
                            data.frame()
  
  # Coefficient names plus manually creating the coefficient names
  end <- length(using_mixed_multi$model)
  coef_names <- c(names(using_mixed_multi$model)[-c(1, end - 1, end)])
  
  last_var <- which(coef_names == "hh_size")
  coef_names <- coef_names[4:last_var]

  coefficient <- coef_names[c]
  me <- effects(using_mixed_multi,
                covariate = coefficient,
                data = mlogit_cpc_v1) %>%
    na.omit() %>%
    apply(., 2, mean)
  
  me_matrix[c, 2:4] <- me
  
}

me_matrix[,1] <- coef_names
colnames(me_matrix) <- c("Covariate", "Cash", "Check", "Electronic")

# I'm leaving this function here instead of putting it in the "Function" section 
# for pedalogical purposes. 
AME_function <- function(thetas){
  
  temp <- using_mixed_multi
  temp$coefficients <- thetas
  end <- length(using_mixed_multi$model)
  coef_names <- c(names(using_mixed_multi$model)[-c(1, end - 1, end)])
  
  last_var <- which(coef_names == "hh_size")
  coef_names <- coef_names[4:last_var]

  # A quick way to compute the marginal effects for each variable. 
  mixl.ME <- sapply(coef_names,
                    function(x){
                      
                      stats::effects(temp,
                                     covariate = x,
                                     data = mlogit_cpc_v1) %>%
                        # There are some NaN electronic random coefficients
                        na.omit() %>% 
                        apply(., 2, mean)
                      
                    },
                    simplify = F)

  me_vector <- mixl.ME %>%
                unlist() %>%
                data.frame(AME = .) %>%
                data.matrix()
  
  return(me_vector)
  
}

# Computing the gradient of the function.
grad <- jacobian(AME_function,
                 using_mixed_multi$coefficients[1:35])

save(list = c("grad"), 
     file = str_c("ame_numGradient_", Sys.Date(), ".Rdata"))

#==========================================================================================

# Matrix of standard errors which are computed using the delta method
multinomial_me_se <- diag(grad %*% vcov(using_mixed_multi)[1:35, 1:35] %*% t(grad)) %>% 
                      sqrt() %>%
                      matrix(., 
                             ncol = 3, 
                             byrow = T)

me_matrix_se <- data.frame(Covariate = coef_names,
                           Cash.SE = multinomial_me_se[, 1],
                           Check.SE = multinomial_me_se[, 2],
                           Electronic.SE = multinomial_me_se[, 3],
                           stringsAsFactors = F)

marginal_effects_data <- inner_join(me_matrix, 
                                    me_matrix_se, 
                                    by = "Covariate")

# Creating a new marginal effects matrix for each banking status variable 
for (b in c("Cash", "Check", "Electronic")) {
  
  temp_data <- marginal_effects_data %>%
    dplyr::select(Covariate, starts_with(b))
  
  temp_data <- temp_data %>%
    mutate(lower.99ci = temp_data[, 2] - (2.576 * temp_data[,3]),
           upper.99ci = temp_data[, 2] + (2.576 * temp_data[,3]),      
           lower.95ci = temp_data[, 2] - (1.96 * temp_data[,3]),
           upper.95ci = temp_data[, 2] + (1.96 * temp_data[,3]),
           lower.90ci = temp_data[, 2] - (1.64 * temp_data[,3]),
           upper.90ci = temp_data[, 2] + (1.64 * temp_data[,3]))
  
  assign(str_c(b,"_me_data"), temp_data)
  
}


sig_gen <- function(data) {
  
  Variables <- data.frame(Covariate = data$Covariate, 
                          stringsAsFactors = )
  new_data <- data %>%
      dplyr::select(-Covariate) %>%
      round(digits = 3) %>%
      mutate(significance = ifelse((lower.90ci < 0 & upper.90ci < 0) | 
                                   (lower.90ci > 0 & upper.90ci > 0), 
                                   "*", 
                                   "")) %>%
    
      mutate(significance = ifelse((lower.95ci < 0 & upper.95ci < 0) | 
                                   (lower.95ci > 0 & upper.95ci > 0), 
                                   "**", 
                                   significance)) %>%
    
      mutate(significance = ifelse((lower.99ci < 0 & upper.99ci < 0) | 
                                   (lower.99ci > 0 & upper.99ci > 0), 
                                   "***", 
                                   significance)) 
  
  output <- cbind(Variables, new_data)  
  
  return(output)
  
}

marginal_effects <- list(Cash_me_data, 
                         Check_me_data, 
                         Electronic_me_data) %>%
                      lapply(., 
                             sig_gen)
me_frame <- bind_cols(marginal_effects[[1]],
                      marginal_effects[[2]],
                      marginal_effects[[3]]) %>%
              dplyr::select(-c(starts_with("lower"), starts_with("upper"))) %>%
  mutate(Covariate...1 = factor(Covariate...1, 
                                levels = c("transaction_cost", "year2016", "year2017", 
                                           "year2018", "year2019", "age",
                                           "yrs_of_schooling", "gender", "race_black", 
                                           "race_white", "hispanic", "household_income", 
                                           "hh_size", "married", "work_employed"))) %>%
  arrange(Covariate...1)
print(me_frame)

