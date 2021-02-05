#================================================================================================
                                        # PREAMBLE #
#================================================================================================

# Written by: Brian Prescott

# Purpose: Calibrate the model of p2p random matching described in the paper 
# How People Pay Each Other: Data, Theory, and Calibrations by 
# Claire Greene, Brian Prescott, and Oz Shy

# Notes: The model relies on the solver `Rsolnp`. The data used to solve the model is precisely
# the same data used to estimate the mixed multinomial logit. 

#================================================================================================
                                        # LIBRARIES #
#================================================================================================

library(plyr)  # Data manipulation
library(tidyverse) # Data manipulation
library(haven) # Import Stata files
library(Rsolnp) # optimization routine
library(xtable) # LaTeX tables
library(gridExtra) # ordering the calibration results
library(readxl) # reading in excel files (frba colors)
library(ggsci) # coloring for ggplots
library(Hmisc) # for the weighted variance function

#================================================================================================
                                #  FILE PATHS & PARAMETERS # 
#================================================================================================

# Change your wd here
setwd("A:/papers/p2p paper/code/p2p-wp-code/p2p-wp-materials/")

#================================================================================================
                                  # GENERAL FUNCTIONS #
#================================================================================================

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
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
  
  if (numPlots==1) {
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

not.nan <- function(x) { 
  
  return(!is.nan(x)) 
  
}


#================================================================================================
                                    # DATA IMPORT #
#================================================================================================

p2p_cs.df <- read_rds("p2p_df.rds")

#================================================================================================
                                # CALIBRATION FUNCTIONS #
#================================================================================================

# These are the functions used to calibrate the theoretical model's predictions. 

# The optimization objective function. 
# p corresponds to the calibrated values
# q corresponds to the data values

rmse <- function(p, q) {
  
  output <- sqrt(mean((p - q)^2)) 
  
  return(output)
  
}

euclidean_distance <- function(p, q) { 
  
  output <- sqrt(sum((p - q)^2)) 
  
  return(output)
  
}

calibrate_model <- function(calibration_values, metric) {
  
  set.seed(10)
  
  calib_values <- calibration_values
  num_pars <- 6

  # The calibration objective function which finds the parameter values to 
  # minimize the RMSE. 
  # Have to restate it each time to make sure that the new vphi is being 
  # used within it. It is defined more clearly above. Collapse this for ease of use
  calibration_function <- function(par) {
          
          phi <- par
        
          for (i in 1:num_pars) { 
            
            assign(str_c("phi_", i), phi[i])
            
          }
  
          mean_C_calib <- phi_1 * (phi_1 + phi_2 + phi_3 + phi_5) + 
                          phi_2 * (phi_2 + phi_3 + phi_5) + 
                          (phi_3 * phi_5)
          
          mean_K_calib <- phi_2 * (phi_4 + phi_6) + 
                          phi_4 * (phi_5 + phi_6) + 
                          phi_5 * (phi_5 + phi_6) + 
                          phi_6^2
          
          mean_E_calib <- phi_1 * (phi_4 + phi_6) +
                          phi_3 * (phi_3 + phi_4 + phi_6) + 
                          phi_4^2
          
          var_C_calib <- mean_C_calib * (1 - mean_C_calib)
          
          var_K_calib <- mean_C_calib * (1 - mean_K_calib)
          
          var_E_calib <- mean_C_calib * (1 - mean_E_calib)

          mean_C_data <- calib_values[1]
          mean_K_data <- calib_values[2]
          mean_E_data <- calib_values[3]
          var_C_data <- calib_values[4]
          var_K_data <- calib_values[5]
          var_E_data <- calib_values[6]
          
    
          calib_input <- rbind(mean_C_calib, # these are the calibrated values
                               mean_K_calib, 
                               mean_E_calib,
                               var_C_calib, 
                               var_K_calib, 
                               var_E_calib
                               )
          
          data_input <- rbind(mean_C_data, # Data implied value
                              mean_K_data, # Data implied value
                              mean_E_data, # Data implied value
                              var_C_data, # Data implied value
                              var_K_data, # Data implied value
                              var_E_data # Data implied value
                              ) 
          
          # The objective to minimize
          if (metric == "euclidean") {
            
            distance <- euclidean_distance(calib_input, data_input)
            
          } else if (metric == "rmse") {
            
            distance <- rmse(calib_input, data_input)
            
          }
          
          return(distance)
          
        }
  
        # An equality constrained minimization problem
        # Assigning each element in phi an equal initial value
  
        calibration <- Rsolnp::solnp(pars = rep(1 / num_pars, num_pars),
                                     fun = calibration_function,
                                     eqfun = function(phi) {
                                       
                                       return(sum(phi))
                                       
                                       },
                                     eqB = 1,
                                     LB = rep(0.00000001, num_pars),
                                     UB = rep(0.95, num_pars),
                                     # Setting the tolerance parameters
                                     control = list(tol = 1e-16, 
                                                    outer.iter = 1000, 
                                                    inner.iter = 1600))
        
        return(calibration)
        
}      

#================================================================================================
                                      # MODEL CALIBRATION #
#================================================================================================

# Turning off scientific notation
options(scipen = 999)

# Computing the sample moment values
input_calib_values <- p2p_cs.df$payment_choice %>%
                        table() %>%
                        prop.table() %>%
                        c(., 
                          p2p_cs.df %>%
                            mutate(used_cash = as.numeric(payment_choice == "Cash"), 
                                   used_check = as.numeric(payment_choice == "Check"), 
                                   used_electronic = as.numeric(payment_choice == "Electronic")
                                   ) %>% 
                            summarize_at(.tbl = ., 
                                         .vars = vars(starts_with("used_")), 
                                         .funs = ~ var(x = ., na.rm = TRUE))) %>%
                        unlist() %>%
                        as.numeric()

calibrated_p2p_model <- calibrate_model(calibration_values = input_calib_values, 
                                        metric = "euclidean")
print(calibrated_p2p_model$pars)


