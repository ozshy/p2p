#================================================================================================
                                        # PREAMBLE #
#================================================================================================

# Written by: Oz Shy

# Purpose: Estimate the machine learning model presented in the paper How People Pay Each Other:
# Data, Theory, and Calibrations by Claire Greene, Brian Prescott, and Oz Shy.

# Notes: The model relies on the rpart packages.

#=================================================================================================
                                        # LIBRARIES #
#=================================================================================================

# Libraries used for the script: 
library(dplyr)
library(stringr) # for the `str_c()` function
library(readr)
library(haven)
library(purrr)
library(rpart)
library(rpart.plot)
library(partykit) # modifies rpart tree plot

#=================================================================================================
                                        # PARAMETERS #
#=================================================================================================

# Setting the working directory 
setwd("A:/Papers/P2P Paper/data") 

#=================================================================================================
                                        # DATA IMPORT #
#=================================================================================================

# This is just reading in the p2p dataset from 2015 to 2019 for the machine learning model.
p2 <- read_rds("p2p_ML_df.rds")

# p2 = p1
dim(p2)

(num_p2p = nrow(p2))# num trans
(num_resp = length(unique(p2$id)))# num resp who p2p
num_p2p/num_resp # num p2p per resp from 2015 to 2019

# checking dist of Amount (whether log is needed)
hist(subset(p2$Amount, p2$Amount < 1000))
hist(subset(p2$Amount, p2$Amount < 200))
hist(subset(p2$Amount, p2$Amount < 150))
hist(subset(p2$Amount, p2$Amount < 100))
nrow(subset(p2, Amount <= 500))
nrow(subset(p2, Amount <= 300))
nrow(subset(p2, Amount <= 200))
nrow(subset(p2, Amount <= 150))
nrow(subset(p2, Amount <= 100))

#=================================================================================================
                                    # MACHINE LEARNING ANALYSIS #
#=================================================================================================

### Classification tree begins
# defining regression model
# check var
table(p2$Method)
names(subset(p2, select = c(Method, Amount, in_person, Age, HH_size, Work,
                            Marital, HH_income, Gender, Education)))
#
p2$ln.Amount <- log(p2$Amount)
p2$year2015 <- as.numeric(p2$year == 2015)
p2$year2016 <- as.numeric(p2$year == 2016)
p2$year2017 <- as.numeric(p2$year == 2017)
p2$year2018 <- as.numeric(p2$year == 2018)
p2$year2019 <- as.numeric(p2$year == 2019)

method_on_demog = Method ~ Amount + in_person + Age + HH_size + 
                           # year2015 + year2016 + year2017 + year2018 + year2019 + 
                           Work + Marital + HH_income + Gender + Education # +Education #

# Extremely-long tree first, then prune it
method_demog_tree = rpart(method_on_demog, 
                          data = p2, 
                          method = "class", 
                          control = rpart.control(cp = 0.005))

# Select the plot that suits the paper size based on the no of splitsmethod_demog_tree$cptable
plotcp(method_demog_tree)

# Now prune the tree 
# Note, you may have to run it more than once to get a nice tree 
# (first time gets distorted for some reason)
# 0.01234568 - this was Oz's original complexity parameter value. 
set.seed(1955)
method_demog_prune = prune.rpart(method_demog_tree, 
                                 cp = 0.0094) # 3 splits 
prp(method_demog_prune, 
    type = 3, 
    box.palette = list("seagreen3", "lightcoral", "skyblue1"), 
    extra = 100, 
    under = TRUE, 
    tweak = 1.25, 
    varlen = 0, 
    faclen = 0, 
    legend.cex = 1,
    legend.x = 0.85,
    legend.y = 1.05,
    clip.right.labs = FALSE,
    Margin = 0.0, 
    digits = -2) # faclet=0 avoids abvreviations, tweak for char size
