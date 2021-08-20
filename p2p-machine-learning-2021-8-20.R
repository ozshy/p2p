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
setwd("~/Papers/p2p/p2p_coding") 
dir()
#=================================================================================================
                                        # DATA IMPORT #
#=================================================================================================

# This is just reading in the p2p dataset from 2015 to 2019 for the machine learning model.
#p2 <- read_rds("p2p_ML_df.rds")# old data set
p2 <- read_rds("p2p_df_2021-08-19.rds")# revised data set

# p2 = p1
dim(p2)

(num_p2p = nrow(p2))# num trans
(num_resp = length(unique(p2$id)))# num resp who p2p
num_p2p/num_resp # num p2p per resp from 2015 to 2019

names(p2)

# checking dist of Amount (whether log is needed)
p2$Amount = p2$amnt
# hist(subset(p2$Amount, p2$Amount < 1000))
# hist(subset(p2$Amount, p2$Amount < 200))
# hist(subset(p2$Amount, p2$Amount < 150))
# hist(subset(p2$Amount, p2$Amount < 100))
# nrow(subset(p2, Amount <= 500))
# nrow(subset(p2, Amount <= 300))
# nrow(subset(p2, Amount <= 200))
# nrow(subset(p2, Amount <= 150))
# nrow(subset(p2, Amount <= 100))

#=================================================================================================
                                    # MACHINE LEARNING ANALYSIS #
#=================================================================================================

### Classification tree begins
# defining regression model
# check var
p2$Method = p2$payment_choice
table(p2$Method)
#
table(p2$in_person)
#
p2$Age = p2$age
summary(p2$Age)
#
p2$HH_size = p2$hh_size
table(p2$HH_size)
#
table(p2$work_employed)
p2$Work = NA
p2$Work[p2$work_employed==1] = "Employed"
p2$Work[p2$work_employed==0] = "Not_employed"
table(p2$Work)
#
table(p2$married)
p2$Marital = NA
p2$Marital[p2$married==1] = "Married"
p2$Marital[p2$married==0] = "Not_married"
table(p2$Marital)
#
summary(p2$hh_income)# need original
p2$HH_income = p2$hh_income
#
table(p2$gender)
p2$Gender = NA
p2$Gender[p2$gender==1] = "Male"
p2$Gender[p2$gender==0] = "Female"
table(p2$Gender)
#
table(p2$highest_education)# NEED origital hh_income
p2$Education = NA
p2$Education[p2$highest_education <= 8] = "elem_or_less"
p2$Education[p2$highest_education >= 9 & p2$highest_education < 11] = "high_school"
p2$Education[p2$highest_education >= 11 & p2$highest_education < 14] = "assoc_or_college"
p2$Education[p2$highest_education >= 14] = "MA_or_higher"
table(p2$Education)

table(p2$ind_payee)# p2p payee
p2$Payee = NA
p2$Payee[p2$ind_payee==1] = "Goods/services (bus)"
p2$Payee[p2$ind_payee==2] = "Goods/services (non_bus)"
p2$Payee[p2$ind_payee==3] = "Family/friend"
p2$Payee[p2$ind_payee==4] = "Coworker/classmate"
p2$Payee[p2$ind_payee==5] = "Other"
table(p2$Payee)

names(select(p2, contains("context")))
p2$Context = NA
table(p2$context_goodServ)
p2$Context[p2$context_goodServ==1] = "GoodServ"
table(p2$context_repay)
p2$Context[p2$context_repay==1] = "Repay"
table(p2$context_lend)
p2$Context[p2$context_lend==1] = "Lend"
table(p2$context_donation)
p2$Context[p2$context_donation==1] = "Donation"
table(p2$Context)

table(p2$in_person)
p2$In_person = NA
p2$In_person[p2$in_person==1] = "Yes"
p2$In_person[p2$in_person==0] = "No"
table(p2$In_person)

# Define regression tree model
#method_on_demog = Method ~ Amount + In_person + Context + Payee + Age + HH_size + Work + Marital + HH_income + Gender + Education

# we removed + Context + Payee because they never appear as top predictors
method_on_demog = Method ~ Amount + In_person + Age + HH_size + Work + Marital + HH_income + Gender + Education 

# Extremely-long tree first, then prune (actually, snip it!) it
method_demog_tree = rpart(method_on_demog, 
                          data = p2, 
                          method = "class", 
                          control = rpart.control(cp = 0.006))

# plot full tree (before snipping or pruning)
rpart.plot(method_demog_tree, type = 5, extra = 100, legend.x=NA, legend.y=NA, tweak = 1.1, fallen.leaves = FALSE, gap = 0, space = 1, digits = 4, compress = T, ycompress = F)

# Select the plot that suits the paper size based on the no of splitsmethod_demog_tree$cptable
plotcp(method_demog_tree)
method_demog_tree$cptable

# Now prune the tree 
# Note, you may have to run it more than once to get a nice tree 
# (first time gets distorted for some reason)
# 0.01234568 - this was Oz's original complexity parameter value. 
set.seed(1955)
method_demog_prune = prune.rpart(method_demog_tree,                                  cp = 0.006) # => x splits 

rpart.plot(method_demog_prune, type = 5, extra = 100, legend.x=NA, legend.y=NA, tweak = 1.1, fallen.leaves = FALSE, gap = 0, space = 1, digits = 4, compress = T, ycompress = F)
# Use tweak to manage font size 
# Export to PDF from R-Studio, use 4x6 or 5x7 device (better), landscape 

# for the caption of Figure 1 in the paper
nrow(p2)# num of payment obs
length(unique(p2$id))# num of respondents

### Unused code
#method_on_demog = Method ~ Amount + In_person + Context  + Age + HH_size + Work + Marital + HH_income + Gender + Education # payee removed

#method_on_demog = Method ~ Amount + In_person + Payee + Age + HH_size + Work + Marital + HH_income + Gender + Education # context removed

#method_on_demog = Method ~ Amount + In_person + Context  + Age + HH_size + Work + Marital + HH_income + Gender #+ Education # payee and education removed

###
# prp(method_demog_prune, 
#     type = 3, 
#     box.palette = list("seagreen3", "lightcoral", "skyblue1"), 
#     extra = 100, 
#     under = TRUE, 
#     tweak = 1.25, 
#     varlen = 0, 
#     faclen = 0, 
#     legend.cex = 1,
#     legend.x = 0.85,
#     legend.y = 1.05,
#     clip.right.labs = FALSE,
#     Margin = 0.0, 
#     digits = -2) # faclet=0 avoids abbreviations, tweak for char size

###
#
#p2$ln.Amount <- log(p2$Amount)
#p2$year2015 <- as.numeric(p2$year == 2015)
#p2$year2016 <- as.numeric(p2$year == 2016)
#p2$year2017 <- as.numeric(p2$year == 2017)
#p2$year2018 <- as.numeric(p2$year == 2018)
#p2$year2019 <- as.numeric(p2$year == 2019)

###
# snip the tree (instead of pruning)
#small_tree <- snip.rpart(method_demog_tree, toss = 5) # toss is num of splits
# plot the snipped tree
#rpart.plot(small_tree, type = 5, extra = 100, legend.x=NA, legend.y=NA, tweak = 1.1, fallen.leaves = FALSE, gap = 0, space = 1, digits = 4, compress = T, ycompress = F)
