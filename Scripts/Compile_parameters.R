## This script defines all parameters necessary for running the model
## Brief descriptions of data stream sources are commented in
## These inputs are used to calculate NPP from the NPP_models.R and NPP_models_wCO2.R scripts

## Author: AM Willson
## Date modified: 10 May 2022

rm(list = ls())

#### Hydraulic parameters ####

# The hydraulic parameters (b1 and b2) were taken from our extended meta-analysis
# of hydraulic functional traits.

# Load conductivity meta-analysis data
K_data = read.csv('data/full_met_analysis_data_8_Feb.csv')
K_data$Units_K = as.character(K_data$Units_K)

# Convert to consistent units
for(i in 1:nrow(K_data)){
  if(K_data$Units_K[i] == 'kg/m/s/Mpa'){
    # Convert from kg to mmol H2O
    K_data$K[i] = K_data$K[i] * 1000 / 18 * 1000
    K_data$Units_K[i] = 'mmol/m/s/Mpa'
  }else{
    if(K_data$Units_K[i] == 'mol/m/s/Mpa'){
      # Convert from mol to mmol H2O
      K_data$K[i] = K_data$K[i] * 1000
      K_data$Units_K[i] = 'mmol/m/s/Mpa'
    }else{
      print('Unknown units')
      print(K_data$Units_K[i])
    }
  }
}

K_data$Units_K = as.factor(K_data$Units_K)

# Remove P50 > -0.75 consistent with Trugman et al. 2020
K_data = K_data %>%
  mutate(P50, replace(P50, P50 > -0.75, NA)) %>%
  dplyr::select(-P50)

# Rename column with NAs
colnames(K_data)[ncol(K_data)] = 'P50'

# Mean for b1 and b2 parameters
b1.mean.all = mean(K_data$Slope, na.rm = T)
b2.mean.all = mean(K_data$P50, na.rm = T)

# Also, save hydraulic conductivity by growth form
# Subset observations for angiosperms (= trees) for K
trees_K = subset(K_data, K_data$Growth.form.2 == 'angiosperm')
# Subset observations for lianas out of PFT dataframe
lianas_K = subset(K_data, K_data$Growth.form.2 == 'liana')

#### DBH ####

# Tree DBH is from a survey of tropical trees in a second-growth
# forest plot in Guanacaste Costa Rica.
# The distribution of DBHs is available in Supplementary Figure 8

# Mean of tree DBHs from forest survey
mean.hori.data.tree.dbh = 18.19389

# Liana DBH is from an unpublished survey of tropical lianas in
# Guanacaste Costa Rica Conducted by Smith-Martin.
# The distribution of DBHs is available in Supplementary Figure 8

# Low and high liana DBHs from field survey
liana.dbh.low = 1.86
liana.dbh.high = 10.649

# Mean of liana DBHs from field survey
mean.data.liana.dbh = 2.652418

#### Tree allometry ####

# The allometric constants (b1Ht and b2Ht) for the tree growth form
# come from the same survey of trees as for DBH. Both height and
# DBH were used in a nonlinear model to estimate the parameters following
# the form height = b1Ht * DBH^b2Ht. Parameters were estimated
# with nonlinear weighted least-squares following the nls() function
# in stats.

tree.b1Ht = 3.063747
tree.b2Ht = 0.4554272

#### Save ####

save(list = ls(), file = 'data/param.input.RData')
