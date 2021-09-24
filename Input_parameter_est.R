## This script collects and processes all inputs to the tree and liana NPP models
## Brief description of data stream sources are commented in
## These inputs are used to calculate NPP from the models in NPP_models.R

## Author: AM Willson
## Date modified: 09 February 2021

rm(list = ls())

###############################
## Load and prepare the data ##
###############################

## Load conductivity meta-analysis data
K_data = read.csv('full_met_analysis_data_8_Feb.csv')
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

# Subset observations for angiosperms (= trees) for K
trees_K = subset(K_data, K_data$Growth.form.2 == 'angiosperm')
# Subset observations for lianas out of PFT dataframe
lianas_K = subset(K_data, K_data$Growth.form.2 == 'liana')

## Load allometry data from Horizontes, from Chris Smith-Martin (from New Phytologist 2020 paper)
hori.data = read.csv('Horizontes_liana_tree_biomass_Martin.csv')

# Subset out observations on juveniles
hori.data.juv = subset(hori.data, hori.data$MatureJuvenile == 'juvenile')
# Subset out observations on mature individuals
hori.data = subset(hori.data, hori.data$MatureJuvenile == 'mature')

# Subset out mature trees
hori.data.tree = subset(hori.data, hori.data$Growth_Form == 'tree')
# Subset out mature lianas
hori.data.liana = subset(hori.data, hori.data$Growth_Form == 'liana')

## Load in DBH data (for lianas only) from Chris Smith-Martin
dbh.data = read.csv('LianaSurveyMay2018DataForAlyssa.csv')
# Convert units
dbh.data$DiamMay2018_cm = dbh.data$DiamMay2018_mm * 0.1

###############################
## Hydraulic trait estimates ##
###############################

# From meta-analysis

# Summary statistics for K by growth form
k.min.tree = min(trees_K$K, na.rm = T)
k.max.tree = max(trees_K$K, na.rm = T)
tree.mean.K = mean(trees_K$K, na.rm = T)
tree.sd.K = sd(trees_K$K, na.rm = T)

k.min.liana = min(lianas_K$K, na.rm = T)
k.max.liana = max(lianas_K$K, na.rm = T)
liana.mean.K = mean(lianas_K$K, na.rm = T)
liana.sd.K = sd(lianas_K$K, na.rm = T)

# Summary statistics for remaining hydraulic parameters
# by growth form and collectively
b1.mean.tree = mean(trees_K$Slope, na.rm = T)
b2.mean.tree = mean(trees_K$P50, na.rm = T) 

b1.mean.liana = mean(lianas_K$Slope, na.rm = T)
b2.mean.liana = mean(lianas_K$P50, na.rm = T)

b1.mean.all = mean(K_data$Slope, na.rm = T)
b2.mean.all = mean(K_data$P50, na.rm = T)

###################
## DBH estimates ##
###################

# From Horizontes data

# DBH's for mature trees and lianas
hori.data.tree.dbh = hori.data.tree$StemSumDBH_mm
hori.data.liana.dbh = hori.data.liana$StemSumDBH_mm

# Make units consistent with model assumptions
hori.data.tree.dbh = hori.data.tree.dbh * 0.1
hori.data.liana.dbh = hori.data.liana.dbh * 0.1

# Summary statistics for trees
mean.hori.data.tree.dbh = mean(hori.data.tree.dbh, na.rm = T)
min.hori.data.tree.dbh = min(hori.data.tree.dbh, na.rm = T)
max.hori.data.tree.dbh = max(hori.data.tree.dbh, na.rm = T)

# From liana dataset (separate from tree dataset)

# Summary statistics for liana DBH
mean.data.liana.dbh = mean(dbh.data$DiamMay2018_cm)
sd.data.liana.dbh = sd(dbh.data$DiamMay2018_cm)
med.data.liana.dbh = median(dbh.data$DiamMay2018_cm)
low.data.liana.dbh = quantile(dbh.data$DiamMay2018_cm, probs = 0.025)
high.data.liana.dbh = quantile(dbh.data$DiamMay2018_cm, probs = 0.975)

###################################
## Canopy area estimates ##
###################################

# Canopy area for mature trees and lianas
hori.data.tree.can = hori.data.tree$TotalCanopyArea_cm2
hori.data.liana.can = hori.data.liana$TotalCanopyArea_cm2
# Make units consistent with model assumptions
hori.data.tree.can = hori.data.tree.can * 0.0001
hori.data.liana.can = hori.data.liana.can * 0.0001

# Average for a tree-liana pair
hori.data.tot.can = mean(c(mean(hori.data.tree.can, na.rm = T), mean(hori.data.liana.can, na.rm = T)))

#############################
## Height/length estimates ##
#############################

# From Horizontes data

# Height for mature trees
hori.data.tree.len = hori.data.tree$MaxStemLength_cm
# Make units consistent with model assumptions
hori.data.tree.len = hori.data.tree.len * 0.01

# Length for mature lianas
hori.data.liana.len = hori.data.liana$MaxStemLength_cm
# Make units consistent with model assumptions
hori.data.liana.len = hori.data.liana.len * 0.01

# Liana length comapred to tree height

# Comparitive liana length at minimum
min.hori.data.comp.len = min(hori.data.liana.len, na.rm = T) / min(hori.data.tree.len, na.rm = T)
# Comparative liana .length at mean
mean.hori.data.comp.len = mean(hori.data.liana.len, na.rm = T) / mean(hori.data.tree.len, na.rm = T)
# Comparative liana length at median
med.hori.data.comp.len = median(hori.data.liana.len, na.rm = T) / median(hori.data.tree.len, na.rm = T)
# Comparative liana length at maximum
max.hori.data.comp.len = max(hori.data.liana.len, na.rm = T) / max(hori.data.tree.len, na.rm = T)

#########################
## Allometry equations ##
#########################

## Estimating b1Ht and b2Ht with nonlinear model

# Subset the data from Horizontes to match DBH and H for each individual tree
allom_data = cbind(hori.data.tree$Individual_ID, hori.data.tree$StemSumDBH_mm * 0.1, hori.data.tree$MaxStemLength_cm * 0.01)
allom_data = as.data.frame(allom_data)  

# Define DBH-H relationship (given in model)
tree.b2h = function(x,b1Ht,b2Ht){b1Ht * x^b2Ht}

# Perform non-linear parameter estimation
tree.b2h.est = nls(V3~tree.b2h(V2,b1Ht,b2Ht), start=c(b1Ht=3.6,b2Ht=0.5), data = allom_data)
tree.b2h.sum = summary(tree.b2h.est)

# Extract parameter estimates from summary
tree.b1Ht = tree.b2h.sum$parameters[1,1]
tree.b2Ht = tree.b2h.sum$parameters[2,1]

##########################
## Save as RData object ##
##########################

save(list = ls(), file = 'param.input.RData')

