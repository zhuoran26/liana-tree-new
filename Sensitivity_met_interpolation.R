## This script interpolates observed vapor pressure deficit (VPD) and soil water potential (SWP) data
## From Horizontes (CR) and BCI (PA) to capture a wide range of plausible water scarcity conditions
## Horizontes and BCI represent extremes, with Horizontes having high water scarcity and BCI having low water scarcity
## Relative to neotropical forests on average in the 21st century

## Author: AM Willson
## Date modified: 20 August 2020

rm(list = ls())

############################
## Load site-specifc data ##
############################

# Load BCI data with dimensions month x hour
# Processed elsewhere
load('bci_met_mxh.RData')
BCI_SWP = SWP
BCI_VPD = VPD

# Load Horizontes data with dimensions month x hour
# Processed elsewhere
load('horizontes_met_mxh.RData')
H_SWP = SWP
H_VPD = VPD

########################################################
## Interpolate growing season length and VPD severity ##
########################################################
## This section interpolates data to run the month x hour model

# Storage
SWP_interp_100 = matrix(, nrow = 12, ncol = 100)
VPD_interp_100 = array(, dim = c(12, 24, 100))

# Soil water potential interpolation
# Dimensions: month x site
for(i in 1:12){
  SWP_interp_100[i,] = seq(from = BCI_SWP[i], to = H_SWP[i], length.out = 100)
}

# Vapor pressure deficit interpolation
# Dimensions: month x hour x site
for(i in 1:12){
  for(j in 1:24){
    VPD_interp_100[i,j,] = seq(from = BCI_VPD[i,j], to = H_VPD[i,j], length.out = 100)
  }
}

##########
## Save ##
##########

save(SWP_interp_100, VPD_interp_100, file = 'Interpolation_mxh_100.RData')