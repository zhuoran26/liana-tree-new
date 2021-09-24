## This script processes meterological inputs for Horizontes CR
## This site represents a dry tropical seasonal forest
## Required met inputs are vapor pressure deficit (VPD) and soil water potential (SWP)

## Author: AM Willson
## Date modified: 07 July 2020

rm(list = ls())
library(rhdf5)

#########
## VPD ##
#########

## VPD is an hourly data product from ECMWF reanalysis data product
## The input here is the derived product that was processed by D Medvigy

# Get attribute name from file
h5ls(file = 'vpd_output.h5')

# Read data into R
VPD_raw = h5read(file = 'vpd_output.h5', name = 'vpd.kPa')

# Fix time to be in CST
VPD_fix = matrix(, nrow = 12, ncol = 24)
VPD_fix[,19] = VPD_raw[,1]
VPD_fix[,20] = VPD_raw[,2]
VPD_fix[,21] = VPD_raw[,3]
VPD_fix[,22] = VPD_raw[,4]
VPD_fix[,23] = VPD_raw[,5]
VPD_fix[,24] = VPD_raw[,6]
VPD_fix[,1] = VPD_raw[,7]
VPD_fix[,2] = VPD_raw[,8]
VPD_fix[,3] = VPD_raw[,9]
VPD_fix[,4] = VPD_raw[,10]
VPD_fix[,5] = VPD_raw[,11]
VPD_fix[,6] = VPD_raw[,12]
VPD_fix[,7] = VPD_raw[,13]
VPD_fix[,8] = VPD_raw[,14]
VPD_fix[,9] = VPD_raw[,15]
VPD_fix[,10] = VPD_raw[,16]
VPD_fix[,11] = VPD_raw[,17]
VPD_fix[,12] = VPD_raw[,18]
VPD_fix[,13] = VPD_raw[,19]
VPD_fix[,14] = VPD_raw[,20]
VPD_fix[,15] = VPD_raw[,21]
VPD_fix[,16] = VPD_raw[,22]
VPD_fix[,17] = VPD_raw[,23]
VPD_fix[,18] = VPD_raw[,24]

# Dimensions months x hours
hori_VPD_matrix = VPD_fix

# Convert to Pa
hori_VPD_matrix = hori_VPD_matrix * 1000

# VPD = hori_VPD_array
VPD = hori_VPD_matrix

#########
## SWP ##
#########

# Read in SWP (in MPa) data  (from Medvigy et al. 2019)
data_swp = read.csv('SWP_CostaRica.csv', skip = 1)

# Change column names
colnames(data_swp) = c('Month', 'D_290_cm', 'D_270_cm', 'D_250_cm', 
                       'D_230_cm', 'D_210_cm', 'D_190_cm', 'D_170_cm', 
                       'D_150_cm', 'D_130_cm', 'D_110_cm', 'D_95_cm', 
                       'D_85_cm', 'D_75_cm', 'D_65_cm', 'D_55_cm', 
                       'D_45_cm', 'D_35_cm', 'D_25_cm', 'D_15_cm', 'D_5_cm')

## Most water used by woody plants comes in forested settings comes from the top soil layers
## Therefore, we are only intrested in SWP from the top layers
## However, between sites, the layers for which data are available are different
## Therefore, we will consider only the 15 cm depth, which gives a reasonable approximation of the soil water available for the shallow soil profile avialble to plant roots

# Keep only 15 cm depth
data_swp = data_swp[,20]

SWP = data_swp

save(VPD, SWP, file = 'horizontes_met_mxh.RData')
