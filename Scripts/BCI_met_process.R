## This script processes meteorological inputs for Barro Colorado Island PA
## This site represents intermediate dryness, or wetter tropical seasonal forest
## Required met inputs for the model are vapor pressure deficit (VPD) and soil water potential (SWP)

## Author: AM Willson
## Date modified: 07 July 2020

rm(list = ls())
library(lubridate)

#########
## VPD ##
#########

## Data to calculate VPD are publically available from STRI website
## Link: https://biogeodb.stri.si.edu/physical_monitoring/research/barrocolorado

## VPD is calculated from relative humidity (RH) and air temperature (AT)
## Temperature and relative humidity data are for the top of the canopy 
## (i.e. tower extending higher than the top of the canopy)
## I chose top of canopy because we are assuming one layer of leaves, 
## such that all leaves would be exposed to the upper most part of the canopy

# Load in data
# These data are unavailable but the processed data that run the model are available
bci_temp = read.csv('bci_elect_48m_at/bci_lutz48m_at_elect.csv')
bci_rh = read.csv('bci_elect_48m_rh/bci_lutz48m_rh_elect.csv')

# Select only required columns: 
## time, air temperature or relative humidity, QAQC
bci_temp = subset(bci_temp, select = c('datetime','at', 'chk_note'))
bci_rh = subset(bci_rh, select = c('datetime', 'rh', 'chk_note'))

# Remove rows failing QAQC
bci_temp = subset(bci_temp, bci_temp$chk_note != 'bad' 
                  & bci_temp$chk_note != 'doubtful' 
                  & bci_temp$chk_note != 'nc' 
                  & bci_temp$chk_note != 'missing')
bci_rh = subset(bci_rh, bci_rh$chk_note != 'bad' 
                & bci_rh$chk_note != 'doubtful' 
                & bci_rh$chk_note != 'nc' 
                & bci_rh$chk_note != 'missing')

# Remove QAQC columns
bci_temp = subset(bci_temp, select=  c('datetime', 'at'))
bci_rh = subset(bci_rh, select = c('datetime', 'rh'))

# Change datetime column to get dates in consistent format
bci_temp$dt = dmy_hms(bci_temp$datetime)
bci_rh$dt = dmy_hms(bci_rh$datetime)

# Remove datetime columns (keeping formatted date column)
bci_temp = subset(bci_temp, select = c('at', 'dt'))
bci_rh = subset(bci_rh, select = c('rh', 'dt'))

# Join data frames
bci_temp_rh = bci_temp %>%
  full_join(bci_rh, by = 'dt')

# Remove rows that are missing at or rh
bci_temp_rh = subset(bci_temp_rh, 
                     !is.na(bci_temp_rh$at) 
                     & !is.na(bci_temp_rh$rh))

# Calculate SWP and VPD
bci_temp_rh$SVP_hPa = 610.78 * exp(bci_temp_rh$at / (bci_temp_rh$at + 238.3) * 17.2694) / 100
bci_temp_rh$VPD_hPa = bci_temp_rh$SVP_hPa * (1 - bci_temp_rh$rh / 100)
bci_temp_rh$VPD_Pa = bci_temp_rh$VPD_hPa * 100

# Storage: month x hour matrix
bci_vpd_matrix = matrix(, nrow = 12, ncol = 24)

# Find mean VPD in format month x hour 
# (averaging over days of the month & years)
for(i in 1:12){
  for(j in 0:23){
    bci_vpd_matrix[i,j+1] = mean(bci_temp_rh$VPD_Pa[which(as.numeric(substr(bci_temp_rh$dt, 6, 7)) == i 
                                                          & as.numeric(substr(bci_temp_rh$dt, 12, 13)) == j)], na.rm = T)
    print(paste('Done with hour',j,'in month',i))
  }
  print(paste('Done with month', i))
}

# But check NA vals and make sure they correspond to times that don't exist
NA_dim = which(is.na(bci_vpd_matrix), arr.ind = T)

# Rename
VPD = bci_vpd_matrix

#########
## SWP ##
#########

# Read in SWP (in MPa) data  (from Levy-Varon et al. 2019)
# Note that these data are not available but the processed data to run the model are
data_swp = read.csv('SWP_BCI.csv', skip = 1)


# Change column names
colnames(data_swp) = c('Month', 'D_550_cm', 'D_450_cm', 'D_350_cm', 'D_250_cm', 
                       'D_175_cm', 'D_125_cm', 'D_90_cm', 'D_70_cm', 'D_50_cm', 
                       'D_35_cm', 'D_25_cm', 'D_15_cm', 'D_7.5_cm', 'D_2.5_cm')

## Most water used by woody plants comes in forested settings comes from the top soil layers
## Therefore, we are only interested in SWP from the top layers
## However, between sites, the layers for which data are available are different
## Therefore, we will consider only the 15 cm depth, 
## which gives a reasonable approximation of the soil water available 
## for the shallow soil profile avialable to plant roots

# Keep only 15 cm depth
data_swp = data_swp[,13]

SWP = data_swp

# Save output
save(VPD, SWP, file = 'bci_met_mxh.RData')
