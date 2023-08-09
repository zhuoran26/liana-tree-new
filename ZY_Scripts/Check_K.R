#checked what K is used in these codes
library(dplyr)
library(ggplot2)
library(ggsci)
library(pracma)
library(reshape2)
library(cowplot)
library(DescTools)
library(gridExtra)
library(directlabels)
library(raster)
library(geomtextpath)

##################
#### Figure 1 ####
##################

# Clean environment and load new data
rm(list = ls())
data = read.csv('/Users/zhuoranyu/Desktop/liana-tree-test/data/full_met_analysis_data_8_Feb.csv')

head(data)
data$Units_K = as.character(data$Units_K)

# Convert to consistent units
for(i in 1:nrow(data)){
  if(data$Units_K[i] == 'kg/m/s/Mpa'){
    # Convert from kg to mmol H2O
    data$K[i] = data$K[i] * 1000 / 18 * 1000
    data$Units_K[i] = 'mmol/m/s/Mpa'
  }else{
    if(data$Units_K[i] == 'mol/m/s/Mpa'){
      # Convert from mol to mmol H2O
      data$K[i] = data$K[i] * 1000
      data$Units_K[i] = 'mmol/m/s/Mpa'
    }else{
      print('Unknown units')
      print(data$Units_K[i])
    }
  }
}

data$Units_K = as.factor(data$Units_K)

# Make column for mol/m/s/MPa
data$K.mol.m.s.MPa = data$K * 0.001

# Remove P50 > -0.75 consistent with Trugman et al. 2020
data = data %>%
  mutate(P50, replace(P50, P50 > -0.75, NA)) %>%
  dplyr::select(-P50)

# Rename column with NAs
colnames(data)[ncol(data)] = 'P50'

# Number of observations for each growth form for each trait (now including slope)
length(which(data$Growth.form == 'liana' & !is.na(data$K.mol.m.s.MPa)))
length(which(data$Growth.form == 'tree' & !is.na(data$K.mol.m.s.MPa)))
length(which(data$Growth.form == 'liana' & !is.na(data$P50)))
length(which(data$Growth.form == 'tree' & !is.na(data$P50)))
length(which(data$Growth.form == 'liana' & !is.na(data$Slope)))
length(which(data$Growth.form == 'tree' & !is.na(data$Slope)))

liana.dif.K <- subset(data, Growth.form == 'liana')
liana.dif.K.Mpa <- liana.dif.K$K.mol.m.s.MPa

summary(liana.dif.K.Mpa)

#NOTE: the first violin plot used K.mol.m.s.MPa

#Check with K I have used in NPP model

load('/Users/zhuoranyu/Desktop/liana-tree-test/data/param.input.RData')

# Load met drivers
load('/Users/zhuoranyu/Desktop/liana-tree-test/data/bci_met_mxh.RData')

nK = 75000
ks = c(trees_K$K / 10, lianas_K$K)
#what is lianas_K$K here
lianas_K$K


varyk = seq(min(ks), max(ks), length.out = nK)
range(varyk)
#make units the same
liana.dif.K.kpa <-liana.dif.K.Mpa*1000
range(varyk)
range(liana.dif.K.Mpa)
summary(lianas_K$K)
length(lianas_K$K)
summary(liana.dif.K.kpa)
length(liana.dif.K.kpa)

#So the lianas_K$K is all the 51 dif liana values


