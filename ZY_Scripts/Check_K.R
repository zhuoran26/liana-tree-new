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

#modified:
data.original =read.csv('/Users/zhuoranyu/Desktop/liana-tree-test/data/full_met_analysis_data_8_Feb.csv')
###
data$Units_K = as.character(data$Units_K)

# Convert to consistent units
for(i in 1:nrow(data)){
  if(data$Units_K[i] == 'kg/m/s/Mpa'){
    # Convert from kg to mmol H2O
    data$K[i] = data$K[i] * 1000 / 18 * 1000 #divided by 18 is to turn to mol, first *1000 is kg to g, second *1000 is mol to mmol
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
data$K.mol.m.s.MPa = data$K * 0.001 #here change everything to mol!

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

# Violin and boxplots overlaid with jittered points
pl1 = data %>%
  dplyr::select(Growth.form, K.mol.m.s.MPa) %>%
  ggplot(aes(x = Growth.form, y = K.mol.m.s.MPa, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', shape = 3, show.legend = F, size = 1.5) +
  theme_linedraw() +
  xlab('') + ylab((bquote(K['s,max']~(mol~m^-1~s^-1~MPa^-1)))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 51', 'Tree\nn = 103')) +
  ggtitle('Hydraulic conductivity') +
  theme(plot.title = element_text(hjust = 0.5, size = 10), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pl1
#modified

liana.data <- subset(data, Growth.form == 'liana')
liana.k.in.plot <- liana.data$K.mol.m.s.MPa
liana_K_mmol <- liana.data$K

#NOTE: the first violin plot used K.mol.m.s.MPa

#Check with K I have used in NPP model

load('/Users/zhuoranyu/Desktop/liana-tree-test/data/param.input.RData')

# Load met drivers
load('/Users/zhuoranyu/Desktop/liana-tree-test/data/bci_met_mxh.RData')

nK = 75000
ks = c(trees_K$K / 10, lianas_K$K)
#what is lianas_K$K here
lianas_K$K
liana_K_mmol

varyk = seq(min(ks), max(ks), length.out = nK)
range(varyk)

#So the lianas_K$K is all the 51 dif liana values in mmol/m/s/Mpa


