## This file creates Figures 1-4 of the main text of the associated manuscript

## This script requires the following inputs:
    ## 1. full_met_analysis_data_12_July.csv: not given but can be
    ## recreated following the steps detailed in the Methods
    ## 2. param.input.RData: created in Input_parameter_est.R
    ## 3. bci_met_mxh.RData: created in BCI_met_process.R
    ## 4. horizontes_met_mxh.RData: created in Horizontes_met_process.R
    ## 5. Interpolation_mxh_100.RData: created in Sensitivity_met_interpolation.R
    ## 6. tree.NPP.mxh & liana.NPP.mxh: compiled in NPP_models.R
    ## 7. tree.NPP.mxh.ca & liana.NPP.mxh.ca: compiled in NPP_models_wCO2.R

## Author: AM Willson
## Date modified: 24 September 2021

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

##################
#### Figure 1 ####
##################

# Clean environment and load new data
rm(list = ls())
data = read.csv('full_met_analysis_data_8_Feb.csv')

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

# Violin and boxplots overlaid with jittered points
pl1 = data %>%
  dplyr::select(Growth.form, K.mol.m.s.MPa) %>%
  ggplot(aes(x = Growth.form, y = K.mol.m.s.MPa, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', shape = 3, show.legend = F, stroke = 1.2, size = 2) +
  theme_linedraw() +
  xlab('') + ylab(expression(paste(K[s],' (mol/m/s/MPa)'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 51', 'Tree\nn = 103')) +
  ggtitle('Hydraulic conductivity') +
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pl2 = data %>%
  dplyr::select(Growth.form, P50) %>%
  ggplot(aes(x = Growth.form, y = P50, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, stroke = 1.2, size = 2) +
  theme_linedraw() +
  xlab('') + ylab(expression(paste(P[50],' (MPa)'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 40', 'Tree\nn = 60')) +
  ggtitle(expression(P[50])) +
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pl3 = data %>%
  dplyr::select(Growth.form, Slope) %>%
  ggplot(aes(x = Growth.form, y = Slope, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, stroke = 1.2, size = 2) +
  theme_linedraw() +
  xlab('') + ylab('Slope (%/MPa)') +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 8', 'Tree\nn = 13')) +
  ggtitle('Slope of PLC curve') +
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Plot together
extendedga = grid.arrange(pl1, pl2, pl3, nrow = 1)

ggsave(extendedga, filename = 'Plots/extended_hydrotrait.jpeg', width = 11, height = 7, units = 'in')

##################
#### Figure 2 ####
##################

rm(list = ls())

# Load parameters
load('param.input.RData')

# Load met drivers
load('bci_met_mxh.RData')

nK = 75000
ks = c(trees_K$K / 10, lianas_K$K)
varyk = seq(min(ks), max(ks), length.out = nK)

# Interpolate DBH (and remove highest DBH)
hori.data.liana.dbh = hori.data.liana.dbh[hori.data.liana.dbh < 29]
dbhs = seq(min(hori.data.liana.dbh), max(hori.data.liana.dbh), length.out = 50)

## Tree: 90% canopy, mean DBH, BCI

treeouteach = c()
tree_out_r1 = c()

for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = j
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 2 * 100
    frac.tree.al = 0.9
    
    psis = SWP[j]
    
    D = VPD[j,]
    
    if(j == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach)
  if(sumtreeout > 0){
    tree_out_r1[1] = sumtreeout
    tree_out_r1[2] = K
    
    break
  }
}

##Liana: 10% canopy, varied DBH, BCI

lianaouteach = c()
liana_out_r1 = matrix(, nrow = length(dbhs), ncol = 2)

for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 2 * 100
      frac.liana.al = 0.1
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach)
    if(sumlianaout > 0){
      liana_out_r1[i,1] = sumlianaout
      liana_out_r1[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

## Tree: 60% canopy, mean DBH, BCI

treeouteach = c()
tree_out_r2 = c()

for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = j
    
    
    tot.al = 2 * 100
    frac.tree.al = 0.6
    
    psis = SWP[j]
    
    D = VPD[j,]
    
    if(j == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach)
  if(sumtreeout > 0){
    tree_out_r2[1] = sumtreeout
    tree_out_r2[2] = K
    
    break
  }
}

## Liana: 40% canopy, varied DBH, BCI

lianaouteach = c()
liana_out_r2 = matrix(, nrow = length(dbhs), ncol = 2)

for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 2 * 100
      frac.liana.al = 0.4
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach)
    if(sumlianaout > 0){
      liana_out_r2[i,1] = sumlianaout
      liana_out_r2[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

# Load new climate drivers (Horizontes)
load('horizontes_met_mxh.RData')

## Tree: 90% canopy, mean DBH, Horizontes
treeouteach = c()
tree_out_r3 = c()

for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = j
    
    
    tot.al = 2 * 100
    frac.tree.al = 0.9
    
    psis = SWP[j]
    
    D = VPD[j,]
    
    if(j == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach)
  if(sumtreeout > 0){
    tree_out_r3[1] = sumtreeout
    tree_out_r3[2] = K
    
    break
  }
}

## Liana: 10% canopy, varied DBH, Horizontes

lianaouteach = c()
liana_out_r3 = matrix(, nrow = length(dbhs), ncol = 2)

for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 2 * 100
      frac.liana.al = 0.1
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach)
    if(sumlianaout > 0){
      liana_out_r3[i,1] = sumlianaout
      liana_out_r3[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

## Tree: 60% canopy, mean DBH, Horizontes

treeouteach = c()
tree_out_r4 = c()

for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = j
    
    
    tot.al = 2 * 100
    frac.tree.al = 0.6
    
    psis = SWP[j]
    
    D = VPD[j,]
    
    if(j == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach)
  if(sumtreeout > 0){
    tree_out_r4[1] = sumtreeout
    tree_out_r4[2] = K
    
    break
  }
}

## Liana: 40% canopy, varied DBH, Horizontes

lianaouteach = c()
liana_out_r4 = matrix(, nrow = length(dbhs), ncol = 2)

for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 2 * 100
      frac.liana.al = 0.4
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach)
    if(sumlianaout > 0){
      liana_out_r4[i,1] = sumlianaout
      liana_out_r4[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

# Format dataframe: 10% liana canopy at BCI
df1 = cbind(dbhs, liana_out_r1[,2], rep('Liana', times = length(dbhs)))
colnames(df1) = c('DBH', 'Kreq', 'PFT')
df1 = rbind(df1, c(mean.hori.data.tree.dbh, tree_out_r1[2], 'Tree'))
df1 = as.data.frame(df1)
df1$DBH = as.numeric(as.character(df1$DBH))
df1$Kreq = as.numeric(as.character(df1$Kreq))
df1$Kreq = df1$Kreq * 0.001
df1$sap_area = rep(NA, nrow(df1))
for(i in 1:nrow(df1)){
  df1$sap_area[i] = min(((pi*df1$DBH[i]^2) / 4), (2.41 * (df1$DBH[i]/2)^1.97)) #sapwood area in cm^2
}
df1$la = c(rep(0.1 * 2 * 100, nrow(df1) - 1), 0.9 * 2 * 100) #leaf area in m^2
df1$hv = df1$sap_area / df1$la

# Format dataframe: 40% liana canopy at BCI
df2 = cbind(dbhs, liana_out_r2[,2], rep('Liana', times = length(dbhs)))
colnames(df2) = c('DBH', 'Kreq', 'PFT')
df2 = rbind(df2, c(mean.hori.data.tree.dbh, tree_out_r2[2], 'Tree'))
df2 = as.data.frame(df2)
df2$DBH = as.numeric(as.character(df2$DBH))
df2$Kreq = as.numeric(as.character(df2$Kreq))
df2$Kreq = df2$Kreq * 0.001
df2$sap_area = rep(NA, nrow(df2))
for(i in 1:nrow(df2)){
  df2$sap_area[i] = min(((pi*df2$DBH[i]^2) / 4), (2.41 * (df2$DBH[i]/2)^1.97))
}
df2$la = c(rep(2 * 100 * 0.4, nrow(df2) - 1), 2 * 100 * 0.6)
df2$hv = df2$sap_area / df2$la

# Format dataframe: 10% liana canopy at Horizontes
df3 = cbind(dbhs, liana_out_r3[,2], rep('Liana', times = length(dbhs)))
colnames(df3) = c('DBH', 'Kreq', 'PFT')
df3 = rbind(df3, c(mean.hori.data.tree.dbh, tree_out_r3[2], 'Tree'))
df3 = as.data.frame(df3)
df3$DBH = as.numeric(as.character(df3$DBH))
df3$Kreq = as.numeric(as.character(df3$Kreq))
df3$Kreq = df3$Kreq * 0.001
df3$sap_area = rep(NA, nrow(df3))
for(i in 1:nrow(df3)){
  df3$sap_area[i] = min(((pi*df3$DBH[i]^2) / 4), (2.41 * (df3$DBH[i]/2)^1.97))
}
df3$la = c(rep(2 * 100 * 0.1, nrow(df3) - 1), 2 * 100 * 0.9)
df3$hv = df3$sap_area / df3$la

# Format dataframe: 40% liana canopy at Horizontes
df4 = cbind(dbhs, liana_out_r4[,2], rep('Liana', times = length(dbhs)))
colnames(df4) = c('DBH', 'Kreq', 'PFT')
df4 = rbind(df4, c(mean.hori.data.tree.dbh, tree_out_r4[2], 'Tree'))
df4 = as.data.frame(df4)
df4$DBH = as.numeric(as.character(df4$DBH))
df4$Kreq = as.numeric(as.character(df4$Kreq))
df4$Kreq = df4$Kreq * 0.001
df4$sap_area = rep(NA, nrow(df4))
for(i in 1:nrow(df4)){
  df4$sap_area[i] = min(((pi*df4$DBH[i]^2) / 4), (2.41 * (df4$DBH[i]/2)^1.97))
}
df4$la = c(rep(2 * 100 * 0.4, nrow(df4) - 1), 2 * 100 * 0.6)
df4$hv = df4$sap_area / df4$la

# Add ID to each dataframe
df1 = as.data.frame(cbind(df1, rep(1, nrow(df1))))
colnames(df1)[7] = 'df'
df2 = as.data.frame(cbind(df2, rep(2, nrow(df2))))
colnames(df2)[7] = 'df'
df3 = as.data.frame(cbind(df3, rep(3, nrow(df3))))
colnames(df3)[7] = 'df'
df4 = as.data.frame(cbind(df4, rep(4, nrow(df4))))
colnames(df4)[7] = 'df'
# Combine
df = as.data.frame(rbind(df1, df2, df3, df4))

# Plot DBH vs Kreq
p1 = ggplot(df, aes(x = DBH, y = log(Kreq), group = df)) +
  geom_line(data = subset(df, PFT == 'Liana'), aes(color = as.factor(df)), size = 1.2) +
  geom_hline(data = subset(df, PFT == 'Tree'), aes(yintercept = log(Kreq), color = as.factor(df)), size = 1, linetype = 'dashed') +
  xlab('DBH (cm)') + ylab(expression(paste('log(',K[req],')',' (mol/m/s/MPa)'))) +
  #scale_color_viridis_d(name = 'Scenario', end = 0.9, labels = c('1' = '\n10% liana\ntropical moist', '2' = '\n40% liana\ntropical moist\n', '3' = '\n10% liana\ntropical dry\n', '4' = '\n40% liana\ntropical dry\n'), breaks = c('4', '2', '3', '1')) +
  scale_color_manual(name = 'Scenario', 
                     labels = c('1' = '\n10% liana\ntropical moist', '2' = '\n40% liana\ntropical moist\n', '3' = '\n10% liana\ntropical dry\n', '4' = '\n40% liana\ntropical dry\n'), 
                     breaks = c('4', '2', '3', '1'),
                     values = c('4' = '#a6611a', '2' = '#80cdc1', '3' = '#dfc27d', '1' = '#018571')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.text = element_text(size = 26), legend.title = element_text(size = 28, hjust = 0.5), panel.grid = element_blank())
p1

# Plot Huber value vs Kreq
p2 = ggplot(df, aes(x = hv, y = log(Kreq), group = df)) +
  geom_line(data = subset(df, PFT == 'Liana'), aes(color = as.factor(df)), size = 1.5) +
  geom_point(data = subset(df, PFT == 'Tree'), aes(color = as.factor(df)), size = 8, shape = 4, stroke = 2, show.legend = F) +
  xlab(expression(paste('Huber value (c', m^2, '/', m^2, ')'))) + ylab(expression(paste('log(',K[req],')',' (mol/m/s/MPa)'))) +
  scale_color_manual(name = 'Scenario', 
                     labels = c('1' = '\n10% liana\ntropical moist', '2' = '\n40% liana\ntropical moist\n', '3' = '\n10% liana\ntropical dry\n', '4' = '\n40% liana\ntropical dry\n'), 
                     breaks = c('4', '2', '3', '1'),
                     values = c('4' = '#a6611a', '2' = '#80cdc1', '3' = '#dfc27d', '1' = '#018571')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.text = element_text(size = 26), legend.title = element_text(size = 28, hjust = 0.5), panel.grid = element_blank())
p2

# Function to extract legend from ggplot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Save legend to make formatting easier
leg = get_legend(p2)

# Plot together
ga = plot_grid(p1 + theme(legend.position = 'none'), 
               p2 + theme(legend.position = 'none'), 
               leg, rel_widths = c(2, 2, 1), nrow = 1)

ggsave(ga, filename = 'Plots/allom_hydro_comp_combined_log_200_dashed_hv_dbh.jpeg', width = 14, height = 8, units = 'in')

##################
#### Figure 3 ####
##################

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Interpolation_mxh_100.RData')

nsite = 100
nmonth = 12
nK = 75000

nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Interpolated VPD & SWP
VPDs = VPD_interp_100
SWPs = SWP_interp_100

################
## Tree plots ##
################
# Input values of K
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

tree.req = matrix(, nrow = nsite, ncol = nsite)

for(vsite in 1:nsite){
  for(ssite in 1:nsite){
    sens = c()
    for(k in 1:nK){
      for(mo in 1:nmonth){
        nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        month = mo
        
        tree.dbh = mean.hori.data.tree.dbh
        
        tot.al = 2 * 100
        frac.tree.al = 0.6
        
        psis = SWPs[mo,ssite]
        
        D = VPDs[mo,,vsite]
        
        if(mo == 1){
          tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
        }else{
          tree.height = tree.out[5]
        }
        K = Ks[k]
        tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
        sens[mo] = tree.out[1]
      }
      sum = sum(sens)
      if(sum > 0){
        tree.req[vsite,ssite] = K
        break
      }
    }
    print(paste('Passing sites',vsite,'&',ssite))
  }
}

# Save first time to be able to load later (long computation time)
#save(tree.req, file = 'tree.req.100interp.monthly_200.RData')

# Load once completed once and format
load(file = 'tree.req.100interp.monthly_200.RData')
tree.req = tree.req * 0.001
tree_req_melt = melt(tree.req)
colnames(tree_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

# Contour plot
ra1 = ggplot(tree_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv, fill = Ksurv)) +
  geom_raster(interpolate = T, show.legend = F) +
  geom_contour(alpha = 0, aes(color = ..level..), breaks = seq(0, 10, by = 1), show.legend = T) +
  geom_contour(alpha = 0.65, color = 'black', breaks = seq(0, 10, by = 1), show.legend = F) +
  xlab('VPD index') + ylab(expression(paste(Psi,' index'))) +
  scale_color_distiller(name = expression(K[req]), palette = 'PuBu', breaks = seq(0, 10, by = 2)) +
  scale_fill_distiller(name = expression(K[req]), palette = 'PuBu', breaks = seq(0, 10, by = 2)) +
  ggtitle('Tree') +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank()) +
  geom_dl(aes(label = ..level..), stat = 'contour', breaks = seq(0, 10, by = 1), method = list('top.pieces', color = 'black', cex = 1.8, vjust = -0.4))
ra1

#################
## Liana plots ##
#################

# Adjusted to improve speed
ks = c(trees_K$K / 10, lianas_K$K)
nK = 5000
Ks = seq(min(ks), max(ks), length.out = nK)
# This was determined by running with the tree K interpolation and keeping the lowest K required for the first index (wettest) and above
Ks = Ks[Ks > 38700]

liana.req = matrix(, nrow = nsite, ncol = nsite)

for(vsite in 1:nsite){
  for(ssite in 1:nsite){
    sens = c()
    for(k in 1:nK){
      for(mo in 1:12){
        nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        month = mo
        
        dbh = mean.data.liana.dbh
        
        tot.al = 2 * 100
        frac.liana.al = 0.4
        
        psis = SWPs[mo,ssite]
        
        D = VPDs[mo,,vsite]
        
        if(mo == 1){
          liana.length = tree.out[2]
        }else{
          liana.length = lianaout[4]
        }
        K = Ks[k]
        lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
        sens[mo] = lianaout[1]
      }
      sum = sum(sens)
      if(sum > 0){
        liana.req[vsite, ssite] = K
        break
      }
    }
    print(paste('Passing sites',vsite,'&',ssite))
  }
}

# Save and load as before
#save(liana.req, file = 'liana.req.100interp.monthly_200.RData')

load(file = 'liana.req.100interp.monthly_200.RData')
liana.req = liana.req * 0.001
liana_req_melt = melt(liana.req)
colnames(liana_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

# Contour plot
ra2 = ggplot(liana_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv, fill = Ksurv)) +
  geom_raster(interpolate = T, show.legend = T) +
  geom_contour(alpha = 0, aes(color = ..level..), show.legend = F, breaks = seq(30, 210, by = 20)) +
  geom_contour(alpha = 0.65, color = 'black', breaks = seq(30, 210, by = 20), show.legend = F) +
  xlab('VPD index') + ylab(expression(paste(Psi,' index'))) +
  scale_fill_distiller(name = expression(K[req]), palette = 'BuGn', breaks = seq(30, 210, by = 30)) +
  ggtitle('Liana') +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank()) +
  geom_dl(aes(label = ..level..), stat = 'contour', breaks = seq(30, 210, by = 20), method = list('top.pieces', color = 'black', cex = 1.8, vjust = -0.4))
ra2

# Descriptive plot
ra4 = ggplot(liana_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv)) +
  geom_blank() +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_segment(aes(x = 10, y = 7, xend = 100, yend = 7)) +
  geom_segment(aes(x = 10, y = 7, xend = 100, yend = 7), arrow = arrow(length = unit(0.25, 'cm'))) +
  geom_text(aes(x = 50, y = 11, label = 'Increasing VPD'), size = 10) +
  geom_segment(aes(x = 7, y = 10, xend = 7, yend = 100)) +
  geom_segment(aes(x = 7, y = 10, xend = 7, yend = 100), arrow = arrow(length = unit(0.25, 'cm'))) +
  geom_text(aes(x = 12, y = 50, label = paste('Decreasing ~~ Psi')), parse = T, angle = 90, size = 10) +
  geom_label(aes(x = 24, y = 92, label = 'Dry soil'), size = 10, fill = NA) +
  geom_label(aes(x = 80, y = 30, label = 'Dry \natmosphere'), size = 10, fill = NA) +
  geom_label(aes(x = 80, y = 92, label = 'Driest site'), size = 10, fill = NA) +
  geom_segment(aes(x = 8, y = 8, xend = 80, yend = 80)) +
  geom_segment(aes(x = 8, y = 8, xend = 80, yend = 80), arrow = arrow(length = unit(0.25, 'cm'))) +
  geom_text(aes(x = 46, y = 51, label = 'Increasing dryness'), angle = 45, size = 10) +
  theme(plot.margin = unit(c(0, 0, -2, 0), 'cm')) +
  xlab('VPD Index') + ylab(expression(paste(Psi,' Index')))

# Format plots
ra1_fin = ra1 + theme(legend.position = 'right', plot.margin = unit(c(0, 0, 0, 0), 'cm'))
ra2_fin = ra2 + theme(legend.position = 'right', plot.margin = unit(c(-1, 0, 1, 0), 'cm'))

# Plot together
ga = plot_grid(ra4, ra1_fin, ra2_fin, align = 'hv', axis = 'tbrl', nrow = 3)
ggsave(ga, filename = 'Plots/tree_liana_concept_final_rb_200.jpeg', height = 20, width = 8, units = 'in')

##################
#### Figure 4 ####
##################

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Driver_calc/Horizontes/horizontes_met_mxh.RData')
H_VPD = VPD
H_SWP = SWP
rm(VPD, SWP)

# Indexing months and days
nmonth = 12
nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

nK = 75000

# Horizontes future scenarios
H_VPD_future = array(, dim = c(12, 24, 2))
H_VPD_future[,,1] = H_VPD
H_VPD_future[,,2] = H_VPD * 2

#### Tree simulation - present day

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

tree.req_h = c()
sens = c()
for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 400
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 2 * 100
    frac.tree.al = 0.6
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,1]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h[1] = K
    break
  }
}

#### Liana simulation -- present day

liana.req_h = c()
sens = c()
for(k in 1:nK){
  for(mo in 1:12){
    Ca = 400
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = mean.data.liana.dbh
    
    tot.al = 2 * 100
    frac.liana.al = 0.4
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,1]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h[1] = K
    break
  }
}

#### Tree simulation -- future scenario
sens = c()
for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 2 * 100
    frac.tree.al = 0.6
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,2]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h[2] = K
    break
  }
}

#### Liana simulation -- future scenario

sens = c()
for(k in 1:nK){
  for(mo in 1:12){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = mean.data.liana.dbh
    
    tot.al = 2 * 100
    frac.liana.al = 0.4
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,2]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h[2] = K
    break
  }
}

tree.req_h = tree.req_h * 0.001
liana.req_h = liana.req_h * 0.001

# Make workable dataframes
liana_kreq = cbind(c('2000', '2100'), 
                   c(liana.req_h[1], liana.req_h[2]))
colnames(liana_kreq) = c('Time', 'Kreq')

tree_kreq = cbind(c('2000', '2100'),
                  c(tree.req_h[1], tree.req_h[2]))
colnames(tree_kreq) = c('Time', 'Kreq')

# Combine and format dataframes
comb_kreq = rbind(liana_kreq, tree_kreq)
comb_kreq = as.data.frame(comb_kreq)
comb_kreq$growth.form = c('Liana', 'Liana', 
                          'Tree', 'Tree')
comb_kreq$Time = as.factor(comb_kreq$Time)
comb_kreq$Kreq = as.numeric(comb_kreq$Kreq)
comb_kreq$growth.form = as.factor(comb_kreq$growth.form)

# Plot
ggplot(comb_kreq, aes(x = Time, y = Kreq, color = growth.form, group = growth.form)) +
  geom_point(size = 4, show.legend = F) +
  geom_line(size = 1.5, show.legend = F, alpha = 0.7) +
  scale_x_discrete(limits = c('2000', '2100')) +
  xlab('') + ylab(expression(paste(K[req], ' (mol/m/s/MPa)'))) +
  scale_color_npg(name = 'PFT') +
  annotate('segment', x = 0.95, xend = 2.05, y = -30, yend = -30, arrow = arrow(length = unit(0.4, 'cm')), size = 1.5) +
  annotate('text', x = '2000', hjust = -0.23, y = -40, label = 'Drying hyrdoclimate', size = 7) +
  annotate('text', x = '2000', hjust = -1, y = 140, label = 'Liana', size = 7, angle = 0, color = '#E64B35FF', fontface = 2) +
  annotate('text', x = '2000', hjust = -0.45, y = 128, label = expression(bold(paste(Delta, K[req], ' = 47'))), size = 6, angle = 0, color = '#E64B35FF') +
  annotate('text', x = '2000', hjust = -2.55, y = 24, label = 'Tree', size = 7, angle = 0, color = '#4DBBD5FF', fontface = 2) +
  annotate('text', x = '2000', hjust = -1.32, y = 12, label = expression(bold(paste(Delta, K[req], ' = 2'))), size = 6, angle = 0, color = '#4DBBD5FF') +
  coord_cartesian(clip = 'off', ylim = c(-2, 150), xlim = c(1.5, 1.5)) +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26)) +
  theme(plot.margin = unit(c(1, 1, 3, 1), 'lines'))

ggsave('Plots/Figure4.jpeg', width = 6, height = 6)
