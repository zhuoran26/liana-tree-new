## Main figures

## Author: AM Willson
## Date modified: 09 February 2021

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

df1 = cbind(dbhs, liana_out_r1[,2], rep('Liana', times = length(dbhs)))
colnames(df1) = c('DBH', 'Kreq', 'PFT')
df1 = rbind(df1, c(mean.hori.data.tree.dbh, tree_out_r1[2], 'Tree'))
df1 = as.data.frame(df1)
df1$DBH = as.numeric(as.character(df1$DBH))
df1$Kreq = as.numeric(as.character(df1$Kreq))
df1$Kreq = df1$Kreq * 0.001

df2 = cbind(dbhs, liana_out_r2[,2], rep('Liana', times = length(dbhs)))
colnames(df2) = c('DBH', 'Kreq', 'PFT')
df2 = rbind(df2, c(mean.hori.data.tree.dbh, tree_out_r2[2], 'Tree'))
df2 = as.data.frame(df2)
df2$DBH = as.numeric(as.character(df2$DBH))
df2$Kreq = as.numeric(as.character(df2$Kreq))
df2$Kreq = df2$Kreq * 0.001

df3 = cbind(dbhs, liana_out_r3[,2], rep('Liana', times = length(dbhs)))
colnames(df3) = c('DBH', 'Kreq', 'PFT')
df3 = rbind(df3, c(mean.hori.data.tree.dbh, tree_out_r3[2], 'Tree'))
df3 = as.data.frame(df3)
df3$DBH = as.numeric(as.character(df3$DBH))
df3$Kreq = as.numeric(as.character(df3$Kreq))
df3$Kreq = df3$Kreq * 0.001

df4 = cbind(dbhs, liana_out_r4[,2], rep('Liana', times = length(dbhs)))
colnames(df4) = c('DBH', 'Kreq', 'PFT')
df4 = rbind(df4, c(mean.hori.data.tree.dbh, tree_out_r4[2], 'Tree'))
df4 = as.data.frame(df4)
df4$DBH = as.numeric(as.character(df4$DBH))
df4$Kreq = as.numeric(as.character(df4$Kreq))
df4$Kreq = df4$Kreq * 0.001

df1 = as.data.frame(cbind(df1, rep(1, nrow(df1))))
colnames(df1)[4] = 'df'
df2 = as.data.frame(cbind(df2, rep(2, nrow(df2))))
colnames(df2)[4] = 'df'
df3 = as.data.frame(cbind(df3, rep(3, nrow(df3))))
colnames(df3)[4] = 'df'
df4 = as.data.frame(cbind(df4, rep(4, nrow(df4))))
colnames(df4)[4] = 'df'
df = as.data.frame(rbind(df1, df2, df3, df4))

pl_log = ggplot(df, aes(x = DBH, y = log(Kreq), group = df)) +
  geom_line(data = subset(df, PFT == 'Liana'), aes(color = as.factor(df)), size = 1.2) +
  geom_point(data = subset(df, PFT == 'Tree'), aes(color = as.factor(df)), size = 2) +
  xlab('DBH (cm)') + ylab(expression(paste('log(',K[req],')',' (mol/m/s/MPa)'))) +
  scale_color_viridis_d(name = 'Scenario', end = 0.9, labels = c('1' = '\n10% liana\ntropical moist', '2' = '\n40% liana\n tropical moist\n', '3' = '\n10% liana\ntropical dry\n', '4' = '\n40% liana\ntropical dry\n'), breaks = c('4', '2', '3', '1')) +
  facet_grid(~PFT, scales = 'free_x', space = 'free') +
  scale_x_continuous(expand = c(-0.04, 0.5), breaks = c(1:21)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.text = element_text(size = 26), legend.title = element_text(size = 28, hjust = 0.5), strip.text = element_text(size = 30))
pl_log

ggsave(pl_log, filename = 'Plots/allom_hydro_comp_combined_log_200.jpeg', width = 14, height = 8, units = 'in')

pl = ggplot(df, aes(x = DBH, y = log(Kreq), group = df)) +
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

ggsave(pl, filename = 'Plots/allom_hydro_comp_combined_log_200_dashed.jpeg', width = 14, height = 8, units = 'in')

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

VPDs = VPD_interp_100
SWPs = SWP_interp_100

################
## Tree plots ##
################
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

#save(tree.req, file = 'tree.req.100interp.monthly_200.RData')

load(file = 'tree.req.100interp.monthly_200.RData')
tree.req = tree.req * 0.001
tree_req_melt = melt(tree.req)
colnames(tree_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

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

#save(liana.req, file = 'liana.req.100interp.monthly_200.RData')

load(file = 'liana.req.100interp.monthly_200.RData')
liana.req = liana.req * 0.001
liana_req_melt = melt(liana.req)
colnames(liana_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

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

ra1_fin = ra1 + theme(legend.position = 'right', plot.margin = unit(c(0, 0, 0, 0), 'cm'))
ra2_fin = ra2 + theme(legend.position = 'right', plot.margin = unit(c(-1, 0, 1, 0), 'cm'))

ga = plot_grid(ra4, ra1_fin, ra2_fin, align = 'hv', axis = 'tbrl', nrow = 3)
ggsave(ga, filename = 'Plots/tree_liana_concept_final_rb_200.jpeg', height = 20, width = 8, units = 'in')

##################
#### Figure 4 ####
##################

rm(list = ls())

# Load parameters and models
load('param.input.RData')

# Load future Kreq for established scenario
# These were created for plots in the Supplement
load('tree_req_future_established_200_wCO2.RData')
load('liana_req_future_established_200_wCO2.RData')

# Convert Ks to mol/m/s/MPa
tree_K = trees_K$K * 0.001
liana_K = lianas_K$K * 0.001

# Convert Ks to Kw by dividing by 10
tree_K = tree_K / 10
liana_K = liana_K / 10

# Make histograms
pl1 = ggplot() +
  geom_histogram(aes(x = tree_K), bins = 20) +
  geom_vline(aes(xintercept = min(tree.req_bci), color = 'Wettest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(tree.req_bci), color = 'Wettest', linetype = 'Future'), size = 1.2) +
  geom_vline(aes(xintercept = min(tree.req_h), color = 'Driest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(tree.req_h), color = 'Driest', linetype = 'Future'), size = 1.2) +
  xlab('K (mol/m/s/MPa)') + ylab('Frequency') +
  scale_color_discrete(name = 'Hydroclimate') +
  scale_linetype_manual(name = 'Time Period', values = c('Present' = 'solid', 'Future' = 'dashed'), breaks = c('Present', 'Future')) +
  ggtitle('Tree') +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl1

pl2 = ggplot() +
  geom_histogram(aes(x = liana_K), bins = 20) +
  geom_vline(aes(xintercept = min(liana.req_bci), color = 'Wettest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(liana.req_bci), color = 'Wettest', linetype = 'Future'), size = 1.2) +
  geom_vline(aes(xintercept = min(liana.req_h), color = 'Driest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(liana.req_h), color = 'Driest', linetype = 'Future'), size = 1.2) +
  xlab('K (mol/m/s/MPa)') + ylab('Frequency') +
  scale_color_discrete(name = 'Hydroclimate') +
  scale_linetype_manual(name = 'Time Period', values = c('Present' = 'solid', 'Future' = 'dashed'), breaks = c('Present', 'Future')) +
  ggtitle('Liana') +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, 'cm'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl2

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend = get_legend(pl1)

ga = plot_grid(pl1 + theme(legend.position = 'none'), 
               pl2 + theme(legend.position = 'none'),
               legend, nrow = 1, rel_widths = c(1, 1, 0.3))
ga

ggsave(ga, filename = 'Plots/Kreq_Kw_comp_200_wCO2.jpeg', width = 20, height = 12, units = 'in')
