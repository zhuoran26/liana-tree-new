## Sensitivity analysis of liana & tree models
## This script walks through a sensitivity analysis, in which each input & parameter besides the
## hydroclimate inputs are perturbed one at a time to determine the sensitivity of Kreq of each growth
## form to changes in the parameters
## We compute Kreq across a range of input values +/- 50% of the default
## We run the sensitivity analysis 4 times:
## Wet hydroclimate + established
## Dry hydroclimate + established
## Wet hydroclimate + invasion
## Dry hydroclimate + invasion

## Note that the model must be recompiled multiple times in this script to allow static parameters
## from the original model to receive input values
## This contributes most of the length of the script

## Author: AM Willson
## Date modified: 20 March 2021

##########################
#### Set up workspace ####
##########################

rm(list = ls())

library(pracma)
library(dplyr)
library(gt)
library(ggplot2)
library(ggsci)
library(cowplot)

# Load functions and parameters
load('param.input.RData')
# Load hydroclimate data
load('bci_met_mxh.RData')
BCI_SWP = SWP
BCI_VPD = VPD
load('horizontes_met_mxh.RData')
H_SWP = SWP
H_VPD = VPD
rm(SWP, VPD)

# Conductivity loop
nK = 75000
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

# Create tables for outputs
tree_table = as_tibble(matrix(, nrow = 100, ncol = 5))
colnames(tree_table) = c('Parameter', 'Mean', 'Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity')
tree_table$Parameter = as.character(tree_table$Parameter)
tree_table$Mean = as.numeric(as.character(tree_table$Mean))
tree_table$`Hydroclimate Scenario` = as.character(tree_table$`Hydroclimate Scenario`)
tree_table$`Competition Scenario` = as.character(tree_table$`Competition Scenario`)
tree_table$Sensitivity = as.numeric(as.character(tree_table$Sensitivity))

liana_table = as_tibble(matrix(, nrow = 100, 5))
colnames(liana_table) = c('Parameter', 'Mean', 'Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity')
liana_table$Parameter = as.character(liana_table$Parameter)
liana_table$Mean = as.numeric(as.character(liana_table$Mean))
liana_table$`Hydroclimate Scenario` = as.character(liana_table$`Hydroclimate Scenario`)
liana_table$`Competition Scenario` = as.character(liana_table$`Competition Scenario`)
liana_table$Sensitivity = as.numeric(as.character(liana_table$Sensitivity))

# Storage to count the number of simulations
tree_count = 0
liana_count = 0

#############
#### DBH ####
#############

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = mean.hori.data.tree.dbh
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = altered[i]
      
      tot.al = 4 * 100
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
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
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[1] = 'DBH'
tree_table$Mean[1] = mu
tree_table$`Hydroclimate Scenario`[1] = 'BCI'
tree_table$`Competition Scenario`[1] = 'Established'
tree_table$Sensitivity[1] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = mean.hori.data.tree.dbh
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = altered[i]
      
      tot.al = 4 * 100
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
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
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[2] = 'DBH'
tree_table$Mean[2] = mu
tree_table$`Hydroclimate Scenario`[2] = 'Horizontes'
tree_table$`Competition Scenario`[2] = 'Established'
tree_table$Sensitivity[2] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Mean and 50% adjustment of parameter
mu = mean.data.liana.dbh
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = altered[i]
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
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
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[1] = 'DBH'
liana_table$Mean[1] = mu
liana_table$`Hydroclimate Scenario`[1] = 'BCI'
liana_table$`Competition Scenario`[1] = 'Established'
liana_table$Sensitivity[1] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = mean.data.liana.dbh
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = altered[i]
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
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
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[2] = 'DBH'
liana_table$Mean[2] = mu
liana_table$`Hydroclimate Scenario`[2] = 'Horizontes'
liana_table$`Competition Scenario`[2] = 'Established'
liana_table$Sensitivity[2] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = mean.hori.data.tree.dbh
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = altered[i]
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
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
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[3] = 'DBH'
tree_table$Mean[3] = mu
tree_table$`Hydroclimate Scenario`[3] = 'BCI'
tree_table$`Competition Scenario`[3] = 'Invasion'
tree_table$Sensitivity[3] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = mean.hori.data.tree.dbh
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = altered[i]
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
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
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[4] = 'DBH'
tree_table$Mean[4] = mu
tree_table$`Hydroclimate Scenario`[4] = 'Horizontes'
tree_table$`Competition Scenario`[4] = 'Invasion'
tree_table$Sensitivity[4] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = mean.data.liana.dbh
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = altered[i]
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
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
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[3] = 'DBH'
liana_table$Mean[3] = mu
liana_table$`Hydroclimate Scenario`[3] = 'BCI'
liana_table$`Competition Scenario`[3] = 'Invasion'
liana_table$Sensitivity[3] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = mean.data.liana.dbh
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = altered[i]
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
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
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[4] = 'DBH'
liana_table$Mean[4] = mu
liana_table$`Hydroclimate Scenario`[4] = 'Horizontes'
liana_table$`Competition Scenario`[4] = 'Invasion'
liana_table$Sensitivity[4] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

#############
#### SLA ####
#############

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 32
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[5] = 'SLA'
tree_table$Mean[5] = mu
tree_table$`Hydroclimate Scenario`[5] = 'BCI'
tree_table$`Competition Scenario`[5] = 'Established'
tree_table$Sensitivity[5] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 32
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[6] = 'SLA'
tree_table$Mean[6] = mu
tree_table$`Hydroclimate Scenario`[6] = 'Horizontes'
tree_table$`Competition Scenario`[6] = 'Established'
tree_table$Sensitivity[6] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Mean and 50% adjustment of parameter
mu = 32
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[5] = 'SLA'
liana_table$Mean[5] = mu
liana_table$`Hydroclimate Scenario`[5] = 'BCI'
liana_table$`Competition Scenario`[5] = 'Established'
liana_table$Sensitivity[5] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 32
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[6] = 'SLA'
liana_table$Mean[6] = mu
liana_table$`Hydroclimate Scenario`[6] = 'Horizontes'
liana_table$`Competition Scenario`[6] = 'Established'
liana_table$Sensitivity[6] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 32
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[7] = 'SLA'
tree_table$Mean[7] = mu
tree_table$`Hydroclimate Scenario`[7] = 'BCI'
tree_table$`Competition Scenario`[7] = 'Invasion'
tree_table$Sensitivity[7] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 32
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[8] = 'SLA'
tree_table$Mean[8] = mu
tree_table$`Hydroclimate Scenario`[8] = 'Horizontes'
tree_table$`Competition Scenario`[8] = 'Invasion'
tree_table$Sensitivity[8] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 32
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[7] = 'SLA'
liana_table$Mean[7] = mu
liana_table$`Hydroclimate Scenario`[7] = 'BCI'
liana_table$`Competition Scenario`[7] = 'Invasion'
liana_table$Sensitivity[7] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 32
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[8] = 'SLA'
liana_table$Mean[8] = mu
liana_table$`Hydroclimate Scenario`[8] = 'Horizontes'
liana_table$`Competition Scenario`[8] = 'Invasion'
liana_table$Sensitivity[8] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### b1 ####
############

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = b1.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[9] = 'b1'
tree_table$Mean[9] = mu
tree_table$`Hydroclimate Scenario`[9] = 'BCI'
tree_table$`Competition Scenario`[9] = 'Established'
tree_table$Sensitivity[9] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = b1.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[10] = 'b1'
tree_table$Mean[10] = mu
tree_table$`Hydroclimate Scenario`[10] = 'Horizontes'
tree_table$`Competition Scenario`[10] = 'Established'
tree_table$Sensitivity[10] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Mean and 50% adjustment of parameter
mu = b1.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[9] = 'b1'
liana_table$Mean[9] = mu
liana_table$`Hydroclimate Scenario`[9] = 'BCI'
liana_table$`Competition Scenario`[9] = 'Established'
liana_table$Sensitivity[9] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = b1.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[10] = 'b1'
liana_table$Mean[10] = mu
liana_table$`Hydroclimate Scenario`[10] = 'Horizontes'
liana_table$`Competition Scenario`[10] = 'Established'
liana_table$Sensitivity[10] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = b1.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[11] = 'b1'
tree_table$Mean[11] = mu
tree_table$`Hydroclimate Scenario`[11] = 'BCI'
tree_table$`Competition Scenario`[11] = 'Invasion'
tree_table$Sensitivity[11] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = b1.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[12] = 'b1'
tree_table$Mean[12] = mu
tree_table$`Hydroclimate Scenario`[12] = 'Horizontes'
tree_table$`Competition Scenario`[12] = 'Invasion'
tree_table$Sensitivity[12] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = b1.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[11] = 'b1'
liana_table$Mean[11] = mu
liana_table$`Hydroclimate Scenario`[11] = 'BCI'
liana_table$`Competition Scenario`[11] = 'Invasion'
liana_table$Sensitivity[11] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = b1.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[12] = 'b1'
liana_table$Mean[12] = mu
liana_table$`Hydroclimate Scenario`[12] = 'Horizontes'
liana_table$`Competition Scenario`[12] = 'Invasion'
liana_table$Sensitivity[12] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### b2 ####
############

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = b2.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[13] = 'b2'
tree_table$Mean[13] = mu
tree_table$`Hydroclimate Scenario`[13] = 'BCI'
tree_table$`Competition Scenario`[13] = 'Established'
tree_table$Sensitivity[13] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = b2.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[14] = 'b2'
tree_table$Mean[14] = mu
tree_table$`Hydroclimate Scenario`[14] = 'Horizontes'
tree_table$`Competition Scenario`[14] = 'Established'
tree_table$Sensitivity[14] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Mean and 50% adjustment of parameter
mu = b2.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[13] = 'b2'
liana_table$Mean[13] = mu
liana_table$`Hydroclimate Scenario`[13] = 'BCI'
liana_table$`Competition Scenario`[13] = 'Established'
liana_table$Sensitivity[13] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = b2.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[14] = 'b2'
liana_table$Mean[14] = mu
liana_table$`Hydroclimate Scenario`[14] = 'Horizontes'
liana_table$`Competition Scenario`[14] = 'Established'
liana_table$Sensitivity[14] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = b2.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[15] = 'b2'
tree_table$Mean[15] = mu
tree_table$`Hydroclimate Scenario`[15] = 'BCI'
tree_table$`Competition Scenario`[15] = 'Invasion'
tree_table$Sensitivity[15] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = b2.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[16] = 'b2'
tree_table$Mean[16] = mu
tree_table$`Hydroclimate Scenario`[16] = 'Horizontes'
tree_table$`Competition Scenario`[16] = 'Invasion'
tree_table$Sensitivity[16] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = b2.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[15] = 'b2'
liana_table$Mean[15] = mu
liana_table$`Hydroclimate Scenario`[15] = 'BCI'
liana_table$`Competition Scenario`[15] = 'Invasion'
liana_table$Sensitivity[15] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = b2.mean.all
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[16] = 'b2'
liana_table$Mean[16] = mu
liana_table$`Hydroclimate Scenario`[16] = 'Horizontes'
liana_table$`Competition Scenario`[16] = 'Invasion'
liana_table$Sensitivity[16] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### AL ####
############

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 200
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = altered[i]
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[17] = 'AL'
tree_table$Mean[17] = mu
tree_table$`Hydroclimate Scenario`[17] = 'BCI'
tree_table$`Competition Scenario`[17] = 'Established'
tree_table$Sensitivity[17] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 200
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = altered[i]
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[18] = 'AL'
tree_table$Mean[18] = mu
tree_table$`Hydroclimate Scenario`[18] = 'Horizontes'
tree_table$`Competition Scenario`[18] = 'Established'
tree_table$Sensitivity[18] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Mean and 50% adjustment of parameter
mu = 200
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = altered[i]
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[17] = 'AL'
liana_table$Mean[17] = mu
liana_table$`Hydroclimate Scenario`[17] = 'BCI'
liana_table$`Competition Scenario`[17] = 'Established'
liana_table$Sensitivity[17] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 200
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = altered[i]
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[18] = 'AL'
liana_table$Mean[18] = mu
liana_table$`Hydroclimate Scenario`[18] = 'Horizontes'
liana_table$`Competition Scenario`[18] = 'Established'
liana_table$Sensitivity[18] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 200
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = altered[i]
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[19] = 'AL'
tree_table$Mean[19] = mu
tree_table$`Hydroclimate Scenario`[19] = 'BCI'
tree_table$`Competition Scenario`[19] = 'Invasion'
tree_table$Sensitivity[19] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 200
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = altered[i]
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[20] = 'AL'
tree_table$Mean[20] = mu
tree_table$`Hydroclimate Scenario`[20] = 'Horizontes'
tree_table$`Competition Scenario`[20] = 'Invasion'
tree_table$Sensitivity[20] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 200
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = altered[i]
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[19] = 'AL'
liana_table$Mean[19] = mu
liana_table$`Hydroclimate Scenario`[19] = 'BCI'
liana_table$`Competition Scenario`[19] = 'Invasion'
liana_table$Sensitivity[19] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 200
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = altered[i]
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[20] = 'AL'
liana_table$Mean[20] = mu
liana_table$`Hydroclimate Scenario`[20] = 'Horizontes'
liana_table$`Competition Scenario`[20] = 'Invasion'
liana_table$Sensitivity[20] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### Lx ####
############

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = tree.b1Ht * mean.hori.data.tree.dbh^tree.b2Ht
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = altered[i]
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[21] = 'Lx'
tree_table$Mean[21] = mu
tree_table$`Hydroclimate Scenario`[21] = 'BCI'
tree_table$`Competition Scenario`[21] = 'Established'
tree_table$Sensitivity[21] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = tree.b1Ht * mean.hori.data.tree.dbh^tree.b2Ht
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = altered[i]
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[22] = 'Lx'
tree_table$Mean[22] = mu
tree_table$`Hydroclimate Scenario`[22] = 'Horizontes'
tree_table$`Competition Scenario`[22] = 'Established'
tree_table$Sensitivity[22] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Mean and 50% adjustment of parameter
mu = tree.b1Ht * mean.hori.data.tree.dbh^tree.b2Ht
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = altered[i]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[21] = 'Lx'
liana_table$Mean[21] = mu
liana_table$`Hydroclimate Scenario`[21] = 'BCI'
liana_table$`Competition Scenario`[21] = 'Established'
liana_table$Sensitivity[21] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = tree.b1Ht * mean.hori.data.tree.dbh^tree.b2Ht
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = altered[i]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[22] = 'Lx'
liana_table$Mean[22] = mu
liana_table$`Hydroclimate Scenario`[22] = 'Horizontes'
liana_table$`Competition Scenario`[22] = 'Established'
liana_table$Sensitivity[22] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = tree.b1Ht * mean.hori.data.tree.dbh^tree.b2Ht
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = altered[i]
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[23] = 'Lx'
tree_table$Mean[23] = mu
tree_table$`Hydroclimate Scenario`[23] = 'BCI'
tree_table$`Competition Scenario`[23] = 'Invasion'
tree_table$Sensitivity[23] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = tree.b1Ht * mean.hori.data.tree.dbh^tree.b2Ht
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = altered[i]
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[24] = 'Lx'
tree_table$Mean[24] = mu
tree_table$`Hydroclimate Scenario`[24] = 'Horizontes'
tree_table$`Competition Scenario`[24] = 'Invasion'
tree_table$Sensitivity[24] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = tree.b1Ht * mean.hori.data.tree.dbh^tree.b2Ht
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = altered[i]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[23] = 'Lx'
liana_table$Mean[23] = mu
liana_table$`Hydroclimate Scenario`[23] = 'BCI'
liana_table$`Competition Scenario`[23] = 'Invasion'
liana_table$Sensitivity[23] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = tree.b1Ht * mean.hori.data.tree.dbh^tree.b2Ht
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = altered[i]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[24] = 'Lx'
liana_table$Mean[24] = mu
liana_table$`Hydroclimate Scenario`[24] = 'Horizontes'
liana_table$`Competition Scenario`[24] = 'Invasion'
liana_table$Sensitivity[24] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### Vm ####
############

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 50
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = altered[i], SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[25] = 'Vm'
tree_table$Mean[25] = mu
tree_table$`Hydroclimate Scenario`[25] = 'BCI'
tree_table$`Competition Scenario`[25] = 'Established'
tree_table$Sensitivity[25] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 50
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = altered[i], b1 = NULL, b2 = NULL)
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[26] = 'Vm'
tree_table$Mean[26] = mu
tree_table$`Hydroclimate Scenario`[26] = 'Horizontes'
tree_table$`Competition Scenario`[26] = 'Established'
tree_table$Sensitivity[26] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Mean and 50% adjustment of parameter
mu = 50
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = altered[i], SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[25] = 'Vm'
liana_table$Mean[25] = mu
liana_table$`Hydroclimate Scenario`[25] = 'BCI'
liana_table$`Competition Scenario`[25] = 'Established'
liana_table$Sensitivity[25] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 50
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = altered[i], SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[26] = 'Vm'
liana_table$Mean[26] = mu
liana_table$`Hydroclimate Scenario`[26] = 'Horizontes'
liana_table$`Competition Scenario`[26] = 'Established'
liana_table$Sensitivity[26] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 50
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = altered[i], SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[27] = 'Vm'
tree_table$Mean[27] = mu
tree_table$`Hydroclimate Scenario`[27] = 'BCI'
tree_table$`Competition Scenario`[27] = 'Invasion'
tree_table$Sensitivity[27] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 50
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = altered[i], SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[28] = 'Vm'
tree_table$Mean[28] = mu
tree_table$`Hydroclimate Scenario`[28] = 'Horizontes'
tree_table$`Competition Scenario`[28] = 'Invasion'
tree_table$Sensitivity[28] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 50
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = altered[i], SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[27] = 'Vm'
liana_table$Mean[27] = mu
liana_table$`Hydroclimate Scenario`[27] = 'BCI'
liana_table$`Competition Scenario`[27] = 'Invasion'
liana_table$Sensitivity[27] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 50
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = altered[i], SLA = NULL, b1 = NULL, b2 = NULL)
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[28] = 'Vm'
liana_table$Mean[28] = mu
liana_table$`Hydroclimate Scenario`[28] = 'Horizontes'
liana_table$`Competition Scenario`[28] = 'Invasion'
liana_table$Sensitivity[28] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

##################
#### sap.frac ####
##################

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, sap.frac){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = sap.frac # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8*10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 0.622
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, sap.frac = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[29] = 'sap.frac'
tree_table$Mean[29] = mu
tree_table$`Hydroclimate Scenario`[29] = 'BCI'
tree_table$`Competition Scenario`[29] = 'Established'
tree_table$Sensitivity[29] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.622
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, sap.frac = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[30] = 'sap.frac'
tree_table$Mean[30] = mu
tree_table$`Hydroclimate Scenario`[30] = 'Horizontes'
tree_table$`Competition Scenario`[30] = 'Established'
tree_table$Sensitivity[30] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, sap.frac){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = sap.frac # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8*10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 0.622
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, sap.frac = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[29] = 'sap.frac'
liana_table$Mean[29] = mu
liana_table$`Hydroclimate Scenario`[29] = 'BCI'
liana_table$`Competition Scenario`[29] = 'Established'
liana_table$Sensitivity[29] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.622
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, sap.frac = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[30] = 'sap.frac'
liana_table$Mean[30] = mu
liana_table$`Hydroclimate Scenario`[30] = 'Horizontes'
liana_table$`Competition Scenario`[30] = 'Established'
liana_table$Sensitivity[30] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 0.622
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, sap.frac = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[31] = 'sap.frac'
tree_table$Mean[31] = mu
tree_table$`Hydroclimate Scenario`[31] = 'BCI'
tree_table$`Competition Scenario`[31] = 'Invasion'
tree_table$Sensitivity[31] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.622
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, sap.frac = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[32] = 'sap.frac'
tree_table$Mean[32] = mu
tree_table$`Hydroclimate Scenario`[32] = 'Horizontes'
tree_table$`Competition Scenario`[32] = 'Invasion'
tree_table$Sensitivity[32] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 0.622
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, sap.frac = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[31] = 'sap.frac'
liana_table$Mean[31] = mu
liana_table$`Hydroclimate Scenario`[31] = 'BCI'
liana_table$`Competition Scenario`[31] = 'Invasion'
liana_table$Sensitivity[31] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.622
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, sap.frac = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[32] = 'sap.frac'
liana_table$Mean[32] = mu
liana_table$`Hydroclimate Scenario`[32] = 'Horizontes'
liana_table$`Competition Scenario`[32] = 'Invasion'
liana_table$Sensitivity[32] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### Ca ####
############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, Ca){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = Ca # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8*10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 400
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Ca = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[33] = 'Ca'
tree_table$Mean[33] = mu
tree_table$`Hydroclimate Scenario`[33] = 'BCI'
tree_table$`Competition Scenario`[33] = 'Established'
tree_table$Sensitivity[33] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 400
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, Ca = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[34] = 'Ca'
tree_table$Mean[34] = mu
tree_table$`Hydroclimate Scenario`[34] = 'Horizontes'
tree_table$`Competition Scenario`[34] = 'Established'
tree_table$Sensitivity[34] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, Ca){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = Ca # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8*10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 400
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Ca = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[33] = 'Ca'
liana_table$Mean[33] = mu
liana_table$`Hydroclimate Scenario`[33] = 'BCI'
liana_table$`Competition Scenario`[33] = 'Established'
liana_table$Sensitivity[33] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 400
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Ca = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[34] = 'Ca'
liana_table$Mean[34] = mu
liana_table$`Hydroclimate Scenario`[34] = 'Horizontes'
liana_table$`Competition Scenario`[34] = 'Established'
liana_table$Sensitivity[34] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 400
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Ca = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[35] = 'Ca'
tree_table$Mean[35] = mu
tree_table$`Hydroclimate Scenario`[35] = 'BCI'
tree_table$`Competition Scenario`[35] = 'Invasion'
tree_table$Sensitivity[35] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 400
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Ca = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[36] = 'Ca'
tree_table$Mean[36] = mu
tree_table$`Hydroclimate Scenario`[36] = 'Horizontes'
tree_table$`Competition Scenario`[36] = 'Invasion'
tree_table$Sensitivity[36] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 400
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Ca = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[35] = 'Ca'
liana_table$Mean[35] = mu
liana_table$`Hydroclimate Scenario`[35] = 'BCI'
liana_table$`Competition Scenario`[35] = 'Invasion'
liana_table$Sensitivity[35] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 400
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Ca = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[36] = 'Ca'
liana_table$Mean[36] = mu
liana_table$`Hydroclimate Scenario`[36] = 'Horizontes'
liana_table$`Competition Scenario`[36] = 'Invasion'
liana_table$Sensitivity[36] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### rG ####
############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, rG){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = rG # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8*10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 0.3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rG = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[37] = 'rG'
tree_table$Mean[37] = mu
tree_table$`Hydroclimate Scenario`[37] = 'BCI'
tree_table$`Competition Scenario`[37] = 'Established'
tree_table$Sensitivity[37] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, rG = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[38] = 'rG'
tree_table$Mean[38] = mu
tree_table$`Hydroclimate Scenario`[38] = 'Horizontes'
tree_table$`Competition Scenario`[38] = 'Established'
tree_table$Sensitivity[38] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, rG){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = rG # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8*10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 0.3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rG = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[37] = 'rG'
liana_table$Mean[37] = mu
liana_table$`Hydroclimate Scenario`[37] = 'BCI'
liana_table$`Competition Scenario`[37] = 'Established'
liana_table$Sensitivity[37] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rG = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[38] = 'rG'
liana_table$Mean[38] = mu
liana_table$`Hydroclimate Scenario`[38] = 'Horizontes'
liana_table$`Competition Scenario`[38] = 'Established'
liana_table$Sensitivity[38] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 0.3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rG = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[39] = 'rG'
tree_table$Mean[39] = mu
tree_table$`Hydroclimate Scenario`[39] = 'BCI'
tree_table$`Competition Scenario`[39] = 'Invasion'
tree_table$Sensitivity[39] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rG = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[40] = 'rG'
tree_table$Mean[40] = mu
tree_table$`Hydroclimate Scenario`[40] = 'Horizontes'
tree_table$`Competition Scenario`[40] = 'Invasion'
tree_table$Sensitivity[40] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 0.3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rG = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[39] = 'rG'
liana_table$Mean[39] = mu
liana_table$`Hydroclimate Scenario`[39] = 'BCI'
liana_table$`Competition Scenario`[39] = 'Invasion'
liana_table$Sensitivity[39] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rG = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[40] = 'rG'
liana_table$Mean[40] = mu
liana_table$`Hydroclimate Scenario`[40] = 'Horizontes'
liana_table$`Competition Scenario`[40] = 'Invasion'
liana_table$Sensitivity[40] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

###########
#### q ####
###########

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, q){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = q # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8*10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 1.89
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, q = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

# Add to table
tree_table$Parameter[41] = 'q'
tree_table$Mean[41] = mu
tree_table$`Hydroclimate Scenario`[41] = 'BCI'
tree_table$`Competition Scenario`[41] = 'Established'
tree_table$Sensitivity[41] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 1.89
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, q = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[42] = 'q'
tree_table$Mean[42] = mu
tree_table$`Hydroclimate Scenario`[42] = 'Horizontes'
tree_table$`Competition Scenario`[42] = 'Established'
tree_table$Sensitivity[42] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, q){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = q # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8*10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 1.89
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, q = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[41] = 'q'
liana_table$Mean[41] = mu
liana_table$`Hydroclimate Scenario`[41] = 'BCI'
liana_table$`Competition Scenario`[41] = 'Established'
liana_table$Sensitivity[41] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 1.89
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, q = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[42] = 'q'
liana_table$Mean[42] = mu
liana_table$`Hydroclimate Scenario`[42] = 'Horizontes'
liana_table$`Competition Scenario`[42] = 'Established'
liana_table$Sensitivity[42] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 1.89
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, q = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[43] = 'q'
tree_table$Mean[43] = mu
tree_table$`Hydroclimate Scenario`[43] = 'BCI'
tree_table$`Competition Scenario`[43] = 'Invasion'
tree_table$Sensitivity[43] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 1.89
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, q = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[44] = 'q'
tree_table$Mean[44] = mu
tree_table$`Hydroclimate Scenario`[44] = 'Horizontes'
tree_table$`Competition Scenario`[44] = 'Invasion'
tree_table$Sensitivity[44] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 1.89
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, q = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[43] = 'q'
liana_table$Mean[43] = mu
liana_table$`Hydroclimate Scenario`[43] = 'BCI'
liana_table$`Competition Scenario`[43] = 'Invasion'
liana_table$Sensitivity[43] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 1.89
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, q = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[44] = 'q'
liana_table$Mean[44] = mu
liana_table$`Hydroclimate Scenario`[44] = 'Horizontes'
liana_table$`Competition Scenario`[44] = 'Invasion'
liana_table$Sensitivity[44] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

#############
#### rho ####
#############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, rho){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = rho # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8*10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 420
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rho = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[45] = 'rho'
tree_table$Mean[45] = mu
tree_table$`Hydroclimate Scenario`[45] = 'BCI'
tree_table$`Competition Scenario`[45] = 'Established'
tree_table$Sensitivity[45] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 420
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, rho = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[46] = 'rho'
tree_table$Mean[46] = mu
tree_table$`Hydroclimate Scenario`[46] = 'Horizontes'
tree_table$`Competition Scenario`[46] = 'Established'
tree_table$Sensitivity[46] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, rho){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = rho # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8*10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 420
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rho = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[45] = 'rho'
liana_table$Mean[45] = mu
liana_table$`Hydroclimate Scenario`[45] = 'BCI'
liana_table$`Competition Scenario`[45] = 'Established'
liana_table$Sensitivity[45] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 420
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rho = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[46] = 'rho'
liana_table$Mean[46] = mu
liana_table$`Hydroclimate Scenario`[46] = 'Horizontes'
liana_table$`Competition Scenario`[46] = 'Established'
liana_table$Sensitivity[46] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 420
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rho = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[47] = 'rho'
tree_table$Mean[47] = mu
tree_table$`Hydroclimate Scenario`[47] = 'BCI'
tree_table$`Competition Scenario`[47] = 'Invasion'
tree_table$Sensitivity[47] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 420
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rho = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[48] = 'rho'
tree_table$Mean[48] = mu
tree_table$`Hydroclimate Scenario`[48] = 'Horizontes'
tree_table$`Competition Scenario`[48] = 'Invasion'
tree_table$Sensitivity[48] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 420
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rho = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[47] = 'rho'
liana_table$Mean[47] = mu
liana_table$`Hydroclimate Scenario`[47] = 'BCI'
liana_table$`Competition Scenario`[47] = 'Invasion'
liana_table$Sensitivity[47] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 420
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rho = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[48] = 'rho'
liana_table$Mean[48] = mu
liana_table$`Hydroclimate Scenario`[48] = 'Horizontes'
liana_table$`Competition Scenario`[48] = 'Invasion'
liana_table$Sensitivity[48] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

#############
#### SRA ####
#############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, SRA){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = SRA # Specific root area (m^2 / kg C), default
    
    Lr = 1.8*10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 80
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, SRA = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[49] = 'SRA'
tree_table$Mean[49] = mu
tree_table$`Hydroclimate Scenario`[49] = 'BCI'
tree_table$`Competition Scenario`[49] = 'Established'
tree_table$Sensitivity[49] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 80
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, SRA = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[50] = 'SRA'
tree_table$Mean[50] = mu
tree_table$`Hydroclimate Scenario`[50] = 'Horizontes'
tree_table$`Competition Scenario`[50] = 'Established'
tree_table$Sensitivity[50] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, SRA){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = SRA # Specific root area (m^2 / kgC), default
    
    Lr = 1.8*10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 80
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, SRA = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[49] = 'SRA'
liana_table$Mean[49] = mu
liana_table$`Hydroclimate Scenario`[49] = 'BCI'
liana_table$`Competition Scenario`[49] = 'Established'
liana_table$Sensitivity[49] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 80
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, SRA = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[50] = 'SRA'
liana_table$Mean[50] = mu
liana_table$`Hydroclimate Scenario`[50] = 'Horizontes'
liana_table$`Competition Scenario`[50] = 'Established'
liana_table$Sensitivity[50] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 80
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, SRA = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[51] = 'SRA'
tree_table$Mean[51] = mu
tree_table$`Hydroclimate Scenario`[51] = 'BCI'
tree_table$`Competition Scenario`[51] = 'Invasion'
tree_table$Sensitivity[51] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 80
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, SRA = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[52] = 'SRA'
tree_table$Mean[52] = mu
tree_table$`Hydroclimate Scenario`[52] = 'Horizontes'
tree_table$`Competition Scenario`[52] = 'Invasion'
tree_table$Sensitivity[52] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 80
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, SRA = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[51] = 'SRA'
liana_table$Mean[51] = mu
liana_table$`Hydroclimate Scenario`[51] = 'BCI'
liana_table$`Competition Scenario`[51] = 'Invasion'
liana_table$Sensitivity[51] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 80
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, SRA = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[52] = 'SRA'
liana_table$Mean[52] = mu
liana_table$`Hydroclimate Scenario`[52] = 'Horizontes'
liana_table$`Competition Scenario`[52] = 'Invasion'
liana_table$Sensitivity[52] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### Lr ####
############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, Lr){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = Lr # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 1.8 * 10^4
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lr = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[53] = 'Lr'
tree_table$Mean[53] = mu
tree_table$`Hydroclimate Scenario`[53] = 'BCI'
tree_table$`Competition Scenario`[53] = 'Established'
tree_table$Sensitivity[53] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 1.8 * 10^4
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, Lr = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[54] = 'Lr'
tree_table$Mean[54] = mu
tree_table$`Hydroclimate Scenario`[54] = 'Horizontes'
tree_table$`Competition Scenario`[54] = 'Established'
tree_table$Sensitivity[54] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, Lr){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = Lr # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 1.8 * 10^4
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lr = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[53] = 'Lr'
liana_table$Mean[53] = mu
liana_table$`Hydroclimate Scenario`[53] = 'BCI'
liana_table$`Competition Scenario`[53] = 'Established'
liana_table$Sensitivity[53] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 1.8 * 10^4
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lr = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[54] = 'Lr'
liana_table$Mean[54] = mu
liana_table$`Hydroclimate Scenario`[54] = 'Horizontes'
liana_table$`Competition Scenario`[54] = 'Established'
liana_table$Sensitivity[54] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 1.8 * 10^4
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lr = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[55] = 'Lr'
tree_table$Mean[55] = mu
tree_table$`Hydroclimate Scenario`[55] = 'BCI'
tree_table$`Competition Scenario`[55] = 'Invasion'
tree_table$Sensitivity[55] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 1.8 * 10^4
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lr = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[56] = 'Lr'
tree_table$Mean[56] = mu
tree_table$`Hydroclimate Scenario`[56] = 'Horizontes'
tree_table$`Competition Scenario`[56] = 'Invasion'
tree_table$Sensitivity[56] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 1.8 * 10^4
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lr = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[55] = 'Lr'
liana_table$Mean[55] = mu
liana_table$`Hydroclimate Scenario`[55] = 'BCI'
liana_table$`Competition Scenario`[55] = 'Invasion'
liana_table$Sensitivity[55] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 1.8 * 10^4
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lr = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[56] = 'Lr'
liana_table$Mean[56] = mu
liana_table$`Hydroclimate Scenario`[56] = 'Horizontes'
liana_table$`Competition Scenario`[56] = 'Invasion'
liana_table$Sensitivity[56] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### Lp ####
############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, Lp){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8 * 10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = Lp # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 0.1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lp = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[57] = 'Lp'
tree_table$Mean[57] = mu
tree_table$`Hydroclimate Scenario`[57] = 'BCI'
tree_table$`Competition Scenario`[57] = 'Established'
tree_table$Sensitivity[57] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, Lp = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[58] = 'Lp'
tree_table$Mean[58] = mu
tree_table$`Hydroclimate Scenario`[58] = 'Horizontes'
tree_table$`Competition Scenario`[58] = 'Established'
tree_table$Sensitivity[58] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, Lp){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8 * 10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = Lp # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 0.1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lp = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[57] = 'Lp'
liana_table$Mean[57] = mu
liana_table$`Hydroclimate Scenario`[57] = 'BCI'
liana_table$`Competition Scenario`[57] = 'Established'
liana_table$Sensitivity[57] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lp = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[58] = 'Lp'
liana_table$Mean[58] = mu
liana_table$`Hydroclimate Scenario`[58] = 'Horizontes'
liana_table$`Competition Scenario`[58] = 'Established'
liana_table$Sensitivity[58] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 0.1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lp = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[59] = 'Lp'
tree_table$Mean[59] = mu
tree_table$`Hydroclimate Scenario`[59] = 'BCI'
tree_table$`Competition Scenario`[59] = 'Invasion'
tree_table$Sensitivity[59] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lp = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[60] = 'Lp'
tree_table$Mean[60] = mu
tree_table$`Hydroclimate Scenario`[60] = 'Horizontes'
tree_table$`Competition Scenario`[60] = 'Invasion'
tree_table$Sensitivity[60] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 0.1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lp = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[59] = 'Lp'
liana_table$Mean[59] = mu
liana_table$`Hydroclimate Scenario`[59] = 'BCI'
liana_table$`Competition Scenario`[59] = 'Invasion'
liana_table$Sensitivity[59] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, Lp = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[60] = 'Lp'
liana_table$Mean[60] = mu
liana_table$`Hydroclimate Scenario`[60] = 'Horizontes'
liana_table$`Competition Scenario`[60] = 'Invasion'
liana_table$Sensitivity[60] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

###########
#### r ####
###########

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, r){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8 * 10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = r # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 0.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, r = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[61] = 'r'
tree_table$Mean[61] = mu
tree_table$`Hydroclimate Scenario`[61] = 'BCI'
tree_table$`Competition Scenario`[61] = 'Established'
tree_table$Sensitivity[61] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, r = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[62] = 'r'
tree_table$Mean[62] = mu
tree_table$`Hydroclimate Scenario`[62] = 'Horizontes'
tree_table$`Competition Scenario`[62] = 'Established'
tree_table$Sensitivity[62] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, r){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8 * 10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = r # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 0.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, r = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[61] = 'r'
liana_table$Mean[61] = mu
liana_table$`Hydroclimate Scenario`[61] = 'BCI'
liana_table$`Competition Scenario`[61] = 'Established'
liana_table$Sensitivity[61] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, r = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[62] = 'r'
liana_table$Mean[62] = mu
liana_table$`Hydroclimate Scenario`[62] = 'Horizontes'
liana_table$`Competition Scenario`[62] = 'Established'
liana_table$Sensitivity[62] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 0.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, r = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[63] = 'r'
tree_table$Mean[63] = mu
tree_table$`Hydroclimate Scenario`[63] = 'BCI'
tree_table$`Competition Scenario`[63] = 'Invasion'
tree_table$Sensitivity[63] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, r = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[64] = 'r'
tree_table$Mean[64] = mu
tree_table$`Hydroclimate Scenario`[64] = 'Horizontes'
tree_table$`Competition Scenario`[64] = 'Invasion'
tree_table$Sensitivity[64] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 0.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, r = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[63] = 'r'
liana_table$Mean[63] = mu
liana_table$`Hydroclimate Scenario`[63] = 'BCI'
liana_table$`Competition Scenario`[63] = 'Invasion'
liana_table$Sensitivity[63] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, r = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[64] = 'r'
liana_table$Mean[64] = mu
liana_table$`Hydroclimate Scenario`[64] = 'Horizontes'
liana_table$`Competition Scenario`[64] = 'Invasion'
liana_table$Sensitivity[64] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

####################
#### gamma.star ####
####################

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, gamma.star){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8 * 10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = gamma.star # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma.star = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[65] = 'gamma.star'
tree_table$Mean[65] = mu
tree_table$`Hydroclimate Scenario`[65] = 'BCI'
tree_table$`Competition Scenario`[65] = 'Established'
tree_table$Sensitivity[65] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, gamma.star = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[66] = 'gamma.star'
tree_table$Mean[66] = mu
tree_table$`Hydroclimate Scenario`[66] = 'Horizontes'
tree_table$`Competition Scenario`[66] = 'Established'
tree_table$Sensitivity[66] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, gamma.star){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8 * 10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = gamma.star # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma.star = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[65] = 'gamma.star'
liana_table$Mean[65] = mu
liana_table$`Hydroclimate Scenario`[65] = 'BCI'
liana_table$`Competition Scenario`[65] = 'Established'
liana_table$Sensitivity[65] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma.star = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[66] = 'gamma.star'
liana_table$Mean[66] = mu
liana_table$`Hydroclimate Scenario`[66] = 'Horizontes'
liana_table$`Competition Scenario`[66] = 'Established'
liana_table$Sensitivity[66] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma.star = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[67] = 'gamma.star'
tree_table$Mean[67] = mu
tree_table$`Hydroclimate Scenario`[67] = 'BCI'
tree_table$`Competition Scenario`[67] = 'Invasion'
tree_table$Sensitivity[67] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma.star = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[68] = 'gamma.star'
tree_table$Mean[68] = mu
tree_table$`Hydroclimate Scenario`[68] = 'Horizontes'
tree_table$`Competition Scenario`[68] = 'Invasion'
tree_table$Sensitivity[68] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma.star = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[67] = 'gamma.star'
liana_table$Mean[67] = mu
liana_table$`Hydroclimate Scenario`[67] = 'BCI'
liana_table$`Competition Scenario`[67] = 'Invasion'
liana_table$Sensitivity[67] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma.star = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[68] = 'gamma.star'
liana_table$Mean[68] = mu
liana_table$`Hydroclimate Scenario`[68] = 'Horizontes'
liana_table$`Competition Scenario`[68] = 'Invasion'
liana_table$Sensitivity[68] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### km ####
############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, km){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8 * 10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = km #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 300
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, km = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[69] = 'km'
tree_table$Mean[69] = mu
tree_table$`Hydroclimate Scenario`[69] = 'BCI'
tree_table$`Competition Scenario`[69] = 'Established'
tree_table$Sensitivity[69] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 300
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, km = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[70] = 'km'
tree_table$Mean[70] = mu
tree_table$`Hydroclimate Scenario`[70] = 'Horizontes'
tree_table$`Competition Scenario`[70] = 'Established'
tree_table$Sensitivity[70] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, km){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8 * 10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = km # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 300
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, km = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[69] = 'km'
liana_table$Mean[69] = mu
liana_table$`Hydroclimate Scenario`[69] = 'BCI'
liana_table$`Competition Scenario`[69] = 'Established'
liana_table$Sensitivity[69] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 300
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, km = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[70] = 'km'
liana_table$Mean[70] = mu
liana_table$`Hydroclimate Scenario`[70] = 'Horizontes'
liana_table$`Competition Scenario`[70] = 'Established'
liana_table$Sensitivity[70] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 300
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, km = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[71] = 'km'
tree_table$Mean[71] = mu
tree_table$`Hydroclimate Scenario`[71] = 'BCI'
tree_table$`Competition Scenario`[71] = 'Invasion'
tree_table$Sensitivity[71] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 300
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, km = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[72] = 'km'
tree_table$Mean[72] = mu
tree_table$`Hydroclimate Scenario`[72] = 'Horizontes'
tree_table$`Competition Scenario`[72] = 'Invasion'
tree_table$Sensitivity[72] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 300
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, km = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[71] = 'km'
liana_table$Mean[71] = mu
liana_table$`Hydroclimate Scenario`[71] = 'BCI'
liana_table$`Competition Scenario`[71] = 'Invasion'
liana_table$Sensitivity[71] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 300
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, km = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[72] = 'km'
liana_table$Mean[72] = mu
liana_table$`Hydroclimate Scenario`[72] = 'Horizontes'
liana_table$`Competition Scenario`[72] = 'Invasion'
liana_table$Sensitivity[72] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### D0 ####
############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, D0){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8 * 10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = D0 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 350
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, D0 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[73] = 'D0'
tree_table$Mean[73] = mu
tree_table$`Hydroclimate Scenario`[73] = 'BCI'
tree_table$`Competition Scenario`[73] = 'Established'
tree_table$Sensitivity[73] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 350
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, D0 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[74] = 'D0'
tree_table$Mean[74] = mu
tree_table$`Hydroclimate Scenario`[74] = 'Horizontes'
tree_table$`Competition Scenario`[74] = 'Established'
tree_table$Sensitivity[74] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, D0){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8 * 10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = D0 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 350
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, D0 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[73] = 'D0'
liana_table$Mean[73] = mu
liana_table$`Hydroclimate Scenario`[73] = 'BCI'
liana_table$`Competition Scenario`[73] = 'Established'
liana_table$Sensitivity[73] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 350
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, D0 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[74] = 'D0'
liana_table$Mean[74] = mu
liana_table$`Hydroclimate Scenario`[74] = 'Horizontes'
liana_table$`Competition Scenario`[74] = 'Established'
liana_table$Sensitivity[74] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 350
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, D0 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[75] = 'D0'
tree_table$Mean[75] = mu
tree_table$`Hydroclimate Scenario`[75] = 'BCI'
tree_table$`Competition Scenario`[75] = 'Invasion'
tree_table$Sensitivity[75] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 350
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, D0 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[76] = 'D0'
tree_table$Mean[76] = mu
tree_table$`Hydroclimate Scenario`[76] = 'Horizontes'
tree_table$`Competition Scenario`[76] = 'Invasion'
tree_table$Sensitivity[76] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 350
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, D0 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[75] = 'D0'
liana_table$Mean[75] = mu
liana_table$`Hydroclimate Scenario`[75] = 'BCI'
liana_table$`Competition Scenario`[75] = 'Invasion'
liana_table$Sensitivity[75] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 350
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, D0 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[76] = 'D0'
liana_table$Mean[76] = mu
liana_table$`Hydroclimate Scenario`[76] = 'Horizontes'
liana_table$`Competition Scenario`[76] = 'Invasion'
liana_table$Sensitivity[76] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

###############
#### gamma ####
###############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, gamma){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8 * 10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = gamma # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 125
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[77] = 'gamma'
tree_table$Mean[77] = mu
tree_table$`Hydroclimate Scenario`[77] = 'BCI'
tree_table$`Competition Scenario`[77] = 'Established'
tree_table$Sensitivity[77] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 125
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, gamma = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[78] = 'gamma'
tree_table$Mean[78] = mu
tree_table$`Hydroclimate Scenario`[78] = 'Horizontes'
tree_table$`Competition Scenario`[78] = 'Established'
tree_table$Sensitivity[78] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, gamma){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8 * 10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = gamma # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 125
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[77] = 'gamma'
liana_table$Mean[77] = mu
liana_table$`Hydroclimate Scenario`[77] = 'BCI'
liana_table$`Competition Scenario`[77] = 'Established'
liana_table$Sensitivity[77] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 125
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[78] = 'gamma'
liana_table$Mean[78] = mu
liana_table$`Hydroclimate Scenario`[78] = 'Horizontes'
liana_table$`Competition Scenario`[78] = 'Established'
liana_table$Sensitivity[78] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 125
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[79] = 'gamma'
tree_table$Mean[79] = mu
tree_table$`Hydroclimate Scenario`[79] = 'BCI'
tree_table$`Competition Scenario`[79] = 'Invasion'
tree_table$Sensitivity[79] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 125
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[80] = 'gamma'
tree_table$Mean[80] = mu
tree_table$`Hydroclimate Scenario`[80] = 'Horizontes'
tree_table$`Competition Scenario`[80] = 'Invasion'
tree_table$Sensitivity[80] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 125
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[79] = 'gamma'
liana_table$Mean[79] = mu
liana_table$`Hydroclimate Scenario`[79] = 'BCI'
liana_table$`Competition Scenario`[79] = 'Invasion'
liana_table$Sensitivity[79] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 125
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, gamma = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[80] = 'gamma'
liana_table$Mean[80] = mu
liana_table$`Hydroclimate Scenario`[80] = 'Horizontes'
liana_table$`Competition Scenario`[80] = 'Invasion'
liana_table$Sensitivity[80] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### rd ####
############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, rd){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8 * 10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = rd
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rd = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[81] = 'rd'
tree_table$Mean[81] = mu
tree_table$`Hydroclimate Scenario`[81] = 'BCI'
tree_table$`Competition Scenario`[81] = 'Established'
tree_table$Sensitivity[81] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, rd = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[82] = 'rd'
tree_table$Mean[82] = mu
tree_table$`Hydroclimate Scenario`[82] = 'Horizontes'
tree_table$`Competition Scenario`[82] = 'Established'
tree_table$Sensitivity[82] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, rd){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8 * 10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = rd
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rd = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[81] = 'rd'
liana_table$Mean[81] = mu
liana_table$`Hydroclimate Scenario`[81] = 'BCI'
liana_table$`Competition Scenario`[81] = 'Established'
liana_table$Sensitivity[81] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rd = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[82] = 'rd'
liana_table$Mean[82] = mu
liana_table$`Hydroclimate Scenario`[82] = 'Horizontes'
liana_table$`Competition Scenario`[82] = 'Established'
liana_table$Sensitivity[82] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rd = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[83] = 'rd'
tree_table$Mean[83] = mu
tree_table$`Hydroclimate Scenario`[83] = 'BCI'
tree_table$`Competition Scenario`[83] = 'Invasion'
tree_table$Sensitivity[83] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rd = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[84] = 'rd'
tree_table$Mean[84] = mu
tree_table$`Hydroclimate Scenario`[84] = 'Horizontes'
tree_table$`Competition Scenario`[84] = 'Invasion'
tree_table$Sensitivity[84] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rd = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[83] = 'rd'
liana_table$Mean[83] = mu
liana_table$`Hydroclimate Scenario`[83] = 'BCI'
liana_table$`Competition Scenario`[83] = 'Invasion'
liana_table$Sensitivity[83] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rd = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[84] = 'rd'
liana_table$Mean[84] = mu
liana_table$`Hydroclimate Scenario`[84] = 'Horizontes'
liana_table$`Competition Scenario`[84] = 'Invasion'
liana_table$Sensitivity[84] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### rx ####
############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, rx){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8 * 10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = rx
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rx = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[85] = 'rx'
tree_table$Mean[85] = mu
tree_table$`Hydroclimate Scenario`[85] = 'BCI'
tree_table$`Competition Scenario`[85] = 'Established'
tree_table$Sensitivity[85] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, rx = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[86] = 'rx'
tree_table$Mean[86] = mu
tree_table$`Hydroclimate Scenario`[86] = 'Horizontes'
tree_table$`Competition Scenario`[86] = 'Established'
tree_table$Sensitivity[86] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, rx){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8 * 10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = rx
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rx = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[85] = 'rx'
liana_table$Mean[85] = mu
liana_table$`Hydroclimate Scenario`[85] = 'BCI'
liana_table$`Competition Scenario`[85] = 'Established'
liana_table$Sensitivity[85] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rx = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[86] = 'rx'
liana_table$Mean[86] = mu
liana_table$`Hydroclimate Scenario`[86] = 'Horizontes'
liana_table$`Competition Scenario`[86] = 'Established'
liana_table$Sensitivity[86] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rx = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[87] = 'rx'
tree_table$Mean[87] = mu
tree_table$`Hydroclimate Scenario`[87] = 'BCI'
tree_table$`Competition Scenario`[87] = 'Invasion'
tree_table$Sensitivity[87] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rx = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[88] = 'rx'
tree_table$Mean[88] = mu
tree_table$`Hydroclimate Scenario`[88] = 'Horizontes'
tree_table$`Competition Scenario`[88] = 'Invasion'
tree_table$Sensitivity[88] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rx = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[87] = 'rx'
liana_table$Mean[87] = mu
liana_table$`Hydroclimate Scenario`[87] = 'BCI'
liana_table$`Competition Scenario`[87] = 'Invasion'
liana_table$Sensitivity[87] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 3
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rx = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[88] = 'rx'
liana_table$Mean[88] = mu
liana_table$`Hydroclimate Scenario`[88] = 'Horizontes'
liana_table$`Competition Scenario`[88] = 'Invasion'
liana_table$Sensitivity[88] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### rp ####
############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, rp){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8 * 10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = rp
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 11.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rp = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[89] = 'rp'
tree_table$Mean[89] = mu
tree_table$`Hydroclimate Scenario`[89] = 'BCI'
tree_table$`Competition Scenario`[89] = 'Established'
tree_table$Sensitivity[89] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 11.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, rp = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[90] = 'rp'
tree_table$Mean[90] = mu
tree_table$`Hydroclimate Scenario`[90] = 'Horizontes'
tree_table$`Competition Scenario`[90] = 'Established'
tree_table$Sensitivity[90] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, rp){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8 * 10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = rp
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 11.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rp = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[89] = 'rp'
liana_table$Mean[89] = mu
liana_table$`Hydroclimate Scenario`[89] = 'BCI'
liana_table$`Competition Scenario`[89] = 'Established'
liana_table$Sensitivity[89] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 11.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rp = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[90] = 'rp'
liana_table$Mean[90] = mu
liana_table$`Hydroclimate Scenario`[90] = 'Horizontes'
liana_table$`Competition Scenario`[90] = 'Established'
liana_table$Sensitivity[90] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 11.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rp = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[91] = 'rp'
tree_table$Mean[91] = mu
tree_table$`Hydroclimate Scenario`[91] = 'BCI'
tree_table$`Competition Scenario`[91] = 'Invasion'
tree_table$Sensitivity[91] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 11.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rp = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[92] = 'rp'
tree_table$Mean[92] = mu
tree_table$`Hydroclimate Scenario`[92] = 'Horizontes'
tree_table$`Competition Scenario`[92] = 'Invasion'
tree_table$Sensitivity[92] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 11.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rp = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[91] = 'rp'
liana_table$Mean[91] = mu
liana_table$`Hydroclimate Scenario`[91] = 'BCI'
liana_table$`Competition Scenario`[91] = 'Invasion'
liana_table$Sensitivity[91] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 11.5
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, rp = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[92] = 'rp'
liana_table$Mean[92] = mu
liana_table$`Hydroclimate Scenario`[92] = 'Horizontes'
liana_table$`Competition Scenario`[92] = 'Invasion'
liana_table$Sensitivity[92] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

###################
#### stem.turn ####
###################

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, stem.turn){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8 * 10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = stem.turn / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 0.02
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, stem.turn = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[93] = 'stem.turn'
tree_table$Mean[93] = mu
tree_table$`Hydroclimate Scenario`[93] = 'BCI'
tree_table$`Competition Scenario`[93] = 'Established'
tree_table$Sensitivity[93] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.02
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, stem.turn = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[94] = 'stem.turn'
tree_table$Mean[94] = mu
tree_table$`Hydroclimate Scenario`[94] = 'Horizontes'
tree_table$`Competition Scenario`[94] = 'Established'
tree_table$Sensitivity[94] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, stem.turn){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8 * 10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = stem.turn / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 0.1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, stem.turn = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[93] = 'stem.turn'
liana_table$Mean[93] = mu
liana_table$`Hydroclimate Scenario`[93] = 'BCI'
liana_table$`Competition Scenario`[93] = 'Established'
liana_table$Sensitivity[93] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, stem.turn = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[94] = 'stem.turn'
liana_table$Mean[94] = mu
liana_table$`Hydroclimate Scenario`[94] = 'Horizontes'
liana_table$`Competition Scenario`[94] = 'Established'
liana_table$Sensitivity[94] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 0.02
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, stem.turn = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[95] = 'stem.turn'
tree_table$Mean[95] = mu
tree_table$`Hydroclimate Scenario`[95] = 'BCI'
tree_table$`Competition Scenario`[95] = 'Invasion'
tree_table$Sensitivity[95] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.02
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, stem.turn = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[96] = 'stem.turn'
tree_table$Mean[96] = mu
tree_table$`Hydroclimate Scenario`[96] = 'Horizontes'
tree_table$`Competition Scenario`[96] = 'Invasion'
tree_table$Sensitivity[96] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 0.1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, stem.turn = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[95] = 'stem.turn'
liana_table$Mean[95] = mu
liana_table$`Hydroclimate Scenario`[95] = 'BCI'
liana_table$`Competition Scenario`[95] = 'Invasion'
liana_table$Sensitivity[95] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 0.1
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, stem.turn = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[96] = 'stem.turn'
liana_table$Mean[96] = mu
liana_table$`Hydroclimate Scenario`[96] = 'Horizontes'
liana_table$`Competition Scenario`[96] = 'Invasion'
liana_table$Sensitivity[96] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

############
#### a1 ####
############

## Rerun model with sap.frac as an input
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, a1){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8 * 10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = a1 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

## Tree, established, BCI

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, a1 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[97] = 'a1'
tree_table$Mean[97] = mu
tree_table$`Hydroclimate Scenario`[97] = 'BCI'
tree_table$`Competition Scenario`[97] = 'Established'
tree_table$Sensitivity[97] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 100 * 4
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = NULL, b2 = NULL, a1 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[98] = 'a1'
tree_table$Mean[98] = mu
tree_table$`Hydroclimate Scenario`[98] = 'Horizontes'
tree_table$`Competition Scenario`[98] = 'Established'
tree_table$Sensitivity[98] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, established, BCI

# Rerun liana model as above
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, a1){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8 * 10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = a1 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, a1 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[97] = 'a1'
liana_table$Mean[97] = mu
liana_table$`Hydroclimate Scenario`[97] = 'BCI'
liana_table$`Competition Scenario`[97] = 'Established'
liana_table$Sensitivity[97] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, established, Horizontes

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, a1 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[98] = 'a1'
liana_table$Mean[98] = mu
liana_table$`Hydroclimate Scenario`[98] = 'Horizontes'
liana_table$`Competition Scenario`[98] = 'Established'
liana_table$Sensitivity[98] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Tree, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, a1 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[99] = 'a1'
tree_table$Mean[99] = mu
tree_table$`Hydroclimate Scenario`[99] = 'BCI'
tree_table$`Competition Scenario`[99] = 'Invasion'
tree_table$Sensitivity[99] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Tree, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
tree.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, a1 = altered[i])
      sens[mo] = tree.out[1]
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req[i] = K
      break
    }
    tree_count = tree_count + 1
  }
}

# Add to table
tree_table$Parameter[100] = 'a1'
tree_table$Mean[100] = mu
tree_table$`Hydroclimate Scenario`[100] = 'Horizontes'
tree_table$`Competition Scenario`[100] = 'Invasion'
tree_table$Sensitivity[100] = abs(tree.req[3] - tree.req[1]) / tree.req[2]

## Liana, invasion, BCI

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, a1 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[99] = 'a1'
liana_table$Mean[99] = mu
liana_table$`Hydroclimate Scenario`[99] = 'BCI'
liana_table$`Competition Scenario`[99] = 'Invasion'
liana_table$Sensitivity[99] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

## Liana, invasion, Horizontes

# Mean and 50% adjustment of parameter
mu = 30
adjust = mu * 0.5

# +/- 50%
altered = c(mu - adjust, mu, mu + adjust)

# Storage
liana.req = c()

# Run model
for(i in 1:3){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD[mo,]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, SLA = NULL, b1 = NULL, b2 = NULL, a1 = altered[i])
      sens[mo] = lianaout[1]
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req[i] = K
      break
    }
    liana_count = liana_count + 1
  }
}

liana_table$Parameter[100] = 'a1'
liana_table$Mean[100] = mu
liana_table$`Hydroclimate Scenario`[100] = 'Horizontes'
liana_table$`Competition Scenario`[100] = 'Invasion'
liana_table$Sensitivity[100] = abs(liana.req[3] - liana.req[1]) / liana.req[2]

##############
#### Save ####
##############

# Tree number of simulations = 25,326
tree_count
# Liana number of simulatinos = 407,282
liana_count
# Total = 432,608

colnames(tree_table) = c('Parameter', 'Mean', 'Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity')
colnames(liana_table) = c('Parameter', 'Mean', 'Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity')

save(tree_table, liana_table, file = 'sensitivity_raw_tables.RData')

#######################
#### Create Tables ####
#######################

tree_table_1 = tree_table[1:24,]
tree_table_2 = tree_table[25:64,]
tree_table_3 = tree_table[65:100,]

tree_table_1 %>% gt() %>%
  fmt_number(columns = vars('Mean', 'Sensitivity'), decimals = 2) %>%
  tab_header(
    title = md('Tree Parameter Sensitivity')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(60),
    heading.title.font.size = 22) %>%
  cols_width(
    vars('Mean') ~ px(80),
    vars('Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity') ~ px(150),
    vars('Parameter') ~ px(120)) %>%
  cols_align(align = 'center', columns = vars('Parameter', 'Mean', 'Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity')) %>%
  tab_style(
    style = cell_fill(color = 'lightyellow'),
    locations = cells_body(rows = Sensitivity > 0.5)) %>%
  gtsave(filename = 'Plots/tree_sensitivity_1.png')

tree_table_2 %>% gt() %>%
  fmt_number(columns = vars('Mean', 'Sensitivity'), decimals = 2) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(60),
    heading.title.font.size = 22) %>%
  cols_width(
    vars('Mean') ~ px(80),
    vars('Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity') ~ px(150),
    vars('Parameter') ~ px(120)) %>%
  cols_align(align = 'center', columns = vars('Parameter', 'Mean', 'Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity')) %>%
  tab_style(
    style = cell_fill(color = 'lightyellow'),
    locations = cells_body(rows = Sensitivity > 0.5)) %>%
  gtsave(filename = 'Plots/tree_sensitivity_2.png')

tree_table_3 %>% gt() %>%
  fmt_number(columns = vars('Mean', 'Sensitivity'), decimals = 2) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(60),
    heading.title.font.size = 22) %>%
  cols_width(
    vars('Mean') ~ px(80),
    vars('Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity') ~ px(150),
    vars('Parameter') ~ px(120)) %>%
  cols_align(align = 'center', columns = vars('Parameter', 'Mean', 'Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity')) %>%
  tab_style(
    style = cell_fill(color = 'lightyellow'),
    locations = cells_body(rows = Sensitivity > 0.5)) %>%
  gtsave(filename = 'Plots/tree_sensitivity_3.png')


liana_table_1 = liana_table[1:24,]
liana_table_2 = liana_table[25:64,]
liana_table_3 = liana_table[65:100,]

liana_table_1 %>% gt() %>%
  fmt_number(columns = vars('Mean', 'Sensitivity'), decimals = 2) %>%
  tab_header(
    title = md('Liana Parameter Sensitivity')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(60),
    heading.title.font.size = 22) %>%
  cols_width(
    vars('Mean') ~ px(80),
    vars('Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity') ~ px(150),
    vars('Parameter') ~ px(120)) %>%
  cols_align(align = 'center', columns = vars('Parameter', 'Mean', 'Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity')) %>%
  tab_style(
    style = cell_fill(color = 'lightyellow'),
    locations = cells_body(rows = Sensitivity > 0.5)) %>%
  gtsave(filename = 'Plots/liana_sensitivity_1.png')

liana_table_2 %>% gt() %>%
  fmt_number(columns = vars('Mean', 'Sensitivity'), decimals = 2) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(60),
    heading.title.font.size = 22) %>%
  cols_width(
    vars('Mean') ~ px(80),
    vars('Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity') ~ px(150),
    vars('Parameter') ~ px(120)) %>%
  cols_align(align = 'center', columns = vars('Parameter', 'Mean', 'Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity')) %>%
  tab_style(
    style = cell_fill(color = 'lightyellow'),
    locations = cells_body(rows = Sensitivity > 0.5)) %>%
  gtsave(filename = 'Plots/liana_sensitivity_2.png')

liana_table_3 %>% gt() %>%
  fmt_number(columns = vars('Mean', 'Sensitivity'), decimals = 2) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(60),
    heading.title.font.size = 22) %>%
  cols_width(
    vars('Mean') ~ px(80),
    vars('Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity') ~ px(150),
    vars('Parameter') ~ px(120)) %>%
  cols_align(align = 'center', columns = vars('Parameter', 'Mean', 'Hydroclimate Scenario', 'Competition Scenario', 'Sensitivity')) %>%
  tab_style(
    style = cell_fill(color = 'lightyellow'),
    locations = cells_body(rows = Sensitivity > 0.5)) %>%
  gtsave(filename = 'Plots/liana_sensitivity_3.png')

########################
#### Create Figures ####
########################
tree_table$`Hydroclimate Scenario`[which(tree_table$`Hydroclimate Scenario` == 'BCI')] = 'Wettest'
tree_table$`Hydroclimate Scenario`[which(tree_table$`Hydroclimate Scenario` == 'Horizontes')] = 'Driest'

tree_DBH = tree_table %>%
  filter(Parameter == 'DBH')

pl1 = ggplot() +
  geom_bar(aes(x = tree_DBH$`Competition Scenario`, y = tree_DBH$Sensitivity, fill = tree_DBH$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('DBH') +
  ylim(c(0, max(tree_DBH$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl1

tree_Ca = tree_table %>%
  filter(Parameter == 'Ca')

pl2 = ggplot() +
  geom_bar(aes(x = tree_Ca$`Competition Scenario`, y = tree_Ca$Sensitivity, fill = tree_Ca$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('Ca') +
  ylim(c(0, max(tree_DBH$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl2

tree_b2 = tree_table %>%
  filter(Parameter == 'b2')

pl3 = ggplot() +
  geom_bar(aes(x = tree_b2$`Competition Scenario`, y = tree_b2$Sensitivity, fill = tree_b2$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('b2') +
  ylim(c(0, max(tree_DBH$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl3

tree_lx = tree_table %>%
  filter(Parameter == 'Lx')

pl4 = ggplot() +
  geom_bar(aes(x = tree_lx$`Competition Scenario`, y = tree_lx$Sensitivity, fill = tree_lx$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('Lx') +
  ylim(c(0, max(tree_DBH$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl4

tree_al = tree_table %>%
  filter(Parameter == 'AL')

pl5 = ggplot() +
  geom_bar(aes(x = tree_al$`Competition Scenario`, y = tree_al$Sensitivity, fill = tree_al$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('AL') +
  ylim(c(0, max(tree_DBH$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl5

tree_sla = tree_table %>%
  filter(Parameter == 'SLA')

pl6 = ggplot() +
  geom_bar(aes(x = tree_sla$`Competition Scenario`, y = tree_sla$Sensitivity, fill = tree_sla$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('SLA') +
  ylim(c(0, max(tree_DBH$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl6

tree_q = tree_table %>%
  filter(Parameter == 'q')

pl7 = ggplot() +
  geom_bar(aes(x = tree_q$`Competition Scenario`, y= tree_q$Sensitivity, fill = tree_q$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('q') +
  ylim(c(0, max(tree_DBH$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl7

tree_rd = tree_table %>%
  filter(Parameter == 'rd')

pl8 = ggplot() +
  geom_bar(aes(x = tree_rd$`Competition Scenario`, y = tree_rd$Sensitivity, fill = tree_rd$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('rd') +
  ylim(c(0, max(tree_DBH$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl8

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend = get_legend(pl1)

ga = plot_grid(pl1 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')), 
               pl2 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')), 
               pl3 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')), 
               pl4 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')), 
               pl5 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')), 
               pl6 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')),
               pl7 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')),
               pl8 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')), nrow = 4)
ga2 = plot_grid(ga, legend, nrow = 1, rel_widths = c(0.8, 0.2))
ga2

ggsave(ga2, filename = 'Plots/tree_sensitivity_hist.jpeg', height = 20, width = 14, units = 'in')

liana_table$`Hydroclimate Scenario`[which(liana_table$`Hydroclimate Scenario` == 'Horizontes')] = 'Driest'
liana_table$`Hydroclimate Scenario`[which(liana_table$`Hydroclimate Scenario` == 'BCI')] = 'Wettest'

liana_dbh = liana_table %>%
  filter(Parameter == 'DBH')

pl1 = ggplot() +
  geom_bar(aes(x = liana_dbh$`Competition Scenario`, y = liana_dbh$Sensitivity, fill = liana_dbh$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('DBH') +
  ylim(c(0, max(liana_dbh$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl1

liana_ca = liana_table %>%
  filter(Parameter == 'Ca')

pl2 = ggplot() +
  geom_bar(aes(x = liana_ca$`Competition Scenario`, y = liana_ca$Sensitivity, fill = liana_ca$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('Ca') +
  ylim(c(0, max(liana_dbh$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl2

liana_b2 = liana_table %>%
  filter(Parameter == 'b2')

pl3 = ggplot() +
  geom_bar(aes(x = liana_b2$`Competition Scenario`, y = liana_b2$Sensitivity, fill = liana_b2$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('b2') +
  ylim(c(0, max(liana_dbh$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl3

liana_lx = liana_table %>%
  filter(Parameter == 'Lx')

pl4 = ggplot() +
  geom_bar(aes(x = liana_lx$`Competition Scenario`, y = liana_lx$Sensitivity, fill = liana_lx$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('Lx') +
  ylim(c(0, max(liana_dbh$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl4

liana_al = liana_table %>%
  filter(Parameter == 'AL')

pl5 = ggplot() +
  geom_bar(aes(x = liana_al$`Competition Scenario`, y = liana_al$Sensitivity, fill = liana_al$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('AL') +
  ylim(c(0, max(liana_dbh$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl5

liana_sla = liana_table %>%
  filter(Parameter == 'SLA')

pl6 = ggplot() +
  geom_bar(aes(x = liana_sla$`Competition Scenario`, y = liana_sla$Sensitivity, fill = liana_sla$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 labels = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('SLA') +
  ylim(c(0, max(liana_dbh$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl6

liana_q = liana_table %>%
  filter(Parameter == 'q')

pl7 = ggplot() +
  geom_bar(aes(x = liana_q$`Competition Scenario`, y = liana_q$Sensitivity, fill = liana_q$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 label = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('q') +
  ylim(c(0, max(liana_dbh$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl7

liana_rd = liana_table %>%
  filter(Parameter == 'rd')

pl8 = ggplot() +
  geom_bar(aes(x = liana_rd$`Competition Scenario`, y = liana_rd$Sensitivity, fill = liana_rd$`Hydroclimate Scenario`), stat = 'identity', position = 'dodge') +
  xlab('Competition Scenario') + ylab('Sensitivity') +
  scale_fill_npg(name = 'Hydroclimate\nScenario',
                 label = c('Driest' = 'Horizontes', 'Wettest' = 'BCI')) +
  theme_linedraw() +
  ggtitle('rd') +
  ylim(c(0, max(liana_dbh$Sensitivity))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl8

ga = plot_grid(pl1 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')),
               pl2 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')),
               pl3 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')),
               pl4 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')),
               pl5 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')),
               pl6 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')),
               pl7 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')),
               pl8 + theme(legend.position = 'none', plot.margin = unit(c(15, 15, 15, 15), 'pt')), nrow = 4)
ga2 = plot_grid(ga, legend, nrow = 1, rel_widths = c(0.8, 0.2))
ga2

ggsave(ga2, filename = 'Plots/liana_sensitivity_hist.jpeg', height = 20, width = 14, unit = 'in')
