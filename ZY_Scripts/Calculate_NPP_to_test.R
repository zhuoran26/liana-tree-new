#this is a copy of calculating NPP from the original scripts
rm(list = ls())

source('scripts/NPP_models.R')
source('scripts/NPP_models_wCO2.R')

# Load parameters
load('data/param.input.RData')

# Load met drivers
#Question where did the data go?
load('data/bci_met_mxh.RData')

nK = 75000 #question
ks = c(trees_K$K / 10, lianas_K$K)
varyk = seq(min(ks), max(ks), length.out = nK) #make a sequence of k

# Interpolate DBH (and remove highest DBH) #Question where is DBH removed here
dbhs = seq(liana.dbh.low, liana.dbh.high, length.out = 50)

## Tree: 90% canopy, mean DBH, BCI

treeouteach = c()
tree_out_r1 = c()

for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) #days in each months
    month = j
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 2 * 100
    frac.tree.al = 0.9
    
    psis = SWP[j] #soil water potential
    
    D = VPD[j,]
    
    if(j == 1){ #the first month
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5] #question: where does this come from
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, 
                            D_vec = D, tree.k = K, nday = nday, 
                            month = month, tot.al = tot.al, 
                            frac.tree.al = frac.tree.al, tree.height = tree.height)
    #tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, 
    #b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, 
    #Vm = NULL)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach) #make a year
  if(sumtreeout > 0){
    tree_out_r1[1] = sumtreeout
    tree_out_r1[2] = K
    
    break
  }
}

##Liana: 10% canopy, varied DBH, BCI

lianaouteach.org = c()
#liana_out_r1 = matrix(, nrow = length(dbhs), ncol = 2)

    for(j in 1){
      nday = 31
      
      
      dbh = dbhs[1]
      
      tot.al = 2 * 100
      frac.liana.al = 0.1
      
      psis = SWP[1]
      
      D = VPD[1,]#average diurnal cycle
      #question:what is liana length and why is the 4th one for liana, 5 for tree
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[1]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, 
                                nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach.org[j] = liana.out[1]#return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
    }

vec=rep(0,10)
for (i in 1:20){
  vec=vec[1]+1
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
load('data/horizontes_met_mxh.RData')

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
df1$sap_area = rep(NA, nrow(df1))
for(i in 1:nrow(df1)){
  df1$sap_area[i] = min(((pi*df1$DBH[i]^2) / 4), (2.41 * (df1$DBH[i]/2)^1.97)) #sapwood area in cm^2
}
df1$la = c(rep(0.1 * 2 * 100, nrow(df1) - 1), 0.9 * 2 * 100) #leaf area in m^2
df1$hv = df1$sap_area / df1$la

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

df1 = as.data.frame(cbind(df1, rep(1, nrow(df1))))
colnames(df1)[7] = 'df'
df2 = as.data.frame(cbind(df2, rep(2, nrow(df2))))
colnames(df2)[7] = 'df'
df3 = as.data.frame(cbind(df3, rep(3, nrow(df3))))
colnames(df3)[7] = 'df'
df4 = as.data.frame(cbind(df4, rep(4, nrow(df4))))
colnames(df4)[7] = 'df'
df = as.data.frame(rbind(df1, df2, df3, df4))

df_new = df
df_new$ratio = NA

# Dividing liana Kreq by the reference scenario tree Kreq

for(i in 1:nrow(df_new)){
  if(df_new$df[i] == 1 & df_new$PFT[i] == 'Liana'){
    df_new$ratio[i] = df_new$Kreq[i] / df_new$Kreq[51]
  }
  if(df_new$df[i] == 2 & df_new$PFT[i] == 'Liana'){
    df_new$ratio[i] = df_new$Kreq[i] / df_new$Kreq[102]
  }
  if(df_new$df[i] == 3 & df_new$PFT[i] == 'Liana'){
    df_new$ratio[i] = df_new$Kreq[i] / df_new$Kreq[153]
  }
  if(df_new$df[i] == 4 & df_new$PFT[i] == 'Liana'){
    df_new$ratio[i] = df_new$Kreq[i] / df_new$Kreq[204]
  }
}

df_est = df %>%
  subset(df %in% c(2, 4))
