rm(list = ls())
source("/Users/zhuoranyu/Desktop/liana-tree-test/ZY_scripts/NPP_functions.R")
library(tidyverse)

#original code
nK = 75000 
ks = c(trees_K$K / 10, lianas_K$K)
varyk = seq(min(ks), max(ks), length.out = nK) #make a sequence of k

# Interpolate DBH (and remove highest DBH) 
dbhs = seq(liana.dbh.low, liana.dbh.high, length.out = 6)

tot.al = 2 * 100
frac.liana.al = 0.4
tree.dbh = mean.hori.data.tree.dbh
liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)

each_day = c(rep(0,times=6),rep(1, times=12),rep(0,times=6))
year_day_light = c(rep(each_day, times=366))

year_SWP = c(c(rep(SWP[1], times=31*24)),c(rep(SWP[2], times=29*24)),c(rep(SWP[3], times=31*24)),c(rep(SWP[4], times=30*24)),
             c(rep(SWP[5], times=31*24)),c(rep(SWP[6], times=30*24)),c(rep(SWP[7], times=31*24)),c(rep(SWP[8], times=31*24)),
             c(rep(SWP[9], times=30*24)),c(rep(SWP[10], times=31*24)),c(rep(SWP[11], times=30*24)),c(rep(SWP[12], times=31*24)))

year_vpd =vpd_year$vpd_pa


liana_npp_var_hour.1 = rep(0,8784)
#find K req
for (i in 1:length(year_day_light)){
  psis= year_SWP
  day_light = year_day_light
  liana_npp_var_hour.1[i] = liana.NPP.mxh.hour(dbh = dbhs[1], liana.k = 44100, tot.al = tot.al, nday=nday,month=month,
                                                  frac.liana.al = frac.liana.al, liana.length = liana.length,
                                                  day_light=day_light,psis = psis,D_vec=year_vpd)}
sum(liana_npp_var_hour.1)  

dbh_1 <- c(dbhs[1],44100,sum(liana_npp_var_hour.1))

liana_npp_var_hour.2 = rep(0,8784)
#find K req
for (i in 1:length(year_day_light)){
  psis= year_SWP
  day_light = year_day_light
  liana_npp_var_hour.2[i] = liana.NPP.mxh.hour(dbh = dbhs[2], liana.k = 12050, tot.al = tot.al, nday=nday,month=month,
                                               frac.liana.al = frac.liana.al, liana.length = liana.length,
                                               day_light=day_light,psis = psis,D_vec=year_vpd)}
sum(liana_npp_var_hour.2)  

dbh_2 <- c(dbhs[2],12050,sum(liana_npp_var_hour.2))

liana_npp_var_hour.3 = rep(0,8784)
#find K req
for (i in 1:length(year_day_light)){
  psis= year_SWP
  day_light = year_day_light
  liana_npp_var_hour.3[i] = liana.NPP.mxh.hour(dbh = dbhs[3], liana.k = 5600, tot.al = tot.al, nday=nday,month=month,
                                               frac.liana.al = frac.liana.al, liana.length = liana.length,
                                               day_light=day_light,psis = psis,D_vec=year_vpd)}
sum(liana_npp_var_hour.3)  

dbh_3 <- c(dbhs[3],5600,sum(liana_npp_var_hour.3))


liana_npp_var_hour.4 = rep(0,8784)
#find K req
for (i in 1:length(year_day_light)){
  psis= year_SWP
  day_light = year_day_light
  liana_npp_var_hour.4[i] = liana.NPP.mxh.hour(dbh = dbhs[4], liana.k = 3250, tot.al = tot.al, nday=nday,month=month,
                                               frac.liana.al = frac.liana.al, liana.length = liana.length,
                                               day_light=day_light,psis = psis,D_vec=year_vpd)}
sum(liana_npp_var_hour.4)  

dbh_4 <- c(dbhs[4],3250,sum(liana_npp_var_hour.4))

liana_npp_var_hour.5 = rep(0,8784)
#find K req
for (i in 1:length(year_day_light)){
  psis= year_SWP
  day_light = year_day_light
  liana_npp_var_hour.5[i] = liana.NPP.mxh.hour(dbh = dbhs[5], liana.k = 2139, tot.al = tot.al, nday=nday,month=month,
                                               frac.liana.al = frac.liana.al, liana.length = liana.length,
                                               day_light=day_light,psis = psis,D_vec=year_vpd)}
sum(liana_npp_var_hour.5)  

dbh_5 <- c(dbhs[5],2139,sum(liana_npp_var_hour.5))

liana_npp_var_hour.6 = rep(0,8784)
#find K req
for (i in 1:length(year_day_light)){
  psis= year_SWP
  day_light = year_day_light
  liana_npp_var_hour.6[i] = liana.NPP.mxh.hour(dbh = dbhs[6], liana.k = 1530, tot.al = tot.al, nday=nday,month=month,
                                               frac.liana.al = frac.liana.al, liana.length = liana.length,
                                               day_light=day_light,psis = psis,D_vec=year_vpd)}
sum(liana_npp_var_hour.6)  

dbh_6 <- c(dbhs[6],1530,sum(liana_npp_var_hour.6))

npp_hour_all <- rbind(dbh_1,dbh_2,dbh_3,dbh_4,dbh_5,dbh_6)


