#this scrip is used to calculate NPP with different VPD
rm(list = ls())
source("/Users/zhuoranyu/Desktop/liana-tree-test/ZY_scripts/NPP_functions.R")
library(tidyverse)
# Interpolate DBH (and remove highest DBH) 
dbhs = seq(liana.dbh.low, liana.dbh.high, length.out = 50)

dbh = dbhs[25]
tot.al = 2 * 100
frac.liana.al = 0.1

tree.dbh = mean.hori.data.tree.dbh
liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
K=5805.55556

each_day = c(rep(0,times=6),rep(1, times=12),rep(0,times=6))
year_day_light = c(rep(each_day, times=366))


month = c(1:12)
SWP_year <- cbind(SWP,month)
SWP_year <- data.frame(SWP_year)
vpd_year <- left_join(vpd_year,SWP_year)
year_vpd_org = c(c(rep(VPD[1,], times=31)),c(rep(VPD[2,], times=29)),c(rep(VPD[3,], times=31)),c(rep(VPD[4,], times=30)),
             c(rep(VPD[5,], times=31)),c(rep(VPD[6,], times=30)),c(rep(VPD[7,], times=31)),c(rep(VPD[8,], times=31)),
             c(rep(VPD[9,], times=30)),c(rep(VPD[10,], times=31)),c(rep(VPD[11,], times=30)),c(rep(VPD[12,], times=31)))

vpd_year <- cbind(vpd_year,year_vpd_org,year_day_light)
liana_npp_org_VPD = rep(0,times=length(year_day_light))
liana_npp_var_VPD = rep(0,times=length(year_day_light))

for (i in 1:length(year_day_light)){
  psis= vpd_year$SWP
  day_light = vpd_year$year_day_light
  liana_npp_org_VPD[i] = liana.NPP.mxh.hour(dbh = dbh, liana.k = K, tot.al = tot.al, nday=nday,month=month,
                                              frac.liana.al = frac.liana.al, liana.length = liana.length,
                                         day_light=day_light,psis = psis,D_vec=vpd_year$year_vpd_org)
  
  liana_npp_var_VPD[i] = liana.NPP.mxh.hour(dbh = dbh, liana.k = K, tot.al = tot.al, nday=nday,month=month,
                                            frac.liana.al = frac.liana.al, liana.length = liana.length,
                                            day_light=day_light,psis = psis,D_vec=vpd_year$vpd_pa)
  
}

liana.data <- cbind(vpd_year,liana_npp_org_VPD,liana_npp_var_VPD)


liana.monthly <- liana.data %>%
  group_by(month) %>%
  summarise( 
    liana.npp.month.org = sum(liana_npp_org_VPD),
    liana.npp.month.var = sum(liana_npp_var_VPD)
    ) 
 summary(liana.monthly$liana.npp.month.var)
plot(liana.monthly$month,liana.monthly$liana.npp.month.org,ylim =c(-0.5,3))
points(liana.monthly$month,liana.monthly$liana.npp.month.var,col="green")
sum(liana_npp_org_VPD)
sum(liana_npp_var_VPD)
