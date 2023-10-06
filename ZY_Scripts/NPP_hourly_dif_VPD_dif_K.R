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
#K=5805.55556

each_day = c(rep(0,times=6),rep(1, times=12),rep(0,times=6))
year_day_light = c(rep(each_day, times=366))


month = c(1:12)
SWP_year <- cbind(SWP,month)
SWP_year <- data.frame(SWP_year)
vpd_year <- left_join(vpd_year,SWP_year)
year_vpd_org = c(c(rep(vpd_month[1,], times=31)),c(rep(vpd_month[2,], times=29)),c(rep(vpd_month[3,], times=31)),c(rep(vpd_month[4,], times=30)),
             c(rep(vpd_month[5,], times=31)),c(rep(vpd_month[6,], times=30)),c(rep(vpd_month[7,], times=31)),c(rep(vpd_month[8,], times=31)),
             c(rep(vpd_month[9,], times=30)),c(rep(vpd_month[10,], times=31)),c(rep(vpd_month[11,], times=30)),c(rep(vpd_month[12,], times=31)))

vpd_year <- cbind(vpd_year,year_vpd_org,year_day_light)
liana_npp_org_VPD = rep(0,times=length(year_day_light))
liana_npp_var_VPD = rep(0,times=length(year_day_light))

sorted.K <- sort(lianas_K$K, decreasing = FALSE)

sorted.K.mPa <- sorted.K/1000

liana_npp_org_VPD_low_k = rep(0,8784)
liana_npp_var_VPD_low_k = rep(0,8784)
liana_npp_org_VPD_mid_k = rep(0,8784)
liana_npp_var_VPD_mid_k = rep(0,8784)
liana_npp_org_VPD_high_k = rep(0,8784)
liana_npp_var_VPD_high_k = rep(0,8784)
#lowest K
for (i in 1:length(year_day_light)){
  psis= vpd_year$SWP
  day_light = vpd_year$year_day_light
  liana_npp_org_VPD_low_k[i] = liana.NPP.mxh.hour(dbh = dbh, liana.k = sorted.K.mPa[1], tot.al = tot.al, nday=nday,month=month,
                                              frac.liana.al = frac.liana.al, liana.length = liana.length,
                                         day_light=day_light,psis = psis,D_vec=vpd_year$year_vpd_org)
  
  liana_npp_var_VPD_low_k[i] = liana.NPP.mxh.hour(dbh = dbh, liana.k = sorted.K.mPa[1], tot.al = tot.al, nday=nday,month=month,
                                            frac.liana.al = frac.liana.al, liana.length = liana.length,
                                            day_light=day_light,psis = psis,D_vec=vpd_year$vpd_pa)
  liana_npp_org_VPD_mid_k[i] = liana.NPP.mxh.hour(dbh = dbh, liana.k = sorted.K.mPa[25], tot.al = tot.al, nday=nday,month=month,
                                                  frac.liana.al = frac.liana.al, liana.length = liana.length,
                                                  day_light=day_light,psis = psis,D_vec=vpd_year$year_vpd_org)
  
  liana_npp_var_VPD_mid_k[i] = liana.NPP.mxh.hour(dbh = dbh, liana.k = sorted.K.mPa[25], tot.al = tot.al, nday=nday,month=month,
                                                  frac.liana.al = frac.liana.al, liana.length = liana.length,
                                                  day_light=day_light,psis = psis,D_vec=vpd_year$vpd_pa)
  
  liana_npp_org_VPD_high_k[i] = liana.NPP.mxh.hour(dbh = dbh, liana.k = sorted.K.mPa[51], tot.al = tot.al, nday=nday,month=month,
                                                  frac.liana.al = frac.liana.al, liana.length = liana.length,
                                                  day_light=day_light,psis = psis,D_vec=vpd_year$year_vpd_org)
  
  liana_npp_var_VPD_high_k[i] = liana.NPP.mxh.hour(dbh = dbh, liana.k = sorted.K.mPa[51], tot.al = tot.al, nday=nday,month=month,
                                                  frac.liana.al = frac.liana.al, liana.length = liana.length,
                                                  day_light=day_light,psis = psis,D_vec=vpd_year$vpd_pa)
  
}

liana.data <- cbind(vpd_year,liana_npp_org_VPD_low_k,liana_npp_var_VPD_low_k,
                    liana_npp_org_VPD_mid_k,liana_npp_var_VPD_mid_k,
                    liana_npp_org_VPD_high_k,liana_npp_var_VPD_high_k)


liana.monthly <- liana.data %>%
  group_by(month) %>%
  summarise( 
    liana.npp.month.org.low = sum(liana_npp_org_VPD_low_k),
    liana.npp.month.var.low = sum(liana_npp_var_VPD_low_k),
    liana.npp.month.org.mid = sum(liana_npp_org_VPD_mid_k),
    liana.npp.month.var.mid = sum(liana_npp_var_VPD_mid_k),
    liana.npp.month.org.high = sum(liana_npp_org_VPD_high_k),
    liana.npp.month.var.high = sum(liana_npp_var_VPD_high_k)
    ) 

sum(liana.monthly$liana.npp.month.org.low)
sum(liana.monthly$liana.npp.month.var.low)
sum(liana.monthly$liana.npp.month.org.mid)
sum(liana.monthly$liana.npp.month.var.mid)
sum(liana.monthly$liana.npp.month.org.high)
sum(liana.monthly$liana.npp.month.var.high)

summary(liana.monthly$liana.npp.month.org.low)
summary(liana.monthly$liana.npp.month.var.low)
plot(liana.monthly$month,liana.monthly$liana.npp.month.org.low,ylim =c(-1,3),main = "K=28.88889",)
points(liana.monthly$month,liana.monthly$liana.npp.month.var.low,col="green")
legend(0, 2, legend=c("original", "variation"),
       col=c("black", "green"), lty=1:2, cex=0.8)

summary(liana.monthly$liana.npp.month.org.mid)
summary(liana.monthly$liana.npp.month.var.mid)
plot(liana.monthly$month,liana.monthly$liana.npp.month.org.mid,ylim =c(-1,3),main = "K=512.8449",)
points(liana.monthly$month,liana.monthly$liana.npp.month.var.mid,col="green")
legend(0, 2, legend=c("original", "variation"),
       col=c("black", "green"), lty=1:2, cex=0.8)


summary(liana.monthly$liana.npp.month.org.high)
summary(liana.monthly$liana.npp.month.var.high)
plot(liana.monthly$month,liana.monthly$liana.npp.month.org.high,ylim =c(-1,3),main = "K=5805.55556",)
points(liana.monthly$month,liana.monthly$liana.npp.month.var.high,col="green")
legend(0, 2, legend=c("original", "variation"),
       col=c("black", "green"), lty=1:2, cex=0.8)



plot(liana.data$vpd_pa,liana_npp_var_VPD_low_k)
plot(liana.data$vpd_pa,liana_npp_var_VPD_mid_k)
plot(liana.data$vpd_pa,liana_npp_var_VPD_high_k)

df[[text]] <- ifelse(df$season == i, 1, 0)
rm(dt)
dt <- matrix(0, 8784, 2) 
dt <-data.frame(dt)

for (K in sorted.K.mPa[1:2]){
  for (i in 1:length(year_day_light)){
  psis= vpd_year$SWP
  day_light = vpd_year$year_day_light
  liana_npp_org_VPD[i] = liana.NPP.mxh.hour(dbh = dbh, liana.k = K, tot.al = tot.al, nday=nday,month=month,
                                            frac.liana.al = frac.liana.al, liana.length = liana.length,
                                            day_light=day_light,psis = psis,D_vec=vpd_year$year_vpd_org)}
  
  dt <- matrix(0, 8784, 2) 
  dt <-data.frame(dt)
  dt[,K] <- liana_npp_org_VPD[i]
  colnames(dt)<-paste0("K_",K))}}
  
  dt[,K]=output
  #liana.data$liana_npp_org_VPD[text][i] = liana.NPP.mxh.hour(dbh = dbh, liana.k = K, tot.al = tot.al, nday=nday,month=month,
                                           # frac.liana.al = frac.liana.al, liana.length = liana.length,
                                           # day_light=day_light,psis = psis,D_vec=vpd_year$vpd_pa)
  
}


#####codes that did not work
#text <- paste("liana_npp_org_VPD_K", K)
#dt[,paste("K_value", K, sep="")] = rep(0, each=8784)}
#df$liana_npp_var_VPD[[text]] = rep(0,times=length(year_day_light))

