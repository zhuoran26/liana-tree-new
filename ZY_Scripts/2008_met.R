#code to calculate and save VPD from 2008
rm(list = ls())
library(rhdf5)
library(tidyverse)
library(gtools)

setwd("/Users/zhuoranyu/Desktop/liana-tree-test/Data")

h5_list <- list.files(pattern="*.h5")
print(h5_list) #check the list of files

jan.met = "CRonesite_2008JAN.h5"  ## HDF5 file
feb.met = "CRonesite_2008FEB.h5"
mar.met = "CRonesite_2008MAR.h5"
apr.met = "CRonesite_2008APR.h5"
may.met = "CRonesite_2008MAY.h5"
jun.met = "CRonesite_2008JUN.h5"
jul.met = "CRonesite_2008JUL.h5"
aug.met = "CRonesite_2008AUG.h5"
sep.met = "CRonesite_2008SEP.h5"
oct.met = "CRonesite_2008OCT.h5"
nov.met = "CRonesite_2008NOV.h5"
dec.met = "CRonesite_2008DEC.h5"

read.h5 <- function(file){
    sh <- h5read(file, 'sh')
    temp_k <- h5read(file,"tmp")
    pres_pa <- h5read(file,"pres")
    df <-data.frame(sh,temp_k,pres_pa)
    #change temp from K to C
    df$temp_c <- df$temp_k-273.15
    #change pressure from Pa to kPa
    df$pres_kpa <- df$pres_pa/1000
    #calculate vapor pressure (e;Pa) from atmosphere pressure (Pa) and specific humidity (kg/kg)
    df$e <-  df$sh*df$pres_pa/(0.622+0.378*df$sh)
    #calculate saturated  vapor pressure (SVP; Pa) from temperature in K
    df$es <- 6.1078*exp(17.2694*(df$temp_k-273.16)/(df$temp_k-35.86))*100
    #calculate relative humidity
    df$rh <- df$e/df$es*100
    #calculate VPD (Pa)
    df$vpd_pa <- df$es*(1 - df$rh / 100)
    df$vpd_hpa <- df$vpd_pa/100
    return(df)
  }

jan_vpd <- read.h5(jan.met)
feb_vpd <- read.h5(feb.met)
mar_vpd <- read.h5(mar.met)
apr_vpd <- read.h5(apr.met)
may_vpd <- read.h5(may.met)
jun_vpd <- read.h5(jun.met)
jul_vpd <- read.h5(jul.met)
aug_vpd <- read.h5(aug.met)
sep_vpd <- read.h5(sep.met)
oct_vpd <- read.h5(oct.met)
nov_vpd <- read.h5(nov.met)
dec_vpd <- read.h5(dec.met)

jan_vpd$month <- 1
feb_vpd$month <- 2
mar_vpd$month <- 3
apr_vpd$month <- 4
may_vpd$month <- 5
jun_vpd$month <- 6
jul_vpd$month <- 7
aug_vpd$month <- 8
sep_vpd$month <- 9
oct_vpd$month <- 10
nov_vpd$month <- 11
dec_vpd$month <- 12

vpd_year <- smartbind(jan_vpd,feb_vpd,mar_vpd,apr_vpd,may_vpd,
                      jun_vpd,jul_vpd,aug_vpd,sep_vpd,oct_vpd,nov_vpd,dec_vpd)

hours <-rep(c(1:24),times=366)
vpd_year$hour <- hours

#this is a leap year
vpd_day<-matrix(vpd_year$vpd_pa,ncol=24,byrow=TRUE)

vpd_month_mean <- vpd_year %>%
  group_by(month,hour) %>%
  summarise(monthly_mean = mean(vpd_pa))

vpd_month <-matrix(vpd_month_mean$monthly_mean,ncol = 24,byrow = T)

save(vpd_day, file = "vpd_2008_day.RData")
save(vpd_month, file = "vpd_2008_month.RData")
save(vpd_year, file = "vpd_2008_year.RData")

