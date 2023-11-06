#use thus script to calculate VPD
library(rhdf5)
library(pvldcurve)
setwd("/Users/zhuoranyu/Desktop/liana-tree-test/Data")

h5_list <- list.files(pattern="*.h5")
print(h5_list) #check the list of files

fname = "CRonesite_2008JAN.h5"  ## HDF5 file

h5_met <- data.frame("sh"=numeric(),'temp_k'=numeric(),"pres_pa"=numeric())

for (file in h5_list)
  h5_met <- data.frame("sh"=numeric(),'temp_k'=numeric(),"pres_pa"=numeric())

  h5_met$sh<-
  
    for (file in h5_list) {
      # Read the file into a data frame, using the specified delimiter
      sh <- h5read(file, 'sh')
      temp_k <- h5read(file,"tmp")
      pres_pa <- h5read(file,"pres")
      
      df <- c(sh, temp_k, pres_pa)
      gsub(".h5$", file) <-df
    }

  
h5ls(fname) ## lists the contents

## pres is the pressure (Pa)
## sh is the specific humidity (g water / g air)
## tmp is the temperature (K)

sh = h5read(fname,"sh") ## read in a variable
temp_k = h5read(fname,"tmp")
pres_pa = h5read(fname,"pres")

met <- data.frame(sh,temp_k,pres_pa)
#change temp from K to C
met$temp_c <- met$temp_k-273.15
#change pressure from Pa to kPa
met$pres_kpa <- met$pres_pa/1000


#calculate vapor pressure (e;Pa?) from atmosphere pressure (Pa) and specific humidity (kg/kg)
met$e <-  met$sh*met$pres_pa/(0.622+0.378*met$sh)
#calculate saturated  vapor pressure (SVP; Pa) from temperature in K
met$es <- 6.1078*exp(17.2694*(met$temp_k-273.16)/(met$temp_k-35.86))*100
#calculate relative humidity
met$rh <- met$e/met$es*100
#calculate VPD (Pa?)
met$vpd_pa <- met$es*(1 - met$rh / 100)
met$vpd_hpa <- met$vpd_pa/100

#check in another way
calc_VPD = VaporPressureDeficit(met, humidity = "rh",
temperature = "temp_c", atmospheric.pressure = 101325)

head(met$vpd_pa)
hist(met$vpd_pa)
head(calc_VPD$vpd)

vpd_daily<-matrix(met$vpd_pa,ncol=24,byrow=TRUE)
plot(vpd_daily[1,])
vpd_monthly <- colMeans(vpd_daily)
