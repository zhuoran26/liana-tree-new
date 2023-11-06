#check July NPP with different K
rm(list = ls())
#load in parameters and meteorological variables
load('/Users/zhuoranyu/Desktop/liana-tree-test/Data/param.input.RData')
load('/Users/zhuoranyu/Desktop/liana-tree-test/Data/bci_met_mxh.RData')
#load in daily and month variations of VPD

#This NPP model compares how the model variables would be different after
#using daily variations for each day and monthly mean variations for each day within
#the same month
load('/Users/zhuoranyu/Desktop/liana-tree-test/Data/vpd_2008_day.RData')
load('/Users/zhuoranyu/Desktop/liana-tree-test/Data/vpd_2008_month.RData')

#load packages
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

########list part of the constants used for model

#input liana dbh for model sensitivity
#Interpolate DBH (and remove highest DBH) 
dbhs = seq(liana.dbh.low, liana.dbh.high, length.out = 50)

#to make this model simple, only DBH value is chosen
dbh = dbhs[25]
tot.al = 2 * 100 # Leaf area
frac.liana.al = 0.1 #percentage of making leaf area for liana

#set tree dbh and make liana dbh equals to tree height
tree.dbh = mean.hori.data.tree.dbh
liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)


#This model calculates liana NPP (photosynthesis rate and respiration rate) based on 
#meteorological variables and set parameters from previous literature. 
liana.NPP.mxh.day = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, 
                             b2 = NULL,nday, month,tot.al, frac.liana.al, 
                             liana.length, Vm = NULL){
  
  #calculate or input more the constants needed for the model
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values, these are vectors for each hour
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  #the D here is the where daily and monthly model differ 
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD, D_vec is 1:24 VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    
    rG = 0.3 # Growth respiration, default
    q = 1.89 # Fine root:leaf ratio, default
    #if SLA is defined in the function input use the input, or use the default
    #SLA = specific leaf area
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
    
    #the K here is the constant we chose for this model
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis, k is a just one value here
    
    #same with SLA, if there is no input for these parameters, use the default values
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
    
    #same with SLA, if there is no input for these parameters, use the default values
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
    # equations from Trugman et al. (2018)
    T1 = 1 / (Ca - gamma_star)
    T2 = 1 / (Ca + km)
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    #smoothed VPD are used there to calcualte photosynthesis (A) rate
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
    #Here shows how VPD input impact the out put
    A1 = V * ((xc - T1) / (xc - T2)) #here is when we used VPD
    
    # Multiply each respiration component by 3600 seconds (in 1 hr), now we have all the respiration
    #rd, rx, rp, rr are caculated above
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    #make A based on if there is no sunlight or not
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
  
  # Monthly turnover & NPP changed to daily turnover
  stem.turn = 0.1 / 365 #nday[month] #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24 #this is now daily
  rx_turn = rx_d - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24 #this is now daily
  rp_turn = rp_d - rp_lost
  
  # NPP under non-water-stressed conditions, nday[month] changed to 31, what is this now? changed to hourly?
  NPP_prelim = ((1 - rG) * (A1_d * liana.al - rd_d * (liana.al) - rr_d - rx_turn - rp_turn) * 1e-9 * 12)#made change
  
  # Allows leaves to drop if water stress is present, nday[month] changed to 31
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - rG) * (A1_d * liana.al - rd_d * (liana.al) - rr_d - rx_turn - rp_turn) * 1 * 10^-9 * 12)#made change
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 /365 #made change
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    #rx_lost = rx_lost * nday[month]
    rx_turn = rx_d - rx_lost 
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    #rp_lost = rp_lost * nday[month]
    rp_turn = rp_d - rp_lost
    
    # Calculate NPP, rG is just a fraction, remove things like /nday[month]
    NPP = (1 - rG) * (A1_d * liana.al - rd_d * (liana.al) - rr_d - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

#Use Liana_K$K

summary(lianas_K$K)
length(lianas_K$K)

#July (Day)

nday = c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


sorted.K <- sort(lianas_K$K, decreasing = FALSE)

sorted.K.mPa <- sorted.K/1000

july.npp.day = rep(0,51)

for(j in 1:51){
  month=7
  K=sorted.K.mPa[j]
  for(i in 1:31){
      psis = SWP[month] 
      D = vpd_day[(cumsum(nday[0:(month-1)])[month-1])+i,]
      
      liana.out = liana.NPP.mxh.day(dbh = dbh, psis = psis, D_vec = D, liana.k = K, tot.al = tot.al, nday=nday,month=month,
                                    frac.liana.al = frac.liana.al, liana.length = liana.length)
      july.npp.day[j] = july.npp.day[j]+liana.out[1]
    }
  }

july.npp.day


#July (month)
july.npp.mon = rep(0,51)

for(j in 1:51){
  month=7
  K=sorted.K.mPa[j]
  for(i in 1:31){
    psis = SWP[month] 
    D = vpd_month[month,]
    
    liana.out = liana.NPP.mxh.day(dbh = dbh, psis = psis, D_vec = D, liana.k = K, tot.al = tot.al, nday=nday,month=month,
                                  frac.liana.al = frac.liana.al, liana.length = liana.length)
    july.npp.mon[j] = july.npp.mon[j]+liana.out[1]
  }
}

july.npp.mon

dif = july.npp.day - july.npp.mon
dif.percentage = dif/july.npp.day*100

dif
dif.percentage

summary(july.npp.day)
summary(july.npp.mon)
df<-data.frame(K=sorted.K.mPa, day_NPP=july.npp.day,month_NPP=july.npp.mon,dif=dif.percentage)
plot(df$K,df$day_NPP)
points(df$K,df$month_NPP, col="green",pch=4)
legend("topleft", legend=c("day_NPP", "month_NPP"),
       col=c("black", "green"), lty=1:2, cex=0.8)

#look at the strange day
#when k=680.55556



