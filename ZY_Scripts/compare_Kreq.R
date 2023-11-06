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

#compared monthly and daily NPP over a year and specifically with new VPD calculated from sample data for Jan
load('/Users/zhuoranyu/Desktop/liana-tree-test/Data/param.input.RData')
load('/Users/zhuoranyu/Desktop/liana-tree-test/Data/bci_met_mxh.RData')

########all the constants
#make K and dbh
nK = 75000 
ks = c(trees_K$K / 10, lianas_K$K)
varyk = seq(min(ks), max(ks), length.out = nK) #make a sequence of k

# Interpolate DBH (and remove highest DBH) 
dbhs = seq(liana.dbh.low, liana.dbh.high, length.out = 6)

tot.al = 2 * 100
frac.liana.al = 0.4
tree.dbh = mean.hori.data.tree.dbh
liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)

#######


#Modify Monthly Loop

#monthly NPP function (changed rG only and liana length)
liana.NPP.mxh.mon = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, 
                             b2 = NULL, nday, month,tot.al, frac.liana.al, 
                             liana.length, Vm = NULL){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  Ca = 400 # Atmospheric CO2 concentration (ppm)
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
  
  # Monthly Anet
  Anet = ((1 - rG) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  # Monthly turnover & NPP
  stem.turn = 0.1 / 365 * nday[month] #yearly turnover / months
  
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
  rp_lost = rp_lost * nday[month]#define this input
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - rG) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - rG) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 365 * nday[month]
    
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
    
    # Calculate NPP
    NPP = (1 - rG) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}


lianaouteach.mon = rep(0,12)
liana_out_mon_1 = matrix(, nrow = length(dbhs), ncol = 2)

for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      psis = SWP[j]
      D = VPD[j,]
      dbh = dbhs[i]
      K = varyk[k]
      liana.out = liana.NPP.mxh.mon(dbh = dbh, psis = psis, D_vec = D, liana.k = K, 
                                    nday = nday, month = month, tot.al = tot.al, 
                                    frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach.mon[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach.mon)
    if(sumlianaout > 0){
      liana_out_mon_1[i,1] = sumlianaout
      liana_out_mon_1[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

lianaouteach.mon
sum(lianaouteach.mon)

liana_out_mon_1 <-data.frame(liana_out_mon_1)
colnames(liana_out_mon_1) <-c("DBH","K")
liana_out_mon_1$K_req <- liana_out_mon_1$K/1000

