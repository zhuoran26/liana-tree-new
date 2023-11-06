## This script expands the tree- and liana-specific NPP models modified from Trugman et al. 
## This iteration allows met drivers to enter at different temporal scales

## Because VPD varies substantially at a subdaily timestep, VPD is input hourly
## This means that Anet and NPP are calculated hourly

## Because SWP does not vary substantially over small timesteps, it only varies monthly

## With the new hourly timestep, another assumption is made:
## the day is split between 12 hours of light and 12 hours of darkness

## This can be changed over the course of a year if desired, but here,
## it remains constant, which is fairly consistent with the tropical systems we are considering
## Meanwhile, all respiration components occur 24 hours per day

## Author: AT Trugman, modified by AM Willson
## Date modified: 07 July 2020

## Tree model
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, 
                        b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, 
                        Vm = NULL){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  Ca = 400 # Atmospheric CO2 concentration (ppm)
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
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

  # Montly Anet
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  # Monthly turnover & NPP
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
    
    # Calculate NPP
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}








## Liana model
liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, 
                         b2 = NULL, nday, month, tot.al, frac.liana.al, 
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
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  # Monthly turnover & NPP
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
    
    # Calculate NPP
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}
