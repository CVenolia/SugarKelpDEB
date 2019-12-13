#sensitivity analyses using sensFun (Local sensitivity analyses)
library(FME)
library(tidyverse)
library(lubridate)

source("SolveR_R.R")
##### Minerals and Organics Section #####
#RL: Conversion coefficients #CTV: organics (n = matrix of chemical indices)
#CTV: "food N" "food C" Stucture "N reserves" "C reserves" products
#RL:     X_N   X_C   V     E_N   E_C   P
n_O <- matrix(
  + c(0.00, 1.00, 1.00, 0.00, 1.00, 1.00,  #RL: C/C, equals 1 by definition
      + 1.50, 0.50, 1.33, 3.00, 2.00, 1.80,  #RL: H/C, these values show that we consider dry-mass
      + 1.50, 2.50, 1.00, 0.00, 1.00, 0.50,  #RL: O/C
      + 1.00, 0.00, 0.04, 1.00, 0.00, 0.04), nrow=4, ncol=6, byrow = TRUE) #RL: N/C
#CTV: X_N is the struture of NO3- and NH3 averaged
#CTV: X_C is the average of structure of HCO3- and CO2
#CTV: V is the C-mol structure of alginate (Alginic acid: (C6H8O6)n)
#CTV: E_N is NH3
#CTV: should E_C be glucose C6H12O6? Going to try it out above?
#CTV: E_C average structure of Laminarin: c18h32o16 and mannitol c6h14o6
#CTV: We aren't using the P collumn here, so this information is not modified from Romain's models (not sure where he got the original n_O matrix)

#This part is commented out and left in MATLAB formatting because it is currently unnecessary
#CTV: minerals (n = matrix of chemical indices) (define the chemical environment of the indivdual)
#CTV: Carbon dioxide Water dioxygen 'nitrogenous waste'
#RL:      C     H     O     N
#RL: n_M = [ 1     0     0     0;    #RL: C/C, equals 0 or 1
#0     2     0     3;    #RL: H/C
#2     1     2     0;    #RL: O/C
#0     0     0     1];   #RL: N/C

#RL: Specific densities
#RL: Parameters that link moles to grams (wet weight), volumes and energy
#RL: X_N  X_C   V   E_N  E_C   P     dry mass per wet volume
#d_O <- c(1.0, 1.0, 0.1, 0.1, 0.1, 0.1) #RL: g/cm^3, specific densities for organics

#RL: Chemical potentials
#RL:     X_N   X_C   V    E_N  E_C   P
#mu_O <- c(56.3, 62.4, 500, 56.3, 516, 480) * 1000 #RL: J/mol, chemical potentials for organics
#RL: mu_O= [32.6; 62.6]; Tab 4.3
#RL: XN = NH3(g) entropy at 20 degC
#RL: XC = CO2(g) entropy at 20 degC
#RL: V  = Kooijman 2010
#RL: EN = NH3(g) entropy at 20 degC
#RL: EC = Tab 4.2, page 150
#RL: P  = Kooijman 2010

#RL: Molecular weights
#Simpler code in MATLAB, t() is a matrix transpose function
#CTV: organics structure matrix multiplied by the atomic masses (mass in grams of one mole of an element) of C H O N
w_O_step <- t(n_O)*matrix(c(12, 1, 16, 14), nrow=6, ncol=4, byrow= TRUE) #RL: g/mol, molecular weights for organics
w_O <- rowSums(w_O_step) #CTV: this provides g/mol of each of the six "pockets of mass" (i.e. X_N, X_C)

#R: Pack coefficients #CTV: not a whole lot of point to packing this data as dwm for the current code
#HOW TO STORE 3 Vectors as collumns and them combine them in a matrix
#dwm <- matrix(c(d_O, w_O, mu_O), nrow = 6, ncol=3, byrow = FALSE) #R: g/cm^3, g/mol, kJ/mol specific density, molecular weight, chemical potential
w_V = w_O[3] #w_V = dwm[3,2] # g/mol       #molecular weight of structure
w_EN = w_O[4] #w_EN = dwm[4,2] # g/mol      #molecular weight of N reserve
w_EC = w_O[5] #w_EC = dwm[5,2] #g/mol       #molecular weight of C reserve
w_O2 <- 32 #g/mol
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##### Parameters compiled #####
#Note Romain's Kelp MATLAB model from July 2018 has shadowing coefficients that CTV has not included
#Thinking that shading will be more important if we extrapolate from an individual model to a population model for this scenario
params_Lo <- c(#shape coefficient #left here without proper modification (value from trials with Romain), relationship not linear, so using the allometric relationship might be stronger
  #del_M  = 0.13, #no unit
  #elemental coefficient of nitrogen in structure #ASKING ROMAIN IF THIS VALUE SHOULD BE CHANGED
  #n_NV   = 0.04, #no unit
  n_NV   = 0.03, # RL: from what I saw in the literature it could be a bit less
  
  #maximum volume-specific assimilation rate of N before temperature correction #Espinoza and Chapman 1983
  #JENAM = 167.5, #mugNO3 * h-1 g-1 dry weight
  #JENAM = 145, #167.5
  JENAM = 2.7e-4, #0.0002848824, #mol N / molM_V / h
  #maximum surface-specific assimilation rate of N #Espinoza and Chapman 1983
  #K_N = 4.4/1000000, #molNO3/L (converted from micromoles/L)
  K_N = 2.667e-6, #4e-6, #4.4e-5
  
  #max volume-specific carbon dioxide assimilation rate #Longphuirt et al (2013)
  #JCO2M = 7.69, #24.91, #112, #50 #originally Longphuirt et al (2013): micromole C g^-1 DW h^-1
  JCO2M = 0.0075, #0.025 #Now molC/molM_V/h
  #half saturation constant of C uptake #Davison 1987 (not the most reliable source based on context)
  #K_C = 4/1000000, #mol C/L (of Rubisco)
  K_C = 4/1000000, #1e-3,
  #maximum volume-specific carbon assimilation rate #value set to make C not limiting (ideally this gets changed through fitting)
  #JECAM = 10, #molC/molM_V/h
  JECAM = 0.282, #8.33, #10, #was 200 when per day #ONLY NEEDS STANDARDIZATION to structure (by multiplication of B/M_V) if using info that needs calibration
  
  #Lorena's Wu and Mercheck photosynthesis model kept here just in case
  #alpha = 1.9e-1, #(mol PSU micro E m^-2)^-1 #PSU excitement coefficient
  #beta = 5.8e-5, #(mol PSU micro E m^-2)^-1  #PSU inhibition coefficient
  #gamma = 1.46e-4*3600*24, #1/molPSU/d       #PSU relaxation rate (returning to ground state from excited state)
  #delta = 4.8e-4*3600*24, #1/molPSU/d        #PSU recovery rate (inhibited returning to ground state)
  
  #rho_PSU = 0.1, #mol PSU/ mol Mv             #Photosynthetic unit density
  rho_PSU = 0.05, #0.05, #ONLY NEEDS STANDARDIZATION to structure (by multiplication of B/M_V) if using info that needs calibration
  #p_L =  1.9e-1, #dimensionless #binding probability of photons to a Light SU
  b_I = 2.8e-6, #5e-5 (the value from the defendible copy of the thesis)
  #k_ATP = 1.46e-4*3600*24, #mol NADPH/molPSU/h #dissociation rate of ATPsynthase and ferrodoxin-NADP reductase products
  k_I = 0.28, #0.29 (the value from the defensible copy of the thesis)
  #Yield factor of C reserve to NADPH
  y_I_C = 2, #mol NADPH mol C-1
  #Yield factor of C reserve to CO2
  y_CO2_C = 1, #mol CO2 mol C-1
  
  #reserve turnover #Lorena (2010) #value changed though overall model fitting
  #kE = 2.6,  #1/h,
  kE_C = 0.05, #0.02, #5/24, #7 
  kE_N = 0.01, #0.03, #7/24, #7
  #fraction of rejection flux from growth SU incorporated back into i-reserve #Lorena (2010) picked an arbitrary value here, found low sensitivity to this value
  #kappa_Ei = 0.7, #no unit
  kappa_Ei = 0.9,
  #yield of structure on N reserve (percent of N in structure), Lorena published with a G instead of a V #Romain calls typo #value from Lorena (2010)
  #y_EN_V = 0.04, #mol EN/mol M_V
  y_EN_V = 0.04, #0.03,
  #yield of structure on C reserve (percent of C in structure), Lorena published with a G instead of a V #Romain calls typo #value from Lorena (2010)
  #y_EC_V = 1.25, #mol EC/mol M_V
  y_EC_V = 1,
  #specific maintenance costs requiring N before temp correction #value from Lorena (2010)
  #JENM = 0.012, #molEN/molM_V/h
  JENM = 3.2e-05, #0.00015, # 5.9e-04
  #specific maintenance costs requiring C before temp correction #value from Lorena (2010)
  #JECM = 0.054, #molEC/molM_V/h
  JECM = 1.4e-05, #0.0000225, #0.000034, # 0.0034
  
  #Arrhenius temperature #calculated in Romain's Matlab code
  T_A = 6314.3, #Romain's second Arrhenius relationship: 1562, # K
  #Upper boundary of temperature tolerance
  T_H = 13.386 + 273.15, #Romain's second Arrhenius relationship: 26 + 273.15, # K
  #Lower boundary of temperature tolerance
  T_L = 273.15, # K #fixed at 0 by Romain
  #Arrhenius temperature outside T_H
  T_AH = 18702, #Romain's second Arrhenius relationship: 38577, #K
  #Arrhenius temperature outside T_L
  T_AL = 4391.9, #Romain's second Arrhenius relationship: 15572, #K
  #temperature at which rate parameters are given (set at 20 degrees C)
  T_0 = 20 + 273.15) # K
###### Set up NOAA data ####
#NOAA irradiance data set-up: NOAASurfaceIrradiance
NOAA_Irradiance <- read.csv("NOAASurfaceIrradiance.csv", header = TRUE) #Import
NOAA_Irradiance$DateTime <- dmy_hms(NOAA_Irradiance$DateTime, tz = "UTC") #NOAA data in UTC (5 hours ahead)
NOAA_Irradiance <- with_tz(NOAA_Irradiance, "America/New_York") #Convert from UTC to EST
NOAA_Irradiance$DownMinusUp <- NOAA_Irradiance$dswrf-NOAA_Irradiance$uswrf #net shortwave radiation at the surface (W/m^2) is obtained by subtracting the upward short wave flux (uswrf) from the downward flux (dswrf)
#PAR = NSW*PAR_frac*C*exp(-k*z)*3600
#NSW=dswrf-uswrf
#PAR_frac is the fraction of the incident flux that is useable for photosynthesis (Dave used 0.368 from wikipedia) #RESEARCH
#C is a conversion factor = 4.56 umol photons/s/W (Dave took from Wikipedia)
#k is the extinction coefficient (Dave used 0.46 m^-1 from CTD casts) and 3600 converts from s^-1 to h^-1. 
NOAA_Irradiance$PAR <- NOAA_Irradiance$DownMinusUp*0.43*4.56*exp(-0.46*1)*3600
ggplot() + 
  geom_point(data = NOAA_Irradiance, aes(DateTime, PAR))

#######

SenseFunc <- function(params_Lo){
  rates_Lo <- function(t, state, parameters) { #CTV: note that the time, state, and parameters in this function definition should always remain in this format regardless of how these pieces of info are labelled
    
    with(as.list(c(state, parameters)), { 
      
      #set-up equations
      #2.16 Lorena Temp correction function:
      C_T <- exp((T_A/T_0)-(T_A/T_field(t))) * (1+exp((T_AL/T_0)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_0))) * ((1+exp((T_AL/T_field(t))-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_field(t))))^-1) #C_T has no units
      #C_T <- 1
      
      #change max assimilation rate of N03 from per unit weight to per unit structural mass (mol M_V)
      #WATCH OUT: as JENAM is calibrated check whether or not w_EN is still necessary here!!!!!!!!!
      #CTV caught a conversion error here 7/25/19: divide by w_EN instead of multiply
      #JENAM <- (JENAM*10^-6/w_EN)*(B/M_V) #mugNO3 * h-1 g-1 dry weight to mol N / molM_V / h
      J_EN_AM <-  JENAM * C_T #temp correct max assimilation rate of N
      #2.1 Lorena specific assimilation rate of N
      J_EN_A <- J_EN_AM*(N_field(t)/(N_field(t)+K_N))  #J_EN_A units: mol N/molM_V/h
      N <- N_field(t) #just to save output
      
      #CTV caught a conversion error here 7/25/19: w_EC does not belong here because the C units are already in moles
      #JCO2M <- (JCO2M*10^-6)*(B/M_V) #micromole C g^-1 DW h^-1 converted to molC/molM_V/h
      J_CO2_M <- JCO2M * C_T #temp correct max uptake rate of CO2 (photon capture not influenced by temp) J_CO2_M units: molC/molM_V/d
      #2.4 Lorena #specific CO2 uptake flux
      J_CO2 <- J_CO2_M*(CO_2/(CO_2+K_C)) #J_CO2 units: molC/molM_V/h
      
      #2.2 Lorena #specific relaxation rate
      #J_I <- (rho_PSU*I_field(t))/((1/alpha)+(1/(gamma*(1+(beta/alpha))*I_field(t)))+((beta/(gamma*delta))*I_field(t)^2)) #Units: mol C/mol M_V/d
      #equation Romain adapted from Papadakis (2005)
      J_I <- (rho_PSU * I_field(t) * b_I)/(1+I_field(t) * b_I/(k_I * C_T)) #Unit: molNADPH molM_V–1 h–1
      J_EC_AM <- JECAM * C_T #temp correct max assimilation rate of C
      #2.3 Lorena #specific carbon assimilation flux #aka photosynthesis rate
      J_EC_A <- ((1/J_EC_AM)+(1/(J_CO2/y_CO2_C))+(1/(J_I/y_I_C))-(1/((J_CO2/y_CO2_C)+(J_I/y_I_C))))^-1 #Units: mol C/mol M_V/h
      I <- I_field(t) # RL: just to save output
      
      #Oxygen
      J_O <- (J_I * M_V * w_O2) / (B * 4) # g O2/g/h O2 production (1 molO2 for 4 molC of e-/NADPH produced)
      
      #Temp corrections for the more internal parts of the model
      k_E <- c(kE_C*C_T, kE_N*C_T) #temp correct reserve turnover rate #Unit: 1/h
      J_EN_M <- JENM*C_T #temp correct volume specific maintenance cost paid by the N reserve #unit: molEN/molM_V/h
      J_EC_M <-JECM*C_T #temp correct volume specific maintenance cost paid by the C reserve #unit: molEC/molM_V/h
      
      #########################
      #setting up holder vectors for the SolveR_R function
      #CARBON in first slot, NITROGEN in second
      m_E <- c(m_EC, m_EN)      #density of both reserves (mol Ei/molM_V)
      J_EM <- c(J_EC_M, J_EN_M) #volume specific maintenance cost paid by both reserves (molEi/molM_V/h)
      y_EV <- c(y_EC_V, y_EN_V) #yield of structure on both reserves (mol Ei/mol M_V)
      J_VM_C <- J_EC_M/y_EC_V  #rate of carbon flux of structure that pays for maintenance when the catabolic flux is not enough (1/h)
      J_VM_N <- J_EN_M/y_EN_V #rate of nitrogen flux of structure that pays for maintenance when the catabolic flux is not enough (1/h)
      #J_VM <- J_VM_C + J_VM_N #rate of maintenance costs paid from structure when the catabolic flux is not enough (allows r to be negative) (1/h)
      J_VM <- c(J_VM_C, J_VM_N) # RL: not additive here, they are first used separately in Solve_R routine
      #r0 <- 0
      r0 <- 1 # RL: I would not start the regression procedure with 0 (not sure why but it could have consequences in the procedure)
      #The loop to solve for r (specific growth rate)
      Output_loop <- SolveR_R(m_E, k_E, J_EM, y_EV, J_VM, r0)
      
      #Unpacking SolveR_R output
      #RETURNS r, J_EC_loop (2), JEM_loop (2), JVM_loop (2), JER_loop (2), info
      r <- Output_loop[1] #units should be 1/h??
      J_EC_C <- Output_loop[2] #unit: mol EC/molM_V/h
      J_EN_C <- Output_loop[3] #unit: mol EN/molM_V/h
      J_EC_M <- Output_loop[4] #unit: mol EC/molM_V/h
      J_EN_M <- Output_loop[5] #unit: mol EN/molM_V/h
      J_VM_C <- Output_loop[6] #unit: 1/h
      J_VM_N <- Output_loop[7] #unit: 1/h
      J_VM <- J_VM_C+J_VM_N #unit: 1/h
      J_EC_R <- Output_loop[8] #unit: mol EC/molM_V/h
      J_EN_R <- Output_loop[9] #unit: mol EN/molM_V/h
      info <- Output_loop[10] #1 if convergence, 0 if no convergence (no units)
      
      #Allocation to growth
      #2.8 Lorena #specific growth flux for nitrogen 
      J_EN_G <- J_EN_C-J_EN_M #unit: mol EN/molM_V/h
      #2.8 Lorena #specific growth flux for carbon
      J_EC_G <- J_EC_C-J_EC_M #unit: mol EC/molM_V/h
      #CTV wrote this function to keep from getting NaN for J_G when growth is not occuring
      #Necessary because Inf-Inf is NaN, but Inf-0 is Inf
      RealNumbas <- function(x) {
        if(is.infinite(x) == TRUE){
          x <- 0 #replace infinite numbers with zero
        }
        return(x)
      }
      #Lorena 2.9 #gross specific growth flux
      J_G <- ((((J_EN_G/y_EN_V)^-1)+ ((J_EC_G/y_EC_V)^-1))-RealNumbas((J_EN_G/y_EN_V+J_EC_G/y_EC_V)^-1))^-1 #unit: 1/h
      
      #CTV: testing whether or not the two different ways of solving J_G are equivalent
      J_G_test <- r + J_VM #unit: 1/h
      JER_C_test <- max(0 , J_EC_C-J_EC_M-y_EV[1]*(J_G)) #units: mol Ei/molM_V/h
      JER_N_test <- max(0, J_EN_C-J_EN_M-y_EV[2]*(J_G))
      
      #Biomass (structure mass+ C reserve mass + N reserve mass)
      B <- (w_V+m_EN*w_EN+m_EC*w_EC)*M_V #unit: g (as long as the units of M_V are mol M_V)
      #N to C ratio Lorena eq 3.1 (simiplified from (n_NV*M_V+M_EN)/(M_V+M_EC)) #NOTE THE DIFFERNCE BETWEEN M_EN (total mass of N reserve) and m_EN (reserve density)!!
      #NC <- (n_NV+m_EN)/(1+m_EC) #units: mol N/mol C
      #structural mass in grams
      #V <- M_V/d_V #units: cm^3 #would need to find d_V (structural density per cubic cm)
      #DEB classic way of calculating length from mass using a shape coefficient
      #L_w <- (V^(1/3))/del_M  #CTV: length unit depends on what we use to set the shape coefficient, because the shape coefficient is volume^(1/3)/length
      #Allometic relationship between length (cm) and dry weight (g) from Gevaert (2001)
      L_allometric <- (B/0.00387)^(1/1.469) #unit: cm
      
      #rates
      #dynamics of the C reserve
      dm_ECdt <- J_EC_A-J_EC_C+(kappa_Ei*J_EC_R)-(r*m_EC)  #unit: mol C/mol M_V/h
      #dynamics of the N reserve
      dm_ENdt <- J_EN_A-J_EN_C+(kappa_Ei*J_EN_R)-(r*m_EN) #unit: mol N/mol M_V/h
      #dynamics of structural mass
      dM_Vdt <- r*M_V  #units of M_V are mol M_V, so dMVdt units are mol MV/h
      
      #tests added by CTV
      #print(r)
      #print(B)
      
      #output
      #m_EC, m_EN, M_V, J_CO2, J_I, J_EC_A, J_EN_A, J_EC_C, J_EN_C, r, J_EN_M, J_EC_M, J_EN_G, J_EC_G, J_G, J_VM, J_EC_R, J_EN_R, N, CO_2, I, B, NC, V, L_w
      return(list(c(dm_ECdt, dm_ENdt, dM_Vdt), r=r, B=B, C_T=C_T, I=I, N=N, J_CO2=J_CO2, J_I=J_I, J_EC_A=J_EC_A, J_EN_A=J_EN_A, J_EC_C=J_EC_C, J_EN_C=J_EN_C, J_EN_M=J_EN_M, J_EC_M=J_EC_M, J_EN_G=J_EN_G, J_EC_G=J_EC_G, J_G=J_G, J_VM=J_VM, J_EC_R=J_EC_R, JER_C_test=JER_C_test, J_EN_R=J_EN_R, JER_N_test=JER_N_test, CO_2=CO_2, L_allometric=L_allometric, J_G_test=J_G_test, J_O=J_O, info=info))
    })
  }
  ####### State Inititial conditions ############
  #Initial conditions of state variables
  #these values are not coming from any field data or literature information (just estimation by CTV)
  state_Lo <- c(m_EC = 0.1, #mol EC/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
                m_EN = 0.03, #mol EN/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
                M_V = 0.00009) #molM_V        #initial mass of structure
  #try for initial biomass around 0.0035
  #(w_V + 0.03 * w_EN + 0.1 * w_EC) * 0.00009
  #cannot have these as initial conditions in the desolve formatting without the program thinking they are differential equations
  #r = 0, # 1/d,            #initial growth rate
  #B = ((29.89 +0.03 * 17 + 0.1 * 30) * 1) # 1.4853 g starting mass (w_V + m_EN * w_EN + m_EC * w_EC) * M_V)
  #state_Lo_Dec <- c(m_EC = 0.1, #mol EC/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
  #m_EN = 0.1, #mol EN/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
  #M_V = 0.01) #molM_V        #initial mass of structure
  # (w_V + m_EN * w_EN + m_EC * w_EC) * M_V)
  
  state_EspinozaChapman1983_Growth <- c(m_EC = 0.01, #cut to 2 cm length at start B = 0.326g
                                        m_EN = 0.001, #paper only gives length Celeste made up the rest of the numbers based around biomass calculated through the allometric relationship
                                        M_V = 0.015)
  #B = ((29.33 +0.001 * 17 + 0.01 * 29.1) * 0.015)
  
  state_Sjotun1993 <- c(m_EC = .1, #mol EC/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
                        m_EN = 0.3, #mol EN/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
                        M_V = 0.125) #molM_V        #initial mass of structure
  #esting out 4.175 = (w_V + 0.03 * w_EN + 0.1 * w_EC) * 0.125
  #technically should be B = 496.5 = ((29.33 + m_EN * 17 + m_EC * 29.1) * M_V)
  #for Olischlager 2017 doesn't give a starting size
  cond <- c(m_EC = 0.01, #cut to 2 cm length at start B = 0.326g
            m_EN = 0.001, 
            M_V = 0.015)
  #B = ((29.33 +0.001 * 17 + 0.01 * 29.1) * 0.015)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #######Time step to run the model on#######
  #(First number of time step, last number of time step, interval to step)
  times_Lo_Sled1 <- seq(0, 4008, 1) #167 days
  times_Lo_Sled2 <- seq(0, 3336, 1) #139 days
  times_Lo_Dredge1 <- seq(0, 4128, 1) #172 days
  times_Lo_Dredge2 <- seq(0, 3456, 1) #144 days
  times_Lo_Wickford1 <- seq(0, 3312, 1) #138 days
  times_Lo_RomePt1 <- seq(0, 4104, 1) #171 days
  times_Lo_RomePt2 <- seq(0, 3264, 1) #136 days
  
  #Sled Sugar Kelp Growth from lines planted on 12/12/18 and 2/6/19
  #Dredge Sugar Kelp Growth from lines planted on 12/12/18 and 2/6/19
  #Wickford Sugar Kelp Growth from line planted on 12/19/18
  #Rome Point Sugar Kelp Growth from lines planted on 12/20/18 and 2/21/19
  times_Y2_Sled1 <- seq(0, 3408, 1) #142 days #end 5/3/19 #19+31+28+31+30+3
  times_Y2_Sled2 <- seq(0, 2064, 1) #86 days #end 5/3/19 #22+31+30+3
  times_Y2_Dredge1 <- seq(0, 3408, 1) #142 days #end 5/3/19 #19+31+28+31+30+3
  times_Y2_Dredge2 <- seq(0, 2064, 1) #86 days #end 5/3/19 #22+31+30+3
  times_Y2_Wickford1 <- seq(0, 3720, 1) #155 days #end 5/23/19 #12+31+28+31+30+23
  times_Y2_RomePt1 <- seq(0, 3720, 1) #155 days #end 5/24/19 #11+31+28+31+30+24
  times_Y2_RomePt2 <- seq(0, 2208, 1) #92 days #end 5/24/19 #7+31+30+24
  
  
  times_Sjotun1993 <- seq(0, 9408, 1) #392 days
  #########
  return(as.data.frame(ode(y = state_Lo, times = times_Lo_Sled1, func = rates_Lo, parms = params_Lo)))
}

#Sled1
B <- 0.003006
###### N forcing set-up##############
#GSO_N <- read.csv("GSO_Ndata2016_2017.csv", header = TRUE)
#GSO_N$Date <- mdy(GSO_N$Date)
#GSO_N$N_combined <- GSO_N$NO3_NO2+GSO_N$NH4
#GSO_N$N_combined <- GSO_N$N_combined/1000000

#ggplot() + 
#geom_point(data = GSO_N, aes(Date, N_combined), color = "red") +
#geom_point(data = N, aes(Date, combined), color = "black")

WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE) #Import water Q data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
#CTV: I see no change created in the dataset by Romain adding the next line (so it's there for the purpose of his computer being able to run this code)
names(WSA2_Y1)[1] <- "Site" # RL: I had to use this because the first column name was different (probably due to some encoding glitch)
Sled_WSA <- filter(WSA2_Y1, Site == "Sled")
daily <- seq(as.Date("2017-11-1"), as.Date("2018-04-17"), by="days")
N <- Sled_WSA[c("Date","NitrateNitrite_uM", "Ammonia_uM")]
N$combined <- N$NitrateNitrite_uM + N$Ammonia_uM
N$combined <- N$combined/1000000 #convert from micromoles/L to moles/L
#Converted to hourly!
N_field <- approxfun(x = c(161*24, 139*24, 105*24, 0, 28*24, 84*24, 172*24), y = N$combined, method = "linear", rule = 2) #N forcing function

#### OLDER CODE 
#N <- mean(Nitrate$Nitrate_uM_calculatedFromDifference) #molNO3/L #CTV: kept here as an alternative to using splines
#N_interpolation <- data.frame(Date=daily, N_interp=spline(N, method="fmm", xout=daily)$y) #CTV: spline interpolation kept here as an alternative
#N_interpolation$N_interp[N_interpolation$N_interp<7.653510e-07]<- 7.653510e-07
#N_matrix <- as.matrix(cbind(0:167, N_interpolation$N_interp))
#ggplot() + 
#geom_point(data = N_interpolation, aes(Date, N_interp), color = "red") +
#geom_point(data = N, aes(Date, combined), color = "black")
#7 actual data points #need to make 161 approximated winter datapoints

#TESTING  FAKE NUMBER!!
#N_field <- approxfun(x = c(0:4008), y = rep(15/1000000, 4009), method = "linear", rule = 2) #N forcing function
###### C forcing set-up ###########
Sled_DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE) #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

###### NOAA Irradiance forcing set-up ####
#as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00")
NOAA_Irradiance_Sledy1 <-  NOAA_Irradiance$PAR[2438:3774]
I_field <- approxfun(x = seq(from = 0, to = 4008, by = 3), y = NOAA_Irradiance_Sledy1, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
#Import Sled Hobo (cynlinder) of just temp: This is the better temp data
Sled_Y1_hobotemp <- read.csv("Sled_Y1_TempLogger2.csv", header = TRUE)
Sled_Y1_hobotemp$DateTime <- mdy_hms(Sled_Y1_hobotemp$Date_Time) #convert time field
Sled_Y1_hobotemp <- Sled_Y1_hobotemp[14:16049,]
Sled_Y1_hobotemp$Temp_K <- Sled_Y1_hobotemp$Temp_C+273.15 #create collumn with temp in K
SledT_hourly <- ceiling_date(Sled_Y1_hobotemp$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(Sled_Y1_hobotemp$Temp_K, by=list(SledT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_part1 <- AvgTempKbyhr$x[0:334]
Dredge_Y1_hobo <- read.csv("Dredge_Y1_hobo.csv", header = TRUE) 
Dredge_Y1_hobo$DateTime <- mdy_hms(Dredge_Y1_hobo$Date_Time) #convert time field
Dredge_Y1_hobo <- Dredge_Y1_hobo[3:16531,] 
Dredge_Y1_hobo$Temp_K <- Dredge_Y1_hobo$Temp_C+273.15 #create collumn with temp in K
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour")
AvgTempKbyhr4FD <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr4FD <- AvgTempKbyhr4FD[4:4132, ]
fd <- AvgTempKbyhr4FD$x[335:859] #526 data points needed from dredge to replace a weird glitch in the sled temp data
AvgTempKbyhr_part2 <- AvgTempKbyhr$x[860:4009]
T_field <- approxfun(x = c(0:4008), y = c(AvgTempKbyhr_part1, fd, AvgTempKbyhr_part2), method = "linear", rule = 2) #the temp forcing function
T_Sled1_Y1 <- T_field(0:4008)
#########
Sled1_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(Sled1_SA))
############

#Sled2
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE) #Import water Q data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" # RL: I had to use this because the first column name was different (probably due to some encoding glitch)
Sled_WSA <- filter(WSA2_Y1, Site == "Sled")
daily <- seq(as.Date("2017-11-29"), as.Date("2018-04-17"), by="days")
N <- Sled_WSA[c("Date","NitrateNitrite_uM", "Ammonia_uM")]
N$combined <- N$NitrateNitrite_uM + N$Ammonia_uM
N$combined <- N$combined/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(133*24, 111*24, 77*24, -28*24, 0, 56*24, 144*24), y = N$combined, method = "linear", rule = 2) #N forcing function
#N <- mean(Nitrate$Nitrate_uM_calculatedFromDifference) #CTV: mean kept her as an alternative to spline and linear interpolation
#N_interpolation <- data.frame(Date=daily, N_interp=spline(Nitrate, method="fmm", xout=daily)$y)
#N_interpolation$N_interp[N_interpolation$N_interp<0]<- 0.5
#ggplot() + 
#geom_point(data = N_interpolation, aes(Date, N_interp), color = "red") +
#geom_point(data = Nitrate, aes(Date, Nitrate_uM_calculatedFromDifference), color = "black")
#7 actual data points #need to make 161 approximated winter datapoints
#N_matrix <- as.matrix(cbind(0:139, N_interpolation$N_interp))
#N_field <- approxfun(x = N_matrix[,1], y = N_matrix[,2], method = "linear", rule = 2) #Nitrate forcing function

###### TIC forcing set-up ###########
Sled_DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE) #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

###### NOAA Irradiance forcing set-up ####
#hourly <- seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
NOAA_Irradiance_Sledy1_L2 <-  NOAA_Irradiance$PAR[2662:3774]
I_field <- approxfun(x = seq(from = 0, to = 3336, by = 3), y = NOAA_Irradiance_Sledy1_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
Sled_Y1_hobotemp <- read.csv("Sled_Y1_TempLogger2.csv", header = TRUE)
Sled_Y1_hobotemp$DateTime <- mdy_hms(Sled_Y1_hobotemp$Date_Time) #convert time field
Sled_Y1_hobotemp <- Sled_Y1_hobotemp[6:16051,]
Sled_Y1_hobotemp$Temp_K <- Sled_Y1_hobotemp$Temp_C+273.15 #create collumn with temp in K
SledT_hourly <- ceiling_date(Sled_Y1_hobotemp$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(Sled_Y1_hobotemp$Temp_K, by=list(SledT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[677:4011,]
Dredge_Y1_hobo <- read.csv("Dredge_Y1_hobo.csv", header = TRUE) 
Dredge_Y1_hobo$DateTime <- mdy_hms(Dredge_Y1_hobo$Date_Time) #convert time field
Dredge_Y1_hobo <- Dredge_Y1_hobo[3:16531,] 
Dredge_Y1_hobo$Temp_K <- Dredge_Y1_hobo$Temp_C+273.15 #create collumn with temp in K
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour")
AvgTempKbyhr4FD <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr4FD <- AvgTempKbyhr4FD[4:4132, ]
fd <- AvgTempKbyhr4FD$x[858:859] #526 data points needed from dredge to replace a weird glitch in the sled temp data
T_field <- approxfun(x = c(0:3336), y = c(fd, AvgTempKbyhr$x), method = "linear", rule = 2) #the temp forcing function
T_Sled2_Y1 <- T_field(0:3336)
#############
Sled2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(Sled2_SA))
#################

#Dredge1
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE) #Import water Q data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" # RL: I had to use this because the first column name was different (probably due to some encoding glitch)
Dredge_WSA <- filter(WSA2_Y1, Site == "Dredge")
daily <- seq(as.Date("2017-11-1"), as.Date("2018-04-22"), by="days")
N <- Dredge_WSA[c("Date","NitrateNitrite_uM", "Ammonia_uM")]
N$combined <- N$NitrateNitrite_uM + N$Ammonia_uM
N$combined <- N$combined/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(139*24, 161*24, 84*24, 0, 105*24), y = N$combined, method = "linear", rule = 2) #N forcing function
#Nitrate <- Dredge_WSA[c("Date","Nitrate_uM_calculatedFromDifference")]
#Nitrate$Nitrate_uM_calculatedFromDifference <- Nitrate$Nitrate_uM_calculatedFromDifference/1000000 #convert from micromoles/L to moles/L
#N <- mean(Nitrate$Nitrate_uM_calculatedFromDifference)
#N <- mean(Nitrate$Nitrate_uM_calculatedFromDifference) #CTV: mean kept her as an alternative to spline and linear interpolation
#N_interpolation <- data.frame(Date=daily, N_interp=spline(Nitrate, method="fmm", xout=daily)$y)
#N_interpolation$N_interp[N_interpolation$N_interp<0]<- 1.5
#N_interpolation$N_interp[N_interpolation$N_interp>10]<- 10
#ggplot() + 
#geom_point(data = N_interpolation, aes(Date, N_interp), color = "red") +
#geom_point(data = Nitrate, aes(Date, Nitrate_uM_calculatedFromDifference), color = "black")
#5 actual data points #need to make 173 approximated winter datapoints
#N_matrix <- as.matrix(cbind(0:172, N_interpolation$N_interp))
#N_field <- approxfun(x = N_matrix[,1], y = N_matrix[,2], method = "linear", rule = 2) #Nitrate forcing function

###### DIC forcing set-up ###########
DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE) #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

###### NOAA Irradiance forcing set-up ####
#11/1/17 to 2018-04-22 12:00:00
NOAA_Irradiance_Dredgey1 <-  NOAA_Irradiance$PAR[2438:3814]
I_field <- approxfun(x = seq(from = 0, to = 4128, by = 3), y = NOAA_Irradiance_Dredgey1, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
Dredge_Y1_hobo <- read.csv("Dredge_Y1_hobo.csv", header = TRUE) 
Dredge_Y1_hobo$DateTime <- mdy_hms(Dredge_Y1_hobo$Date_Time) #convert time field
Dredge_Y1_hobo <- Dredge_Y1_hobo[3:16531,] 
Dredge_Y1_hobo$Temp_K <- Dredge_Y1_hobo$Temp_C+273.15 #create collumn with temp in K
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[4:4132, ]
T_field <- approxfun(x = c(0:4128), y = AvgTempKbyhr$x, method = "linear", rule = 2) #the temp forcing function
T_Dredge1_Y1 <- T_field(0:4128)
##############
Dredge1_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(Dredge1_SA))
##############

#Dredge2
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE) #Import water Q data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" # RL: I had to use this because the first column name was different (probably due to some encoding glitch)
Dredge_WSA <- filter(WSA2_Y1, Site == "Dredge")
daily <- seq(as.Date("2017-11-29"), as.Date("2018-04-22"), by="days")
N <- Dredge_WSA[c("Date","NitrateNitrite_uM", "Ammonia_uM")]
N$combined <- N$NitrateNitrite_uM + N$Ammonia_uM
N$combined <- N$combined/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(111*24, 133*24, 56*24, -28*24, 77*24), y = N$combined, method = "linear", rule = 2) #N forcing function
#Nitrate <- Dredge_WSA[c("Date","Nitrate_uM_calculatedFromDifference")]
#Nitrate$Nitrate_uM_calculatedFromDifference <- Nitrate$Nitrate_uM_calculatedFromDifference/1000000 #convert from micromoles/L to moles/L
#N <- mean(Nitrate$Nitrate_uM_calculatedFromDifference) #CTV: kept here as an alternative to linear and spline interpolation
#N_interpolation <- data.frame(Date=daily, N_interp=spline(Nitrate, method="fmm", xout=daily)$y)
#N_interpolation$N_interp[N_interpolation$N_interp<0]<- 1.5
#N_interpolation$N_interp[N_interpolation$N_interp>10]<- 10
#ggplot() + 
#geom_point(data = N_interpolation, aes(Date, N_interp), color = "red") +
#geom_point(data = Nitrate, aes(Date, Nitrate_uM_calculatedFromDifference), color = "black")
#5 actual data points #need to make 173 approximated winter datapoints
#N_matrix <- as.matrix(cbind(0:144, N_interpolation$N_interp))
#N_field <- approxfun(x = N_matrix[,1], y = N_matrix[,2], method = "linear", rule = 2) #Nitrate forcing function

###### TIC forcing set-up ###########
DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE) #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

###### NOAA Irradiance forcing set-up ####
#hourly <- seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
NOAA_Irradiance_Dredgey1_L2 <-  NOAA_Irradiance$PAR[2662:3814]
I_field <- approxfun(x = seq(from = 0, to = 3456, by = 3), y = NOAA_Irradiance_Dredgey1_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
Dredge_Y1_hobo <- read.csv("Dredge_Y1_hobo.csv", header = TRUE) 
Dredge_Y1_hobo$DateTime <- mdy_hms(Dredge_Y1_hobo$Date_Time) #convert time field
Dredge_Y1_hobo <- Dredge_Y1_hobo[3:16531,] #this is when I extimate the hobo gets in the water (cut 2 points in beginning)
Dredge_Y1_hobo$Temp_K <- Dredge_Y1_hobo$Temp_C+273.15 #create collumn with temp in K
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_sub <- AvgTempKbyhr[676:4132,]
T_field <- approxfun(x = c(0:3456), y = c(AvgTempKbyhr_sub$x), method = "linear", rule = 2) #the temp forcing function
T_Dredge2_Y1 <- T_field(0:3456)
##############
Dredge2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(Dredge2_SA))
##############

#Wickford1
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
#CHRP Nitrate data for comparison
#CHRP_Nitrate <- read.csv("CHRP_Nitrate.csv", header = TRUE)
#CHRP_Nitrate$Date <- mdy(CHRP_Nitrate$Date)
#CHRP_Nitrate_1 <- filter(CHRP_Nitrate, Station == "1")
#CHRP_Nitrate_1$NO3 <- as.numeric(as.character(CHRP_Nitrate_1$NO3))
#ggplot() + 
#geom_point(data = CHRP_Nitrate_1, aes(Date, NO3), color = "blue") +
#geom_point(data = N_interpolation, aes(Date, N_interp), color = "red") +
#geom_point(data = Nitrate, aes(Date, Nitrate_uM_calculatedFromDifference), color = "black")
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE) #Import water Q data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" # RL: I had to use this because the first column name was different (probably due to some encoding glitch)
Wickford_WSA <- filter(WSA2_Y1, Site == "Wickford")
daily <- seq(as.Date("2017-12-4"), as.Date("2018-04-21"), by="days")
N <- Wickford_WSA[c("Date","NitrateNitrite_uM", "Ammonia_uM")]
N$combined <- N$NitrateNitrite_uM + N$Ammonia_uM
N$combined <- N$combined/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(115*24, 0, 81*24, 59*24, 38*24, 138*24), y = c(N$combined[1], N$combined[3:7]), method = "linear", rule = 2) #N forcing function
#Nitrate <- Wickford_WSA[c("Date","Nitrate_uM_calculatedFromDifference")]
#Nitrate$Nitrate_uM_calculatedFromDifference <- Nitrate$Nitrate_uM_calculatedFromDifference/1000000 #convert from micromoles/L to moles/L
#N <- mean(Nitrate$Nitrate_uM_calculatedFromDifference) #CTV: mean kept here as an alternative to spline and linear interpolation
#N_interpolation <- data.frame(Date=daily, N_interp=spline(Nitrate, method="fmm", xout=daily)$y)
#N_interpolation$N_interp[N_interpolation$N_interp<0.2]<- 0.3
#ggplot() + 
#geom_point(data = N_interpolation, aes(Date, N_interp), color = "red") +
#geom_point(data = Nitrate, aes(Date, Nitrate_uM_calculatedFromDifference), color = "black")
#5 actual data points #need to make 132 approximated winter datapoints (more data points outside time range)
#N_matrix <- as.matrix(cbind(0:138, N_interpolation$N_interp))
#N_field <- approxfun(x = N_matrix[,1], y = N_matrix[,2], method = "linear", rule = 2) #Nitrate forcing function

###### TIC forcing set-up ###########
#BrentonPoint_Segarra2002CarbonData
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE) #Import lit TCO2 data
names(Segarra2002Carbon)[1] <- "Date" # RL: had to add this line (probably due to some encoding glitch in loaded file)
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date)
mean(Segarra2002Carbon$TCO2_micromolPERkg)
#ggplot() + 
#geom_point(data = Segarra2002Carbon, aes(Date, TCO2_micromolPERkg))
CO_2 <- 1956.143/1000000 #(mol CO2/L)

###### NOAA Irradiance forcing set-up ####
#hourly <- seq(as_datetime("2017-12-4 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
NOAA_Irradiance_Wickfordy1 <-  NOAA_Irradiance$PAR[2702:3806]
I_field <- approxfun(x = seq(from = 0, to = 3312, by = 3), y = NOAA_Irradiance_Wickfordy1, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
#Import Wickford Hobo data
#the overall shape of the two Wickford temp files in practically identical, so I'm using Wickford_Y1_hobo.csv and not Wickford_Y1_hobotemp.csv
Wickford_Y1_hobo <- read.csv("Wickford_Y1_hobo.csv", header = TRUE) 
Wickford_Y1_hobo$DateTime <- mdy_hms(Wickford_Y1_hobo$DateTime) #convert time field
#2017-12-04 15:30:00 first real temp measurement
#2018-04-21 11:30:00 last real temp measurement
Wickford_Y1_hobo <- Wickford_Y1_hobo[607:13839,] #subset based on the above date/times
Wickford_Y1_hobo$Temp_K <- Wickford_Y1_hobo$Temp_C+273.15 #create collumn with temp in K
WickfordT_hourly <- ceiling_date(Wickford_Y1_hobo$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(Wickford_Y1_hobo$Temp_K, by=list(WickfordT_hourly), mean) #calculate average hourly temp
fd <- AvgTempKbyhr[1:4,]
T_field <- approxfun(x = c(0:3312), y = c(fd$x, AvgTempKbyhr$x), method = "linear", rule = 2) #the temp forcing function
T_Wickford1_Y1 <- T_field(0:3312)
##############
Wickford1_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(Wickford1_SA))
##############

#RomePt1
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE) #Import water Q data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" # RL: I had to use this because the first column name was different (probably due to some encoding glitch)
RomePt_WSA <- filter(WSA2_Y1, Site == "Rome Point")
daily <- seq(as.Date("2017-11-1"), as.Date("2018-04-21"), by="days")
N <- RomePt_WSA[c("Date","NitrateNitrite_uM", "Ammonia_uM")]
N$combined <- N$NitrateNitrite_uM + N$Ammonia_uM
N$combined <- N$combined/1000000 #convert from micromoles/L to moles/Lmean
N_field <- approxfun(x = c(114*24, 148*24, 0, 71*24, 92*24, 171*24), y = c(N$combined[1:3], N$combined[5:7]), method = "linear", rule = 2) #N forcing function
#Nitrate <- RomePt_WSA[c("Date","Nitrate_uM_calculatedFromDifference")]
#Nitrate$Nitrate_uM_calculatedFromDifference <- Nitrate$Nitrate_uM_calculatedFromDifference/1000000 #convert from micromoles/L to moles/L
#N <- mean(Nitrate$Nitrate_uM_calculatedFromDifference) #CTV: mean kept here as an alternative to spline and linear interpolation
#N_interpolation <- data.frame(Date=daily, N_interp=spline(Nitrate, method="fmm", xout=daily)$y)
#N_interpolation$N_interp[N_interpolation$N_interp<0.2]<- 0.5
#ggplot() + 
#geom_point(data = N_interpolation, aes(Date, N_interp), color = "red") +
#geom_point(data = Nitrate, aes(Date, Nitrate_uM_calculatedFromDifference), color = "black")
#5 actual data points #need to make 173 approximated winter datapoints (more data points outside time range)
#N_matrix <- as.matrix(cbind(0:171, N_interpolation$N_interp))
#N_field <- approxfun(x = N_matrix[,1], y = N_matrix[,2], method = "linear", rule = 2) #Nitrate forcing function


###### TIC forcing set-up ###########
#BrentonPoint_Segarra2002CarbonData
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE) #Import lit TCO2 data
names(Segarra2002Carbon)[1] <- "Date" # RL: had to add this line (probably due to some encoding glitch in loaded file)
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date)
mean(Segarra2002Carbon$TCO2_micromolPERkg)
#ggplot() + 
#geom_point(data = Segarra2002Carbon, aes(Date, TCO2_micromolPERkg))
CO_2 <- 1956.143/1000000 #(mol CO2/L)

###### NOAA Irradiance forcing set-up ####
#hourly <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
NOAA_Irradiance_RomePty1 <-  NOAA_Irradiance$PAR[2438:3806]
I_field <- approxfun(x = seq(from = 0, to = 4104, by = 3), y = NOAA_Irradiance_RomePty1, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
#RomePoint_Y1_hobotemp
RomePoint_Y1_hobotemp <- read.csv("RomePoint_Y1_hobotemp.csv", header = TRUE)
RomePoint_Y1_hobotemp$DateTime <- mdy_hms(RomePoint_Y1_hobotemp$DateTime) #convert time field
#2017-11-01 13:15:00 start
#2018-04-21 14:00:00 end
RomePoint_Y1_hobotemp <- RomePoint_Y1_hobotemp[6:16425,] #subset based on the above date/times
RomePoint_Y1_hobotemp$Temp_K <- RomePoint_Y1_hobotemp$Temp_C+273.15 #create collumn with temp in K
RomePointT_hourly <- ceiling_date(RomePoint_Y1_hobotemp$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(RomePoint_Y1_hobotemp$Temp_K, by=list(RomePointT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[1:4103,]
fd <- AvgTempKbyhr[1:2,]
T_field <- approxfun(x = c(0:4104), y = c(fd$x, AvgTempKbyhr$x), method = "linear", rule = 2) #the temp forcing function
T_RomePt1_Y1 <- T_field(0:4104)
##############
RomePt1_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(RomePt1_SA))
##############

#RomePt2
B <- 0.003006 #1.6545 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE) #Import water Q data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" # RL: I had to use this because the first column name was different (probably due to some encoding glitch)
RomePt_WSA <- filter(WSA2_Y1, Site == "Rome Point")
daily <- seq(as.Date("2017-12-6"), as.Date("2018-04-21"), by="days")
N <- RomePt_WSA[c("Date","NitrateNitrite_uM", "Ammonia_uM")]
N$combined <- N$NitrateNitrite_uM + N$Ammonia_uM
N$combined <- N$combined/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(79*24, 113*24, -25*24, 57*24, 136*24), y = c(N$combined[1:2], N$combined[5:7]), method = "linear", rule = 2) #N forcing function
#Nitrate <- RomePt_WSA[c("Date","Nitrate_uM_calculatedFromDifference")]
#Nitrate$Nitrate_uM_calculatedFromDifference <- Nitrate$Nitrate_uM_calculatedFromDifference/1000000 #convert from micromoles/L to moles/L
#N <- mean(Nitrate$Nitrate_uM_calculatedFromDifference) #CTV: mean kept here as an alternative to spline and linear interpolation
#N_interpolation <- data.frame(Date=daily, N_interp=spline(Nitrate, method="fmm", xout=daily)$y)
#N_interpolation$N_interp[N_interpolation$N_interp<0.2]<- 0.5
#ggplot() + 
#geom_point(data = N_interpolation, aes(Date, N_interp), color = "red") +
#geom_point(data = Nitrate, aes(Date, Nitrate_uM_calculatedFromDifference), color = "black")
#5 actual data points #need to make 132 approximated winter datapoints (more data points outside time range)
#N_matrix <- as.matrix(cbind(0:136, N_interpolation$N_interp))
#N_field <- approxfun(x = N_matrix[,1], y = N_matrix[,2], method = "linear", rule = 2) #Nitrate forcing function

###### TIC forcing set-up ###########
#BrentonPoint_Segarra2002CarbonData
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE) #Import lit TCO2 data
names(Segarra2002Carbon)[1] <- "Date" # RL: had to add this line (probably due to some encoding glitch in loaded file)
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date)
mean(Segarra2002Carbon$TCO2_micromolPERkg)
#ggplot() + 
#geom_point(data = Segarra2002Carbon, aes(Date, TCO2_micromolPERkg))
CO_2 <- 1956.143/1000000 #(mol CO2/L)

###### NOAA Irradiance forcing set-up ####
#hourly <- seq(as_datetime("2017-12-6 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
NOAA_Irradiance_RomePty1_L2 <-  NOAA_Irradiance$PAR[2718:3806]
I_field <- approxfun(x = seq(from = 0, to = 3264, by = 3), y = NOAA_Irradiance_RomePty1_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
RomePoint_Y1_hobotemp <- read.csv("RomePoint_Y1_hobotemp.csv", header = TRUE)
RomePoint_Y1_hobotemp$DateTime <- mdy_hms(RomePoint_Y1_hobotemp$DateTime) #convert time field
#2017-12-06 00:00:00 start
#2018-04-21 14:00:00 end
RomePoint_Y1_hobotemp <- RomePoint_Y1_hobotemp[3313:16425,] #subset based on the above date/times
RomePoint_Y1_hobotemp$Temp_K <- RomePoint_Y1_hobotemp$Temp_C+273.15 #create collumn with temp in K
RomePointT_hourly <- ceiling_date(RomePoint_Y1_hobotemp$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(RomePoint_Y1_hobotemp$Temp_K, by=list(RomePointT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[13:3277, ]
T_field <- approxfun(x = c(0:3264), y = c(AvgTempKbyhr$x), method = "linear", rule = 2) #the temp forcing function
T_RomePt2_Y1 <- T_field(0:3264)
##############
RomePt2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(RomePt2_SA))
##############

#Sled1_Y2
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE) #Import water Q data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #for RL
Sled_WSA <- filter(WSA_Y2, Site == "Moonstone Sled")
daily <- seq(as.Date("2018-12-12"), as.Date("2019-05-03"), by="days")
key <- as.data.frame(daily)
Sled_WSA$Ncombined <- Sled_WSA$NO3NO2_µM + Sled_WSA$NH4_µM
Sled_WSA$Nmol_combined <- Sled_WSA$Ncombined/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(1*24, 57*24, 93*24, 124*24, 142*24, 163*24), y = c(Sled_WSA$Nmol_combined), method = "linear", rule = 2) #N forcing function

###### C forcing set-up ###########
Sled_DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE) #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

###### NOAA Irradiance forcing set-up ####
#hourly <- seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
NOAA_Irradiance_Sledy2 <-  NOAA_Irradiance$PAR[5686:6822]
I_field <- approxfun(x = seq(from = 0, to = 3408, by = 3), y = NOAA_Irradiance_Sledy2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
Sled_Y2_Hobo <- read.csv("Sled_Y2_HoboLightTemp.csv", header = TRUE)
Sled_Y2_Hobo$DateTime <- mdy_hms(Sled_Y2_Hobo$DateTime)
Sled_Y2_Hobo$Temp_K <- Sled_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
SledY2T_hourly <- ceiling_date(Sled_Y2_Hobo$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(Sled_Y2_Hobo$Temp_K, by=list(SledY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[2:3385,]
fd <- rep(285, 25)
T_field <- approxfun(x = c(0:3408), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Sled1_Y2 <- T_field(0:3408)
##############
Sled1_Y2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(Sled1_Y2_SA))
##############

#Sled2_Y2
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE) #Import water Q data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #for RL
Sled_WSA <- filter(WSA_Y2, Site == "Moonstone Sled")
daily <- seq(as.Date("2019-02-06"), as.Date("2019-05-03"), by="days")
key <- as.data.frame(daily)
Sled_WSA$Ncombined <- Sled_WSA$NO3NO2_µM + Sled_WSA$NH4_µM
Sled_WSA$Nmol_combined <- Sled_WSA$Ncombined/1000000 #convert from micromoles/L to moles/L
Sled_WSA <- Sled_WSA[2:6,]
N_field <- approxfun(x = c(1*24, 37*24, 68*24, 86*24, 107*24), y = c(Sled_WSA$Nmol_combined), method = "linear", rule = 2) #N forcing function

###### C forcing set-up ###########
Sled_DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE) #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

###### NOAA Irradiance forcing set-up ####
#hourly <- seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
NOAA_Irradiance_Sledy2_L2 <-  NOAA_Irradiance$PAR[6134:6822]
I_field <- approxfun(x = seq(from = 0, to = 2064, by = 3), y = NOAA_Irradiance_Sledy2_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
Sled_Y2_Hobo <- read.csv("Sled_Y2_HoboLightTemp.csv", header = TRUE)
Sled_Y2_Hobo$DateTime <- mdy_hms(Sled_Y2_Hobo$DateTime)
Sled_Y2_Hobo$Temp_K <- Sled_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
SledY2T_hourly <- ceiling_date(Sled_Y2_Hobo$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(Sled_Y2_Hobo$Temp_K, by=list(SledY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[1346:3385,]
fd <- rep(285, 25)
T_field <- approxfun(x = c(0:2064), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Sled2_Y2 <- T_field(0:2064)
##############
Sled2_Y2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(Sled2_Y2_SA))
##############

#Dredge1_Y2
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE) #Import water Q data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #for RL
Dredge_WSA <- filter(WSA_Y2, Site == "Moonstone Dredge")
daily <- seq(as.Date("2018-12-12"), as.Date("2019-05-03"), by="days")
key <- as.data.frame(daily)
Dredge_WSA$Ncombined <- Dredge_WSA$NO3NO2_µM + Dredge_WSA$NH4_µM
Dredge_WSA$Nmol_combined <- Dredge_WSA$Ncombined/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(1*24, 93*24, 124*24, 142*24, 163*24), y = c(Dredge_WSA$Nmol_combined), method = "linear", rule = 2) #N forcing function

###### C forcing set-up ###########
Sled_DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE) #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

###### NOAA Irradiance forcing set-up ####
#hourly <- seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
NOAA_Irradiance_Dredgey2 <-  NOAA_Irradiance$PAR[5686:6822]
I_field <- approxfun(x = seq(from = 0, to = 3408, by = 3), y = NOAA_Irradiance_Dredgey2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
Dredge_Y2_Hobo <- read.csv("Dredge_Y2_HoboTempLight.csv", header = TRUE)
Dredge_Y2_Hobo$DateTime <- mdy_hms(Dredge_Y2_Hobo$DateTime)
Dredge_Y2_Hobo$Temp_K <- Dredge_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
DredgeY2T_hourly <- ceiling_date(Dredge_Y2_Hobo$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(Dredge_Y2_Hobo$Temp_K, by=list(DredgeY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[2:3384,] #3385
fd <- rep(285, 26)
T_field <- approxfun(x = c(0:3408), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Dredge1_Y2 <- T_field(0:3408)
##############
Dredge1_Y2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(Dredge1_Y2_SA))
##############

#Dredge2_Y2
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE) #Import water Q data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #for RL
Dredge_WSA <- filter(WSA_Y2, Site == "Moonstone Dredge")
daily <- seq(as.Date("2019-02-06"), as.Date("2019-05-03"), by="days")
key <- as.data.frame(daily)
Dredge_WSA$Ncombined <- Dredge_WSA$NO3NO2_µM + Dredge_WSA$NH4_µM
Dredge_WSA$Nmol_combined <- Dredge_WSA$Ncombined/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(-25*24, 37*24, 68*24, 86*24, 107*24), y = c(Dredge_WSA$Nmol_combined), method = "linear", rule = 2) #N forcing function

###### C forcing set-up ###########
Sled_DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE) #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

###### NOAA Irradiance forcing set-up ####
#hourly <- seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
NOAA_Irradiance_Dredgey2_L2 <-  NOAA_Irradiance$PAR[6134:6822]
I_field <- approxfun(x = seq(from = 0, to = 2064, by = 3), y = NOAA_Irradiance_Dredgey2_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
Dredge_Y2_Hobo <- read.csv("Dredge_Y2_HoboTempLight.csv", header = TRUE)
Dredge_Y2_Hobo$DateTime <- mdy_hms(Dredge_Y2_Hobo$DateTime)
Dredge_Y2_Hobo$Temp_K <- Dredge_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
DredgeY2T_hourly <- ceiling_date(Dredge_Y2_Hobo$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(Dredge_Y2_Hobo$Temp_K, by=list(DredgeY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[1346:3385,]
fd <- rep(285, 25)
T_field <- approxfun(x = c(0:2064), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Dredge2_Y2 <- T_field(0:2064)
##############
Dredge2_Y2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(Dredge2_Y2_SA))
##############

#Wickford1_Y2
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE) #Import water Q data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #for RL
Wickford_WSA <- filter(WSA_Y2, Site == "Wickford")

daily <- seq(as.Date("2018-12-19"), as.Date("2019-05-23"), by="days")
key <- as.data.frame(daily)
Wickford_WSA$Ncombined <- Wickford_WSA$NO3NO2_µM + Wickford_WSA$NH4_µM
Wickford_WSA$Nmol_combined <- Wickford_WSA$Ncombined/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(1*24, 55*24, 85*24, 156*24), y = c(Wickford_WSA$Nmol_combined), method = "linear", rule = 2) #N forcing function

###### C forcing set-up ###########
#BrentonPoint_Segarra2002CarbonData
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE) #Import lit TCO2 data
names(Segarra2002Carbon)[1] <- "Date" # RL: had to add this line (probably due to some encoding glitch in loaded file)
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date)
mean(Segarra2002Carbon$TCO2_micromolPERkg)
#ggplot() + 
#geom_point(data = Segarra2002Carbon, aes(Date, TCO2_micromolPERkg))
CO_2 <- 1956.143/1000000 #(mol CO2/L)

###### NOAA Irradiance forcing set-up ####
#hourly <- seq(as_datetime("2018-12-19 12:00:00"), as_datetime("2019-05-23 12:00:00"), by="hour")
NOAA_Irradiance_Wickfordy2 <-  NOAA_Irradiance$PAR[5742:6982]
I_field <- approxfun(x = seq(from = 0, to = 3720, by = 3), y = NOAA_Irradiance_Wickfordy2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
hourly <- seq(as_datetime("2018-12-19 12:00:00"), as_datetime("2019-05-23 12:00:00"), by="hour")
key <- as.data.frame(hourly)

Wickford_Y2_Hobo <- read.csv("Wickford_Y2_HoboLightTemp.csv", header = TRUE)
Wickford_Y2_Hobo$DateTime <- mdy_hms(Wickford_Y2_Hobo$DateTime)
Wickford_Y2_Hobo$Temp_K <- Wickford_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
WickfordY2T_hourly <- ceiling_date(Wickford_Y2_Hobo$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(Wickford_Y2_Hobo$Temp_K, by=list(WickfordY2T_hourly), mean) #calculate average hourly temp
fd <- rep(278, 4)
AvgTempKbyhr_sub <- AvgTempKbyhr[4:3716,]
fd2 <- rep(287, 4)
T_field <- approxfun(x = c(0:3720), y = c(fd, AvgTempKbyhr_sub$x, fd2), method = "linear", rule = 2) #the temp forcing function
T_Wickford1_Y2 <- T_field(0:3720)
##############
Wickford1_Y2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(Wickford1_Y2_SA))
##############

#RomePt1_Y2
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE) #Import water Q data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #for RL
RomePt_WSA <- filter(WSA_Y2, Site == "Rome Point")
daily <- seq(as.Date("2018-12-20"), as.Date("2019-05-24"), by="days")
key <- as.data.frame(daily)
RomePt_WSA$Ncombined <- RomePt_WSA$NO3NO2_µM + RomePt_WSA$NH4_µM
RomePt_WSA$Nmol_combined <- RomePt_WSA$Ncombined/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(1*24, 64*24, 84*24, 155*24), y = c(RomePt_WSA$Nmol_combined), method = "linear", rule = 2) #N forcing function

###### C forcing set-up ###########
#BrentonPoint_Segarra2002CarbonData
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE) #Import lit TCO2 data
names(Segarra2002Carbon)[1] <- "Date" # RL: had to add this line (probably due to some encoding glitch in loaded file)
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date)
mean(Segarra2002Carbon$TCO2_micromolPERkg)
#ggplot() + 
#geom_point(data = Segarra2002Carbon, aes(Date, TCO2_micromolPERkg))
CO_2 <- 1956.143/1000000 #(mol CO2/L)

###### NOAA Irradiance forcing set-up ####
#hourly <- seq(as_datetime("2018-12-20 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
NOAA_Irradiance_RomePty2 <-  NOAA_Irradiance$PAR[5750:6990]
I_field <- approxfun(x = seq(from = 0, to = 3720, by = 3), y = NOAA_Irradiance_RomePty2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
hourly <- seq(as_datetime("2018-12-20 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
key <- as.data.frame(hourly)

RomePt_Y2_Hobo <- read.csv("RomePt_Y2_HoboTempLight.csv", header = TRUE)
RomePt_Y2_Hobo$DateTime <- mdy_hm(RomePt_Y2_Hobo$DateTime)
RomePt_Y2_Hobo$Temp_K <- RomePt_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
RomePtY2T_hourly <- ceiling_date(RomePt_Y2_Hobo$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(RomePt_Y2_Hobo$Temp_K, by=list(RomePtY2T_hourly), mean) #calculate average hourly temp

fd <- rep(280, 4) #4
AvgTempKbyhr_sub <- AvgTempKbyhr[28:2414,] #2387
#Using Wickford temp to fill in th gap in the Rome Pt temp
Wickford_Y2_Hobo <- read.csv("Wickford_Y2_HoboLightTemp.csv", header = TRUE)
Wickford_Y2_Hobo$DateTime <- mdy_hms(Wickford_Y2_Hobo$DateTime)
Wickford_Y2_Hobo$Temp_K <- Wickford_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
WickfordY2T_hourly <- ceiling_date(Wickford_Y2_Hobo$DateTime, unit = "hour")
AvgTempKbyhr_Wickford <- aggregate(Wickford_Y2_Hobo$Temp_K, by=list(WickfordY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_W <- AvgTempKbyhr_Wickford[2415:3716,] #1302
fd2 <- rep(287, 28)
T_field <- approxfun(x = c(0:3720), y = c(fd, AvgTempKbyhr_sub$x, AvgTempKbyhr_W$x, fd2), method = "linear", rule = 2) #the temp forcing function
T_RomePt1_Y2 <- T_field(0:3720)
##############
RomePt1_Y2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(RomePt1_Y2_SA))
##############

#RomePt2_Y2
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE) #Import water Q data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #for RL
RomePt_WSA <- filter(WSA_Y2, Site == "Rome Point")
daily <- seq(as.Date("2019-2-21"), as.Date("2019-05-24"), by="days")
key <- as.data.frame(daily)
RomePt_WSA$Ncombined <- RomePt_WSA$NO3NO2_µM + RomePt_WSA$NH4_µM
RomePt_WSA$Nmol_combined <- RomePt_WSA$Ncombined/1000000 #convert from micromoles/L to moles/L
RomePt_WSA <- RomePt_WSA[2:4,]
N_field <- approxfun(x = c(1*24, 21*24, 92*24), y = c(RomePt_WSA$Nmol_combined), method = "linear", rule = 2) #N forcing function

###### C forcing set-up ###########
#BrentonPoint_Segarra2002CarbonData
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE) #Import lit TCO2 data
names(Segarra2002Carbon)[1] <- "Date" # RL: had to add this line (probably due to some encoding glitch in loaded file)
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date)
mean(Segarra2002Carbon$TCO2_micromolPERkg)
#ggplot() + 
#geom_point(data = Segarra2002Carbon, aes(Date, TCO2_micromolPERkg))
CO_2 <- 1956.143/1000000 #(mol CO2/L)

###### NOAA Irradiance forcing set-up ####
#hourly <- seq(as_datetime("2019-2-21 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
NOAA_Irradiance_RomePty2_L2 <-  NOAA_Irradiance$PAR[6254:6990]
I_field <- approxfun(x = seq(from = 0, to = 2208, by = 3), y = NOAA_Irradiance_RomePty2_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
hourly <- seq(as_datetime("2019-2-21 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
key <- as.data.frame(hourly)

RomePt_Y2_Hobo <- read.csv("RomePt_Y2_HoboTempLight.csv", header = TRUE)
RomePt_Y2_Hobo$DateTime <- mdy_hm(RomePt_Y2_Hobo$DateTime)
RomePt_Y2_Hobo$Temp_K <- RomePt_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
RomePtY2T_hourly <- ceiling_date(RomePt_Y2_Hobo$DateTime, unit = "hour")
AvgTempKbyhr <- aggregate(RomePt_Y2_Hobo$Temp_K, by=list(RomePtY2T_hourly), mean) #calculate average hourly temp

AvgTempKbyhr_sub <- AvgTempKbyhr[1536:2414,] #879
#Using Wickford temp to fill in th gap in the Rome Pt temp
Wickford_Y2_Hobo <- read.csv("Wickford_Y2_HoboLightTemp.csv", header = TRUE)
Wickford_Y2_Hobo$DateTime <- mdy_hms(Wickford_Y2_Hobo$DateTime)
Wickford_Y2_Hobo$Temp_K <- Wickford_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
WickfordY2T_hourly <- ceiling_date(Wickford_Y2_Hobo$DateTime, unit = "hour")
AvgTempKbyhr_Wickford <- aggregate(Wickford_Y2_Hobo$Temp_K, by=list(WickfordY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_W <- AvgTempKbyhr_Wickford[2415:3716,] #1302
fd2 <- rep(287, 28)
T_field <- approxfun(x = c(0:2208), y = c(AvgTempKbyhr_sub$x, AvgTempKbyhr_W$x, fd2), method = "linear", rule = 2) #the temp forcing function
T_RomePt2_Y2 <- T_field(0:2208)
##############
RomePt2_Y2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
plot(summary(RomePt2_Y2_SA))
##############

Compile_SA <- rbind(Sled1_SA, Sled2_SA, Dredge1_SA, Dredge2_SA, Wickford1_SA, RomePt1_SA, RomePt2_SA, Sled1_Y2_SA, Sled2_Y2_SA, Dredge1_Y2_SA, Dredge2_Y2_SA, Wickford1_Y2_SA, RomePt1_Y2_SA, RomePt2_Y2_SA)
summary(Compile_SA)
plot(summary(Compile_SA))
SA_SUM <- summary(Compile_SA)

SA_SUM[ "Parameter" ] <- rownames(SA_SUM)
SA_SUM <- subset(SA_SUM, Parameter != "n_NV")
SA_SUM$Parameter[SA_SUM$Parameter=="JENAM"] <- "J_EN_Am"
SA_SUM$Parameter[SA_SUM$Parameter=="JCO2M"] <- "J_CO2_m"
SA_SUM$Parameter[SA_SUM$Parameter=="JECAM"] <- "J_EC_Am"
SA_SUM$Parameter[SA_SUM$Parameter=="rho_PSU"] <- "ρ_PSU"
SA_SUM$Parameter[SA_SUM$Parameter=="y_CO2_C"] <- "y_CO2_C"
SA_SUM$Parameter[SA_SUM$Parameter=="JENM"] <- "J_EN_M"
SA_SUM$Parameter[SA_SUM$Parameter=="JECM"] <- "J_EC_M"
#SA_SUM$Parameter[SA_SUM$Parameter=="JENAM"] <-  "Max. volume spec. N assimilation"
#SA_SUM$Parameter[SA_SUM$Parameter=="K_N"] <- "Half-saturation concentration for NO3- and NH4+ uptake"
#SA_SUM$Parameter[SA_SUM$Parameter=="JCO2M"] <- "Max. volume spec. CO2 uptake rate"
#SA_SUM$Parameter[SA_SUM$Parameter=="K_C"] <- "Half-saturation concentration for CO2 uptake"
#SA_SUM$Parameter[SA_SUM$Parameter=="JECAM"] <- "Max. volume spec. C assimilation"
#SA_SUM$Parameter[SA_SUM$Parameter=="rho_PSU"] <- "Photosynthetic unit (PSU) density"
#SA_SUM$Parameter[SA_SUM$Parameter=="b_I"] <- "Binding probability of a photon to a free light SU"
#SA_SUM$Parameter[SA_SUM$Parameter=="k_I"] <- "Dissociation rate of releasing ATP and NADPH+"
#SA_SUM$Parameter[SA_SUM$Parameter=="y_I_C"] <- "Yield factor of C reserve to NADPH"
#SA_SUM$Parameter[SA_SUM$Parameter=="y_CO2_C"] <- "Yield factor of C reserve to CO2"
#SA_SUM$Parameter[SA_SUM$Parameter=="kE_C"] <- "C reserve turnover rate"
#SA_SUM$Parameter[SA_SUM$Parameter=="kE_N"] <- "N reserve turnover rate"
#SA_SUM$Parameter[SA_SUM$Parameter=="kappa_Ei"] <- "Fraction of rejection flux incorporated back in i-reserve"
#SA_SUM$Parameter[SA_SUM$Parameter=="y_EN_V"] <- "Yield factor of N reserve to structure"
#SA_SUM$Parameter[SA_SUM$Parameter=="y_EC_V"] <- "Yield factor of C reserve to structure"
#SA_SUM$Parameter[SA_SUM$Parameter=="JENM"] <- "Volume spec. maintenance cost paid by N reserve"
#SA_SUM$Parameter[SA_SUM$Parameter=="JECM"] <- "Volumes pecific maintenance cost paid by C reserve"
#SA_SUM$Parameter[SA_SUM$Parameter=="T_A"] <- "Arrhenius temperature"
#SA_SUM$Parameter[SA_SUM$Parameter=="T_H"] <- "Upper boundary of temperature tolerance"
#SA_SUM$Parameter[SA_SUM$Parameter=="T_L"] <- "Lower boundary of temperature tolerance"
#SA_SUM$Parameter[SA_SUM$Parameter=="T_AH"] <- "Arrhenius temperature outside T_H"
#SA_SUM$Parameter[SA_SUM$Parameter=="T_AL"] <- "Arrhenius temperature outside T_L"
#SA_SUM$Parameter[SA_SUM$Parameter=="T_0"] <- "Reference temperature"

SA_SUM$Parameter <- factor(SA_SUM$Parameter, as.character(SA_SUM$Parameter)) #to keep ggplot from automatically alphabetizing

ggplot(data = SA_SUM) +
  geom_point(mapping = aes(x = L1, y = Parameter), size = 3) +
  geom_abline(mapping = aes(slope = 0, intercept = 0)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "L1 (sensitivity function)", y = "Parameter Name")


ident <- collin(test) #Inf collinearity when trying to include any large number of the parameters

ident2 <- collin(test, parset = c("JCO2M", "K_C", "kE_C", "kE_N", "kappa_Ei", "JENM", "JECM"))


