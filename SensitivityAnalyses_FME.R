#sensitivity analyses using sensFun (Local sensitivity analyses)
library(FME)
library(tidyverse)
library(lubridate)

source("SolveR_R.R")
##### Minerals and Organics Section #####
#Conversion coefficients, organics (n = matrix of chemical indices)
# "food N" "food C" Stucture "N reserves" "C reserves" products
#     X_N   X_C      V    E_N    E_C    P
n_O <- matrix(
  + c(0.00, 1.00, 1.00, 0.00, 1.00, 1.00,  #C/C, equals 1 by definition
      + 0.00, 0.50, 1.33, 0.00, 2.00, 1.80,  #H/C, these values show that we consider dry-mass
      + 3.00, 2.50, 1.00, 2.50, 1.00, 0.50,  #O/C
      + 1.00, 0.00, 0.04, 1.00, 0.00, 0.04), nrow=4, ncol=6, byrow = TRUE) #N/C
#V is the C-mol structure of alginate (Alginic acid: (C6H8O6)n)
#E_N is N03- and N02- averaged
#E_C is glucose C6H12O6 (Laminarin: c18h32o16 and mannitol c6h14o6)
#We aren't using the X_N, X_C, or P collumn here

#Molecular weights
#t() is a matrix transpose function
#organics structure matrix multiplied by the atomic masses (mass in grams of one mole of an element) of C H O N
w_O_step <- t(n_O)*matrix(c(12, 1, 16, 14), nrow=6, ncol=4, byrow= TRUE) #g/mol, molecular weights for organics
w_O <- rowSums(w_O_step) #this provides g/mol of each of the six "pockets of mass" (i.e. X_N, X_C)

#define molecular weights
w_V <- w_O[3]  # g/mol       #molecular weight of structure
w_EN <- w_O[4]  # g/mol      #molecular weight of N reserve
w_EC <- w_O[5]  #g/mol       #molecular weight of C reserve
w_O2 <- 32 #g/mol
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##### Parameters compiled #####
params_Lo <- c(#maximum volume-specific assimilation rate of N before temperature correction
  JENAM = 1.5e-4, #mol N / molM_V / h
  #maximum surface-specific assimilation rate of N
  K_N = 2.5e-6, #molNO3 and NO2/L
  #max volume-specific carbon dioxide assimilation rate
  JCO2M = 0.0075, #molC/molM_V/h
  #half saturation constant of C uptake
  K_C = 4e-7, #mol DIC/L
  #maximum volume-specific carbon assimilation rate
  JECAM = 0.282, #molC/molM_V/h
  #Photosynthetic unit density
  rho_PSU = 0.5, #mol PSU/ mol Mv
  #binding probability of photons to a Light SU
  b_I = 0.5, #dimensionless
  #Specific photon arrival cross section
  alpha = 1, #m^2 mol PSU–1
  #dissociation rate
  k_I = 0.075, #molγ molPS–1 h–1
  #Yield factor of C reserve to photon
  y_I_C = 10, #mol γ mol C-1
  #Yield factor of C reserve to CO2
  y_CO2_C = 1, #mol CO2 mol C-1
  #Yield factor of photon to O2
  y_LO2 = 0.125, #molO2 molγ –1
  #reserve turnover
  kE_C = 0.05,
  kE_N = 0.01, 
  #fraction of rejection flux from growth SU incorporated back into i-reserve
  kappa_Ei = 0.9, #dimensionless
  #yield of structure on N reserve (percent of N in structure)
  y_EN_V = 0.04, #mol N/mol M_V
  #yield of structure on C reserve (percent of C in structure)
  y_EC_V = 1, #mol C/mol M_V
  #specific maintenance costs requiring N before temp correction
  JENM = 4e-6, #mol N/molM_V/h
  #specific maintenance costs requiring C before temp correction
  JECM = 1e-6, #mol C/molM_V/h
  #Arrhenius temperature
  T_A = 6314.3, # K
  #Upper boundary of temperature tolerance
  T_H = 13.386 + 273.15, # K
  #Lower boundary of temperature tolerance
  T_L = 273.15, # K
  #Arrhenius temperature outside T_H
  T_AH = 18702, #K
  #Arrhenius temperature outside T_L
  T_AL = 4391.9, #K
  #temperature at which rate parameters are given
  T_0 = 20 + 273.15) # K
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### Set up NOAA data (for Irradiance forcing) ####
#NOAA irradiance data set-up: NOAASurfaceIrradiance
NOAA_Irradiance <- read.csv("NOAASurfaceIrradiance.csv", header = TRUE, fileEncoding="UTF-8-BOM")
NOAA_Irradiance$DateTime <- dmy_hms(NOAA_Irradiance$DateTime, tz = "UTC") #NOAA data in UTC (5 hours ahead)
NOAA_Irradiance <- with_tz(NOAA_Irradiance, "America/New_York") #Convert from UTC to EST
NOAA_Irradiance$DownMinusUp <- NOAA_Irradiance$dswrf-NOAA_Irradiance$uswrf #net shortwave radiation at the surface (W/m^2) is obtained by subtracting the upward short wave flux (uswrf) from the downward flux (dswrf)
#PAR = NSW*PAR_frac*C*exp(-k*z)*3600
#NSW=dswrf-uswrf
#PAR_frac is the fraction of the incident flux that is useable for photosynthesis
#C is a conversion factor = 4.56 umol photons/s/W
#k is the extinction coefficient
#3600 converts from s^-1 to h^-1
#1e-6 converts from micomoles to moles
NOAA_Irradiance$PAR <- NOAA_Irradiance$DownMinusUp*0.43*4.56*exp(-0.46*1)*3600*1e-6

#######

SenseFunc <- function(params_Lo){
  rates_Lo <- function(t, state, parameters) { #note that the time, state, and parameters in this function definition should always remain in this format regardless of how these pieces of info are labelled
    
    with(as.list(c(state, parameters)), { 
      
      #set-up equations (to calculate values necessary for the differential equations)
      #Temperature correction:
      C_T <- exp((T_A/T_0)-(T_A/T_field(t))) * (1+exp((T_AL/T_0)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_0))) * ((1+exp((T_AL/T_field(t))-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_field(t))))^-1)
      
      J_EN_AM <-  JENAM * C_T #temperature correction of max assimilation rate of N
      #Specific assimilation rate of N
      J_EN_A <- J_EN_AM*(N_field(t)/(N_field(t)+K_N))  #J_EN_A units: mol N/molM_V/h
      N <- N_field(t) #just to save N forcing output
      
      J_CO2_M <- JCO2M * C_T #Temperature correction max uptake rate of CO2 (photon capture not influenced by temp) J_CO2_M units: molC/molM_V/d
      #Specific CO2 uptake flux
      J_CO2 <- J_CO2_M*(CO_2/(CO_2+K_C)) #J_CO2 units: molC/molM_V/h
      
      #Specific relaxation rate
      J_I <- (rho_PSU * I_field(t) * b_I * alpha)/(1+I_field(t) * b_I * alpha/(k_I * C_T)) #Unit: molγ molM_V–1 h–1
      J_EC_AM <- JECAM * C_T #temperature correction max assimilation rate of C
      #Specific carbon assimilation flux
      J_EC_A <- ((1/J_EC_AM)+(1/(J_CO2/y_CO2_C))+(1/(J_I/y_I_C))-(1/((J_CO2/y_CO2_C)+(J_I/y_I_C))))^-1 #Units: mol C/mol M_V/h
      I <- I_field(t) #just to save I forcing output
      
      #rate of oxygen production
      J_O <- J_I * (M_V/W) * w_O2 * y_LO2 # g O2/g/h O2 production
      
      #Temperature corrections for the more internal parts of the model
      k_E <- c(kE_C * C_T, kE_N * C_T) #reserve turnover rate #Unit: 1/h
      J_EN_M <- JENM * C_T #volume specific maintenance cost paid by the N reserve #unit: molEN/molM_V/h
      J_EC_M <-JECM * C_T #volume specific maintenance cost paid by the C reserve #unit: molEC/molM_V/h
      
      #########################
      #setting up holder vectors for the SolveR_R function
      #CARBON in first slot, NITROGEN in second
      m_E <- c(m_EC, m_EN)      #density of both reserves (mol Ei/molM_V)
      J_EM <- c(J_EC_M, J_EN_M) #volume specific maintenance cost paid by both reserves (molEi/molM_V/h)
      y_EV <- c(y_EC_V, y_EN_V) #yield of structure on both reserves (mol Ei/mol M_V)
      J_VM_C <- J_EC_M/y_EC_V  #rate of carbon flux of structure that pays for maintenance when the catabolic flux is not enough (1/h)
      J_VM_N <- J_EN_M/y_EN_V #rate of nitrogen flux of structure that pays for maintenance when the catabolic flux is not enough (1/h)
      J_VM <- c(J_VM_C, J_VM_N) #rate of maintenance costs paid from structure when the catabolic flux is not enough (allows r to be negative) (1/h)
      r0 <- 0.01 # RL: I would not start the regression procedure with 0 (not sure why but it could have consequences in the procedure)
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
      info <- Output_loop[10] #1 if convergence, 0 if no convergence
      
      #Allocation to growth
      #Specific growth flux for nitrogen 
      J_EN_G <- J_EN_C-J_EN_M #unit: mol EN/molM_V/h
      #Specific growth flux for carbon
      J_EC_G <- J_EC_C-J_EC_M #unit: mol EC/molM_V/h
      #this function keeps J_G from being NaN when growth is not occuring
      #Necessary because Inf-Inf is NaN, but Inf-0 is Inf
      RealNumbas <- function(x) {
        if(is.infinite(x) == TRUE){
          x <- 0 #replace infinite numbers with zero
        }
        return(x)
      }
      #Gross specific growth flux
      J_G <- ((((J_EN_G/y_EN_V)^-1)+ ((J_EC_G/y_EC_V)^-1))-RealNumbas((J_EN_G/y_EN_V+J_EC_G/y_EC_V)^-1))^-1 #unit: 1/h
      
      #Biomass (structure mass+ C reserve mass + N reserve mass)
      W <- (w_V+m_EN*w_EN+m_EC*w_EC)*M_V #unit: g (as long as the units of M_V are mol M_V)
      #Allometic relationship between length (cm) and dry weight (g) from Gevaert (2001)
      L_allometric <- (W/0.00387)^(1/1.469) #unit: cm
      
      #rates
      #dynamics of the C reserve
      dm_ECdt <- J_EC_A-J_EC_C+(kappa_Ei*J_EC_R)-(r*m_EC)  #unit: mol C/mol M_V/h
      #dynamics of the N reserve
      dm_ENdt <- J_EN_A-J_EN_C+(kappa_Ei*J_EN_R)-(r*m_EN) #unit: mol N/mol M_V/h
      #dynamics of structural mass
      dM_Vdt <- r*M_V  #units of M_V are mol M_V, so dMVdt units are mol MV/h
      
      
      #output
      return(list(c(dm_ECdt, dm_ENdt, dM_Vdt), r=r, W=W, C_T=C_T, I=I, N=N, J_CO2=J_CO2, J_I=J_I, J_EC_A=J_EC_A, J_EN_A=J_EN_A, J_EC_C=J_EC_C, J_EN_C=J_EN_C, J_EN_M=J_EN_M, J_EC_M=J_EC_M, J_EN_G=J_EN_G, J_EC_G=J_EC_G, J_G=J_G, J_VM=J_VM, J_EC_R=J_EC_R, J_EN_R=J_EN_R, CO_2=CO_2, L_allometric=L_allometric, J_O=J_O, info=info))
    })
  }
  ####### State Inititial conditions ############
  #Initial conditions of state variables
  #these values are not coming from any field data or literature information, estimated
  state_Lo <- c(m_EC = 0.002, #0.1, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
                m_EN = 0.01, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
                M_V = 0.05/(w_V+0.01*w_EN+0.002*w_EC)) #molM_V #initial mass of structure
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #######Time step to run the model on#######
  #(First number of time step, last number of time step, interval to step)
  times_Lo_Sled1 <- seq(0, 4008, 1) #167 days stepped hourly
  times_Lo_Sled2 <- seq(0, 3336, 1) #139 days stepped hourly
  times_Lo_Dredge1 <- seq(0, 4128, 1) #172 days stepped hourly
  times_Lo_Dredge2 <- seq(0, 3456, 1) #144 days stepped hourly
  times_Lo_Wickford1 <- seq(0, 3312, 1) #138 days stepped hourly
  times_Lo_RomePt1 <- seq(0, 4104, 1) #171 days stepped hourly
  times_Lo_RomePt2 <- seq(0, 3264, 1) #136 days stepped hourly
  
  times_Y2_Sled1 <- seq(0, 3408, 1) #142 days stepped hourly
  times_Y2_Sled2 <- seq(0, 2064, 1) #86 days stepped hourly
  times_Y2_Dredge1 <- seq(0, 3408, 1) #142 days stepped hourly
  times_Y2_Dredge2 <- seq(0, 2064, 1) #86 days stepped hourly
  times_Y2_Wickford1 <- seq(0, 3720, 1) #155 days stepped hourly
  times_Y2_RomePt1 <- seq(0, 3720, 1) #155 days stepped hourly
  times_Y2_RomePt2 <- seq(0, 2208, 1) #92 days stepped hourly
  
  #########
  return(as.data.frame(ode(y = state_Lo, times = times_Lo_Sled1, func = rates_Lo, parms = params_Lo)))
}

SenseFunc2 <- function(params_Lo){
  rates_Lo <- function(t, state, parameters) { #note that the time, state, and parameters in this function definition should always remain in this format regardless of how these pieces of info are labelled
    
    with(as.list(c(state, parameters)), { 
      
      #set-up equations (to calculate values necessary for the differential equations)
      #Temperature correction:
      C_T <- exp((T_A/T_0)-(T_A/T_field(t))) * (1+exp((T_AL/T_0)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_0))) * ((1+exp((T_AL/T_field(t))-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_field(t))))^-1)
      
      J_EN_AM <-  JENAM * C_T #temperature correction of max assimilation rate of N
      #Specific assimilation rate of N
      J_EN_A <- J_EN_AM*(N_field(t)/(N_field(t)+K_N))  #J_EN_A units: mol N/molM_V/h
      N <- N_field(t) #just to save N forcing output
      
      J_CO2_M <- JCO2M * C_T #Temperature correction max uptake rate of CO2 (photon capture not influenced by temp) J_CO2_M units: molC/molM_V/d
      #Specific CO2 uptake flux
      J_CO2 <- J_CO2_M*(CO_2/(CO_2+K_C)) #J_CO2 units: molC/molM_V/h
      
      #Specific relaxation rate
      J_I <- (rho_PSU * I_field(t) * b_I * alpha)/(1+I_field(t) * b_I * alpha/(k_I * C_T)) #Unit: molγ molM_V–1 h–1
      J_EC_AM <- JECAM * C_T #temperature correction max assimilation rate of C
      #Specific carbon assimilation flux
      J_EC_A <- ((1/J_EC_AM)+(1/(J_CO2/y_CO2_C))+(1/(J_I/y_I_C))-(1/((J_CO2/y_CO2_C)+(J_I/y_I_C))))^-1 #Units: mol C/mol M_V/h
      I <- I_field(t) #just to save I forcing output
      
      #rate of oxygen production
      J_O <- J_I * (M_V/W) * w_O2 * y_LO2 # g O2/g/h O2 production
      
      #Temperature corrections for the more internal parts of the model
      k_E <- c(kE_C * C_T, kE_N * C_T) #reserve turnover rate #Unit: 1/h
      J_EN_M <- JENM * C_T #volume specific maintenance cost paid by the N reserve #unit: molEN/molM_V/h
      J_EC_M <-JECM * C_T #volume specific maintenance cost paid by the C reserve #unit: molEC/molM_V/h
      
      #########################
      #setting up holder vectors for the SolveR_R function
      #CARBON in first slot, NITROGEN in second
      m_E <- c(m_EC, m_EN)      #density of both reserves (mol Ei/molM_V)
      J_EM <- c(J_EC_M, J_EN_M) #volume specific maintenance cost paid by both reserves (molEi/molM_V/h)
      y_EV <- c(y_EC_V, y_EN_V) #yield of structure on both reserves (mol Ei/mol M_V)
      J_VM_C <- J_EC_M/y_EC_V  #rate of carbon flux of structure that pays for maintenance when the catabolic flux is not enough (1/h)
      J_VM_N <- J_EN_M/y_EN_V #rate of nitrogen flux of structure that pays for maintenance when the catabolic flux is not enough (1/h)
      J_VM <- c(J_VM_C, J_VM_N) #rate of maintenance costs paid from structure when the catabolic flux is not enough (allows r to be negative) (1/h)
      r0 <- 0.01 # RL: I would not start the regression procedure with 0 (not sure why but it could have consequences in the procedure)
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
      info <- Output_loop[10] #1 if convergence, 0 if no convergence
      
      #Allocation to growth
      #Specific growth flux for nitrogen 
      J_EN_G <- J_EN_C-J_EN_M #unit: mol EN/molM_V/h
      #Specific growth flux for carbon
      J_EC_G <- J_EC_C-J_EC_M #unit: mol EC/molM_V/h
      #this function keeps J_G from being NaN when growth is not occuring
      #Necessary because Inf-Inf is NaN, but Inf-0 is Inf
      RealNumbas <- function(x) {
        if(is.infinite(x) == TRUE){
          x <- 0 #replace infinite numbers with zero
        }
        return(x)
      }
      #Gross specific growth flux
      J_G <- ((((J_EN_G/y_EN_V)^-1)+ ((J_EC_G/y_EC_V)^-1))-RealNumbas((J_EN_G/y_EN_V+J_EC_G/y_EC_V)^-1))^-1 #unit: 1/h
      
      #Biomass (structure mass+ C reserve mass + N reserve mass)
      W <- (w_V+m_EN*w_EN+m_EC*w_EC)*M_V #unit: g (as long as the units of M_V are mol M_V)
      #Allometic relationship between length (cm) and dry weight (g) from Gevaert (2001)
      L_allometric <- (W/0.00387)^(1/1.469) #unit: cm
      
      #rates
      #dynamics of the C reserve
      dm_ECdt <- J_EC_A-J_EC_C+(kappa_Ei*J_EC_R)-(r*m_EC)  #unit: mol C/mol M_V/h
      #dynamics of the N reserve
      dm_ENdt <- J_EN_A-J_EN_C+(kappa_Ei*J_EN_R)-(r*m_EN) #unit: mol N/mol M_V/h
      #dynamics of structural mass
      dM_Vdt <- r*M_V  #units of M_V are mol M_V, so dMVdt units are mol MV/h
      
      
      #output
      return(list(c(dm_ECdt, dm_ENdt, dM_Vdt), r=r, W=W, C_T=C_T, I=I, N=N, J_CO2=J_CO2, J_I=J_I, J_EC_A=J_EC_A, J_EN_A=J_EN_A, J_EC_C=J_EC_C, J_EN_C=J_EN_C, J_EN_M=J_EN_M, J_EC_M=J_EC_M, J_EN_G=J_EN_G, J_EC_G=J_EC_G, J_G=J_G, J_VM=J_VM, J_EC_R=J_EC_R, J_EN_R=J_EN_R, CO_2=CO_2, L_allometric=L_allometric, J_O=J_O, info=info))
    })
  }
  ####### State Inititial conditions ############
  #Initial conditions of state variables
  state_LoY2 <- c(m_EC = 0.01, #0.9 #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
                  m_EN = 0.09, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
                  M_V = 0.05/(w_V+0.09*w_EN+0.01*w_EC)) #molM_V #initial mass of structure
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #######Time step to run the model on#######
  #(First number of time step, last number of time step, interval to step)
  times_Lo_Sled1 <- seq(0, 4008, 1) #167 days stepped hourly
  times_Lo_Sled2 <- seq(0, 3336, 1) #139 days stepped hourly
  times_Lo_Dredge1 <- seq(0, 4128, 1) #172 days stepped hourly
  times_Lo_Dredge2 <- seq(0, 3456, 1) #144 days stepped hourly
  times_Lo_Wickford1 <- seq(0, 3312, 1) #138 days stepped hourly
  times_Lo_RomePt1 <- seq(0, 4104, 1) #171 days stepped hourly
  times_Lo_RomePt2 <- seq(0, 3264, 1) #136 days stepped hourly
  
  times_Y2_Sled1 <- seq(0, 3408, 1) #142 days stepped hourly
  times_Y2_Sled2 <- seq(0, 2064, 1) #86 days stepped hourly
  times_Y2_Dredge1 <- seq(0, 3408, 1) #142 days stepped hourly
  times_Y2_Dredge2 <- seq(0, 2064, 1) #86 days stepped hourly
  times_Y2_Wickford1 <- seq(0, 3720, 1) #155 days stepped hourly
  times_Y2_RomePt1 <- seq(0, 3720, 1) #155 days stepped hourly
  times_Y2_RomePt2 <- seq(0, 2208, 1) #92 days stepped hourly
  
  #########
  return(as.data.frame(ode(y = state_LoY2, times = times_Lo_Sled1, func = rates_Lo, parms = params_Lo)))
}

#Sled1
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water Q data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" #only necessary for some computers running this code
Sled_WSA <- filter(WSA2_Y1, Site == "Sled") #filter for Sled site
daily <- seq(as.Date("2017-11-1"), as.Date("2018-04-17"), by="days") #days kelp in water for this site
N <- Sled_WSA[c("Date","NitrateNitrite_uM")] #subset
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000
#Converted to hourly by multiply by 24
N_field <- approxfun(x = c(161*24, 139*24, 105*24, 0, 28*24, 84*24, 172*24), y = N$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
Sled_DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_Sledy1 <-  NOAA_Irradiance$PAR[2438:3774] # subset based on as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00")
I_field <- approxfun(x = seq(from = 0, to = 4008, by = 3), y = NOAA_Irradiance_Sledy1, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
#Import Sled Hobo (cynlinder) of just temp
Sled_Y1_hobotemp <- read.csv("Sled_Y1_TempLogger2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Sled_Y1_hobotemp$DateTime <- mdy_hms(Sled_Y1_hobotemp$Date_Time) #convert time field
Sled_Y1_hobotemp <- Sled_Y1_hobotemp[14:16049,] #subset
Sled_Y1_hobotemp$Temp_K <- Sled_Y1_hobotemp$Temp_C+273.15 #create collumn with temp in K
SledT_hourly <- ceiling_date(Sled_Y1_hobotemp$DateTime, unit = "hour") #set the values to aggregate around
AvgTempKbyhr <- aggregate(Sled_Y1_hobotemp$Temp_K, by=list(SledT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_part1 <- AvgTempKbyhr$x[0:334] #subset
Dredge_Y1_hobo <- read.csv("Dredge_Y1_hobo.csv", header = TRUE, fileEncoding="UTF-8-BOM") #importing neighboring temp file to replace corrupted section
Dredge_Y1_hobo$DateTime <- mdy_hms(Dredge_Y1_hobo$Date_Time) #convert time field
Dredge_Y1_hobo <- Dredge_Y1_hobo[3:16531,] #subset
Dredge_Y1_hobo$Temp_K <- Dredge_Y1_hobo$Temp_C+273.15 #create collumn with temp in K
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour") #set the values to aggregate around
AvgTempKbyhr4FD <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr4FD <- AvgTempKbyhr4FD[4:4132, ] #subset
fd <- AvgTempKbyhr4FD$x[335:859] #526 data points needed from dredge to replace a weird glitch in the sled temp data
AvgTempKbyhr_part2 <- AvgTempKbyhr$x[860:4009] #later part of original temp file
T_field <- approxfun(x = c(0:4008), y = c(AvgTempKbyhr_part1, fd, AvgTempKbyhr_part2), method = "linear", rule = 2) #the temp forcing function
T_Sled1_Y1 <- T_field(0:4008) #saving the forcing this way for ease of later visualization
#########
Sled1_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(Sled1_SA))
############

#Sled2
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water Q data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" #only necessary on some computers for the code to run
Sled_WSA <- filter(WSA2_Y1, Site == "Sled") #filter for site
daily <- seq(as.Date("2017-11-29"), as.Date("2018-04-17"), by="days") #aquaculture season for this site
N <- Sled_WSA[c("Date","NitrateNitrite_uM")] #create a dataframe with just the relevant collums
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/L
#multiplying by 24 to set as hourly
N_field <- approxfun(x = c(133*24, 111*24, 77*24, -28*24, 0, 56*24, 144*24), y = N$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
Sled_DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_Sledy1_L2 <-  NOAA_Irradiance$PAR[2662:3774] #subset by seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3336, by = 3), y = NOAA_Irradiance_Sledy1_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #############
Sled_Y1_hobotemp <- read.csv("Sled_Y1_TempLogger2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Sled_Y1_hobotemp$DateTime <- mdy_hms(Sled_Y1_hobotemp$Date_Time) #convert time field
Sled_Y1_hobotemp <- Sled_Y1_hobotemp[6:16051,] #subset
Sled_Y1_hobotemp$Temp_K <- Sled_Y1_hobotemp$Temp_C+273.15 #create collumn with temp in K
SledT_hourly <- ceiling_date(Sled_Y1_hobotemp$DateTime, unit = "hour") #set values to aggregate around
AvgTempKbyhr <- aggregate(Sled_Y1_hobotemp$Temp_K, by=list(SledT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[677:4011,] #subset
Dredge_Y1_hobo <- read.csv("Dredge_Y1_hobo.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import nearby site temperature data to fix an error in the Sled data
Dredge_Y1_hobo$DateTime <- mdy_hms(Dredge_Y1_hobo$Date_Time) #convert time field
Dredge_Y1_hobo <- Dredge_Y1_hobo[3:16531,] #subset
Dredge_Y1_hobo$Temp_K <- Dredge_Y1_hobo$Temp_C+273.15 #create collumn with temp in K
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour") #set values to aggregate around
AvgTempKbyhr4FD <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr4FD <- AvgTempKbyhr4FD[4:4132, ] #subset
fd <- AvgTempKbyhr4FD$x[858:859] #526 data points needed from dredge to replace a weird glitch in the sled temp data
T_field <- approxfun(x = c(0:3336), y = c(fd, AvgTempKbyhr$x), method = "linear", rule = 2) #the temp forcing function
T_Sled2_Y1 <- T_field(0:3336) #for later ease in plotting the forcing
#############
Sled2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(Sled2_SA))
#################

#Dredge1
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" #only necessary on some computers to make the code run
Dredge_WSA <- filter(WSA2_Y1, Site == "Dredge") #filter by site
daily <- seq(as.Date("2017-11-1"), as.Date("2018-04-22"), by="days") #date range relevant for this site
N <- Dredge_WSA[c("Date","NitrateNitrite_uM")] #new dataframe with the relevant collumns
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/L
#multipling by 24 to take from daily to hourly
N_field <- approxfun(x = c(139*24, 161*24, 84*24, 0, 105*24), y = N$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000
###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_Dredgey1 <-  NOAA_Irradiance$PAR[2438:3814] #subset by 11/1/17 to 2018-04-22 12:00:00
I_field <- approxfun(x = seq(from = 0, to = 4128, by = 3), y = NOAA_Irradiance_Dredgey1, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
Dredge_Y1_hobo <- read.csv("Dredge_Y1_hobo.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Dredge_Y1_hobo$DateTime <- mdy_hms(Dredge_Y1_hobo$Date_Time) #convert time field
Dredge_Y1_hobo <- Dredge_Y1_hobo[3:16531,] #subset
Dredge_Y1_hobo$Temp_K <- Dredge_Y1_hobo$Temp_C+273.15 #create collumn with temp in K
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour") #set values to aggregate around
AvgTempKbyhr <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[4:4132, ] #subset
T_field <- approxfun(x = c(0:4128), y = AvgTempKbyhr$x, method = "linear", rule = 2) #the temp forcing function
T_Dredge1_Y1 <- T_field(0:4128) #for ease in later plotting of the forcing
##############
Dredge1_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(Dredge1_SA))
##############

#Dredge2
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" #only necessary for some computers to run the following code
Dredge_WSA <- filter(WSA2_Y1, Site == "Dredge") #Filter by site
daily <- seq(as.Date("2017-11-29"), as.Date("2018-04-22"), by="days") #day range relevant to this site
N <- Dredge_WSA[c("Date","NitrateNitrite_uM")] #create new dataframe with just the relevant collumns
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/L
#multiplied by 24 to convert from daily to hourly
N_field <- approxfun(x = c(111*24, 133*24, 56*24, -28*24, 77*24), y = N$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000
###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_Dredgey1_L2 <-  NOAA_Irradiance$PAR[2662:3814] #subset by seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3456, by = 3), y = NOAA_Irradiance_Dredgey1_L2, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
Dredge_Y1_hobo <- read.csv("Dredge_Y1_hobo.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Dredge_Y1_hobo$DateTime <- mdy_hms(Dredge_Y1_hobo$Date_Time) #convert time field
Dredge_Y1_hobo <- Dredge_Y1_hobo[3:16531,] #cut 2 points in beginning, logger not yet in water
Dredge_Y1_hobo$Temp_K <- Dredge_Y1_hobo$Temp_C+273.15 #create collumn with temp in K
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour") #determine dates to aggregate around
AvgTempKbyhr <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_sub <- AvgTempKbyhr[676:4132,] #subset
T_field <- approxfun(x = c(0:3456), y = c(AvgTempKbyhr_sub$x), method = "linear", rule = 2) #the temp forcing function
T_Dredge2_Y1 <- T_field(0:3456) #for ease of later plotting the temperature forcing
##############
Dredge2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(Dredge2_SA))
##############

#Wickford1
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" #only necessary to run this code on some computers
Wickford_WSA <- filter(WSA2_Y1, Site == "Wickford") #filter by site
daily <- seq(as.Date("2017-12-4"), as.Date("2018-04-21"), by="days") #daily range relevant for this site
N <- Wickford_WSA[c("Date","NitrateNitrite_uM")] #new data frame with the relevant collumns
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/L
#multiplied by 24 to take the values from their daily to hourly positions
N_field <- approxfun(x = c(115*24, 0, 81*24, 59*24, 38*24, 138*24), y = c(N$NitrateNitrite_uM[1], N$NitrateNitrite_uM[3:7]), method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import literature DIC data
names(Segarra2002Carbon)[1] <- "Date" #only necessary from some computers to run this code
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date) #convert date collumn
CO_2 <- mean(Segarra2002Carbon$TCO2_micromolPERkg)/1000000 #(mol CO2/L)
###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_Wickfordy1 <-  NOAA_Irradiance$PAR[2702:3806] #subset by seq(as_datetime("2017-12-4 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3312, by = 3), y = NOAA_Irradiance_Wickfordy1, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
Wickford_Y1_hobo <- read.csv("Wickford_Y1_hobo.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Wickford Hobo data
Wickford_Y1_hobo$DateTime <- mdy_hms(Wickford_Y1_hobo$DateTime) #convert time field
Wickford_Y1_hobo <- Wickford_Y1_hobo[607:13839,] #subset based on 2017-12-04 15:30:00 to 2018-04-21 11:30:00 
Wickford_Y1_hobo$Temp_K <- Wickford_Y1_hobo$Temp_C+273.15 #create collumn with temp in K
WickfordT_hourly <- ceiling_date(Wickford_Y1_hobo$DateTime, unit = "hour") #determine values to aggregate around
AvgTempKbyhr <- aggregate(Wickford_Y1_hobo$Temp_K, by=list(WickfordT_hourly), mean) #calculate average hourly temp
fd <- AvgTempKbyhr[1:4,] #a few replacement data points at the front of the forcing
T_field <- approxfun(x = c(0:3312), y = c(fd$x, AvgTempKbyhr$x), method = "linear", rule = 2) #the temp forcing function
T_Wickford1_Y1 <- T_field(0:3312) #For ease of plotting the temp forcing
##############
Wickford1_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(Wickford1_SA))
##############

#RomePt1
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" #only necessary for some computers to run this code
RomePt_WSA <- filter(WSA2_Y1, Site == "Rome Point") #filter by site
daily <- seq(as.Date("2017-11-1"), as.Date("2018-04-21"), by="days") #days kelp in the water at this site
N <- RomePt_WSA[c("Date","NitrateNitrite_uM")] #new dataframe with relevant collumns
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/Lmean
#multiplied by 24 to take from daily to hourly
N_field <- approxfun(x = c(114*24, 148*24, 0, 71*24, 92*24, 171*24), y = c(N$NitrateNitrite_uM[1:3], N$NitrateNitrite_uM[5:7]), method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import literature DIC data
names(Segarra2002Carbon)[1] <- "Date" #only necessary for some computers to write this code
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date) #convert date collumn for easier use
CO_2 <- mean(Segarra2002Carbon$TCO2_micromolPERkg)/1000000 #(mol CO2/L)

###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_RomePty1 <-  NOAA_Irradiance$PAR[2438:3806] #subset by seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 4104, by = 3), y = NOAA_Irradiance_RomePty1, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
RomePoint_Y1_hobotemp <- read.csv("RomePoint_Y1_hobotemp.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
RomePoint_Y1_hobotemp$DateTime <- mdy_hms(RomePoint_Y1_hobotemp$DateTime) #convert time field
RomePoint_Y1_hobotemp <- RomePoint_Y1_hobotemp[6:16425,] #subset based on 2017-11-01 13:15:00 start and 2018-04-21 14:00:00 end
RomePoint_Y1_hobotemp$Temp_K <- RomePoint_Y1_hobotemp$Temp_C+273.15 #create collumn with temp in K
RomePointT_hourly <- ceiling_date(RomePoint_Y1_hobotemp$DateTime, unit = "hour") #determine dates to aggregate around
AvgTempKbyhr <- aggregate(RomePoint_Y1_hobotemp$Temp_K, by=list(RomePointT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[1:4103,] #subset
fd <- AvgTempKbyhr[1:2,] #two points of simulated data
T_field <- approxfun(x = c(0:4104), y = c(fd$x, AvgTempKbyhr$x), method = "linear", rule = 2) #the temp forcing function
T_RomePt1_Y1 <- T_field(0:4104) #for ease in plotting the temperature forcing
##############
RomePt1_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(RomePt1_SA))
##############

#RomePt2
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" #only necessary for some computers to run this code
RomePt_WSA <- filter(WSA2_Y1, Site == "Rome Point") #filter by site
daily <- seq(as.Date("2017-12-6"), as.Date("2018-04-21"), by="days") #days kelp in the water at this site
N <- RomePt_WSA[c("Date","NitrateNitrite_uM")] #new dataframe with relevant collumns
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/L
#multiplied by 24 to convert from daily positioning to hourly positioning
N_field <- approxfun(x = c(79*24, 113*24, -25*24, 57*24, 136*24), y = c(N$NitrateNitrite_uM[1:2], N$NitrateNitrite_uM[5:7]), method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import literature DIC data
names(Segarra2002Carbon)[1] <- "Date" #only necessary for this code to run in some computers
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date) #convert date collumn
CO_2 <- mean(Segarra2002Carbon$TCO2_micromolPERkg)/1000000 #(mol CO2/L)
###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_RomePty1_L2 <-  NOAA_Irradiance$PAR[2718:3806] #subset by seq(as_datetime("2017-12-6 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3264, by = 3), y = NOAA_Irradiance_RomePty1_L2, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
RomePoint_Y1_hobotemp <- read.csv("RomePoint_Y1_hobotemp.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
RomePoint_Y1_hobotemp$DateTime <- mdy_hms(RomePoint_Y1_hobotemp$DateTime) #convert time field
RomePoint_Y1_hobotemp <- RomePoint_Y1_hobotemp[3313:16425,] #subset based on 2017-12-06 00:00:00 - 2018-04-21 14:00:00
RomePoint_Y1_hobotemp$Temp_K <- RomePoint_Y1_hobotemp$Temp_C+273.15 #create collumn with temp in K
RomePointT_hourly <- ceiling_date(RomePoint_Y1_hobotemp$DateTime, unit = "hour") #detertime times to aggregate around
AvgTempKbyhr <- aggregate(RomePoint_Y1_hobotemp$Temp_K, by=list(RomePointT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[13:3277, ] #subset
T_field <- approxfun(x = c(0:3264), y = c(AvgTempKbyhr$x), method = "linear", rule = 2) #the temp forcing function
T_RomePt2_Y1 <- T_field(0:3264) #for ease in later plotting
##############
RomePt2_SA <- sensFun(SenseFunc, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(RomePt2_SA))
##############

#Sled1_Y2
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #only necessary to run this code on some computers
Sled_WSA <- filter(WSA_Y2, Site == "Moonstone Sled") #filter by site
daily <- seq(as.Date("2018-12-12"), as.Date("2019-05-03"), by="days") #date range kelp in water
Sled_WSA$NO3NO2_µM <- Sled_WSA$NO3NO2_µM/1000000 #convert from micromoles/L to moles/L
#multiplied by 24 to go from daily to hourly
N_field <- approxfun(x = c(1*24, 57*24, 93*24, 124*24, 142*24, 163*24), y = c(Sled_WSA$NO3NO2_µM), method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
Sled_DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
CO_2 <- CO_2/1000000 #need units to match K_C (molDIC/L)
###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_Sledy2 <-  NOAA_Irradiance$PAR[5686:6822] #subset by seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3408, by = 3), y = NOAA_Irradiance_Sledy2, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
Sled_Y2_Hobo <- read.csv("Sled_Y2_HoboLightTemp.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Sled_Y2_Hobo$DateTime <- mdy_hms(Sled_Y2_Hobo$DateTime) #convert date time field
Sled_Y2_Hobo$Temp_K <- Sled_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
SledY2T_hourly <- ceiling_date(Sled_Y2_Hobo$DateTime, unit = "hour") #determine times to aggregate around
AvgTempKbyhr <- aggregate(Sled_Y2_Hobo$Temp_K, by=list(SledY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[2:3385,] #subset
fd <- rep(285, 25) #small section of replacement
T_field <- approxfun(x = c(0:3408), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Sled1_Y2 <- T_field(0:3408) #for ease in plotting the temperature forcing
##############
Sled1_Y2_SA <- sensFun(SenseFunc2, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(Sled1_Y2_SA))
##############

#Sled2_Y2
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #only necessary for some computers to run this code
Sled_WSA <- filter(WSA_Y2, Site == "Moonstone Sled") #filter by site
daily <- seq(as.Date("2019-02-06"), as.Date("2019-05-03"), by="days") #days the kelp was in the water
Sled_WSA$NO3NO2_µM <- Sled_WSA$NO3NO2_µM/1000000 #convert from micromoles/L to moles/L
Sled_WSA <- Sled_WSA[2:6,] #remove the point before the relevant time range
N_field <- approxfun(x = c(1*24, 37*24, 68*24, 86*24, 107*24), y = c(Sled_WSA$NO3NO2_µM), method = "linear", rule = 2) #N forcing function

###### DIC forcing set-up ###########
Sled_DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
CO_2 <- CO_2/1000000 #need units to match K_C (molDIC/L)
###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_Sledy2_L2 <-  NOAA_Irradiance$PAR[6134:6822] #subset by seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 2064, by = 3), y = NOAA_Irradiance_Sledy2_L2, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
Sled_Y2_Hobo <- read.csv("Sled_Y2_HoboLightTemp.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Sled_Y2_Hobo$DateTime <- mdy_hms(Sled_Y2_Hobo$DateTime) #convert date time field
Sled_Y2_Hobo$Temp_K <- Sled_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
SledY2T_hourly <- ceiling_date(Sled_Y2_Hobo$DateTime, unit = "hour") #determine times to aggregate around
AvgTempKbyhr <- aggregate(Sled_Y2_Hobo$Temp_K, by=list(SledY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[1346:3385,] #subset
fd <- rep(285, 25) #small data replacement
T_field <- approxfun(x = c(0:2064), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Sled2_Y2 <- T_field(0:2064) #for ease in later plotting
##############
Sled2_Y2_SA <- sensFun(SenseFunc2, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(Sled2_Y2_SA))
##############

#Dredge1_Y2
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #only necessary for some computers to run this code
Dredge_WSA <- filter(WSA_Y2, Site == "Moonstone Dredge") #filter by site
daily <- seq(as.Date("2018-12-12"), as.Date("2019-05-03"), by="days") #growth date range
Dredge_WSA$NO3NO2_µM <- Dredge_WSA$NO3NO2_µM/1000000 #convert from micromoles/L to moles/L
#multiply by 24 to go from daily to hourly
N_field <- approxfun(x = c(1*24, 93*24, 124*24, 142*24, 163*24), y = c(Dredge_WSA$NO3NO2_µM), method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
Sled_DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
CO_2 <- CO_2/1000000 #need units to match K_C (molDIC/L)
###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_Dredgey2 <-  NOAA_Irradiance$PAR[5686:6822] #subset by seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3408, by = 3), y = NOAA_Irradiance_Dredgey2, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
Dredge_Y2_Hobo <- read.csv("Dredge_Y2_HoboTempLight.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Dredge_Y2_Hobo$DateTime <- mdy_hms(Dredge_Y2_Hobo$DateTime) #convert date time field
Dredge_Y2_Hobo$Temp_K <- Dredge_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
DredgeY2T_hourly <- ceiling_date(Dredge_Y2_Hobo$DateTime, unit = "hour") #determine what times to aggregate around
AvgTempKbyhr <- aggregate(Dredge_Y2_Hobo$Temp_K, by=list(DredgeY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[2:3384,] #subset
fd <- rep(285, 26) #small amount of replacement data
T_field <- approxfun(x = c(0:3408), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Dredge1_Y2 <- T_field(0:3408) #for ease of plotting
##############
Dredge1_Y2_SA <- sensFun(SenseFunc2, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(Dredge1_Y2_SA))
##############

#Dredge2_Y2
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #only necessary for some computers to run this code
Dredge_WSA <- filter(WSA_Y2, Site == "Moonstone Dredge") #filter by site
daily <- seq(as.Date("2019-02-06"), as.Date("2019-05-03"), by="days") #date range kelp at farm
Dredge_WSA$NO3NO2_µM <- Dredge_WSA$NO3NO2_µM/1000000 #convert from micromoles/L to moles/L
#multiply by 24 to switch from daily to hourly
N_field <- approxfun(x = c(-25*24, 37*24, 68*24, 86*24, 107*24), y = c(Dredge_WSA$NO3NO2_µM), method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
Sled_DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Ninigret DIC data
CO_2 <- mean(Sled_DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
CO_2 <- CO_2/1000000 #need units to match K_C (molDIC/L)
###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_Dredgey2_L2 <-  NOAA_Irradiance$PAR[6134:6822] #subset by seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 2064, by = 3), y = NOAA_Irradiance_Dredgey2_L2, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
Dredge_Y2_Hobo <- read.csv("Dredge_Y2_HoboTempLight.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Dredge_Y2_Hobo$DateTime <- mdy_hms(Dredge_Y2_Hobo$DateTime) #convert date time field
Dredge_Y2_Hobo$Temp_K <- Dredge_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
DredgeY2T_hourly <- ceiling_date(Dredge_Y2_Hobo$DateTime, unit = "hour") #determine values to aggregate around
AvgTempKbyhr <- aggregate(Dredge_Y2_Hobo$Temp_K, by=list(DredgeY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[1346:3385,] #subset
fd <- rep(285, 25) #estimation to fill in gap in data
T_field <- approxfun(x = c(0:2064), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Dredge2_Y2 <- T_field(0:2064) #for ease in plotting
##############
Dredge2_Y2_SA <- sensFun(SenseFunc2, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(Dredge2_Y2_SA))
##############

#Wickford1_Y2
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #only necessary for some computers to run this code
Wickford_WSA <- filter(WSA_Y2, Site == "Wickford") #filter by site
daily <- seq(as.Date("2018-12-19"), as.Date("2019-05-23"), by="days") #the date range for the kelp in the field
Wickford_WSA$NO3NO2_µM <- Wickford_WSA$NO3NO2_µM/1000000 #convert from micromoles/L to moles/L
#multiply by 24 to convert daily values to hourly values
N_field <- approxfun(x = c(1*24, 55*24, 85*24, 156*24), y = c(Wickford_WSA$NO3NO2_µM), method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import lit TCO2 data
names(Segarra2002Carbon)[1] <- "Date" #Only necessary for running the code on some computer
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date) #time field conversion
CO_2 <- mean(Segarra2002Carbon$TCO2_micromolPERkg)/1000000 #(mol CO2/L)
###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_Wickfordy2 <-  NOAA_Irradiance$PAR[5742:6982] #subset by seq(as_datetime("2018-12-19 12:00:00"), as_datetime("2019-05-23 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3720, by = 3), y = NOAA_Irradiance_Wickfordy2, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
Wickford_Y2_Hobo <- read.csv("Wickford_Y2_HoboLightTemp.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Wickford_Y2_Hobo$DateTime <- mdy_hms(Wickford_Y2_Hobo$DateTime) #convert date time field
Wickford_Y2_Hobo$Temp_K <- Wickford_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
WickfordY2T_hourly <- ceiling_date(Wickford_Y2_Hobo$DateTime, unit = "hour") #determine the times to aggregate around
AvgTempKbyhr <- aggregate(Wickford_Y2_Hobo$Temp_K, by=list(WickfordY2T_hourly), mean) #calculate average hourly temp
fd <- rep(278, 4) #replacement data
AvgTempKbyhr_sub <- AvgTempKbyhr[4:3716,] #subset
fd2 <- rep(287, 4) #second small section of replacement data
T_field <- approxfun(x = c(0:3720), y = c(fd, AvgTempKbyhr_sub$x, fd2), method = "linear", rule = 2) #the temp forcing function
T_Wickford1_Y2 <- T_field(0:3720) #for ease of plotting
##############
Wickford1_Y2_SA <- sensFun(SenseFunc2, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(Wickford1_Y2_SA))
##############

#RomePt1_Y2
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #only necessary for some computers to run this code
RomePt_WSA <- filter(WSA_Y2, Site == "Rome Point") #filter by site
daily <- seq(as.Date("2018-12-20"), as.Date("2019-05-24"), by="days") #the date range for kelp on the farm
RomePt_WSA$NO3NO2_µM <- RomePt_WSA$NO3NO2_µM/1000000 #convert from micromoles/L to moles/L
#multiplied by 24 to take daily values to hourly values
N_field <- approxfun(x = c(1*24, 64*24, 84*24, 155*24), y = c(RomePt_WSA$NO3NO2_µM), method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import literature DIC data
names(Segarra2002Carbon)[1] <- "Date" #Only necessary to run this code on some computers
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date) #convert date field
CO_2 <- mean(Segarra2002Carbon$TCO2_micromolPERkg)/1000000 #(mol CO2/L)
###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_RomePty2 <-  NOAA_Irradiance$PAR[5750:6990] #subset by seq(as_datetime("2018-12-20 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3720, by = 3), y = NOAA_Irradiance_RomePty2, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
RomePt_Y2_Hobo <- read.csv("RomePt_Y2_HoboTempLight.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
RomePt_Y2_Hobo$DateTime <- mdy_hm(RomePt_Y2_Hobo$DateTime) #convert date time field
RomePt_Y2_Hobo$Temp_K <- RomePt_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
RomePtY2T_hourly <- ceiling_date(RomePt_Y2_Hobo$DateTime, unit = "hour") #determine the times to aggregate around
AvgTempKbyhr <- aggregate(RomePt_Y2_Hobo$Temp_K, by=list(RomePtY2T_hourly), mean) #calculate average hourly temp
fd <- rep(280, 4) #small bit of replacement data
AvgTempKbyhr_sub <- AvgTempKbyhr[28:2414,] #subset
#Using Wickford temp to fill in th gap in the Rome Pt temp
Wickford_Y2_Hobo <- read.csv("Wickford_Y2_HoboLightTemp.csv", header = TRUE, fileEncoding="UTF-8-BOM")
Wickford_Y2_Hobo$DateTime <- mdy_hms(Wickford_Y2_Hobo$DateTime)
Wickford_Y2_Hobo$Temp_K <- Wickford_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
WickfordY2T_hourly <- ceiling_date(Wickford_Y2_Hobo$DateTime, unit = "hour")
AvgTempKbyhr_Wickford <- aggregate(Wickford_Y2_Hobo$Temp_K, by=list(WickfordY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_W <- AvgTempKbyhr_Wickford[2415:3716,]
fd2 <- rep(287, 28)
T_field <- approxfun(x = c(0:3720), y = c(fd, AvgTempKbyhr_sub$x, AvgTempKbyhr_W$x, fd2), method = "linear", rule = 2) #the temp forcing function
T_RomePt1_Y2 <- T_field(0:3720) #for ease of plotting
##############
RomePt1_Y2_SA <- sensFun(SenseFunc2, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(RomePt1_Y2_SA))
##############

#RomePt2_Y2
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)
###### N forcing set-up##############
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #only necessary for some computers to run this code
RomePt_WSA <- filter(WSA_Y2, Site == "Rome Point") #filter by site
daily <- seq(as.Date("2019-2-21"), as.Date("2019-05-24"), by="days") #date range for the kelp on the farm
RomePt_WSA$NO3NO2_µM <- RomePt_WSA$NO3NO2_µM/1000000 #convert from micromoles/L to moles/L
RomePt_WSA <- RomePt_WSA[2:4,] #subset based on what data is relevant here
#multiplying by 24 to converty from daily values to hourly values
N_field <- approxfun(x = c(1*24, 21*24, 92*24), y = c(RomePt_WSA$NO3NO2_µM), method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import lit TCO2 data
names(Segarra2002Carbon)[1] <- "Date" #Only necessary to run this code on some computers
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date) #convert date values
CO_2 <- mean(Segarra2002Carbon$TCO2_micromolPERkg)/1000000 #(mol CO2/L)
###### NOAA Irradiance forcing set-up ####
NOAA_Irradiance_RomePty2_L2 <-  NOAA_Irradiance$PAR[6254:6990] #subset by seq(as_datetime("2019-2-21 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 2208, by = 3), y = NOAA_Irradiance_RomePty2_L2, method = "linear", rule = 2) #irradiance forcing function
###### Temp forcing set-Up #############
RomePt_Y2_Hobo <- read.csv("RomePt_Y2_HoboTempLight.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
RomePt_Y2_Hobo$DateTime <- mdy_hm(RomePt_Y2_Hobo$DateTime) #convert date time
RomePt_Y2_Hobo$Temp_K <- RomePt_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
RomePtY2T_hourly <- ceiling_date(RomePt_Y2_Hobo$DateTime, unit = "hour") #determine what times to aggregate around
AvgTempKbyhr <- aggregate(RomePt_Y2_Hobo$Temp_K, by=list(RomePtY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_sub <- AvgTempKbyhr[1536:2414,] #subset
#Using Wickford temp to fill in th gap in the Rome Pt temp
Wickford_Y2_Hobo <- read.csv("Wickford_Y2_HoboLightTemp.csv", header = TRUE, fileEncoding="UTF-8-BOM")
Wickford_Y2_Hobo$DateTime <- mdy_hms(Wickford_Y2_Hobo$DateTime)
Wickford_Y2_Hobo$Temp_K <- Wickford_Y2_Hobo$Temp_C+273.15 #create collumn with temp in K
WickfordY2T_hourly <- ceiling_date(Wickford_Y2_Hobo$DateTime, unit = "hour")
AvgTempKbyhr_Wickford <- aggregate(Wickford_Y2_Hobo$Temp_K, by=list(WickfordY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_W <- AvgTempKbyhr_Wickford[2415:3716,]
fd2 <- rep(287, 28)
T_field <- approxfun(x = c(0:2208), y = c(AvgTempKbyhr_sub$x, AvgTempKbyhr_W$x, fd2), method = "linear", rule = 2) #the temp forcing function
T_RomePt2_Y2 <- T_field(0:2208) #for ease in plotting
##############
RomePt2_Y2_SA <- sensFun(SenseFunc2, params_Lo, sensvar = c("m_EC", "m_EN", "M_V"))
#plot(summary(RomePt2_Y2_SA))
##############

Compile_SA <- rbind(Sled1_SA, Sled2_SA, Dredge1_SA, Dredge2_SA, Wickford1_SA, RomePt1_SA, RomePt2_SA, Sled1_Y2_SA, Sled2_Y2_SA, Dredge1_Y2_SA, Dredge2_Y2_SA, Wickford1_Y2_SA, RomePt1_Y2_SA, RomePt2_Y2_SA)
summary(Compile_SA)
#plot(summary(Compile_SA))
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
#SA_SUM$Parameter[SA_SUM$Parameter=="k_I"] <- "Dissociation rate"
#SA_SUM$Parameter[SA_SUM$Parameter=="y_I_C"] <- "Yield factor of C reserve to photons"
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

#Figure 10
ggplot(data = SA_SUM) +
  geom_point(mapping = aes(x = L1, y = Parameter), size = 3) +
  geom_abline(mapping = aes(slope = 0, intercept = 0)) +
  xlim(0, 8000) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "L1 (sensitivity function)", y = "Parameter Name")

