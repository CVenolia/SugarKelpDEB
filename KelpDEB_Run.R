####################################################################################################################
#Run file of Sugar Kelp model using data from Sled line 1 and 2 and Dredge line 1 and 2 (two sites at Point Judith Pond, Moonstone/Cedar Island Oysters Farm)
#and 1 line at Wickford and 2 lines at Rome Point (both sites in Narragansett Bay)
#Sled = Pt Judith Pond N
#Dredge = Pt Judith Pond S
#Wickford = Narragansett Bay N
#Rome Point = Narragansett Bay S

#File created by Celeste Venolia in March 2018-December 2019, please contact celestevenolia@gmail.com with questions

####################################################################################################################

#Reminder to set working directory to location of data

#Import libraries
library(deSolve)
library(tidyverse)
library(lubridate)
library(gridExtra)
library(gdata)
library(Metrics)
library(pracma)

#Required for model runs
source("SolveR_R.R")
source("KelpDEB_model.R")
#Required for Calibration Code
source("N_uptake_Calibration.R")
source("Photosynthesis_Calibration.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#############
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#YEAR 1
#
#Setting up the forcing functions with field data for Sled line 1
B <- 0.003006 #inital biomass for conversions (cannot put in initial conditions)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Sled1 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Sled line 2
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Sled2 <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Dredge line 1
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Dredge1 <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Dredge line 2
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
################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Dredge2 <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Wickford line 1
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
################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Wickford1 <- ode(y = state_Lo, t = times_Lo_Wickford1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Rome Point line 1
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
################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_RomePt1 <- ode(y = state_Lo, t = times_Lo_RomePt1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Rome Point line 2
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
################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_RomePt2 <- ode(y = state_Lo, t = times_Lo_RomePt2, func = rates_Lo, parms = params_Lo)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#YEAR 2 Kelp data
#Setting up the forcing functions with field data for Sled line 1 (y2)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Sled1_Y2 <- ode(y = state_Lo, t = times_Y2_Sled1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Sled line 2 (y2)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Sled2_Y2 <- ode(y = state_Lo, t = times_Y2_Sled2, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Dredge line 1 (y2)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Dredge1_Y2 <- ode(y = state_Lo, t = times_Y2_Dredge1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Dredge line 2 (y2)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Dredge2_Y2 <- ode(y = state_Lo, t = times_Y2_Dredge2, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Wickford line 1 (y2)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Wickford1_Y2 <- ode(y = state_Lo, t = times_Y2_Wickford1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Rome Point line 1 (y2)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_RomePt1_Y2 <- ode(y = state_Lo, t = times_Y2_RomePt1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Rome Point line 2 (y2)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_RomePt2_Y2 <- ode(y = state_Lo, t = times_Y2_RomePt2, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#END FIELD DATA, START LITERATURE DATA FOR CALIBRATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#N uptake first (skipping Subandar et al. (1993), because their N concentrations are not well documented)
#Setting up the forcing functions for Espinoza and Chapman (1983) nitrogen uptake #9C
###### N forcing set-up##############
Nmax <- 73.1221719/1000000 #M

###### Temp forcing set-Up #############
#9 C is 282.15 K
T_dat <- 9 #C (conversion in Nuptake function to K)
###################################
#Model run (the differential equation solver)
sol_EspinozaChapman1983_N_9 <- Nuptake(params_Lo, T_dat, Nmax, w_EN)
sol_EspinozaChapman1983_N_9 <- as.data.frame(sol_EspinozaChapman1983_N_9)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions for Espinoza and Chapman (1983) nitrogen uptake #18C
###### N forcing set-up##############
Nmax <- 76.9543147/1000000 #M

###### Temp forcing set-Up #############
#18 C is 282.15 K
T_dat <- 18 #C (conversion in Nuptake function to K)
#################################
#Model run (the differential equation solver)
sol_EspinozaChapman1983_N_18 <- Nuptake(params_Lo, T_dat, Nmax, w_EN)
sol_EspinozaChapman1983_N_18 <- as.data.frame(sol_EspinozaChapman1983_N_18)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Photosynthesis related
#forcings
##### Temp forcing set-up ####
T_dat <- 14 #C (maintained for entire experiment
###### I forcing set-up #####
I_max <- 3233205 #micromol photons m-2 h-1
##############
sol_Johansson2002 <- Photosynthesis(params_Lo, cond, w_V, w_EN, w_EC, w_O2, T_dat, I_max)
sol_Johansson2002 <- as.data.frame(sol_Johansson2002)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###### Convert DeSolve solutions into data frame for broader plotting use ####
  ##### Year 1 #####
sol_Sled1 <- as.data.frame(sol_Sled1)
sol_Sled2 <- as.data.frame(sol_Sled2)
sol_Dredge1 <- as.data.frame(sol_Dredge1)
sol_Dredge2 <- as.data.frame(sol_Dredge2)
sol_Wickford1 <- as.data.frame(sol_Wickford1)
sol_RomePt1 <- as.data.frame(sol_RomePt1)
sol_RomePt2 <- as.data.frame(sol_RomePt2)

sol_Sled1$Date <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
sol_Sled2$Date <- seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
sol_Dredge1$Date <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
sol_Dredge2$Date <-seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
sol_Wickford1$Date <- seq(as_datetime("2017-12-4 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
sol_RomePt1$Date <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
sol_RomePt2$Date  <- seq(as_datetime("2017-12-6 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")

sol_Sled1$Temp_C <- T_Sled1_Y1 - 273.15
sol_Sled2$Temp_C <- T_Sled2_Y1 - 273.15
sol_Dredge1$Temp_C <- T_Dredge1_Y1 - 273.15
sol_Dredge2$Temp_C <- T_Dredge2_Y1 - 273.15
sol_Wickford1$Temp_C <- T_Wickford1_Y1 - 273.15
sol_RomePt1$Temp_C <- T_RomePt1_Y1 - 273.15
sol_RomePt2$Temp_C <- T_RomePt2_Y1 - 273.15

sol_Sled1$source <- "sol_Sled1"
sol_Sled2$source <- "sol_Sled2"
sol_Dredge1$source <- "sol_Dredge1"
sol_Dredge2$source <- "sol_Dredge2"
sol_Wickford1$source <- "sol_Wickford1"
sol_RomePt1$source <- "sol_RomePt1"
sol_RomePt2$source  <- "sol_RomePt2"

#combine all Y1 field data into one dataframe
sol_all <- rbind(sol_Dredge1, sol_Dredge2, sol_RomePt1, sol_RomePt2, sol_Sled1, sol_Sled2, sol_Wickford1)
#sol_all <- combine(sol_Dredge1, sol_Dredge2, sol_RomePt1, sol_RomePt2, sol_Sled1, sol_Sled2, sol_Wickford1)

#sol_all$Date <- as.Date(sol_all$Date)
sol_all$source <- as.character(sol_all$source)
sol_all$source[sol_all$source == "sol_Dredge1"] <- "Point Judith Pond S 1"
sol_all$source[sol_all$source == "sol_Dredge2"] <- "Point Judith Pond S 2"
sol_all$source[sol_all$source == "sol_Sled1"] <- "Point Judith Pond N 1"
sol_all$source[sol_all$source == "sol_Sled2"] <- "Point Judith Pond N 2"
sol_all$source[sol_all$source == "sol_Wickford1"] <- "Narragansett Bay N 1"
sol_all$source[sol_all$source == "sol_RomePt1"] <- "Narragansett Bay S 1"
sol_all$source[sol_all$source == "sol_RomePt2"] <- "Narragansett Bay S 2"

  ###### Year 2 #####
sol_Sled1_Y2 <- as.data.frame(sol_Sled1_Y2)
sol_Sled2_Y2 <- as.data.frame(sol_Sled2_Y2)
sol_Dredge1_Y2 <- as.data.frame(sol_Dredge1_Y2)
sol_Dredge2_Y2 <- as.data.frame(sol_Dredge2_Y2)
sol_Wickford1_Y2 <- as.data.frame(sol_Wickford1_Y2)
sol_RomePt1_Y2 <- as.data.frame(sol_RomePt1_Y2)
sol_RomePt2_Y2 <- as.data.frame(sol_RomePt2_Y2)

sol_Sled1_Y2$Date <- seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour") 
sol_Sled2_Y2$Date <- seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
sol_Dredge1_Y2$Date <- seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
sol_Dredge2_Y2$Date <- seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
sol_Wickford1_Y2$Date <- seq(as_datetime("2018-12-19 12:00:00"), as_datetime("2019-05-23 12:00:00"), by="hour")
sol_RomePt1_Y2$Date <- seq(as_datetime("2018-12-20 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
sol_RomePt2_Y2$Date <- seq(as_datetime("2019-2-21 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")

sol_Sled1_Y2$Temp_C <- T_Sled1_Y2 - 273.15
sol_Sled2_Y2$Temp_C <- T_Sled2_Y2 - 273.15
sol_Dredge1_Y2$Temp_C <- T_Dredge1_Y2 - 273.15
sol_Dredge2_Y2$Temp_C <- T_Dredge2_Y2 - 273.15
sol_Wickford1_Y2$Temp_C <- T_Wickford1_Y2 - 273.15
sol_RomePt1_Y2$Temp_C <- T_RomePt1_Y2 - 273.15
sol_RomePt2_Y2$Temp_C <- T_RomePt2_Y2 - 273.15

sol_Sled1_Y2$source <- "sol_Sled1_Y2"
sol_Sled2_Y2$source <- "sol_Sled2_Y2"
sol_Dredge1_Y2$source <- "sol_Dredge1_Y2"
sol_Dredge2_Y2$source <- "sol_Dredge2_Y2"
sol_Wickford1_Y2$source <- "sol_Wickford1_Y2"
sol_RomePt1_Y2$source <- "sol_RomePt1_Y2"
sol_RomePt2_Y2$source  <- "sol_RomePt2_Y2"

#combine all Y2 field data into one dataframe
sol_all_Y2 <- rbind(sol_Dredge1_Y2, sol_Dredge2_Y2, sol_RomePt1_Y2, sol_RomePt2_Y2, sol_Sled1_Y2, sol_Sled2_Y2, sol_Wickford1_Y2)
#sol_all_Y2 <- combine(sol_Dredge1_Y2, sol_Dredge2_Y2, sol_RomePt1_Y2, sol_RomePt2_Y2, sol_Sled1_Y2, sol_Sled2_Y2, sol_Wickford1_Y2)

#sol_all_Y2$Date <- as.Date(sol_all_Y2$Date)
sol_all_Y2$source <- as.character(sol_all_Y2$source)
sol_all_Y2$source[sol_all_Y2$source == "sol_Dredge1_Y2"] <- "Point Judith Pond S 1"
sol_all_Y2$source[sol_all_Y2$source == "sol_Dredge2_Y2"] <- "Point Judith Pond S 2"
sol_all_Y2$source[sol_all_Y2$source == "sol_Sled1_Y2"] <- "Point Judith Pond N 1"
sol_all_Y2$source[sol_all_Y2$source == "sol_Sled2_Y2"] <- "Point Judith Pond N 2"
sol_all_Y2$source[sol_all_Y2$source == "sol_Wickford1_Y2"] <- "Narragansett Bay N 1"
sol_all_Y2$source[sol_all_Y2$source == "sol_RomePt1_Y2"] <- "Narragansett Bay S 1"
sol_all_Y2$source[sol_all_Y2$source == "sol_RomePt2_Y2"] <- "Narragansett Bay S 2"

##### Parameter Plots for field data #####
    #######YEAR 1 ####
# Simulations only
ggplot() + 
  geom_point(data = sol_all, aes(Date, L_allometric, shape = source)) +
  scale_shape_manual(values=1:nlevels(sol_all$source)) +
  labs(x= "Date", y = "Blade length (cm)") +
  ggtitle("Year 1 Kelp growth model output")

plot_CO_2 <- ggplot(data = sol_all, aes(Date, CO_2, color = source), position = "jitter") + 
  geom_point() +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  labs(x= "Date", y = "mol DIC L-1") +
  ggtitle("A")
plot_N <- ggplot(data = sol_all, aes(Date, N, color = source)) + 
  geom_point() +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  theme_bw() +
  theme(legend.position="none") + 
  labs(x= "Date", y = "mol NO3- and NH4+ L-1") +
  ggtitle("B")
plot_C_T <- ggplot(data = sol_all, aes(Date, C_T, color = source)) + 
  geom_point() +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.title = element_blank()) +
  labs(x= "Date", y = "Temperature Correction (no unit)")
  #ggtitle("C_T")
plot_I <- ggplot(data = sol_all, aes(Date, I, color = source)) + 
  geom_point() +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  theme_bw() +
  labs(x= "Date", y = "μE m-2 h-1") +
  ggtitle("Irradiance")
grid.arrange(plot_CO_2, plot_N, plot_C_T, plot_I, ncol=2)

grid.arrange(plot_CO_2, plot_N, ncol=2)

plot_J_CO2 <- ggplot(data = sol_all, aes(Date, J_CO2, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "molCO2 /molM_V/h") +
  ggtitle("J_CO2")
plot_J_I <- ggplot(data = sol_all, aes(Date, J_I, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "molNADPH/molV/h") +
  ggtitle("J_I")
plot_J_EC_A <- ggplot(data = sol_all, aes(Date, J_EC_A, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol C/mol M_V/h") +
  ggtitle("J_EC_A")
plot_J_EN_A <- ggplot(data = sol_all, aes(Date, J_EN_A, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol N/molM_V/h") +
  ggtitle("J_EN_A")
grid.arrange(plot_J_CO2, plot_J_I, plot_J_EC_A, plot_J_EN_A, ncol=2)

plot_J_EC_C <- ggplot(data = sol_all, aes(Date, J_EC_C, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EC/molM_V/h") +
  ggtitle("J_EC_C")
plot_J_EN_C <- ggplot(data = sol_all, aes(Date, J_EN_C, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EN/molM_V/h") +
  ggtitle("J_EN_C")
plot_J_EC_M <- ggplot(data = sol_all, aes(Date, J_EC_M, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EC/molM_V/h") +
  ggtitle("J_EC_M")
plot_J_EN_M <- ggplot(data = sol_all, aes(Date, J_EN_M, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EN/molM_V/h") +
  ggtitle("J_EN_M")
grid.arrange(plot_J_EC_C, plot_J_EN_C, plot_J_EC_M, plot_J_EN_M, ncol=2)

plot_J_EC_G <- ggplot(data = sol_all, aes(Date, J_EC_G, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EC/molM_V/h") +
  ggtitle("J_EC_G")
plot_J_EN_G <- ggplot(data = sol_all, aes(Date, J_EN_G, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EN/molM_V/h") +
  ggtitle("J_EN_G")
plot_J_G <- ggplot(data = sol_all, aes(Date, J_G, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "1/h") +
  ggtitle("J_G")
plot_r <- ggplot(data = sol_all, aes(Date, r, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "r (1/h)") +
  ggtitle("r")
grid.arrange(plot_J_EC_G, plot_J_EN_G, plot_J_G, plot_r, ncol=2)

plot_J_EC_R <- ggplot(data = sol_all, aes(Date, J_EC_R, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EC/molM_V/h") +
  ggtitle("J_EC_R")
plot_J_EN_R <- ggplot(data = sol_all, aes(Date, J_EN_R, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EN/molM_V/h") +
  ggtitle("J_EN_R")
plot_J_VM <- ggplot(data = sol_all, aes(Date, J_VM, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "1/h") +
  ggtitle("J_VM")
plot_NC <- ggplot(data = sol_all, aes(Date, NC, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol N/mol C") +
  ggtitle("N:C ratio")
grid.arrange(plot_J_EC_R, plot_J_EN_R, plot_J_VM, plot_NC, ncol=2)

plot_m_EC <- ggplot(data = sol_all, aes(Date, m_EC, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "m_EC (mol C/mol M_V)") +
  ggtitle("m_EC")
plot_m_EN <- ggplot(data = sol_all, aes(Date, m_EN, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "m_EN (mol N/mol M_V)") +
  ggtitle("m_EN")
plot_M_V <- ggplot(data = sol_all, aes(Date, M_V, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "M_V (mol)") +
  ggtitle("M_V")
plot_B <- ggplot(data = sol_all, aes(Date, B, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "B (g)") +
  ggtitle("B")
grid.arrange(plot_m_EC, plot_m_EN, plot_M_V, plot_B, ncol=2)

#Oxygen
ggplot(data = sol_all) + 
  geom_point(mapping = aes(x = Date, y = J_O)) +
  labs(x= "Date", y = "mgO2/g/d") +
  ggtitle("Rate of Oxygen Production")

    #######YEAR 2 #####
ggplot() + 
  geom_point(data = sol_all_Y2, aes(Date, L_allometric, shape = source)) +
  scale_shape_manual(values=1:nlevels(sol_all_Y2$source)) +
  labs(x= "Date", y = "Blade length (cm)") +
  ggtitle("Year 1 Kelp growth in comparison to model output")

plot_CO_2 <- ggplot(data = sol_all_Y2, aes(Date, CO_2, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "molDIC/L") +
  ggtitle("DIC forcing")
plot_N <- ggplot(data = sol_all_Y2, aes(Date, N, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "molNO3/L") +
  ggtitle("Nitrate forcing")
plot_C_T <- ggplot(data = sol_all_Y2, aes(Date, C_T, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "C_T (no unit)") +
  ggtitle("C_T")
plot_I <- ggplot(data = sol_all_Y2, aes(Date, I, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "μE m-2 d-1") +
  ggtitle("Irradiance")
grid.arrange(plot_CO_2, plot_N, plot_C_T, plot_I, ncol=2)

plot_J_CO2 <- ggplot(data = sol_all_Y2, aes(Date, J_CO2, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "molC/molM_V/d") +
  ggtitle("J_CO2")
plot_J_I <- ggplot(data = sol_all_Y2, aes(Date, J_I, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol C/mol M_V/d") +
  ggtitle("J_I")
plot_J_EC_A <- ggplot(data = sol_all_Y2, aes(Date, J_EC_A, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol C/mol M_V/d") +
  ggtitle("J_EC_A")
plot_J_EN_A <- ggplot(data = sol_all_Y2, aes(Date, J_EN_A, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "molNO3/molM_V/d") +
  ggtitle("J_EN_A")
grid.arrange(plot_J_CO2, plot_J_I, plot_J_EC_A, plot_J_EN_A, ncol=2)

plot_J_EC_C <- ggplot(data = sol_all_Y2, aes(Date, J_EC_C, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EC/molM_V/d") +
  ggtitle("J_EC_C")
plot_J_EN_C <- ggplot(data = sol_all_Y2, aes(Date, J_EN_C, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EN/molM_V/d") +
  ggtitle("J_EN_C")
plot_J_EC_M <- ggplot(data = sol_all_Y2, aes(Date, J_EC_M, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EC/molM_V/d") +
  ggtitle("J_EC_M")
plot_J_EN_M <- ggplot(data = sol_all_Y2, aes(Date, J_EN_M, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EN/molM_V/d") +
  ggtitle("J_EN_M")
grid.arrange(plot_J_EC_C, plot_J_EN_C, plot_J_EC_M, plot_J_EN_M, ncol=2)

plot_J_EC_G <- ggplot(data = sol_all_Y2, aes(Date, J_EC_G, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EC/molM_V/d") +
  ggtitle("J_EC_G")
plot_J_EN_G <- ggplot(data = sol_all_Y2, aes(Date, J_EN_G, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EN/molM_V/d") +
  ggtitle("J_EN_G")
plot_J_G <- ggplot(data = sol_all_Y2, aes(Date, J_G, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "1/d") +
  ggtitle("J_G")
plot_r <- ggplot(data = sol_all_Y2, aes(Date, r, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "r (1/d)") +
  ggtitle("r")
grid.arrange(plot_J_EC_G, plot_J_EN_G, plot_J_G, plot_r, ncol=2)

plot_J_EC_R <- ggplot(data = sol_all_Y2, aes(Date, J_EC_R, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EC/molM_V/d") +
  ggtitle("J_EC_R")
plot_J_EN_R <- ggplot(data = sol_all_Y2, aes(Date, J_EN_R, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol EN/molM_V/d") +
  ggtitle("J_EN_R")
plot_J_VM <- ggplot(data = sol_all_Y2, aes(Date, J_VM, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "1/d") +
  ggtitle("J_VM")
plot_NC <- ggplot(data = sol_all_Y2, aes(Date, NC, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "mol N/mol C") +
  ggtitle("N:C ratio")
grid.arrange(plot_J_EC_R, plot_J_EN_R, plot_J_VM, plot_NC, ncol=2)

plot_m_EC <- ggplot(data = sol_all_Y2, aes(Date, m_EC, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "m_EC (mol C/mol M_V)") +
  ggtitle("m_EC")
plot_m_EN <- ggplot(data = sol_all_Y2, aes(Date, m_EN, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "m_EN (mol N/mol M_V)") +
  ggtitle("m_EN")
plot_M_V <- ggplot(data = sol_all_Y2, aes(Date, M_V, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "M_V (mol)") +
  ggtitle("M_V")
plot_B <- ggplot(data = sol_all_Y2, aes(Date, B, color = source)) + 
  geom_point() +
  scale_color_viridis_d() +
  labs(x= "Date", y = "B (g)") +
  ggtitle("B")
grid.arrange(plot_m_EC, plot_m_EN, plot_M_V, plot_B, ncol=2)

#Oxygen
ggplot(data = sol_all_Y2) + 
  geom_point(mapping = aes(x = Date, y = J_O)) +
  labs(x= "Date", y = "mgO2/g/d") +
  ggtitle("Rate of Oxygen Production")
    #######Combined Figures for Pub ####
plot_T_Y1 <- ggplot(data = sol_all, aes(Date, Temp_C, color = source)) + 
  geom_point() +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  ylim(-2, 25) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2017-2018)", y = "Temperature (°C)") +
  ggtitle("A")
plot_T_Y2 <- ggplot(data = sol_all_Y2, aes(Date, Temp_C, color = source)) + 
  geom_point() +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  ylim(-2, 25) +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2018-2019)", y = "Temperature (°C)") +
  ggtitle("B")
plot_C_T <- ggplot(data = sol_all, aes(Date, C_T, color = source)) + 
  geom_point() +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2017-2018)", y = "Temperature Correction (no unit)") +
  ggtitle("C")
plot_C_T_Y2 <- ggplot(data = sol_all_Y2, aes(Date, C_T, color = source)) + 
  geom_point() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2018-2019)", y = "Temperature Correction (no unit)") +
  ggtitle("D")
grid.arrange(plot_T_Y1, plot_T_Y2, plot_C_T, plot_C_T_Y2, ncol=2)

plot_N <- ggplot(data = sol_all, aes(Date, N, color = source)) + 
  geom_line(size = 2) +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  ylim(0, 1.25e-05) +
  labs(x= "Date (2017-2018)", y = bquote('mol' ~NO[3]^{"-"}~ 'and' ~NH[4]^{"+"}~ 'L'^"-1")) +
  ggtitle("A")
plot_N_Y2 <- ggplot(data = sol_all_Y2, aes(Date, N, color = source)) + 
  geom_line(size = 2) +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  ylim(0, 1.25e-05) +
  labs(x= "Date (2018-2019)", y = bquote('mol' ~NO[3]^{"-"}~ 'and' ~NH[4]^{"+"}~ 'L'^"-1")) +
  ggtitle("B")
grid.arrange(plot_N, plot_N_Y2, ncol=2)

plot_I_NBN <- ggplot() + 
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, I), color = "gray0") +
  #geom_point(data = Ibyhr_Wickford, aes(Group.1, x), color = "gray60") +
  theme_bw() +
  ylim(0, 4500000) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  xlab("Date (2017-2018)") +
  ylab(bquote('μE m'^"-2"*' h'^"-1")) + #"μE m-2 d-1"
  ggtitle("Narragansett Bay N")
RomePt_Irradiance <- read.csv("RomePt_Calibrated_Irradiance.csv", header = TRUE)
names(RomePt_Irradiance)[2] <- "Date" # RL: I had to use this because the first column name was different (probably due to some encoding glitch)
RomePt_Irradiance$datetime <- paste(RomePt_Irradiance$Date, RomePt_Irradiance$Time)
RomePt_Irradiance$datetime <- dmy_hms(RomePt_Irradiance$datetime)
R_hourly <- ceiling_date(RomePt_Irradiance$datetime, unit = "hour")
Ibyhr_RomePt <- aggregate(RomePt_Irradiance$Irradiance, by=list(R_hourly), sum)

plot_I_NBS <- ggplot() + 
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, I), color = "gray0") +
  #all replaced with Wickford light data
  #geom_point(data = Ibyhr_RomePt, aes(Group.1, x), color = "gray60") +
  theme_bw() +
  ylim(0, 4500000) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2017-2018)", y = bquote('μE m'^"-2"*' h'^"-1")) +
  ggtitle("Narragansett Bay S")
plot_I_PJN <- ggplot() + 
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, I), color = "gray0") +
  #geom_point(data = Ibyhr_Sled2_sub, aes(Group.1, x), color = "gray60") +
  #geom_point(data = Ibyhr_Sled1, aes(Group.1, x), color = "gray60") +
  #geom_point(data = Ibyhr_Sled3, aes(Group.1, x), color = "gray60") +
  theme_bw() +
  ylim(0, 4500000) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2017-2018)", y = bquote('μE m'^"-2"*' h'^"-1")) +
  ggtitle("Point Judith Pond N")
plot_I_PJS <- ggplot() + 
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, I), color = "gray0") +
  #geom_point(data = Ibyhr_Dredge, aes(Group.1, x), color = "gray60") +
  theme_bw() +
  ylim(0, 4500000) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2017-2018)", y = bquote('μE m'^"-2"*' h'^"-1")) +
  ggtitle("Point Judith Pond S")

plot_I_NBN_Y2 <- ggplot() + 
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, I), color = "gray0") +
  #geom_point(data = Ibyhr_WickfordY2L1_subset, aes(Group.1, x), color = "gray60") +
  ylim(0, 4500000) +
  xlim(as.POSIXct(c("2018-11-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2018-2019)", y = bquote('μE m'^"-2"*' h'^"-1")) +
  ggtitle("Narragansett Bay N")
plot_I_NBS_Y2 <- ggplot() + 
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, I), color = "gray0") +
  #logger flooded
  ylim(0, 4500000) +
  xlim(as.POSIXct(c("2018-11-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2018-2019)", y = bquote('μE m'^"-2"*' h'^"-1")) +
  ggtitle("Narragansett Bay S")
plot_I_PJN_Y2 <- ggplot() + 
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, I), color = "gray0") +
  #geom_point(data = Ibyhr_SledY2L1, aes(Group.1, x), color = "gray60") +
  ylim(0, 4500000) +
  xlim(as.POSIXct(c("2018-11-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2018-2019)", y = bquote('μE m'^"-2"*' h'^"-1")) +
  ggtitle("Point Judith Pond N")
plot_I_PJS_Y2 <- ggplot() + 
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, I), color = "gray0") +
  #geom_point(data = Ibyhr_DredgeY2L1, aes(Group.1, x), color = "gray60") +
  ylim(0, 4500000) +
  xlim(as.POSIXct(c("2018-11-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2018-2019)", y = bquote('μE m'^"-2"*' h'^"-1")) +
  ggtitle("Point Judith Pond S")
grid.arrange(plot_I_NBN, plot_I_NBS, plot_I_PJN, plot_I_PJS, plot_I_NBN_Y2, plot_I_NBS_Y2, plot_I_PJN_Y2, plot_I_PJS_Y2, ncol=4)

plot_J_I_NB <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, J_I, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, J_I, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = bquote('Spec. relaxation (mol NADPH mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("A)")
plot_J_I_PJ <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_I, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_I, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = bquote('Spec. relaxation (mol NADPH mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("B)")
plot_J_I_NB_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, J_I, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, J_I, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = bquote('Spec. relaxation (mol NADPH mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("C)")
plot_J_I_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, J_I, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, J_I, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = bquote('Spec. relaxation (mol NADPH mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("D)")

grid.arrange(plot_J_I_NB, plot_J_I_PJ, plot_J_I_NB_Y2, plot_J_I_PJ_Y2, ncol=2)

plot_J_EC_A_NB <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, J_EC_A, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, J_EC_A, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = bquote('C assimilation (mol C mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("A)")
plot_J_EC_A_PJ <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_EC_A, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_EC_A, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = bquote('C assimilation (mol C mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("B)")
plot_J_EC_A_NB_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, J_EC_A, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, J_EC_A, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = bquote('C assimilation (mol C mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("C)")
plot_J_EC_A_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, J_EC_A, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, J_EC_A, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = bquote('C assimilation (mol C mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("D)")
grid.arrange(plot_J_EC_A_NB, plot_J_EC_A_PJ, plot_J_EC_A_NB_Y2, plot_J_EC_A_PJ_Y2, ncol=2)

plot_J_EN_A_NB <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, J_EN_A, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, J_EN_A, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 3.1e-04) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol N mol Mv-1 h-1") +
  ggtitle("Narragansett Bay N and S")
plot_J_EN_A_PJ <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_EN_A, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_EN_A, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 3.1e-04) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol N mol Mv-1 h-1") +
  ggtitle("Point Judith Pond N and S")
plot_J_EN_A_NB_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, J_EN_A, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, J_EN_A, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 3.1e-04) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol N mol Mv-1 h-1") +
  ggtitle("Narragansett Bay N and S")
plot_J_EN_A_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, J_EN_A, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, J_EN_A, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 3.1e-04) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol N mol Mv-1 h-1") +
  ggtitle("Point Judith Pond N and S")
grid.arrange(plot_J_EN_A_NB, plot_J_EN_A_PJ, plot_J_EN_A_NB_Y2, plot_J_EN_A_PJ_Y2, ncol=2)

plot_J_EN_C_NB <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, J_EN_C, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, J_EN_C, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 5.1e-04) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol EN mol Mv-1 h-1") +
  ggtitle("Narragansett Bay N and S")
plot_J_EN_C_PJ <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_EN_C, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_EN_C, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 5.1e-04) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol EN mol Mv-1 h-1") +
  ggtitle("Point Judith Pond N and S")
plot_J_EN_C_NB_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, J_EN_C, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, J_EN_C, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 5.1e-04) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol EN mol Mv-1 h-1") +
  ggtitle("Narragansett Bay N and S")
plot_J_EN_C_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, J_EN_C, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, J_EN_C, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 5.1e-04) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol EN mol Mv-1 h-1") +
  ggtitle("Point Judith Pond N and S")
grid.arrange(plot_J_EN_C_NB, plot_J_EN_C_PJ, plot_J_EN_C_NB_Y2, plot_J_EN_C_PJ_Y2, ncol=2)

plot_J_EC_C_NB <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, J_EC_C, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, J_EC_C, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 0.027) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol EC mol Mv-1 h-1") +
  ggtitle("Narragansett Bay N and S")
plot_J_EC_C_PJ <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_EC_C, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_EC_C, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 0.027) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol EC mol Mv-1 h-1") +
  ggtitle("Point Judith Pond N and S")
plot_J_EC_C_NB_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, J_EC_C, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, J_EC_C, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 0.027) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol EC mol Mv-1 h-1") +
  ggtitle("Narragansett Bay N and S")
plot_J_EC_C_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, J_EC_C, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, J_EC_C, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 0.027) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol EC mol Mv-1 h-1") +
  ggtitle("Point Judith Pond N and S")
grid.arrange(plot_J_EC_C_NB, plot_J_EC_C_PJ, plot_J_EC_C_NB_Y2, plot_J_EC_C_PJ_Y2, ncol=2)

plot_J_EC_G_NB <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, J_EC_G, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, J_EC_G, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 0.027) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol EC mol Mv-1 h-1") +
  ggtitle("Narragansett Bay N and S")
plot_J_EC_G_PJ <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_EC_G, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_EC_G, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 0.027) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol EC mol Mv-1 h-1") +
  ggtitle("Point Judith Pond N and S")
plot_J_EC_G_NB_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, J_EC_G, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, J_EC_G, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 0.027) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol EC mol Mv-1 h-1") +
  ggtitle("Narragansett Bay N and S")
plot_J_EC_G_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, J_EC_G, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, J_EC_G, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 0.027) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol EC mol Mv-1 h-1") +
  ggtitle("Point Judith Pond N and S")
grid.arrange(plot_J_EC_G_NB, plot_J_EC_G_PJ, plot_J_EC_G_NB_Y2, plot_J_EC_G_PJ_Y2, ncol=2)

plot_J_EN_G_NB <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, J_EN_G, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, J_EN_G, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 5.1e-04) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol EN mol Mv-1 h-1") +
  ggtitle("Narragansett Bay N and S")
plot_J_EN_G_PJ <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_EN_G, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_EN_G, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 5.1e-04) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol EN mol Mv-1 h-1") +
  ggtitle("Point Judith Pond N and S")
plot_J_EN_G_NB_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, J_EN_G, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, J_EN_G, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 5.1e-04) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol EN mol Mv-1 h-1") +
  ggtitle("Narragansett Bay N and S")
plot_J_EN_G_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, J_EN_G, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, J_EN_G, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 5.1e-04) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol EN mol Mv-1 h-1") +
  ggtitle("Point Judith Pond N and S")
grid.arrange(plot_J_EN_G_NB, plot_J_EN_G_PJ, plot_J_EN_G_NB_Y2, plot_J_EN_G_PJ_Y2, ncol=2)

plot_r_NB <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, r, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, r, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 0.0047) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "r (h-1)") +
  ggtitle("Narragansett Bay N and S")
plot_r_PJ <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, r, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, r, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 0.0047) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "r (h-1)") +
  ggtitle("Point Judith Pond N and S")
plot_r_NB_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, r, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, r, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 0.0047) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "r (h-1)") +
  ggtitle("Narragansett Bay N and S")
plot_r_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, r, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, r, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 0.0047) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "r (h-1)") +
  ggtitle("Point Judith Pond N and S")
grid.arrange(plot_r_NB, plot_r_PJ, plot_r_NB_Y2, plot_r_PJ_Y2, ncol=2)

plot_J_EC_R_NB <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, J_EC_R, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, J_EC_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 0.035) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = bquote('Rejected C (mol EC mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("A)")
plot_J_EC_R_PJ <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_EC_R, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_EC_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 0.035) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = bquote('Rejected C (mol EC mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("B)")
plot_J_EC_R_NB_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, J_EC_R, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, J_EC_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 0.035) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = bquote('Rejected C (mol EC mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("C)")
plot_J_EC_R_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, J_EC_R, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, J_EC_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 0.035) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = bquote('Rejected C (mol EC mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("D)")
grid.arrange(plot_J_EC_R_NB, plot_J_EC_R_PJ, plot_J_EC_R_NB_Y2, plot_J_EC_R_PJ_Y2, ncol=2)

plot_J_EN_R_NB <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, J_EN_R, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, J_EN_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 3.2e-04) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = bquote('Rejected N (mol EN mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("A)")
plot_J_EN_R_PJ <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_EN_R, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_EN_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 3.2e-04) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = bquote('Rejected N (mol EN mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("B)")
plot_J_EN_R_NB_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, J_EN_R, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, J_EN_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 3.2e-04) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = bquote('Rejected N (mol EN mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("C)")
plot_J_EN_R_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, J_EN_R, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, J_EN_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 3.2e-04) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = bquote('Rejected N (mol EN mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("D)")
grid.arrange(plot_J_EN_R_NB, plot_J_EN_R_PJ, plot_J_EN_R_NB_Y2, plot_J_EN_R_PJ_Y2, ncol=2)

plot_m_EC_NB <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, m_EC, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, m_EC, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 0.45) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol C mol Mv-1") +
  ggtitle("Narragansett Bay N and S")
plot_m_EC_PJ <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, m_EC, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, m_EC, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 0.45) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol C mol Mv-1") +
  ggtitle("Point Judith Pond N and S")
plot_m_EC_NB_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, m_EC, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, m_EC, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 0.45) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol C mol Mv-1") +
  ggtitle("Narragansett Bay N and S")
plot_m_EC_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, m_EC, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, m_EC, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 0.45) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol C mol Mv-1") +
  ggtitle("Point Judith Pond N and S")
grid.arrange(plot_m_EC_NB, plot_m_EC_PJ, plot_m_EC_NB_Y2, plot_m_EC_PJ_Y2, ncol=2)

plot_m_EN_NB <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, m_EN, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, m_EN, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 0.06) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol N mol Mv-1") +
  ggtitle("Narragansett Bay N and S")
plot_m_EN_PJ <- ggplot() +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, m_EN, color = source)) +
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, m_EN, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0, 0.06) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = "mol N mol Mv-1") +
  ggtitle("Point Judith Pond N and S")
plot_m_EN_NB_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, m_EN, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, m_EN, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 0.06) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol N mol Mv-1") +
  ggtitle("Narragansett Bay N and S")
plot_m_EN_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, m_EN, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, m_EN, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  ylim(0, 0.06) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = "mol N mol Mv-1") +
  ggtitle("Point Judith Pond N and S")
grid.arrange(plot_m_EN_NB, plot_m_EN_PJ, plot_m_EN_NB_Y2, plot_m_EN_PJ_Y2, ncol=2)

##### Kelp Field Data for Comparison, filtered####
  ####Year 1####
KelpY1 <- read.csv("Year1kelpdata.csv", header = TRUE)
names(KelpY1)[2] <- "Site"
KelpY1 <- filter(KelpY1, Site != "Fox Island")
KelpY1$Date <- mdy(KelpY1$SamplingDate)
KelpY1$SiteLine <- paste(KelpY1$Site, KelpY1$Line)
KelpY1 <- filter(KelpY1, SiteLine != "Narragansett Bay N 2")

#Plot in format to compare field data with model output

plot_data <- ggplot() + 
  geom_point(data = KelpY1, aes(Date, Length, color = SiteLine), position = "jitter") +
  geom_smooth(data = KelpY1, aes(Date, Length, lty = SiteLine, color = SiteLine), se = FALSE) +
  theme_bw() +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed","solid", "dashed", "solid", "dashed")) +
  theme(legend.key.size =  unit(0.4, "in")) +
  labs(x= "Date", y = "Blade length (cm)") +
  ggtitle("A") +
  theme(legend.title = element_blank())

plot_model <- ggplot() + 
  geom_point(data = KelpY1, aes(Date, Length, color = SiteLine), position = "jitter") +
  geom_smooth(data = sol_all, aes(Date, L_allometric, lty = source, color = source), se = FALSE) +
  #scale_color_manual(values = c("gray0", "gray0", "gray60", "gray60", "gray30", "gray30", "gray77")) +
  #scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed","solid", "dashed", "solid")) +
  theme_bw() +
  labs(x= "Date", y = "Blade length (cm)") +
  ggtitle("B")
  #theme(legend.position="none") #CUTTING LEGEND BY CROPPING

grid.arrange(plot_data, plot_model, ncol=2)



#Sled Sugar Kelp Growth from lines planted on 11/1/17 and 11/29/17
#Dredge Sugar Kelp Growth from lines planted on 11/1/17 and 11/29/17
#Wickford Sugar Kelp Growth from line planted on 12/4/17
#Rome Point Sugar Kelp Growth from lines planted on 11/1/17 and 12/6/17

  ####Year 2####
KelpY2 <- read.csv("Year2kelpdata.csv", header = TRUE)
names(KelpY2)[2] <- "Site"
KelpY2 <- filter(KelpY2, Site != "Fox Island")
KelpY2$Date <- mdy(KelpY2$SamplingDate)
KelpY2$SiteLine <- paste(KelpY2$Site, KelpY2$Line)
KelpY2 <- filter(KelpY2, SiteLine != "Narragansett Bay N 2")

#Plot in format to compare field data with model output 
plot_dataY2 <- ggplot() + 
  geom_point(data = KelpY2, aes(Date, Length, color = SiteLine), position = "jitter") +
  geom_smooth(data = KelpY2, aes(Date, Length, lty = SiteLine, color = SiteLine), se = FALSE) +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed","solid", "dashed", "solid", "dashed")) +
  theme_bw() +
  theme(legend.key.size =  unit(0.4, "in")) +
  theme(legend.title = element_blank()) +
  labs(x= "Date", y = "Blade length (cm)") +
  ggtitle("A")

plot_modelY2 <- ggplot() + 
  geom_smooth(data = sol_all_Y2, aes(Date, L_allometric, lty = source, color = source)) +
  scale_color_manual(values = c("gray0", "gray0", "gray60", "gray60", "gray30", "gray30", "gray77")) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed","solid", "dashed", "solid")) +
  theme_bw() +
  labs(x= "Date", y = "Blade length (cm)") +
  ggtitle("B")

grid.arrange(plot_dataY2, plot_modelY2, ncol=2)

#Sled Sugar Kelp Growth from lines planted on 12/12/18 and 2/6/19
#Dredge Sugar Kelp Growth from lines planted on 12/12/18 and 2/6/19
#Wickford Sugar Kelp Growth from line planted on 12/19/18
#Rome Point Sugar Kelp Growth from lines planted on 12/20/18 and 2/21/19
  ####y1+y2 combined for publication #####
#model on field data
KelpY1 <- read.csv("Year1kelpdata.csv", header = TRUE)
names(KelpY1)[2] <- "Site"
KelpY1 <- filter(KelpY1, Site != "Fox Island")
KelpY1$Date <- mdy(KelpY1$SamplingDate)
KelpY1$SiteLine <- paste(KelpY1$Site, KelpY1$Line)
KelpY1 <- filter(KelpY1, SiteLine != "Narragansett Bay N 2")
KelpY2 <- read.csv("Year2kelpdata.csv", header = TRUE)
names(KelpY2)[2] <- "Site"
KelpY2 <- filter(KelpY2, Site != "Fox Island")
KelpY2$Date <- mdy(KelpY2$SamplingDate)
KelpY2$SiteLine <- paste(KelpY2$Site, KelpY2$Line)
KelpY2 <- filter(KelpY2, SiteLine != "Narragansett Bay N 2")


NBN1_meandat <- KelpY1[KelpY1$SiteLine == "Narragansett Bay N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length))
NBN1_meandat$Date <- as.POSIXct(NBN1_meandat$Date)
NBN1_meandat_sub <- NBN1_meandat[2:6,]
erNBN1 <- merge(NBN1_meandat_sub, sol_all[sol_all$source == "Narragansett Bay N 1",], all.x = TRUE)
NBN1_rmse <- rmse(erNBN1$mean_length, erNBN1$L_allometric)
NBN1_rmse <- round(NBN1_rmse, 2) #THIS RIGHT??
  
NBN1_meandat$Date <- as.POSIXct(NBN1_meandat$Date)
NBN1 <- ggplot() + 
  geom_point(data = NBN1_meandat, aes(Date, mean_length)) +
  geom_errorbar(NBN1_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", NBN1_rmse)) +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  scale_color_manual(values = c("gray0")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2017-2018)", y = "Blade length (cm)") +
  ggtitle("Narragansett Bay N") +
  theme(legend.position="none")

NBS1_meandat <- KelpY1[KelpY1$SiteLine == "Narragansett Bay S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
NBS1_meandat$Date <- as.POSIXct(NBS1_meandat$Date)
NBS1_meandat_sub <- NBS1_meandat[2:6,]
erNBS1 <- merge(NBS1_meandat_sub, sol_all[sol_all$source == "Narragansett Bay S 1",], all.x = TRUE)
NBS1_rmse <- rmse(erNBS1$mean_length, erNBS1$L_allometric)

NBS2_meandat <- KelpY1[KelpY1$SiteLine == "Narragansett Bay S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
NBS2_meandat$Date <- as.POSIXct(NBS2_meandat$Date)
NBS2_meandat_sub <- NBS2_meandat[2:6,]
erNBS2 <- merge(NBS2_meandat_sub, sol_all[sol_all$source == "Narragansett Bay S 2",], all.x = TRUE)
NBS2_rmse <- rmse(erNBS2$mean_length, erNBS2$L_allometric)

NBS1_2 <- ggplot() + 
  geom_point(data = NBS1_meandat, aes(Date, mean_length)) +
  geom_errorbar(NBS1_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", NBS1_rmse)) +
  geom_point(data = NBS2_meandat, aes(Date, mean_length), color = "gray50") +
  geom_errorbar(NBS2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all[sol_all$source == "Narragansett Bay S 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", NBS2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2017-2018)", y = "Blade length (cm)") +
  ggtitle("Narragansett Bay S") +
  theme(legend.position="none")

PJS1_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJS1_meandat$Date <- as.POSIXct(PJS1_meandat$Date)
PJS1_meandat_sub <- PJS1_meandat[2:6,]
erPJS1 <- merge(PJS1_meandat_sub, sol_all[sol_all$source == "Point Judith Pond S 1",], all.x = TRUE)
PJS1_rmse <- rmse(erPJS1$mean_length, erPJS1$L_allometric)

PJS2_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJS2_meandat$Date <- as.POSIXct(PJS2_meandat$Date)
PJS2_meandat_sub <- PJS2_meandat[2:6,]
erPJS2 <- merge(PJS2_meandat_sub, sol_all[sol_all$source == "Point Judith Pond S 2",], all.x = TRUE)
PJS2_rmse <- rmse(erPJS2$mean_length, erPJS2$L_allometric)

PJS1_2 <- ggplot() + 
  geom_point(data = PJS1_meandat, aes(Date, mean_length)) +
  geom_errorbar(PJS1_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-01 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", PJS1_rmse)) +
  geom_point(data = PJS2_meandat, aes(Date, mean_length), color = "gray50") +
  geom_errorbar(PJS2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all[sol_all$source == "Point Judith Pond S 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-01 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", PJS2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2017-2018)", y = "Blade length (cm)") +
  ggtitle("Point Judith Pond S") +
  theme(legend.position="none")

PJN1_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJN1_meandat$Date <- as.POSIXct(PJN1_meandat$Date)
PJN1_meandat_sub <- PJN1_meandat[2:6,]
erPJN1 <- merge(PJN1_meandat_sub, sol_all[sol_all$source == "Point Judith Pond N 1",], all.x = TRUE)
PJN1_rmse <- rmse(erPJN1$mean_length, erPJN1$L_allometric)

PJN2_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond N 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJN2_meandat$Date <- as.POSIXct(PJN2_meandat$Date)
PJN2_meandat_sub <- PJN2_meandat[2:6,]
erPJN2 <- merge(PJN2_meandat_sub, sol_all[sol_all$source == "Point Judith Pond N 2",], all.x = TRUE)
PJN2_rmse <- rmse(erPJN2$mean_length, erPJN2$L_allometric)


PJN1_2 <- ggplot() + 
  geom_point(data = PJN1_meandat, aes(Date, mean_length)) +
  geom_errorbar(PJN1_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", PJN1_rmse)) +
  geom_point(data = PJN2_meandat, aes(Date, mean_length), color ="gray50") +
  geom_errorbar(PJN2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all[sol_all$source == "Point Judith Pond N 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", PJN2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2017-2018)", y = "Blade length (cm)") +
  ggtitle("Point Judith Pond N") +
  theme(legend.position="none")

NBN1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Narragansett Bay N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
NBN1_Y2_meandat$Date <- as.POSIXct(NBN1_Y2_meandat$Date)
NBN1_Y2_meandat_sub <- NBN1_Y2_meandat[2:5,]
erNBN1_Y2 <- merge(NBN1_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], all.x = TRUE)
NBN1_Y2_rmse <- rmse(erNBN1_Y2$mean_length, erNBN1_Y2$L_allometric)

NBN1_Y2 <- ggplot() + 
  geom_point(data = NBN1_Y2_meandat, aes(Date, mean_length)) +
  geom_errorbar(NBN1_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", NBN1_Y2_rmse)) +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2018-12-01 12:00:00", "2019-05-30 12:00:00"))) +
  scale_color_manual(values = c("gray0")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2018-2019)", y = "Blade length (cm)") +
  ggtitle("Narragansett Bay N") +
  theme(legend.position="none")

NBS1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Narragansett Bay S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
NBS1_Y2_meandat$Date <- as.POSIXct(NBS1_Y2_meandat$Date)
NBS1_Y2_meandat_sub <- NBS1_Y2_meandat[2:4,]
erNBS1_Y2 <- merge(NBS1_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], all.x = TRUE)
NBS1_Y2_rmse <- rmse(erNBS1_Y2$mean_length, erNBS1_Y2$L_allometric)

NBS2_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Narragansett Bay S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
NBS2_Y2_meandat$Date <- as.POSIXct(NBS2_Y2_meandat$Date)
NBS2_Y2_meandat_sub <- NBS2_Y2_meandat[2:3,]
erNBS2_Y2 <- merge(NBS2_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 2",], all.x = TRUE)
NBS2_Y2_rmse <- rmse(erNBS2_Y2$mean_length, erNBS2_Y2$L_allometric)

NBS1_2_Y2 <- ggplot() + 
  geom_point(data = NBS1_Y2_meandat, aes(Date, mean_length)) +
  geom_errorbar(NBS1_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", NBS1_Y2_rmse)) +
  geom_point(data = NBS2_Y2_meandat, aes(Date, mean_length), color = "gray50") +
  geom_errorbar(NBS2_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", NBS2_Y2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2018-12-01 12:00:00", "2019-05-30 12:00:00"))) +
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2018-2019)", y = "Blade length (cm)") +
  ggtitle("Narragansett Bay S") +
  theme(legend.position="none")

PJS1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJS1_Y2_meandat$Date <- as.POSIXct(PJS1_Y2_meandat$Date)
PJS1_Y2_meandat_sub <- PJS1_Y2_meandat[2:6,]
erPJS1_Y2 <- merge(PJS1_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], all.x = TRUE)
PJS1_Y2_rmse <- rmse(erPJS1_Y2$mean_length, erPJS1_Y2$L_allometric)

PJS2_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJS2_Y2_meandat$Date <- as.POSIXct(PJS2_Y2_meandat$Date)
PJS2_Y2_meandat_sub <- PJS2_Y2_meandat[2:4,]
erPJS2_Y2 <- merge(PJS2_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 2",], all.x = TRUE)
PJS2_Y2_rmse <- rmse(erPJS2_Y2$mean_length, erPJS2_Y2$L_allometric)

PJS1_2_Y2 <- ggplot() + 
  geom_point(data = PJS1_Y2_meandat, aes(Date, mean_length)) +
  geom_errorbar(PJS1_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", PJS1_Y2_rmse)) +
  geom_point(data = PJS2_Y2_meandat, aes(Date, mean_length), color = "gray50") +
  geom_errorbar(PJS2_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", PJS2_Y2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2018-12-01 12:00:00", "2019-05-30 12:00:00"))) + 
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2018-2019)", y = "Blade length (cm)") +
  ggtitle("Point Judith Pond S") +
  theme(legend.position="none")

PJN1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJN1_Y2_meandat$Date <- as.POSIXct(PJN1_Y2_meandat$Date)
PJN1_Y2_meandat_sub <- PJN1_Y2_meandat[2:6,]
erPJN1_Y2 <- merge(PJN1_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], all.x = TRUE)
PJN1_Y2_rmse <- rmse(erPJN1_Y2$mean_length, erPJN1_Y2$L_allometric)

PJN2_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond N 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJN2_Y2_meandat$Date <- as.POSIXct(PJN2_Y2_meandat$Date)
PJN2_Y2_meandat_sub <- PJN2_Y2_meandat[2:4,]
erPJN2_Y2 <- merge(PJN2_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 2",], all.x = TRUE)
PJN2_Y2_rmse <- rmse(erPJN2_Y2$mean_length, erPJN2_Y2$L_allometric)

PJN1_2_Y2 <- ggplot() + 
  geom_point(data = PJN1_Y2_meandat, aes(Date, mean_length)) +
  geom_errorbar(PJN1_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", PJN1_Y2_rmse)) +
  geom_point(data = PJN2_Y2_meandat, aes(Date, mean_length), color = "gray50") +
  geom_errorbar(PJN2_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", PJN2_Y2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2018-12-01 12:00:00", "2019-05-30 12:00:00"))) +
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2018-2019)", y = "Blade length (cm)") +
  ggtitle("Point Judith Pond N") +
  theme(legend.position="none")

grid.arrange(NBN1, NBS1_2, PJN1_2, PJS1_2, NBN1_Y2, NBS1_2_Y2, PJN1_2_Y2, PJS1_2_Y2, ncol=4)

#aiming to create a figure like Romain's figure 6
er_all_y1 <- rbind(erNBN1, erNBS1, erNBS2, erPJS1, erPJS2, erPJN1, erPJN2)
er_all_y2 <- rbind(erNBN1_Y2, erNBS1_Y2, erNBS2_Y2, erPJN1_Y2, erPJN2_Y2, erPJS1_Y2, erPJS2_Y2)

y1ObsSim <- ggplot() +
  geom_point(data = er_all_y1, aes(mean_length, L_allometric, color = source)) +
  geom_errorbarh(er_all_y1, mapping = aes(xmin = mean_length-sd_length, xmax = mean_length+sd_length, y = L_allometric, color = source)) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_grey() +
  xlim(0, 250) +
  ylim(0, 250) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(legend.title = element_blank()) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Observed length (cm)", y = "Simulated length (cm)")

y2ObsSim <- ggplot() +
  geom_point(data = er_all_y2, aes(mean_length, L_allometric, color = source)) +
  geom_errorbarh(er_all_y2, mapping = aes(xmin = mean_length-sd_length, xmax = mean_length+sd_length, y = L_allometric, color = source)) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_grey() +
  xlim(0, 250) +
  ylim(0, 250) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(legend.title = element_blank()) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Observed length (cm)", y = "Simulated length (cm)")

grid.arrange(y1ObsSim, y2ObsSim, ncol=2)

#Austin wants a version with just end points:
er_ends_y1 <- rbind(erNBN1[5,], erNBS1[5,], erNBS2[5,], erPJS1[5,], erPJS2[5,], erPJN1[5,], erPJN2[5,]) 
er_ends_y2 <- rbind(erNBN1_Y2[4,], erNBS1_Y2[3,], erNBS2_Y2[2,], erPJN1_Y2[5,], erPJN2_Y2[3,], erPJS1_Y2[5,], erPJS2_Y2[3,])

y1ObsSim_ends <- ggplot() +
  geom_point(data = er_ends_y1, aes(mean_length, L_allometric, color = source), size = 2) +
  geom_errorbarh(er_ends_y1, mapping = aes(xmin = mean_length-sd_length, xmax = mean_length+sd_length, y = L_allometric, color = source)) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  xlim(0, 250) +
  ylim(0, 250) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(legend.title = element_blank()) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  ggtitle("A") +
  labs(x= "Observed length (cm)", y = "Simulated length (cm)")

y2ObsSim_ends <- ggplot() +
  geom_point(data = er_ends_y2, aes(mean_length, L_allometric, color = source), size = 2) +
  geom_errorbarh(er_ends_y2, mapping = aes(xmin = mean_length-sd_length, xmax = mean_length+sd_length, y = L_allometric, color = source)) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  xlim(0, 250) +
  ylim(0, 250) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(legend.title = element_blank()) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  ggtitle("B") +
  labs(x= "Observed length (cm)", y = "Simulated length (cm)")

grid.arrange(y1ObsSim_ends, y2ObsSim_ends, ncol=2)

########
##### Literature data for comparison/Calibration ####
  ######## Nitrate uptake ####
#Espinoza and Chapman (1983) and Ahn et al. (1998)
EC1983_9C_Nuptake_Fundy <- read.csv("EspinozaChapman1983_Nuptake_9C_Fundy.csv", header = TRUE)
EC1983_9C_Nuptake_StM <- read.csv("EspinozaChapman1983_Nuptake_9C_StMargaretsBay.csv", header = TRUE)
EC1983_18C_Nuptake_Fundy <- read.csv("EspinozaChapman1983_Nuptake_18C_Fundy.csv", header = TRUE)
EC1983_18C_Nuptake_StM <- read.csv("EspinozaChapman1983_Nuptake_18C_StMargaretsBay.csv", header = TRUE)

#conversions
EC1983_9C_Nuptake_StM$ResidualNitrateConcentration <- EC1983_9C_Nuptake_StM$ResidualNitrateConcentration/1000000 #microM to M
EC1983_9C_Nuptake_StM$NuptakeRate <- EC1983_9C_Nuptake_StM$NuptakeRate/1000000/w_EN #convert micro g N gDW–1 h–1 mol N gDW–1 h–1
EC1983_18C_Nuptake_StM$ResidualNitrateConcentration <- EC1983_18C_Nuptake_StM$ResidualNitrateConcentration/1000000 #microM to M
EC1983_18C_Nuptake_StM$NuptakeRate <- EC1983_18C_Nuptake_StM$NuptakeRate/1000000/w_EN
  
ggplot() +
  geom_line(data = sol_EspinozaChapman1983_N_9, mapping = aes(x = N, y = J_EN_A, color = "Model of Espinoza and Chapman (1983) at 9°C")) +
  geom_line(data = sol_EspinozaChapman1983_N_18, mapping = aes(x = N, y = J_EN_A, color = "Model of Espinoza and Chapman (1983) at 18°C")) +
  #geom_point(data = EC1983_9C_Nuptake_Fundy, mapping = aes(x = ResidualNitrateConcentration, y = NuptakeRate, color="Espinoza and Chapman (1983), Bay of Fundy sample, 9°C")) +
  geom_point(data = EC1983_9C_Nuptake_StM, mapping = aes(x = ResidualNitrateConcentration, y = NuptakeRate, color="Espinoza and Chapman (1983), St. Margaret's Bay, 9°C"), size=3) +
  #geom_point(data = EC1983_18C_Nuptake_Fundy, mapping = aes(x = ResidualNitrateConcentration, y = NuptakeRate, color="Espinoza and Chapman (1983), Bay of Fundy sample, 18°C")) +
  geom_point(data = EC1983_18C_Nuptake_StM, mapping = aes(x = ResidualNitrateConcentration, y = NuptakeRate, color="Espinoza and Chapman (1983), St. Margaret's Bay, 18°C"), size=3) +
  xlim(0, 8e-05) +
  scale_color_manual(values = c("gray60", "gray0", "gray60", "gray0")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= bquote('Nitrate Concentration (mol' ~NO[3]^{"-"}~ 'g DW'^"-1"*' h'^"-1"*')'), y = bquote('Nitrate uptake (mol N g DW'^"-1"*' h'^"-1"*')'))

#Error calculations
EC1983_9C_Nuptake_StM$N <- EC1983_9C_Nuptake_StM$ResidualNitrateConcentration
EC1983_9C_Nuptake_StM$N <- round(EC1983_9C_Nuptake_StM$N, digits = 7)
er9 <- merge(EC1983_9C_Nuptake_StM, sol_EspinozaChapman1983_N_9, all.x = TRUE)
er9$J_EN_A[3] <- 3.321263e-06 #fixing mystery NA
er9$J_EN_A[6] <-  5.821903e-06 #fixing mystery NA
er9$J_EN_A[14] <-  9.696641e-06 #fixing mystery NA
rmse(er9$NuptakeRate, er9$J_EN_A) #2.7e-4 and 2.667e-6 ->1.42595e-06

EC1983_18C_Nuptake_StM$N <- EC1983_18C_Nuptake_StM$ResidualNitrateConcentration
EC1983_18C_Nuptake_StM$N <- round(EC1983_18C_Nuptake_StM$N, digits = 7)
er18 <- merge(EC1983_18C_Nuptake_StM, sol_EspinozaChapman1983_N_18, all.x = TRUE)
er18$J_EN_A[14] <- 9.116292e-06 #fixing mystery NA
rmse(er18$NuptakeRate, er18$J_EN_A) #2.7e-4 and 2.667e-6 -> 9.727939e-07

  ####### Photosynthesis related ####
#Johansson2002
Johansson2002 <- read.csv("Johansson2002.csv", header = TRUE)
#conversions
Johansson2002$Irradiance <- Johansson2002$Irradiance*3600 #micromol photons m-2 s-1 to micromol photons m-2 h-1
Johansson2002$O2production <- Johansson2002$O2production/1e+6*32/1000*3600 #micromol O2 kg DW-1 s-1 to g O2/g/h
Johansson2002$O2productionSHIFT <- Johansson2002$O2production + 0.001720976 #0.00172097569728

#color = "Model of Johansson and Snoeijs (2002)"
#color = "Johansson and Snoeijs (2002)"
ggplot(data = Johansson2002) +
  geom_line(data = sol_Johansson2002, mapping = aes(x = I, y = J_O)) +
  geom_point(mapping = aes(x = Irradiance, y = O2productionSHIFT), size = 3) +
  scale_color_grey() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= bquote('Irradiance (μE m'^"-2"*' d'^"-1"*')'), y = bquote('Oxygen production (g' ~O[2]~ 'g DW'^"-1"*' h'^"-1"*')'))

#error calculations
Johansson2002$I <- round(Johansson2002$Irradiance, digits = 0)
erPhoto <- merge(Johansson2002, sol_Johansson2002, all.x = TRUE)
rmse(erPhoto$O2productionSHIFT, erPhoto$J_O)

