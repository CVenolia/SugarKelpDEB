####################################################################################################################
#Run file of Sugar Kelp model from Venolia et al (in press)
#Site names here begin with names other than those used in the manuscript
#Sled = Pt Judith Pond N
#Dredge = Pt Judith Pond S
#Wickford = Narragansett Bay N
#Rome Point = Narragansett Bay S
#File created by Celeste Venolia in March 2018-December 2019

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
#Conversion coefficients, organics (n = matrix of chemical indices)
# "food N" "food C" Stucture "N reserves" "C reserves" products
#     X_N   X_C      V    E_N    E_C    P
n_O <- matrix(
    + c(0.00, 1.00, 1.00, 0.00, 1.00, 1.00,  #C/C, equals 1 by definition
      + 1.50, 0.50, 1.33, 0.00, 2.00, 1.80,  #H/C, these values show that we consider dry-mass
      + 1.50, 2.50, 1.00, 3.00, 1.00, 0.50,  #O/C
      + 1.00, 0.00, 0.04, 1.00, 0.00, 0.04), nrow=4, ncol=6, byrow = TRUE) #N/C
#V is the C-mol structure of alginate (Alginic acid: (C6H8O6)n)
#E_N is N03-
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
               JENAM = 2.5e-4, #1.2e-4, #mol N / molM_V / h
               #maximum surface-specific assimilation rate of N
               K_N = 2.667e-6, #molNO3 and NO2/L
               #max volume-specific carbon dioxide assimilation rate
               JCO2M = 0.0075, #molC/molM_V/h
               #half saturation constant of C uptake
               K_C = 4e-06, #mol DIC/L
               #maximum volume-specific carbon assimilation rate
               JECAM = 0.282, #molC/molM_V/h
               #Photosynthetic unit density
               rho_PSU = 0.9, #mol PSU/ mol Mv
               #binding probability of photons to a Light SU
               b_I = 0.55, #dimensionless
               #Specific photon arrival cross section
               alpha = 1, #m^2 mol PSU–1
               #dissociation rate
               k_I = 0.065, #molγ molPS–1 h–1
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
               JENM = 4.2e-07, #mol N/molM_V/h
               #specific maintenance costs requiring C before temp correction
               JECM = 1e-07, #mol C/molM_V/h
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
####### State Inititial conditions ############
#Initial conditions of state variables
#these values are not coming from any field data or literature information, estimated
state_Lo <- c(m_EC = 0.9, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
              m_EN = 0.03, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
              M_V = 4.503984e-05) #molM_V #initial mass of structure
#0.0039 = (29.89+0.9*62+0.03*30)*M_V
#0.4
#4.503984e-05
#8.914286e-05

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

#############
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#YEAR 1
#
#Setting up the forcing functions with field data for Sled line 1
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Sled1 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Sled line 2
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Sled2 <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Dredge line 1
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Dredge1 <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Dredge line 2
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Dredge2 <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Wickford line 1
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Wickford1 <- ode(y = state_Lo, t = times_Lo_Wickford1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Rome Point line 1
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_RomePt1 <- ode(y = state_Lo, t = times_Lo_RomePt1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Rome Point line 2
W <- 0.0039 #1.6545 #inital biomass for conversions (cannot put in initial conditions)
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
################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_RomePt2 <- ode(y = state_Lo, t = times_Lo_RomePt2, func = rates_Lo, parms = params_Lo)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#YEAR 2 Kelp data
#Setting up the forcing functions with field data for Sled line 1 (y2)
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Sled1_Y2 <- ode(y = state_Lo, t = times_Y2_Sled1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Sled line 2 (y2)
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Sled2_Y2 <- ode(y = state_Lo, t = times_Y2_Sled2, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Dredge line 1 (y2)
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Dredge1_Y2 <- ode(y = state_Lo, t = times_Y2_Dredge1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Dredge line 2 (y2)
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Dredge2_Y2 <- ode(y = state_Lo, t = times_Y2_Dredge2, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Wickford line 1 (y2)
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Wickford1_Y2 <- ode(y = state_Lo, t = times_Y2_Wickford1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Rome Point line 1 (y2)
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_RomePt1_Y2 <- ode(y = state_Lo, t = times_Y2_RomePt1, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for Rome Point line 2 (y2)
W <- 0.0039 #inital biomass for conversions (cannot put in initial conditions)
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
#####################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_RomePt2_Y2 <- ode(y = state_Lo, t = times_Y2_RomePt2, func = rates_Lo, parms = params_Lo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#END FIELD DATA, START LITERATURE DATA FOR CALIBRATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions for Espinoza and Chapman (1983) nitrogen uptake #9C
###### N forcing set-up##############
Nmax <- 73.1221719/1000000 #M

###### Temperature forcing set-Up #############
#9 C is 282.15 K
T_dat <- 9 #C (conversion in Nuptake function to K)
###################################
#Model run (the differential equation solver)
sol_EspinozaChapman1983_N_9 <- Nuptake(params_Lo, T_dat, Nmax, w_EN) #function from N_uptake_Calibration.R code
sol_EspinozaChapman1983_N_9 <- as.data.frame(sol_EspinozaChapman1983_N_9) #conversion to dataframe for later use
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions for Espinoza and Chapman (1983) nitrogen uptake #18C
###### N forcing set-up##############
Nmax <- 76.9543147/1000000 #M

###### Temperature forcing set-Up #############
#18 C is 282.15 K
T_dat <- 18 #C (conversion in Nuptake function to K)
#################################
#Model run (the differential equation solver)
sol_EspinozaChapman1983_N_18 <- Nuptake(params_Lo, T_dat, Nmax, w_EN) #function from N_uptake_Calibration.R code
sol_EspinozaChapman1983_N_18 <- as.data.frame(sol_EspinozaChapman1983_N_18) #conversion to dataframe for later use
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Photosynthesis model calibration
##### Temp forcing set-up ####
T_dat <- 14 #C (maintained for entire experiment
###### I forcing set-up #####
I_max <- 3233205 #micromol photons m-2 h-1
##############
sol_Johansson2002 <- Photosynthesis(params_Lo, state_Lo, w_V, w_EN, w_EC, w_O2, T_dat, I_max) #function from Photosynthesis_Calibration.R
sol_Johansson2002 <- as.data.frame(sol_Johansson2002) #conversion to dataframe for later use
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###### Convert DeSolve solutions into data frame for broader plotting use ####
  ##### Year 1 #####
#conversions to dataframes
sol_Sled1 <- as.data.frame(sol_Sled1)
sol_Sled2 <- as.data.frame(sol_Sled2)
sol_Dredge1 <- as.data.frame(sol_Dredge1)
sol_Dredge2 <- as.data.frame(sol_Dredge2)
sol_Wickford1 <- as.data.frame(sol_Wickford1)
sol_RomePt1 <- as.data.frame(sol_RomePt1)
sol_RomePt2 <- as.data.frame(sol_RomePt2)

#addition of a date variable
sol_Sled1$Date <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
sol_Sled2$Date <- seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
sol_Dredge1$Date <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
sol_Dredge2$Date <-seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
sol_Wickford1$Date <- seq(as_datetime("2017-12-4 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
sol_RomePt1$Date <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
sol_RomePt2$Date  <- seq(as_datetime("2017-12-6 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")

#conversion back to Celsius from Kelvin
sol_Sled1$Temp_C <- T_Sled1_Y1 - 273.15
sol_Sled2$Temp_C <- T_Sled2_Y1 - 273.15
sol_Dredge1$Temp_C <- T_Dredge1_Y1 - 273.15
sol_Dredge2$Temp_C <- T_Dredge2_Y1 - 273.15
sol_Wickford1$Temp_C <- T_Wickford1_Y1 - 273.15
sol_RomePt1$Temp_C <- T_RomePt1_Y1 - 273.15
sol_RomePt2$Temp_C <- T_RomePt2_Y1 - 273.15

#create source collumn to prepare for binding all these dataframes together
sol_Sled1$source <- "Point Judith Pond N 1"
sol_Sled2$source <- "Point Judith Pond N 2"
sol_Dredge1$source <- "Point Judith Pond S 1"
sol_Dredge2$source <- "Point Judith Pond S 2"
sol_Wickford1$source <- "Narragansett Bay N 1"
sol_RomePt1$source <- "Narragansett Bay S 1"
sol_RomePt2$source  <- "Narragansett Bay S 2"

#combine all Y1 field data into one dataframe
sol_all <- rbind(sol_Dredge1, sol_Dredge2, sol_RomePt1, sol_RomePt2, sol_Sled1, sol_Sled2, sol_Wickford1)

  ##### Year 2 #####
#conversions to dataframes
sol_Sled1_Y2 <- as.data.frame(sol_Sled1_Y2)
sol_Sled2_Y2 <- as.data.frame(sol_Sled2_Y2)
sol_Dredge1_Y2 <- as.data.frame(sol_Dredge1_Y2)
sol_Dredge2_Y2 <- as.data.frame(sol_Dredge2_Y2)
sol_Wickford1_Y2 <- as.data.frame(sol_Wickford1_Y2)
sol_RomePt1_Y2 <- as.data.frame(sol_RomePt1_Y2)
sol_RomePt2_Y2 <- as.data.frame(sol_RomePt2_Y2)

#addition of a date variable
sol_Sled1_Y2$Date <- seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour") 
sol_Sled2_Y2$Date <- seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
sol_Dredge1_Y2$Date <- seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
sol_Dredge2_Y2$Date <- seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
sol_Wickford1_Y2$Date <- seq(as_datetime("2018-12-19 12:00:00"), as_datetime("2019-05-23 12:00:00"), by="hour")
sol_RomePt1_Y2$Date <- seq(as_datetime("2018-12-20 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
sol_RomePt2_Y2$Date <- seq(as_datetime("2019-2-21 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")

#conversion back to Celsius from Kelvin
sol_Sled1_Y2$Temp_C <- T_Sled1_Y2 - 273.15
sol_Sled2_Y2$Temp_C <- T_Sled2_Y2 - 273.15
sol_Dredge1_Y2$Temp_C <- T_Dredge1_Y2 - 273.15
sol_Dredge2_Y2$Temp_C <- T_Dredge2_Y2 - 273.15
sol_Wickford1_Y2$Temp_C <- T_Wickford1_Y2 - 273.15
sol_RomePt1_Y2$Temp_C <- T_RomePt1_Y2 - 273.15
sol_RomePt2_Y2$Temp_C <- T_RomePt2_Y2 - 273.15

#create source collumn to prepare for binding all these dataframes together
sol_Sled1_Y2$source <- "Point Judith Pond N 1"
sol_Sled2_Y2$source <- "Point Judith Pond N 2"
sol_Dredge1_Y2$source <- "Point Judith Pond S 1"
sol_Dredge2_Y2$source <- "Point Judith Pond S 2"
sol_Wickford1_Y2$source <- "Narragansett Bay N 1"
sol_RomePt1_Y2$source <- "Narragansett Bay S 1"
sol_RomePt2_Y2$source  <- "Narragansett Bay S 2"

#combine all Y2 field data into one dataframe
sol_all_Y2 <- rbind(sol_Dredge1_Y2, sol_Dredge2_Y2, sol_RomePt1_Y2, sol_RomePt2_Y2, sol_Sled1_Y2, sol_Sled2_Y2, sol_Wickford1_Y2)

##### Model Plots #####
#Figure 6
#Temperature y1 plot
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
#Temperature y2 plot
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
#N forcing y1 plot
plot_N <- ggplot(data = sol_all, aes(Date, N, color = source)) + 
  geom_line(size = 2) +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  ylim(0, 1e-05) +
  labs(x= "Date (2017-2018)", y = bquote('mol' ~NO[3]^{"-"}~ 'and' ~NO[2]^{"-"}~ 'L'^"-1")) +
  ggtitle("C")
#N forcing y2 plot
plot_N_Y2 <- ggplot(data = sol_all_Y2, aes(Date, N, color = source)) + 
  geom_line(size = 2) +
  scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  ylim(0, 1e-05) +
  labs(x= "Date (2018-2019)", y = bquote('mol' ~NO[3]^{"-"}~ 'and'~NO[2]^{"-"}~ 'L'^"-1")) +
  ggtitle("D")
grid.arrange(plot_T_Y1, plot_T_Y2, plot_N, plot_N_Y2, ncol=2) #gridded plot

#Figure 3: combining all irradiance forcings
plot_I_NBN <- ggplot() + 
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, I), color = "gray0") +
  theme_bw() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  xlab("Date (2017-2018)") +
  ylab(bquote('E m'^"-2"*' h'^"-1")) + #"E m-2 d-1"
  ggtitle("Narragansett Bay N")

plot_I_NBS <- ggplot() + 
  geom_point(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, I), color = "gray0") +
  theme_bw() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2017-2018)", y = bquote('E m'^"-2"*' h'^"-1")) +
  ggtitle("Narragansett Bay S")

plot_I_PJN <- ggplot() + 
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, I), color = "gray0") +
  theme_bw() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2017-2018)", y = bquote('E m'^"-2"*' h'^"-1")) +
  ggtitle("Point Judith Pond N")

plot_I_PJS <- ggplot() + 
  geom_point(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, I), color = "gray0") +
  theme_bw() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2017-2018)", y = bquote('E m'^"-2"*' h'^"-1")) +
  ggtitle("Point Judith Pond S")

plot_I_NBN_Y2 <- ggplot() + 
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, I), color = "gray0") +
  xlim(as.POSIXct(c("2018-11-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2018-2019)", y = bquote('E m'^"-2"*' h'^"-1")) +
  ggtitle("Narragansett Bay N")

plot_I_NBS_Y2 <- ggplot() + 
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, I), color = "gray0") +
  xlim(as.POSIXct(c("2018-11-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2018-2019)", y = bquote('E m'^"-2"*' h'^"-1")) +
  ggtitle("Narragansett Bay S")

plot_I_PJN_Y2 <- ggplot() + 
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, I), color = "gray0") +
  xlim(as.POSIXct(c("2018-11-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2018-2019)", y = bquote('E m'^"-2"*' h'^"-1")) +
  ggtitle("Point Judith Pond N")

plot_I_PJS_Y2 <- ggplot() + 
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, I), color = "gray0") +
  xlim(as.POSIXct(c("2018-11-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2018-2019)", y = bquote('E m'^"-2"*' h'^"-1")) +
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

#Figure 7
plot_J_EC_R_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, J_EC_R, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, J_EC_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = bquote('Rejected C (mol C mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("A)")
plot_J_EN_R_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, J_EN_R, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, J_EN_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = bquote('Rejected N (mol N mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("B)")
grid.arrange(plot_J_EC_R_PJ_Y2, plot_J_EN_R_PJ_Y2, ncol=2)

#Figure 8
plot_J_I_PJ_Y2 <- ggplot() +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, J_I, color = source)) +
  geom_point(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, J_I, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2018-10-30 23:00:00", "2019-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2018-2019)", y = bquote('Spec. relaxation (mol E mol Mv'^"-1"*' h'^"-1"*')')) +
  ggtitle("A)")
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
  ggtitle("B)")
grid.arrange(plot_J_I_PJ_Y2, plot_J_EC_A_PJ_Y2, ncol=2)

##### Kelp Field Data Comparison plot (Figure 9) ####
#import field data
KelpY1 <- read.csv("Year1kelpdata.csv", header = TRUE, fileEncoding="UTF-8-BOM")
names(KelpY1)[2] <- "Site"
KelpY1 <- filter(KelpY1, Site != "Fox Island")
KelpY1$Date <- mdy(KelpY1$SamplingDate)
KelpY1$SiteLine <- paste(KelpY1$Site, KelpY1$Line)
KelpY1 <- filter(KelpY1, SiteLine != "Narragansett Bay N 2")
KelpY2 <- read.csv("Year2kelpdata.csv", header = TRUE, fileEncoding="UTF-8-BOM")
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
NBN1_rmse <- round(NBN1_rmse, 2)()
  
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

########
##### Literature data for comparison/Calibration ####
  ######## Nitrate uptake ####
#Espinoza and Chapman (1983) and Ahn et al. (1998)
EC1983_9C_Nuptake_StM <- read.csv("EspinozaChapman1983_Nuptake_9C_StMargaretsBay.csv", header = TRUE, fileEncoding="UTF-8-BOM")
EC1983_18C_Nuptake_StM <- read.csv("EspinozaChapman1983_Nuptake_18C_StMargaretsBay.csv", header = TRUE, fileEncoding="UTF-8-BOM")

#conversions 9C
EC1983_9C_Nuptake_StM$N <- EC1983_9C_Nuptake_StM$ResidualNitrateConcentration
EC1983_9C_Nuptake_StM$N <- round(EC1983_9C_Nuptake_StM$N, digits = 2)
EC1983_9C_Nuptake_StM$N <- EC1983_9C_Nuptake_StM$N/1000000 #microM to M
EC1983_9C_Nuptake_StM$NuptakeRate <- EC1983_9C_Nuptake_StM$NuptakeRate/1000000/w_EN #convert micro g N gDW–1 h–1 to mol N gDW–1 h–1
#conversions 18C
EC1983_18C_Nuptake_StM$N <- EC1983_18C_Nuptake_StM$ResidualNitrateConcentration
EC1983_18C_Nuptake_StM$N <- round(EC1983_18C_Nuptake_StM$N, digits = 2)
EC1983_18C_Nuptake_StM$N <- EC1983_18C_Nuptake_StM$N/1000000 #microM to M
EC1983_18C_Nuptake_StM$NuptakeRate <- EC1983_18C_Nuptake_StM$NuptakeRate/1000000/w_EN

#testing rounding
sol_EspinozaChapman1983_N_9$N <- round(sol_EspinozaChapman1983_N_9$N*1000000, digits = 3)/1000000
sol_EspinozaChapman1983_N_18$N <- round(sol_EspinozaChapman1983_N_18$N*1000000, digits = 3)/1000000

N_calibration <- ggplot() +
  geom_line(data = sol_EspinozaChapman1983_N_9, mapping = aes(x = N, y = J_EN_A, color = "Model of Espinoza and Chapman (1983) at 9°C")) +
  geom_line(data = sol_EspinozaChapman1983_N_18, mapping = aes(x = N, y = J_EN_A, color = "Model of Espinoza and Chapman (1983) at 18°C")) +
  geom_point(data = EC1983_9C_Nuptake_StM, mapping = aes(x = N, y = NuptakeRate, color="Espinoza and Chapman (1983), St. Margaret's Bay, 9°C"), size=3) +
  geom_point(data = EC1983_18C_Nuptake_StM, mapping = aes(x = N, y = NuptakeRate, color="Espinoza and Chapman (1983), St. Margaret's Bay, 18°C"), size=3) +
  xlim(0, 8e-05) +
  scale_color_manual(values = c("gray60", "gray0", "gray60", "gray0")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= bquote('Nitrate Concentration (mol' ~NO[3]^{"-"}~ 'L'^"-1"*')'), y = bquote('N uptake (mol' ~NO[3]^{"-"}~ 'and' ~NO[2]^{"-"}~ 'g DW'^"-1"*' h'^"-1"*')')) +
  ggtitle('A')

#Error calculations
er9 <- merge(EC1983_9C_Nuptake_StM, sol_EspinozaChapman1983_N_9, all.x = TRUE)
rmse(er9$NuptakeRate, er9$J_EN_A) #3.683799e-07

er18 <- merge(EC1983_18C_Nuptake_StM, sol_EspinozaChapman1983_N_18, all.x = TRUE)
rmse(er18$NuptakeRate, er18$J_EN_A) #2.606024e-07

  ######## Photosynthesis related ####
#Johansson2002
Johansson2002 <- read.csv("Johansson2002.csv", header = TRUE, fileEncoding="UTF-8-BOM")
#conversions
Johansson2002$Irradiance <- Johansson2002$Irradiance*3600*1e-6 #micromol photons m-2 s-1 to mol photons m-2 h-1
Johansson2002$O2production <- Johansson2002$O2production/1e+6*32/1000*3600 #micromol O2 kg DW-1 s-1 to g O2/g/h
Johansson2002$O2productionSHIFT <- Johansson2002$O2production + 0.001720976 #from net to gross

Photosynthesis_calibration <- ggplot(data = Johansson2002) +
  geom_line(data = sol_Johansson2002, mapping = aes(x = I, y = J_O)) +
  geom_point(mapping = aes(x = Irradiance, y = O2productionSHIFT), size = 3) +
  scale_color_grey() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= bquote('Irradiance (E m'^"-2"*' d'^"-1"*')'), y = bquote('Oxygen production (g' ~O[2]~ 'g DW'^"-1"*' h'^"-1"*')')) +
  ggtitle('B')

#error calculations
Johansson2002$I <- round(Johansson2002$Irradiance, digits = 6)
sol_Johansson2002$I <- round(sol_Johansson2002$I, digits = 6)
erPhoto <- merge(Johansson2002, sol_Johansson2002, all.x = TRUE)
rmse(erPhoto$O2productionSHIFT, erPhoto$J_O)

  ######## Combine calibration plot (Figure 5) #######
#Figure 5
grid.arrange(N_calibration, Photosynthesis_calibration, ncol=2)
