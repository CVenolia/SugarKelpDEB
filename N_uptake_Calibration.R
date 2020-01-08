#########################################################################
#Code for calibration of N uptake parameters for the sugar kelp DEB model
#for Venolia et al (in press) created on 7/12/19
#########################################################################

#function including just the parts of the sugar kelp DEB model related to N uptake
Nuptake <- function(params_Lo, T_dat, Nmax, w_EN){
  #set-up relevant parameters needed
  params_Lo <- as.list(c(params_Lo))
  T_A <- params_Lo$T_A
  T_0 <- params_Lo$T_0
  T_AL <- params_Lo$T_AL
  T_L <- params_Lo$T_L
  T_H <- params_Lo$T_H
  T_AH <- params_Lo$T_AH
  JENAM <- params_Lo$JENAM
  K_N <- params_Lo$K_N
  W <- 0.019 #These should be set based on whatever information you have on literature or field data being used
  M_V <- 0.0003418 #in Espinoza and Chapman (1983), blade length starts at around 3 cm, which allows for an estmiation of B
  #0.019 = (29.89+0.4*62+0.03*30)*M_V
  
  T_dat <- T_dat + 273.15 #if converstion to Kelvin is needed
  #temperature correction
  C_T <- exp((T_A/T_0)-(T_A/T_dat)) * (1+exp((T_AL/T_0)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_0))) * ((1+exp((T_AL/T_dat)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_dat)))^-1)
  N <- seq(0, Nmax, 1e-08) #array from 0 to Nmax stepped by an interval of 1e-08
  
  #this first conversion is taking this parameter from model units and putting it in the units of the calibration data
  JENAM <- JENAM * M_V/W #convert from mol N / molM_V / h to mol N / g DW / hour
  J_EN_AM <-  JENAM * C_T #temperature correct max assimilation rate of N03
  #Specific assimilation rate of N
  J_EN_A <- J_EN_AM*(N/(N+K_N))
  return(cbind(N, J_EN_A))  #this is what the function outputs
}

