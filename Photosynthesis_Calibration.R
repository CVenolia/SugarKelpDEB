#####################################################################
#Photosythesis parameter calibration for sugar kelp DEB model
#Created 7/12/19 for Venolia et al (in press)
#####################################################################

#Function code to run only the photosynthesis related part of the sugar kelp DEB model
Photosynthesis <- function(params_Lo, state_Lo, w_V, w_EN, w_EC, w_O2, T_dat, I_max){
  #set-up relevant parameters in a usable form
  params_Lo <- as.list(params_Lo)
  T_A <- params_Lo$T_A
  T_0 <- params_Lo$T_0
  T_AL <- params_Lo$T_AL
  T_L <- params_Lo$T_L
  T_H <- params_Lo$T_H
  T_AH <- params_Lo$T_AH
  rho_PSU <- params_Lo$rho_PSU
  b_I <- params_Lo$b_I
  alpha <- params_Lo$alpha
  k_I <- params_Lo$k_I
  y_LO2 <- params_Lo$y_LO2
  state_Lo <- as.list(state_Lo)
  m_EC <- state_Lo$m_EC
  m_EN <- state_Lo$m_EN
  M_V <- state_Lo$M_V
  
  T_dat <- T_dat + 273.15 #if converstion to Kelvin is needed
  #temperature correction
  C_T <- exp((T_A/T_0)-(T_A/T_dat)) * (1+exp((T_AL/T_0)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_0))) * ((1+exp((T_AL/T_dat)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_dat)))^-1)
  #Dry weight Biomass calculation
  W <- (w_V + m_EN * w_EN + m_EC * w_EC) * M_V  # g
  
  #Irradiance steps
  I <- seq(0, I_max, 1)
  I <- I*1e-6 #micromol to mol
  #Specific relaxation rate
  J_I <- (rho_PSU * I * b_I * alpha)/(1 + (I * b_I * alpha)/(k_I * C_T)) #Unit: molγ molM_V–1 h–1
  #rate of oxygen production
  J_O <- J_I * (M_V/W) * w_O2 * y_LO2 # g O2/g/h O2
  
  return(cbind(I, J_O))
}
