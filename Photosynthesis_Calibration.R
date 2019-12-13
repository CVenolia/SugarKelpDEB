#Conversion of Romain's photosynthesis.m code
#Celeste Venolia on 7/12/19

Photosynthesis <- function(params_Lo, cond, w_V, w_EN, w_EC, w_O2, T_dat, I_max){
  params_Lo <- as.list(params_Lo)
  T_A <- params_Lo$T_A
  T_0 <- params_Lo$T_0
  T_AL <- params_Lo$T_AL
  T_L <- params_Lo$T_L
  T_H <- params_Lo$T_H
  T_AH <- params_Lo$T_AH
  rho_PSU <- params_Lo$rho_PSU #0.05
  #rho_PSU <- 0.05 #these commented out values are the "tight fit"
  b_I <- params_Lo$b_I #9e-3
  #b_I <- 2.8e-6
  k_I <- params_Lo$k_I #0.29
  #k_I <- 0.28
  condinit <- as.list(cond)
  m_EC <- condinit$m_EC
  m_EN <- condinit$m_EN
  M_V <- condinit$M_V
  
  T_dat <- T_dat + 273.15 #if converstion to Kelvin is needed
  C_T <- exp((T_A/T_0)-(T_A/T_dat)) * (1+exp((T_AL/T_0)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_0))) * ((1+exp((T_AL/T_dat)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_dat)))^-1) #C_T has no units
  B <- (w_V + m_EN * w_EN + m_EC * w_EC) * M_V  # g
  
  I <- seq(0, I_max, 1)
  J_I <- (rho_PSU * I * b_I)/(1 + I * b_I/(k_I * C_T)) #molNADPH molM_V–1 h–1
  J_O <- (J_I * M_V * w_O2) / (B * 4) # g O2/g/h O2
  
  #setting this up for comparison to mgO2/m^2/h units
  #J_I molNADPH molM_V–1 h–1
  #J_O <- J_I * M_V * w_O2/ (B * 4) # mg O2/g biomass/h O2 production (1 molO2 for 4 molC of e-/NADPH produced)
  #J_O <- J_I * w_O2 * 1000/ 4 #mg O2/molM_V/h
  #d_V <- 1 #the specific structural mass (d_V typically equals 1 g/cm^3)
  #J_O <- (J_I * w_O2 * 1000/ 4)*(d_V / w_V)*2/10000 #mgO2/m^2/h
  #J_O (mgO2/cm^2/h) = 2 * J_O (mgO2/cm^3/h)
  
  return(cbind(I, J_O))
}
