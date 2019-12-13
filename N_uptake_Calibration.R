#Conversion of Romain's Nuptake.m code
#Celeste Venolia on 7/12/19

Nuptake <- function(params_Lo, T_dat, Nmax, w_EN){
  params_Lo <- as.list(c(params_Lo))
  T_A <- params_Lo$T_A
  T_0 <- params_Lo$T_0
  T_AL <- params_Lo$T_AL
  T_L <- params_Lo$T_L
  T_H <- params_Lo$T_H
  T_AH <- params_Lo$T_AH
  JENAM <- params_Lo$JENAM
  K_N <- params_Lo$K_N
  B <- 0.334
  M_V <- 0.01
  
  T_dat <- T_dat + 273.15 #if converstion to Kelvin is needed
  C_T <- exp((T_A/T_0)-(T_A/T_dat)) * (1+exp((T_AL/T_0)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_0))) * ((1+exp((T_AL/T_dat)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_dat)))^-1) #C_T has no units
  N <- seq(0, Nmax, 0.00000001) #array from 0 to Nmax stepped by an interval of 0.01
  
  JENAM <- JENAM * M_V/B #convert from mol N / molM_V / h to mol N / g DW / hour
  J_EN_AM <-  JENAM * C_T #temp correct max assimilation rate of N03
  #2.1 Lorena specific assimilation rate of N
  J_EN_A <- J_EN_AM*(N/(N+K_N)) #this is what the function outputs
  return(cbind(N, J_EN_A))
}

