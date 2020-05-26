########################################################################################################
#This code was created by Celeste Venolia in August 2018 for Venolia et al. (in press)
#deSolve structured sugar kelp DEB model
########################################################################################################

#Import libraries
library(deSolve)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#The deSolve structured model
#note that the time, state, and parameters in this function definition should always remain in this format
#regardless of how these pieces of info are labelled
rates_Lo <- function(t, state, parameters) { 
  
  with(as.list(c(state, parameters)), { 
    
    #set-up equations (to calculate values necessary for the differential equations)
    #Temperature correction:
    C_T <- exp((T_A/T_0)-(T_A/T_field(t))) * (1+exp((T_AL/T_0)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_0))) * ((1+exp((T_AL/T_field(t))-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_field(t))))^-1)
    
    J_EN_AM <-  JENAM * C_T #temperature correction of max assimilation rate of N
    #Specific assimilation rate of N
    J_EN_A <- J_EN_AM*(N_field(t)/(N_field(t)+K_N))  #J_EN_A units: mol N/molV/h
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


