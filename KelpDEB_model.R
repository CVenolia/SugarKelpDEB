########################################################################################################
#This code was created by Celeste Venolia in August 2018 based on  
##### Romain Lavaud's Kelp model made in MATLAB
#If there is an "RL:" in the front of comment text this means this is commentary from Romain Lavaud
#The model for Sugar Kelp is adapted from the Lorena et al. 2011 microalgae model
#It utilizes a function currently housed in the file 'SolveR_R.R'
########################################################################################################

#Import libraries
library(deSolve)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#The deSolve structured model
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
    #PUT BACK NC=NC when relevant!
    #m_EC, m_EN, M_V, J_CO2, J_I, J_EC_A, J_EN_A, J_EC_C, J_EN_C, r, J_EN_M, J_EC_M, J_EN_G, J_EC_G, J_G, J_VM, J_EC_R, J_EN_R, N, CO_2, I, B, NC, V, L_w
    return(list(c(dm_ECdt, dm_ENdt, dM_Vdt), r=r, B=B, C_T=C_T, I=I, N=N, J_CO2=J_CO2, J_I=J_I, J_EC_A=J_EC_A, J_EN_A=J_EN_A, J_EC_C=J_EC_C, J_EN_C=J_EN_C, J_EN_M=J_EN_M, J_EC_M=J_EC_M, J_EN_G=J_EN_G, J_EC_G=J_EC_G, J_G=J_G, J_VM=J_VM, J_EC_R=J_EC_R, JER_C_test=JER_C_test, J_EN_R=J_EN_R, JER_N_test=JER_N_test, CO_2=CO_2, L_allometric=L_allometric, J_G_test=J_G_test, J_O=J_O, info=info))
  })
}


