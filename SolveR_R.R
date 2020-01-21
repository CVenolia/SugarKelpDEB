###################################################################################################
#R conversion of the DEB tool alga sgr2.m function (DEBtool package – scripts freely available (https://www.bio.vu.nl/thb/deb/deblab/debtool/DEBtool_M/manual/index.html))
#Created August 2018 for Venolia et al. (in press)
#This function is used in the R file KelpDEB_Run.R
###################################################################################################

SolveR_R <- function(m_E, k_E, J_EM, y_EV, J_VM, r0){
#statements
r <- r0 #set the r equal to zero the first time the loop runs
M <- m_E / y_EV #help quantities #unit: mol i/mol Ei
sM <- sum(M)  #sum of M necessary because m_E and y_EV are both vector objects with 2 elements #unit: mol i/mol Ei
i <- 0 #set intial number of interations to zero (the loop increases i by 1 with each run)
n <- 20 #max of interations was set at 20
info <- 1 #If there is no convergence info is set to 0
f <- 1 #allows the while loop to begin #initiate norm; make sure that iteration procedure is started

#this while loop continues to run while the number of iterations is less than 20 and for the function f: f^2 > 1e-15
while(f^2 > 1e-15 && i < n) { #test norm
  #structure specific catabolic fluxes J_EC_loop <- m_E * (k_E - r) #units: mol i/molM_V/h
  J_EC_loop_C <- m_E[1] * (k_E[1] - r) #units: mol i/molM_V/h
  J_EC_loop_N <- m_E[2] * (k_E[2] - r) #units: mol i/molM_V/h
  J_EC_loop <- c(J_EC_loop_C, J_EC_loop_N)
  #Is the catabolic flux large enough to fill maintenance requirements
  JEM_loop_C <- max(1e-8, min(J_EM[1], J_EC_loop[1])) #unit: molEC/molM_V/h or mol C/molM_V/h depending on what quanity is smaller
  JEM_loop_N <- max(1e-8, min(J_EM[2], J_EC_loop[2])) #unit: molEN/molM_V/h or mol N/molM_V/h depending on what quanity is smaller
  JEM_loop <- c(JEM_loop_C, JEM_loop_N)
  #relevant if structure is required to fill maintenance requirements
  JVM_loop <- J_VM * (1 - JEM_loop / J_EM) #unit: 1/h
  #growth fluxes
  R_C <- max(1e-6, (J_EC_loop[1]-JEM_loop[1])/y_EV[1]) #unit: 1/h
  R_N <- max(1e-6, (J_EC_loop[2]-JEM_loop[2])/y_EV[2]) #unit: 1/h
  R <- c(R_C, R_N) 
  sR <- R_C + R_N #unit: 1/h
  s <- sum(1/R)-1/sR #units: h
  f <- r + sum(JVM_loop)-(1/s) #norm #units: 1/h
  df <- 1+(sum(M/(R^2))-sM/(sR^2))/(s^2)  #d/dr f
  r <- r-(f/df) #new step in NR-procedure #units: 1/h
  i <- i+1 #mark the increase in iteration
}

#triggered to provide no convergence message when 20 interations result in no convergence
if(i==n){
  info <- 0  # "no convergence"
  text <- paste("no convergence of SolveR_R in ", num2str(n, 0), " steps")
  print(text)
  break
}

#reject fluxes
JER_loop_C <- max(0 , J_EC_loop_C-JEM_loop_C-y_EV[1]*(r+sum(JVM_loop))) #units: mol Ei/molM_V/h
JER_loop_N <- max(0, J_EC_loop_N-JEM_loop_N-y_EV[2]*(r+sum(JVM_loop))) #units: mol Ei/molM_V/h
JER_loop <- c(JER_loop_C, JER_loop_N)

return(c(r, J_EC_loop, JEM_loop, JVM_loop, JER_loop, info)) #function to return these values
}

