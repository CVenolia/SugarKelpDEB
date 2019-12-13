###################################################################################################
#Celeste's conversion of the DEB tool alga function modified by Romain Lavaud to be called SolveR
#Created August 2018
#"R:" comments are from Romain
#This function is used in the R file KelpDEBCTV.R
###################################################################################################

SolveR_R <- function(m_E, k_E, J_EM, y_EV, J_VM, r0){
#statements
r <- r0 #set the r equal to zero the first time the loop runs
M <- m_E / y_EV #R: help quantities #unit: mol i/mol Ei
sM <- sum(M)  #sum of M necessary because m_E and y_EV are both vector objects with 2 elements #unit: mol i/mol Ei
i <- 0 #set intial number of interations to zero (the loop increases i by 1 with each run)
n <- 20 #max of interations was set at 20
info <- 1 #If there is no convergence info is set to 0, I have this set to output with the graphs
f <- 1 #CTV: allows the while loop to begin #R: initiate norm; make sure that iteration procedure is started

#this while loop continues to run while the number of iterations is less than 20 and for the function f: f^2 > 1e-15
while(f^2 > 1e-15 && i < n) { #R: test norm
  #2.5 Lorena #structure specific catabolic flux
  #J_EC_loop <- m_E * (k_E - r) #units: mol i/molM_V/h
  J_EC_loop_C <- m_E[1] * (k_E[1] - r) #units: mol i/molM_V/h
  J_EC_loop_N <- m_E[2] * (k_E[2] - r) #units: mol i/molM_V/h
  J_EC_loop <- c(J_EC_loop_C, J_EC_loop_N)
  #print(paste0("J_EC_loop is ", J_EC_loop))
  #I believe this comes from Lorena 2.11
  JEM_loop_C <- max(1e-6, min(J_EM[1], J_EC_loop[1])) #unit: molEC/molM_V/h or mol C/molM_V/h depending on what quanity is smaller
  JEM_loop_N <- max(1e-6, min(J_EM[2], J_EC_loop[2])) #unit: molEN/molM_V/h or mol N/molM_V/h depending on what quanity is smaller
  JEM_loop <- c(JEM_loop_C, JEM_loop_N)
  #print(paste0("JEM_loop is ", JEM_loop))
  #This doesn't quite match the rest of Lorena 2.11, but it is along the same vein
  JVM_loop <- J_VM * (1 - JEM_loop / J_EM) #unit: 1/h
  #Above this line in the while statement was written in Matlab in a for loop
  
  #inside the first parantheses for Lorena equation 2.12 for Carbon
  R_C <- max(1e-6, (J_EC_loop[1]-JEM_loop[1])/y_EV[1]) #unit: 1/h
  #inside the first parantheses for Lorena equation 2.12 for Nitrogen
  R_N <- max(1e-6, (J_EC_loop[2]-JEM_loop[2])/y_EV[2]) #unit: 1/h
  R <- c(R_C, R_N) 
  #similar question to above about adding N and C units together
  sR <- R_C + R_N #unit: 1/h
  #the part of Lorena equation 2.12 inside the square brackets
  s <- sum(1/R)-1/sR #units: h
  #Lorena equation 2.12
  f <- r + sum(JVM_loop)-(1/s) #R: norm #units: 1/h
  #print(paste0("f is ", f^2))
  df <- 1+(sum(M/(R^2))-sM/(sR^2))/(s^2) #R: d/dr f #units: (mol i/mol Ei/(1/h)^2)/d^2
  r <- r-(f/df) #R: new step in NR-procedure #units: 1/h (f/df units don't quite cancel cleanly)
  #print(r)
  #print(paste("f is", f))
  #print(paste("df is", df))
  #print(paste("J_EM", JEM_loop[1]))
  #print(paste("J_EC_loop", J_EC_loop[1]))
  #print(paste("cat-main =", J_EC_loop[1]-JEM_loop[1]))
  i <- i+1 #mark the increase in iteration
  #testing <- paste("iteration number", num2str(i, 0), "r estimation is", num2str(r, 4)) #added by CTV to observe the outcome of SolveR
  #print(testing)
}

if(i==n){
  info <- 0  # "no convergence"
  text <- paste("no convergence of SolveR_R in ", num2str(n, 0), " steps")
  print(text)
  break
}

#2-vector with rejected flux of reserve #Romain was going to think more about this part!: y_EV*(r+sum(JVM_loop)
JER_loop_C <- max(0 , J_EC_loop_C-JEM_loop_C-y_EV[1]*(r+sum(JVM_loop))) #units: mol Ei/molM_V/h
JER_loop_N <- max(0, J_EC_loop_N-JEM_loop_N-y_EV[2]*(r+sum(JVM_loop))) #units: mol Ei/molM_V/h
JER_loop <- c(JER_loop_C, JER_loop_N)

return(c(r, J_EC_loop, JEM_loop, JVM_loop, JER_loop, info)) #object
}

