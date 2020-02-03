# SugarKelpDEB
This repository contains code associated with the manuscript Venolia et al. (in press).
Any questions about this code or data should be directed towards Celeste Venolia (celestevenolia@gmail.com) 
and Austin Humphries (humphries@uri.edu).

**Code files key:**

  *KelpDEB_Run.R*: The runfile for the Kelp DEB model which contains most of the plot codes.
  
  *KelpDEB_model.R*: The Kelp DEB model turned into a function using the R package deSolve.
  
  *SolveR_R.R*: The Newton-Raphson solver for the net specific growth rate (r) in the Kelp DEB model.
  
  *Photosynthesis_Calibration.R*: Function code used to calibrate the photosynthesis parameters for the Kelp DEB model.
  
  *N_uptake_Calibration.R*: Function code used to calibrate the N parameters for the Kelp DEB model.
  
  *SensitivityAnalyses_FME.R*: Function code and run code for sensitivity analysis of the Kelp DEB parameters using the R package FME.
  
**Field data key:**

*Dredge_Y1_hobo.csv*: Temperature data for Dredge (Pt Judith Pond S) in 2017-2018

*Dredge_Y2_HoboTempLight.csv*: Temperature data for Dredge (Pt Judith Pond S) in 2018-2019

*RomePoint_Y1_hobotemp.csv*: Temperature data for Rome Pt (Narragansett Bay S) in 2017-2018

*RomePt_Y2_HoboTempLight.csv*: Temperature data for Rome Pt (Narragansett Bay S) in 2018-2019

*Sled_Y1_TempLogger2.csv*: Temperature data for Sled (Pt Judith Pond N) in 2017-2018

*Sled_Y2_HoboLightTemp.csv*: Temperature data for Sled (Pt Judith Pond N) in 2018-2019

*WaterSampleAnalysis2Y1.csv*: Nitrate and nitrite concentration from all sites in 2017-2018

*WaterSamplesY2.csv*: Nitrate and nitrite concentration from all sites in 2018-2019

*Wickford_Y1_hobo.csv*: Temperature data for Wickford (Narragansett Bay N) in 2017-2018

*Wickford_Y2_HoboLightTemp.csv*: Temperature data for Wickford (Narragansett Bay N) in 2018-2019

*Year1kelpdata.csv*: Kelp blade length growth from all sites in 2017-2018

*Year2kelpdata.csv*: Kelp blade length growth from all sites in 2018-2019

**Literature data key:** for references see Venolia et al. (in press)

 *BrentonPoint_Segarra2002CarbonData.csv*: Dissolved inorganic carbon concentration used for Narragansett Bay sites (Segarra, 2002)
 
 *EspinozaChapman1983_Nuptake_9C_StMargaretsBay.csv*: Nitrate uptake data for kelp from St Margarets Bay with an experimental temperature of 9 degrees C from Espinoza and Chapman (1983)
 
 *EspinozaChapman1983_Nuptake_18C_StMargaretsBay.csv*: Nitrate uptake data for kelp from St Margarets Bay with an experimental temperature of 18 degrees C from Espinoza and Chapman (1983)
 
 *Johansson2002.csv*: Oxygen production rate for kelp from Johansson and Snoeijs (2002)
 
 *Ninigret_EPA_DIC.csv*: Dissolved inorganic carbon concentration used for Pt Judith sites (J. Grear, unpublished data)
 
 *NOAASurfaceIrradiance.csv*: radiative forcing from the North American Regional Reanalysis (Mesinger et al., 2006) used to estimate photosynthetically active radiation
 
 
