This repository contains the scripts necessary to semi-automatically obtain N2O concentrations of discrete samples analysed with the open-loop method.
Key steps:
1. Import raw data and produce injections map based on Li-COR remarks (Script: "Map injections.R")
2. Correct injection map based on excell with lab anotations (Manual step)
3. Import raw data and corrected injection maps, automatically integrate peaks and calculate N2O concentrations. Produce integration plots and associated quality parameters of injections. 
