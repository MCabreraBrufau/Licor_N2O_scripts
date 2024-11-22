#Peaks to ppm 

# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Nov 2024"
# https://github.com/MCabreraBrufau/Licor_N2O_scripts
# ---

#Description: this script uses integrated_injections files produced in the Raw_to_peaks_LicorN2O.R script and calculates ppm for each peak based on the calibration curve and volume injected. It outputs ppm data for each peak and for each sample

#This script uses the calibration curve obtained in Licor_cal_n2o.R : 
#i.e. #ppm = ((1/203.05)*(peaksum-5.8236))/ml_injection

intercept<- 5.8236
slope<- (1/203.05)


# ---- Directories ----

#Root
dropbox_root <- "C:/Users/Miguel/Dropbox" # You have to make sure this is pointing to the write folder on your local machine

#Data folders
folder_raw <- paste0(dropbox_root,"/Licor_N2O/Rawdata") #contains unedited files downloaded from licor

folder_mapinjections<- paste0(dropbox_root,"/Licor_N2O/Map_injections") #Contains csvs with startstop times of injections and their corresponding labels, corrections should be made manually when needed (editting the csvs and re-saving with "corrected_" prefix)

folder_plots<-  paste0(dropbox_root,"/Licor_N2O/Integration plots") #One pdf per dayofinjections (auto-name from rawfile name), plots of each injection sequence (baseline correction & integration)

folder_results<- paste0(dropbox_root,"/Licor_N2O/Results_ppm")



# ---- packages & functions ----
library(tidyverse)
library(readxl)
library(lubridate)
# library(ptw)
# library(pracma)
library(stringr)
library(ggpmisc)


#Get extracted data
integratedfiles<- list.files(path = folder_results, pattern = "^integrated_injections_")
ppmfiles<- list.files(path = folder_results, pattern = "^ppm_samples_")

#Select integratedfiles without ppm data
integratedtoppm<- gsub(".csv","",gsub("integrated_injections_","",integratedfiles[
  !gsub(".csv","",gsub("integrated_injections_","",integratedfiles))%in%gsub(".csv","",gsub("ppm_samples_","",ppmfiles))]))#  integrated files "rawcode" without corresponding ppmfiles "rawcode"
  

for (i in integratedtoppm){
  
  #Load integrated peaks of integratedfile i
  int<- read.csv(paste0(folder_results,"/","integrated_injections_",i,".csv"))
  
  #Calculate ppm of each peak
  peak_ppm<- int %>% 
    separate(peak_id, into = c("sample", "ml_injected","peak_no"), sep = "_",remove = F) %>% 
    mutate(ml_injected=as.numeric(ml_injected), n2o_ppm=(slope*(peaksum-intercept))/ml_injected) %>% 
    select(dayofanalysis, sample, ml_injected, peak_id, n2o_ppm, unixtime_ofmax) %>% 
    mutate(datetime=as.POSIXct(unixtime_ofmax))
  
  #Save ppm of peaks
  write.csv(peak_ppm, file = paste0(folder_results, "/","ppm_peaks_",i,".csv"), row.names = F)
  
  sample_ppm<- peak_ppm %>% 
    group_by(dayofanalysis,sample, ml_injected) %>% 
    summarise(n2o_ppm_avg=mean(n2o_ppm,na.rm=T),
              n2o_ppm_sd=sd(n2o_ppm,na.rm=T),
              n2o_ppm_cv=n2o_ppm_sd/n2o_ppm_avg,
              nobs=n(),
              datetime=as.POSIXct(min(unixtime_ofmax)),
              .groups="keep") %>% ungroup()
  
  write.csv(sample_ppm, file = paste0(folder_results, "/","ppm_samples_",i,".csv"), row.names = F)
  
}


