#Injection peaks to ppm#

# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Nov 2024"
# https://github.com/MCabreraBrufau/Licor_N2O_scripts
# ---


#Cosas para hacer: 

#Implementar collección de linea base: "lab_baseline" y "loop_baseline" en protocolo de lab, extraer para cada dia estos datos, guardar valores medios, desviacion, nobs. 

#Implementar en el protocolo dos inyecciones de estandar interno: 6ppm-[1,2,3,..]_1ml y aire-[1,2,3,..]_1ml. Para calidad, podremos utilizar 6ppm_1ml como drift, y lab_1ml vs lab_baseline como estandar interno. 2 puntos de calibracion diarios 6ppm y aire ambiental (autocalibracion)

#crear un script separado para obtener map_injections (label, start, stop) para cada rawfile #HECHO#, hacerlo editable, luego importar desde este script los datos con los label corregidos/eliminados/anadidos


#crear loop para cada rawfile, comprobando primero si los datos estan extraidos ya (no ejecutar para todos los rawfile, sino solo para los que no tengan ya los datos extraidos, en base a un metadata_file ) 


#Description: this script takes raw-files from Li-COR 7820 containing discrete injections, corrected injection_sequences (with label, start and stop) and calculates N2O concentration (in ppm) along with signal-to-noise ratio for each injection. It also generates inspection plots (baseline correction & integration) and stores the results in csv format. 

#This script uses the calibration curve obtained in Licor_cal_n2o.R : i.e. #ppm = ((1/203.05)*(peaksum-5.8236))/ml


# ---- Directories ----

#Root
dropbox_root <- "C:/Users/Miguel/Dropbox" # You have to make sure this is pointing to the write folder on your local machine

#Data folders
folder_raw <- paste0(dropbox_root,"/Licor_N2O/Rawdata") #contains unedited files downloaded from licor

folder_mapinjections<- paste0(dropbox_root,"/Licor_N2O/Map_injections") #Contains csvs with startstop times of injections and their corresponding labels, corrections should be made manually when needed (editting the csvs and re-saving with "corrected_" prefix)

folder_plots<-  paste0(dropbox_root,"/Licor_N2O/Integration plots") #One pdf per dayofinjections (auto-name from rawfile name), plots of each injection sequence (baseline correction & integration)

folder_results<- paste0(dropbox_root,"/Licor_N2O/Results_ppm")#One csv per dayofinjections will be created (auto-name from rawfile name), with individual peak parameters (Original file, exetainer_id, label, peak_id, unixtime of max, Area, SNR, ppmN2O)



# ---- packages & functions ----
library(tidyverse)
library(readxl)
library(lubridate)
library(ptw)
library(pracma)
library(stringr)
library(ggpmisc)

#Load repository functions
repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}



###1. Check data to integrate####

#Get rawfiles
rawfiles<- list.files(path = folder_raw, pattern = ".data")

#Get corrected maps of injections
mapscorrect<- list.files(path = folder_mapinjections, pattern = "^corrected_map_injection_")

#Get extracted data
ppmfiles<- list.files(path = folder_results, pattern = "^ppm_injections_day_")


#Select code of rawfiles with corresponding mapscorrect but without ppmfiles
rawtointegrate<- gsub(".data","",rawfiles[
  gsub(".data","",rawfiles)%in%gsub(".csv","",gsub("corrected_map_injection_","",mapscorrect))& #Match raw with maps
    !gsub(".data","",rawfiles)%in%gsub(".csv","",gsub("ppm_injections_day_","",ppmfiles)) #Not match raw with ppmfiles
])



###2. Integration loop####

#Loop over raw files
for (i in rawtointegrate){
  
  message(paste("Integrating peaks from",i))
  
  #Import data from rawfile
  raw_data<- read_Licor_n2o(paste0(folder_raw,"/",i,".data"))
  
  #Import corrected map of injections
  mapinj<- read.csv(paste0(folder_mapinjections,"/","corrected_map_injection_",i,".csv")) %>% 
    filter(!is.na(label_correct)) %>% 
    filter(label_correct!="") %>% 
    mutate(date=as.Date.character(date,format="%d/%m/%Y"))
  
  #loop over different labels of rawfile i
  for (inj in mapinj$label_correct){
   
    #Get Tstart and Tend
    
    #Unixstart, Tstart_correct from mapinj in unix time format
    unixstart<- as.numeric(as.POSIXct(paste(mapinj[mapinj$label_correct==inj,]$date,mapinj[mapinj$label_correct==inj,]$Tstart_correct), tz = "UTC"))
    
    #Unixend, Tend_correct from mapinj in unix time format
    unixend<- as.numeric(as.POSIXct(paste(mapinj[mapinj$label_correct==inj,]$date,mapinj[mapinj$label_correct==inj,]$Tend_correct), tz = "UTC"))
    
    
    #Subset data from inj 
    inj_data<- raw_data[between(raw_data$unixtime, unixstart,unixend),]  
   
    
    #IF section for baseline process
    
    #calculate descriptive statistics and save them in a different csv 
    
    #ELSEIF section for actual injections
    
    #Detect and integrate peaks, plot results, calculate  baseline SD within label for Signal to Noise ratio
    #Save data in folder_results following the scheme: "ppm_injections_day_" [rawfile without ".data"
     
    
  }
}




