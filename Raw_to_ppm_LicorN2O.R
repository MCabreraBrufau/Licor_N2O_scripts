#Injection peaks to ppm#

# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Nov 2024"
# https://github.com/MCabreraBrufau/Licor_N2O_scripts
# ---


#Cosas para hacer: 

#Implementar collección de linea base: "lab_baseline" y "loop_baseline" en protocolo de lab, extraer para cada dia estos datos, guardar valores medios, desviacion, nobs. 

#Implementar en el protocolo dos inyecciones de estandar interno: 6ppm_1ml y lab_1ml. Para calidad, podremos utilizar 6ppm_1ml como drift, y lab_1ml vs lab_baseline como estandar interno. 2 puntos de calibracion diarios 6ppm y aire ambiental (autocalibracion)

#crear un script separado para obtener injection_sequence (label, start, stop) para cada rawfile, hacerlo editable, luego importar desde este script los datos con los label corregidos/eliminados/anadidos

#crear loop para cada rawfile, comprobando primero si los datos estan extraidos ya (no ejecutar para todos los rawfile, sino solo para los que no tengan ya los datos extraidos, en base a un metadata_file ) 


#Description: this script takes raw-files from Li-COR 7820 containing discrete injections, corrected injection_sequences (with label, start and stop) and calculates N2O concentration (in ppm) along with signal-to-noise ratio for each injection. It also generates inspection plots (baseline correction & integration) and stores the results in csv format. 

#This script uses the calibration curve obtained in Licor_cal_n2o.R : i.e. #ppm = ((1/203.05)*(peaksum-5.8236))/ml


# ---- Directories ----

#Root
dropbox_root <- "C:/Users/Miguel/Dropbox" # You have to make sure this is pointing to the write folder on your local machine

#Data folders
folder_raw <- paste0(dropbox_root,"/Licor_N2O/Rawdata") #contains unedited files downloaded from licor

folder_injections<- paste0(dropbox_root,"/Licor_N2O/Injection sequences") #Contains csvs with startstop times of injections and their corresponding labels, corrections should be made manually when needed (editting the csvs and re-saving with "_correct" sufix)

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


#Modified function Licor N2o
read_Licor_n2o <- function(file){
  message(paste0("reading ",file))
  data_raw <- read_lines(file)
  prefix <- substr(data_raw, start = 1,stop = 5) # isolate first 5 characters for each line
  
  # find line corresponding to headers
  headers <- unlist(strsplit(data_raw[which(prefix == "DATAH")], "\t"))
  units <- unlist(strsplit(data_raw[which(prefix == "DATAU")], "\t"))
  
  data <- read.delim(file, sep = "\t", header = F, skip = which(prefix == "DATAU"), na.strings = "nan")
  names(data) <- headers
  
  my_data <- data.frame(date = data$DATE,
                        UTCtime = data$TIME,
                        unixtime = data$SECONDS,
                        H2O = data$H2O,
                        N2O = data$N2O,
                        # CO2 = data$CO2,
                        # CH4 = data$CH4/1000, #ppm
                        Press = data$CAVITY_P,
                        label = data$REMARK)
  
  return(my_data)
}



