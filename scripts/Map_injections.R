#Map injections rawdata


# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Nov 2024"
# https://github.com/MCabreraBrufau/Licor_N2O_scripts
# ---


#Description: This script creates a map_injection csv for every rawdata file. These are stored in the Map_injections folder and will be manually edited to correct label text, adapt time of labels or add labels that were not written in the data at the time of collection. 


# ---- Directories ----

#Root
dropbox_root <- "C:/Users/Miguel/Dropbox" # You have to make sure this is pointing to the write folder on your local machine

#Data folders
folder_raw <- paste0(dropbox_root,"/Licor_N2O/Rawdata") #contains unedited files downloaded from licor

folder_mapinjections<- paste0(dropbox_root,"/Licor_N2O/Map_injections") #Where csv files with the injection details will be saved. 

# ---- packages & functions ----
library(tidyverse)

#Import functions of repo 
repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


#Check rawdata files and map_injection files

#raw_files are named as: paste("TG20-01377-", year-month-day, "-T",hourstartexport, ".data")
#maps_done will be named as raw_map_injection_rawfile.csv, the rawfile part will not end in .data

#List maps that are already created in folder_mapinjections
maps_done<- list.files(path=folder_mapinjections, pattern = "^raw_map_injection")

#List raw files present in folder_raw
raw_files<- list.files(path = folder_raw, pattern =  "^TG20-01377")


#Get raw files without corresponding map injection
raw_files_withoutmap<- raw_files[!gsub(".data", "",raw_files)%in%gsub(".csv", "", gsub(pattern = "raw_map_injection_","",maps_done))]


#Collect Tstart Tend and labels for all unique remarks of every raw_file_withoutmap
#Save these details in csv files named "raw_map_injection_"[rawfilename without".data"].csv 
for (i in raw_files_withoutmap){
  
  a<- read_Licor_n2o(paste0(folder_raw,"/",i))%>% 
    group_by(label) %>% 
    summarise(date=first(date),Tstart=first(UTCtime), Tend =last(UTCtime)) %>% 
    arrange(Tstart) %>% 
    mutate(rawfile=i) %>% 
    select(date, Tstart, Tend, label, rawfile)
  
  write.csv(a,file = paste0(folder_mapinjections,"/raw_map_injection_", gsub(".data","",i), ".csv"),row.names = F)
  
}

rm(a,i, raw_files,raw_files_withoutmap,maps_done)
