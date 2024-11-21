#Injection peaks to ppm#

# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Nov 2024"
# https://github.com/MCabreraBrufau/Licor_N2O_scripts
# ---


#Cosas para hacer: 

#Implementar collección de linea base: "lab_baseline" y "loop_baseline"  extraer para cada dia estos datos, guardar valores medios, desviacion, nobs. 

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
ppmfiles<- list.files(path = folder_results, pattern = "^ppm_of_injections_")


#Select code of rawfiles with corresponding mapscorrect but without ppmfiles
rawtointegrate<- gsub(".data","",rawfiles[
  gsub(".data","",rawfiles)%in%gsub(".csv","",gsub("corrected_map_injection_","",mapscorrect))& #Match raw with maps
    !gsub(".data","",rawfiles)%in%gsub(".csv","",gsub("ppm_of_injections_","",ppmfiles)) #Not match raw with ppmfiles
])





###2. Integration loop####

#Loop over raw files
for (i in rawtointegrate){
  
  message(paste("Integrating peaks from",i))
  
  #Get date of analysis 
  dayofanalysis<-  str_extract(pattern = "[0-9]{4}-[0-9]{2}-[0-9]{2}", string = i)
  
  #Import data from rawfile
  raw_data<- read_Licor_n2o(paste0(folder_raw,"/",i,".data"))
  
  #Import corrected map of injections
  mapinj<- read.csv(paste0(folder_mapinjections,"/","corrected_map_injection_",i,".csv")) %>% 
    filter(!is.na(label_correct)) %>% 
    filter(label_correct!="") %>% 
    mutate(date=as.Date.character(date,format="%d/%m/%Y"))
  
  #Create tables where baseline and injections will be saved
  
  #Initialize data frame for injections
  A<- data.frame(
    dayofanalysis=character(),
    label = character(),
    peak_id = character(),
    peaksum = double(),
    peakmax = double(),
    unixtime_ofmax = double(),
    raw_peaksum = double(),
    peakSNR = double())
  
  #Initialize data frame for baselines
  B<- data.frame(
    dayofanalysis=character(),
    label = character(),
    base_avg = double(),
    base_sd = double(),
    base_cv = double(),
    base_n = integer(),
    stringsAsFactors = FALSE
  )

  #Initialize list of plots to save integration plots
  plotspeak <- list()
  
  
    
  #loop over different labels of rawfile i
  for (inj in mapinj$label_correct){
   
    #Unixstart, Tstart_correct from mapinj in unix time format
    unixstart<- as.numeric(as.POSIXct(paste(mapinj[mapinj$label_correct==inj,]$date,mapinj[mapinj$label_correct==inj,]$Tstart_correct), tz = "UTC"))
    
    #Unixend, Tend_correct from mapinj in unix time format
    unixend<- as.numeric(as.POSIXct(paste(mapinj[mapinj$label_correct==inj,]$date,mapinj[mapinj$label_correct==inj,]$Tend_correct), tz = "UTC"))
    
    
    #Subset data from injection sequence inj 
    inj_data<- raw_data[between(raw_data$unixtime, unixstart,unixend),]  
   
    #Make sure whole inj_data has the correct label inj
    inj_data$label<- inj
    
    ######2.1. Baselines#####
    if (grepl("baseline", inj)){
      print(paste0('Baseline recording: ',inj))
      
      #calculate descriptive statistics for baseline
      b<- inj_data %>% 
        summarise(base_date=dayofanalysis,
                  label=inj,
                  base_avg= mean(n2o,na.rm = T), 
                  base_sd= sd(n2o,na.rm=T),
                  base_cv=base_sd/base_avg,
                  base_n= n())
      
      #Add baseline statistics to baseline table
      B<- rbind(B,b)
    } 
    
      ###2.2. Injections#####
    else {
      print(paste0("Injection sample: ", inj))
      
      #Detect and integrate peaks, plot results, calculate  baseline SD within label for Signal to Noise ratio

      ##____Base-correction#####
      #Base-correct injection sequence, using asymetric least-square. 
      inj_data<-inj_data %>% 
        mutate(n2o_bc=baseline.corr(n2o,lambda=1e5, p=0.0001))
      
      ##____Peak-max detection#####
      
      #Find local maxima in sequence and add max_id (label_1,label_2,...) : 
      #Criteria for local maximum:
      # at least 1 increase before and 1 decrease afterwards
      # minimum peak height to be detected is > 1/5 of maximum point in all remark
      # at leas 5 points between localmaxima
      
      inj_data <- inj_data %>%
        mutate(is_localmaxn2o = ifelse(row_number() %in% findpeaks(n2o_bc, 
                                                                   minpeakheight = max(n2o_bc,na.rm = T)/5, 
                                                                   nups=1, ndowns=1,
                                                                   minpeakdistance = 5)[, 2], TRUE, FALSE)) %>%
        mutate(peak_id = ifelse(is_localmaxn2o, paste0(label,"_",cumsum(is_localmaxn2o)), NA)) %>%  #Add unique code for local maxima 
        ungroup()
      
      ##____Peak-window selection#####
      #Consider peakwindow as max height + 4 leading and 7 trailing points. (i.e. peak width == 9points), 
      
      inj_data <- inj_data %>%
        mutate(peak_id = map_chr(row_number(), function(idx) {
          seq_start <- max(1, idx - 7)
          seq_end <- min(n(), idx + 4)
          
          # Check for peak_id in the surrounding window
          surrounding_codes <- peak_id[seq(seq_start, seq_end)]
          
          # Return the peak_id if it's available, otherwise return NA
          if (any(!is.na(surrounding_codes))) {
            return(first(na.omit(surrounding_codes)))  # Use the first valid peak_id found
          } else {
            return(NA)
          }
        }))
      
      
      ##____Peak integration#####
      
      #Get baseline noise to outside the peak areas
      baseline_sd<-inj_data %>% 
        filter(is.na(peak_id)) %>%
        summarise(baseline_sd=sd(n2o_bc, na.rm=T)) %>%  ungroup()%>% pull(baseline_sd)
      
      #Summarise each peak_id (peaksum, peakmax, unixtimeofmax, raw_peaksum, peakSNR)
      integrated<- inj_data %>% 
        filter(!is.na(peak_id)) %>% #keep only data of peaks
        group_by(label, peak_id) %>% 
        summarise(peaksum=sum(n2o_bc),
                  peakmax=max(n2o_bc), 
                  unixtime_ofmax=unixtime[n2o_bc==max(n2o_bc)],
                  raw_peaksum=sum(n2o),.groups = "keep") %>%
        mutate(dayofanalysis=dayofanalysis,
               peakSNR=peaksum/(3*baseline_sd)) %>% 
        ungroup()
      
      
      avg_peaksum<- mean(integrated$peaksum)
      sd_peaksum<- sd(integrated$peaksum)
      
        ###____Create integration plots#####
      
      p<-ggplot()+
        geom_point(data=subset(inj_data,!is.na(peak_id)), aes(x=as.POSIXct(unixtime),y=n2o_bc,col="2_peaks"))+
        geom_point(data = integrated, aes(x=as.POSIXct(unixtime_ofmax), y=peaksum, col="3_integrated"))+
        geom_line(data = inj_data, aes(x=as.POSIXct(unixtime), y=n2o_bc, col="1_base-corrected"))+
        geom_line(data = inj_data, aes(x=as.POSIXct(unixtime), y=n2o, col="0_raw"), linetype = 2)+
        scale_y_continuous(name="signal")+
        scale_x_datetime(name="Licor time")+
        labs(col="")+
        ggtitle(paste0(dayofanalysis,", injection: ",inj))+
        theme_bw()+
        # Add label for average peaksum value
        # geom_text(data=integrated, aes(x = as.POSIXct(min(unixtime_ofmax)-50), 
        #               y = min(peaksum)*0.8, 
        #               label = paste("Avg: ", round(avg_peaksum, 2), " ± ", round(sd_peaksum, 2), " (CV= ",round(sd_peaksum/avg_peaksum,2),")" )), color = "black", hjust = 0, 
        #           vjust = 1, 
        #           size = 4, 
        #           fontface = "italic")
        # 
      annotate("text",x = as.POSIXct(min(integrated$unixtime_ofmax)-50), 
               y = min(integrated$peaksum)*0.8, 
               label = paste ("Avg: ", round(avg_peaksum, 2), " ± ", round(sd_peaksum, 2), " (CV= ",round(sd_peaksum/avg_peaksum,2),")" ), color = "black", hjust = 0, 
      vjust = 1, 
      size = 4, 
      fontface = "italic")
      
      # Store each plot in the list
      plotspeak[[inj]] <- p
      
      #Add integrations of inj to injections table
      A<-rbind(A,integrated)
      
      
    }
    
    
  } 
  
  
  #Save baseline statistics of rawfile i 
  write.csv(B,file = paste0(folder_results,"/", "baselines_", i, ".csv"),row.names = F)
  
  #Save concentrations of injections for rawfile i   
  write.csv(A,file = paste0(folder_results,"/", "ppm_of_injections_", i, ".csv"),row.names = F)
  
  #Save plots of integrations: use i for naming convention of pdf
  
  #plot every injection sequence and their integrals: 
  pdf(file = paste0(folder_plots,"/Integrations_",i,".pdf"))  # Open PDF device
  
  # Loop through the list of plots and print each plot
  for (plot_name in names(plotspeak)) {
    print(plotspeak[[plot_name]])
  }
  
  dev.off()  # Close the PDF device
  
  
  
  
}

