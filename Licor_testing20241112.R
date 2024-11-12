##Licor N2O peak process calibration##


# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Oct 2024"
# https://github.com/MCabreraBrufau/Licor_N2O_scripts
# ---

#Description: 
# This script contains the necessary steps to obtain the calibration factor of Licor 7820 for N2O for discrete sample injections. Quality checks for the processing are also included.

#Basic steps from rawlicorfiles to peak-areas:
  # 1. Import raw-data (using function addapted from camille's scripts)
  # 2. Correct labels: either adding remarks based on time of injection and removing remarks that are not correct.
  # 3. Get starstop times of injections per label
    #FOR every injection sequence in startstop:
      # 4. Do baseline correction 
      # 5. Detect peaks
      # 6. Integrate peaks
      # 7. Plot results (original signal, base-corrected, peaks detected, windows of peaks & integrated peaks)
  # 8. Decide on linear range of injection volumes    
  # 9. Calculate calibration factor and quality indexes


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


# ---- Directories ----

#Root
dropbox_root <- "C:/Users/Miguel/Dropbox" # You have to make sure this is pointing to the write folder on your local machine

#Rawdata
datapath_licor <- paste0(dropbox_root,"/Licor_N2O/Rawdata")


# ---- 1.Import ----

#Import data from cal day
calfile<- paste0(datapath_licor,"/TG20-01377-2024-11-12T070000.data")

a<- read_Licor_n2o(calfile)


# ---- 2.Correct labels ----

#Wrong remarks until 10:28:30, delete them
a[between(a$unixtime,
          a[a$UTCtime=="08:55:00",]$unixtime,
          a[a$UTCtime=="10:28:30",]$unixtime),]$label<- ""


#Substitute 1.2ppm-b_1ml with 1.2ppm-b-i_1ml 
a[a$label=="1.2ppm-b_1ml",]$label<- "1.2ppm-b-i_1ml"

#Remove weird peak    
a[a$label=="6ppm-b_1ml",]$label<-""

#Some wrong injections:
a[between(a$unixtime,
          a[a$UTCtime=="11:33:20",]$unixtime,
          a[a$UTCtime=="11:38:20",]$unixtime),]$label<- ""


unique(a$label)  

# ---- 3.Get startstop times----

#Keep all injection sequences identified by label and mark their start and end
startend_inj<- a %>% 
  filter(label!="") %>% 
  group_by(label) %>%
  summarise(
    start_inj = first(unixtime),
    end_inj = last(unixtime),
    .groups = 'drop'  # To ungroup after summarization
  )


#Object for saving plots
plotspeak <- list()

for (idinj in startend_inj$label){
  
  # Subset full dataset for expanded start and end of injection sequence "idinj"
  startinj<- startend_inj[startend_inj$label==idinj,]$start_inj
  endinj<- startend_inj[startend_inj$label==idinj,]$end_inj
  
  dat<- a[between(a$unixtime, startinj,endinj),]  
  # -----1. Base-correct----
  #Base-correct injection sequence, using asymetric least-square. 
  dat_bc<-dat %>% 
    mutate(N2Obc=baseline.corr(N2O,lambda=1e5, p=0.0001))
  

  # -----2. Peak-max detection ----
  #Find local maxima in sequence and add max_id (label_1,label_2,...) : 
  #Criteria for local maximum:
  # at least 1 increase before and 1 decrease afterwards
  # minimum peak height to be detected is > 1/5 of maximum point in all remark
  # at leas 5 points between localmaxima
  
  dat_peakid <- dat_bc %>%
    mutate(is_localmaxn2o = ifelse(row_number() %in% findpeaks(N2Obc, 
                                                               minpeakheight = max(N2Obc,na.rm = T)/5, 
                                                               nups=1, ndowns=1,
                                                               minpeakdistance = 5)[, 2], TRUE, FALSE)) %>%
    mutate(peak_id = ifelse(is_localmaxn2o, paste0(label,"_",cumsum(is_localmaxn2o)), NA)) %>%  #Add unique code for local maxima 
    ungroup()
  
  # -----3. Define peak width ----
  #Identify peakwindow: 
  #Consider peakwindow as max height + 4 leading and 7 trailing points. (i.e. peak width == 9points), 
  
  dat_peakwindow <- dat_peakid %>%
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
  
  # -----4. Integrate peaks ----
  #Summarise each peak_id (sum, max, unixtime)
  dat_integrated<- dat_peakwindow %>% 
    filter(!is.na(peak_id)) %>% #keep only data of peaks
    group_by(label, peak_id) %>% 
    summarise(peaksum=sum(N2Obc), peakmax=max(N2Obc), unixtime_ofmax=unixtime[N2Obc==max(N2Obc)],
              peaksum_nobasec=sum(N2O)) %>% 
    ungroup()
  
  # -----5. Plot results ----
  #Create a plot for inspection:
  p<- ggplot(dat_peakwindow, aes(x=unixtime, y=N2Obc))+
    geom_point(aes(col=peak_id))+
    geom_line(aes(col="1_base-correct"))+
    geom_point(data = dat_integrated, aes(x=unixtime_ofmax, y=peaksum, col="integral"))+
    geom_line(data = dat_peakwindow, aes(x=unixtime, y=N2O, col="0_raw"))+
    geom_point(data = dat_peakwindow, aes(x=unixtime, y=N2O, col="0_raw"))+
    ggtitle(idinj)
  
  # Store each plot in the list
  plotspeak[[idinj]] <- p

  #Add results to objects all_peaks and all_data
  if (idinj==startend_inj$label[1]) {
    all_peaks<-dat_integrated 
    all_data<-dat_peakwindow
  }else{
    all_peaks<-rbind(all_peaks,dat_integrated)
    all_data<-rbind(all_data,dat_peakwindow)
  }
  
  #Remove intermediate objects in last run of loop
  if (idinj==startend_inj$label[length(startend_inj$label)]){
    rm(startinj,endinj,dat,dat_bc,dat_integrated,dat_peakid,dat_peakwindow,p,idinj)
  }
}




#plot every injection sequence and their integrals: 
pdf(file = paste0(datapath_licor,"/plots_20241112_testing.pdf"))  # Open PDF device

# Loop through the list of plots and print each plot
for (plot_name in names(plotspeak)) {
  print(plotspeak[[plot_name]])
}

dev.off()  # Close the PDF device




ppmpeaks<- all_peaks %>% 
  mutate(n2o=((1/203.05)*(peaksum-5.8236))/1) %>% 
  filter(!peak_id%in%c("1.2ppm-e-m_1ml_5","6ppm-e-m_1ml_4","6ppm-e-m_1ml_5"))


ppmpeaks %>% 
group_by(label) %>% 
summarise(avg=mean(n2o), SD=sd(n2o), CV=(SD/avg)*100, nobs=n()) 

ppmpeaks %>% 
  group_by(label) %>% 
  summarise(avg=mean(n2o), SD=sd(n2o), CV=(SD/avg)*100, nobs=n()) %>% 
  ggplot(aes(x=0, y=CV))+
  geom_boxplot()+
  geom_point()+
  ggtitle("CV % basado en 13 muestras de 5-7 injecciones")
  


ppm_comparisons<-
ppmpeaks %>% 
  filter(!grepl("^Pu",label)) %>% 
  filter(!peak_id%in%c("6ppm-e-m_1ml_4","6ppm-e-m_1ml_5","1.2ppm-e-m_1ml_5")) %>% #remove last inections, exetainer not full
  separate(label, into = c("concentration","origin","user"),sep = "-",remove = F) %>% 
  mutate(user=gsub("_1ml","",user)) %>% 
  mutate(user=case_when(user=="m"~"Miguel",
                        user=="i"~"Ilenia"),
         origin=case_when(origin=="b"~"Bolsa",
                          origin=="e"~"Exetainer")) %>% 
  group_by(concentration) %>% 
  mutate(avg_n2o=mean(n2o))

ggplot(ppm_comparisons, aes(x=origin, y=n2o, col=user))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~concentration, scales="free")


ggplot(ppm_comparisons, aes(x=origin, y=100*(n2o-avg_n2o)/avg_n2o, col=user))+
  geom_boxplot()+
  facet_wrap(~concentration, scales="free")+
  scale_y_continuous(name="Percent difference from mean")


  ggplot(aes(x=label, y=(100*(n2o-mean(n2o))/mean(n2o)), fill = label))+
  geom_boxplot()+
  geom_point()+
  # geom_label(aes(label=peak_id))+
  scale_y_continuous("Difference from mean (%)")


ppmpeaks %>% 
  filter(grepl("^6ppm-",label)) %>% 
  filter(!peak_id%in%c("6ppm-e-m_1ml_4","6ppm-e-m_1ml_5")) %>% #remove last inection, exetainer not full
  ggplot(aes(x=label, y=100*((n2o-mean(n2o))/mean(n2o))))+
  geom_boxplot()+
  geom_point()+
  # geom_label(aes(label=peak_id))
  scale_y_continuous("Difference from mean (%)")



#Purga: es suficiente con 3 veces x 0.2ml 
ppmpeaks %>% 
  filter(grepl("^Pu",label)) %>% 
  filter(!label%in%c("Pu_aire_lab")) %>% 
  ggplot(aes(x=label, y=n2o, fill = label))+
  geom_boxplot()+
  geom_point()
  


  #Comparison peaks aire lab vs aire lab continuous:
  
  ppmpeaks %>% 
    filter(grepl("^Pu",label)) %>% 
    filter(label%in%c("Pu_aire_lab")) %>% 
    ggplot(aes(x=label, y=n2o))+
    geom_boxplot()+
    geom_point()+
    geom_hline(yintercept = 0.34)+ #Concentracion lab a ojo del LiCor antes de encender botella aire.
    scale_y_continuous("N2O (ppm)", limits = c(0,0.4))
  
  
  ppm_comparisons %>% 
    group_by(concentration, user, origin) %>% 
    summarise(avg=mean(n2o), SD=sd(n2o), CV=(SD/avg)*100, nobs=n())   
  
  ppm_comparisons %>% 
    mutate(n2o_norm=n2o/avg_n2o) %>% 
    ggplot(aes(x=origin, y= n2o_norm, col=user))+
    geom_boxplot()+
    geom_point()

  model <- aov(n2o/avg_n2o ~ user * origin, data = ppm_comparisons)
  summary(model)# origin is only slightly signnificant
  
  