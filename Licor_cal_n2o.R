##Licor N2O peak process##


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
  # 9. Calculate calibration factor and quality indeces


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

#Plots
plotfolder <- paste0(dropbox_root, "/Licor_N2O/Plots")



# ---- 1.Import ----

#Import data from cal day
calfile<- paste0(datapath_licor,"/TG20-01377-2024-10-31T080000_fullday.data")

a<- read_Licor_n2o(calfile)


# ---- 2.Correct labels ----

#6ppm_06ml without remark in licor: injected between 10:55 and 10:59. Add remark:
a[between(a$unixtime,
          a[a$UTCtime=="10:55:00",]$unixtime,
          a[a$UTCtime=="10:59:00",]$unixtime),]$label<- "6ppm_06ml"

#cal6ppm_20ml remark without injections, remmove remark
a[a$label=="cal6ppm_20ml",]$label<- ""

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
plots <- list()


for (idinj in startend_inj$label){
  
  # Subset full dataset for expanded start and end of injection sequence "idinj"
  startinj<- startend_inj[startend_inj$label==idinj,]$start_inj
  endinj<- startend_inj[startend_inj$label==idinj,]$end_inj
  
  dat<- a[between(a$unixtime, startinj,endinj),]  
  # -----1. Base-correct----
  #Base-correct injection sequence
  dat_bc<-dat %>% 
    mutate(N2Obc=baseline.corr(N2O,lambda=1e2, p=0.01))
  
  
  # -----2. Peak-max detection ----
  #Find local maxima in sequence and add max_id (label_1,label_2,...) : 
  #Criteria for local maximum:
  # at least 1 increase before and 1 decrease afterwards
  # minimum peak height to be detected is > percentil 0.9
  # at leas 5 points between localmaxima
  
  dat_peakid <- dat_bc %>%
    mutate(is_localmaxn2o = ifelse(row_number() %in% findpeaks(N2Obc, 
                                                               minpeakheight = quantile(N2Obc, 0.90), 
                                                               nups=1, ndowns=1,
                                                               minpeakdistance = 5)[, 2], TRUE, FALSE)) %>%
    mutate(peak_id = ifelse(is_localmaxn2o, paste0(label,"_",cumsum(is_localmaxn2o)), NA)) %>%  #Add unique code for local maxima 
    ungroup()
  
  # -----3. Define peak width ----
  #Identify peakwindow: 
  #Consider peakwindow as max height + 4 leading and 4 trailing points. (i.e. peak width == 9points), 
  
  dat_peakwindow <- dat_peakid %>%
    mutate(peak_id = map_chr(row_number(), function(idx) {
      seq_start <- max(1, idx - 4)
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
    summarise(peaksum=sum(N2Obc), peakmax=max(N2Obc), unixtime_ofmax=unixtime[N2Obc==max(N2Obc)]) %>% 
    ungroup()
  
  # -----5. Plot results ----
  #Create a plot for inspection:
  p<- ggplot(dat_peakwindow, aes(x=unixtime, y=N2Obc))+
    geom_point(aes(col=peak_id))+
    geom_line()+
    geom_point(data = dat_integrated, aes(x=unixtime_ofmax, y=peaksum, col="integral"))+
    ggtitle(idinj)
  
  # Store each plot in the list
  plots[[idinj]] <- p
  
  #Add results to object all_peaks
  if (idinj==startend_inj$label[1]) {
    all_peaks<-dat_integrated 
  }else{
    all_peaks<-rbind(all_peaks,dat_integrated)
  }
  
  #Remove intermediate objects in last run of loop
  if (idinj==startend_inj$label[length(startend_inj$label)]){
    rm(startinj,endinj,dat,dat_bc,dat_integrated,dat_peakid,dat_peakwindow)
  }
}


#plot every injection sequence and their integrals: 
for (plot_name in names(plots)) {
  print(plots[[plot_name]])
}


all_peaks %>% 
  filter(grepl("^6ppm",label)) %>% 
  filter(!peak_id%in%c("6ppm_03ml_1","6ppm_05ml_6")) %>% 
  mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10) %>% 
  filter(vol<=1) %>%   
  ggplot( aes(x=vol, y=peaksum))+
  geom_point()+
  geom_smooth(aes(col=vol<11),method = "lm", se=F)


peaks_linearity<- 
  all_peaks %>% 
  filter(grepl("^6ppm",label)) %>% 
  filter(!peak_id%in%c("6ppm_03ml_1","6ppm_05ml_6")) %>% 
  mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10)


peaks_precission<-
  all_peaks %>% 
  filter(grepl("^Prec",label)) %>% 
  mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml"))))

peaks_precission %>% 
  summarise(avg=mean(peaksum), SD=sd(peaksum), CV=(SD/avg)*100)


peaks_linearity %>% 
  group_by(vol) %>% 
  summarise(avg=mean(peaksum), SD=sd(peaksum), CV=(SD/avg)*100) %>% 
  filter(vol<2) %>% 
  ggplot(aes(x=vol, y=CV))+
  geom_point()



cal_peaks<-all_peaks %>% 
  filter(grepl("^cal", label)) %>% 
  separate(label, into = c("conc", "vol"),remove = F, sep = "_") %>% 
  mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10) %>%
  mutate(ppm=case_when(conc=="cal6ppm"~6,
                       conc=="cal03ppm"~0.3,
                       conc=="cal1.2ppm"~1.2,
                       TRUE~NA_real_)) %>% 
  mutate(mlppm=ppm*vol) %>% 
  filter(peaksum>5) # REMOVE missidentified peaks from injection cal03ppm_01ml

ggplot(cal_peaks,aes(x=mlppm, y=peaksum))+
  geom_smooth(method = "lm", se=F, col="black")+
  geom_point(aes(col=conc, shape=vol>1))+
  stat_poly_eq(
    aes(label = ..eq.label..,), 
    formula = y~x, 
    parse = TRUE, 
    size = 5
  )


ggplot(subset(cal_peaks, vol<=1),aes(x=mlppm, y=peaksum))+
  geom_smooth(formula = y~x,method = "lm", se=T,col="black")+
  geom_point(aes(shape=conc, col= "0.1-1ml"))+
  stat_poly_eq(
    aes(label = ..eq.label..,), 
    formula = y~x, 
    parse = TRUE, 
    size = 5
  ) +
  geom_point(data=subset(cal_peaks, vol>1), aes(x=mlppm, y=peaksum, col=">1ml"))


ggplot(subset(cal_peaks, vol<=1),aes(x=mlppm, y=peaksum))+
  geom_smooth(method = "lm", se=T,col="black")+
  geom_point(aes(shape=conc, col= "0.1-1ml"))+
  geom_point(data=subset(cal_peaks, vol>1), aes(x=mlppm, y=peaksum, col=">1ml"))+
  stat_poly_eq(
    aes(label = ..eq.label..,), 
    formula = y~x, 
    parse = TRUE, 
    size = 5
  ) +
  geom_abline(intercept = 0, slope = 200, col="green")+
  facet_wrap(.~conc, scales="free")

ggplot(cal_peaks,aes(x=vol, y=peaksum,col=conc))+
  geom_smooth(method = "lm", se=F)+
  geom_point(aes(col=conc, shape=vol>1))+
  facet_wrap(.~conc, scales="free")




peaks_ini<- all_peaks %>% 
  filter(grepl("^P_05", label)) %>% 
  filter(peaksum>20)%>% 
  mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10) 

ggplot(cal6ppm, aes(x=vol, y=peaksum, col="cal"))+
  geom_point()+
  geom_point(data=peaks_linearity, aes(x=vol, y=peaksum, col="linearity"))+
  geom_point(data=peaks_ini, aes(x=vol, y=peaksum, colour = "bottle"))


head(cal6ppm)
head(peaks_linearity)
