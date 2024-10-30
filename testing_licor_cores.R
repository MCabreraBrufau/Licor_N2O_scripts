
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Oct 2024"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
##Working script for decission making into how to process the peaks obtain with the Valencia Li-COR 
#We will use one of the last rawfiles created as an example to build the processing pipeline

#Steps: 
#1 BASELINE correction of "spectra": using function baseline.corr from ptw package (tunning of parameters after visual inspection)
#2 Detect peaks within remarked spectra (not always same number of peaks with a single remark, 2-5). Function findpeaks from pracma package
#3 Define width of peaks (generally well constrained, but no unique width):
    #Option 1: fixed width (as data is already baseline-corrected, adding a bit of baseline should not be an issue)
    #Option 2: define peak width based on derivative (of smoothed data?)
#4 Integrate area of peaks



#Generally, much more clear signal for CH4 than for CO2 (much bigger difference between baseline and peaks). Peaks of CH4 and CO2 are simultaneous (make sure!!), so we will define the peaks based on the CH4 signal alone, then pull the CO2 values corresponding to CH4 peak. For the baseline we cannot use the same approach, as CO2 baseline is much more noisy, with some Humps in between and throughout the peaks. 


rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- packages ----
library(tidyverse)
library(readxl)
library(lubridate)
library(ptw)
library(pracma)
# library(zoo)
# library(ggplot2)
# library(grid)
# library(egg)
# library(goFlux)
# require(dplyr)
# require(purrr)
# require(msm)
# require(data.table)
# require(tools)
# 
# require(pbapply)

#Browse restore4cs-scripts functions:
repo_root <- "C:/Users/Miguel/Documents/GitHub/restore4cs-scripts"
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
files.sources

#Load functions needed: 
source(files.sources[grep("read_Licor",files.sources)])



# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
datapath_licor <- paste0(dropbox_root,"/GHG/RAW data/RAW Data Licor-7810")


# load incubation map (run get_exact_incubation_times_Licor.R to update it)
map_incubations<- read.csv(file = paste0(datapath_licor, "/map_incubations.csv"))


# ---- Import ----
# We will eventually need a complete and exclusive list of all cores injection remarks
#For the moment, use an example (day with only core incubations)

example_file<- paste0(dropbox_root, "/GHG/RAW data/RAW Data Licor-7810/S3-CA/TG10-01275-2024-04-26T060000.data.txt")


#Import using camille's function
a<- read_Licor(example_file)


a %>% filter(label!="") %>% 
  group_by(label) %>% 
  summarise(nobs=n()) %>% 
  ggplot(aes(x=1,y=nobs))+
  geom_boxplot()+
  geom_point()

#Typically, 3-4 peaks during a ~60s remark duration (40-120s)
#max intensity of peaks should be reliably among the top 0.8 percentile

a %>% 
  mutate(datetime=ymd_hms(paste(date,UTCtime))) %>% 
  filter(label%in%unique(a$label)[2:12])%>% 
  ggplot(aes(x=datetime))+
  geom_point(aes(y=CH4*1000,col="ch4"))+
  geom_line(aes(y=CH4*1000,col="ch4"))+
  geom_point(aes(y=CO2*4,col="co2"))+
  geom_line(aes(y=CO2*4,col="co2"))


#We can use the findpeaks within the pracma package to find peaks and filter setting the minimum peak height as greater than the 0.8 percentile.
#We will have to manually explore peak identification performance within our dataset


#---isolated examples----
library(pracma)

# Create a sample time series
time_series_data <- a %>%   filter(label==unique(a$label)[3]) %>% pull(CH4)

# Find peaks
peaks <- findpeaks(time_series_data,minpeakheight = quantile(time_series_data, 0.85))

# Print peaks
print(peaks)

# Visualize the data and detected peaks
plot(time_series_data, type = 'b', main = "Peak Detection Example", xlab = "Index", ylab = "Value")
points(peaks[,2], peaks[,1], col = "red", pch = 19)  # Mark peaks in red





#Other options: 

# The zoo package provides methods for dealing with irregular time series data. While it doesn’t directly detect peaks, you can apply smoothing functions and then use which.max() or custom logic to identify peaks.
library(zoo)

data <- c(1, 3, 7, 1, 2, 5, 8, 6, 4, 9, 3, 2)
smoothed <- rollmean(data, k = 3, fill = NA, align = "center")
peaks <- which(diff(sign(diff(smoothed))) == -2) + 1  # Local maxima



# 5. Signal Package
# The signal package provides tools for signal processing, including peak detection through convolution or filtering techniques.
# Example 
library(signal)

data <- c(1, 3, 7, 1, 2, 5, 8, 6, 4, 9, 3, 2)
smoothed <- sgolayfilt(data, p = 3, n = 11)  # Savitzky-Golay filter
peaks <- findpeaks(smoothed)






#1. Baseline correction: 

#function ptw::baseline.corr(testb[, 2])
library(ptw)
# asysm: Trend estimation with asymmetric least squares 
baseline.corr(y = )# this function uses the asysm function to provide a baseline-corrected chromatogram

asysm(y,#data to correct, either as vector or as matrix with each chromatogram as a row
      lambda = 1e7, #smoothing parameter
      p=0.001, #asymetry parameter
      eps=1e-8, #numerical precission for convergence
      maxit = 25)#max number of iterations (warns when no convergence reached)


#Default seems to work relatively well for CH4 when RT windows are supplied
a%>%
  filter(label==unique(a$label)[14])%>%
  ggplot(aes(x=unixtime))+
  geom_line(aes(y=CO2-400.027, col="OG-2"))+
  geom_line(aes(y=baseline.corr(CO2,lambda=1e2, p=0.01), col="bs-corr"))


a%>%
  filter(label==unique(a$label)[14])%>%
  ggplot(aes(x=unixtime))+
  geom_line(aes(y=CH4-2.027, col="OG-2"))+
  geom_line(aes(y=baseline.corr(CH4,lambda=1e2, p=0.01), col="bs-corr"))



a%>%
  filter(label==unique(a$label)[14])%>%
  ggplot(aes(x=unixtime))+
  geom_line(aes(y=c(NA,diff(CO2)), col="OG-2"))

#_________####
# ---- Example pipeline 1----

# -----1. Base-correct----

#Base-correct each "spectra" with a label that is not "" 

#ISSUE: if peak is very near the end of label, baseline correction does not work properly.
#TOFIX: for each label, add 5 seconds before and after to the data 

#ISSUE: for CO2 with very high concentrations, lower-than-baseline end of peaks, causes problems for baseline correction.


a_basec<- a %>% 
  filter(label!="") %>% 
  group_by(label) %>% 
  mutate(CH4bc=baseline.corr(CH4,lambda=1e2, p=0.01),
         CO2bc=baseline.corr(CO2,lambda=1e3, p=0.1))



a_basec%>%
  filter(label==unique(a_basec$label)[15])%>%
  ggplot(aes(x=unixtime))+
  geom_line(aes(y=CH4-2, col="OG-2"))+
  geom_line(aes(y=CH4bc, col="bs-corr"))


a_basec%>%
  filter(label==unique(a_basec$label)[4])%>%
  ggplot(aes(x=unixtime))+
  geom_line(aes(y=CO2-400, col="OG-2"))+
  geom_line(aes(y=CO2bc, col="bs-corr"))


# -----2. Peak-max detection ----


# Find peaks in base-corrected data
#Criteria that defines a peak-maximum:
    # at least 1 increase before and 1 decrease afterwards
    # minimum peak height to be detected is > percentil 0.9

#Flag each local maxima with a unique code (label_1, label_2,....)

a_basec_peakid <- a_basec %>%
  group_by(label) %>%
  mutate(is_localmaxch4 = ifelse(row_number() %in% findpeaks(CH4bc, minpeakheight = quantile(CH4bc, 0.90), nups=1, ndowns=1,minpeakdistance = 5)[, 2], TRUE, FALSE)) %>%
  mutate(unique_code = ifelse(is_localmaxch4, paste0(label,"_",cumsum(is_localmaxch4)), NA)) %>%  #Add unique code for local maxima 
  ungroup()




a_basec_peakid %>% 
  filter(label%in%unique(.$label)[2]) %>% 
  ggplot(aes(x=unixtime, y=CH4bc))+
  geom_point(aes(col=is_localmaxch4))+
  geom_line()

# -----3. Define peak width ----


#OptionA: fixed width. 
#Consider peak as max height + 4 leading and 4 trailing points. (i.e. peak width == 9points), 


a_basec_peakwindow <- a_basec_peakid %>%
  mutate(unique_code = map_chr(row_number(), function(idx) {
    seq_start <- max(1, idx - 4)
    seq_end <- min(n(), idx + 4)
    
    # Check for unique codes in the surrounding window
    surrounding_codes <- unique_code[seq(seq_start, seq_end)]
    
    # Return the unique code if it's available, otherwise return NA
    if (any(!is.na(surrounding_codes))) {
      return(first(na.omit(surrounding_codes)))  # Use the first valid unique code found
    } else {
      return(NA)
    }
  }))


a_basec_peakwindow %>% 
  filter(label%in%unique(.$label)[1:11]) %>% 
  ggplot(aes(x=unixtime, y=CH4bc))+
  geom_point(aes(col=!is.na(unique_code)))+
  geom_line()







# -----4. Integrate peaks ----


a_integrated<- a_basec_peakwindow %>% 
  filter(!is.na(unique_code)) %>% 
  group_by(date, label, unique_code) %>% 
  summarise(peaksum=sum(CH4bc), peakmax=max(CH4bc)) %>% 
  ungroup()


a_integrated_qual<- a_integrated %>% 
  group_by(label) %>% 
  mutate(peaksum_relativetomean=peaksum/mean(peaksum),
         peakmax_relativetomean=peakmax/mean(peakmax))


ggplot(a_integrated_qual, aes(x=label))+
  geom_point(aes( y=peaksum_relativetomean,col="asum"))+
  geom_point(aes( y=peakmax_relativetomean,col="amax"))
  


ggplot(a_integrated, aes(x=peaksum, y=peakmax, col=label))+
  geom_point()




#_________####

# ---- Example in loop 2----

#process each label individually then join the results
# -----0. Expand data injection----

#Keep all injection sequences identified by label and mark their start and end (plus 5 seconds)
startend_inj<- a %>% 
  filter(label!="") %>% 
  group_by(label) %>%
  summarise(
    start_inj = first(unixtime)-5,
    end_inj = last(unixtime)+5,
    .groups = 'drop'  # To ungroup after summarization
  )


#There is one injection sequence that has a weird peak at the end and does not correspond to any injection.
#We will crop this one injection
startend_inj[startend_inj$label=="S3-CA-A2-C2-TFM",]$end_inj<- 1714123810


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
    mutate(CH4bc=baseline.corr(CH4,lambda=1e2, p=0.01))
    
  
  # -----2. Peak-max detection ----
  #Find local maxima in sequence and add max_id (label_1,label_2,...) : 
    #Criteria for local maximum:
          # at least 1 increase before and 1 decrease afterwards
          # minimum peak height to be detected is > percentil 0.9
          # at leas 5 points between localmaxima
  
  dat_peakid <- dat_bc %>%
    mutate(is_localmaxch4 = ifelse(row_number() %in% findpeaks(CH4bc, 
                                                               minpeakheight = quantile(CH4bc, 0.90), 
                                                               nups=1, ndowns=1,
                                                               minpeakdistance = 5)[, 2], TRUE, FALSE)) %>%
    mutate(peak_id = ifelse(is_localmaxch4, paste0(label,"_",cumsum(is_localmaxch4)), NA)) %>%  #Add unique code for local maxima 
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
    summarise(peaksum=sum(CH4bc), peakmax=max(CH4bc), unixtime_ofmax=unixtime[CH4bc==max(CH4bc)]) %>% 
    ungroup()
  
  # -----5. Plot results ----
  #Create a plot for inspection:
  p<- ggplot(dat_peakwindow, aes(x=unixtime, y=CH4bc))+
    geom_point(aes(col=peak_id))+
    geom_line()+
    geom_point(data = dat_integrated, aes(x=unixtime_ofmax, y=peaksum, col="integral"))
  
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
 


