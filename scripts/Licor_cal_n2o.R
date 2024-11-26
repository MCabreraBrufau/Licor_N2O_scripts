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
datapath_licor <- paste0(dropbox_root,"/Licor_N2O/Calibration")


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
pdf(file = paste0(datapath_licor,"/plots_integration.pdf"))  # Open PDF device

# Loop through the list of plots and print each plot
for (plot_name in names(plotspeak)) {
  print(plotspeak[[plot_name]])
}

dev.off()  # Close the PDF device



#Subset injections with properly diluted concentrations to be used for calibration (they start with "cal")
cal_peaks<-all_peaks %>% 
  filter(grepl("^cal", label)) %>% 
  separate(label, into = c("conc", "vol"),remove = F, sep = "_") %>% 
  mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10) %>%
  mutate(ppm=case_when(conc=="cal6ppm"~6,
                       conc=="cal03ppm"~0.3,
                       conc=="cal1.2ppm"~1.2,
                       TRUE~NA_real_)) %>% 
  mutate(mlppm=ppm*vol)


#Examine effect of injected volume on linearity

#All 3 dilutions together
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

#Per concentration sequence
ggplot(subset(cal_peaks, vol<=1),aes(x=vol, y=peaksum))+
  geom_smooth(method = "lm", se=T,col="black",formula = y~x)+
  geom_point(aes(shape=conc, col= "0.1-1ml"))+
  geom_point(data=subset(cal_peaks, vol>1), aes(x=vol, y=peaksum, col=">1ml"))+
  facet_wrap(.~conc, scales="free")+
  theme_bw()
  

#Linear range for injection volumes of 01-1ml (higher volumes cause underestimation of peaks due to diffussion, probably too wide peaks for integration window, baseline takes too long to recover).

#DECISSION: use 1ml injections for all samples (unless saturation is detected (very unlikely), in those cases reduce volume of injection)


#CALIBRATION FACTOR #----
cal_peaks %>% 
  filter(vol<=1) %>% #Remove peaks with more than 1ml injection
  filter(peak_id!="cal1.2ppm_01ml_1") %>% #Remove 1 outlier from cal1.2ppm injection sequence
ggplot(aes(x=mlppm, y=peaksum))+
  geom_smooth(formula = y~x,method = "lm", se=T, aes(weight = 1/mlppm,col="weighted"))+
  geom_smooth(formula = y~x,method = "lm", se=T, aes(col="unweighted"))+
  geom_point(aes(shape=conc))

#There is no apparent difference between performing a weighted vs unweighted calibration 

#DECISSION: use unweighted linear regression for calibration curve. We will use all data from the 3 dilutions (0.3ppm, 1.2ppm and 6ppm) to obtain a single calibration curve (even if there are small differences between their slopes and intercepts)

calibration_data<- cal_peaks %>% 
  filter(vol<=1) %>% #Remove peaks with more than 1ml injection
  filter(peak_id!="cal1.2ppm_01ml_1")#Remove 1 outlier from cal1.2ppm injection sequence


cal_model<-  lm(data=calibration_data, formula= peaksum~mlppm)

summary(cal_model)

#mlppm= (1/203.05)*(peaksum-5.8236)

#ppm = ((1/203.05)*(peaksum-5.8236))/ml


pdf(file=paste0(datapath_licor,"/calibration curve.pdf"))

ggplot(calibration_data, aes(x=mlppm, y=peaksum))+
  geom_smooth(method = "lm", formula= y~x, col="black")+
  geom_point(aes(col=conc))+
  scale_x_continuous(name="ml*ppm N2O", breaks = seq(0,6,by=1))+
  scale_y_continuous(name="Peak-area")+
  labs(col="Dilution")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep="~~~~")),
    formula = y~x, 
    parse = TRUE, 
    size = 5,coef.digits = 4,rr.digits = 4
  )+
  theme_bw()

dev.off()


#PRECISSION##-----
#Generally good precission for calibration data (n=5 for most) with more than 0.1ml: CV < 4%
calibration_data %>% 
  group_by(label,vol) %>% 
  summarise(avg=mean(peaksum), SD=sd(peaksum), CV=(SD/avg)*100,nobs=n()) %>% 
  filter(vol<5) %>% 
  ggplot(aes(x=vol, y=CV))+
  geom_point()+
  geom_label(aes(label = paste(label,", n=",nobs)))


#Good precission for 15 injections of 0.5ml 6ppm N2O: CV < 2%
all_peaks %>% 
  filter(grepl("^Prec",label)) %>% 
  mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10) %>% 
  group_by(label,vol) %>% 
  summarise(avg=mean(peaksum), SD=sd(peaksum), CV=(SD/avg)*100, nobs=n())




#SENSITIVITY####

#Injecting ambient lab air produces a very clear signal, even with just 0.5ml (half of what we plan to inject)

plotspeak[["aire_05ml"]]

#With our current calibration this sample has a N2O concentration of: 0.318 +- 0.01 ppm N2O

all_peaks %>% 
  filter(label=="aire_05ml") %>% 
  mutate(N2Oppm=((1/203.05)*(peaksum-5.8236))/0.5) %>% 
  summarise(mean(N2Oppm), sd(N2Oppm), cv=sd(N2Oppm)/mean(N2Oppm))



#Limit of detection: based on the quality of the calibration curve. 
#The general approach is LOD = 3 times residual standard error / slope

# Extract the slope, intercept, and their standard errors
intercept <- coef(cal_model)[1]
slope <- coef(cal_model)[2]
se_intercept <- summary(cal_model)$coefficients[1, 2]  # Standard error of the intercept
se_slope <- summary(cal_model)$coefficients[2, 2]  # Standard error of the slope
rse <- summary(cal_model)$sigma  # Residual standard error

# Calculate LOD using both intercept and slope
LOD <- (3 * sqrt(se_intercept^2 + se_slope^2 + rse^2)) / abs(slope)

# Print the LOD: 0.1137
print(LOD)



#Calibration parameters of different calibration dilutions: 
#6ppm N2O dilution:
#intercept: 16.2
#slope: 201
#LOD: 0.1562984 mlppm


#1.2ppm N2O dilution:
#intercept: 4.8
#slope: 197
#LOD: 0.06210055  mlppm


#0.3ppm N2O dilution: 
#intercept: 2.5
#slope: 232
#LOD: 0.03922452  mlppm


#Signal-to-noise ratio based on signal noise: the SD of base-corrected signal when there is no peak

noise_baseline<- all_data %>% 
  filter(label!="") %>% 
  filter(is.na(peak_id)) %>% 
  group_by(label) %>% 
  summarise(noise_sd=sd(N2Obc))

signaltonoise<- all_peaks %>% 
  merge.data.frame(noise_baseline, by = "label") %>% 
  mutate(SNR=peaksum/noise_sd)

#Gennerally, signal-to-noise ratio is >>100  
summary(signaltonoise$SNR)


#Signal-to-noise ratio for near-ambient injections (aire & 0.3ppm dilution): SNR >100, >300 for ambient air
signaltonoise %>% 
  filter(grepl("^cal03|aire", label)) %>% 
  ggplot(aes(x=label, y=SNR))+
  geom_boxplot()




#Observed vs expected 

cal_peaks %>% 
  # filter(grepl("^cal",label)) %>% 
  # mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10) %>% 
  mutate(N2Oppm=((1/203.05)*(peaksum-5.8236))/vol) %>% 
  filter(vol<=1) %>% 
  ggplot( aes(x=vol, y=N2Oppm, col=factor(ppm)))+
  geom_point()+
  geom_abline(intercept=c(0.3,1.2,6), slope=0)+
  scale_x_continuous(name = "Volume of injection (ml)", breaks = seq(0,1,by=0.1))+
  scale_y_continuous(name = "Measured concentration (ppm N2O)", breaks = c(0,0.3,1.2,6))


cal_peaks %>% 
  # filter(grepl("^cal",label)) %>% 
  # mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10) %>% 
  mutate(N2Oppm=((1/203.05)*(peaksum-5.8236))/vol) %>% 
  filter(vol<=1) %>% 
  ggplot( aes(x=ppm, y=N2Oppm, col=factor(ppm)))+
  geom_point()+
  geom_abline(intercept=0, slope=1)+
  scale_y_continuous(name = "Measured concentration (ppm N2O)", breaks = c(0,0.3,1.2,6))+
  scale_x_continuous(name = "Expected concentration (ppm N2O)", breaks = c(0,0.3,1.2,6))



calibration_data %>% 
  # filter(grepl("^cal",label)) %>% 
  # mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10) %>% 
  mutate(N2Oppm=((1/203.05)*(peaksum-5.8236))/vol) %>% 
  filter(vol<=1) %>% 
  ggplot( aes(x=vol, y=N2Oppm/ppm, col=factor(ppm)))+
  geom_boxplot(aes(group = paste(vol, ppm)))+
  geom_abline(intercept=c(1), slope=0)+
  scale_x_continuous(name = "Volume of injection (ml)", breaks = seq(0,1,by=0.1))+
  scale_y_continuous(name = "Measured/expected (ppm:ppm N2O)", breaks = seq(0.9,1.1,by=0.05))







#___________________####

#Miscelanea####


#Baseline correction does not perform well when air supply is running short (for these two injections air flow from bottle dropped a bit)

ggplot(all_peaks, aes(x=peaksum, y=peaksum_nobasec))+
  geom_smooth(formula = y~x,method = "lm", se=T,col="black")+
  geom_point()+
  stat_poly_eq(
    aes(label = ..eq.label..,), 
    formula = y~x, 
    parse = TRUE, 
    size = 5
  ) 


#baseline correction works well for volumes less than 1ml, when more than that is injected, the peak tail spreads too long, and its difficult to discriminate for baseline and peak-window

#Volumes injected should be between 0.2 and 1ml. higher volumes cause problems of diffusion and impact the baseline. The method is sensitive enough to produce accurate results CV of 2.5% with 0.5ml of ambient concentrations. 



all_peaks %>% 
  filter(grepl("^6ppm",label)) %>% 
  # filter(!peak_id%in%c("6ppm_03ml_1","6ppm_05ml_6")) %>% 
  mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10) %>% 
  filter(vol<=4) %>%   
  ggplot( aes(x=vol, y=peaksum))+
  geom_point()+
  geom_smooth(aes(col=vol<=1),method = "lm", se=F)+
  geom_smooth(aes(col="all"),method = "lm", se=F)


peaks_linearity<- 
  all_peaks %>% 
  filter(grepl("^6ppm",label)) %>% 
  # filter(!peak_id%in%c("6ppm_03ml_1","6ppm_05ml_6")) %>% 
  mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10)






all_peaks %>% 
  mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10) %>% 
  filter(peak_id!="cal1.2ppm_01ml_1")%>% 
  filter(grepl("^cal", label)) %>% 
  group_by(label,vol) %>% 
  summarise(avg=mean(peaksum), SD=sd(peaksum), CV=(SD/avg)*100) %>% 
  filter(vol<5) %>% 
  ggplot(aes(x=vol, y=CV))+
  geom_point()+
  geom_label(aes(label = label))



#Coeficient of variation is aproximately 2.5% sligthly higher for volumes <0.2ml. Even for ambient concentrations (see aire_0.5ml results)

# I would keep is it worth it to perform baseline correction? no apparent differences for 
cal_peaks<-all_peaks %>% 
  filter(grepl("^cal", label)) %>% 
  separate(label, into = c("conc", "vol"),remove = F, sep = "_") %>% 
  mutate(vol=as.numeric(gsub("ml","",str_extract(label, "[0-9]{2}ml")))/10) %>%
  mutate(ppm=case_when(conc=="cal6ppm"~6,
                       conc=="cal03ppm"~0.3,
                       conc=="cal1.2ppm"~1.2,
                       TRUE~NA_real_)) %>% 
  mutate(mlppm=ppm*vol)


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


# Compare weighted vs unweighted linear fits for calibration data with volumes of 0.1-1ml 
cal_data<- subset(cal_peaks, vol<=1) %>% filter(peak_id!="cal1.2ppm_01ml_1")

ggplot(cal_data,aes(x=mlppm, y=peaksum))+
  geom_smooth(formula = y~x,method = "lm", se=T, aes(weight = 1/mlppm,col="weighted"))+
  geom_smooth(formula = y~x,method = "lm", se=T, aes(col="unweighted"))+
  geom_point(aes(shape=conc))+
  stat_poly_eq(
    aes(label = ..eq.label..), 
    formula = y~x, 
    parse = TRUE, 
    size = 5
  )

#weighted model
# Fit weighted linear regression model
weighted_model <- lm(cal_data$peaksum ~ cal_data$mlppm, weights = 1/cal_data$mlppm)
# Fit unweighted linear regression model
unweighted_model <- lm(cal_data$peaksum ~ cal_data$mlppm)

# Summary of the weighted model
summary(weighted_model)


# Summary of the unweighted model
summary(unweighted_model)


anova(weighted_model,unweighted_model)


plot(unweighted_model)

plot(weighted_model)






ggplot(cal_data,aes(x=vol, y=peaksum))+
  geom_smooth(method = "lm", se=T,col="black",formula = y~x)+
  geom_point(aes(shape=conc, col= "0.1-1ml"))+
  geom_point(data=subset(cal_peaks, vol>1), aes(x=vol, y=peaksum, col=">1ml"))+
  stat_poly_eq(
    aes(label = ..eq.label..,), 
    formula = y~x, 
    parse = TRUE, 
    size = 5
  ) +
  geom_abline(intercept = 5.82, slope = 203, col="green")+
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
