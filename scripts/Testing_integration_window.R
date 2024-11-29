#Test effect of integration window on precission and differences by operator



# ---- Directories ----

#Root
dropbox_root <- "C:/Users/Miguel/Dropbox" # You have to make sure this is pointing to the write folder on your local machine

#Data folders
folder_raw <- paste0(dropbox_root,"/Licor_N2O/Rawdata") #contains unedited files downloaded from licor

folder_mapinjections<- paste0(dropbox_root,"/Licor_N2O/Map_injections") #Contains csvs with startstop times of injections and their corresponding labels, corrections should be made manually when needed (editting the csvs and re-saving with "corrected_" prefix)

folder_plots<-  paste0(dropbox_root,"/Licor_N2O/Testing integration window") #One pdf per dayofinjections (auto-name from rawfile name), plots of each injection sequence (baseline correction & integration)

folder_results<- paste0(dropbox_root,"/Licor_N2O/Testing integration window")#One csv per dayofinjections will be created (auto-name from rawfile name), with individual peak parameters (label, peak_id, peaksum, peakmax, unixtime_ofmax, raw_peaksum, dayofanalysis, SNR)

folder_results_done<- paste0(dropbox_root,"/Licor_N2O/Results_ppm")

# ---- packages & functions ----
library(tidyverse)
library(readxl)
library(lubridate)
library(ptw)
library(pracma)
library(stringr)
library(ggpmisc)
library(broom)

#Load repository functions
repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}



####1. Inspect lags and decide windows to test#####
#Collect all peaks timings:

peakfiles<- list.files(path = folder_results_done, pattern = "^integrated_")

peaksdata<- data.frame()
for(i in peakfiles){
  
  a<- read.csv(file = paste0(folder_results_done,"/", i)) 
  peaksdata<- rbind(peaksdata, a)
}

#Get lag data (time between injections of same label)
lagdata<- peaksdata %>%
  arrange(unixtime_ofmax) %>%
  group_by(label,dayofanalysis) %>%
  mutate(
    lagpeak = unixtime_ofmax - lag(unixtime_ofmax),  # Calculate lag with previous row
    first_instance = peak_id,  # Value of the first instance
    second_instance = lag(peak_id)     # Value of the second instance
  ) %>%
  filter(!is.na(lagpeak)) %>%  # Remove NA rows where lag cannot be calculated
  ungroup() %>% 
  select(dayofanalysis, label, peak_id, first_instance, second_instance, lagpeak)

#Summarise the distances between peaks of the same remark for calibration and samples
lagdata %>% 
  filter(grepl("cal|S4", label)) %>% 
  arrange(lagpeak) %>% 
  summarise(avg_lag=mean(lagpeak),
            sd_lag=sd(lagpeak),
            min_lag=min(lagpeak),
            max_lag=max(lagpeak))

lagdata %>% 
  filter(grepl("cal|S4", label)) %>% 
ggplot( aes(x=dayofanalysis, y=lagpeak))+
  geom_boxplot()+
  geom_point()+
  scale_y_continuous(limits = c(0,100), breaks = seq(0,20,by=2))

#We could use wider integration windows for Restore4Cs samples, but the calibration injections are the limiting factor: 
#The widest peak integration that we can do is 17 seconds before starting to overlap between peaks. 
#exploring the integration plots, we can delay slightly the start of the integration window to 3 seconds before maxpeak

#Test widths of window with the following: 
#secondsbefore=3
#secondsafter= c(7, 10, 11, 12, 13, 14)

#We will have some minor overlap for a couple of samples (peaklag< (3+1+14)):
lagdata %>% 
  filter(grepl("cal|S4", label)) %>% 
  filter(lagpeak<=18)



#We have to test:
  #effect on precission (CV of same label)
  #effect on calibration factor and recovery 

#To do this we will process with != windows data from the calibration curve and from the precission tests of 20241125

#Modify integration loop to output only label, peak_id, unixtimeofmax, int7,int8,int9, int10, 


#Get rawfiles
rawfiles<- list.files(path = folder_raw, pattern = ".data")

#Get corrected maps of injections
mapscorrect<- list.files(path = folder_mapinjections, pattern = "^corrected_map_injection_")


#Select code of rawfiles with corresponding mapscorrect but without integratedfiles
rawtointegrate<- gsub(".data","",rawfiles[
  gsub(".data","",rawfiles)%in%gsub(".csv","",gsub("corrected_map_injection_","",mapscorrect)) #Match raw with maps
])

#Choose only days of calibration and of testing Ilenia vs Miguel standards 6ppm & aire
rawtointegrate<- rawtointegrate[grepl("2024-10-31|2024-11-25", rawtointegrate)]



###2. Integration loop####
seconds_aftertotest<-c(7,8,9,10,11,12,13,14) #max is 14s after peak

for (s in seconds_aftertotest){

#Loop over raw files
for (i in rawtointegrate){
  
  message(paste("Integrating peaks from",i, "with ",s,"-second window" ))
  
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
      #Consider peakwindow as max height + 4 leading and 7 trailing points. (i.e. peak width == 12points), 
      
      inj_data <- inj_data %>%
        mutate(peak_id = map_chr(row_number(), function(idx) {
          #For each row, search for a non-na peak_id, look up to 7 seconds before and 4 seconds after the row i. Then assing the value of peak_id to the row i.
          #This results in the spread of the value of "peak_id" of the local maximum to secondsbefore_max seconds before and to secondsafter_max seconds after each identified maximum. 
          secondsbefore_max<- 3
          secondsafter_max<- s
          
          # Check for peak_id in the window:
          surrounding_codes <- peak_id[seq(max(1, idx - secondsafter_max), min(n(), idx + secondsbefore_max))]  
          
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
        summarise(intwind=s, 
                  peaksum=sum(n2o_bc),
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
        ggtitle(paste0(dayofanalysis,", injection: ",inj, " and ",s, "-s window"))+
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
  write.csv(A,file = paste0(folder_results,"/", "integrated_injections_",s,"seconds_", i, ".csv"),row.names = F)
  
  #Save plots of integrations: use i for naming convention of pdf
  print(paste0("Plotting integrations of day: ", i,"with",s,"seconds"))
  #plot every injection sequence and their integrals: 
  pdf(file = paste0(folder_plots,"/Integrations_",s,"seconds_",i,".pdf"))  # Open PDF device
  
  # Loop through the list of plots and print each plot
  for (plot_name in names(plotspeak)) {
    print(plotspeak[[plot_name]])
  }
  
  dev.off()  # Close the PDF device
  
}

}


rm(s,unixstart,unixend,seconds_aftertotest,sd_peaksum,rawtointegrate,raw_data,plot_name,i,inj,a,A,b,B,integrated, inj_data, mapinj,p,plotspeak,avg_peaksum,baseline_sd,dayofanalysis)
#Import and compare 

####3. Import integrations####
peakfiles<- list.files(path = folder_results, pattern = "^integrated_")

peaksdata<- data.frame()
for(i in peakfiles){
  
  a<- read.csv(file = paste0(folder_results,"/", i)) 
  peaksdata<- rbind(peaksdata, a)
}

rm(a)


#####3.1. Check precission####
peaksdataf<- peaksdata %>% 
  select(dayofanalysis, label, peak_id, intwind, peaksum) %>% 
  filter(grepl("^cal|^M|^I", label)) %>% #get only calibration data and tests of operators
  filter(!grepl("_0.1$", label))#Remove 0.1ml injections 

cvs<- peaksdataf %>% 
  mutate(sampletype=case_when(grepl("^cal", label)~"cal",
                              grepl("^M|^I", label)~"samples")) %>% 
  group_by(label, intwind,sampletype) %>% 
  summarise(avg_sum=mean(peaksum),
            sd_sum=sd(peaksum),
            cv_sum=100*(sd_sum/avg_sum)) %>% 
  separate(label, into = c("gas", "vol"), sep = "_", remove = F)


cvs %>% 
  ggplot(aes(x=intwind, y=cv_sum,col=sampletype))+
  geom_point()+
  geom_line(aes(group = label))+
  geom_boxplot(aes(group = intwind))+
  facet_wrap(.~sampletype)+
  ggtitle("Precission as a function of integration window")
#Precission improves with wider integration windows for samples (glass syringe), not for calibration (plastic syringe)





#####3.2. Fit calcurves#####

cals<- cvs %>% 
  filter(sampletype=="cal") %>% 
  ungroup() %>% 
  mutate(gas=case_when(gas=="cal03ppm"~"cal0.3ppm",TRUE~gas)) %>% 
  mutate(ppm=parse_number(gas), vol=parse_number(vol),mlppm=ppm*vol)



#All 3 dilutions together as a function of integration window
ggplot(subset(cals, vol<=1),aes(x=mlppm, y=avg_sum))+
  geom_smooth(formula = y~x,method = "lm", se=T,col="black")+
  geom_point(aes(shape=gas, col= "0.1-1ml"))+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep="~~~~")),
    formula = y~x, 
    parse = TRUE, 
    size = 5,coef.digits = 4,rr.digits = 4
  )+
  geom_point(data=subset(cals, vol>1), aes(x=mlppm, y=avg_sum, col=">1ml"))+
  facet_wrap(.~intwind)
#Linear range remains, slight difference in cal factors as a function of integration window 



#Obtain 1 calibration curve per integration window

calfactors<- cals %>% 
  filter(vol<=1) %>% # Remove injections with large volumes.
  group_by(intwind) %>%  # Group by the integration window
  do({
    model <- lm(avg_sum ~ mlppm, data = .)  # Fit a linear model within each group
    tidy_model <- tidy(model)  # Extract model coefficients (including intercept)
    glance_model <- glance(model)  # Extract model summary statistics (including R-squared)
    tidy_model %>%
      bind_cols(glance_model %>% select(r.squared))  # Bind R-squared with the coefficients
  }) %>%
  select(intwind, term, estimate,r.squared) %>% 
  mutate(term=case_when(term=="(Intercept)"~"intercept",
                        term=="mlppm"~"slope")) %>% 
  pivot_wider(names_from = term,values_from = estimate) # Select relevant columns 



#CALIBRATION using glass-syringe injections

#Plot:
cvs %>% 
  filter(sampletype=="samples"&grepl("cal",label)) %>% 
  mutate(ppm=6, vol=parse_number(vol), mlppm=ppm*vol) %>% 
  ggplot(aes(x=mlppm, y=avg_sum))+
  geom_smooth(formula = y~x,method = "lm", se=T,col="black")+
  geom_point(aes(shape=gas, col= "0.1-1ml"))+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep="~~~~")),
    formula = y~x, 
    parse = TRUE, 
    size = 5,coef.digits = 4,rr.digits = 4
  )+
  # geom_point(data=subset(cals, vol>1), aes(x=mlppm, y=avg_sum, col=">1ml"))+
  facet_wrap(.~intwind)

#Calfactors:
calfactorglass<- cvs %>% 
  filter(sampletype=="samples"&grepl("cal",label)) %>% 
  mutate(ppm=6, vol=parse_number(vol), mlppm=ppm*vol) %>% 
  group_by(intwind) %>%  # Group by the integration window
  do({
    model <- lm(avg_sum ~ mlppm, data = .)  # Fit a linear model within each group
    tidy_model <- tidy(model)  # Extract model coefficients (including intercept)
    glance_model <- glance(model)  # Extract model summary statistics (including R-squared)
    tidy_model %>%
      bind_cols(glance_model %>% select(r.squared))  # Bind R-squared with the coefficients
  }) %>%
  select(intwind, term, estimate,r.squared) %>% 
  mutate(term=case_when(term=="(Intercept)"~"intercept",
                        term=="mlppm"~"slope")) %>% 
  pivot_wider(names_from = term,values_from = estimate) # Select relevant columns 


#####3.3. Check recovery####


#get baseline of day of tests 20241125
base_ppm<- read.csv(paste0(folder_results,"/baselines_TG20-01377-2024-11-25T060000.csv")) %>% 
  filter(label=="lab-1_baseline") %>% 
  mutate(base_ppm=base_avg/1000) %>% 
  pull(base_ppm)



#######__a. plastic calibrations#####
#Calculate recovery of all samples using plastic calibrations
recovery_plastical<- peaksdataf %>% 
  filter(!grepl("cal",label)) %>% 
  merge.data.frame(calfactors, by="intwind") %>% 
  separate(label, into = c("gas", "vol"), sep = "_", remove = F) %>%
  mutate(vol=parse_number(vol)) %>% 
  mutate(n2o_ppm= ((1/slope)*(peaksum-intercept))/vol) %>% 
  select(label, intwind, n2o_ppm) %>% 
group_by(label, intwind) %>% 
  mutate(operator=case_when(grepl("I",label)~"Ilenia",
                            grepl("M",label)~"Miguel"),
         gas=case_when(grepl("air", label)~"Lab Air",
                       grepl("6ppm", label)~"6ppm")) %>% 
  mutate(reference_n2o=case_when(grepl("air", label)~base_ppm,
                                 grepl("6ppm",label)~6),
         recovery=100*n2o_ppm/reference_n2o)

#Plot recovery per int window and for the 2 types of gas
ggplot(recovery_plastical, aes(x=intwind, y=recovery, col=operator))+
  geom_point()+
  # geom_line(aes(group=))
  geom_violin(draw_quantiles = c(0.5),aes(group=paste(intwind,operator)))+
  facet_wrap(.~gas)
#Very little increase in recovery with increasing integration window 


recovery_plastical %>% 
  group_by(intwind, operator,label,gas) %>% 
  summarise(avg_recovery=mean(recovery)) %>% 
  ggplot(aes(x=intwind, y=avg_recovery, col=operator))+
  geom_point()+
  geom_boxplot(aes(group=paste(intwind,operator)))+
  geom_line(aes(group = label))+
  geom_boxplot(aes(group=paste(intwind,operator)))+
  facet_wrap(.~gas)+
  scale_y_continuous(breaks=seq(92,100,by=0.5))+
  ggtitle("Average recovery as a function of integration window")


#Change in recovery (with respect to 7s window)

recovery_plastical %>% 
  group_by(intwind, operator,label,gas) %>% 
  summarise(avg_recovery=mean(recovery)) %>% 
  ungroup() %>% 
  group_by(label) %>% 
  mutate(recovery_7s=case_when(intwind==7~avg_recovery, intwind!=7~NA_real_)) %>% 
  mutate(recovery_7s=mean(recovery_7s,na.rm = T)) %>% 
  mutate(increase_recovery=avg_recovery-recovery_7s) %>% 
  ggplot(aes(x=intwind, y=increase_recovery, col=operator))+
  geom_point()+
  geom_boxplot(aes(group=paste(intwind,operator)))+
  geom_line(aes(group = label))+
  geom_boxplot(aes(group=paste(intwind,operator)))+
  facet_wrap(.~gas)+
  ggtitle("Increase in recovery (with respect to 7swindow)")



#####DECISSION 13s window####

#Increasing the integration window improves precission for normal samples (13 seconds reduces the CV of injections from an average of 1.6% to 1% )
#Increasing the integration window also improves recovery (13 seconds improves 1% recovery of air samples for ilenia, from 96.2 to 97.2)



cvs %>% 
  group_by(label) %>% 
  mutate(cv_7s=case_when(intwind==7~cv_sum, intwind!=7~NA_real_)) %>% 
  mutate(cv_7s=mean(cv_7s,na.rm = T)) %>% 
  mutate(increase_CV=cv_sum-cv_7s) %>% 
  ggplot(aes(x=intwind, y=increase_CV,col=sampletype))+
  geom_point()+
  geom_line(aes(group = label))+
  geom_boxplot(aes(group = intwind))+
  facet_wrap(.~sampletype)+
  ggtitle("Precission as a function of integration window")








#Precission improves with wider integration windows for samples (glass syringe), not for calibration (plastic syringe)








########__b. glass cal (DONT USE)#####
#Calculate recovery of all samples using glass calibrations
recovery_glass<- peaksdataf %>% 
  filter(!grepl("cal",label)) %>% 
  merge.data.frame(calfactorglass, by="intwind") %>% 
  separate(label, into = c("gas", "vol"), sep = "_", remove = F) %>%
  mutate(vol=parse_number(vol)) %>% 
  mutate(n2o_ppm= ((1/slope)*(peaksum-intercept))/vol) %>% 
  select(label, intwind, n2o_ppm) %>% 
  group_by(label, intwind) %>% 
  mutate(operator=case_when(grepl("I",label)~"Ilenia",
                            grepl("M",label)~"Miguel"),
         gas=case_when(grepl("air", label)~"Lab Air",
                       grepl("6ppm", label)~"6ppm")) %>% 
  mutate(reference_n2o=case_when(grepl("air", label)~base_ppm,
                                 grepl("6ppm",label)~6),
         recovery=100*n2o_ppm/reference_n2o)

#Plot recovery per int window and for the 2 types of gas
ggplot(recovery_glass, aes(x=intwind, y=recovery, col=operator))+
  geom_point()+
  # geom_line(aes(group=))
  geom_violin(draw_quantiles = c(0.5),aes(group=paste(intwind,operator)))+
  facet_wrap(.~gas)





