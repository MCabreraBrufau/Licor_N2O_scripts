#Quality check N2O


# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Nov 2024"
# https://github.com/MCabreraBrufau/Licor_N2O_scripts
# ---


#Description: this script uses ppm_peak files and baselines files produced in the Raw_to_peaks_LicorN2O.R and Peaks_to_ppm_LicorN2O.R scripts to compare the concentrations retrieved by the open loop method with nominal standard gas concentration (6ppm) and with ambient lab-air concentration (obtained in baseline measurements)



# ---- Directories ----

#Root
dropbox_root <- "C:/Users/Miguel/Dropbox" # You have to make sure this is pointing to the write folder on your local machine

#Data folders

folder_results<- paste0(dropbox_root,"/Licor_N2O/Results_ppm")#Contains baseline files and ppm_peak files.


#Collect all baseline lab measurements

baselinefiles<- list.files(path = folder_results, pattern = "^baselines_")
ppmpeakfiles<- list.files(path = folder_results, pattern = "^ppm_peaks_")
ppmsamplefiles<- list.files(path = folder_results, pattern = "^ppm_samples_")

baselinedata<- data.frame()
for(i in baselinefiles){
  
  a<- read.csv(file = paste0(folder_results,"/", i)) %>% 
    filter(grepl("^lab", label))

  baselinedata<- rbind(baselinedata, a)
  }

#Collect all injection measurements of aire and 6ppm 
injdata<- data.frame()
for(i in ppmpeakfiles){
  
  a<- read.csv(file = paste0(folder_results,"/", i))
    # filter(grepl("^aire|^6ppm", sample))
  
  injdata<- rbind(injdata, a)
}

sampledata<- data.frame()
for(i in ppmsamplefiles){
  
  a<- read.csv(file = paste0(folder_results,"/", i))
  # filter(grepl("^aire|^6ppm", sample))
  
  sampledata<- rbind(sampledata, a)
}


sampledata$dayofanalysis<- as.POSIXct(sampledata$dayofanalysis,format = "%Y-%m-%d")




#CV of R4Cs samples from exetainers with 3 and 4 injections :
sampledata %>% 
  filter(grepl("^S", sample)) %>% 
ggplot(aes(x=as.character(nobs), y=n2o_ppm_cv*100, group = as.character(nobs)))+
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  scale_x_discrete(name="Injections per sample")+
  scale_y_continuous(name = "CV %")+
  ggtitle("Precission of R4Cs Samples (from exetainers)")+
  theme_bw()

#Summary precission CV
sampledata %>% 
  filter(grepl("^S", sample)) %>% 
  group_by(nobs) %>% 
  summarise(n=n(),
            avg_cv=mean(n2o_ppm_cv*100),
            sd_cv=sd(n2o_ppm_cv*100),
            CI_lower = mean(n2o_ppm_cv*100) - qnorm(0.975) * (sd(n2o_ppm_cv*100) / sqrt(n())),
            CI_upper = mean(n2o_ppm_cv*100) + qnorm(0.975) * (sd(n2o_ppm_cv*100) / sqrt(n())))

#Summary precission SD ppm
sampledata %>% 
  filter(grepl("^S", sample)) %>% 
  group_by(nobs) %>% 
  summarise(n=n(),
            avg_sd=mean(n2o_ppm_sd*1000),
            sd_sd=sd(n2o_ppm_sd*1000),
            CI_lower = mean(n2o_ppm_sd*1000) - qnorm(0.975) * (sd(n2o_ppm_sd*1000) / sqrt(n())),
            CI_upper = mean(n2o_ppm_sd*1000) + qnorm(0.975) * (sd(n2o_ppm_sd*1000) / sqrt(n())))



#Accuracy (%retrieval)

#Sampling ambient lab air from exetainers
sampledata %>% 
  filter(dayofanalysis==as.POSIXct("2024-11-25",format = "%Y-%m-%d")) %>% 
  filter(!grepl("bot|cal|Dif",sample)) %>% 
  separate(sample, into = c("operator", "origin","replicate"), remove = F) %>% 
  mutate(retrieval=case_when(origin=="6ppm"~(n2o_ppm_avg/6)*100,
                             origin=="air"~(n2o_ppm_avg/0.343)*100)) %>%
  ggplot(aes(x=operator,y=retrieval, col=operator))+
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  geom_hline(yintercept = 100)+
  facet_wrap(.~origin)+
  theme_bw()

#Sampling ambient lab air from exetainers
sampledata %>% 
  filter(dayofanalysis==as.POSIXct("2024-11-25",format = "%Y-%m-%d")) %>% 
  filter(!grepl("bot|cal|Dif",sample)) %>% 
  separate(sample, into = c("operator", "origin","replicate"), remove = F) %>% 
  mutate(retrieval=case_when(origin=="6ppm"~(n2o_ppm_avg/6)*100,
                             origin=="air"~(n2o_ppm_avg/0.343)*100)) %>%
  ggplot(aes(x=origin,y=retrieval))+
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  # geom_hline(yintercept = 100)+
  facet_wrap(.~origin, scales = "free")+
  theme_bw()



# Summary %retrieval
sampledata %>% 
  filter(dayofanalysis==as.POSIXct("2024-11-25",format = "%Y-%m-%d")) %>% 
  filter(!grepl("bot|cal|Dif",sample)) %>% 
  separate(sample, into = c("operator", "origin","replicate"), remove = F) %>% 
  mutate(retrieval=case_when(origin=="6ppm"~(n2o_ppm_avg/6)*100,
                             origin=="air"~(n2o_ppm_avg/0.343)*100)) %>% #Ave concentration baseline-1 (when exetainers were filled)
  group_by(operator, origin) %>%
  summarise(n=n(),
            avg_retrieval=mean(retrieval),
            sd_retrieval=sd(retrieval),
            CI_lower = mean(retrieval) - qnorm(0.975) * (sd(retrieval) / sqrt(n())),
            CI_upper = mean(retrieval) + qnorm(0.975) * (sd(retrieval) / sqrt(n())))



#Sampling ambient lab air from exetainers
injdata %>% 
  filter(dayofanalysis==as.POSIXct("2024-11-25",format = "%Y-%m-%d")) %>% 
  filter(!grepl("bot|cal|Dif",sample)) %>% 
  separate(sample, into = c("operator", "origin","replicate"), remove = F) %>% 
  mutate(retrieval=case_when(origin=="6ppm"~(n2o_ppm/6)*100,
                             origin=="air"~(n2o_ppm/0.343)*100)) %>%
  # filter(origin=="6ppm") %>% 
  ggplot(aes(x=operator,y=retrieval,col=operator))+
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  geom_hline(yintercept = 100)+
  facet_grid(.~origin)+
  theme_bw()+
  scale_y_continuous(name="% Recovery of method ")


injdata %>% 
  filter(dayofanalysis==as.POSIXct("2024-11-25",format = "%Y-%m-%d")) %>% 
  filter(!grepl("bot|cal|Dif",sample)) %>% 
  separate(sample, into = c("operator", "origin","replicate"), remove = F) %>% 
  mutate(retrieval=case_when(origin=="6ppm"~(n2o_ppm/6)*100,
                             origin=="air"~(n2o_ppm/0.343)*100)) %>%
  group_by(operator, origin) %>%
  summarise(n=n(),
            avg_retrieval=mean(retrieval),
            sd_retrieval=sd(retrieval),
            CI_lower = mean(retrieval) - qnorm(0.975) * (sd(retrieval) / sqrt(n())),
            CI_upper = mean(retrieval) + qnorm(0.975) * (sd(retrieval) / sqrt(n())))



#Difference in behaviour between operators
injdata %>% 
  filter(dayofanalysis==as.POSIXct("2024-11-25",format = "%Y-%m-%d")) %>% 
  filter(!grepl("bot|cal|Dif",sample)) %>% 
  separate(sample, into = c("operator", "origin","replicate"), remove = F) %>% 
  arrange(unixtime_ofmax) %>% 
  group_by(operator, sample) %>% 
  mutate(
    lagpeak = unixtime_ofmax - lag(unixtime_ofmax),  # Calculate lag with previous row
    first_instance = peak_id,  # Value of the first instance
    second_instance = lag(peak_id)     # Value of the second instance
  ) %>%
  filter(!is.na(lagpeak)) %>%  # Remove NA rows where lag cannot be calculated
  ungroup() %>% 
  ggplot(aes(x=operator, y=lagpeak))+
  geom_boxplot()+
  geom_jitter()+
  scale_y_continuous(breaks = seq(20,120,by=10))
  
  
  




#Example concentrations Samples R4Cs
injdata %>% 
  filter(grepl("^S|air", sample)) %>% 
  ggplot(aes(x=sample, y=n2o_ppm))+
  geom_point()+
  scale_x_discrete(name="Field samples + Lab air")+
  scale_y_continuous(name = "N2O (ppm)")+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )





#Diffusion test:
injdata %>% 
  filter(dayofanalysis==as.POSIXct("2024-11-25",format = "%Y-%m-%d")) %>% 
  filter(grepl("Dif",sample)) %>% 
  separate(sample, into=c("diftest", "conc", "difusion_time"), remove = F) %>% 
  mutate(difusion_time=parse_number(difusion_time)) %>% 
  ggplot(aes(x=difusion_time,y=n2o_ppm/6*100, col=as.factor(difusion_time)))+
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  labs(col="Diffusion seconds",x="Diffusion seconds",y="Recovery")+
  ggtitle("Losses by difusion (6ppm from bottle)")+
  theme_bw()



ggplot(sampledata, aes(x=n2o_ppm_avg, y=n2o_ppm_cv*100, group = nobs))+
  geom_jitter(width = 0.1)



#Calculate average lab air baseline across days:
lab_air_mean<- mean(baselinedata$base_avg)/1000
lab_air_sd<- sd(baselinedata$base_avg)/1000



injdata %>% 
  filter(grepl("^aire", sample)) %>% 
  ggplot(aes(x=sample, y=n2o_ppm,col=sample))+
  geom_point()+
  geom_hline(yintercept = lab_air_mean)+
  geom_hline(yintercept = lab_air_mean+lab_air_sd, linetype=2)+
  geom_hline(yintercept = lab_air_mean-lab_air_sd, linetype=2)+
  facet_wrap(.~dayofanalysis)+
  ggtitle("Lab air injections vs mean baseline")



injdata %>% 
  filter(grepl("^aire", sample)) %>% 
  ggplot(aes(x=sample, y=n2o_ppm/lab_air_mean,col=sample))+
  geom_point()+
  geom_hline(yintercept = lab_air_mean/lab_air_mean)+
  geom_hline(yintercept = (lab_air_mean+lab_air_sd)/lab_air_mean, linetype=2)+
  geom_hline(yintercept = (lab_air_mean-lab_air_sd)/lab_air_mean, linetype=2)+
  facet_wrap(.~dayofanalysis)+
  ggtitle("%recovery Lab air injections vs mean baseline")



injdata %>% 
  filter(grepl("^aire", sample)) %>% 
  ggplot(aes(x="All_inections", y=n2o_ppm/lab_air_mean))+
  geom_boxplot()+
  geom_point(aes(col=sample))+
  geom_hline(yintercept = lab_air_mean/lab_air_mean)+
  geom_hline(yintercept = (lab_air_mean+lab_air_sd)/lab_air_mean, linetype=2)+
  geom_hline(yintercept = (lab_air_mean-lab_air_sd)/lab_air_mean, linetype=2)+
  facet_wrap(.~dayofanalysis)+
  ggtitle("%recovery Lab air injections vs mean baseline")




injdata %>% 
  filter(grepl("^6ppm-[0-9]", sample)) %>% 
  ggplot(aes(x=sample, y=n2o_ppm,col=sample))+
  geom_point()+
  geom_hline(yintercept = 6)+
  facet_wrap(.~dayofanalysis)+
  ggtitle("Standard injections vs nominal N2O")



injdata %>% 
  filter(grepl("^6ppm-[0-9]", sample)) %>% 
  ggplot(aes(x=sample, y=n2o_ppm/6,col=sample))+
  geom_point()+
  geom_hline(yintercept = 1)+
  facet_wrap(.~dayofanalysis)+
  ggtitle("%recovery Standard injections")



injdata %>% 
  filter(grepl("^6ppm-[0-9]", sample)) %>% 
  ggplot(aes(x="All_inections", y=n2o_ppm/6))+
  geom_boxplot()+
  geom_point(aes(col=sample))+
  geom_hline(yintercept = 1)+
  facet_wrap(.~dayofanalysis)+
  ggtitle("%recovery Standard injections")


injdata %>% 
  filter(grepl("^6ppm-[0-9]|^aire", sample)) %>% 
  mutate(n2o_esperado=case_when(grepl("^6ppm-[0-9]",sample)~6,
                                grepl("^aire",sample)~lab_air_mean)) %>%
    ggplot(aes(x=sample, y=n2o_ppm/n2o_esperado, col=dayofanalysis))+
  geom_point()+
  geom_hline(yintercept = 1)

injdata %>% 
  filter(grepl("^6ppm-[0-9]|^aire", sample)) %>% 
  mutate(n2o_esperado=case_when(grepl("^6ppm-[0-9]",sample)~6,
                                grepl("^aire",sample)~lab_air_mean)) %>%
  mutate(recovery=n2o_ppm/n2o_esperado*100, 
         injtype=gsub("-[0-9]","",sample)) %>%
  group_by(injtype, dayofanalysis) %>% 
  summarise(mean_recovery=mean(recovery), sd_recovery=sd(recovery))

#%recovery is ~90% for standardgas and ~95% for air

#Expected recoveries ~100 for both types of injections. No purge problem, as recovery for air is as bad as for standard.

#Can be related with:
# 1.calibration curve used != 1ml-syringe
# 2.N2O reacts when stored in plastic syringes
# 3.User error (user in Calcurve is different than user in standard injections) Some injections are good for air recovery but others are below

#Different spread between injections part of the same injection sequence (CV<5%) and those of same gas origin but separated throughout the day


