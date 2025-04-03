###
###
#' Script to prepare data for hibernation meta-analysis
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())



##
##### libraries #####
##

#renv::init(s)
#renv::restore()

library(openxlsx)
library(metafor)
library(tidyverse)
library(ggokabeito)

#renv::snapshot()

##
##### data #####
##

##
## study specific data
##
data  <- read.xlsx("./01_data/Estivation_dataset_2_0.xlsx",
                   colNames=T,
                   sheet = 1)
head(data)
summary(data)

#####

##
##### Descriptive summary of final data and checks to catch errors, typos, etc #####
##
summary(data$Euthermia_M)
summary(data$Euthermia_SD)

table(data$Class)
table(data$Class_2)
table(data$Thermoregulation)
table(data$Biomarker_cat_2)

#####

##
##
##### Complete SD if missing based on SE #####
##
##
data$Euthermia_SD <- ifelse(is.na(data$Euthermia_SD), 
                            data$Euthermia_SE * sqrt(data$Euthermia_N), 
                            data$Euthermia_SD)

data$Hib_SD <- ifelse(is.na(data$Hib_SD), 
                      data$Hib_SE * sqrt(data$Hib_N), 
                      data$Hib_SD)

data$Arousal_SD <- ifelse(is.na(data$Arousal_SD), 
                          data$Arousal_SE * sqrt(data$Arousal_N), 
                          data$Arousal_SD)

####


##
##### Calculating meta-analysis y variables #####
##

## SMD
data <- as.data.frame(escalc(measure="SMDH", 
                             n1i=Hib_N, 
                             m1i=Hib_M, 
                             sd1i=Hib_SD,
                             
                             n2i=Euthermia_N, 
                             m2i=Euthermia_M, 
                             sd2i= Euthermia_SD, 
                             
                             data=data,
                             var.names=c("SMDH","SMDH.sv"), 
                             add.measure=F,
                             append=TRUE))

## SMD
data <- as.data.frame(escalc(measure="SMDH",
                             n1i=Arousal_N, 
                             m1i=Arousal_M, 
                             sd1i=Arousal_SD,
                             
                             n2i=Euthermia_N, 
                             m2i=Euthermia_M, 
                             sd2i= Euthermia_SD, 
                             
                             data=data,
                             var.names=c("SMDH_arousal","SMDH.sv_arousal"), 
                             add.measure=F,
                             append=TRUE))

## checking very high SMD values
data %>% 
  arrange(desc(SMD)) %>% 
  head()

df_check <- data %>% 
  filter(!is.na(SMDH)) %>% 
  arrange(desc(SMDH)) %>%
  filter(SMDH < -4 | SMDH > 4)

write.csv(x = df_check, file = './01_data/data_check/GC_aestivation_papers_to_check.csv')


##
##
##### Save full table #####
saveRDS(object = data, file = "./01_data/processed_RDS_data_files/04_OS_metaanalysis_estivation_data.RDS")

#####













