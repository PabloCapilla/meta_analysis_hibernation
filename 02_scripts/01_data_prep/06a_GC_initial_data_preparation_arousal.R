###
###
#' Script to prepare data for hibernation meta-analysis (hibernation effect size)
#'
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list = ls())


##
##### libraries #####
##

library(openxlsx)
library(metafor)
library(tidyverse)
library(ggokabeito)

##
##### data #####
##

##
## study specific data
##
data <- read.xlsx(
  "./01_data/raw_data/Hibernation_dataset_4_3_arousal.xlsx",
  colNames = T,
  sheet = 3
)
head(data)
summary(data)

#####

##
##### Descriptive summary of final data and checks to catch errors, typos, etc #####
##
summary(data$Arousal_M)
summary(data$Arousal_SD)

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
data$Euthermia_SD <- ifelse(
  is.na(data$Euthermia_SD),
  data$Euthermia_SE * sqrt(data$Euthermia_N),
  data$Euthermia_SD
)

data$Arousal_SD <- ifelse(
  is.na(data$Arousal_SD),
  data$Hib_SE * sqrt(data$Arousal_N),
  data$Arousal_SD
)

####

##
##### Calculating meta-analysis y variables #####
##

## SMD
data <- as.data.frame(escalc(
  measure = "SMDH",
  n1i = Arousal_N,
  m1i = Arousal_M,
  sd1i = Arousal_SD,

  n2i = Euthermia_N,
  m2i = Euthermia_M,
  sd2i = Euthermia_SD,

  data = data,
  var.names = c("SMDH", "SMDH.sv"),
  add.measure = F,
  append = TRUE
))


## checking very high SMD values
data %>%
  arrange(desc(SMDH)) %>%
  head()

## list of papers to check sd values - checked by PB
df_check <- data %>%
  filter(!is.na(SMDH)) %>%
  arrange(desc(SMDH)) %>%
  filter(SMDH < -4 | SMDH > 4)

##
##
##### Save full table #####
saveRDS(
  object = data |> filter(!is.na(SMDH)),
  file = "./01_data/processed_RDS_data_files/06_GC_metaanalysis_arousal_full_data.RDS"
)
#####
