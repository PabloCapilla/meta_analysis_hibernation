###
###
#' Script to prepare data for hibernation meta-analysis
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
  "./01_data/Hibernation_dataset_4_2_clean.xlsx",
  colNames = T,
  sheet = 3
)
head(data)
summary(data)

#####

##
##### remove duplicated entries (due to how arousal data was collected) #####
##
data <- data |>
  filter(!is.na(Hib_M)) |>
  group_by(Study_ID, Tissue, Hib_M) |>
  filter(row_number() == 1)

#####

##
##### Descriptive summary of final data and checks to catch errors, typos, etc #####
##
summary(data$Euthermia_M)
summary(data$Euthermia_SD)

table(data$Class)
table(data$Thermoregulation)
table(data$Biomarker_category)

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

data$Hib_SD <- ifelse(
  is.na(data$Hib_SD),
  data$Hib_SE * sqrt(data$Hib_N),
  data$Hib_SD
)

data$Arousal_SD <- ifelse(
  is.na(data$Arousal_SD),
  data$Arousal_SE * sqrt(data$Arousal_N),
  data$Arousal_SD
)

####

##
##### Calculating meta-analysis y variables #####
##
data$Hib_M
data$Euthermia_M

## SMD
data <- as.data.frame(escalc(
  measure = "SMDH",
  n1i = Hib_N,
  m1i = Hib_M,
  sd1i = Hib_SD,

  n2i = Euthermia_N,
  m2i = Euthermia_M,
  sd2i = Euthermia_SD,

  data = data,
  var.names = c("SMDH", "SMDH.sv"),
  add.measure = F,
  append = TRUE
))

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
  var.names = c("SMDH_arousal", "SMDH.sv_arousal"),
  add.measure = F,
  append = TRUE
))

## checking very high SMD values
data %>%
  arrange(desc(SMDH)) %>%
  head()

## list of papers to check sd values - sent to Pablo B
df_check <- data %>%
  filter(!is.na(SMDH)) %>%
  arrange(desc(SMDH)) %>%
  filter(SMDH < -4 | SMDH > 4)
write.csv(
  x = df_check,
  file = './01_data/data_check/GC_papers_check_high_low_var.csv'
)

## list of papers with repeated Eu_M values per study and tissue - sent to Pablo B to be checked
df_eu_rep <- data |>
  filter(!is.na(Euthermia_M)) |>
  group_by(Study_ID, Tissue, Euthermia_M) |>
  count() |>
  arrange(desc(n)) |>
  filter(n > 1)

write.csv(
  x = df_eu_rep,
  file = './01_data/data_check/GC_papers_check_repeated_Eu_M_SD.csv'
)

##
##
##### Save full table #####
saveRDS(
  object = data,
  file = "./01_data/processed_RDS_data_files/03_GC_metaanalysis_full_data.RDS"
)

#####
