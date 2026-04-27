##
##
## script to compile list of all (unique) studies included
##
##

# libraries
library(openxlsx)
library(tidyverse)

# read dataset
data1 <- readRDS("./01_data/processed_RDS_data_files/01_OS_metaanalysis_full_data.RDS")
data2 <- readRDS("./01_data/processed_RDS_data_files/02_IM_metaanalysis_full_data.RDS")
data3 <- readRDS("./01_data/processed_RDS_data_files/03_GC_metaanalysis_full_data.RDS")
data4 <- readRDS("./01_data/processed_RDS_data_files/04_OS_metaanalysis_estivation_data.RDS")

df_studies <- rbind(select(data1, Study_ID, Reference), 
                    select(data2, Study_ID, Reference = Study), 
                    select(data3, Study_ID, Reference), 
                    select(data4, Study_ID, Reference))

df_studies_unique <- df_studies %>% 
  group_by(Study_ID) %>% 
  filter(row_number() == 1)

head(df_studies_unique)

write.csv(x = df_studies_unique, file = './01_data/LIST UNIQUE STUDIES INCLUDED.csv')
