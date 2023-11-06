## This script adds species name to the TRY output
## This will allow for manual filtering by species name as described in the methods of the manuscript
## The filtered dataframe was then processed
## The output of this file is the final matrix of species-averaged trait values in TRY

## This script requires the following inputs:
    ## 1. output_9oct.csv: not given but created following the steps from TRY_analysis_1.R
    ## 2. 7754.txt: not given but can be downloaded from the TRY database

## Author: AM Willson
## Date modified: 16 March 2021

library(tidyverse)
library(data.table)

rm(list = ls())

# Processed TRY data from full_TRY_analysis.R
# Note that this file is not available but can be made available upon request
data = read.csv('output_9oct.csv')

# Remove column that doesn't mean anything
data = data[,-1]

# Load data with species name and ID numbers
# Note: this file is not available but can be downloaded from TRY (global species names)
spp.name = fread('TRY_analysis/Data/7754.txt')
# Keep only species name and ID number columns and rename
name.id = spp.name[,c(6:7)]
colnames(name.id) = c('Species.ID', 'Name')

# Match the species names to the rest of the data
data = data %>%
  full_join(name.id, by = 'Species.ID')

# Subset to only keep tree and liana species to speed computation
data_final = data %>%
  filter(Growth.form %in% c('tree', 'liana'))

# Keep only the variables from TRY that were used in the final analysis/shown in paper
# The following 2 functions were redone at ~1M row intervals to create 10 CSVs
data_after_talargada = data_final %>%
  filter(!is.nan(Leaf.lifespan..m.) | 
           !is.nan(SLA..mm2.mg.) | 
           !is.nan(Narea..g.m2.) | 
           !is.nan(Aarea..mmol.CO2.m2.s) |
           !is.nan(SSD..g.cm3.) |
           !is.nan(SSC..kg.m.s.MPa.) |
           !is.nan(SRL..m.g.) |
           !is.nan(FRD..mm.) |
           !is.nan(MC....) |
           !is.nan(WVD..um.) |
           !is.nan(Dens..1.mm.2.) |
           !is.nan(RRD..m.) |
           !is.nan(P50..MPa.) |
           !is.nan(Nmass..mg.g.) |
           !is.nan(Pmass..mg.g.) |
           !is.nan(LA..cm.2.)) %>%
  filter(Species.ID > 73094)

write.csv(, file = 'TRY_output_with_species_10.csv')

# Load in all 10 CSVs
data1 = read.csv('TRY_output_with_species.csv')
data2 = read.csv('TRY_output_with_species_2.csv')
data3 = read.csv('TRY_output_with_species_3.csv')
data4 = read.csv('TRY_output_with_species_4.csv')
data5 = read.csv('TRY_output_with_species_5.csv')
data6 = read.csv('TRY_output_with_species_6.csv')
data7 = read.csv('TRY_output_with_species_7.csv')
data8 = read.csv('TRY_output_with_species_8.csv')
data9 = read.csv('TRY_output_with_species_9.csv')
data10 = read.csv('TRY_output_with_species_10.csv')

# Clean working environment
# Too much stuff right now to complete the next step
rm(data,
   data_after_lkaempferi,
   data_after_pacuminatum,
   data_after_psylvestris,
   data_after_qilex,
   data_after_spolyphylla,
   data_after_talargada,
   data_final,
   name.id,
   spp.name)

# Combine all dataframes
data_filtered = rbind(data1, data2, data3, data4, data5,
                      data6, data7, data8, data9, data10)

# Remove parital dataframes
rm(data1, data2, data3, data4, data5,
   data6, data7, data8, data9, data10)

# Fix multiple row problem from original dataset
unique_id = unique(data_filtered$Species.ID)
count = 0

for(i in unique_id){
  count = count + 1
  sub = data_filtered[which(data_filtered$Species.ID == i),]
  row = sub[1,]
  if(i == unique_id[1]){
    data_final = row
  }else{
    data_final = rbind(data_final, row)
  }
  print(count / length(unique_id))
}

# Save final dataframe
write.csv(data_final, 'data/filtered_TRY_analysis_16-03-21.csv')