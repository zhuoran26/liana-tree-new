## This script assesses the sensitivity of our meta-analysis to differences
## in the geographic extent of tree and liana observations.
## Specifically, we subset the meta-analysis observations to only literature 
## measuring both trees and lianas in the same location.
## We compute the Mann-Whitney and Glass' Delta tests on this data subset.

## This analysis is not shown in the manuscript.

## Author: AM Willson
## Date modified: 19 January 2022

rm(list = ls())

library(tidyverse)
library(gt)
library(effectsize)

data = read.csv('geographic_met_analysis_10_Jan_22.csv')

data$Units_K = as.character(data$Units_K)

# Convert to consistent units
for(i in 1:nrow(data)){
  if(data$Units_K[i] == 'kg/m/s/Mpa'){
    # Convert from kg to mmol H2O
    data$K[i] = data$K[i] * 1000 / 18 * 1000
    data$Units_K[i] = 'mmol/m/s/Mpa'
  }else{
    if(data$Units_K[i] == 'mol/m/s/Mpa'){
      # Convert from mol to mmol H2O
      data$K[i] = data$K[i] * 1000
      data$Units_K[i] = 'mmol/m/s/Mpa'
    }else{
      print('Unknown units')
      print(data$Units_K[i])
    }
  }
}

data$Units_K = as.factor(data$Units_K)

# Make column for mol/m/s/MPa
data$K.mol.m.s.MPa = data$K * 0.001

# Remove P50 > -0.75 consistent with Trugman et al. 2020
data = data %>%
  mutate(P50, replace(P50, P50 > -0.75, NA)) %>%
  dplyr::select(-P50)

# Rename column with NAs
colnames(data)[ncol(data)] = 'P50'

######################
## Table 1: U tests ##
######################

sig_tests = as_tibble(matrix(, nrow = 3, ncol = 7))
colnames(sig_tests) = c('Trait', 'Tree mean', 'Liana mean', 'n tree', 'n liana', 'Test Statistic', 'p-value')

sig_tests$Trait[1] = 'Stem-specific hydraulic conductivity (mol/m/s/MPa)'
sig_tests$`Tree mean`[1] = mean(data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[1] = mean(data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[1] = length(data$K.mol.m.s.MPa[which(data$Growth.form == 'tree' & !is.na(data$K.mol.m.s.MPa))])
sig_tests$`n liana`[1] = length(data$K.mol.m.s.MPa[which(data$Growth.form == 'liana' & !is.na(data$K.mol.m.s.MPa))])
sig_tests$`Test Statistic`[1] = wilcox.test(data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')],  
                                            data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[1] = wilcox.test(data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')], 
                                     data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')])$p.value

sig_tests$Trait[2] = 'P<sub>50</sub> (MPa)'
sig_tests$`Tree mean`[2] = mean(data$P50[which(data$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[2] = mean(data$P50[which(data$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[2] = length(data$P50[which(data$Growth.form == 'tree' & !is.na(data$P50))])
sig_tests$`n liana`[2] = length(data$P50[which(data$Growth.form == 'liana' & !is.na(data$P50))])
sig_tests$`Test Statistic`[2] = wilcox.test(data$P50[which(data$Growth.form == 'tree')],
                                            data$P50[which(data$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[2] = wilcox.test(data$P50[which(data$Growth.form == 'tree')],
                                     data$P50[which(data$Growth.form == 'liana')])$p.value

sig_tests$Trait[3] = 'Slope of PLC curve (%/MPa)'
sig_tests$`Tree mean`[3] = mean(data$Slope[which(data$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[3] = mean(data$Slope[which(data$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[3] = length(data$Slope[which(data$Growth.form == 'tree' & !is.na(data$Slope))])
sig_tests$`n liana`[3] = length(data$Slope[which(data$Growth.form == 'liana' & !is.na(data$Slope))])
sig_tests$`Test Statistic`[3] = wilcox.test(data$Slope[which(data$Growth.form == 'tree')],
                                            data$Slope[which(data$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[3] = wilcox.test(data$Slope[which(data$Growth.form == 'tree')],
                                     data$Slope[which(data$Growth.form == 'liana')])$p.value

sig_tests$`p-value`
sig_tests$`p-value` = c('2.08 x 10<sup>-3</sup>',
                        '1.75 x 10<sup>-1</sup>',
                        '1.85 x 10<sup>-1</sup>')
sig_tests %>% gt() %>%
  fmt_markdown(columns = c('Trait', 'p-value')) %>%
  fmt_number(columns = c('Tree mean', 'Liana mean'), decimals = 2) %>%
  fmt_number(columns = c('n tree', 'n liana'), decimals = 0) %>%
  fmt_number(columns = c('Test Statistic'), decimals = 0) %>%
  cols_label('Tree mean' = html('&mu;<sub><sub>tree</sub></sub>')) %>%
  cols_label('Liana mean' = html('&mu;<sub><sub>liana</sub></sub>')) %>%
  cols_label('n tree' = html('n<sub><sub>tree</sub></sub>')) %>%
  cols_label('n liana' = html('n<sub><sub>liana</sub></sub>')) %>%
  tab_header(
    title = md('Mann-Whitney U Tests for extended meta-analysis')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(90),
    heading.title.font.size = 22) %>%
  tab_style(
    style = cell_fill(color = 'lightyellow'),
    locations = cells_body(rows = 1)) %>%
  cols_align(align = 'center',
             columns = c('Tree mean', 'Liana mean', 'n tree', 'n liana', 'Test Statistic', 'p-value')) %>%
  cols_align(align = 'left',
             columns = 'Trait') %>%
  gtsave(filename = 'Plots/extended_mann-whitney_tests_geographic.png')

##########################
## Table 2: Effect size ##
##########################

effect_size = as_tibble(matrix(, nrow = 3, ncol = 4))
colnames(effect_size) = c('Trait', 'Glass Delta', 'Lower CI', 'Upper CI')

effect_size$Trait[1] = 'Stem-specific hydraulic conductivity (mol/m/s/MPa)'
effect_size$`Glass Delta`[1] = glass_delta(data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')], data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[1] = glass_delta(data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')], data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[1] = glass_delta(data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')], data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')])$CI_high

effect_size$Trait[2] = 'P<sub>50</sub> (MPa)'
effect_size$`Glass Delta`[2] = glass_delta(data$P50[which(data$Growth.form == 'liana')], data$P50[which(data$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[2] = glass_delta(data$P50[which(data$Growth.form == 'liana')], data$P50[which(data$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[2] = glass_delta(data$P50[which(data$Growth.form == 'liana')], data$P50[which(data$Growth.form == 'tree')])$CI_high

effect_size$Trait[3] = 'Slope of PLC curve (%/MPa)'
effect_size$`Glass Delta`[3] = glass_delta(data$Slope[which(data$Growth.form == 'liana')], data$Slope[which(data$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[3] = glass_delta(data$Slope[which(data$Growth.form == 'liana')], data$Slope[which(data$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[3] = glass_delta(data$Slope[which(data$Growth.form == 'liana')], data$Slope[which(data$Growth.form == 'tree')])$CI_high

effect_size %>% gt() %>%
  fmt_markdown(columns = c('Trait')) %>%
  fmt_number(columns = c('Glass Delta', 'Lower CI', 'Upper CI'), decimals = 2) %>%
  cols_label('Glass Delta' = html("Glass' &Delta;")) %>%
  tab_header(
    title = md('Effect size for extended meta-analysis')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(90),
    heading.title.font.size = 22) %>%
  tab_style(
    style = cell_fill(color = 'lightyellow'),
    locations = cells_body(rows = 1)) %>%
  cols_align(align = 'center',
             columns = c('Glass Delta', 'Lower CI', 'Upper CI')) %>%
  cols_align(align = 'left',
             columns = c('Trait')) %>%
  gtsave(filename = 'Plots/extended_effect_size_tests_geographic.png')
