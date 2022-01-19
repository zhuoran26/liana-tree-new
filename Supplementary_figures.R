## This script provides the code for all the supplementary figures
## including Extended Data Figures & Supplementary Figures

## This script requires the following inputs:
    ## 1. Supplement_parameters.csv: spreadsheet created in Excel to 
    ## create Extended Data Table 4
    ## 2. Sources_2.csv: spreadsheet created in Excel to create
    ## Extended Data Table 3
    ## 3. full_met_analysis_data_8_Feb.csv: not given but can be
    ## recreated following the steps detailed in the Methods
    ## 4. param.input.RData: created in Input_parameter_est.R
    ## 5. bci_met_mxh.RData: created in BCI_met_process.R
    ## 6. horizontes_met_mxh.RData: created in Horizontes_met_process.R
    ## 7. filtered_TRY_analysis_16-03-21.csv: created in TRY_analysis1.R and 
    ## TRY_analysis2.R and filtered by hand
    ## 8. interpolation_mxh_100.RData: created in Sensitivity_met_interpolation.R
    ## 9.  tree.NPP.mxh & liana.NPP.mxh: compiled in NPP_models.R
    ## 10. tree.NPP.mxh.ca & liana.NPP.mxh.ca: compiled in NPP_models_wCO2.R

## Author: AM Willson
## Date modified: 11 October 2021

library(ggplot2)
library(ggsci)
library(dplyr)
library(gridExtra)
library(grid)
library(gt)
library(effectsize)
library(pracma)
library(reshape2)
library(cowplot)
library(DescTools)

#####################################
#### Part 1: Inputs & parameters ####
#####################################

## gt() version of Supplementary Table 4

rm(list = ls())

## Parameters

# Load table from CSV
params = read.csv('Supplement_parameters.csv')

# Make character objects
params$Name = as.character(params$Name)
params$Definition = as.character(params$Definition)
params$Value = as.character(params$Value)
params$Units = as.character(params$Units)
params$Source = as.character(params$Source)
params$Tree.or.liana.function = as.character(params$Tree.or.liana.function)

# Make table
params = as_tibble(params)

# Split table into two 
# (for formatting in Word)
params_a = params[1:10,]
params_b = params[11:20,]

params_a %>%
  gt() %>%
  fmt_markdown(columns = c('Name', 'Definition', 'Value', 'Units', 'Source', 'Tree.or.liana.function')) %>%
  tab_header(
    title = md('Changed Model Parameters')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(90),
    heading.title.font.size = 22) %>%
  cols_label('Tree.or.liana.function' = 'Tree or Liana Function') %>%
  cols_width('Units' ~ px(150)) %>%
  cols_width('Value' ~ px(200)) %>%
  cols_width('Definition' ~ px(200)) %>%
  cols_align(align = 'center') %>%
  gtsave(filename = 'Plots/SupplementaryTable4a.png')

params_b %>%
  gt() %>%
  fmt_markdown(columns = c('Name', 'Definition', 'Value', 'Units', 'Source', 'Tree.or.liana.function')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(90),
    heading.title.font.size = 22) %>%
  cols_label('Tree.or.liana.function' = 'Tree or Liana Function') %>%
  cols_width('Units' ~ px(150)) %>%
  cols_width('Value' ~ px(200)) %>%
  cols_width('Definition' ~ px(200)) %>%
  cols_align(align = 'center') %>%
  gtsave(filename = 'Plots/SupplementaryTable4b.png')

## gt() version of Supplementary Table 3

## Data sources
source = read.csv('Sources_2.csv')

source$Source = as.character(source$Source)
source$K = as.character(source$K)
source$P50 = as.character(source$P50)
source$Slope = as.character(source$Slope)

# Split into two
# (for formatting in Word)
source_a = source[1:8,]
source_b = source[9:15,]

source_a %>%
  gt() %>%
  fmt_markdown(columns = c('Source')) %>%
  tab_header(
    title = md('Extended meta-analysis data sources')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(80),
    heading.title.font.size = 22) %>%
  cols_label('K' = html('K<sub><sub>s</sub></sub>')) %>%
  cols_label('P50' = html('P<sub><sub>50</sub></sub>')) %>%
  cols_label('n.liana' = html('n<sub><sub>liana</sub></sub>')) %>%
  cols_label('n.tree' = html('n<sub><sub>tree</sub></sub>')) %>%
  cols_align(align = 'center',
             columns = c('K', 'P50', 'Slope', 'n.liana', 'n.tree')) %>%
  cols_align(align = 'left',
             columns = c('Source')) %>%
  cols_width('Source' ~ px(500)) %>%
  gtsave(filename = 'Plots/SupplementaryTable3a.png')

source_b %>%
  gt() %>%
  fmt_markdown(columns = c('Source')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(80),
    heading.title.font.size = 22) %>%
  cols_label('K' = html('K<sub><sub>s</sub></sub>')) %>%
  cols_label('P50' = html('P<sub><sub>50</sub></sub>')) %>%
  cols_label('n.liana' = html('n<sub><sub>liana</sub></sub>')) %>%
  cols_label('n.tree' = html('n<sub><sub>tree</sub></sub>')) %>%
  cols_align(align = 'center',
             columns = c('K', 'P50', 'Slope', 'n.liana', 'n.tree')) %>%
  cols_align(align = 'left',
             columns = c('Source')) %>%
  cols_width('Source' ~ px(500)) %>%
  gtsave(filename = 'Plots/SupplementaryTable3b.png')

############################################
#### Part 2: Hydraulic trait trade-offs ####
############################################

## Supplementary Figure 4

# Clean up working environment for next part
rm(list = ls())

# Load in new data
hydro_data = read.csv('full_met_analysis_data_8_Feb.csv')

# Make units easier to work with
hydro_data$Units_K = as.character(as.factor(hydro_data$Units_K))

# Convert units for K measurements
for(i in 1:nrow(hydro_data)){
  if(hydro_data$Units_K[i] == 'kg/m/s/Mpa'){
    hydro_data$K[i] = hydro_data$K[i] * 1000 / 18 * 1000
    hydro_data$Units_K[i] = 'mmol/m/s/Mpa'
  }else{
    if(hydro_data$Units_K[i] == 'mol/m/s/Mpa'){
      hydro_data$K[i] = hydro_data$K[i] * 1000
      hydro_data$Units_K[i] = 'mmol/m/s/Mpa'
    }else{
      if(hydro_data$Units_K[i] == ''){
        next
      }else{
        print('Error with units')
      }
    }
  }
}

# Remove P50 > -0.75 consistent with Trugman et al. 2020
hydro_data = hydro_data %>%
  mutate(P50, replace(P50, P50 > -0.75, NA)) %>%
  dplyr::select(-P50)

# Rename column with NAs
colnames(hydro_data)[ncol(hydro_data)] = 'P50'

K = hydro_data %>%
  filter(!is.na(K)) %>%
  dplyr::select(Growth.form, K)

P50 = hydro_data %>%
  filter(!is.na(P50)) %>%
  dplyr::select(Growth.form, P50)

Slope = hydro_data %>%
  filter(!is.na(Slope)) %>%
  dplyr::select(Growth.form, Slope)

# Log transform variables
# Notice that P50 is now positive, but it does not change the fundamental relationship
hydro_data$log_P50 = log(-hydro_data$P50)
hydro_data$log_K = log(hydro_data$K)
hydro_data$log_slope = log(hydro_data$Slope)

# Summarize for both growth forms
all_lm = summary(lm(hydro_data$log_P50 ~ hydro_data$log_K))
all_lm$r.squared

cor = cor.test(hydro_data$log_K, hydro_data$log_P50, method = 's', exact = F)
cor$p.value

# Trees only
tree_hydro = hydro_data %>%
  filter(Growth.form == 'tree')

tree_lm = summary(lm(tree_hydro$log_P50 ~ tree_hydro$log_K))
tree_lm$r.squared

tree_cor = cor.test(tree_hydro$log_K, tree_hydro$log_P50, method = 's', exact = F)
tree_cor$p.value

# Lianas only
liana_hydro = hydro_data %>%
  filter(Growth.form == 'liana')

liana_lm = summary(lm(liana_hydro$log_P50 ~ liana_hydro$log_K))
liana_lm$r.squared

liana_cor = cor.test(liana_hydro$log_K, liana_hydro$log_P50, method = 's', exact = F)
liana_cor$p.value

pl1 = ggplot() +
  geom_point(aes(x = hydro_data$log_K, y = hydro_data$log_P50, color = hydro_data$Growth.form), size = 2) +
  geom_smooth(aes(x = hydro_data$log_K, y = hydro_data$log_P50), color = 'black', method = 'lm', se = F, size = 1.2) +
  geom_smooth(aes(x = tree_hydro$log_K, y = tree_hydro$log_P50, color = tree_hydro$Growth.form), method = 'lm', se = F, size = 1.2) +
  geom_smooth(aes(x = liana_hydro$log_K, y = liana_hydro$log_P50, color = liana_hydro$Growth.form), method = 'lm', se = F, size = 1.2) +
  xlab(expression(paste('log(K) (mmol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) + ylab(expression(paste('log(-',P[50],') (MPa)'))) +
  ggtitle('Hydraulic efficiency vs. safety') +
  theme_linedraw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_npg(labels = c('tree' = 'Tree', 'liana' = 'Liana')) +
  labs(color = 'Growth form') +
  scale_y_continuous(breaks = c(0, 1))
pl1

# Summarize for both growth forms
all_lm = summary(lm(hydro_data$log_slope ~ hydro_data$log_K))
all_lm$r.squared

cor = cor.test(hydro_data$log_K, hydro_data$log_slope, method = 's', exact = F)
cor$p.value

# Trees only
tree_lm = summary(lm(tree_hydro$log_slope ~ tree_hydro$log_K))
tree_lm$r.squared

tree_cor = cor.test(tree_hydro$log_K, tree_hydro$log_slope, method = 's', exact = F)
tree_cor$p.value

# Lianas only
liana_lm = summary(lm(liana_hydro$log_slope ~ liana_hydro$log_K))
liana_lm$r.squared

liana_cor = cor.test(liana_hydro$log_K, liana_hydro$log_slope, method = 's', exact = F)
liana_cor$p.value

pl2 = ggplot() +
  geom_point(aes(x = hydro_data$log_K, y = hydro_data$log_slope, color = hydro_data$Growth.form), size = 2) +
  geom_smooth(aes(x = hydro_data$log_K, y = hydro_data$log_slope), color = 'black', method = 'lm', se = F, size = 1.2) +
  geom_smooth(aes(x = tree_hydro$log_K, y = tree_hydro$log_slope, color = tree_hydro$Growth.form), method = 'lm', se = F, size = 1.2) +
  geom_smooth(aes(x = liana_hydro$log_K, y = liana_hydro$log_slope, color = liana_hydro$Growth.form), method = 'lm', se = F, size = 1.2) +
  xlab(expression(paste('log(K) (mmol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) + ylab(expression(paste('log(Slope) (% MP', a^-1, ')'))) +
  ggtitle('Hydraulic efficiency vs. slope') +
  theme_linedraw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_npg(labels = c('tree' = 'Tree', 'liana' = 'Liana')) +
  labs(color = 'Growth form')
pl2

# Summarize for both growth forms
all_lm = summary(lm(hydro_data$log_slope ~ hydro_data$log_P50))
all_lm$r.squared

cor = cor.test(hydro_data$log_P50, hydro_data$log_slope, method = 's', exact = F)
cor$p.value

# Trees only
tree_lm = summary(lm(tree_hydro$log_slope ~ tree_hydro$log_P50))
tree_lm$r.squared

tree_cor = cor.test(tree_hydro$log_P50, tree_hydro$log_slope, method = 's', exact = F)
tree_cor$p.value

# Lianas only
liana_lm = summary(lm(liana_hydro$log_slope ~ liana_hydro$log_P50))
liana_lm$r.squared

liana_cor = cor.test(liana_hydro$log_P50, liana_hydro$log_slope, method = 's', exact = F)
liana_cor$p.value

pl3 = ggplot() +
  geom_point(aes(x = hydro_data$log_P50, y = hydro_data$log_slope, color = hydro_data$Growth.form), size = 2) +
  geom_smooth(aes(x = hydro_data$log_P50, y = hydro_data$log_slope), color = 'black', method = 'lm', se = F, size = 1.2) +
  geom_smooth(aes(x = tree_hydro$log_P50, y = tree_hydro$log_slope, color = tree_hydro$Growth.form), method = 'lm', se = F, size = 1.2) +
  geom_smooth(aes(x = liana_hydro$log_P50, y = liana_hydro$log_slope, color = liana_hydro$Growth.form), method = 'lm', se = F, size = 1.2) +
  xlab(expression(paste('log(-',P[50],') (MPa)'))) + ylab(expression(paste('log(Slope) (% MP', a^-1, ')'))) +
  ggtitle('Hydraulic safety vs. slope') +
  theme_linedraw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_npg(labels = c('tree' = 'Tree', 'liana' = 'Liana')) +
  labs(color = 'Growth form')
pl3

# Function to extract legend from ggplot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Get legend
leg = get_legend(pl1)

# Plot together
ga = plot_grid(pl1 + theme(legend.position = 'none'), 
               pl2 + theme(legend.position = 'none'), 
               pl3 + theme(legend.position = 'none'),
               leg, nrow = 2, ncol = 2, 
               rel_widths = c(1, 1),
               labels = c('A', 'B', 'C', ''),
               label_size = 30)
ga

ggsave(ga, filename = 'Plots/SupplementaryFigure4.jpeg', height = 13, width = 13, unit = 'in')

###################################
#### Part 3: DBH distributions ####
###################################

## Supplementary Figure 5

rm(list = ls())

# Load parameters & input
load('param.input.RData')

# Plot tree DBH
pl1 = ggplot() +
  geom_histogram(aes(x = hori.data.tree.dbh), bins = 10, color = 'black', fill = 'lightgray') +
  theme_linedraw() +
  xlab('DBH (cm)') + ylab('Frequency') +
  ggtitle('Tree') +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), plot.title = element_text(size = 30, hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pl1

# Plot liana DBH
pl2 = ggplot() +
  geom_histogram(aes(x = dbh.data$DiamMay2018_cm), bins = 10, color = 'black', fill = 'lightgray') +
  theme_linedraw() +
  xlab('DBH (cm)') + ylab('Frequency') +
  ggtitle('Liana') +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), plot.title = element_text(size = 30, hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pl2

# Plot together
ga = plot_grid(pl1, pl2,
               nrow = 1,
               rel_widths = 1,
               labels = c('A', 'B'),
               label_size = 30)
ga

ggsave(ga, filename = 'Plots/SupplementaryFigure5.jpeg', width = 20, height = 10, units = 'in')

#################################
#### Part 4: climate drivers ####
#################################

## Supplementary Figure 6

rm(list = ls())

# Load climate data
load('bci_met_mxh.RData')
BCI_VPD = VPD
BCI_SWP = SWP
load('horizontes_met_mxh.RData')
H_VPD = VPD
H_SWP = SWP

# SWP for each site
pl1 = ggplot() +
  geom_line(aes(x = c(1:12), y = BCI_SWP, color = 'BCI'), size = 1.2) +
  geom_line(aes(x = c(1:12), y = H_SWP, color = 'Horizontes'), size = 1.2) +
  xlab('Month') + ylab(expression(paste(Psi,' (MPa)'))) +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  theme_linedraw() +
  scale_color_discrete(name = 'Site') +
  theme(axis.title = element_text(size = 28), 
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pl1

# Average VPD over hours
BCI_VPD_avg = apply(BCI_VPD, 1, mean)
H_VPD_avg = apply(H_VPD, 1, mean)

# Plot average VPD for each site
pl2 = ggplot() +
  geom_line(aes(x = c(1:12), y = BCI_VPD_avg, color = 'BCI'), size = 1.2) +
  geom_line(aes(x = c(1:12), y = H_VPD_avg, color = 'Horizontes'), size = 1.2) +
  xlab('Month') + ylab('VPD (Pa)') +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  theme_linedraw() +
  scale_color_discrete(name = 'Site') +
  theme(axis.title = element_text(size = 28), 
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pl2

# Melt hourly VPD data
BCI_VPD_melt = melt(BCI_VPD)
colnames(BCI_VPD_melt) = c('Month', 'Hour', 'VPD')
BCI_VPD_melt$Hour = BCI_VPD_melt$Hour - 1
H_VPD_melt = melt(H_VPD)
colnames(H_VPD_melt) = c('Month', 'Hour', 'VPD')
H_VPD_melt$Hour = H_VPD_melt$Hour - 1

# Plot hourly VPD for BCI
pl3 = ggplot() +
  geom_line(aes(x = BCI_VPD_melt$Month, y = BCI_VPD_melt$VPD, color = BCI_VPD_melt$Hour, group = BCI_VPD_melt$Hour)) +
  xlab('Month') + ylab('VPD (Pa)') +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  theme_linedraw() +
  scale_color_continuous(name = 'Hour', type = 'viridis') +
  ggtitle('Barro Colorado Island, Panama') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pl3

# Plot hourly VPD for Horizontes
pl4 = ggplot() +
  geom_line(aes(x = H_VPD_melt$Month, y = H_VPD_melt$VPD, color = H_VPD_melt$Hour, group = H_VPD_melt$Hour)) +
  xlab('Month') + ylab('VPD (Pa)') +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  theme_linedraw() +
  scale_color_continuous(name = 'Hour', type = 'viridis', breaks = c(0, 12, 23)) +
  ggtitle('Horizontes, Costa Rica') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pl4

# Plot together
ga1 = plot_grid(pl1 + theme(legend.position = 'none'), 
                pl2, nrow = 1,
                rel_widths = c(1, 1.4),
                labels = c('A', 'B'),
                label_size = 30)
ga1

ga2 = plot_grid(pl3 + theme(legend.position = 'none'), 
                pl4, nrow = 1,
                rel_widths = c(1, 1.2),
                labels = c('C', 'D'),
                label_size = 30)
ga2

ggsave(ga1, filename = 'Plots/SupplementaryFigure6a.jpeg', height = 8, width = 16, units = 'in')
ggsave(ga2, filename = 'Plots/SupplementaryFigure6b.jpeg', height = 8, width = 16, units = 'in')

#####################
#### Part 5: TRY ####
#####################

## Supplementary Figures 7-9

rm(list = ls())

## Load data

# Read in TRY data and clean the data frame
TRY = read.csv('filtered_TRY_analysis_16-03-21.csv')
TRY = TRY[,-1]
colnames(TRY) = c('Species.ID', 'Growth.form', 'LL.months', 
                  'SLA.mm2.mg', 'Narea.g.m2', 'Aarea.mmolC02.m2.s', 
                  'SSD.g.cm3', 'K.kg.m.s.MPa', 'SRL.m.g', 'FRD.mm', 
                  'CRD.mm', 'MC.%', 'Kroot.cm2.cm.s.MPa10.6', 'WVD.um', 
                  'Dens.1.mm2', 'H.m', 'LA.cm2', 'SM.mg', 'WVEL.um', 
                  'RRD.m', 'P50.MPa', 'Nmass.mg.g', 'Pmass.mg.g', 
                  'Crown.depth.m', 'Species')

TRY$K.mmol.m.s.MPa = TRY$K.kg.m.s.MPa * 1000 / 18 * 1000

# Create mol/m/s/MPa column
TRY$K.mol.m.s.MPa = TRY$K.mmol.m.s.MPa * 0.001

#########################
## Plot 1: Leaf traits ##
#########################

## Supplementary Figure 8

# Number of species-level observations
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$LL.months)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$LL.months)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$SLA.mm2.mg)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$SLA.mm2.mg)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$Narea.g.m2)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$Narea.g.m2)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$Aarea.mmolC02.m2.s)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$Aarea.mmolC02.m2.s)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$Nmass.mg.g)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$Nmass.mg.g)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$Pmass.mg.g)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$Pmass.mg.g)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$LA.cm2)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$LA.cm2)))

# Violin and boxplots overlaid with jittered points
pla = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, LL.months) %>%
  ggplot(aes(x = Growth.form, y = LL.months, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', shape = 3, show.legend = F, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab('Leaf lifespan (months)') +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 20', 'Tree\nn = 500')) +
  theme(axis.title.y = element_text(size = 15.5), axis.text = element_text(size = 15))

plb = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, SLA.mm2.mg) %>%
  ggplot(aes(x = Growth.form, y = SLA.mm2.mg, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste('SLA (m',m^2,' m', g^-1, ')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 82', 'Tree\nn = 1207')) +
  theme(axis.title.y = element_text(size = 15.5), axis.text = element_text(size = 15))

plc = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Narea.g.m2) %>%
  ggplot(aes(x = Growth.form, y = Narea.g.m2, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(N[area],' (g ',m^-2,')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 138', 'Tree\nn = 1716')) +
  theme(axis.title.y = element_text(size = 15.5), axis.text = element_text(size = 15))

pld = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Aarea.mmolC02.m2.s) %>%
  ggplot(aes(x = Growth.form, y = Aarea.mmolC02.m2.s, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(A[area],' (mmol C',O[2],' ',m^-2,' ', s^-1, ')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 91', 'Tree\nn = 1427')) +
  theme(axis.title.y = element_text(size = 14), axis.text = element_text(size = 15))

ple = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Nmass.mg.g) %>%
  ggplot(aes(x = Growth.form, y = Nmass.mg.g, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(N[mass],' (mg ', g^-1, ')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 265', 'Tree\nn = 5264')) +
  theme(axis.title.y = element_text(size = 15.5), axis.text = element_text(size = 15))

plf = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Pmass.mg.g) %>%
  ggplot(aes(x = Growth.form, y = Pmass.mg.g, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(P[mass],' (mg ', g^-1, ')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 163', 'Tree\nn = 3089')) +
  theme(axis.title.y = element_text(size = 15.5), axis.text = element_text(size = 15))

plg = TRY %>%
  filter(Growth.form %in% c('tree', 'liana')) %>%
  dplyr::select(Growth.form, LA.cm2) %>%
  ggplot(aes(x = Growth.form, y = LA.cm2, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste('Leaf area (c',m^2,')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 561', 'Tree\nn = 4029')) +
  theme(axis.title.y = element_text(size = 15.5), axis.text = element_text(size = 15))

ga_leaf = plot_grid(pla, plb, plc, pld, ple, plf, plg, 
                       nrow = 2, ncol = 4, 
                       labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G'),
                       label_size = 22)

ggsave(ga_leaf, filename = 'Plots/SupplementaryFigure8.jpeg', width = 12, height = 7, units = 'in')

#########################
## Plot 2: Stem traits ##
#########################

## Supplementary Table 7

# Number of species-level observations
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$SSD.g.cm3)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$SSD.g.cm3)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$WVD.um)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$WVD.um)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$Dens.1.mm2)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$Dens.1.mm2)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$K.mol.m.s.MPa)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$K.mol.m.s.MPa)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$P50.MPa)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$P50.MPa)))

pla = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, SSD.g.cm3) %>%
  ggplot(aes(x = Growth.form, y = SSD.g.cm3, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste('Stem specific density (g c',m^-3,')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 33', 'Tree\nn = 7682')) +
  theme(axis.title.y = element_text(size = 11), axis.text = element_text(size = 15))

plb = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, WVD.um) %>%
  ggplot(aes(x = Growth.form, y = WVD.um, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste('Vessel diameter (',mu,'m)'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 65', 'Tree\nn = 626')) +
  theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 15))

plc = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Dens.1.mm2) %>%
  ggplot(aes(x = Growth.form, y = Dens.1.mm2, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste('Vessel density (m',m^-2,')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 118', 'Tree\nn = 2009')) +
  theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 15))

pld = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, K.mol.m.s.MPa) %>%
  ggplot(aes(x = Growth.form, y = K.mol.m.s.MPa, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste('Conductivity (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 18', 'Tree\nn = 400')) +
  theme(axis.title.y = element_text(size = 10.5), axis.text = element_text(size = 15))

ple = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, P50.MPa) %>%
  ggplot(aes(x = Growth.form, y = P50.MPa, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(P[50],' (MPa)'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 4', 'Tree\nn = 211')) +
  theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 15))

ga_stem = plot_grid(pla, plb, plc, pld, ple, 
                    nrow = 2, ncol = 3, 
                    labels = c('A', 'B', 'C', 'D', 'E'),
                    label_size = 22)

ggsave(ga_stem, filename = 'Plots/SupplementaryFigure7.jpeg', width = 10, height = 7, units = 'in')

#########################
## Plot 3: Root traits ##
#########################

## Supplementary Table 9

# Number of species-level observations
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$SRL.m.g)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$SRL.m.g)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$FRD.mm)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$FRD.mm)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$`MC.%`)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$`MC.%`)))
length(which(TRY$Growth.form == 'liana' & !is.na(TRY$RRD.m)))
length(which(TRY$Growth.form == 'tree' & !is.na(TRY$RRD.m)))

pla = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, SRL.m.g) %>%
  ggplot(aes(x = Growth.form, y = SRL.m.g, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste('Specific root length (m ', g^-1, ')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 15', 'Tree\nn = 179')) +
  theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 15))

plb = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, FRD.mm) %>%
  ggplot(aes(x = Growth.form, y = FRD.mm, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab('Fine root diameter (mm)') +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 13', 'Tree\nn = 169')) +
  theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 15))

plc = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, `MC.%`) %>%
  ggplot(aes(x = Growth.form, y = `MC.%`, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab('Mycorrhizal colonization (%)') +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 7', 'Tree\nn = 55')) +
  theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 15))

pld = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, RRD.m) %>%
  ggplot(aes(x = Growth.form, y = RRD.m, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab('Rooting depth (m)') +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 6', 'Tree\nn = 264')) +
  theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 15))

ga_root = plot_grid(pla, plb, plc, pld, 
                    nrow = 2, ncol = 2,
                    labels = c('A', 'B', 'C', 'D'),
                    label_size = 22)

ggsave(ga_root, filename = 'Plots/SupplementaryFigure9.jpeg', width = 10, height = 7, units = 'in')

######################
## Table 1: U tests ##
######################

## gt() version of Supplementary Table 5

# Storage & formatting
sig_tests = as_tibble(matrix(, nrow = 16, ncol = 7))
colnames(sig_tests) = c('Trait', 'Tree mean', 'Liana mean', 'n tree', 'n liana', 'Test Statistic', 'p-value')

# Fill in table for wood density
sig_tests$Trait[1] = 'Stem specific density (g/cm<sup>3</sup>)'
sig_tests$`Tree mean`[1] = mean(TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[1] = mean(TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[1] = length(TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree' & !is.na(TRY$SSD.g.cm3))])
sig_tests$`n liana`[1] = length(TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana' & !is.na(TRY$SSD.g.cm3))])
sig_tests$`Test Statistic`[1] = wilcox.test(TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')],
                                            TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[1] = wilcox.test(TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')],
                                     TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for vessel diameter
sig_tests$Trait[2] = 'Vessel diameter (&mu;m)'
sig_tests$`Tree mean`[2] = mean(TRY$WVD.um[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[2] = mean(TRY$WVD.um[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[2] = length(TRY$WVD.um[which(TRY$Growth.form == 'tree' & !is.na(TRY$WVD.um))])
sig_tests$`n liana`[2] = length(TRY$WVD.um[which(TRY$Growth.form == 'liana' & !is.na(TRY$WVD.um))])
sig_tests$`Test Statistic`[2] = wilcox.test(TRY$WVD.um[which(TRY$Growth.form == 'tree')],
                                            TRY$WVD.um[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[2] = wilcox.test(TRY$WVD.um[which(TRY$Growth.form == 'tree')],
                                     TRY$WVD.um[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for vessel density
sig_tests$Trait[3] = 'Vessel density (1/mm<sup>2</sup>)'
sig_tests$`Tree mean`[3] = mean(TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[3] = mean(TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[3] = length(TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree' & !is.na(TRY$Dens.1.mm2))])
sig_tests$`n liana`[3] = length(TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana' & !is.na(TRY$Dens.1.mm2))])
sig_tests$`Test Statistic`[3] = wilcox.test(TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')],
                                            TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[3] = wilcox.test(TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')],
                                     TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for K
sig_tests$Trait[4] = 'Stem specific hydraulic conductivity (mol/m/s/MPa)'
sig_tests$`Tree mean`[4] = mean(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[4] = mean(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[4] = length(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree' & !is.na(TRY$K.mol.m.s.MPa))])
sig_tests$`n liana`[4] = length(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana' & !is.na(TRY$K.mol.m.s.MPa))])
sig_tests$`Test Statistic`[4] = wilcox.test(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')],
                                            TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[4] = wilcox.test(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')],
                                     TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for P50
sig_tests$Trait[5] = 'P<sub>50</sub> (MPa)'
sig_tests$`Tree mean`[5] = mean(TRY$P50.MPa[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[5] = mean(TRY$P50.MPa[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[5] = length(TRY$P50.MPa[which(TRY$Growth.form == 'tree' & !is.na(TRY$P50.MPa))])
sig_tests$`n liana`[5] = length(TRY$P50.MPa[which(TRY$Growth.form == 'liana' & !is.na(TRY$P50.MPa))])
sig_tests$`Test Statistic`[5] = wilcox.test(TRY$P50.MPa[which(TRY$Growth.form == 'tree')],
                                            TRY$P50.MPa[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[5] = wilcox.test(TRY$P50.MPa[which(TRY$Growth.form == 'tree')],
                                     TRY$P50.MPa[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for leaf lifespan
sig_tests$Trait[6] = 'Leaf lifespan (months)'
sig_tests$`Tree mean`[6] = mean(TRY$LL.months[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[6] = mean(TRY$LL.months[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[6] = length(TRY$LL.months[which(TRY$Growth.form == 'tree' & !is.na(TRY$LL.months))])
sig_tests$`n liana`[6] = length(TRY$LL.months[which(TRY$Growth.form == 'liana' & !is.na(TRY$LL.months))])
sig_tests$`Test Statistic`[6] = wilcox.test(TRY$LL.months[which(TRY$Growth.form == 'tree')],
                                            TRY$LL.months[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[6] = wilcox.test(TRY$LL.months[which(TRY$Growth.form == 'tree')],
                                     TRY$LL.months[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for SLA
sig_tests$Trait[7] = 'Specific leaf area (mm<sup>2</sup>/mg)'
sig_tests$`Tree mean`[7] = mean(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[7] = mean(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[7] = length(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree' & !is.na(TRY$SLA.mm2.mg))])
sig_tests$`n liana`[7] = length(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana' & !is.na(TRY$SLA.mm2.mg))])
sig_tests$`Test Statistic`[7] = wilcox.test(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')],
                                            TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[7] = wilcox.test(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')],
                                     TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for leaf N
sig_tests$Trait[8] = 'Area-based leaf nitrogen (g/m<sup>2</sup>)'
sig_tests$`Tree mean`[8] = mean(TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[8] = mean(TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[8] = length(TRY$Narea.g.m2[which(TRY$Growth.form == 'tree' & !is.na(TRY$Narea.g.m2))])
sig_tests$`n liana`[8] = length(TRY$Narea.g.m2[which(TRY$Growth.form == 'liana' & !is.na(TRY$Narea.g.m2))])
sig_tests$`Test Statistic`[8] = wilcox.test(TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')],
                                            TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[8] = wilcox.test(TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')],
                                     TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for A
sig_tests$Trait[9] = 'Area-based photosynthetic rate (mmol CO<sub>2</sub>/m<sup>2</sup>/s)'
sig_tests$`Tree mean`[9] = mean(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[9] = mean(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[9] = length(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree' & !is.na(TRY$Aarea.mmolC02.m2.s))])
sig_tests$`n liana`[9] = length(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana' & !is.na(TRY$Aarea.mmolC02.m2.s))])
sig_tests$`Test Statistic`[9] = wilcox.test(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')],
                                            TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[9] = wilcox.test(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')],
                                     TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for leaf N
sig_tests$Trait[10] = 'Mass-based leaf nitrogen (mg/g)'
sig_tests$`Tree mean`[10] = mean(TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[10] = mean(TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[10] = length(TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree' & !is.na(TRY$Nmass.mg.g))])
sig_tests$`n liana`[10] = length(TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana' & !is.na(TRY$Nmass.mg.g))])
sig_tests$`Test Statistic`[10] = wilcox.test(TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')],
                                             TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[10] = wilcox.test(TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')],
                                      TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for leaf P
sig_tests$Trait[11] = 'Mass-based leaf phosophorus (mg/g)'
sig_tests$`Tree mean`[11] = mean(TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[11] = mean(TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[11] = length(TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree' & !is.na(TRY$Pmass.mg.g))])
sig_tests$`n liana`[11] = length(TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana' & !is.na(TRY$Pmass.mg.g))])
sig_tests$`Test Statistic`[11] = wilcox.test(TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')],
                                             TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[11] = wilcox.test(TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')],
                                      TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for leaf area
sig_tests$Trait[12] = 'Leaf area (cm<sup>2</sup>)'
sig_tests$`Tree mean`[12] = mean(TRY$LA.cm2[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[12] = mean(TRY$LA.cm2[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[12] = length(TRY$LA.cm2[which(TRY$Growth.form == 'tree' & !is.na(TRY$LA.cm2))])
sig_tests$`n liana`[12] = length(TRY$LA.cm2[which(TRY$Growth.form == 'liana' & !is.na(TRY$LA.cm2))])
sig_tests$`Test Statistic`[12] = wilcox.test(TRY$LA.cm2[which(TRY$Growth.form == 'tree')], 
                                             TRY$LA.cm2[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[12] = wilcox.test(TRY$LA.cm2[which(TRY$Growth.form == 'tree')],
                                      TRY$LA.cm2[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for SRL
sig_tests$Trait[13] = 'Specific root length (m/g)'
sig_tests$`Tree mean`[13] = mean(TRY$SRL.m.g[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[13] = mean(TRY$SRL.m.g[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[13] = length(TRY$SRL.m.g[which(TRY$Growth.form == 'tree' & !is.na(TRY$SRL.m.g))])
sig_tests$`n liana`[13] = length(TRY$SRL.m.g[which(TRY$Growth.form == 'liana' & !is.na(TRY$SRL.m.g))])
sig_tests$`Test Statistic`[13] = wilcox.test(TRY$SRL.m.g[which(TRY$Growth.form == 'tree')],
                                             TRY$SRL.m.g[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[13] = wilcox.test(TRY$SRL.m.g[which(TRY$Growth.form == 'tree')],
                                      TRY$SRL.m.g[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for root diameter
sig_tests$Trait[14] = 'Fine root diameter (mm)'
sig_tests$`Tree mean`[14] = mean(TRY$FRD.mm[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[14] = mean(TRY$FRD.mm[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[14] = length(TRY$FRD.mm[which(TRY$Growth.form == 'tree' & !is.na(TRY$FRD.mm))])
sig_tests$`n liana`[14] = length(TRY$FRD.mm[which(TRY$Growth.form == 'liana' & !is.na(TRY$FRD.mm))])
sig_tests$`Test Statistic`[14] = wilcox.test(TRY$FRD.mm[which(TRY$Growth.form == 'tree')],
                                             TRY$FRD.mm[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[14] = wilcox.test(TRY$FRD.mm[which(TRY$Growth.form == 'tree')],
                                      TRY$FRD.mm[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for mycorrhizal colonization
sig_tests$Trait[15] = 'Mycorrhizal colonization (%)'
sig_tests$`Tree mean`[15] = mean(TRY$`MC.%`[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[15] = mean(TRY$`MC.%`[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[15] = length(TRY$`MC.%`[which(TRY$Growth.form == 'tree' & !is.na(TRY$`MC.%`))])
sig_tests$`n liana`[15] = length(TRY$`MC.%`[which(TRY$Growth.form == 'liana' & !is.na(TRY$`MC.%`))])
sig_tests$`Test Statistic`[15] = wilcox.test(TRY$`MC.%`[which(TRY$Growth.form == 'tree')],
                                             TRY$`MC.%`[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[15] = wilcox.test(TRY$`MC.%`[which(TRY$Growth.form == 'tree')],
                                      TRY$`MC.%`[which(TRY$Growth.form == 'liana')])$p.value

# Fill in table for rooting depth
sig_tests$Trait[16] = 'Rooting depth (m)'
sig_tests$`Tree mean`[16] = mean(TRY$RRD.m[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[16] = mean(TRY$RRD.m[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[16] = length(TRY$RRD.m[which(TRY$Growth.form == 'tree' & !is.na(TRY$RRD.m))])
sig_tests$`n liana`[16] = length(TRY$RRD.m[which(TRY$Growth.form == 'liana' & !is.na(TRY$RRD.m))])
sig_tests$`Test Statistic`[16] = wilcox.test(TRY$RRD.m[which(TRY$Growth.form == 'tree')],
                                             TRY$RRD.m[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[16] = wilcox.test(TRY$RRD.m[which(TRY$Growth.form == 'tree')],
                                      TRY$RRD.m[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$`p-value`
sig_tests$`p-value` = c('3.75 x 10<sup>-3</sup>',
                        '< 1.0 x 10<sup>-5</sup>',
                        '7.61 x 10<sup>-1</sup>',
                        '8.11 x 10<sup>-5</sup>',
                        '5.73 x 10<sup>-3</sup>',
                        '8.96 x 10<sup>-2</sup>',
                        '< 1.0 x 10<sup>-5</sup>',
                        '< 1.0 x 10<sup>-5</sup>',
                        '7.55 x 10<sup>-5</sup>',
                        '2.03 x 10<sup>-2</sup>',
                        '< 1.0 x 10<sup>-5</sup>',
                        '1.80 x 10<sup>-3</sup>',
                        '6.46 x 10<sup>-1</sup>',
                        '2.84 x 10<sup>-1</sup>',
                        '8.60 x 10<sup>-1</sup>',
                        '7.07 x 10<sup>-1</sup>')
sig_tests %>% 
  gt() %>%
  fmt_markdown(columns = vars('Trait', 'p-value')) %>%
  fmt_number(columns = vars('Tree mean', 'Liana mean'), decimals = 2) %>%
  fmt_number(columns = vars('n tree', 'n liana'), decimals = 0) %>%
  fmt_number(columns = vars('Test Statistic'), decimals = 0) %>%
  cols_label('Tree mean' = html('&mu;<sub><sub>tree</sub></sub>')) %>%
  cols_label('Liana mean' = html('&mu;<sub><sub>liana</sub></sub>')) %>%
  cols_label('n tree' = html('n<sub><sub>tree</sub></sub>')) %>%
  cols_label('n liana' = html('n<sub><sub>liana</sub></sub>')) %>%
  tab_header(
    title = md('Mann-Whitney U Tests for TRY traits')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(90),
    heading.title.font.size = 22) %>%
  tab_style(
    style = cell_fill(color = 'lightyellow'),
    locations = cells_body(rows = c(1:2, 4:5, 7:12))) %>%
  cols_align(align = 'center',
             columns = c('Tree mean', 'Liana mean', 'n tree', 'n liana', 'Test Statistic', 'p-value')) %>%
  cols_align(align = 'left',
             columns = 'Trait') %>%
  gtsave(filename = 'Plots/SupplementaryTable5.png')

###########################
## Table 2: Effect sizes ##
###########################

## gt() version of Supplementary Table 6

## Note that Glass' Delta is computed using either the liana
## or tree growth form as the reference group

# Storage & formatting
effect_size = as_tibble(matrix(, nrow = 16, ncol = 7))
colnames(effect_size) = c('Trait', 'Glass &Delta<sub>T</sub>', 'Lower CI<sub>T</sub>', 'Upper CI<sub>T</sub>', 'Glass &Delta<sub>L</sub>', 'Lower CI<sub>L</sub>', 'Upper CI<sub>L</sub>')

# Fill in table for wood density
effect_size$Trait[1] = 'Stem specific density (g/cm<sup>3</sup>)'
effect_size$`Glass &Delta<sub>T</sub>`[1] = glass_delta(TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')], 
                                                        TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[1] = glass_delta(TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')], 
                                                    TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[1] = glass_delta(TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')], 
                                                    TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[1] = glass_delta(TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')], 
                                                        TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[1] = glass_delta(TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')], 
                                                    TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[1] = glass_delta(TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')], 
                                                    TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for vessel diameter
effect_size$Trait[2] = 'Vessel diameter (&mu;m)'
effect_size$`Glass &Delta<sub>T</sub>`[2] = glass_delta(TRY$WVD.um[which(TRY$Growth.form == 'liana')], 
                                                        TRY$WVD.um[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[2] = glass_delta(TRY$WVD.um[which(TRY$Growth.form == 'liana')], 
                                                    TRY$WVD.um[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[2] = glass_delta(TRY$WVD.um[which(TRY$Growth.form == 'liana')], 
                                                    TRY$WVD.um[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[2] = glass_delta(TRY$WVD.um[which(TRY$Growth.form == 'tree')], 
                                                        TRY$WVD.um[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[2] = glass_delta(TRY$WVD.um[which(TRY$Growth.form == 'tree')], 
                                                    TRY$WVD.um[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[2] = glass_delta(TRY$WVD.um[which(TRY$Growth.form == 'tree')], 
                                                    TRY$WVD.um[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for vessel density
effect_size$Trait[3] = 'Vessel density (1/mm<sup>2</sup>)'
effect_size$`Glass &Delta<sub>T</sub>`[3] = glass_delta(TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')], 
                                                        TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[3] = glass_delta(TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')], 
                                                    TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[3] = glass_delta(TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')], 
                                                    TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[3] = glass_delta(TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')], 
                                                        TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[3] = glass_delta(TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')], 
                                                    TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[3] = glass_delta(TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')], 
                                                    TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for K
effect_size$Trait[4] = 'Stem specific hydraulic conductivity (mol/m/s/MPa)'
effect_size$`Glass &Delta<sub>T</sub>`[4] = glass_delta(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')], 
                                                        TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[4] = glass_delta(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')], 
                                                    TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[4] = glass_delta(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')], 
                                                    TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[4] = glass_delta(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')], 
                                                        TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[4] = glass_delta(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')], 
                                                    TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[4] = glass_delta(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')], 
                                                    TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for P50
effect_size$Trait[5] = 'P<sub>50</sub> (MPa)'
effect_size$`Glass &Delta<sub>T</sub>`[5] = glass_delta(TRY$P50.MPa[which(TRY$Growth.form == 'liana')], 
                                                        TRY$P50.MPa[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[5] = glass_delta(TRY$P50.MPa[which(TRY$Growth.form == 'liana')], 
                                                    TRY$P50.MPa[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[5] = glass_delta(TRY$P50.MPa[which(TRY$Growth.form == 'liana')], 
                                                    TRY$P50.MPa[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[5] = glass_delta(TRY$P50.MPa[which(TRY$Growth.form == 'tree')], 
                                                        TRY$P50.MPa[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[5] = glass_delta(TRY$P50.MPa[which(TRY$Growth.form == 'tree')], 
                                                    TRY$P50.MPa[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[5] = glass_delta(TRY$P50.MPa[which(TRY$Growth.form == 'tree')], 
                                                    TRY$P50.MPa[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for leaf lifespan
effect_size$Trait[6] = 'Leaf lifespan (months)'
effect_size$`Glass &Delta<sub>T</sub>`[6] = glass_delta(TRY$LL.months[which(TRY$Growth.form == 'liana')], 
                                                        TRY$LL.months[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[6] = glass_delta(TRY$LL.months[which(TRY$Growth.form == 'liana')], 
                                                    TRY$LL.months[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[6] = glass_delta(TRY$LL.months[which(TRY$Growth.form == 'liana')], 
                                                    TRY$LL.months[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[6] = glass_delta(TRY$LL.months[which(TRY$Growth.form == 'tree')], 
                                                        TRY$LL.months[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[6] = glass_delta(TRY$LL.months[which(TRY$Growth.form == 'tree')], 
                                                    TRY$LL.months[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[6] = glass_delta(TRY$LL.months[which(TRY$Growth.form == 'tree')], 
                                                    TRY$LL.months[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for SLA
effect_size$Trait[7] = 'Specific leaf area (mm<sup>2</sup>/mg)'
effect_size$`Glass &Delta<sub>T</sub>`[7] = glass_delta(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')], 
                                                        TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[7] = glass_delta(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')], 
                                                    TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[7] = glass_delta(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')], 
                                                    TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[7] = glass_delta(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')], 
                                                        TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[7] = glass_delta(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')], 
                                                    TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[7] = glass_delta(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')], 
                                                    TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for leaf N
effect_size$Trait[8] = 'Area-based leaf nitrogen (g/m<sup>2</sup>)'
effect_size$`Glass &Delta<sub>T</sub>`[8] = glass_delta(TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')], 
                                                        TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[8] = glass_delta(TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')], 
                                                    TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[8] = glass_delta(TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')], 
                                                    TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[8] = glass_delta(TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')], 
                                                        TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[8] = glass_delta(TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')], 
                                                    TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[8] = glass_delta(TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')], 
                                                    TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for A
effect_size$Trait[9] = 'Area-based photosynthetic rate (mmol CO<sub>2</sub>/m<sup>2</sup>/s)'
effect_size$`Glass &Delta<sub>T</sub>`[9] = glass_delta(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')], 
                                                        TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[9] = glass_delta(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')], 
                                                    TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[9] = glass_delta(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')], 
                                                    TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[9] = glass_delta(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')], 
                                                        TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[9] = glass_delta(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')], 
                                                    TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[9] = glass_delta(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')], 
                                                    TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for leaf N
effect_size$Trait[10] = 'Mass-based leaf nitrogen (mg/g)'
effect_size$`Glass &Delta<sub>T</sub>`[10] = glass_delta(TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')], 
                                                         TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[10] = glass_delta(TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')], 
                                                     TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[10] = glass_delta(TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')], 
                                                     TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[10] = glass_delta(TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')], 
                                                         TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[10] = glass_delta(TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')], 
                                                     TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[10] = glass_delta(TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')], 
                                                     TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for leaf P
effect_size$Trait[11] = 'Mass-based leaf phosphorus (mg/g)'
effect_size$`Glass &Delta<sub>T</sub>`[11] = glass_delta(TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')], 
                                                         TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[11] = glass_delta(TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')], 
                                                     TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[11] = glass_delta(TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')],
                                                     TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[11] = glass_delta(TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')], 
                                                         TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[11] = glass_delta(TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')], 
                                                     TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[11] = glass_delta(TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')], 
                                                     TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for leaf area
effect_size$Trait[12] = 'Leaf area (cm<sup>2</sup>)'
effect_size$`Glass &Delta<sub>T</sub>`[12] = glass_delta(TRY$LA.cm2[which(TRY$Growth.form == 'liana')], 
                                                         TRY$LA.cm2[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[12] = glass_delta(TRY$LA.cm2[which(TRY$Growth.form == 'liana')], 
                                                     TRY$LA.cm2[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[12] = glass_delta(TRY$LA.cm2[which(TRY$Growth.form == 'liana')], 
                                                     TRY$LA.cm2[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[12] = glass_delta(TRY$LA.cm2[which(TRY$Growth.form == 'tree')], 
                                                         TRY$LA.cm2[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[12] = glass_delta(TRY$LA.cm2[which(TRY$Growth.form == 'tree')], 
                                                     TRY$LA.cm2[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[12] = glass_delta(TRY$LA.cm2[which(TRY$Growth.form == 'tree')], 
                                                     TRY$LA.cm2[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for SRL
effect_size$Trait[13] = 'Specific root length (m/g)'
effect_size$`Glass &Delta<sub>T</sub>`[13] = glass_delta(TRY$SRL.m.g[which(TRY$Growth.form == 'liana')], 
                                                         TRY$SRL.m.g[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[13] = glass_delta(TRY$SRL.m.g[which(TRY$Growth.form == 'liana')], 
                                                     TRY$SRL.m.g[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[13] = glass_delta(TRY$SRL.m.g[which(TRY$Growth.form == 'liana')], 
                                                     TRY$SRL.m.g[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[13] = glass_delta(TRY$SRL.m.g[which(TRY$Growth.form == 'tree')], 
                                                         TRY$SRL.m.g[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[13] = glass_delta(TRY$SRL.m.g[which(TRY$Growth.form == 'tree')], 
                                                     TRY$SRL.m.g[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[13] = glass_delta(TRY$SRL.m.g[which(TRY$Growth.form == 'tree')], 
                                                     TRY$SRL.m.g[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for root diameter
effect_size$Trait[14] = 'Fine root diameter (mm)'
effect_size$`Glass &Delta<sub>T</sub>`[14] = glass_delta(TRY$FRD.mm[which(TRY$Growth.form == 'liana')], 
                                                         TRY$FRD.mm[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[14] = glass_delta(TRY$FRD.mm[which(TRY$Growth.form == 'liana')], 
                                                     TRY$FRD.mm[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[14] = glass_delta(TRY$FRD.mm[which(TRY$Growth.form == 'liana')], 
                                                     TRY$FRD.mm[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[14] = glass_delta(TRY$FRD.mm[which(TRY$Growth.form == 'tree')], 
                                                         TRY$FRD.mm[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[14] = glass_delta(TRY$FRD.mm[which(TRY$Growth.form == 'tree')], 
                                                     TRY$FRD.mm[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[14] = glass_delta(TRY$FRD.mm[which(TRY$Growth.form == 'tree')], 
                                                     TRY$FRD.mm[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for mycorrhizal colonization
effect_size$Trait[15] = 'Mycorrhizal colonization (%)'
effect_size$`Glass &Delta<sub>T</sub>`[15] = glass_delta(TRY$`MC.%`[which(TRY$Growth.form == 'liana')], 
                                                         TRY$`MC.%`[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[15] = glass_delta(TRY$`MC.%`[which(TRY$Growth.form == 'liana')], 
                                                     TRY$`MC.%`[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[15] = glass_delta(TRY$`MC.%`[which(TRY$Growth.form == 'liana')], 
                                                     TRY$`MC.%`[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[15] = glass_delta(TRY$`MC.%`[which(TRY$Growth.form == 'tree')], 
                                                         TRY$`MC.%`[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[15] = glass_delta(TRY$`MC.%`[which(TRY$Growth.form == 'tree')], 
                                                     TRY$`MC.%`[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[15] = glass_delta(TRY$`MC.%`[which(TRY$Growth.form == 'tree')], 
                                                     TRY$`MC.%`[which(TRY$Growth.form == 'liana')])$CI_high

# Fill in table for rooting depth
effect_size$Trait[16] = 'Rooting depth (m)'
effect_size$`Glass &Delta<sub>T</sub>`[16] = glass_delta(TRY$RRD.m[which(TRY$Growth.form == 'liana')], 
                                                         TRY$RRD.m[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI<sub>T</sub>`[16] = glass_delta(TRY$RRD.m[which(TRY$Growth.form == 'liana')], 
                                                     TRY$RRD.m[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI<sub>T</sub>`[16] = glass_delta(TRY$RRD.m[which(TRY$Growth.form == 'liana')], 
                                                     TRY$RRD.m[which(TRY$Growth.form == 'tree')])$CI_high
effect_size$`Glass &Delta<sub>L</sub>`[16] = glass_delta(TRY$RRD.m[which(TRY$Growth.form == 'tree')], 
                                                         TRY$RRD.m[which(TRY$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI<sub>L</sub>`[16] = glass_delta(TRY$RRD.m[which(TRY$Growth.form == 'tree')], 
                                                     TRY$RRD.m[which(TRY$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI<sub>L</sub>`[16] = glass_delta(TRY$RRD.m[which(TRY$Growth.form == 'tree')], 
                                                     TRY$RRD.m[which(TRY$Growth.form == 'liana')])$CI_high

effect_size %>% gt() %>%
  fmt_markdown(columns = c('Trait')) %>%
  fmt_number(columns = c('Glass &Delta<sub>T</sub>', 'Lower CI<sub>T</sub>', 'Upper CI<sub>T</sub>',
                         'Glass &Delta<sub>L</sub>', 'Lower CI<sub>L</sub>', 'Upper CI<sub>L</sub>'), decimals = 2) %>%
  cols_label('Glass &Delta<sub>T</sub>' = html("Glass' &Delta;<sub>T</sub>")) %>%
  cols_label('Lower CI<sub>T</sub>' = html('Lower CI<sub>T</sub>')) %>%
  cols_label('Upper CI<sub>T</sub>' = html('Upper CI<sub>T</sub>')) %>%
  cols_label('Glass &Delta<sub>L</sub>' = html("Glass' &Delta;<sub>L</sub>")) %>%
  cols_label('Lower CI<sub>L</sub>' = html('Lower CI<sub>L</sub>')) %>%
  cols_label('Upper CI<sub>L</sub>' = html('Upper CI<sub>L</sub>')) %>%
  tab_header(
    title = md('Effect Size for TRY traits')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(90),
    heading.title.font.size = 22) %>%
  tab_style(
    style = cell_fill(color = 'lightyellow'),
    locations = cells_body(rows = c(1:2, 4:6, 7:9, 11:12))) %>%
  cols_align(align = 'center',
             columns = c('Glass &Delta<sub>T</sub>', 'Lower CI<sub>T</sub>', 'Upper CI<sub>T</sub>',
                         'Glass &Delta<sub>L</sub>', 'Lower CI<sub>L</sub>', 'Upper CI<sub>L</sub>')) %>%
  cols_align(align = 'left',
             columns = 'Trait') %>%
  gtsave(filename = 'Plots/SupplementaryTable6.png')

########################################
#### Part 6: Extended Meta-analysis ####
########################################

## Supplementary Tables 1 & 2

rm(list = ls())

data = read.csv('full_met_analysis_data_8_Feb.csv')

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

## gt() version of Supplementary Table 1

# Storage & formatting
sig_tests = as_tibble(matrix(, nrow = 3, ncol = 7))
colnames(sig_tests) = c('Trait', 'Tree mean', 'Liana mean', 'n tree', 'n liana', 'Test Statistic', 'p-value')

# Fill in table for K
sig_tests$Trait[1] = 'Stem-specific hydraulic conductivity (mol/m/s/MPa)'
sig_tests$`Tree mean`[1] = mean(data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[1] = mean(data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[1] = length(data$K.mol.m.s.MPa[which(data$Growth.form == 'tree' & !is.na(data$K.mol.m.s.MPa))])
sig_tests$`n liana`[1] = length(data$K.mol.m.s.MPa[which(data$Growth.form == 'liana' & !is.na(data$K.mol.m.s.MPa))])
sig_tests$`Test Statistic`[1] = wilcox.test(data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')],  
                                            data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[1] = wilcox.test(data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')], 
                                     data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')])$p.value

# Fill in table for P50
sig_tests$Trait[2] = 'P<sub>50</sub> (MPa)'
sig_tests$`Tree mean`[2] = mean(data$P50[which(data$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[2] = mean(data$P50[which(data$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[2] = length(data$P50[which(data$Growth.form == 'tree' & !is.na(data$P50))])
sig_tests$`n liana`[2] = length(data$P50[which(data$Growth.form == 'liana' & !is.na(data$P50))])
sig_tests$`Test Statistic`[2] = wilcox.test(data$P50[which(data$Growth.form == 'tree')],
                                            data$P50[which(data$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[2] = wilcox.test(data$P50[which(data$Growth.form == 'tree')],
                                     data$P50[which(data$Growth.form == 'liana')])$p.value

# Fill in table for Slope
sig_tests$Trait[3] = 'Slope of PLC curve (%/MPa)'
sig_tests$`Tree mean`[3] = mean(data$Slope[which(data$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[3] = mean(data$Slope[which(data$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[3] = length(data$Slope[which(data$Growth.form == 'tree' & !is.na(data$Slope))])
sig_tests$`n liana`[3] = length(data$Slope[which(data$Growth.form == 'liana' & !is.na(data$Slope))])
sig_tests$`Test Statistic`[3] = wilcox.test(data$Slope[which(data$Growth.form == 'tree')],
                                            data$Slope[which(data$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[3] = wilcox.test(data$Slope[which(data$Growth.form == 'tree')],
                                     data$Slope[which(data$Growth.form == 'liana')])$p.value

# Adjust p-value formatting
sig_tests$`p-value`
sig_tests$`p-value` = c('< 1.0 x 10<sup>-5</sup>',
                        '1.29 x 10<sup>-1</sup>',
                        '1.85 x 10<sup>-1</sup>')
sig_tests %>% gt() %>%
  fmt_markdown(columns = vars('Trait', 'p-value')) %>%
  fmt_number(columns = vars('Tree mean', 'Liana mean'), decimals = 2) %>%
  fmt_number(columns = vars('n tree', 'n liana'), decimals = 0) %>%
  fmt_number(columns = vars('Test Statistic'), decimals = 0) %>%
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
  gtsave(filename = 'Plots/SupplementaryTable1.png')

##########################
## Table 2: Effect size ##
##########################

## gt() version of Supplementary Table 2

# Storage & formatting
effect_size = as_tibble(matrix(, nrow = 3, ncol = 7))
colnames(effect_size) = c('Trait', 
                          'Glass Delta_T', 'Lower CI_T', 'Upper CI_T', 
                          'Glass Delta_L', 'Lower CI_L', 'Upper CI_L')

# Fill in table for K
effect_size$Trait[1] = 'Stem-specific hydraulic conductivity (mol/m/s/MPa)'
effect_size$`Glass Delta_T`[1] = glass_delta(data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')], 
                                             data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI_T`[1] = glass_delta(data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')], 
                                          data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI_T`[1] = glass_delta(data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')], 
                                          data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')])$CI_high
effect_size$`Glass Delta_L`[1] = glass_delta(data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')], 
                                             data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI_L`[1] = glass_delta(data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')], 
                                          data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI_L`[1] = glass_delta(data$K.mol.m.s.MPa[which(data$Growth.form == 'tree')], 
                                          data$K.mol.m.s.MPa[which(data$Growth.form == 'liana')])$CI_high

# Fill in table for P50
effect_size$Trait[2] = 'P<sub>50</sub> (MPa)'
effect_size$`Glass Delta_T`[2] = glass_delta(data$P50[which(data$Growth.form == 'liana')], 
                                             data$P50[which(data$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI_T`[2] = glass_delta(data$P50[which(data$Growth.form == 'liana')], 
                                          data$P50[which(data$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI_T`[2] = glass_delta(data$P50[which(data$Growth.form == 'liana')], 
                                          data$P50[which(data$Growth.form == 'tree')])$CI_high
effect_size$`Glass Delta_L`[2] = glass_delta(data$P50[which(data$Growth.form == 'tree')], 
                                             data$P50[which(data$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI_L`[2] = glass_delta(data$P50[which(data$Growth.form == 'tree')], 
                                          data$P50[which(data$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI_L`[2] = glass_delta(data$P50[which(data$Growth.form == 'tree')], 
                                          data$P50[which(data$Growth.form == 'liana')])$CI_high

# Fill in table for Slope
effect_size$Trait[3] = 'Slope of PLC curve (%/MPa)'
effect_size$`Glass Delta_T`[3] = glass_delta(data$Slope[which(data$Growth.form == 'liana')], 
                                             data$Slope[which(data$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI_T`[3] = glass_delta(data$Slope[which(data$Growth.form == 'liana')], 
                                          data$Slope[which(data$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI_T`[3] = glass_delta(data$Slope[which(data$Growth.form == 'liana')], 
                                          data$Slope[which(data$Growth.form == 'tree')])$CI_high
effect_size$`Glass Delta_L`[3] = glass_delta(data$Slope[which(data$Growth.form == 'tree')], 
                                             data$Slope[which(data$Growth.form == 'liana')])$Glass_delta
effect_size$`Lower CI_L`[3] = glass_delta(data$Slope[which(data$Growth.form == 'tree')], 
                                          data$Slope[which(data$Growth.form == 'liana')])$CI_low
effect_size$`Upper CI_L`[3] = glass_delta(data$Slope[which(data$Growth.form == 'tree')], 
                                          data$Slope[which(data$Growth.form == 'liana')])$CI_high

effect_size %>% gt() %>%
  fmt_markdown(columns = c('Trait')) %>%
  fmt_number(columns = c('Glass Delta_T', 'Lower CI_T', 'Upper CI_T',
                         'Glass Delta_L', 'Lower CI_L', 'Upper CI_L'), decimals = 2) %>%
  cols_label('Glass Delta_T' = html("Glass' &Delta;<sub>T</sub>")) %>%
  cols_label('Lower CI_T' = html('Lower CI<sub>T</sub>')) %>%
  cols_label('Upper CI_T' = html('Upper CI<sub>T</sub>')) %>%
  cols_label('Glass Delta_L' = html("Glass' &Delta;<sub>L</sub>")) %>%
  cols_label('Lower CI_L' = html('Lower CI<sub>L</sub>')) %>%
  cols_label('Upper CI_L' = html('Upper CI<sub>L</sub>')) %>%
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
             columns = c('Glass Delta_T', 'Lower CI_T', 'Upper CI_T',
                         'Glass Delta_L', 'Lower CI_L', 'Upper CI_L')) %>%
  cols_align(align = 'left',
             columns = c('Trait')) %>%
  gtsave(filename = 'Plots/SupplementaryTable2.png')

######################
#### Part 7: Kreq ####
######################

## Supplementary Figure 1

## Before running this section, you must compile tree.NPP.mxh
## & liana.NPP.mxh from NPP_models.R

rm(list = ls())

# Load parameters
load('param.input.RData')

# Load met drivers
load('bci_met_mxh.RData')

# Interpolate K on log scale and back-transform
log_k = c(log(trees_K$K / 10), log(lianas_K$K))
log_varyk = seq(min(log_k), max(log_k), length.out = 500)
varyk = exp(log_varyk)

# Interpolate DBH (and remove highest DBH)
hori.data.liana.dbh = hori.data.liana.dbh[hori.data.liana.dbh < 29]
dbhs = c(2, 6)

# Storage
lianaouteach = c()
liana_out_mat = matrix(, nrow = length(dbhs), ncol = length(varyk))

for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 2 * 100
      frac.liana.al = 0.4
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * mean.hori.data.tree.dbh^tree.b2Ht) * 1.25
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      liana_out_mat[i,k] = liana.out[1]
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

# Format
rownames(liana_out_mat) = dbhs
colnames(liana_out_mat) = varyk

# Format
liana_melt = melt(liana_out_mat)
colnames(liana_melt) = c('DBH', 'Ks', 'NPP')
liana_melt$Ks = liana_melt$Ks * 0.001

# Plot reference scenarios
plot = liana_melt %>%
  ggplot(aes(x = Ks, y = NPP, group = as.factor(DBH), color = as.factor(DBH))) +
  geom_line(size = 1.2) +
  xlab(expression(paste(K[w],' (mol ', m^-1, ' ', s^-1, ' MP', a^-1,')'))) + ylab(expression(paste('NPP (kg C yea', r^-1,')'))) +
  geom_vline(aes(xintercept = liana_melt$Ks[which(liana_melt$DBH == 2) & liana_melt$NPP == min(abs(liana_melt$NPP[which(liana_melt$DBH == 2)]))], linetype = 'Kreq')) +
  geom_vline(aes(xintercept = liana_melt$Ks[which(liana_melt$DBH == 6) & liana_melt$NPP == min(abs(liana_melt$NPP[which(liana_melt$DBH == 6)]))], linetype = 'Kreq')) +
  scale_linetype_manual(name = element_blank(), values = c('Kreq' = 'dashed'), labels = expression(K[req])) +
  scale_color_npg(name = 'DBH (cm)') +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 25), legend.text = element_text(size = 25), legend.title = element_text(size = 28)) +
  guides(color = guide_legend(reverse = T)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1, 1, 1.2, 1), 'cm'))

# Inset zoomed in version
inset = liana_melt %>%
  ggplot(aes(x = Ks, y = NPP, group = as.factor(DBH), color = as.factor(DBH))) +
  geom_line(size = 1.2, show.legend = F) +
  xlab('') + ylab('') +
  geom_vline(aes(xintercept = liana_melt$Ks[which(liana_melt$DBH == 2) & liana_melt$NPP == min(abs(liana_melt$NPP[which(liana_melt$DBH == 2)]))]), linetype = 'dashed', show.legend = F) +
  geom_vline(aes(xintercept = liana_melt$Ks[which(liana_melt$DBH == 6) & liana_melt$NPP == min(abs(liana_melt$NPP[which(liana_melt$DBH == 6)]))]), linetype = 'dashed', show.legend = F) +
  scale_color_npg() +
  theme(axis.title = element_text(size = 22), axis.text = element_text(size = 22), legend.text = element_text(size = 22), legend.title = element_text(size = 22)) +
  coord_cartesian(xlim = c(0, 500), ylim = c(-5, 5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Create plot with inset
plot.w.inset = ggdraw() +
  draw_plot(plot) +
  draw_plot(inset, x = 0.35, y = 0.25, width = 0.4, height = 0.4) +
  theme(axis.text = element_text(size = 18))

ggsave(plot = plot.w.inset, filename = 'Plots/SupplementaryFigure1.jpeg', width = 12, height = 7, units = 'in')

####################################################
#### Part 8: Invading liana climate sensitivity ####
####################################################

## Supplementary Figures 10

## Before running this section, you must compile tree.NPP.mxh
## & liana.NPP.mxh from NPP_models.R

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Interpolation_mxh_100.RData')

# Indices
nsite = 100
nmonth = 12
nK = 75000

# Number of days in each month
nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Rename matrices
VPDs = VPD_interp_100
SWPs = SWP_interp_100

# Possible values of K
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

# Storage
tree.req = matrix(, nrow = nsite, ncol = nsite)

# Run model for invading liana scenario
for(vsite in 1:nsite){
  for(ssite in 1:nsite){
    sens = c()
    for(k in 1:nK){
      for(mo in 1:nmonth){
        nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        month = mo
        
        tree.dbh = mean.hori.data.tree.dbh
        
        tot.al = 2 * 100
        frac.tree.al = 0.9
        
        psis = SWPs[mo,ssite]
        
        D = VPDs[mo,,vsite]
        
        if(mo == 1){
          tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
        }else{
          tree.height = tree.out[5]
        }
        K = Ks[k]
        tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
        sens[mo] = tree.out[1]
      }
      sum = sum(sens)
      if(sum > 0){
        tree.req[vsite,ssite] = K
        break
      }
    }
    print(paste('Passing sites',vsite,'&',ssite))
  }
}

# First time running, save output (long computation time)
#save(tree.req, file = 'tree.req.100interp.monthly_invasion_200.RData')

# After first time, load output for further formatting
load(file = 'tree.req.100interp.monthly_invasion_200.RData')
# Format
tree.req = tree.req * 0.001
tree_req_melt = melt(tree.req)
colnames(tree_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

# Contour plot
ra1 = ggplot(tree_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv, fill = Ksurv)) +
  geom_raster(interpolate = T, show.legend = F) +
  geom_contour(alpha = 0, aes(color = ..level..), breaks = seq(2, 12, by = 1), show.legend = T) +
  geom_contour(alpha = 0.65, color = 'black', breaks = seq(2, 12, by = 1), show.legend = F) +
  xlab('VPD index') + ylab(expression(paste(Psi,' index'))) +
  scale_color_distiller(name = expression(K[req]), palette = 'PuBu', breaks = seq(4, 12, by = 1)) +
  scale_fill_distiller(name = expression(K[req]), palette = 'PuBu', breaks = seq(4, 12, by = 1)) +
  ggtitle('Tree') +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank()) +
  geom_dl(aes(label = ..level..), stat = 'contour', breaks = seq(2, 12, by = 1), method = list('top.pieces', color = 'black', cex = 1.8, vjust = -0.4))
ra1

# Truncate possible K values to reduce computation time for lianas
# Cut-off was determined by running the model for the wettest site and
# choosing the lowest value of K = Kreq within our simulations
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)
Ks = Ks[Ks > 17550]

# Storage
liana.req = matrix(, nrow = nsite, ncol = nsite)

# Run model for invading liana scenario
for(vsite in 1:nsite){
  for(ssite in 1:nsite){
    sens = c()
    for(k in 1:nK){
      for(mo in 1:12){
        nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        month = mo
        
        dbh = 2
        
        tot.al = 2 * 100
        frac.liana.al = 0.1
        
        psis = SWPs[mo,ssite]
        
        D = VPDs[mo,,vsite]
        
        if(mo == 1){
          liana.length = tree.out[2]
        }else{
          liana.length = lianaout[4]
        }
        K = Ks[k]
        lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
        sens[mo] = lianaout[1]
      }
      sum = sum(sens)
      if(sum > 0){
        liana.req[vsite, ssite] = K
        break
      }
    }
    print(paste('Passing sites',vsite,'&',ssite))
  }
}

# Save and load as described above
#save(liana.req, file = 'liana.req.100interp.monthly_invasion_200.RData')

load(file = 'liana.req.100interp.monthly_invasion_200.RData')
# Format
liana.req = liana.req * 0.001
liana_req_melt = melt(liana.req)
colnames(liana_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

# Countour plot
ra2 = ggplot(liana_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv, fill = Ksurv)) +
  geom_raster(interpolate = T, show.legend = T) +
  geom_contour(alpha = 0, aes(color = ..level..), show.legend = F, breaks = seq(10, 90, by = 10)) +
  geom_contour(alpha = 0.65, color = 'black', breaks = seq(10, 90, by = 10), show.legend = F) +
  xlab('VPD index') + ylab(expression(paste(Psi,' index'))) +
  scale_fill_distiller(name = expression(K[req]), palette = 'BuGn', breaks = seq(10, 90, by = 10)) +
  ggtitle('Liana') +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank()) +
  geom_dl(aes(label = ..level..), stat = 'contour', breaks = seq(10, 90, by = 10), method = list('top.pieces', color = 'black', cex = 1.8, vjust = -0.4))
ra2

# Plot together
ga_fin = plot_grid(ra1, ra2, 
                   nrow = 1, rel_widths = 1,
                   labels = c('A', 'B'),
                   label_size = 30)

ggsave(ga_fin, filename = 'Plots/SupplementaryFigure10.jpeg', height = 12, width = 20, units = 'in')

##################################
#### Part 9: Future scenarios ####
##################################

## Supplementary Figure 2

## Before running this scenario, you must compile tree.NPP.mxh.ca 
## & liana.NPP.mxh.ca from NPP_models_wCO2.R

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('bci_met_mxh.RData')
BCI_VPD = VPD
BCI_SWP = SWP
load('horizontes_met_mxh.RData')
H_VPD = VPD
H_SWP = SWP
rm(VPD, SWP)

# Indexing months and days
nmonth = 12
nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# More indexing
vpdsite = 11
nK = 75000

# Sequence of multipliers for VPD
# (from 1x present day to 2x present day)
mult = seq(1, 2, by = 0.1)

# BCI future scenarios
BCI_VPD_future = array(, dim = c(12, 24, 11))
BCI_VPD_future[,,1] = BCI_VPD

for(i in 2:11){
  BCI_VPD_future[,,i] = BCI_VPD * mult[i]
}

# Horizontes future scenarios
H_VPD_future = array(, dim = c(12, 24, 11))
H_VPD_future[,,1] = H_VPD

for(i in 2:11){
  H_VPD_future[,,i] = H_VPD * mult[i]
}

#####################################
## Tree runs, established scenario ##
#####################################

# Make input K values
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI

# Storage
tree.req_bci = c()

count = 0

# Run future established liana scenario at BCI (tropical moist forest)
for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:nmonth){
      # Specify Ca = 550
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 2 * 100
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD_future[mo,,vsite]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
      sens[mo] = tree.out[1]
      count = count + 1
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req_bci[vsite] = K
      break
    }
  }
  print(vsite)
}

tree.req_bci = tree.req_bci * 0.001

## Horizontes

# Storage
tree.req_h = c()

count = 0

# Run future established liana scenario at Horizontes (tropical dry forest)
for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:nmonth){
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 2 * 100
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD_future[mo,,vsite]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
      sens[mo] = tree.out[1]
      count = count + 1
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req_h[vsite] = K
      break
    }
  }
  print(vsite)
}

tree.req_h = tree.req_h * 0.001

pl1 = ggplot() +
  geom_point(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Tree') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pl1

######################################
## Liana runs, established scenario ##
######################################

# Possible K values
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI

# Storage
liana.req_bci = c()

count = 0

# Run future established liana scenario at BCI (tropical moist forest)
for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      # Specify Ca = 550
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 2 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD_future[mo,,vsite]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
      sens[mo] = lianaout[1]
      count = count + 1
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req_bci[vsite] = K
      break
    }
  }
  print(vsite)
}

liana.req_bci = liana.req_bci * 0.001

## Horizontes

# Storage
liana.req_h = c()

count = 0

# Run future established liana scenario at Horizontes (tropical dry forest)
for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      # Specify Ca = 550
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 2 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD_future[mo,,vsite]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
      sens[mo] = lianaout[1]
      count = count + 1
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req_h[vsite] = K
      break
    }
  }
  print(vsite)
}

liana.req_h = liana.req_h * 0.001

pl2 = ggplot() +
  geom_point(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Liana') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pl2

# Function to extract legend from ggplot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Extract legend
legend = get_legend(pl2)

# Plot together
ga = plot_grid(pl1 + theme(legend.position = 'none'), 
               pl2 + theme(legend.position = 'none'), 
               legend, nrow = 1, 
               rel_widths = c(2.3, 2.3, 0.7),
               labels = c('A', 'B', ''),
               label_size = 30, hjust = -1)

ggsave(ga, filename = 'Plots/SupplementaryFigure2a.jpeg', width = 16, height = 8, units = 'in')

##################################
## Tree runs, invasion scenario ##
##################################

# Input K values
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI

# Storage
tree.req_bci = c()

count = 0

# Run future invading liana scenario for BCI (tropical moist forest)
for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:nmonth){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 2 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD_future[mo,,vsite]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
      sens[mo] = tree.out[1]
      count = count + 1
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req_bci[vsite] = K
      break
    }
  }
  print(vsite)
}

tree.req_bci = tree.req_bci * 0.001

## Horizontes

# Storage
tree.req_h = c()

count = 0

# Run future invading liana scenario at Horizontes (tropical dry forest)
for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:nmonth){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 2 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD_future[mo,,vsite]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
      sens[mo] = tree.out[1]
      count = count + 1
    }
    sum = sum(sens)
    if(sum > 0){
      tree.req_h[vsite] = K
      break
    }
  }
  print(vsite)
}

tree.req_h = tree.req_h * 0.001

pl1 = ggplot() +
  geom_point(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Tree') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pl1

###################################
## Liana runs, invasion scenario ##
###################################

# Input K values
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI

# Storage
liana.req_bci = c()

count = 0

# Run future invading liana scenario at BCI (tropical moist forest)
for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 2 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD_future[mo,,vsite]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
      sens[mo] = lianaout[1]
      count = count + 1
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req_bci[vsite] = K
      break
    }
  }
  print(vsite)
}

liana.req_bci = liana.req_bci * 0.001

## Horizontes

# Storage
liana.req_h = c()

count = 0

# Run future invading liana scenario at Horizontes (tropical dry forest)
for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 2 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD_future[mo,,vsite]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
      sens[mo] = lianaout[1]
      count = count + 1
    }
    sum = sum(sens)
    if(sum > 0){
      liana.req_h[vsite] = K
      break
    }
  }
  print(vsite)
}

liana.req_h = liana.req_h * 0.001

pl2 = ggplot() +
  geom_point(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Liana') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pl2

# Get legend from function above
legend = get_legend(pl1)

# Plot together
ga = plot_grid(pl1 + theme(legend.position = 'none'), 
               pl2 + theme(legend.position = 'none'),
               legend, nrow = 1, 
               rel_widths = c(2.3, 2.3, 0.7),
               labels = c('C', 'D', ''),
               label_size = 30, hjust = -1.5)

ggsave(ga, filename = 'Plots/SupplementaryFigure2b.jpeg', width = 16, height = 8, units = 'in')

##############################
#### Part 10: future Kreq ####
##############################

## Supplementary Figure 3

## Before running this section, you must compile
## tree.NPP.mxh.ca & liana.NPP.mxh.ca from NPP_models_wCO2.R

## These simulations correspond to the invading liana scenario
## in the present and future (2100)

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Driver_calc/Horizontes/horizontes_met_mxh.RData')
H_VPD = VPD
H_SWP = SWP
rm(VPD, SWP)

# Indexing months and days
nmonth = 12
nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

nK = 75000

# Horizontes VPD future scenarios
H_VPD_future = array(, dim = c(12, 24, 2))
H_VPD_future[,,1] = H_VPD
H_VPD_future[,,2] = H_VPD * 2

#### Tree simulation - present day

# Input K values
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

# Storage
tree.req_h = c()
sens = c()

for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 400
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 2 * 100
    frac.tree.al = 0.9
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,1]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h[1] = K
    break
  }
}

#### Liana simulation -- present day

# Storage
liana.req_h = c()
sens = c()

for(k in 1:nK){
  for(mo in 1:12){
    Ca = 400
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = 2
    
    tot.al = 2 * 100
    frac.liana.al = 0.1
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,1]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h[1] = K
    break
  }
}

#### Tree simulation -- future scenario

# Storage
sens = c()

for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 2 * 100
    frac.tree.al = 0.9
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,2]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h[2] = K
    break
  }
}

#### Liana simulation -- future scenario

# Storage
sens = c()

for(k in 1:nK){
  for(mo in 1:12){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = 2
    
    tot.al = 2 * 100
    frac.liana.al = 0.1
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,2]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h[2] = K
    break
  }
}

tree.req_h = tree.req_h * 0.001
liana.req_h = liana.req_h * 0.001

# Make workable dataframes
liana_kreq = cbind(c('2000', '2100'), 
                   c(liana.req_h[1], liana.req_h[2]))
colnames(liana_kreq) = c('Time', 'Kreq')

tree_kreq = cbind(c('2000', '2100'),
                  c(tree.req_h[1], tree.req_h[2]))
colnames(tree_kreq) = c('Time', 'Kreq')

# Combine and format dataframes
comb_kreq = rbind(liana_kreq, tree_kreq)
comb_kreq = as.data.frame(comb_kreq)
comb_kreq$growth.form = c('Liana', 'Liana', 
                          'Tree', 'Tree')
comb_kreq$Time = as.factor(comb_kreq$Time)
comb_kreq$Kreq = as.numeric(comb_kreq$Kreq)
comb_kreq$growth.form = as.factor(comb_kreq$growth.form)

# Plot
ggplot(comb_kreq, aes(x = Time, y = Kreq, color = growth.form, group = growth.form)) +
  geom_point(size = 4, show.legend = F) +
  geom_line(size = 1.5, show.legend = F, alpha = 0.7) +
  scale_x_discrete(limits = c('2000', '2100')) +
  xlab('') + ylab(expression(paste(K[req], ' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  scale_color_npg(name = 'PFT') +
  annotate('segment', x = 0.95, xend = 2.05, y = -20, yend = -20, arrow = arrow(length = unit(0.4, 'cm')), size = 1.5) +
  annotate('text', x = '2000', hjust = -0.23, y = -24, label = 'Drying hyrdoclimate', size = 7) +
  annotate('text', x = '2000', hjust = -1, y = 65, label = 'Liana', size = 7, angle = 0, color = '#E64B35FF', fontface = 2) +
  annotate('text', x = '2000', hjust = -0.45, y = 58, label = expression(bold(paste(Delta, K[req], ' = 21'))), size = 6, angle = 0, color = '#E64B35FF') +
  annotate('text', x = '2000', hjust = -2.55, y = 18, label = 'Tree', size = 7, angle = 0, color = '#4DBBD5FF', fontface = 2) +
  annotate('text', x = '2000', hjust = -1.32, y = 11, label = expression(bold(paste(Delta, K[req], ' = 3'))), size = 6, angle = 0, color = '#4DBBD5FF') +
  coord_cartesian(clip = 'off', ylim = c(-2, 75), xlim = c(1.5, 1.5)) +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26)) +
  theme(plot.margin = unit(c(1, 1, 3, 1), 'lines'))

ggsave('SupplementaryFigure3.jpeg', width = 6, height = 6)

##########################################
#### Part 11: Canopy area sensitivity ####
##########################################

## Supplementary Figures 11-15

##########################
## Main Figure 2 150 m2 ##
##########################

## Supplementary Figure 11

## Before running this section, you must compile
## tree.NPP.mxh & liana.NPP.mxh from NPP_models.R

rm(list = ls())

# Load parameters
load('param.input.RData')

# Load met drivers
load('bci_met_mxh.RData')

# Input K values
nK = 75000
ks = c(trees_K$K / 10, lianas_K$K)
varyk = seq(min(ks), max(ks), length.out = nK)

# Interpolate DBH (and remove highest DBH)
hori.data.liana.dbh = hori.data.liana.dbh[hori.data.liana.dbh < 29]
dbhs = seq(min(hori.data.liana.dbh), max(hori.data.liana.dbh), length.out = 50)

## Tree: 90% canopy, mean DBH, BCI

# Storage
treeouteach = c()
tree_out_r1 = c()

# Run present day with 90% tree canopy at BCI (tropical moist forest)
for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = j
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 1.5 * 100
    frac.tree.al = 0.9
    
    psis = SWP[j]
    
    D = VPD[j,]
    
    if(j == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach)
  if(sumtreeout > 0){
    tree_out_r1[1] = sumtreeout
    tree_out_r1[2] = K
    
    break
  }
}

##Liana: 10% canopy, varied DBH, BCI

# Storage
lianaouteach = c()
liana_out_r1 = matrix(, nrow = length(dbhs), ncol = 2)

# Run present day with 10% liana canopy at BCI (tropical moist forest)
for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 1.5 * 100
      frac.liana.al = 0.1
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach)
    if(sumlianaout > 0){
      liana_out_r1[i,1] = sumlianaout
      liana_out_r1[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

## Tree: 60% canopy, mean DBH, BCI

# Storage
treeouteach = c()
tree_out_r2 = c()

# Run present day with 60% tree canopy at BCI (tropical moist forest)
for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = j
    
    tot.al = 1.5 * 100
    frac.tree.al = 0.6
    
    psis = SWP[j]
    
    D = VPD[j,]
    
    if(j == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach)
  if(sumtreeout > 0){
    tree_out_r2[1] = sumtreeout
    tree_out_r2[2] = K
    
    break
  }
}

## Liana: 40% canopy, varied DBH, BCI

# Storage 
lianaouteach = c()
liana_out_r2 = matrix(, nrow = length(dbhs), ncol = 2)

# Run present day with 40% liana canopy at BCI (tropical moist forest)
for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 1.5 * 100
      frac.liana.al = 0.4
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach)
    if(sumlianaout > 0){
      liana_out_r2[i,1] = sumlianaout
      liana_out_r2[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

# Load new climate drivers (Horizontes)
load('horizontes_met_mxh.RData')

## Tree: 90% canopy, mean DBH, Horizontes

# Storage
treeouteach = c()
tree_out_r3 = c()

# Run present day 90% tree canopy at Horizontes (tropical dry forest)
for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = j
    
    
    tot.al = 1.5 * 100
    frac.tree.al = 0.9
    
    psis = SWP[j]
    
    D = VPD[j,]
    
    if(j == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach)
  if(sumtreeout > 0){
    tree_out_r3[1] = sumtreeout
    tree_out_r3[2] = K
    
    break
  }
}

## Liana: 10% canopy, varied DBH, Horizontes

# Storage
lianaouteach = c()
liana_out_r3 = matrix(, nrow = length(dbhs), ncol = 2)

# Run present day 10% liana canopy at Horizontes (tropical dry forest)
for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 1.5 * 100
      frac.liana.al = 0.1
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach)
    if(sumlianaout > 0){
      liana_out_r3[i,1] = sumlianaout
      liana_out_r3[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

## Tree: 60% canopy, mean DBH, Horizontes

# Storage
treeouteach = c()
tree_out_r4 = c()

# Run present day 60% tree canopy at Horizontes (tropical dry forest)
for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = j
    
    
    tot.al = 1.5 * 100
    frac.tree.al = 0.6
    
    psis = SWP[j]
    
    D = VPD[j,]
    
    if(j == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach)
  if(sumtreeout > 0){
    tree_out_r4[1] = sumtreeout
    tree_out_r4[2] = K
    
    break
  }
}

## Liana: 40% canopy, varied DBH, Horizontes

# Storage
lianaouteach = c()
liana_out_r4 = matrix(, nrow = length(dbhs), ncol = 2)

# Run present day 40% liana canopy at Horizontes (tropical dry forest)
for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 1.5 * 100
      frac.liana.al = 0.4
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach)
    if(sumlianaout > 0){
      liana_out_r4[i,1] = sumlianaout
      liana_out_r4[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

# Formatting
df1 = cbind(dbhs, liana_out_r1[,2], rep('Liana', times = length(dbhs)))
colnames(df1) = c('DBH', 'Kreq', 'PFT')
df1 = rbind(df1, c(mean.hori.data.tree.dbh, tree_out_r1[2], 'Tree'))
df1 = as.data.frame(df1)
df1$DBH = as.numeric(as.character(df1$DBH))
df1$Kreq = as.numeric(as.character(df1$Kreq))
df1$Kreq = df1$Kreq * 0.001

df2 = cbind(dbhs, liana_out_r2[,2], rep('Liana', times = length(dbhs)))
colnames(df2) = c('DBH', 'Kreq', 'PFT')
df2 = rbind(df2, c(mean.hori.data.tree.dbh, tree_out_r2[2], 'Tree'))
df2 = as.data.frame(df2)
df2$DBH = as.numeric(as.character(df2$DBH))
df2$Kreq = as.numeric(as.character(df2$Kreq))
df2$Kreq = df2$Kreq * 0.001

df3 = cbind(dbhs, liana_out_r3[,2], rep('Liana', times = length(dbhs)))
colnames(df3) = c('DBH', 'Kreq', 'PFT')
df3 = rbind(df3, c(mean.hori.data.tree.dbh, tree_out_r3[2], 'Tree'))
df3 = as.data.frame(df3)
df3$DBH = as.numeric(as.character(df3$DBH))
df3$Kreq = as.numeric(as.character(df3$Kreq))
df3$Kreq = df3$Kreq * 0.001

df4 = cbind(dbhs, liana_out_r4[,2], rep('Liana', times = length(dbhs)))
colnames(df4) = c('DBH', 'Kreq', 'PFT')
df4 = rbind(df4, c(mean.hori.data.tree.dbh, tree_out_r4[2], 'Tree'))
df4 = as.data.frame(df4)
df4$DBH = as.numeric(as.character(df4$DBH))
df4$Kreq = as.numeric(as.character(df4$Kreq))
df4$Kreq = df4$Kreq * 0.001

# Add index for simulation scenario
df1 = as.data.frame(cbind(df1, rep(1, nrow(df1))))
colnames(df1)[4] = 'df'
df2 = as.data.frame(cbind(df2, rep(2, nrow(df2))))
colnames(df2)[4] = 'df'
df3 = as.data.frame(cbind(df3, rep(3, nrow(df3))))
colnames(df3)[4] = 'df'
df4 = as.data.frame(cbind(df4, rep(4, nrow(df4))))
colnames(df4)[4] = 'df'
# Combine data frames
df = as.data.frame(rbind(df1, df2, df3, df4))

pl = ggplot(df, aes(x = DBH, y = log(Kreq), group = df)) +
  geom_line(data = subset(df, PFT == 'Liana'), aes(color = as.factor(df)), size = 1.2) +
  geom_hline(data = subset(df, PFT == 'Tree'), aes(yintercept = log(Kreq), color = as.factor(df)), size = 1, linetype = 'dashed') +
  xlab('DBH (cm)') + ylab(expression(paste('log(',K[req],')',' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  #scale_color_viridis_d(name = 'Scenario', end = 0.9, labels = c('1' = '\n10% liana\ntropical moist', '2' = '\n40% liana\ntropical moist\n', '3' = '\n10% liana\ntropical dry\n', '4' = '\n40% liana\ntropical dry\n'), breaks = c('4', '2', '3', '1')) +
  scale_color_manual(name = 'Scenario', 
                     labels = c('1' = '\n10% liana\ntropical moist', '2' = '\n40% liana\ntropical moist\n', '3' = '\n10% liana\ntropical dry\n', '4' = '\n40% liana\ntropical dry\n'), 
                     breaks = c('4', '2', '3', '1'),
                     values = c('4' = '#a6611a', '2' = '#80cdc1', '3' = '#dfc27d', '1' = '#018571')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.text = element_text(size = 26), legend.title = element_text(size = 28, hjust = 0.5), panel.grid = element_blank())

pl = plot_grid(pl, nrow = 1,
               labels = c('A'),
               label_size = 30)

ggsave(pl, filename = 'Plots/SupplementaryFigure11a.jpeg', width = 14, height = 8, units = 'in')

##########################
## Main Figure 2 400 m2 ##
##########################

## Supplementary Figure 11

## Before running this section, you must compile
## tree.NPP.mxh & liana.NPP.mxh from NPP_models.R

rm(list = ls())

# Load parameters
load('param.input.RData')

# Load met drivers
load('bci_met_mxh.RData')

# Input K values
nK = 75000
ks = c(trees_K$K / 10, lianas_K$K)
varyk = seq(min(ks), max(ks), length.out = nK)

# Interpolate DBH (and remove highest DBH)
hori.data.liana.dbh = hori.data.liana.dbh[hori.data.liana.dbh < 29]
dbhs = seq(min(hori.data.liana.dbh), max(hori.data.liana.dbh), length.out = 50)

## Tree: 90% canopy, mean DBH, BCI

# Storage
treeouteach = c()
tree_out_r1 = c()

# Run present day 90% tree canopy at BCI (tropical moist forest)
for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = j
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 4 * 100
    frac.tree.al = 0.9
    
    psis = SWP[j]
    
    D = VPD[j,]
    
    if(j == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach)
  if(sumtreeout > 0){
    tree_out_r1[1] = sumtreeout
    tree_out_r1[2] = K
    
    break
  }
}

##Liana: 10% canopy, varied DBH, BCI

# Storage
lianaouteach = c()
liana_out_r1 = matrix(, nrow = length(dbhs), ncol = 2)

# Run present day 10% liana canopy at BCI (tropical moist forest)
for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach)
    if(sumlianaout > 0){
      liana_out_r1[i,1] = sumlianaout
      liana_out_r1[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

## Tree: 60% canopy, mean DBH, BCI

# Storage
treeouteach = c()
tree_out_r2 = c()

# Run present day 60% tree canopy at BCI (tropical moist forest)
for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = j
    
    tot.al = 4 * 100
    frac.tree.al = 0.6
    
    psis = SWP[j]
    
    D = VPD[j,]
    
    if(j == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach)
  if(sumtreeout > 0){
    tree_out_r2[1] = sumtreeout
    tree_out_r2[2] = K
    
    break
  }
}

## Liana: 40% canopy, varied DBH, BCI

# Storage
lianaouteach = c()
liana_out_r2 = matrix(, nrow = length(dbhs), ncol = 2)

# Run present day 40% canopy at BCI (tropical moist forest)
for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach)
    if(sumlianaout > 0){
      liana_out_r2[i,1] = sumlianaout
      liana_out_r2[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

# Load new climate drivers (Horizontes)
load('horizontes_met_mxh.RData')

## Tree: 90% canopy, mean DBH, Horizontes

# Storage
treeouteach = c()
tree_out_r3 = c()

# Run present day 90% tree canopy at Horizontes (tropical dry forest)
for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = j
    
    
    tot.al = 4 * 100
    frac.tree.al = 0.9
    
    psis = SWP[j]
    
    D = VPD[j,]
    
    if(j == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach)
  if(sumtreeout > 0){
    tree_out_r3[1] = sumtreeout
    tree_out_r3[2] = K
    
    break
  }
}

## Liana: 10% canopy, varied DBH, Horizontes

# Storage
lianaouteach = c()
liana_out_r3 = matrix(, nrow = length(dbhs), ncol = 2)

# Run present day 10% liana canopy at Horizontes (tropical dry forest)
for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach)
    if(sumlianaout > 0){
      liana_out_r3[i,1] = sumlianaout
      liana_out_r3[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

## Tree: 60% canopy, mean DBH, Horizontes

# Storage
treeouteach = c()
tree_out_r4 = c()

# Run present day 60% tree canopy at Horizontes (tropical dry forest)
for(k in 1:length(varyk)){
  for(j in 1:12){
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = j
    
    
    tot.al = 4 * 100
    frac.tree.al = 0.6
    
    psis = SWP[j]
    
    D = VPD[j,]
    
    if(j == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    
    K = varyk[k]
    tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height)
    treeouteach[j] = tree.out[1]
  }
  sumtreeout = sum(treeouteach)
  if(sumtreeout > 0){
    tree_out_r4[1] = sumtreeout
    tree_out_r4[2] = K
    
    break
  }
}

## Liana: 40% canopy, varied DBH, Horizontes

# Storage
lianaouteach = c()
liana_out_r4 = matrix(, nrow = length(dbhs), ncol = 2)

# Run present day 40% liana canopy at Horizontes (tropical dry forest)
for(i in 1:length(dbhs)){
  for(k in 1:length(varyk)){
    for(j in 1:12){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = j
      
      dbh = dbhs[i]
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = SWP[j]
      
      D = VPD[j,]
      
      if(j == 1){
        liana.length = (tree.b1Ht * tree.dbh^tree.b2Ht)
      }else{
        liana.length = liana.out[4]
      }
      K = varyk[k]
      liana.out = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length)
      lianaouteach[j] = liana.out[1]
    }
    sumlianaout = sum(lianaouteach)
    if(sumlianaout > 0){
      liana_out_r4[i,1] = sumlianaout
      liana_out_r4[i,2] = K
      
      break
    }
  }
  print(paste('Done with DBH',i,'of',length(dbhs)))
}

# Formatting
df1 = cbind(dbhs, liana_out_r1[,2], rep('Liana', times = length(dbhs)))
colnames(df1) = c('DBH', 'Kreq', 'PFT')
df1 = rbind(df1, c(mean.hori.data.tree.dbh, tree_out_r1[2], 'Tree'))
df1 = as.data.frame(df1)
df1$DBH = as.numeric(as.character(df1$DBH))
df1$Kreq = as.numeric(as.character(df1$Kreq))
df1$Kreq = df1$Kreq * 0.001

df2 = cbind(dbhs, liana_out_r2[,2], rep('Liana', times = length(dbhs)))
colnames(df2) = c('DBH', 'Kreq', 'PFT')
df2 = rbind(df2, c(mean.hori.data.tree.dbh, tree_out_r2[2], 'Tree'))
df2 = as.data.frame(df2)
df2$DBH = as.numeric(as.character(df2$DBH))
df2$Kreq = as.numeric(as.character(df2$Kreq))
df2$Kreq = df2$Kreq * 0.001

df3 = cbind(dbhs, liana_out_r3[,2], rep('Liana', times = length(dbhs)))
colnames(df3) = c('DBH', 'Kreq', 'PFT')
df3 = rbind(df3, c(mean.hori.data.tree.dbh, tree_out_r3[2], 'Tree'))
df3 = as.data.frame(df3)
df3$DBH = as.numeric(as.character(df3$DBH))
df3$Kreq = as.numeric(as.character(df3$Kreq))
df3$Kreq = df3$Kreq * 0.001

df4 = cbind(dbhs, liana_out_r4[,2], rep('Liana', times = length(dbhs)))
colnames(df4) = c('DBH', 'Kreq', 'PFT')
df4 = rbind(df4, c(mean.hori.data.tree.dbh, tree_out_r4[2], 'Tree'))
df4 = as.data.frame(df4)
df4$DBH = as.numeric(as.character(df4$DBH))
df4$Kreq = as.numeric(as.character(df4$Kreq))
df4$Kreq = df4$Kreq * 0.001

# Add index for each scenario
df1 = as.data.frame(cbind(df1, rep(1, nrow(df1))))
colnames(df1)[4] = 'df'
df2 = as.data.frame(cbind(df2, rep(2, nrow(df2))))
colnames(df2)[4] = 'df'
df3 = as.data.frame(cbind(df3, rep(3, nrow(df3))))
colnames(df3)[4] = 'df'
df4 = as.data.frame(cbind(df4, rep(4, nrow(df4))))
colnames(df4)[4] = 'df'
# Combine dataframes
df = as.data.frame(rbind(df1, df2, df3, df4))

# Plot
pl = ggplot(df, aes(x = DBH, y = log(Kreq), group = df)) +
  geom_line(data = subset(df, PFT == 'Liana'), aes(color = as.factor(df)), size = 1.2) +
  geom_hline(data = subset(df, PFT == 'Tree'), aes(yintercept = log(Kreq), color = as.factor(df)), size = 1, linetype = 'dashed') +
  xlab('DBH (cm)') + ylab(expression(paste('log(',K[req],')',' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  #scale_color_viridis_d(name = 'Scenario', end = 0.9, labels = c('1' = '\n10% liana\ntropical moist', '2' = '\n40% liana\ntropical moist\n', '3' = '\n10% liana\ntropical dry\n', '4' = '\n40% liana\ntropical dry\n'), breaks = c('4', '2', '3', '1')) +
  scale_color_manual(name = 'Scenario', 
                     labels = c('1' = '\n10% liana\ntropical moist', '2' = '\n40% liana\ntropical moist\n', '3' = '\n10% liana\ntropical dry\n', '4' = '\n40% liana\ntropical dry\n'), 
                     breaks = c('4', '2', '3', '1'),
                     values = c('4' = '#a6611a', '2' = '#80cdc1', '3' = '#dfc27d', '1' = '#018571')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.text = element_text(size = 26), legend.title = element_text(size = 28, hjust = 0.5), panel.grid = element_blank())

pl = plot_grid(pl,
               nrow = 1,
               labels = c('B'),
               label_size = 30)

ggsave(pl, filename = 'Plots/SupplementaryFigure11b.jpeg', width = 14, height = 8, units = 'in')

#########################
## Interpolation plots ##
#########################

## Supplementary Figure 12

## Before running this section, you must compile
## tree.NPP.mxh & liana.NPP.mxh from NPP_models.R

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Interpolation_mxh_100.RData')

# Indexing
nsite = 100
nmonth = 12
nK = 75000

# Number of days in each month
nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Renaming climate inputs
VPDs = VPD_interp_100
SWPs = SWP_interp_100

#########################
## Established: 150 m2 ##
#########################

# Input K values
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

# Storage
tree.req = matrix(, nrow = nsite, ncol = nsite)

# Run present day established scenario with 150 m2 canopy
for(vsite in 1:nsite){
  for(ssite in 1:nsite){
    sens = c()
    for(k in 1:nK){
      for(mo in 1:nmonth){
        nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        month = mo
        
        tree.dbh = mean.hori.data.tree.dbh
        
        tot.al = 1.5 * 100
        frac.tree.al = 0.6
        
        psis = SWPs[mo,ssite]
        
        D = VPDs[mo,,vsite]
        
        if(mo == 1){
          tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
        }else{
          tree.height = tree.out[5]
        }
        K = Ks[k]
        tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
        sens[mo] = tree.out[1]
      }
      sum = sum(sens)
      if(sum > 0){
        tree.req[vsite,ssite] = K
        break
      }
    }
    print(paste('Passing sites',vsite,'&',ssite))
  }
}

# Save output due to high computation time
#save(tree.req, file = 'tree.req.100interp.monthly_150.RData')

# Load once this has been run once
load(file = 'tree.req.100interp.monthly_150.RData')
# Process
tree.req = tree.req * 0.001
tree_req_melt = melt(tree.req)
colnames(tree_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

# Contour plot
ra1 = ggplot(tree_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv, fill = Ksurv)) +
  geom_raster(interpolate = T, show.legend = F) +
  geom_contour(alpha = 0, aes(color = ..level..), breaks = seq(0, 10, by = 1), show.legend = T) +
  geom_contour(alpha = 0.65, color = 'black', breaks = seq(0, 10, by = 1), show.legend = F) +
  xlab('VPD index') + ylab(expression(paste(Psi,' index'))) +
  scale_color_distiller(name = expression(K[req]), palette = 'PuBu', breaks = seq(0, 10, by = 1)) +
  scale_fill_distiller(name = expression(K[req]), palette = 'PuBu', breaks = seq(0, 10, by = 1)) +
  ggtitle('Tree') +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank()) +
  geom_dl(aes(label = ..level..), stat = 'contour', breaks = seq(0, 10, by = 1), method = list('top.pieces', color = 'black', cex = 1.8, vjust = -0.4))
ra1

# Make coarser input K values to improve computation time
ks = c(trees_K$K / 10, lianas_K$K)
nK = 5000
Ks = seq(min(ks), max(ks), length.out = nK)

# Storage
liana.req = matrix(, nrow = nsite, ncol = nsite)

# Run present day established scenario with 150 m2 canopy
for(vsite in 1:nsite){
  for(ssite in 1:nsite){
    sens = c()
    for(k in 1:nK){
      for(mo in 1:12){
        nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        month = mo
        
        dbh = mean.data.liana.dbh
        
        tot.al = 1.5 * 100
        frac.liana.al = 0.4
        
        psis = SWPs[mo,ssite]
        
        D = VPDs[mo,,vsite]
        
        if(mo == 1){
          liana.length = tree.out[2]
        }else{
          liana.length = lianaout[4]
        }
        K = Ks[k]
        lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
        sens[mo] = lianaout[1]
      }
      sum = sum(sens)
      if(sum > 0){
        liana.req[vsite, ssite] = K
        break
      }
    }
    print(paste('Passing sites',vsite,'&',ssite))
  }
}

# Save and load as above
#save(liana.req, file = 'liana.req.100interp.monthly_150.RData')

load(file = 'liana.req.100interp.monthly_150.RData')
# Processing
liana.req = liana.req * 0.001
liana_req_melt = melt(liana.req)
colnames(liana_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

# Contour plot
ra2 = ggplot(liana_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv, fill = Ksurv)) +
  geom_raster(interpolate = T, show.legend = T) +
  geom_contour(alpha = 0, aes(color = ..level..), show.legend = F, breaks = seq(30, 210, by = 20)) +
  geom_contour(alpha = 0.65, color = 'black', breaks = seq(30, 210, by = 20), show.legend = F) +
  xlab('VPD index') + ylab(expression(paste(Psi,' index'))) +
  scale_fill_distiller(name = expression(K[req]), palette = 'BuGn', breaks = seq(30, 210, by = 30)) +
  ggtitle('Liana') +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank()) +
  geom_dl(aes(label = ..level..), stat = 'contour', breaks = seq(30, 210, by = 20), method = list('top.pieces', color = 'black', cex = 1.8, vjust = -0.4))
ra2

# Format plots
ra1_fin = ra1 + theme(legend.position = 'bottom', plot.margin = unit(c(0, 0, 0, 0), 'cm'))
ra2_fin = ra2 + theme(legend.position = 'bottom', plot.margin = unit(c(-1, 0, 1, 0), 'cm'))

# Plot together
ga = plot_grid(ra1_fin, ra2_fin, 
               nrow = 1,
               align = 'hv', axis = 'tbrl',
               labels = c('A', 'B'),
               label_size = 30)

ggsave(ga, filename = 'Plots/SupplementaryFigure12a.jpeg', height = 8, width = 14, units = 'in')

#########################
## Established: 400 m2 ##
#########################

## Before running this section, you must compile
## tree.NPP.mxh & liana.NPP.mxh from NPP_models.R

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Interpolation_mxh_100.RData')

# Indexing
nsite = 100
nmonth = 12
nK = 75000

# Number of days in each month
nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Rename climate inputs
VPDs = VPD_interp_100
SWPs = SWP_interp_100

# Input K values
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

# Storage
tree.req = matrix(, nrow = nsite, ncol = nsite)

# Run present day established scenario with 400 m2 canopy
for(vsite in 1:nsite){
  for(ssite in 1:nsite){
    sens = c()
    for(k in 1:nK){
      for(mo in 1:nmonth){
        nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        month = mo
        
        tree.dbh = mean.hori.data.tree.dbh
        
        tot.al = 4 * 100
        frac.tree.al = 0.6
        
        psis = SWPs[mo,ssite]
        
        D = VPDs[mo,,vsite]
        
        if(mo == 1){
          tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
        }else{
          tree.height = tree.out[5]
        }
        K = Ks[k]
        tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
        sens[mo] = tree.out[1]
      }
      sum = sum(sens)
      if(sum > 0){
        tree.req[vsite,ssite] = K
        break
      }
    }
    print(paste('Passing sites',vsite,'&',ssite))
  }
}

# Save and load as above
#save(tree.req, file = 'tree.req.100interp.monthly.RData')

load(file = 'tree.req.100interp.monthly.RData')
# Process
tree.req = tree.req * 0.001
tree_req_melt = melt(tree.req)
colnames(tree_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

ra1 = ggplot(tree_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv, fill = Ksurv)) +
  geom_raster(interpolate = T, show.legend = F) +
  geom_contour(alpha = 0, aes(color = ..level..), breaks = seq(0, 10, by = 1), show.legend = T) +
  geom_contour(alpha = 0.65, color = 'black', breaks = seq(0, 10, by = 1), show.legend = F) +
  xlab('VPD index') + ylab(expression(paste(Psi,' index'))) +
  scale_color_distiller(name = expression(K[req]), palette = 'PuBu', breaks = seq(0, 10, by = 2)) +
  scale_fill_distiller(name = expression(K[req]), palette = 'PuBu', breaks = seq(0, 10, by = 2)) +
  ggtitle('Tree') +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank()) +
  geom_dl(aes(label = ..level..), stat = 'contour', breaks = seq(0, 10, by = 1), method = list('top.pieces', color = 'black', cex = 1.8, vjust = -0.4))
ra1

# Make coarser input K values and truncate
# Truncation was based on running the model with the wettest hydroclimate
# scenario and using a K value lower than Kreq for that scenario
ks = c(trees_K$K / 10, lianas_K$K)
nK = 5000
Ks = seq(min(ks), max(ks), length.out = nK)
Ks = Ks[Ks > 75925]

# Storage
liana.req = matrix(, nrow = nsite, ncol = nsite)

# Run present day established scenario with 400 m2 canopy
for(vsite in 1:nsite){
  for(ssite in 1:nsite){
    sens = c()
    for(k in 1:nK){
      for(mo in 1:12){
        nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        month = mo
        
        dbh = mean.data.liana.dbh
        
        tot.al = 4 * 100
        frac.liana.al = 0.4
        
        psis = SWPs[mo,ssite]
        
        D = VPDs[mo,,vsite]
        
        if(mo == 1){
          liana.length = tree.out[2]
        }else{
          liana.length = lianaout[4]
        }
        K = Ks[k]
        lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
        sens[mo] = lianaout[1]
      }
      sum = sum(sens)
      if(sum > 0){
        liana.req[vsite, ssite] = K
        break
      }
    }
    print(paste('Passing sites',vsite,'&',ssite))
  }
}

# Save and load as before
#save(liana.req, file = 'liana.req.100interp.monthly.RData')

load(file = 'liana.req.100interp.monthly.RData')
# Processing
liana.req = liana.req * 0.001
liana_req_melt = melt(liana.req)
colnames(liana_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

# Contour plot
ra2 = ggplot(liana_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv, fill = Ksurv)) +
  geom_raster(interpolate = T, show.legend = T) +
  geom_contour(alpha = 0, aes(color = ..level..), show.legend = F, breaks = seq(30, 210, by = 20)) +
  geom_contour(alpha = 0.65, color = 'black', breaks = seq(30, 210, by = 20), show.legend = F) +
  xlab('VPD index') + ylab(expression(paste(Psi,' index'))) +
  scale_fill_distiller(name = expression(K[req]), palette = 'BuGn', breaks = seq(30, 210, by = 60)) +
  ggtitle('Liana') +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank()) +
  geom_dl(aes(label = ..level..), stat = 'contour', breaks = seq(30, 210, by = 20), method = list('top.pieces', color = 'black', cex = 1.8, vjust = -0.4))
ra2

# Format plots
ra1_fin = ra1 + theme(legend.position = 'bottom', plot.margin = unit(c(0, 0, 0, 0), 'cm'))
ra2_fin = ra2 + theme(legend.position = 'bottom', plot.margin = unit(c(-1, 0, 1, 0), 'cm'))

# Plot together
ga = plot_grid(ra1_fin, ra2_fin, 
               nrow = 1,
               align = 'hv', axis = 'tbrl',
               labels = c('C', 'D'),
               label_size = 30)

ggsave(ga, filename = 'Plots/SupplementaryFigure12b.jpeg', height = 8, width = 14, units = 'in')

###############################
## Invasion scenario: 150 m2 ##
###############################

## Supplementary Figure 13

## Before running this section, you must compile
## tree.NPP.mxh & liana.NPP.mxh from NPP_models.R

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Interpolation_mxh_100.RData')

# Indexing
nsite = 100
nmonth = 12
nK = 75000

# Number of days in each month
nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Rename climate inputs
VPDs = VPD_interp_100
SWPs = SWP_interp_100

# Input K values
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

# Storage
tree.req = matrix(, nrow = nsite, ncol = nsite)

# Run present day invading liana scenario with 150 m2 canopy
for(vsite in 1:nsite){
  for(ssite in 1:nsite){
    sens = c()
    for(k in 1:nK){
      for(mo in 1:nmonth){
        nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        month = mo
        
        tree.dbh = mean.hori.data.tree.dbh
        
        tot.al = 1.5 * 100
        frac.tree.al = 0.9
        
        psis = SWPs[mo,ssite]
        
        D = VPDs[mo,,vsite]
        
        if(mo == 1){
          tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
        }else{
          tree.height = tree.out[5]
        }
        K = Ks[k]
        tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
        sens[mo] = tree.out[1]
      }
      sum = sum(sens)
      if(sum > 0){
        tree.req[vsite,ssite] = K
        break
      }
    }
    print(paste('Passing sites',vsite,'&',ssite))
  }
}

# Save and load as above
#save(tree.req, file = 'tree.req.100interp.monthly_invasion_150.RData')

load(file = 'tree.req.100interp.monthly_invasion_150.RData')
# Processing
tree.req = tree.req * 0.001
tree_req_melt = melt(tree.req)
colnames(tree_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

# Contour plot
ra1 = ggplot(tree_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv, fill = Ksurv)) +
  geom_raster(interpolate = T, show.legend = F) +
  geom_contour(alpha = 0, aes(color = ..level..), breaks = seq(2, 12, by = 1), show.legend = T) +
  geom_contour(alpha = 0.65, color = 'black', breaks = seq(2, 12, by = 1), show.legend = F) +
  xlab('VPD index') + ylab(expression(paste(Psi,' index'))) +
  scale_color_distiller(name = expression(K[req]), palette = 'PuBu', breaks = seq(2, 12, by = 1)) +
  scale_fill_distiller(name = expression(K[req]), palette = 'PuBu', breaks = seq(2, 12, by = 1)) +
  ggtitle('Tree') +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank()) +
  geom_dl(aes(label = ..level..), stat = 'contour', breaks = seq(2, 12, by = 1), method = list('top.pieces', color = 'black', cex = 1.8, vjust = -0.4))
ra1

# Input K values truncated as above
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)
Ks = Ks[Ks > 13525]

# Storage
liana.req = matrix(, nrow = nsite, ncol = nsite)

# Run present day invasion scenario with 150 m2 canopy
for(vsite in 1:nsite){
  for(ssite in 1:nsite){
    sens = c()
    for(k in 1:nK){
      for(mo in 1:12){
        nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        month = mo
        
        dbh = 2
        
        tot.al = 1.5 * 100
        frac.liana.al = 0.1
        
        psis = SWPs[mo,ssite]
        
        D = VPDs[mo,,vsite]
        
        if(mo == 1){
          liana.length = tree.out[2]
        }else{
          liana.length = lianaout[4]
        }
        K = Ks[k]
        lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
        sens[mo] = lianaout[1]
      }
      sum = sum(sens)
      if(sum > 0){
        liana.req[vsite, ssite] = K
        break
      }
    }
    print(paste('Passing sites',vsite,'&',ssite))
  }
}

# Save and load as above
#save(liana.req, file = 'liana.req.100interp.monthly_invasion_150.RData')

load(file = 'liana.req.100interp.monthly_invasion_150.RData')
# Processing
liana.req = liana.req * 0.001
liana_req_melt = melt(liana.req)
colnames(liana_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

# Contour plot
ra2 = ggplot(liana_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv, fill = Ksurv)) +
  geom_raster(interpolate = T, show.legend = T) +
  geom_contour(alpha = 0, aes(color = ..level..), show.legend = F, breaks = seq(10, 90, by = 5)) +
  geom_contour(alpha = 0.65, color = 'black', breaks = seq(10, 90, by = 5), show.legend = F) +
  xlab('VPD index') + ylab(expression(paste(Psi,' index'))) +
  scale_fill_distiller(name = expression(K[req]), palette = 'BuGn', breaks = seq(10, 90, by = 10)) +
  ggtitle('Liana') +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank()) +
  geom_dl(aes(label = ..level..), stat = 'contour', breaks = seq(10, 90, by = 5), method = list('top.pieces', color = 'black', cex = 1.8, vjust = -0.4))
ra2

# Plot together
ga_fin = plot_grid(ra1, ra2, 
                   nrow = 1, 
                   rel_widths = 1,
                   labels = c('A', 'B'),
                   label_size = 30)

ggsave(ga_fin, filename = 'Plots/SupplementaryFigure13a.jpeg', height = 8, width = 14, units = 'in')

###############################
## Invasion scenario: 400 m2 ##
###############################

## Before running this section, you must compile
## tree.NPP.mxh & liana.NPP.mxh from NPP_models.R

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Interpolation_mxh_100.RData')

# Indexing
nsite = 100
nmonth = 12
nK = 75000

# Number of days in each month
nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Rename climate inputs
VPDs = VPD_interp_100
SWPs = SWP_interp_100

# Input K values
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

# Storage
tree.req = matrix(, nrow = nsite, ncol = nsite)

# Run present day invasion scenario with 400 m2 canopy
for(vsite in 1:nsite){
  for(ssite in 1:nsite){
    sens = c()
    for(k in 1:nK){
      for(mo in 1:nmonth){
        nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        month = mo
        
        tree.dbh = mean.hori.data.tree.dbh
        
        tot.al = 4 * 100
        frac.tree.al = 0.9
        
        psis = SWPs[mo,ssite]
        
        D = VPDs[mo,,vsite]
        
        if(mo == 1){
          tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
        }else{
          tree.height = tree.out[5]
        }
        K = Ks[k]
        tree.out = tree.NPP.mxh(tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
        sens[mo] = tree.out[1]
      }
      sum = sum(sens)
      if(sum > 0){
        tree.req[vsite,ssite] = K
        break
      }
    }
    print(paste('Passing sites',vsite,'&',ssite))
  }
}

# Save and load as above
#save(tree.req, file = 'tree.req.100interp.monthly_invasion.RData')

load(file = 'tree.req.100interp.monthly_invasion.RData')
# Processing
tree.req = tree.req * 0.001
tree_req_melt = melt(tree.req)
colnames(tree_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

# Contour plot
ra1 = ggplot(tree_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv, fill = Ksurv)) +
  geom_raster(interpolate = T, show.legend = F) +
  geom_contour(alpha = 0, aes(color = ..level..), breaks = seq(2, 12, by = 2), show.legend = T) +
  geom_contour(alpha = 0.65, color = 'black', breaks = seq(2, 12, by = 2), show.legend = F) +
  xlab('VPD index') + ylab(expression(paste(Psi,' index'))) +
  scale_color_distiller(name = expression(K[req]), palette = 'PuBu', breaks = seq(4, 12, by = 2)) +
  scale_fill_distiller(name = expression(K[req]), palette = 'PuBu', breaks = seq(4, 12, by = 2)) +
  ggtitle('Tree') +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank()) +
  geom_dl(aes(label = ..level..), stat = 'contour', breaks = seq(2, 12, by = 2), method = list('top.pieces', color = 'black', cex = 1.8, vjust = -0.4))
ra1

# Input K values truncated as above
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)
Ks = Ks[Ks > 33727]

# Storage
liana.req = matrix(, nrow = nsite, ncol = nsite)

# Run present day invasion scenario with 400 m2 canopy
for(vsite in 1:nsite){
  for(ssite in 1:nsite){
    sens = c()
    for(k in 1:nK){
      for(mo in 1:12){
        nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        month = mo
        
        dbh = 2
        
        tot.al = 4 * 100
        frac.liana.al = 0.1
        
        psis = SWPs[mo,ssite]
        
        D = VPDs[mo,,vsite]
        
        if(mo == 1){
          liana.length = tree.out[2]
        }else{
          liana.length = lianaout[4]
        }
        K = Ks[k]
        lianaout = liana.NPP.mxh(dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
        sens[mo] = lianaout[1]
      }
      sum = sum(sens)
      if(sum > 0){
        liana.req[vsite, ssite] = K
        break
      }
    }
    print(paste('Passing sites',vsite,'&',ssite))
  }
}

# Save and load as above
#save(liana.req, file = 'liana.req.100interp.monthly_invasion.RData')

load(file = 'liana.req.100interp.monthly_invasion.RData')
# Processing
liana.req = liana.req * 0.001
liana_req_melt = melt(liana.req)
colnames(liana_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

# Contour plot
ra2 = ggplot(liana_req_melt, aes(x = VPD_site, y = SWP_site, z = Ksurv, fill = Ksurv)) +
  geom_raster(interpolate = T, show.legend = T) +
  geom_contour(alpha = 0, aes(color = ..level..), show.legend = F, breaks = seq(10, 90, by = 10)) +
  geom_contour(alpha = 0.65, color = 'black', breaks = seq(10, 90, by = 10), show.legend = F) +
  xlab('VPD index') + ylab(expression(paste(Psi,' index'))) +
  scale_fill_distiller(name = expression(K[req]), palette = 'BuGn', breaks = seq(10, 90, by = 20)) +
  ggtitle('Liana') +
  coord_fixed(ratio = 1, ylim = c(0, 103)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), plot.title = element_text(size = 30, hjust = 0.5), legend.title = element_text(size = 28), legend.text = element_text(size = 26), axis.text = element_text(size = 26)) +
  theme(panel.border = element_blank(), legend.position = 'bottom', panel.grid.minor = element_blank()) +
  geom_dl(aes(label = ..level..), stat = 'contour', breaks = seq(10, 90, by = 10), method = list('top.pieces', color = 'black', cex = 1.8, vjust = -0.4))
ra2

# Plot together
ga_fin = plot_grid(ra1, ra2, 
                   nrow = 1, 
                   rel_widths = 1,
                   labels = c('C', 'D'),
                   label_size = 30)

ggsave(ga_fin, filename = 'Plots/SupplementaryFigure13b.jpeg', height = 8, width = 14, units = 'in')

#####################################
## Future Kreq: 150 m2 established ##
#####################################

## Supplementary Figure 14

## Before running this section, you must compile
## tree.NPP.mxh.ca & liana.NPP.mxh.ca from NPP_models_wCO2.R

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Driver_calc/Horizontes/horizontes_met_mxh.RData')
H_VPD = VPD
H_SWP = SWP
rm(VPD, SWP)

# Indexing months and days
nmonth = 12
nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

nK = 75000

# Horizontes future scenarios
H_VPD_future = array(, dim = c(12, 24, 2))
H_VPD_future[,,1] = H_VPD
H_VPD_future[,,2] = H_VPD * 2

#### Tree simulation - present day

# Input K values
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

# Storage
tree.req_h = c()
sens = c()

# Run present day established liana scenario with 150 m2 canopy
for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 400
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 1.5 * 100
    frac.tree.al = 0.6
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,1]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h[1] = K
    break
  }
}

#### Liana simulation -- present day

# Storage
liana.req_h = c()
sens = c()

# Run present day established liana scenario with 150 m2 canopy
for(k in 1:nK){
  for(mo in 1:12){
    Ca = 400
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = mean.data.liana.dbh
    
    tot.al = 1.5 * 100
    frac.liana.al = 0.4
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,1]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h[1] = K
    break
  }
}

#### Tree simulation -- future scenario

# Storage
sens = c()

# Run future established liana scenario with 150 m2 canopy
for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 1.5 * 100
    frac.tree.al = 0.6
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,2]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h[2] = K
    break
  }
}

#### Liana simulation -- future scenario

# Storage
sens = c()

# Run future established liana scenario with 150 m2 canopy
for(k in 1:nK){
  for(mo in 1:12){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = mean.data.liana.dbh
    
    tot.al = 1.5 * 100
    frac.liana.al = 0.4
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,2]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h[2] = K
    break
  }
}

tree.req_h = tree.req_h * 0.001
liana.req_h = liana.req_h * 0.001

# Make workable dataframes
liana_kreq = cbind(c('2000', '2100'), 
                   c(liana.req_h[1], liana.req_h[2]))
colnames(liana_kreq) = c('Time', 'Kreq')

tree_kreq = cbind(c('2000', '2100'),
                  c(tree.req_h[1], tree.req_h[2]))
colnames(tree_kreq) = c('Time', 'Kreq')

# Combine and format dataframes
comb_kreq = rbind(liana_kreq, tree_kreq)
comb_kreq = as.data.frame(comb_kreq)
comb_kreq$growth.form = c('Liana', 'Liana', 
                          'Tree', 'Tree')
comb_kreq$Time = as.factor(comb_kreq$Time)
comb_kreq$Kreq = as.numeric(comb_kreq$Kreq)
comb_kreq$growth.form = as.factor(comb_kreq$growth.form)

# Plot
p1 = ggplot(comb_kreq, aes(x = Time, y = Kreq, color = growth.form, group = growth.form)) +
  geom_point(size = 4, show.legend = F) +
  geom_line(size = 1.5, show.legend = F, alpha = 0.7) +
  scale_x_discrete(limits = c('2000', '2100')) +
  xlab('') + ylab(expression(paste(K[req], ' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  scale_color_npg(name = 'PFT') +
  ggtitle(expression(paste('Total leaf area = 150', m^2))) +
  annotate('segment', x = 0.95, xend = 2.05, y = -25, yend = -25, arrow = arrow(length = unit(0.4, 'cm')), size = 1.5) +
  annotate('text', x = '2000', hjust = -0.23, y = -32, label = 'Drying hyrdoclimate', size = 7) +
  annotate('text', x = '2000', hjust = -1, y = 110, label = 'Liana', size = 7, angle = 0, color = '#E64B35FF', fontface = 2) +
  annotate('text', x = '2000', hjust = -0.45, y = 98, label = expression(bold(paste(Delta, K[req], ' = 35'))), size = 6, angle = 0, color = '#E64B35FF') +
  annotate('text', x = '2000', hjust = -2.55, y = 22, label = 'Tree', size = 7, angle = 0, color = '#4DBBD5FF', fontface = 2) +
  annotate('text', x = '2000', hjust = -1.32, y = 10, label = expression(bold(paste(Delta, K[req], ' = 1'))), size = 6, angle = 0, color = '#4DBBD5FF') +
  coord_cartesian(clip = 'off', ylim = c(-2, 120), xlim = c(1.5, 1.5)) +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26)) +
  theme(plot.title = element_text(size = 28, hjust = 0.5)) +
  theme(plot.margin = unit(c(1, 1, 3, 1), 'lines'))
p1

############################################
## Kreq-Kw comparison: 400 m2 established ##
############################################

## This section should be run in combination with the previous section

# Storage
tree.req_h_2 = c()
sens = c()

# Run present day established scenario with 400 m2 canopy
for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 400
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 4 * 100
    frac.tree.al = 0.6
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,1]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h_2[1] = K
    break
  }
}

#### Liana simulation -- present day

# Storage
liana.req_h_2 = c()
sens = c()

# Run present day established liana scenario with 400 m2 canopy
for(k in 1:nK){
  for(mo in 1:12){
    Ca = 400
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = mean.data.liana.dbh
    
    tot.al = 4 * 100
    frac.liana.al = 0.4
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,1]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h_2[1] = K
    break
  }
}

#### Tree simulation -- future scenario

# Storage
sens = c()

# Run future established liana scenario with 400 m2 canopy
for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 4 * 100
    frac.tree.al = 0.6
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,2]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h_2[2] = K
    break
  }
}

#### Liana simulation -- future scenario

# Storage
sens = c()

# Run future established liana scenario with 400 m2 canopy area
for(k in 1:nK){
  for(mo in 1:12){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = mean.data.liana.dbh
    
    tot.al = 4 * 100
    frac.liana.al = 0.4
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,2]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h_2[2] = K
    break
  }
}

tree.req_h_2 = tree.req_h_2 * 0.001
liana.req_h_2 = liana.req_h_2 * 0.001

# Make workable dataframes
liana_kreq = cbind(c('2000', '2100'), 
                   c(liana.req_h_2[1], liana.req_h_2[2]))
colnames(liana_kreq) = c('Time', 'Kreq')

tree_kreq = cbind(c('2000', '2100'),
                  c(tree.req_h_2[1], tree.req_h_2[2]))
colnames(tree_kreq) = c('Time', 'Kreq')

# Combine and format dataframes
comb_kreq_2 = rbind(liana_kreq, tree_kreq)
comb_kreq_2 = as.data.frame(comb_kreq_2)
comb_kreq_2$growth.form = c('Liana', 'Liana', 
                            'Tree', 'Tree')
comb_kreq_2$Time = as.factor(comb_kreq_2$Time)
comb_kreq_2$Kreq = as.numeric(comb_kreq_2$Kreq)
comb_kreq_2$growth.form = as.factor(comb_kreq_2$growth.form)

# Plot
p2 = ggplot(comb_kreq_2, aes(x = Time, y = Kreq, color = growth.form, group = growth.form)) +
  geom_point(size = 4, show.legend = F) +
  geom_line(size = 1.5, show.legend = F, alpha = 0.7) +
  scale_x_discrete(limits = c('2000', '2100')) +
  xlab('') + ylab(expression(paste(K[req], ' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  ggtitle(expression(paste('Total leaf area = 400', m^2))) +
  scale_color_npg(name = 'PFT') +
  annotate('segment', x = 0.95, xend = 2.05, y = -70, yend = -70, arrow = arrow(length = unit(0.4, 'cm')), size = 1.5) +
  annotate('text', x = '2000', hjust = -0.23, y = -90, label = 'Drying hyrdoclimate', size = 7) +
  annotate('text', x = '2000', hjust = -1, y = 285, label = 'Liana', size = 7, angle = 0, color = '#E64B35FF', fontface = 2) +
  annotate('text', x = '2000', hjust = -0.45, y = 255, label = expression(bold(paste(Delta, K[req], ' = 92'))), size = 6, angle = 0, color = '#E64B35FF') +
  annotate('text', x = '2000', hjust = -2.55, y = 55, label = 'Tree', size = 7, angle = 0, color = '#4DBBD5FF', fontface = 2) +
  annotate('text', x = '2000', hjust = -1.32, y = 25, label = expression(bold(paste(Delta, K[req], ' = 4'))), size = 6, angle = 0, color = '#4DBBD5FF') +
  coord_cartesian(clip = 'off', ylim = c(-2, 320), xlim = c(1.5, 1.5)) +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26)) +
  theme(plot.title = element_text(size = 28, hjust = 0.5)) +
  theme(plot.margin = unit(c(1, 1, 3, 1), 'lines'))
p2

# Plot together
ga = plot_grid(p1, p2,
               labels = c('A', 'B'),
               label_size = 30)

ggsave(plot = ga, 'Plots/SupplementaryFigure14.jpeg', width = 12, height = 6)

#########################################
## Kreq-Kw comparison: 150 m2 invasion ##
#########################################

## Before running this section, you must compile
## tree.NPP.mxh.ca & liana.NPP.mxh.ca from NPP_models_wCO2.R

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Driver_calc/Horizontes/horizontes_met_mxh.RData')
H_VPD = VPD
H_SWP = SWP
rm(VPD, SWP)

# Indexing months and days
nmonth = 12
nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

nK = 75000

# Horizontes future scenarios
H_VPD_future = array(, dim = c(12, 24, 2))
H_VPD_future[,,1] = H_VPD
H_VPD_future[,,2] = H_VPD * 2

# Input K values
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

# Storage
tree.req_h = c()
sens = c()

# Run present day invading liana scenario with 150 m2 canopy
for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 400
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 1.5 * 100
    frac.tree.al = 0.9
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,1]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h[1] = K
    break
  }
}

#### Liana simulation -- present day

# Storage
liana.req_h = c()
sens = c()

# Run present day invading liana scenario with 150 m2 canopy
for(k in 1:nK){
  for(mo in 1:12){
    Ca = 400
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = 2
    
    tot.al = 1.5 * 100
    frac.liana.al = 0.1
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,1]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h[1] = K
    break
  }
}

#### Tree simulation -- future scenario

# Storage
sens = c()

# Run future invading liana scenario with 150 m2 canopy
for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 1.5 * 100
    frac.tree.al = 0.9
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,2]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h[2] = K
    break
  }
}

#### Liana simulation -- future scenario

# Storage
sens = c()

# Run future invading liana scenario with 150 m2 canopy
for(k in 1:nK){
  for(mo in 1:12){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = 2
    
    tot.al = 1.5 * 100
    frac.liana.al = 0.1
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,2]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h[2] = K
    break
  }
}

tree.req_h = tree.req_h * 0.001
liana.req_h = liana.req_h * 0.001

# Make workable dataframes
liana_kreq = cbind(c('2000', '2100'), 
                   c(liana.req_h[1], liana.req_h[2]))
colnames(liana_kreq) = c('Time', 'Kreq')

tree_kreq = cbind(c('2000', '2100'),
                  c(tree.req_h[1], tree.req_h[2]))
colnames(tree_kreq) = c('Time', 'Kreq')

# Combine and format dataframes
comb_kreq = rbind(liana_kreq, tree_kreq)
comb_kreq = as.data.frame(comb_kreq)
comb_kreq$growth.form = c('Liana', 'Liana', 
                          'Tree', 'Tree')
comb_kreq$Time = as.factor(comb_kreq$Time)
comb_kreq$Kreq = as.numeric(comb_kreq$Kreq)
comb_kreq$growth.form = as.factor(comb_kreq$growth.form)

# Plot
p1 = ggplot(comb_kreq, aes(x = Time, y = Kreq, color = growth.form, group = growth.form)) +
  geom_point(size = 4, show.legend = F) +
  geom_line(size = 1.5, show.legend = F, alpha = 0.7) +
  scale_x_discrete(limits = c('2000', '2100')) +
  xlab('') + ylab(expression(paste(K[req], ' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  scale_color_npg(name = 'PFT') +
  ggtitle(expression(paste('Total leaf area = 150', m^2))) +
  annotate('segment', x = 0.95, xend = 2.05, y = -18, yend = -18, arrow = arrow(length = unit(0.4, 'cm')), size = 1.5) +
  annotate('text', x = '2000', hjust = -0.23, y = -22, label = 'Drying hyrdoclimate', size = 7) +
  annotate('text', x = '2000', hjust = -1, y = 52, label = 'Liana', size = 7, angle = 0, color = '#E64B35FF', fontface = 2) +
  annotate('text', x = '2000', hjust = -0.45, y = 45, label = expression(bold(paste(Delta, K[req], ' = 16'))), size = 6, angle = 0, color = '#E64B35FF') +
  annotate('text', x = '2000', hjust = -2.55, y = 15, label = 'Tree', size = 7, angle = 0, color = '#4DBBD5FF', fontface = 2) +
  annotate('text', x = '2000', hjust = -1.32, y = 9, label = expression(bold(paste(Delta, K[req], ' = 2'))), size = 6, angle = 0, color = '#4DBBD5FF') +
  coord_cartesian(clip = 'off', ylim = c(-2, 60), xlim = c(1.5, 1.5)) +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26)) +
  theme(plot.title = element_text(size = 28, hjust = 0.5)) +
  theme(plot.margin = unit(c(1, 1, 3, 1), 'lines'))
p1

#########################################
## Kreq-Kw comparison: 400 m2 invasion ##
#########################################

## This section should be run in combination with the previous section

# Storage
tree.req_h_2 = c()
sens = c()

# Run present day invading liana scenario with 400 m2 canopy
for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 400
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 4 * 100
    frac.tree.al = 0.9
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,1]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h_2[1] = K
    break
  }
}

#### Liana simulation -- present day

# Storage
liana.req_h_2 = c()
sens = c()

# Run present day invading liana scenario with 400 m2 canopy
for(k in 1:nK){
  for(mo in 1:12){
    Ca = 400
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = 2
    
    tot.al = 4 * 100
    frac.liana.al = 0.1
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,1]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h_2[1] = K
    break
  }
}

#### Tree simulation -- future scenario

# Storage
sens = c()

# Run future invading liana scenario with 400 m2 canopy
for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 4 * 100
    frac.tree.al = 0.9
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,2]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h_2[2] = K
    break
  }
}

#### Liana simulation -- future scenario

# Storage
sens = c()

# Run future invading liana scenario with 400 m2 canopy
for(k in 1:nK){
  for(mo in 1:12){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = 2
    
    tot.al = 4 * 100
    frac.liana.al = 0.1
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,2]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h_2[2] = K
    break
  }
}

tree.req_h_2 = tree.req_h_2 * 0.001
liana.req_h_2 = liana.req_h_2 * 0.001

# Make workable dataframes
liana_kreq = cbind(c('2000', '2100'), 
                   c(liana.req_h_2[1], liana.req_h_2[2]))
colnames(liana_kreq) = c('Time', 'Kreq')

tree_kreq = cbind(c('2000', '2100'),
                  c(tree.req_h_2[1], tree.req_h_2[2]))
colnames(tree_kreq) = c('Time', 'Kreq')

# Combine and format dataframes
comb_kreq_2 = rbind(liana_kreq, tree_kreq)
comb_kreq_2 = as.data.frame(comb_kreq_2)
comb_kreq_2$growth.form = c('Liana', 'Liana', 
                            'Tree', 'Tree')
comb_kreq_2$Time = as.factor(comb_kreq_2$Time)
comb_kreq_2$Kreq = as.numeric(comb_kreq_2$Kreq)
comb_kreq_2$growth.form = as.factor(comb_kreq_2$growth.form)

# Plot
p2 = ggplot(comb_kreq_2, aes(x = Time, y = Kreq, color = growth.form, group = growth.form)) +
  geom_point(size = 4, show.legend = F) +
  geom_line(size = 1.5, show.legend = F, alpha = 0.7) +
  scale_x_discrete(limits = c('2000', '2100')) +
  xlab('') + ylab(expression(paste(K[req], ' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  ggtitle(expression(paste('Total leaf area = 400', m^2))) +
  scale_color_npg(name = 'PFT') +
  annotate('segment', x = 0.95, xend = 2.05, y = -35, yend = -35, arrow = arrow(length = unit(0.4, 'cm')), size = 1.5) +
  annotate('text', x = '2000', hjust = -0.23, y = -45, label = 'Drying hyrdoclimate', size = 7) +
  annotate('text', x = '2000', hjust = -1, y = 130, label = 'Liana', size = 7, angle = 0, color = '#E64B35FF', fontface = 2) +
  annotate('text', x = '2000', hjust = -0.45, y = 115, label = expression(bold(paste(Delta, K[req], ' = 41'))), size = 6, angle = 0, color = '#E64B35FF') +
  annotate('text', x = '2000', hjust = -2.55, y = 37, label = 'Tree', size = 7, angle = 0, color = '#4DBBD5FF', fontface = 2) +
  annotate('text', x = '2000', hjust = -1.32, y = 22, label = expression(bold(paste(Delta, K[req], ' = 5'))), size = 6, angle = 0, color = '#4DBBD5FF') +
  coord_cartesian(clip = 'off', ylim = c(-2, 140), xlim = c(1.5, 1.5)) +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26)) +
  theme(plot.title = element_text(size = 28, hjust = 0.5)) +
  theme(plot.margin = unit(c(1, 1, 3, 1), 'lines'))
p2

# Plot together
ga = plot_grid(p1, p2,
               labels = c('A', 'B'),
               label_size = 30)

ggsave(plot = ga, 'Plots/SupplementaryFigure15.jpeg', width = 12, height = 6)

########################
#### P50 Adaptation ####
########################

## Supplementary Figure 16

rm(list = ls())

#### Use empirical data to make P50-Slope scenarios

# Load in new data
hydro_data = read.csv('full_met_analysis_data_8_Feb.csv')

# Make units easier to work with
hydro_data$Units_K = as.character(as.factor(hydro_data$Units_K))

# Convert units for K measurements
for(i in 1:nrow(hydro_data)){
  if(hydro_data$Units_K[i] == 'kg/m/s/Mpa'){
    hydro_data$K[i] = hydro_data$K[i] * 1000 / 18 * 1000
    hydro_data$Units_K[i] = 'mmol/m/s/Mpa'
  }else{
    if(hydro_data$Units_K[i] == 'mol/m/s/Mpa'){
      hydro_data$K[i] = hydro_data$K[i] * 1000
      hydro_data$Units_K[i] = 'mmol/m/s/Mpa'
    }else{
      if(hydro_data$Units_K[i] == ''){
        next
      }else{
        print('Error with units')
      }
    }
  }
}

# Remove P50 > -0.75 consistent with Trugman et al. 2020
hydro_data = hydro_data %>%
  mutate(P50, replace(P50, P50 > -0.75, NA)) %>%
  dplyr::select(-P50)

# Rename column with NAs
colnames(hydro_data)[ncol(hydro_data)] = 'P50'

# Make sure only one unit remains for each
unique(hydro_data$Units_K)
unique(hydro_data$Units_P50)
unique(hydro_data$Units_Slope)

# Log transform variables
# Notice that P50 is now positive, but it does not change the fundamental relationship
hydro_data$log_P50 = log(-hydro_data$P50)
hydro_data$log_K = log(hydro_data$K)
hydro_data$log_slope = log(hydro_data$Slope)

# Summarize for both growth forms
all_lm = summary(lm(hydro_data$log_slope ~ hydro_data$log_P50))

#### Use coefficients from linear regression to develop new P50-Slope
#### combinations for lianas

# liana P50 = -2.25
slope_2.25 = all_lm$coefficients[1] + all_lm$coefficients[2] * log(2.25)
slope_2.25 = 10^slope_2.25

# Ensure the point falls within the range of observations
ggplot() +
  geom_point(aes(x = hydro_data$P50, y = hydro_data$Slope), alpha = 0.5) +
  geom_point(aes(x = -2.25, y = slope_2.25), color = 'maroon')

# liana P50 = -2.5
slope_2.5 = all_lm$coefficients[1] + all_lm$coefficients[2] * log(2.5)
slope_2.5 = 10^slope_2.5

# Ensure the point falls within the range of observations
ggplot() +
  geom_point(aes(x = hydro_data$P50, y = hydro_data$Slope), alpha = 0.5) +
  geom_point(aes(x = -2.5, y = slope_2.5), color = 'maroon')

# liana P50 = -3
slope_3 = all_lm$coefficients[1] + all_lm$coefficients[2] * log(3)
slope_3 = 10^slope_3

# Ensure the point falls within the range of observations
ggplot() +
  geom_point(aes(x = hydro_data$P50, y = hydro_data$Slope), alpha = 0.5) +
  geom_point(aes(x = -3, y = slope_3), color = 'maroon')

# Input parameters
load('param.input.RData')

# Met drivers
load('Driver_calc/BCI/bci_met_mxh.RData')
BCI_VPD = VPD
BCI_SWP = SWP
load('Driver_calc/Horizontes/horizontes_met_mxh.RData')
H_VPD = VPD
H_SWP = SWP
rm(VPD, SWP)

# Indexing months and days
nmonth = 12
nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

vpdsite = 11
nK = 75000

# Sequence of multipliers
mult = seq(1, 2, by = 0.1)

# BCI future scenarios
BCI_VPD_future = array(, dim = c(12, 24, 11))
BCI_VPD_future[,,1] = BCI_VPD

for(i in 2:11){
  BCI_VPD_future[,,i] = BCI_VPD * mult[i]
}

# Horizontes future scenarios
H_VPD_future = array(, dim = c(12, 24, 11))
H_VPD_future[,,1] = H_VPD

for(i in 2:11){
  H_VPD_future[,,i] = H_VPD * mult[i]
}

#### Tree - established scenario, no P50 adaptation

# Load saved simulation results from Supplementary figures.R
load('tree_req_future_established_200_wCO2.RData')

#### Tree - established scenario, P50 adaptation = -2.25

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI

sens = c()
for(k in 1:nK){
  for(mo in 1:nmonth){
    # Specify Ca = 550
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 2 * 100
    frac.tree.al = 0.6
    
    psis = BCI_SWP[mo]
    
    D = BCI_VPD_future[mo,,11]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    # Specify b1 & b2
    b1 = slope_2.25
    b2 = -2.25
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = b1, b2 = b2)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_bci_adap = K
    break
  }
}

tree.req_bci_adap = tree.req_bci_adap * 0.001

## Horizontes

for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 2 * 100
    frac.tree.al = 0.6
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,11]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    # Specify b1 & b2
    b1 = slope_2.25
    b2 = -2.25
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = b1, b2 = b2)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h_adap = K
    break
  }
}

tree.req_h_adap = tree.req_h_adap * 0.001

#### Liana - established scenario, no P50 adaptation

load('liana_req_future_established_200_wCO2.RData')

#### Liana - established scenario, P50 adaptation = -2.25

## BCI

for(k in 1:nK){
  for(mo in 1:12){
    # Specify Ca = 550
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = mean.data.liana.dbh
    
    tot.al = 2 * 100
    frac.liana.al = 0.4
    
    psis = BCI_SWP[mo]
    
    D = BCI_VPD_future[mo,,11]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    # Specify b1 & b2
    b1 = slope_2.25
    b2 = -2.25
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, b1 = b1, b2 = b2)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_bci_adap = K
    break
  }
}

liana.req_bci_adap = liana.req_bci_adap * 0.001

## Horizontes

for(k in 1:nK){
  for(mo in 1:12){
    # Specify Ca = 550
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = mean.data.liana.dbh
    
    tot.al = 2 * 100
    frac.liana.al = 0.4
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,11]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    # Specify b1 & b2
    b1 = slope_2.25
    b2 = -2.25
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, b1 = b1, b2 = b2)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h_adap = K
    break
  }
}

liana.req_h_adap = liana.req_h_adap * 0.001

#### Tree - established scenario, P50 adaptation = -2.5

## BCI

sens = c()
for(k in 1:nK){
  for(mo in 1:nmonth){
    # Specify Ca = 550
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 2 * 100
    frac.tree.al = 0.6
    
    psis = BCI_SWP[mo]
    
    D = BCI_VPD_future[mo,,11]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    # Specify b1 & b2
    b1 = slope_2.5
    b2 = -2.5
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = b1, b2 = b2)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_bci_adap_2.5 = K
    break
  }
}

tree.req_bci_adap_2.5 = tree.req_bci_adap_2.5 * 0.001

## Horizontes

for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 2 * 100
    frac.tree.al = 0.6
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,11]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    # Specify b1 & b2
    b1 = slope_2.5
    b2 = -2.5
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = b1, b2 = b2)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h_adap_2.5 = K
    break
  }
}

tree.req_h_adap_2.5 = tree.req_h_adap_2.5 * 0.001

#### Liana - established scenario, P50 adaptation = -2.5

## BCI

for(k in 1:nK){
  for(mo in 1:12){
    # Specify Ca = 550
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = mean.data.liana.dbh
    
    tot.al = 2 * 100
    frac.liana.al = 0.4
    
    psis = BCI_SWP[mo]
    
    D = BCI_VPD_future[mo,,11]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    # Specify b1 & b2
    b1 = slope_2.5
    b2 = -2.5
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, b1 = b1, b2 = b2)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_bci_adap_2.5 = K
    break
  }
}

liana.req_bci_adap_2.5 = liana.req_bci_adap_2.5 * 0.001

## Horizontes

for(k in 1:nK){
  for(mo in 1:12){
    # Specify Ca = 550
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = mean.data.liana.dbh
    
    tot.al = 2 * 100
    frac.liana.al = 0.4
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,11]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    # Specify b1 & b2
    b1 = slope_2.5
    b2 = -2.5
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, b1 = b1, b2 = b2)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h_adap_2.5 = K
    break
  }
}

liana.req_h_adap_2.5 = liana.req_h_adap_2.5 * 0.001

#### Tree - established scenario, P50 adaptation = -3

## BCI

sens = c()
for(k in 1:nK){
  for(mo in 1:nmonth){
    # Specify Ca = 550
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 2 * 100
    frac.tree.al = 0.6
    
    psis = BCI_SWP[mo]
    
    D = BCI_VPD_future[mo,,11]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    # Specify b1 & b2
    b1 = slope_3
    b2 = -3
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = b1, b2 = b2)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_bci_adap_3 = K
    break
  }
}

tree.req_bci_adap_3 = tree.req_bci_adap_3 * 0.001

## Horizontes

for(k in 1:nK){
  for(mo in 1:nmonth){
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    tree.dbh = mean.hori.data.tree.dbh
    
    tot.al = 2 * 100
    frac.tree.al = 0.6
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,11]
    
    if(mo == 1){
      tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
    }else{
      tree.height = tree.out[5]
    }
    K = Ks[k]
    # Specify b1 & b2
    b1 = slope_3
    b2 = -3
    tree.out = tree.NPP.mxh.ca(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL, b1 = b1, b2 = b2)
    sens[mo] = tree.out[1]
  }
  sum = sum(sens)
  if(sum > 0){
    tree.req_h_adap_3 = K
    break
  }
}

tree.req_h_adap_3 = tree.req_h_adap_3 * 0.001

# Plot with & without P50 adaptation for tree alone
pl1 = ggplot() +
  geom_point(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.2) +
  geom_point(aes(x = mult[11], y = tree.req_bci_adap, color = 'Wettest'), shape = 17, size = 4, show.legend = F) +
  geom_point(aes(x = mult[11], y = tree.req_h_adap, color = 'Driest'), shape = 17, size = 4, show.legend = F) +
  geom_point(aes(x = mult[11], y = tree.req_bci_adap_2.5, color = 'Wettest'), shape = 18, size = 5, show.legend = F) + 
  geom_point(aes(x = mult[11], y = tree.req_h_adap_2.5, color = 'Driest'), shape = 18, size = 5, show.legend = F) + 
  geom_point(aes(x = mult[11], y = tree.req_bci_adap_3, color = 'Wettest'), shape = 15, size = 4, show.legend = F) +
  geom_point(aes(x = mult[11], y = tree.req_h_adap_3, color = 'Driest'), shape = 15, size = 4, show.legend = F) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Tree') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pl1

#### Liana - established scenario, P50 adaptation = -3

## BCI

for(k in 1:nK){
  for(mo in 1:12){
    # Specify Ca = 550
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = mean.data.liana.dbh
    
    tot.al = 2 * 100
    frac.liana.al = 0.4
    
    psis = BCI_SWP[mo]
    
    D = BCI_VPD_future[mo,,11]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    # Specify b1 & b2
    b1 = slope_3
    b2 = -3
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, b1 = b1, b2 = b2)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_bci_adap_3 = K
    break
  }
}

liana.req_bci_adap_3 = liana.req_bci_adap_3 * 0.001

## Horizontes

for(k in 1:nK){
  for(mo in 1:12){
    # Specify Ca = 550
    Ca = 550
    nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month = mo
    
    dbh = mean.data.liana.dbh
    
    tot.al = 2 * 100
    frac.liana.al = 0.4
    
    psis = H_SWP[mo]
    
    D = H_VPD_future[mo,,11]
    
    if(mo == 1){
      liana.length = tree.out[2]
    }else{
      liana.length = lianaout[4]
    }
    K = Ks[k]
    # Specify b1 & b2
    b1 = slope_3
    b2 = -3
    lianaout = liana.NPP.mxh.ca(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL, b1 = b1, b2 = b2)
    sens[mo] = lianaout[1]
  }
  sum = sum(sens)
  if(sum > 0){
    liana.req_h_adap_3 = K
    break
  }
}

liana.req_h_adap_3 = liana.req_h_adap_3 * 0.001

pl2 = ggplot() +
  geom_point(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.2) +
  geom_point(aes(x = mult[11], y = liana.req_bci_adap, color = 'Wettest'), shape = 17, size = 4, show.legend = F) +
  geom_point(aes(x = mult[11], y = liana.req_h_adap, color = 'Driest'), shape = 17, size = 4, show.legend = F) +
  geom_point(aes(x = mult[11], y = liana.req_bci_adap_2.5, color = 'Wettest'), shape = 18, size = 5, show.legend = F) +
  geom_point(aes(x = mult[11], y = liana.req_h_adap_2.5, color = 'Driest'), shape = 18, size = 5, show.legend = F) +
  geom_point(aes(x = mult[11], y = liana.req_bci_adap_3, color = 'Wettest'), shape = 15, size = 4, show.legend = F) +
  geom_point(aes(x = mult[11], y = liana.req_h_adap_3, color = 'Driest'), shape = 15, size = 4, show.legend = F) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol ', m^-1, ' ', s^-1, ' MP', a^-1, ')'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Liana') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title = element_text(size = 28), 
        axis.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pl2

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend = get_legend(pl2)

ga = plot_grid(pl1 + theme(legend.position = 'none'), 
               pl2 + theme(legend.position = 'none'), 
               legend, nrow = 1, 
               rel_widths = c(2.3, 2.3, 0.7),
               labels = c('A', 'B', ''),
               label_size = 30)

ggsave(ga, filename = 'Plots/SupplementaryFigure16.jpeg', width = 16, height = 8, units = 'in')
