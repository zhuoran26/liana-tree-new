## Supplementary figures
## This script provides the code for all the supplementary figures

## Author: AM Willson
## Date modified: 30 March 2021

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
params_a = params[1:10,]
params_b = params[11:20,]

params_a %>%
  gt() %>%
  fmt_markdown(columns = vars('Name', 'Definition', 'Value', 'Units', 'Source', 'Tree.or.liana.function')) %>%
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
  gtsave(filename = 'Plots/Supplement_summary_parameters.png')

params_b %>%
  gt() %>%
  fmt_markdown(columns = vars('Name', 'Definition', 'Value', 'Units', 'Source', 'Tree.or.liana.function')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(90),
    heading.title.font.size = 22) %>%
  cols_label('Tree.or.liana.function' = 'Tree or Liana Function') %>%
  cols_width('Units' ~ px(150)) %>%
  cols_width('Value' ~ px(200)) %>%
  cols_width('Definition' ~ px(200)) %>%
  cols_align(align = 'center') %>%
  gtsave(filename = 'Plots/Supplement_summary_parameter_b.png')

## Data sources
source = read.csv('Sources_2.csv')

source$Source = as.character(source$Source)
source$K = as.character(source$K)
source$P50 = as.character(source$P50)
source$Slope = as.character(source$Slope)

# Split into two
source_a = source[1:8,]
source_b = source[9:15,]

source_a %>%
  gt() %>%
  fmt_markdown(columns = vars('Source')) %>%
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
             columns = vars('K', 'P50', 'Slope', 'n.liana', 'n.tree')) %>%
  cols_align(align = 'left',
             columns = vars('Source')) %>%
  cols_width('Source' ~ px(500)) %>%
  gtsave(filename = 'Plots/Supplement_sources.png')

source_b %>%
  gt() %>%
  fmt_markdown(columns = vars('Source')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(80),
    heading.title.font.size = 22) %>%
  cols_label('K' = html('K<sub><sub>s</sub></sub>')) %>%
  cols_label('P50' = html('P<sub><sub>50</sub></sub>')) %>%
  cols_label('n.liana' = html('n<sub><sub>liana</sub></sub>')) %>%
  cols_label('n.tree' = html('n<sub><sub>tree</sub></sub>')) %>%
  cols_align(align = 'center',
             columns = vars('K', 'P50', 'Slope', 'n.liana', 'n.tree')) %>%
  cols_align(align = 'left',
             columns = vars('Source')) %>%
  cols_width('Source' ~ px(500)) %>%
  gtsave(filename = 'Plots/Supplement_sources_b.png')

############################################
#### Part 2: Hydraulic trait trade-offs ####
############################################

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

# Make sure only one unit remains for each
unique(hydro_data$Units_K)
unique(hydro_data$Units_P50)
unique(hydro_data$Units_Slope)

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
  xlab('log(K) (mmol/m/s/MPa)') + ylab(expression(paste('log(-',P[50],') (MPa)'))) +
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

ggsave(pl1, filename = 'Plots/K-P50_tradeoff.jpeg', height = 8, width = 10, unit = 'in')

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
  xlab('log(K) (mmol/m/s/MPa)') + ylab('log(Slope) (%/MPa)') +
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

ggsave(pl1, filename = 'Plots/K-Slope_tradeoff.jpeg', height = 8, width = 10, unit = 'in')

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
  xlab(expression(paste('log(-',P[50],') (MPa)'))) + ylab('log(Slope) (%/MPa)') +
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

ggsave(pl1, filename = 'Plots/P50-Slope_tradeoff.jpeg', height = 8, width = 10, unit = 'in')

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

leg = get_legend(pl1)

ga = plot_grid(pl1 + theme(legend.position = 'none'), 
               pl2 + theme(legend.position = 'none'), 
               pl3 + theme(legend.position = 'none'),
               leg, rel_widths = c(1, 1), nrow = 2, ncol = 2)
ga

ggsave(ga, filename = 'Plots/hydraulic_tradeoffs.jpeg', height = 13, width = 13, unit = 'in')

###################################
#### Part 3: DBH distributions ####
###################################

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

ga = plot_grid(pl1, pl2, rel_widths = 1, nrow = 1)
ga

ggsave(ga, filename = 'Plots/DBH_distributions.jpeg', width = 20, height = 10, units = 'in')

#################################
#### Part 4: climate drivers ####
#################################

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

ggsave(pl1, filename = 'Plots/SWP.jpeg', height = 10, width = 10, units = 'in')

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

ggsave(pl2, filename = 'Plots/VPD_average.jpeg', width = 10, height = 10, units = 'in')

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

ga1 = plot_grid(pl1 + theme(legend.position = 'none'), 
                pl2, rel_widths = c(1, 1.4), nrow = 1)
ga1

ga2 = plot_grid(pl3 + theme(legend.position = 'none'), 
                pl4, rel_widths = c(1, 1.2), nrow = 1)
ga2

ggsave(ga1, filename = 'Plots/VPD_SWP_monthly.jpeg', height = 8, width = 16, units = 'in')
ggsave(ga2, filename = 'Plots/VPD_hourly.jpeg', height = 8, width = 16, units = 'in')

#####################
#### Part 5: TRY ####
#####################

rm(list = ls())

## Load data

# Read in TRY data and clean the data frame
TRY = read.csv('filtered_TRY_analysis_16-03-21.csv')
TRY = TRY[,-1]
colnames(TRY) = c('Species.ID', 'Growth.form', 'LL.months', 'SLA.mm2.mg', 'Narea.g.m2', 'Aarea.mmolC02.m2.s', 'SSD.g.cm3', 'K.kg.m.s.MPa', 'SRL.m.g', 'FRD.mm', 'CRD.mm', 'MC.%', 'Kroot.cm2.cm.s.MPa10.6', 'WVD.um', 'Dens.1.mm2', 'H.m', 'LA.cm2', 'SM.mg', 'WVEL.um', 'RRD.m', 'P50.MPa', 'Nmass.mg.g', 'Pmass.mg.g', 'Crown.depth.m', 'Species')
TRY$K.mmol.m.s.MPa = TRY$K.kg.m.s.MPa * 1000 / 18 * 1000

# Create mol/m/s/MPa column
TRY$K.mol.m.s.MPa = TRY$K.mmol.m.s.MPa * 0.001

#########################
## Plot 1: Leaf traits ##
#########################

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

plb = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, SLA.mm2.mg) %>%
  ggplot(aes(x = Growth.form, y = SLA.mm2.mg, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste('SLA (m',m^2,'/mg)'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 82', 'Tree\nn = 1207')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

plc = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Narea.g.m2) %>%
  ggplot(aes(x = Growth.form, y = Narea.g.m2, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(N[area],' (g/',m^2,')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 138', 'Tree\nn = 1716')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

pld = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Aarea.mmolC02.m2.s) %>%
  ggplot(aes(x = Growth.form, y = Aarea.mmolC02.m2.s, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(A[area],' (mmol C',O[2],'/',m^2,'/s)'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 91', 'Tree\nn = 1427')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

ple = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Nmass.mg.g) %>%
  ggplot(aes(x = Growth.form, y = Nmass.mg.g, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(N[mass],' (mg/g)'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 265', 'Tree\nn = 5264')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

plf = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Pmass.mg.g) %>%
  ggplot(aes(x = Growth.form, y = Pmass.mg.g, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(P[mass],' (mg/g)'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 163', 'Tree\nn = 3089')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

ga_leaf = grid.arrange(pla, plb, plc, pld, ple, plf, plg, nrow = 2, ncol = 4, 
                       top = textGrob('Leaf Traits', gp = gpar(fontsize = 20)))

ggsave(ga_leaf, filename = 'Plots/try_leaftraits.jpeg', width = 12, height = 7, units = 'in')

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 50))

plb = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, SLA.mm2.mg) %>%
  ggplot(aes(x = Growth.form, y = SLA.mm2.mg, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste('SLA (m',m^2,'/mg)'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 82', 'Tree\nn = 1207')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 80))

plc = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Narea.g.m2) %>%
  ggplot(aes(x = Growth.form, y = Narea.g.m2, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(N[area],' (g/',m^2,')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 138', 'Tree\nn = 1716')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 5))

pld = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Aarea.mmolC02.m2.s) %>%
  ggplot(aes(x = Growth.form, y = Aarea.mmolC02.m2.s, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(A[area],' (mmol C',O[2],'/',m^2,'/s)'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 91', 'Tree\nn = 1427')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 30))

ple = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Nmass.mg.g) %>%
  ggplot(aes(x = Growth.form, y = Nmass.mg.g, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(N[mass],' (mg/g)'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 265', 'Tree\nn = 5624')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 60))

plf = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, Pmass.mg.g) %>%
  ggplot(aes(x = Growth.form, y = Pmass.mg.g, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste(P[mass],' (mg/g)'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 163', 'Tree\nn = 3089')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 7.5))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 150))

ga_leaf = grid.arrange(pla, plb, plc, pld, ple, plf, plg, nrow = 2, ncol = 4, 
                       top = textGrob('Leaf Traits', gp = gpar(fontsize = 20)))

ggsave(ga_leaf, filename = 'Plots/try_leaftraits_trimmed_y.jpeg', width = 12, height = 7, units = 'in')

#########################
## Plot 2: Stem traits ##
#########################

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
  xlab('') + ylab(expression(paste('Stem specific density (g/c',m^3,')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 33', 'Tree\nn = 7682')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

pld = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, K.mol.m.s.MPa) %>%
  ggplot(aes(x = Growth.form, y = K.mol.m.s.MPa, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab('Conductivity (mol/m/s/MPa)') +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 18', 'Tree\nn = 400')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

ga_stem = grid.arrange(pla, plb, plc, pld, ple, nrow = 2, ncol = 3, 
                       top = textGrob('Stem & Hydraulic Traits', gp = gpar(fontsize = 20)))

ggsave(ga_stem, filename = 'Plots/try_stemtraits.jpeg', width = 10, height = 7, units = 'in')

pla = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, SSD.g.cm3) %>%
  ggplot(aes(x = Growth.form, y = SSD.g.cm3, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab(expression(paste('Stem specific density (g/c',m^3,')'))) +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 33', 'Tree\nn = 7682')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 1.25))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 300))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 250))

pld = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, K.mol.m.s.MPa) %>%
  ggplot(aes(x = Growth.form, y = K.mol.m.s.MPa, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab('Conductivity (mol/m/s/MPa)') +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 18', 'Tree\nn = 400')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 1000))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(-6, 0))

ga_stem = grid.arrange(pla, plb, plc, pld, ple, nrow = 2, ncol = 3, 
                       top = textGrob('Stem & Hydraulic Traits', gp = gpar(fontsize = 20)))

ggsave(ga_stem, filename = 'Plots/try_stemtraits_trimmed_y.jpeg', width = 10, height = 7, units = 'in')

#########################
## Plot 3: Root traits ##
#########################

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
  xlab('') + ylab('Specific root length (m/g)') +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 15', 'Tree\nn = 179')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

ga_root = grid.arrange(pla, plb, plc, pld, nrow = 2, ncol = 2,
                       top = textGrob('Root Traits', gp = gpar(fontsize = 20)))

ggsave(ga_root, filename = 'Plots/try_roottraits.jpeg', width = 10, height = 7, units = 'in')

pla = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  dplyr::select(Growth.form, SRL.m.g) %>%
  ggplot(aes(x = Growth.form, y = SRL.m.g, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab('Specific root length (m/g)') +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 15', 'Tree\nn = 179')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 250))

plb = TRY %>%
  filter(Growth.form %in% c('tree','liana')) %>%
  select(Growth.form, FRD.mm) %>%
  ggplot(aes(x = Growth.form, y = FRD.mm, fill = Growth.form)) +
  geom_violin(alpha = 0.5, show.legend = F) +
  stat_summary(fun = median, geom = 'point', show.legend = F, shape = 3, size = 2, stroke = 1.2) +
  theme_minimal() +
  xlab('') + ylab('Fine root diameter (mm)') +
  scale_fill_npg(name = 'PFT') +
  scale_x_discrete(labels = c('Liana\nn = 13', 'Tree\nn = 169')) +
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 1.25))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15))

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
  theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 10))

ga_root = grid.arrange(pla, plb, plc, pld, nrow = 2, ncol = 2,
                       top = textGrob('Root Traits', gp = gpar(fontsize = 20)))

ggsave(ga_root, filename = 'Plots/try_roottraits_trimmed_y.jpeg', width = 10, height = 7, units = 'in')

######################
## Table 1: U tests ##
######################

sig_tests = as_tibble(matrix(, nrow = 16, ncol = 7))
colnames(sig_tests) = c('Trait', 'Tree mean', 'Liana mean', 'n tree', 'n liana', 'Test Statistic', 'p-value')

sig_tests$Trait[1] = 'Stem specific density (g/cm<sup>3</sup>)'
sig_tests$`Tree mean`[1] = mean(TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[1] = mean(TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[1] = length(TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree' & !is.na(TRY$SSD.g.cm3))])
sig_tests$`n liana`[1] = length(TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana' & !is.na(TRY$SSD.g.cm3))])
sig_tests$`Test Statistic`[1] = wilcox.test(TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')],
                                            TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[1] = wilcox.test(TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')],
                                     TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[2] = 'Vessel diameter (&mu;m)'
sig_tests$`Tree mean`[2] = mean(TRY$WVD.um[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[2] = mean(TRY$WVD.um[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[2] = length(TRY$WVD.um[which(TRY$Growth.form == 'tree' & !is.na(TRY$WVD.um))])
sig_tests$`n liana`[2] = length(TRY$WVD.um[which(TRY$Growth.form == 'liana' & !is.na(TRY$WVD.um))])
sig_tests$`Test Statistic`[2] = wilcox.test(TRY$WVD.um[which(TRY$Growth.form == 'tree')],
                                            TRY$WVD.um[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[2] = wilcox.test(TRY$WVD.um[which(TRY$Growth.form == 'tree')],
                                     TRY$WVD.um[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[3] = 'Vessel density (1/mm<sup>2</sup>)'
sig_tests$`Tree mean`[3] = mean(TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[3] = mean(TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[3] = length(TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree' & !is.na(TRY$Dens.1.mm2))])
sig_tests$`n liana`[3] = length(TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana' & !is.na(TRY$Dens.1.mm2))])
sig_tests$`Test Statistic`[3] = wilcox.test(TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')],
                                            TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[3] = wilcox.test(TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')],
                                     TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[4] = 'Stem specific hydraulic conductivity (mol/m/s/MPa)'
sig_tests$`Tree mean`[4] = mean(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[4] = mean(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[4] = length(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree' & !is.na(TRY$K.mol.m.s.MPa))])
sig_tests$`n liana`[4] = length(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana' & !is.na(TRY$K.mol.m.s.MPa))])
sig_tests$`Test Statistic`[4] = wilcox.test(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')],
                                            TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[4] = wilcox.test(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')],
                                     TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[5] = 'P<sub>50</sub> (MPa)'
sig_tests$`Tree mean`[5] = mean(TRY$P50.MPa[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[5] = mean(TRY$P50.MPa[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[5] = length(TRY$P50.MPa[which(TRY$Growth.form == 'tree' & !is.na(TRY$P50.MPa))])
sig_tests$`n liana`[5] = length(TRY$P50.MPa[which(TRY$Growth.form == 'liana' & !is.na(TRY$P50.MPa))])
sig_tests$`Test Statistic`[5] = wilcox.test(TRY$P50.MPa[which(TRY$Growth.form == 'tree')],
                                            TRY$P50.MPa[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[5] = wilcox.test(TRY$P50.MPa[which(TRY$Growth.form == 'tree')],
                                     TRY$P50.MPa[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[6] = 'Leaf lifespan (months)'
sig_tests$`Tree mean`[6] = mean(TRY$LL.months[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[6] = mean(TRY$LL.months[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[6] = length(TRY$LL.months[which(TRY$Growth.form == 'tree' & !is.na(TRY$LL.months))])
sig_tests$`n liana`[6] = length(TRY$LL.months[which(TRY$Growth.form == 'liana' & !is.na(TRY$LL.months))])
sig_tests$`Test Statistic`[6] = wilcox.test(TRY$LL.months[which(TRY$Growth.form == 'tree')],
                                            TRY$LL.months[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[6] = wilcox.test(TRY$LL.months[which(TRY$Growth.form == 'tree')],
                                     TRY$LL.months[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[7] = 'Specific leaf area (mm<sup>2</sup>/mg)'
sig_tests$`Tree mean`[7] = mean(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[7] = mean(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[7] = length(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree' & !is.na(TRY$SLA.mm2.mg))])
sig_tests$`n liana`[7] = length(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana' & !is.na(TRY$SLA.mm2.mg))])
sig_tests$`Test Statistic`[7] = wilcox.test(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')],
                                            TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[7] = wilcox.test(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')],
                                     TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[8] = 'Area-based leaf nitrogen (g/m<sup>2</sup>)'
sig_tests$`Tree mean`[8] = mean(TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[8] = mean(TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[8] = length(TRY$Narea.g.m2[which(TRY$Growth.form == 'tree' & !is.na(TRY$Narea.g.m2))])
sig_tests$`n liana`[8] = length(TRY$Narea.g.m2[which(TRY$Growth.form == 'liana' & !is.na(TRY$Narea.g.m2))])
sig_tests$`Test Statistic`[8] = wilcox.test(TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')],
                                            TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[8] = wilcox.test(TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')],
                                     TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[9] = 'Area-based photosynthetic rate (mmol CO<sub>2</sub>/m<sup>2</sup>/s)'
sig_tests$`Tree mean`[9] = mean(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[9] = mean(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[9] = length(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree' & !is.na(TRY$Aarea.mmolC02.m2.s))])
sig_tests$`n liana`[9] = length(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana' & !is.na(TRY$Aarea.mmolC02.m2.s))])
sig_tests$`Test Statistic`[9] = wilcox.test(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')],
                                            TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[9] = wilcox.test(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')],
                                     TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[10] = 'Mass-based leaf nitrogen (mg/g)'
sig_tests$`Tree mean`[10] = mean(TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[10] = mean(TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[10] = length(TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree' & !is.na(TRY$Nmass.mg.g))])
sig_tests$`n liana`[10] = length(TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana' & !is.na(TRY$Nmass.mg.g))])
sig_tests$`Test Statistic`[10] = wilcox.test(TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')],
                                             TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[10] = wilcox.test(TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')],
                                      TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[11] = 'Mass-based leaf phosophorus (mg/g)'
sig_tests$`Tree mean`[11] = mean(TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[11] = mean(TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[11] = length(TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree' & !is.na(TRY$Pmass.mg.g))])
sig_tests$`n liana`[11] = length(TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana' & !is.na(TRY$Pmass.mg.g))])
sig_tests$`Test Statistic`[11] = wilcox.test(TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')],
                                             TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[11] = wilcox.test(TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')],
                                      TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[12] = 'Leaf area (cm<sup>2</sup>)'
sig_tests$`Tree mean`[12] = mean(TRY$LA.cm2[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[12] = mean(TRY$LA.cm2[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[12] = length(TRY$LA.cm2[which(TRY$Growth.form == 'tree' & !is.na(TRY$LA.cm2))])
sig_tests$`n liana`[12] = length(TRY$LA.cm2[which(TRY$Growth.form == 'liana' & !is.na(TRY$LA.cm2))])
sig_tests$`Test Statistic`[12] = wilcox.test(TRY$LA.cm2[which(TRY$Growth.form == 'tree')], 
                                             TRY$LA.cm2[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[12] = wilcox.test(TRY$LA.cm2[which(TRY$Growth.form == 'tree')],
                                      TRY$LA.cm2[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[13] = 'Specific root length (m/g)'
sig_tests$`Tree mean`[13] = mean(TRY$SRL.m.g[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[13] = mean(TRY$SRL.m.g[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[13] = length(TRY$SRL.m.g[which(TRY$Growth.form == 'tree' & !is.na(TRY$SRL.m.g))])
sig_tests$`n liana`[13] = length(TRY$SRL.m.g[which(TRY$Growth.form == 'liana' & !is.na(TRY$SRL.m.g))])
sig_tests$`Test Statistic`[13] = wilcox.test(TRY$SRL.m.g[which(TRY$Growth.form == 'tree')],
                                             TRY$SRL.m.g[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[13] = wilcox.test(TRY$SRL.m.g[which(TRY$Growth.form == 'tree')],
                                      TRY$SRL.m.g[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[14] = 'Fine root diameter (mm)'
sig_tests$`Tree mean`[14] = mean(TRY$FRD.mm[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[14] = mean(TRY$FRD.mm[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[14] = length(TRY$FRD.mm[which(TRY$Growth.form == 'tree' & !is.na(TRY$FRD.mm))])
sig_tests$`n liana`[14] = length(TRY$FRD.mm[which(TRY$Growth.form == 'liana' & !is.na(TRY$FRD.mm))])
sig_tests$`Test Statistic`[14] = wilcox.test(TRY$FRD.mm[which(TRY$Growth.form == 'tree')],
                                             TRY$FRD.mm[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[14] = wilcox.test(TRY$FRD.mm[which(TRY$Growth.form == 'tree')],
                                      TRY$FRD.mm[which(TRY$Growth.form == 'liana')])$p.value

sig_tests$Trait[15] = 'Mycorrhizal colonization (%)'
sig_tests$`Tree mean`[15] = mean(TRY$`MC.%`[which(TRY$Growth.form == 'tree')], na.rm = T)
sig_tests$`Liana mean`[15] = mean(TRY$`MC.%`[which(TRY$Growth.form == 'liana')], na.rm = T)
sig_tests$`n tree`[15] = length(TRY$`MC.%`[which(TRY$Growth.form == 'tree' & !is.na(TRY$`MC.%`))])
sig_tests$`n liana`[15] = length(TRY$`MC.%`[which(TRY$Growth.form == 'liana' & !is.na(TRY$`MC.%`))])
sig_tests$`Test Statistic`[15] = wilcox.test(TRY$`MC.%`[which(TRY$Growth.form == 'tree')],
                                             TRY$`MC.%`[which(TRY$Growth.form == 'liana')])$statistic
sig_tests$`p-value`[15] = wilcox.test(TRY$`MC.%`[which(TRY$Growth.form == 'tree')],
                                      TRY$`MC.%`[which(TRY$Growth.form == 'liana')])$p.value

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
  gtsave(filename = 'Plots/TRY_mann-whitney_tests.png')

###########################
## Table 2: Effect sizes ##
###########################

effect_size = as_tibble(matrix(, nrow = 16, ncol = 4))
colnames(effect_size) = c('Trait', 'Glass Delta', 'Lower CI', 'Upper CI')

effect_size$Trait[1] = 'Stem specific density (g/cm<sup>3</sup>)'
effect_size$`Glass Delta`[1] = glass_delta(TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')], TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[1] = glass_delta(TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')], TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[1] = glass_delta(TRY$SSD.g.cm3[which(TRY$Growth.form == 'liana')], TRY$SSD.g.cm3[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[2] = 'Vessel diameter (&mu;m)'
effect_size$`Glass Delta`[2] = glass_delta(TRY$WVD.um[which(TRY$Growth.form == 'liana')], TRY$WVD.um[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[2] = glass_delta(TRY$WVD.um[which(TRY$Growth.form == 'liana')], TRY$WVD.um[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[2] = glass_delta(TRY$WVD.um[which(TRY$Growth.form == 'liana')], TRY$WVD.um[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[3] = 'Vessel density (1/mm<sup>2</sup>)'
effect_size$`Glass Delta`[3] = glass_delta(TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')], TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[3] = glass_delta(TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')], TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[3] = glass_delta(TRY$Dens.1.mm2[which(TRY$Growth.form == 'liana')], TRY$Dens.1.mm2[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[4] = 'Stem specific hydraulic conductivity (mol/m/s/MPa)'
effect_size$`Glass Delta`[4] = glass_delta(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')], TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[4] = glass_delta(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')], TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[4] = glass_delta(TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'liana')], TRY$K.mol.m.s.MPa[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[5] = 'P<sub>50</sub> (MPa)'
effect_size$`Glass Delta`[5] = glass_delta(TRY$P50.MPa[which(TRY$Growth.form == 'liana')], TRY$P50.MPa[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[5] = glass_delta(TRY$P50.MPa[which(TRY$Growth.form == 'liana')], TRY$P50.MPa[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[5] = glass_delta(TRY$P50.MPa[which(TRY$Growth.form == 'liana')], TRY$P50.MPa[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[6] = 'Leaf lifespan (months)'
effect_size$`Glass Delta`[6] = glass_delta(TRY$LL.months[which(TRY$Growth.form == 'liana')], TRY$LL.months[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[6] = glass_delta(TRY$LL.months[which(TRY$Growth.form == 'liana')], TRY$LL.months[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[6] = glass_delta(TRY$LL.months[which(TRY$Growth.form == 'liana')], TRY$LL.months[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[7] = 'Specific leaf area (mm<sup>2</sup>/mg)'
effect_size$`Glass Delta`[7] = glass_delta(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')], TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[7] = glass_delta(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')], TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[7] = glass_delta(TRY$SLA.mm2.mg[which(TRY$Growth.form == 'liana')], TRY$SLA.mm2.mg[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[8] = 'Area-based leaf nitrogen (g/m<sup>2</sup>)'
effect_size$`Glass Delta`[8] = glass_delta(TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')], TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[8] = glass_delta(TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')], TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[8] = glass_delta(TRY$Narea.g.m2[which(TRY$Growth.form == 'liana')], TRY$Narea.g.m2[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[9] = 'Area-based photosynthetic rate (mmol CO<sub>2</sub>/m<sup>2</sup>/s)'
effect_size$`Glass Delta`[9] = glass_delta(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')], TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[9] = glass_delta(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')], TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[9] = glass_delta(TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'liana')], TRY$Aarea.mmolC02.m2.s[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[10] = 'Mass-based leaf nitrogen (mg/g)'
effect_size$`Glass Delta`[10] = glass_delta(TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')], TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[10] = glass_delta(TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')], TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[10] = glass_delta(TRY$Nmass.mg.g[which(TRY$Growth.form == 'liana')], TRY$Nmass.mg.g[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[11] = 'Mass-based leaf phosphorus (mg/g)'
effect_size$`Glass Delta`[11] = glass_delta(TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')], TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[11] = glass_delta(TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')], TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[11] = glass_delta(TRY$Pmass.mg.g[which(TRY$Growth.form == 'liana')], TRY$Pmass.mg.g[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[12] = 'Leaf area (cm<sup>2</sup>)'
effect_size$`Glass Delta`[12] = glass_delta(TRY$LA.cm2[which(TRY$Growth.form == 'liana')], TRY$LA.cm2[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[12] = glass_delta(TRY$LA.cm2[which(TRY$Growth.form == 'liana')], TRY$LA.cm2[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[12] = glass_delta(TRY$LA.cm2[which(TRY$Growth.form == 'liana')], TRY$LA.cm2[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[13] = 'Specific root length (m/g)'
effect_size$`Glass Delta`[13] = glass_delta(TRY$SRL.m.g[which(TRY$Growth.form == 'liana')], TRY$SRL.m.g[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[13] = glass_delta(TRY$SRL.m.g[which(TRY$Growth.form == 'liana')], TRY$SRL.m.g[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[13] = glass_delta(TRY$SRL.m.g[which(TRY$Growth.form == 'liana')], TRY$SRL.m.g[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[14] = 'Fine root diameter (mm)'
effect_size$`Glass Delta`[14] = glass_delta(TRY$FRD.mm[which(TRY$Growth.form == 'liana')], TRY$FRD.mm[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[14] = glass_delta(TRY$FRD.mm[which(TRY$Growth.form == 'liana')], TRY$FRD.mm[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[14] = glass_delta(TRY$FRD.mm[which(TRY$Growth.form == 'liana')], TRY$FRD.mm[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[15] = 'Mycorrhizal colonization (%)'
effect_size$`Glass Delta`[15] = glass_delta(TRY$`MC.%`[which(TRY$Growth.form == 'liana')], TRY$`MC.%`[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[15] = glass_delta(TRY$`MC.%`[which(TRY$Growth.form == 'liana')], TRY$`MC.%`[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[15] = glass_delta(TRY$`MC.%`[which(TRY$Growth.form == 'liana')], TRY$`MC.%`[which(TRY$Growth.form == 'tree')])$CI_high

effect_size$Trait[16] = 'Rooting depth (m)'
effect_size$`Glass Delta`[16] = glass_delta(TRY$RRD.m[which(TRY$Growth.form == 'liana')], TRY$RRD.m[which(TRY$Growth.form == 'tree')])$Glass_delta
effect_size$`Lower CI`[16] = glass_delta(TRY$RRD.m[which(TRY$Growth.form == 'liana')], TRY$RRD.m[which(TRY$Growth.form == 'tree')])$CI_low
effect_size$`Upper CI`[16] = glass_delta(TRY$RRD.m[which(TRY$Growth.form == 'liana')], TRY$RRD.m[which(TRY$Growth.form == 'tree')])$CI_high

effect_size %>% gt() %>%
  fmt_markdown(columns = vars('Trait')) %>%
  fmt_number(columns = vars('Glass Delta', 'Lower CI', 'Upper CI'), decimals = 2) %>%
  cols_label('Glass Delta' = html("Glass' &Delta;")) %>%
  tab_header(
    title = md('Effect Size for TRY traits')) %>%
  tab_options(
    column_labels.font.size = 20,
    table.width = pct(90),
    heading.title.font.size = 22) %>%
  tab_style(
    style = cell_fill(color = 'lightyellow'),
    locations = cells_body(rows = c(1:2, 4:5, 7:8, 11:12))) %>%
  cols_align(align = 'center',
             columns = c('Glass Delta', 'Lower CI', 'Upper CI')) %>%
  cols_align(align = 'left',
             columns = 'Trait') %>%
  gtsave(filename = 'Plots/TRY_effect_size_tests.png')

########################################
#### Part 6: Extended Meta-analysis ####
########################################

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
  gtsave(filename = 'Plots/extended_mann-whitney_tests.png')

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
  fmt_markdown(columns = vars('Trait')) %>%
  fmt_number(columns = vars('Glass Delta', 'Lower CI', 'Upper CI'), decimals = 2) %>%
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
  gtsave(filename = 'Plots/extended_effect_size_tests.png')

######################
#### Part 7: Kreq ####
######################

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

rownames(liana_out_mat) = dbhs
colnames(liana_out_mat) = varyk

liana_melt = melt(liana_out_mat)
colnames(liana_melt) = c('DBH', 'Ks', 'NPP')
liana_melt$Ks = liana_melt$Ks * 0.001

plot = liana_melt %>%
  ggplot(aes(x = Ks, y = NPP, group = as.factor(DBH), color = as.factor(DBH))) +
  geom_line(size = 1.2) +
  xlab(expression(K[s],' (mol/m/s/MPa)')) + ylab('NPP (kg C/year)') +
  geom_vline(aes(xintercept = liana_melt$Ks[which(liana_melt$DBH == 2) & liana_melt$NPP == min(abs(liana_melt$NPP[which(liana_melt$DBH == 2)]))], linetype = 'Kreq')) +
  geom_vline(aes(xintercept = liana_melt$Ks[which(liana_melt$DBH == 6) & liana_melt$NPP == min(abs(liana_melt$NPP[which(liana_melt$DBH == 6)]))], linetype = 'Kreq')) +
  scale_linetype_manual(name = element_blank(), values = c('Kreq' = 'dashed'), labels = expression(K[req])) +
  scale_color_npg(name = 'DBH (cm)') +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 25), legend.text = element_text(size = 25), legend.title = element_text(size = 28)) +
  guides(color = guide_legend(reverse = T)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

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

plot.w.inset = ggdraw() +
  draw_plot(plot) +
  draw_plot(inset, x = 0.4, y = 0.15, width = 0.4, height = 0.4) +
  theme(axis.text = element_text(size = 18))

ggsave(plot = plot.w.inset, filename = 'Plots/Kreq_explain.jpeg', width = 12, height = 7, units = 'in')

####################################################
#### Part 8: Invading liana climate sensitivity ####
####################################################

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Interpolation_mxh_100.RData')

nsite = 100
nmonth = 12
nK = 75000

nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

VPDs = VPD_interp_100
SWPs = SWP_interp_100

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

tree.req = matrix(, nrow = nsite, ncol = nsite)

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

#save(tree.req, file = 'tree.req.100interp.monthly_invasion_200.RData')

load(file = 'tree.req.100interp.monthly_invasion_200.RData')
tree.req = tree.req * 0.001
tree_req_melt = melt(tree.req)
colnames(tree_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

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

ggsave(plot = ra1, file = 'Plots/tree.sens.ksurv.raster_final_invasion_200.jpeg')

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)
Ks = Ks[Ks > 17550]

liana.req = matrix(, nrow = nsite, ncol = nsite)

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

#save(liana.req, file = 'liana.req.100interp.monthly_invasion_200.RData')

load(file = 'liana.req.100interp.monthly_invasion_200.RData')
liana.req = liana.req * 0.001
liana_req_melt = melt(liana.req)
colnames(liana_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

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

ggsave(ra2, file = 'Plots/liana.sens.ksurv.raster_final_invasion_200.jpeg')

ga_fin = plot_grid(ra1, ra2, nrow = 1, rel_widths = 1)
ggsave(ga_fin, filename = 'Plots/tree_liana_relsens_concept_final_invasion_200.jpeg', height = 12, width = 20, units = 'in')

##################################
#### Part 9: Future scenarios ####
##################################

rm(list = ls())

# Input parameters
load('param.input.RData')

# Rerun model with Ca as an input (instead of fixed in model)
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, Ca){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = Ca # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8*10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, Ca){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = Ca # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8*10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

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

#####################################
## Tree runs, established scenario ##
#####################################

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI
tree.req_bci = c()

count = 0

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
      tree.out = tree.NPP.mxh(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
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
tree.req_h = c()

count = 0

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
      tree.out = tree.NPP.mxh(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
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

save(tree.req_bci, tree.req_h, file = 'tree_req_future_established_200_wCO2.RData')

pl1 = ggplot() +
  geom_point(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol/m/s/MPa)'))) +
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

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI
liana.req_bci = c()

count = 0

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
      lianaout = liana.NPP.mxh(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
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
liana.req_h = c()

count = 0

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
      lianaout = liana.NPP.mxh(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
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

save(liana.req_bci, liana.req_h, file = 'liana_req_future_established_200_wCO2.RData')

pl2 = ggplot() +
  geom_point(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol/m/s/MPa)'))) +
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

ga = grid.arrange(pl1 + theme(legend.position = 'none'), 
                  pl2 + theme(legend.position = 'none'), 
                  legend, nrow = 1, widths = c(2.3, 2.3, 0.7))

ggsave(ga, filename = 'Plots/both_established_fugure_200_wCO2.jpeg', width = 16, height = 8, units = 'in')

##################################
## Tree runs, invasion scenario ##
##################################

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI
tree.req_bci = c()

count = 0

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
      tree.out = tree.NPP.mxh(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
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
tree.req_h = c()

count = 0

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
      tree.out = tree.NPP.mxh(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
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

save(tree.req_bci, tree.req_h, file = 'tree_req_future_invasion_200_wCO2.RData')

pl1 = ggplot() +
  geom_point(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol/m/s/MPa)'))) +
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
ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI
liana.req_bci = c()

count = 0

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
      lianaout = liana.NPP.mxh(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
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
liana.req_h = c()

count = 0

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
      lianaout = liana.NPP.mxh(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
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

save(liana.req_bci, liana.req_h, file = 'liana_req_future_invasion_200_wCO2.RData')

pl2 = ggplot() +
  geom_point(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol/m/s/MPa)'))) +
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

legend = get_legend(pl1)

ga = grid.arrange(pl1 + theme(legend.position = 'none'), 
                  pl2 + theme(legend.position = 'none'),
                  legend, nrow = 1, widths = c(2.3, 2.3, 0.7))

ggsave(ga, filename = 'Plots/both_invasion_fugure_200_wCO2.jpeg', width = 16, height = 8, units = 'in')

#####################################
#### Part 10: Kreq-Kw comparison ####
#####################################

rm(list = ls())

# Load parameters and models
load('param.input.RData')

# Load future Kreq for invasion scenario
# These were created for plots in the Supplement
load('tree_req_future_invasion_200_wCO2.RData')
load('liana_req_future_invasion_200_wCO2.RData')

# Convert Ks to mol/m/s/MPa
tree_K = trees_K$K * 0.001
liana_K = lianas_K$K * 0.001

# Convert Ks to Kw by dividing by 10
tree_K = tree_K / 10
liana_K = liana_K / 10

# Make histograms
pl1 = ggplot() +
  geom_histogram(aes(x = tree_K), bins = 20) +
  geom_vline(aes(xintercept = min(tree.req_bci), color = 'Wettest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(tree.req_bci), color = 'Wettest', linetype = 'Future'), size = 1.2) +
  geom_vline(aes(xintercept = min(tree.req_h), color = 'Driest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(tree.req_h), color = 'Driest', linetype = 'Future'), size = 1.2) +
  xlab('K (mol/m/s/MPa)') + ylab('Frequency') +
  scale_color_discrete(name = 'Hydroclimate') +
  scale_linetype_manual(name = 'Time Period', values = c('Present' = 'solid', 'Future' = 'dashed'), breaks = c('Present', 'Future')) +
  ggtitle('Tree') +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl1

pl2 = ggplot() +
  geom_histogram(aes(x = liana_K), bins = 20) +
  geom_vline(aes(xintercept = min(liana.req_bci), color = 'Wettest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(liana.req_bci), color = 'Wettest', linetype = 'Future'), size = 1.2) +
  geom_vline(aes(xintercept = min(liana.req_h), color = 'Driest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(liana.req_h), color = 'Driest', linetype = 'Future'), size = 1.2) +
  xlab('K (mol/m/s/MPa)') + ylab('Frequency') +
  scale_color_discrete(name = 'Hydroclimate') +
  scale_linetype_manual(name = 'Time Period', values = c('Present' = 'solid', 'Future' = 'dashed'), breaks = c('Present', 'Future')) +
  ggtitle('Liana') +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, 'cm'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl2

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend = get_legend(pl1)

ga = plot_grid(pl1 + theme(legend.position = 'none'), 
               pl2 + theme(legend.position = 'none'),
               legend, nrow = 1, rel_widths = c(1, 1, 0.5))
ga

ggsave(ga, filename = 'Plots/Kreq_Kw_comp_invasion_200_wCO2.jpeg', width = 16, height = 10, units = 'in')

##########################################
#### Part 11: Canopy area sensitivity ####
##########################################

##########################
## Main Figure 2 150 m2 ##
##########################

rm(list = ls())

# Load parameters
load('param.input.RData')

# Load met drivers
load('bci_met_mxh.RData')

nK = 75000
ks = c(trees_K$K / 10, lianas_K$K)
varyk = seq(min(ks), max(ks), length.out = nK)

# Interpolate DBH (and remove highest DBH)
hori.data.liana.dbh = hori.data.liana.dbh[hori.data.liana.dbh < 29]
dbhs = seq(min(hori.data.liana.dbh), max(hori.data.liana.dbh), length.out = 50)

## Tree: 90% canopy, mean DBH, BCI

treeouteach = c()
tree_out_r1 = c()

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

lianaouteach = c()
liana_out_r1 = matrix(, nrow = length(dbhs), ncol = 2)

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

treeouteach = c()
tree_out_r2 = c()

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

lianaouteach = c()
liana_out_r2 = matrix(, nrow = length(dbhs), ncol = 2)

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
treeouteach = c()
tree_out_r3 = c()

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

lianaouteach = c()
liana_out_r3 = matrix(, nrow = length(dbhs), ncol = 2)

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

treeouteach = c()
tree_out_r4 = c()

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

lianaouteach = c()
liana_out_r4 = matrix(, nrow = length(dbhs), ncol = 2)

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

df1 = as.data.frame(cbind(df1, rep(1, nrow(df1))))
colnames(df1)[4] = 'df'
df2 = as.data.frame(cbind(df2, rep(2, nrow(df2))))
colnames(df2)[4] = 'df'
df3 = as.data.frame(cbind(df3, rep(3, nrow(df3))))
colnames(df3)[4] = 'df'
df4 = as.data.frame(cbind(df4, rep(4, nrow(df4))))
colnames(df4)[4] = 'df'
df = as.data.frame(rbind(df1, df2, df3, df4))

pl = ggplot(df, aes(x = DBH, y = log(Kreq), group = df)) +
  geom_line(data = subset(df, PFT == 'Liana'), aes(color = as.factor(df)), size = 1.2) +
  geom_hline(data = subset(df, PFT == 'Tree'), aes(yintercept = log(Kreq), color = as.factor(df)), size = 1, linetype = 'dashed') +
  xlab('DBH (cm)') + ylab(expression(paste('log(',K[req],')',' (mol/m/s/MPa)'))) +
  #scale_color_viridis_d(name = 'Scenario', end = 0.9, labels = c('1' = '\n10% liana\ntropical moist', '2' = '\n40% liana\ntropical moist\n', '3' = '\n10% liana\ntropical dry\n', '4' = '\n40% liana\ntropical dry\n'), breaks = c('4', '2', '3', '1')) +
  scale_color_manual(name = 'Scenario', 
                     labels = c('1' = '\n10% liana\ntropical moist', '2' = '\n40% liana\ntropical moist\n', '3' = '\n10% liana\ntropical dry\n', '4' = '\n40% liana\ntropical dry\n'), 
                     breaks = c('4', '2', '3', '1'),
                     values = c('4' = '#a6611a', '2' = '#80cdc1', '3' = '#dfc27d', '1' = '#018571')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.text = element_text(size = 26), legend.title = element_text(size = 28, hjust = 0.5), panel.grid = element_blank())

ggsave(pl, filename = 'Plots/allom_hydro_comp_combined_log_150_dashed.jpeg', width = 14, height = 8, units = 'in')

##########################
## Main Figure 2 400 m2 ##
##########################

rm(list = ls())

# Load parameters
load('param.input.RData')

# Load met drivers
load('bci_met_mxh.RData')

nK = 75000
ks = c(trees_K$K / 10, lianas_K$K)
varyk = seq(min(ks), max(ks), length.out = nK)

# Interpolate DBH (and remove highest DBH)
hori.data.liana.dbh = hori.data.liana.dbh[hori.data.liana.dbh < 29]
dbhs = seq(min(hori.data.liana.dbh), max(hori.data.liana.dbh), length.out = 50)

## Tree: 90% canopy, mean DBH, BCI

treeouteach = c()
tree_out_r1 = c()

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

lianaouteach = c()
liana_out_r1 = matrix(, nrow = length(dbhs), ncol = 2)

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

treeouteach = c()
tree_out_r2 = c()

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

lianaouteach = c()
liana_out_r2 = matrix(, nrow = length(dbhs), ncol = 2)

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
treeouteach = c()
tree_out_r3 = c()

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

lianaouteach = c()
liana_out_r3 = matrix(, nrow = length(dbhs), ncol = 2)

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

treeouteach = c()
tree_out_r4 = c()

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

lianaouteach = c()
liana_out_r4 = matrix(, nrow = length(dbhs), ncol = 2)

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

df1 = as.data.frame(cbind(df1, rep(1, nrow(df1))))
colnames(df1)[4] = 'df'
df2 = as.data.frame(cbind(df2, rep(2, nrow(df2))))
colnames(df2)[4] = 'df'
df3 = as.data.frame(cbind(df3, rep(3, nrow(df3))))
colnames(df3)[4] = 'df'
df4 = as.data.frame(cbind(df4, rep(4, nrow(df4))))
colnames(df4)[4] = 'df'
df = as.data.frame(rbind(df1, df2, df3, df4))

pl = ggplot(df, aes(x = DBH, y = log(Kreq), group = df)) +
  geom_line(data = subset(df, PFT == 'Liana'), aes(color = as.factor(df)), size = 1.2) +
  geom_hline(data = subset(df, PFT == 'Tree'), aes(yintercept = log(Kreq), color = as.factor(df)), size = 1, linetype = 'dashed') +
  xlab('DBH (cm)') + ylab(expression(paste('log(',K[req],')',' (mol/m/s/MPa)'))) +
  scale_color_manual(name = 'Scenario', 
                     labels = c('1' = '\n10% liana\ntropical moist', '2' = '\n40% liana\ntropical moist\n', '3' = '\n10% liana\ntropical dry\n', '4' = '\n40% liana\ntropical dry\n'), 
                     breaks = c('4', '2', '3', '1'),
                     values = c('4' = '#a6611a', '2' = '#80cdc1', '3' = '#dfc27d', '1' = '#018571')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.text = element_text(size = 26), legend.title = element_text(size = 28, hjust = 0.5), panel.grid = element_blank())

ggsave(pl, filename = 'Plots/allom_hydro_comp_combined_log_400_dashed.jpeg', width = 14, height = 8, units = 'in')

#########################
## Interpolation plots ##
#########################

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Interpolation_mxh_100.RData')

nsite = 100
nmonth = 12
nK = 75000

nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

VPDs = VPD_interp_100
SWPs = SWP_interp_100

#########################
## Established: 150 m2 ##
#########################

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

tree.req = matrix(, nrow = nsite, ncol = nsite)

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

#save(tree.req, file = 'tree.req.100interp.monthly_150.RData')

load(file = 'tree.req.100interp.monthly_150.RData')
tree.req = tree.req * 0.001
tree_req_melt = melt(tree.req)
colnames(tree_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

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

ks = c(trees_K$K / 10, lianas_K$K)
nK = 5000
Ks = seq(min(ks), max(ks), length.out = nK)

liana.req = matrix(, nrow = nsite, ncol = nsite)

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

#save(liana.req, file = 'liana.req.100interp.monthly_150.RData')

load(file = 'liana.req.100interp.monthly_150.RData')
liana.req = liana.req * 0.001
liana_req_melt = melt(liana.req)
colnames(liana_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

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

ra1_fin = ra1 + theme(legend.position = 'bottom', plot.margin = unit(c(0, 0, 0, 0), 'cm'))
ra2_fin = ra2 + theme(legend.position = 'bottom', plot.margin = unit(c(-1, 0, 1, 0), 'cm'))

ga = plot_grid(ra1_fin, ra2_fin, align = 'hv', axis = 'tbrl', nrow = 1)
ggsave(ga, filename = 'Plots/tree_liana_concept_final_rb_150.jpeg', height = 8, width = 14, units = 'in')

#########################
## Established: 400 m2 ##
#########################

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Interpolation_mxh_100.RData')

nsite = 100
nmonth = 12
nK = 75000

nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

VPDs = VPD_interp_100
SWPs = SWP_interp_100

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

tree.req = matrix(, nrow = nsite, ncol = nsite)

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

#save(tree.req, file = 'tree.req.100interp.monthly.RData')

load(file = 'tree.req.100interp.monthly.RData')
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

ks = c(trees_K$K / 10, lianas_K$K)
nK = 5000
Ks = seq(min(ks), max(ks), length.out = nK)
Ks = Ks[Ks > 75925]

liana.req = matrix(, nrow = nsite, ncol = nsite)

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

#save(liana.req, file = 'liana.req.100interp.monthly.RData')

load(file = 'liana.req.100interp.monthly.RData')
liana.req = liana.req * 0.001
liana_req_melt = melt(liana.req)
colnames(liana_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

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

ra1_fin = ra1 + theme(legend.position = 'bottom', plot.margin = unit(c(0, 0, 0, 0), 'cm'))
ra2_fin = ra2 + theme(legend.position = 'bottom', plot.margin = unit(c(-1, 0, 1, 0), 'cm'))

ga = plot_grid(ra1_fin, ra2_fin, align = 'hv', axis = 'tbrl', nrow = 1)
ggsave(ga, filename = 'Plots/tree_liana_concept_final_rb_400.jpeg', height = 8, width = 14, units = 'in')

###############################
## Invasion scenario: 150 m2 ##
###############################

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Interpolation_mxh_100.RData')

nsite = 100
nmonth = 12
nK = 75000

nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

VPDs = VPD_interp_100
SWPs = SWP_interp_100

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

tree.req = matrix(, nrow = nsite, ncol = nsite)

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

#save(tree.req, file = 'tree.req.100interp.monthly_invasion_150.RData')

load(file = 'tree.req.100interp.monthly_invasion_150.RData')
tree.req = tree.req * 0.001
tree_req_melt = melt(tree.req)
colnames(tree_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

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

ggsave(plot = ra1, file = 'Plots/tree.sens.ksurv.raster_final_invasion_150.jpeg')

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)
Ks = Ks[Ks > 13525]

liana.req = matrix(, nrow = nsite, ncol = nsite)

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

#save(liana.req, file = 'liana.req.100interp.monthly_invasion_150.RData')

load(file = 'liana.req.100interp.monthly_invasion_150.RData')
liana.req = liana.req * 0.001
liana_req_melt = melt(liana.req)
colnames(liana_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

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

ggsave(ra2, file = 'Plots/liana.sens.ksurv.raster_final_invasion_150.jpeg')

ga_fin = plot_grid(ra1, ra2, nrow = 1, rel_widths = 1)
ggsave(ga_fin, filename = 'Plots/tree_liana_relsens_concept_final_invasion_150.jpeg', height = 8, width = 14, units = 'in')

###############################
## Invasion scenario: 400 m2 ##
###############################

rm(list = ls())

# Input parameters
load('param.input.RData')

# Met drivers
load('Interpolation_mxh_100.RData')

nsite = 100
nmonth = 12
nK = 75000

nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

VPDs = VPD_interp_100
SWPs = SWP_interp_100

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

tree.req = matrix(, nrow = nsite, ncol = nsite)

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

#save(tree.req, file = 'tree.req.100interp.monthly_invasion.RData')

load(file = 'tree.req.100interp.monthly_invasion.RData')
tree.req = tree.req * 0.001
tree_req_melt = melt(tree.req)
colnames(tree_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

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

ggsave(plot = ra1, file = 'Plots/tree.sens.ksurv.raster_final_invasion_400.jpeg')

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)
Ks = Ks[Ks > 33727]

liana.req = matrix(, nrow = nsite, ncol = nsite)

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


#save(liana.req, file = 'liana.req.100interp.monthly_invasion.RData')

load(file = 'liana.req.100interp.monthly_invasion.RData')
liana.req = liana.req * 0.001
liana_req_melt = melt(liana.req)
colnames(liana_req_melt) = c('VPD_site', 'SWP_site', 'Ksurv')

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

ggsave(ra2, file = 'Plots/liana.sens.ksurv.raster_final_invasion_400.jpeg')

ga_fin = plot_grid(ra1, ra2, nrow = 1, rel_widths = 1)
ggsave(ga_fin, filename = 'Plots/tree_liana_relsens_concept_final_invasion_400.jpeg', height = 8, width = 14, units = 'in')

############################################
## Kreq-Kw comparison: 150 m2 established ##
############################################

rm(list = ls())

# Input parameters
load('param.input.RData')

# Rerun model with Ca as an input (instead of fixed in model)
tree.NPP.mxh = function(tree.dbh, psis, D_vec, tree.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.tree.al, tree.height, Vm = NULL, Ca){
  
  tree.tot.area = ((pi*tree.dbh^2) / 4) * 0.0001 # Total stem cross-sectional area (cm^2)
  
  tree.sap_frac = 0.622 # Fraction of total cross-sectional area that is sapwood (unitless)
  
  tree.ax =  min(tree.tot.area, (2.41 * (tree.dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectional area, Reyes-García et al. 2012 (m^2)
  
  tree.al = tot.al * frac.tree.al # Leaf area (m^2)
  
  Ca = Ca # Atmospheric CO2 concentration (ppm)
  
  #tree.stem.biom = (pi * (tree.dbh / 2)^2 * (100 * tree.height)) * tree.ssd
  
  #starX = (tree.stem.biom * 0.001) * tree.sap_frac
  starX = 0.0815 * tree.dbh^2.5 * tree.sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  #rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  # Extract one day's worth of VPD values and make moving average with one previous step
  D = movavg(D_vec, 2, type = 's')
  
  # Loop through 1 day
  # Note that because of how VPD got indexed 1 = midnight, 2 = 1 AM, etc
  for(hour in 1:24){
    
    rG = 0.3 # Tree growth respiration, default
    
    q = 1.89 # Fine root:leaf area ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # If specified in inputs
    }else{
      SLA = 32 # Specific leaf area (m^2 / kg C), default, supported by TRY analysis for both trees & lianas
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kg C), default
    
    Lr = 1.8*10^4 # Equivalent path length for roots, m, default
    
    Lx = tree.height # Stem length
    
    Lp = 0.1 # Petiole length, m, default
    
    Xstar = starX
    
    ar = tree.al / SLA * q * SRA # Fine root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg C)
    
    X = tree.ax * rho * Lx / 2 # Actional functional xylem biomass (kg), default equation
    
    Stem.biom = tree.tot.area * rho * Lx / 2 # Total stem biomass (kg), same as X but with total cross-sectional area, not xylem area
    
    Leaf.biom = (1 / (SLA/2)) * tree.al # Total leaf biomass (kg)
    
    ap = median(tree.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m2), default equation
    
    Kmax = tree.k # Maximum hydraulic conductivity (mmol/m/s/MPa), from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b1 = b1 # If specified as an input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Parameter for logistic function, from meta-analysis
    }else{
      b2 = b2 # If specified as an input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum hydraulic conductivity (mmol/m/s), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of hydraulic conductivity (mmol/m/s), default equation
    
    if(is.null(Vm)){
      Vm = 50 # maximum carboxylation rate (umol/m2/s), default value
    }else{
      Vm = Vm # If specified as an input
    }
    
    r = 0.5 # Leaf daytime respiration rate (umol/m2/s)
    
    gamma_star = 30 # CO2 compesnation point (ppm)
    
    km = 300 #Michaelis constant (ppm)
    
    D0 = 350 # Empirical coefficient (Pa)
    
    a1 = 30 # Empirical coefficient (Pa)
    
    gamma = 125 # CO2 compensation point (ppm)
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / tree.ax)) * (Lx * ar + tree.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default
    
    A = ((-Z * phi_max) / g) - (tree.al * V * D2) + (tree.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * tree.al * V * D2) - (tree.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Sum the same way
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_m - rp_m) * 1e-9 * 12
  
  stem.turn = 0.02 / 12 # yearly turnover / months in a year
  
  Lx_lost = Lx * stem.turn # fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = (1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    tree.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  # If NPP is still negative, repeat with higher stem turnover rate
  
  if(NPP < 0){
    stem.turn = 0.0324 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (tree.dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = ((1 - (rG/nday[month])) * (A1_m * tree.al - rd_m * (tree.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  }
  return(c(NPP, Lx, tree.al, stem.turn, Lx_turn, A1))
}

liana.NPP.mxh = function(dbh, psis, D_vec, liana.k, SLA = NULL, b1 = NULL, b2 = NULL, nday, month, tot.al, frac.liana.al, liana.length, Vm = NULL, Ca){
  
  tot.area = ((pi*dbh^2) / 4) * 0.0001 # Total cross-sectional stem area (cm^2)
  
  sap_frac = 0.622 # Fraction of stem that is sapwood, default
  
  liana.ax =  min(tot.area, (2.41 * (dbh/2)^1.97) * 0.0001) # Functional xylem cross-sectiional area, from Reyes-García et al. 2012
  
  liana.al = tot.al * frac.liana.al # Leaf area, from inputs (m^2)
  
  Ca = Ca # Atmospheric CO2 concentration (ppm)
  
  #liana.stem.biom = (pi * ((dbh / 100) / 2)^2 * (liana.length)) * rho
  
  #liana.starX = liana.stem.biom * sap_frac
  liana.starX = 0.0815 * dbh^2.5 * sap_frac # Optimal xylem dry C weight (kg), default equation
  
  # Storage for hourly values
  rG_h = c()
  rd_h = c()
  rr_h = c()
  rx_h = c()
  rp_h = c()
  A1_h = c()
  
  D = movavg(D_vec, 2, type = 's') # Smoothed, hourly VPD
  
  # Loop through 1 day
  # Note hour remark above
  for(hour in 1:24){
    rG = 0.3 # Growth respiration, default
    
    q = 1.89 # Fine root:leaf ratio, default
    
    if(!is.null(SLA)){
      SLA = SLA # SLA if it is included in inputs
    }else{
      SLA = 32 # Default SLA (m^2 / kgC)
    }
    
    rho = 420 # Wood density (kg / m^3)
    
    SRA = 80 # Specific root area (m^2 / kgC), default
    
    Lr = 1.8*10^4 # Equivalent root path length (m), default
    
    Lx = liana.length # Stem length (m)
    
    Lp = 0.1 # Petiole length (m), default
    
    Xstar = liana.starX # Optimal xylem dry C weight, same as above
    
    ar = liana.al / SLA * q * SRA # Root area (m^2), default equation
    
    Br = ar / SRA # Root biomass (kg)
    
    X = liana.ax * rho * Lx / 2 # Actual functional xylem biomass (kg), default equation
    
    Stem.biom = tot.area * rho * Lx / 2 # Stem biomass (kg), same as X but for whole cross-sectional area
    
    Leaf.biom = (1 / (SLA/2)) * liana.al # Leaf biomass (kg), LMA * AL
    
    ap = median(liana.al / SLA / (5 * rho / 3) / Lp) # Petiole area (m^2), default equation
    
    Kmax = liana.k # Hydraulic conductivity, from meta-analysis
    
    if(is.null(b1)){
      b1 = b1.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b1 = b1 # Or input
    }
    
    if(is.null(b2)){
      b2 = b2.mean.all # Logistic function parameter, from meta-analysis
    }else{
      b2 = b2 # Or input
    }
    
    phi_max = (Kmax * log(exp(-b1 * b2) + 1)) / b1 # Integral of maximum conductivity (mmol/m/s/MPa), default equation
    
    phi_s = (Kmax * log(exp(b1 * psis - b1 * b2) + 1)) / b1 # Integral of conductivity (mmol/m/s/MPa), default equation
    
    if(is.null(Vm)){
      Vm = 50 # Maximum carboxylation rate (umol/m2/s), default
    }else{
      Vm = Vm # Specified input rate
    }
    
    r = 0.5 # Daytime leaf respiration rate (umol/m2/s), default
    
    gamma_star = 30 # CO2 compensation point (ppm), default
    
    km = 300 # Michaelis constant (ppm), default
    
    D0 = 350 # Empirical coefficient (Pa), default
    
    a1 = 30 # Empirical coefficient (Pa), default
    
    gamma = 125 # CO2 compensation point (ppm), default
    
    # All equations as default
    
    T1 = 1 / (Ca - gamma_star)
    
    T2 = 1 / (Ca + km)
    
    g = a1 / ((Ca - gamma) * (1 + (D[hour] / D0)))
    
    V = Vm * ((Ca - gamma_star) / (Ca + km))
    
    Z = (ap / Lp) / (1 + ((ap / Lp) * (Lx / liana.ax)) * (Lx * ar + liana.ax * Lr) / (Lx * ar))
    
    D2 = D[hour] / 101 * 1.6
    
    # All equations as default 
    
    A = ((-Z * phi_max) / g) - (liana.al * V * D2) + (liana.al * r * D2)
    
    B = (Z * phi_s) + ((Z * T2 * phi_max) / g) + (T1 * liana.al * V * D2) - (liana.al * r * T2 * D2)
    
    C = -Z * T2 * phi_s
    
    xc = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    rd = 1
    
    rx = 3
    
    rx = rx * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    
    rp = 11.5
    
    rp = rp * (Lx * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    
    rr = 0.0949 * Br / (122 * 24 * 3600) * 10^9
    
    A1 = V * ((xc - T1) / (xc - T2))
    
    # Multiply each respiration component by 3600 seconds (in 1 hr)
    rd_h[hour] = rd * 3600
    rr_h[hour] = rr * 3600
    rx_h[hour] = rx * 3600
    rp_h[hour] = rp * 3600
    
    # If it's night,
    if(hour %in% c(1:6,19:24)){
      # GPP = 0
      A1_h[hour] = 0
    }
    
    # If it's day
    if(hour %in% c(7:18)){
      # multiply instantaneous GPP by 3600 seconds (in 1 hr)
      A1_h[hour] = A1 * 3600
    }
  }
  
  # Sum respiration for each hour in one day
  # Equivalent to umol/h * 1 h = umol
  rd_d = sum(rd_h)
  rr_d = sum(rr_h)
  rx_d = sum(rx_h)
  rp_d = sum(rp_h)
  
  # Sum GPP the same way
  A1_d = sum(A1_h)
  
  # Repeat for monthly scale
  rd_m = rd_d * nday[month]
  rr_m = rr_d * nday[month]
  rx_m = rx_d * nday[month]
  rp_m = rp_d * nday[month]
  A1_m = A1_d * nday[month]
  
  Anet = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_m - rp_m) * 1e-9 * 12)
  
  stem.turn = 0.1 / 12 #yearly turnover / months
  
  Lx_lost = Lx * stem.turn #fraction of stem lost
  Lx_turn = Lx - Lx_lost
  
  # Calculates new xylem respiration by subtracting the proportion lost
  rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
  rx_lost = rx_lost * 3600
  rx_lost = rx_lost * 24
  rx_lost = rx_lost * nday[month]
  rx_turn = rx_m - rx_lost
  
  # Calculates new phloem respiration by subtracting the proportion lost
  rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
  rp_lost = rp_lost * 3600
  rp_lost = rp_lost * 24
  rp_lost = rp_lost * nday[month]
  rp_turn = rp_m - rp_lost
  
  # NPP under non-water-stressed conditions
  NPP_prelim = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12)
  
  # Allows leaves to drop if water stress is present
  if(NPP_prelim < 0){
    liana.al = 0
    NPP = ((1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1 * 10^-9 * 12)
  }else{
    NPP = NPP_prelim
  }
  
  if(NPP < 0){
    stem.turn = 0.162 / 12
    
    Lx_lost = Lx * stem.turn
    Lx_turn = Lx - Lx_lost
    
    rx_lost = rx * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6 * (X / Xstar)
    rx_lost = rx_lost * 3600
    rx_lost = rx_lost * 24
    rx_lost = rx_lost * nday[month]
    rx_turn = rx_m - rx_lost
    
    rp_lost = rp * (Lx_lost * pi * (dbh / 100)) / (122 * 24 * 3600) * 10^6
    rp_lost = rp_lost * 3600
    rp_lost = rp_lost * 24
    rp_lost = rp_lost * nday[month]
    rp_turn = rp_m - rp_lost
    
    NPP = (1 - (rG/nday[month])) * (A1_m * liana.al - rd_m * (liana.al) - rr_m - rx_turn - rp_turn) * 1e-9 * 12
  }
  return(c(NPP, liana.al, stem.turn, Lx_turn, A1))
}

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

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI
tree.req_bci = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:nmonth){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 1.5 * 100
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD_future[mo,,vsite]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(Ca = 550, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
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
tree.req_h = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:nmonth){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 1.5 * 100
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD_future[mo,,vsite]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
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

save(tree.req_bci, tree.req_h, file = 'tree_req_future_established_150_wCO2.RData')

pl1 = ggplot() +
  geom_point(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Tree') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl1

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI
liana.req_bci = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 1.5 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD_future[mo,,vsite]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
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
liana.req_h = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 1.5 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD_future[mo,,vsite]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
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

save(liana.req_bci, liana.req_h, file = 'liana_req_future_established_150_wCO2.RData')

pl2 = ggplot() +
  geom_point(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Liana') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl2

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend = get_legend(pl2)

ga = grid.arrange(pl1 + theme(legend.position = 'none'), 
                  pl2 + theme(legend.position = 'none'), 
                  legend, nrow = 1, widths = c(2.3, 2.3, 0.4))

ggsave(ga, filename = 'Plots/both_established_fugure_150_wCO2.jpeg', width = 24, height = 10, units = 'in')

# Convert Ks to mol/m/s/MPa
tree_K = trees_K$K * 0.001
liana_K = lianas_K$K * 0.001

# Convert Ks to Kw by dividing by 10
tree_K = tree_K / 10
liana_K = liana_K / 10

# Make histograms
pl1 = ggplot() +
  geom_histogram(aes(x = tree_K), bins = 20) +
  geom_vline(aes(xintercept = min(tree.req_bci), color = 'Wettest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(tree.req_bci), color = 'Wettest', linetype = 'Future'), size = 1.2) +
  geom_vline(aes(xintercept = min(tree.req_h), color = 'Driest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(tree.req_h), color = 'Driest', linetype = 'Future'), size = 1.2) +
  xlab('K (mol/m/s/MPa)') + ylab('Frequency') +
  scale_color_discrete(name = 'Hydroclimate') +
  scale_linetype_manual(name = 'Time Period', values = c('Present' = 'solid', 'Future' = 'dashed'), breaks = c('Present', 'Future')) +
  ggtitle('Tree') +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl1

pl2 = ggplot() +
  geom_histogram(aes(x = liana_K), bins = 20) +
  geom_vline(aes(xintercept = min(liana.req_bci), color = 'Wettest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(liana.req_bci), color = 'Wettest', linetype = 'Future'), size = 1.2) +
  geom_vline(aes(xintercept = min(liana.req_h), color = 'Driest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(liana.req_h), color = 'Driest', linetype = 'Future'), size = 1.2) +
  xlab('K (mol/m/s/MPa)') + ylab('Frequency') +
  scale_color_discrete(name = 'Hydroclimate') +
  scale_linetype_manual(name = 'Time Period', values = c('Present' = 'solid', 'Future' = 'dashed'), breaks = c('Present', 'Future')) +
  ggtitle('Liana') +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, 'cm'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl2

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend = get_legend(pl1)

ga = plot_grid(pl1 + theme(legend.position = 'none'), 
               pl2 + theme(legend.position = 'none'),
               legend, nrow = 1, rel_widths = c(1, 1, 0.3))
ga

ggsave(ga, filename = 'Plots/Kreq_Kw_comp_150_wCO2.jpeg', width = 20, height = 12, units = 'in')

############################################
## Kreq-Kw comparison: 400 m2 established ##
############################################

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

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI
tree.req_bci = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:nmonth){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.6
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD_future[mo,,vsite]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
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
tree.req_h = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    # Specify Ca = 550 for future
    Ca = 550
    for(mo in 1:nmonth){
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.6
      
      psis = H_SWP[mo]
      
      D = H_VPD_future[mo,,vsite]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
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

save(tree.req_bci, tree.req_h, file = 'tree_req_future_established_400_wCO2.RData')

pl1 = ggplot() +
  geom_point(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Tree') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl1

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI
liana.req_bci = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD_future[mo,,vsite]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
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
liana.req_h = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = mean.data.liana.dbh
      
      tot.al = 4 * 100
      frac.liana.al = 0.4
      
      psis = H_SWP[mo]
      
      D = H_VPD_future[mo,,vsite]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
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

save(liana.req_bci, liana.req_h, file = 'liana_req_future_established_400_wCO2.RData')

pl2 = ggplot() +
  geom_point(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Liana') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl2

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend = get_legend(pl2)

ga = grid.arrange(pl1 + theme(legend.position = 'none'), 
                  pl2 + theme(legend.position = 'none'), 
                  legend, nrow = 1, widths = c(2.3, 2.3, 0.4))

ggsave(ga, filename = 'Plots/both_established_fugure_400_wCO2.jpeg', width = 24, height = 10, units = 'in')

# Convert Ks to mol/m/s/MPa
tree_K = trees_K$K * 0.001
liana_K = lianas_K$K * 0.001

# Convert Ks to Kw by dividing by 10
tree_K = tree_K / 10
liana_K = liana_K / 10

# Make histograms
pl1 = ggplot() +
  geom_histogram(aes(x = tree_K), bins = 20) +
  geom_vline(aes(xintercept = min(tree.req_bci), color = 'Wettest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(tree.req_bci), color = 'Wettest', linetype = 'Future'), size = 1.2) +
  geom_vline(aes(xintercept = min(tree.req_h), color = 'Driest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(tree.req_h), color = 'Driest', linetype = 'Future'), size = 1.2) +
  xlab('K (mol/m/s/MPa)') + ylab('Frequency') +
  scale_color_discrete(name = 'Hydroclimate') +
  scale_linetype_manual(name = 'Time Period', values = c('Present' = 'solid', 'Future' = 'dashed'), breaks = c('Present', 'Future')) +
  ggtitle('Tree') +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl1

pl2 = ggplot() +
  geom_histogram(aes(x = liana_K), bins = 20) +
  geom_vline(aes(xintercept = min(liana.req_bci), color = 'Wettest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(liana.req_bci), color = 'Wettest', linetype = 'Future'), size = 1.2) +
  geom_vline(aes(xintercept = min(liana.req_h), color = 'Driest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(liana.req_h), color = 'Driest', linetype = 'Future'), size = 1.2) +
  xlab('K (mol/m/s/MPa)') + ylab('Frequency') +
  scale_color_discrete(name = 'Hydroclimate') +
  scale_linetype_manual(name = 'Time Period', values = c('Present' = 'solid', 'Future' = 'dashed'), breaks = c('Present', 'Future')) +
  ggtitle('Liana') +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, 'cm'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl2

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend = get_legend(pl1)

ga = plot_grid(pl1 + theme(legend.position = 'none'), 
               pl2 + theme(legend.position = 'none'),
               legend, nrow = 1, rel_widths = c(1, 1, 0.3))
ga

ggsave(ga, filename = 'Plots/Kreq_Kw_comp_400_wCO2.jpeg', width = 20, height = 12, units = 'in')

#########################################
## Kreq-Kw comparison: 150 m2 invasion ##
#########################################

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI
tree.req_bci = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:nmonth){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 1.5 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD_future[mo,,vsite]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
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
tree.req_h = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:nmonth){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 1.5 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD_future[mo,,vsite]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
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

save(tree.req_bci, tree.req_h, file = 'tree_req_future_invasion_150_wCO2.RData')

pl1 = ggplot() +
  geom_point(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Tree') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl1

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI
liana.req_bci = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 1.5 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD_future[mo,,vsite]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
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
liana.req_h = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 1.5 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD_future[mo,,vsite]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
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

save(liana.req_bci, liana.req_h, file = 'liana_req_future_invasion_150_wCO2.RData')

pl2 = ggplot() +
  geom_point(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Liana') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl2

legend = get_legend(pl1)

ga = grid.arrange(pl1 + theme(legend.position = 'none'), 
                  pl2 + theme(legend.position = 'none'),
                  legend, nrow = 1, widths = c(2.3, 2.3, 0.4))

ggsave(ga, filename = 'Plots/both_invasion_fugure_150_wCO2.jpeg', width = 24, height = 10, units = 'in')

# Convert Ks to mol/m/s/MPa
tree_K = trees_K$K * 0.001
liana_K = lianas_K$K * 0.001

# Convert Ks to Kw by dividing by 10
tree_K = tree_K / 10
liana_K = liana_K / 10

# Make histograms
pl1 = ggplot() +
  geom_histogram(aes(x = tree_K), bins = 20) +
  geom_vline(aes(xintercept = min(tree.req_bci), color = 'Wettest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(tree.req_bci), color = 'Wettest', linetype = 'Future'), size = 1.2) +
  geom_vline(aes(xintercept = min(tree.req_h), color = 'Driest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(tree.req_h), color = 'Driest', linetype = 'Future'), size = 1.2) +
  xlab('K (mol/m/s/MPa)') + ylab('Frequency') +
  scale_color_discrete(name = 'Hydroclimate') +
  scale_linetype_manual(name = 'Time Period', values = c('Present' = 'solid', 'Future' = 'dashed'), breaks = c('Present', 'Future')) +
  ggtitle('Tree') +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl1

pl2 = ggplot() +
  geom_histogram(aes(x = liana_K), bins = 20) +
  geom_vline(aes(xintercept = min(liana.req_bci), color = 'Wettest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(liana.req_bci), color = 'Wettest', linetype = 'Future'), size = 1.2) +
  geom_vline(aes(xintercept = min(liana.req_h), color = 'Driest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(liana.req_h), color = 'Driest', linetype = 'Future'), size = 1.2) +
  xlab('K (mol/m/s/MPa)') + ylab('Frequency') +
  scale_color_discrete(name = 'Hydroclimate') +
  scale_linetype_manual(name = 'Time Period', values = c('Present' = 'solid', 'Future' = 'dashed'), breaks = c('Present', 'Future')) +
  ggtitle('Liana') +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, 'cm'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl2

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend = get_legend(pl1)

ga = plot_grid(pl1 + theme(legend.position = 'none'), 
               pl2 + theme(legend.position = 'none'),
               legend, nrow = 1, rel_widths = c(1, 1, 0.3))
ga

ggsave(ga, filename = 'Plots/Kreq_Kw_comp_invasion_150_wCO2.jpeg', width = 20, height = 12, units = 'in')

#########################################
## Kreq-Kw comparison: 400 m2 invasion ##
#########################################

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI
tree.req_bci = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:nmonth){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD_future[mo,,vsite]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
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
tree.req_h = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:nmonth){
      # Specify Ca = 550  for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      tree.dbh = mean.hori.data.tree.dbh
      
      tot.al = 4 * 100
      frac.tree.al = 0.9
      
      psis = H_SWP[mo]
      
      D = H_VPD_future[mo,,vsite]
      
      if(mo == 1){
        tree.height = tree.b1Ht * tree.dbh^tree.b2Ht
      }else{
        tree.height = tree.out[5]
      }
      K = Ks[k]
      tree.out = tree.NPP.mxh(Ca = Ca, tree.dbh = tree.dbh, psis = psis, D_vec = D, tree.k = K, nday = nday, month = month, tot.al = tot.al, frac.tree.al = frac.tree.al, tree.height = tree.height, Vm = NULL)
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

save(tree.req_bci, tree.req_h, file = 'tree_req_future_invasion_400_wCO2.RData')

pl1 = ggplot() +
  geom_point(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = tree.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = tree.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Tree') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl1

ks = c(trees_K$K / 10, lianas_K$K)
Ks = seq(min(ks), max(ks), length.out = nK)

## BCI
liana.req_bci = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = BCI_SWP[mo]
      
      D = BCI_VPD_future[mo,,vsite]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
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
liana.req_h = c()

count = 0

for(vsite in 1:vpdsite){
  sens = c()
  for(k in 1:nK){
    for(mo in 1:12){
      # Specify Ca = 550 for future
      Ca = 550
      nday = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      month = mo
      
      dbh = 2
      
      tot.al = 4 * 100
      frac.liana.al = 0.1
      
      psis = H_SWP[mo]
      
      D = H_VPD_future[mo,,vsite]
      
      if(mo == 1){
        liana.length = tree.out[2]
      }else{
        liana.length = lianaout[4]
      }
      K = Ks[k]
      lianaout = liana.NPP.mxh(Ca = Ca, dbh = dbh, psis = psis, D_vec = D, liana.k = K, SLA = NULL, b1 = NULL, b2 = NULL, nday = nday, month = month, tot.al = tot.al, frac.liana.al = frac.liana.al, liana.length = liana.length, Vm = NULL)
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

save(liana.req_bci, liana.req_h, file = 'liana_req_future_invasion_400_wCO2.RData')

pl2 = ggplot() +
  geom_point(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.3) +
  geom_point(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.3) +
  geom_line(aes(x = mult, y = liana.req_bci, color = 'Wettest'), size = 1.2) +
  geom_line(aes(x = mult, y = liana.req_h, color = 'Driest'), size = 1.2) +
  xlab('Increase from present') + ylab(expression(paste(K[req],' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  scale_color_npg(name = 'Site') +
  ggtitle('Liana') +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl2

legend = get_legend(pl1)

ga = grid.arrange(pl1 + theme(legend.position = 'none'), 
                  pl2 + theme(legend.position = 'none'),
                  legend, nrow = 1, widths = c(2.3, 2.3, 0.4))

ggsave(ga, filename = 'Plots/both_invasion_fugure_400_wCO2.jpeg', width = 24, height = 10, units = 'in')

# Convert Ks to mol/m/s/MPa
tree_K = trees_K$K * 0.001
liana_K = lianas_K$K * 0.001

# Convert Ks to Kw by dividing by 10
tree_K = tree_K / 10
liana_K = liana_K / 10

# Make histograms
pl1 = ggplot() +
  geom_histogram(aes(x = tree_K), bins = 20) +
  geom_vline(aes(xintercept = min(tree.req_bci), color = 'Wettest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(tree.req_bci), color = 'Wettest', linetype = 'Future'), size = 1.2) +
  geom_vline(aes(xintercept = min(tree.req_h), color = 'Driest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(tree.req_h), color = 'Driest', linetype = 'Future'), size = 1.2) +
  xlab('K (mol/m/s/MPa)') + ylab('Frequency') +
  scale_color_discrete(name = 'Hydroclimate') +
  scale_linetype_manual(name = 'Time Period', values = c('Present' = 'solid', 'Future' = 'dashed'), breaks = c('Present', 'Future')) +
  ggtitle('Tree') +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl1

pl2 = ggplot() +
  geom_histogram(aes(x = liana_K), bins = 20) +
  geom_vline(aes(xintercept = min(liana.req_bci), color = 'Wettest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(liana.req_bci), color = 'Wettest', linetype = 'Future'), size = 1.2) +
  geom_vline(aes(xintercept = min(liana.req_h), color = 'Driest', linetype = 'Present'), size = 1.2) +
  geom_vline(aes(xintercept = max(liana.req_h), color = 'Driest', linetype = 'Future'), size = 1.2) +
  xlab('K (mol/m/s/MPa)') + ylab('Frequency') +
  scale_color_discrete(name = 'Hydroclimate') +
  scale_linetype_manual(name = 'Time Period', values = c('Present' = 'solid', 'Future' = 'dashed'), breaks = c('Present', 'Future')) +
  ggtitle('Liana') +
  theme_linedraw() +
  theme(legend.key.height = unit(1.5, 'cm'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))
pl2

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend = get_legend(pl1)

ga = plot_grid(pl1 + theme(legend.position = 'none'), 
               pl2 + theme(legend.position = 'none'),
               legend, nrow = 1, rel_widths = c(1, 1, 0.3))
ga

ggsave(ga, filename = 'Plots/Kreq_Kw_comp_invasion_400_wCO2.jpeg', width = 20, height = 12, units = 'in')
