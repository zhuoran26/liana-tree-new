## This script explores relationships between site dryness variables
## and measured Kmax in the literature

## Dry season length is taken for the papers associated with Kmax observations
## If the site does not have a dry season, dry season length = 0

## Kmax values are those from the meta analysis

## The spreadsheet was compiled by matching site to species-level observation
## for each study and compiling with K and P50 data from those studies
## If the same species was measured in multiple papers, study-level data
## were used
## If more than one site was used in the same paper, species-level observations
## were matched to the site to the author's best ability

## Author: AM Willson
## Date modified: 24 September 2021

rm(list = ls())
library(ggplot2)
setwd('~/Google Drive 2/Medvigy_SP20/liana_model/liana_model/')

# Read in compiled data frame
data = read.csv('hydraulic.trait_clim_geo.csv')

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
data$K_molmsMPa = data$K * 0.001

##################
## K - Latitude ##
##################

ggplot(data, aes(x = Latitude, y = K_molmsMPa)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Latitude') + ylab(expression(paste(K[s], ' (mol/m/s/MPa)'))) +
  annotate('text', x = -32, y = 6000, label = 'y = 413 + 6.52x', size = 5) +
  annotate('text', x = -35.5, y = 5600, label = expression(paste(R[adj]^2, ' = 0.01')), size = 5) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26))

summary(lm(data$K_molmsMPa ~ data$Latitude))

ggplot(data, aes(x = Latitude, y = K_molmsMPa, group = Growth.form, color = Growth.form)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Latitude') + ylab(expression(paste(K[s], ' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  scale_color_discrete(name = 'Growth form') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))

summary(lm(data$K[which(data$Growth.form == 'tree')] ~ data$Latitude[which(data$Growth.form == 'tree')]))
summary(lm(data$K[which(data$Growth.form == 'liana')] ~ data$Latitude[which(data$Growth.form == 'liana')]))

###################
## K - Longitude ##
###################

ggplot(data, aes(x = Longitude, y = K_molmsMPa)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Longitude') + ylab(expression(paste(K[s], ' (mol/m/s/MPa)'))) +
  annotate('text', x = -50, y = 6000, label = 'y = 496 - 1.54x', size = 5) +
  annotate('text', x = -62, y = 5600, label = expression(paste(R[adj]^2, ' = 0.03')), size = 5) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26))

summary(lm(data$K_molmsMPa ~ data$Longitude))

ggplot(data, aes(x = Longitude, y = K_molmsMPa, group = Growth.form, color = Growth.form)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Longitude') + ylab(expression(paste(K[s], ' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  scale_color_discrete(name = 'Growth form') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))

summary(lm(data$K[which(data$Growth.form == 'tree')] ~ data$Longitude[which(data$Growth.form == 'tree')]))
summary(lm(data$K[which(data$Growth.form == 'liana')] ~ data$Longitude[which(data$Growth.form == 'liana')]))

##################
## K - Altitude ##
##################

ggplot(data, aes(x = Altitude, y = K_molmsMPa)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Altitude (m)') + ylab(expression(paste(K[s], ' (mol/m/s/MPa)'))) +
  annotate('text', x = 310, y = 6100, label = 'y = 641 - 0.39x', size = 5) +
  annotate('text', x = 210, y = 5700, label = expression(paste(R[adj]^2, ' = 0.03')), size = 5) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26))

summary(lm(data$K_molmsMPa ~ data$Altitude))

ggplot(data, aes(x = Altitude, y = K_molmsMPa, group = Growth.form, color = Growth.form)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Altitude (m)') + ylab(expression(paste(K[s], ' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  scale_color_discrete(name = 'Growth form') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))

summary(lm(data$K[which(data$Growth.form == 'tree')] ~ data$Altitude[which(data$Growth.form == 'tree')]))
summary(lm(data$K[which(data$Growth.form == 'liana')] ~ data$Altitude[which(data$Growth.form == 'liana')]))

###########################
## K - dry season length ##
###########################

ggplot(data, aes(x = Dry_len, y = K_molmsMPa)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Dry season\nlength (months)') + ylab(expression(paste(K[s], ' (mol/m/s/MPa)'))) +
  annotate('text', x = 1, y = 6000, label = 'y = 312 + 39.2x', size = 5) +
  annotate('text', x = 0.62, y = 5600, label = expression(paste(R[adj]^2, ' < 0.01')), size = 5) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26))

summary(lm(data$K_molmsMPa ~ data$Dry_len))
summary(aov(data$K ~ data$Dry_len))

ggplot(data, aes(x = Dry_len, y = K, group = Growth.form, color = Growth.form)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Dry season\nlength (months)') + ylab(expression(paste(K[s], ' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  scale_color_discrete(name = 'Growth form') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))

summary(lm(data$K[which(data$Growth.form == 'tree')] ~ data$Dry_len[which(data$Growth.form == 'tree')]))
summary(lm(data$K[which(data$Growth.form == 'liana')] ~ data$Dry_len[which(data$Growth.form == 'liana')]))

summary(aov(data$K[which(data$Growth.form == 'tree')] ~ data$Dry_len[which(data$Growth.form == 'tree')]))
summary(aov(data$K[which(data$Growth.form == 'liana')] ~ data$Dry_len[which(data$Growth.form == 'liana')]))

######################
## K - field season ##
######################

data %>%
  filter(Meas_szn != '') %>%
  ggplot(aes(x = Meas_szn, y = K_molmsMPa, group = Meas_szn)) + 
  geom_boxplot() +
  scale_x_discrete(labels = c('Dry', 'Wet', 'Both')) +
  xlab('Season of measurement') + ylab(expression(paste(K[s], ' (mol/m/s/MPa)'))) +
  annotate('text', x = 'Dry', hjust = -1.25, y = 3500, label = 'ANOVA p < 0.01', size = 5) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26))

summary(aov(data$K ~ data$Meas_szn))
t.test(data$K[which(data$Meas_szn == 'Dry')], data$K[which(data$Meas_szn == 'Wet')])

data %>%
  filter(Meas_szn != '') %>%
  ggplot(aes(x = Meas_szn, y = K_molmsMPa, group = Meas_szn)) + 
  geom_boxplot() + facet_wrap(~Growth.form, labeller = labeller(Growth.form = c('liana' = 'Liana', 'tree' = 'Tree'))) +
  scale_x_discrete(labels = c('Dry', 'Wet', 'Both')) +
  xlab('Season of measurement') + ylab(expression(paste(K[s], ' (mol/m/s/MPa)'))) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26, angle = 45), strip.text = element_text(size = 30))

summary(aov(data$K[which(data$Growth.form == 'tree')] ~ data$Meas_szn[which(data$Growth.form == 'tree')]))
summary(aov(data$K[which(data$Growth.form == 'liana')] ~ data$Meas_szn[which(data$Growth.form == 'liana')]))

####################
## P50 - Latitude ##
####################

ggplot(data, aes(x = Latitude, y = P50)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Latitude') + ylab(expression(paste(P[50], ' (MPa)'))) +
  coord_cartesian(xlim = c(0, 30)) +
  annotate('text', x = 4, y = 0, label = 'y = -1.51 - 0.02x', size = 5) +
  annotate('text', x = 2.3, y = -0.3, label = expression(paste(R[adj]^2, ' = 0.01')), size = 5) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26))

summary(lm(data$P50 ~ data$Latitude))

ggplot(data, aes(x = Latitude, y = P50, group = Growth.form, color = Growth.form)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Latitude') + ylab(expression(paste(P[50], ' (MPa)'))) +
  coord_cartesian(xlim = c(0, 30)) +
  theme_linedraw() +
  scale_color_discrete(name = 'Growth form') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))

summary(lm(data$P50[which(data$Growth.form == 'tree')] ~ data$Latitude[which(data$Growth.form == 'tree')]))
summary(lm(data$P50[which(data$Growth.form == 'liana')] ~ data$Latitude[which(data$Growth.form == 'liana')]))

#####################
## P50 - Longitude ##
#####################

ggplot(data, aes(x = Longitude, y = P50)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Longitude') + ylab(expression(paste(P[50], ' (MPa)'))) +
  annotate('text', x = 130, y = 0, label = 'y = -1.78 - 0.001x', size = 5) +
  annotate('text', x = 147, y = -0.3, label = expression(paste(R[adj]^2, ' = 0.01')), size = 5) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26))

summary(lm(data$P50 ~ data$Longitude))

ggplot(data, aes(x = Longitude, y = P50, group = Growth.form, color = Growth.form)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Longitude') + ylab(expression(paste(P[50], ' (MPa)'))) +
  theme_linedraw() +
  scale_color_discrete(name = 'Growth form') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))

summary(lm(data$P50[which(data$Growth.form == 'tree')] ~ data$Longitude[which(data$Growth.form == 'tree')]))
summary(lm(data$P50[which(data$Growth.form == 'liana')] ~ data$Longitude[which(data$Growth.form == 'liana')]))

####################
## P50 - Altitude ##
####################

ggplot(data, aes(x = Altitude, y = P50)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Altitude (m)') + ylab(expression(paste(P[50], ' (MPa)'))) +
  coord_cartesian(xlim = c(0, 1000)) +
  annotate('text', x = 825, y = 0, label = 'y = -1.61 - 0.0004x', size = 5) +
  annotate('text', x = 910, y = -0.3, label = expression(paste(R[adj]^2, ' = 0.01')), size = 5) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26))

summary(lm(data$P50 ~ data$Altitude))

ggplot(data, aes(x = Altitude, y = P50, group = Growth.form, color = Growth.form)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Altitude (m)') + ylab(expression(paste(P[50], ' (MPa)'))) +
  coord_cartesian(xlim = c(0, 1000)) +
  theme_linedraw() +
  scale_color_discrete(name = 'Growth form') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))

summary(lm(data$P50[which(data$Growth.form == 'tree')] ~ data$Altitude[which(data$Growth.form == 'tree')]))
summary(lm(data$P50[which(data$Growth.form == 'liana')] ~ data$Altitude[which(data$Growth.form == 'liana')]))

#############################
## P50 - dry season length ##
#############################

ggplot(data, aes(x = Dry_len, y = P50)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Dry season\nlength (months)') + ylab(expression(paste(P[50], ' (MPa)'))) +
  coord_cartesian(xlim = c(3, 6)) +
  annotate('text', x = 3.42, y = 0, label = 'y = -1.13 - 0.13x', size = 5) +
  annotate('text', x = 3.25, y = -0.3, label = expression(paste(R[adj]^2, ' = 0.01')), size = 5) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26))

summary(lm(data$P50 ~ data$Dry_len))

ggplot(data, aes(x = Dry_len, y = P50, group = Growth.form, color = Growth.form)) + 
  geom_point() + geom_smooth(method = 'lm', se = F) +
  xlab('Dry season\nlength (months)') + ylab(expression(paste(P[50], ' (MPa)'))) +
  coord_cartesian(xlim = c(3, 6)) +
  theme_linedraw() +
  scale_color_discrete(name = 'Growth form') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), legend.title = element_text(size = 28), legend.text = element_text(size = 26))

summary(lm(data$P50[which(data$Growth.form == 'tree')] ~ data$Dry_len[which(data$Growth.form == 'tree')]))
summary(lm(data$P50[which(data$Growth.form == 'liana')] ~ data$Dry_len[which(data$Growth.form == 'liana')]))

summary(aov(data$P50[which(data$Growth.form == 'tree')] ~ data$Dry_len[which(data$Growth.form == 'tree')]))
summary(aov(data$P50[which(data$Growth.form == 'liana')] ~ data$Dry_len[which(data$Growth.form == 'liana')]))

########################
## P50 - field season ##
########################

data %>%
  filter(Meas_szn != '') %>%
  filter(Meas_szn != 'WetDry') %>%
  ggplot(aes(x = Meas_szn, y = P50, group = Meas_szn)) + 
  geom_boxplot() +
  annotate('text', x = 'Dry', hjust = -1.2, y = 0, label = 't = 2.38, p = 0.02', size = 5) +
  xlab('Season of measurement') + ylab(expression(paste(P[50], ' (MPa)'))) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26))

summary(aov(data$P50 ~ data$Meas_szn))
t.test(data$P50[which(data$Meas_szn == 'Dry')], data$P50[which(data$Meas_szn == 'Wet')])

data %>%
  filter(Meas_szn != '') %>%
  filter(Meas_szn != 'WetDry') %>%
  ggplot(aes(x = Meas_szn, y = P50, group = Meas_szn)) + 
  geom_boxplot() + facet_wrap(~Growth.form, labeller = labeller(Growth.form = c('liana' = 'Liana', 'tree' = 'Tree'))) +
  xlab('Season of measurement') + ylab(expression(paste(P[50], ' (MPa)'))) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 26), strip.text = element_text(size = 30))

summary(aov(data$P50[which(data$Growth.form == 'tree')] ~ data$Meas_szn[which(data$Growth.form == 'tree')]))
summary(aov(data$P50[which(data$Growth.form == 'liana')] ~ data$Meas_szn[which(data$Growth.form == 'liana')]))
t.test(data$P50[which(data$Meas_szn == 'Dry' & data$Growth.form == 'liana')], data$P50[which(data$Meas_szn == 'Wet' & data$Growth.form == 'liana')])
t.test(data$P50[which(data$Meas_szn == 'Dry' & data$Growth.form == 'tree')], data$P50[which(data$Meas_szn == 'Wet' & data$Growth.form == 'tree')])
