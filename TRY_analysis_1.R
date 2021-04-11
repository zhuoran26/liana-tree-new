## This script processes the TRY data for each trait
## Each text file contains data for one trait and all species with observations
## The output is a matrix with species-averaged trait values and growth form attribute
## This file must be used in TRY_analysis_2.R before being used for analysis

## Author: AM Willson
## Date modified: 28 January 2021

##################################################
## Plant growth form and species numbers (7754) ##
##################################################

library(data.table)

#Read in dataframe from text file
sp.ll=fread('7754.txt')
#Subset only necessary columns
sp.ll=subset(sp.ll,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.ll=subset(sp.ll,TraitID!=0)

#See the species and growth forms
spp=subset(sp.ll,sp.ll$TraitID==42)
spp = subset(spp, !is.na(spp$OrigValueStr))
#Sort the growth forms from most to least frequently appearing
speciestable=sort(table(spp$OrigValueStr),decreasing=T)
#View the list of growth forms
View(speciestable)

#Make a list of all the unique species ID numbers
uniquespecies=unique(spp$AccSpeciesID)
#Make the list into an integer for use in the for loop
uniquespecies=unlist(uniquespecies)
#Make a matrix for the for loop that has the number of rows corresponding to the maximum species ID number
output=matrix(,nrow=max(spp$AccSpeciesID),ncol=2)

count = 0

for(i in uniquespecies){
  sp = subset(spp, AccSpeciesID == i)
  pfts = sp$OrigValueStr
  if(length(grep('tree', pfts, ignore.case = T)) > 0 | length(grep('liana', pfts, ignore.case = T)) > 0){
    if(length(grep('tree', pfts, ignore.case = T)) > length(grep('liana', pfts, ignore.case = T))){
      type = 'tree'
    }
    if(length(grep('tree', pfts, ignore.case = T)) < length(grep('liana', pfts, ignore.case = T))){
      type = 'liana'
    }
  }else{
    if('T' %in% pfts | 'TS' %in% pfts | 'L' %in% pfts){
      if(length(which(pfts %in% c('T', 'TS'))) > length(which(pfts == 'L'))){
        type = 'tree'
      }
      if(length(which(pfts %in% c('T', 'TS'))) < length(which(pfts == 'L'))){
        type = 'liana'
      }
    }else{
      type = 'non-woody'
    }
  }
  output[i,] = c(i, type)
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

colnames=c('Species ID','Growth form')
colnames(output)=colnames

write.csv(output, 'growth_form_8-10-20.csv')

#############################################
## Leaf lifespan and species number (7754) ##
#############################################

#View the units of the leaf lifespan characteristic
unique(sp.ll$OrigUnitStr)
#Subset for only leaf lifespan observations
units=subset(sp.ll,sp.ll$TraitID==12)

#Check the measurements of leaf lifespan
unique(units$OriglName)
#Make a table of the frequency of the measurements of leaf lifespan
unit=sort(table(units$OriglName),decreasing=T)
#View frequency table
View(unit)
#Remove inconsistent measurements
units=subset(units,units$OriglName!='annual turnover of foliage (fraction of total leaf biomass lost per annum)')

#Check the units of this subset
unique(units$OrigUnitStr)
#Make a table of the frequency of the units of leaf lifespan
unit=sort(table(units$OrigUnitStr),decreasing=T)
#View the frequency table
View(unit)
#Take out nonsense units
units=subset(units,units$OrigUnitStr!='')

#Make the units consistent
for(i in 1:length(units$OrigUnitStr)){
  #Index the unit for each observation
  iunit=units$OrigUnitStr[i]
  #If the units are in months, then make the unit say "months"
  if(iunit%in%c('months','month','mo')){
    units$OrigUnitStr[i]='months'
    #If the unit is in years, then make the unit say "years"
  }else{
    if(iunit%in%c('years','yrs','yr','yr-1','year')){
      units$OrigUnitStr[i]='years'
      #If the unit is in days, then make the unit say "days"
    }else{
      if(iunit%in%c('days')){
        units$OrigUnitStr[i]='days'
        #Catch-all to make sure nothing is missed
      }else{
        units$OrigUnitStr[i]='NA'
      }
    }
  }
}

#Check to make sure there are no NA
unique(units$OrigUnitStr)

#Now, there are some values of leaf lifespan that are not numeric
for(i in 1:length(units$DatasetID)){
  if(units$OrigValueStr[i]%in%c('<=1yr','>=1yr')){
    units$OrigValueStr[i]=1
  }else{
    if(units$OrigValueStr[i]=='>24'){
      units$OrigValueStr[i]=24
    }else{
      if(units$OrigValueStr[i]=='0-2'){
        units$OrigValueStr[i]=1
      }else{
        if(units$OrigValueStr[i]=='12-24'){
          units$OrigValueStr[i]=18
        }else{
          if(units$OrigValueStr[i]=='2-6'){
            units$OrigValueStr[i]=4
          }else{
            if(units$OrigValueStr[i]=='3-6'){
              units$OrigValueStr[i]=4.5
            }else{
              if(units$OrigValueStr[i]=='6-12'){
                units$OrigValueStr[i]=9
              }else{
                if(units$OrigValueStr[i]=='7-12'){
                  units$OrigValueStr[i]=9.5
                }else{
                  units$OrigValueStr[i]=units$OrigValueStr[i]
                }
              }
            }
          }
        }
      }
    }
  }
}

#Make LL numeric
units$OrigValueStr=as.numeric(as.character(units$OrigValueStr))

#Convert all leaf lifespan values to months
for(i in 1:length(units$DatasetID)){
  #Index the units for each observation
  iunit=units$OrigUnitStr[i]
  #Index the leaf lifespan value for each observation
  ill=units$OrigValueStr[i]
  #If the unit is in months, keep the same value
  if(iunit=='months'){
    units$OrigValueStr[i]=ill
    #If the unit is in years, multiply the original value by 12
  }else{
    if(iunit=='years'){
      units$OrigValueStr[i]=ill*12
      #If the unit is in days, divide the original value by 30.42
    }else{
      if(iunit=='days'){
        units$OrigValueStr[i]=ill/30.42
        #Catch-all to make sure nothing is missed
      }else{
        units$OrigValueStr[i]=NA
      }
    }
  }
}

#Check for NAs
subset(units,units$OrigValueStr=='NA')

#Add a column to the output
output=cbind(output,NA)

count = 0
#Add leaf lifespan measurements to output
for(i in uniquespecies){
  #The species subset by species ID number
  sp=subset(units,units$AccSpeciesID==i)
  #Mean leaf lifespan for each species
  ll=mean(sp$OrigValueStr, na.rm = T)
  #Put in output created for plant growth form and species number
  output[i,3]=ll
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Check to make sure for one species that the leaf lifespan was computed correctly
llcheck=subset(sp.ll,sp.ll$AccSpeciesID==294)

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`Leaf lifespan (m)` = as.numeric(as.character(df$`Leaf lifespan (m)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`Leaf lifespan (m)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`Leaf lifespan (m)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Leaf lifespan (m)`))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Leaf lifespan (m)`)) +
  coord_cartesian(ylim = c(0, 100))

###################################
## SLA and species number (7766) ##
###################################

#Read in the dataframe
sp.sla=fread('7766.txt')
#Subset only necessary columns
sp.sla=subset(sp.sla,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.sla=subset(sp.sla,TraitID!=0)

#Look at the measurements of SLA
unique(sp.sla$OriglName)
#Make a table of the frequency of the measurements of SLA
unit=sort(table(sp.sla$OriglName),decreasing=T)
#View frequency table
View(unit)
#Remove inconsistent measurements
sp.sla=subset(sp.sla,sp.sla$OriglName!='SLA_max')
sp.sla=subset(sp.sla,sp.sla$OriglName!='SLA_min')
sp.sla=subset(sp.sla,sp.sla$OriglName!='SLA Max (mm\xb2/mg)')
sp.sla=subset(sp.sla,sp.sla$OriglName!='SLA Min (mm\xb2/mg)')

#Look at the units
unique(sp.sla$OrigUnitStr)
#Make a table of the frequency of the units of SLA
unit=sort(table(sp.sla$OrigUnitStr),decreasing=T)
#View the frequency table
View(unit)

#Make the units consistent
for(i in 1:length(sp.sla$OrigUnitStr)){
  if(sp.sla$OrigUnitStr[i]%in%c('mm2/mg','mm2 / mg','mm2 mg-1','mm2*mg-1','mm^2/mg','mm\xb2/mg')){
    sp.sla$OrigUnitStr[i]='mm2/mg'
  }else{
    if(sp.sla$OrigUnitStr[i]%in%c('cm2/g','cm2/g (n.r.)','cm2 g-1','cm2 / g','cm-2 g')){
      sp.sla$OrigUnitStr[i]='cm2/g'
    }else{
      if(sp.sla$OrigUnitStr[i]%in%c('m2/kg','m2 kg-1','m2/kg leaf DM','m\xb2 Kg-1')){
        sp.sla$OrigUnitStr[i]='m2/kg'
      }else{
        if(sp.sla$OrigUnitStr[i]%in%c('g/cm2')){
          sp.sla$OrigUnitStr[i]='g/cm2'
        }else{
          if(sp.sla$OrigUnitStr[i]%in%c('gDM m-2','g/m2')){
            sp.sla$OrigUnitStr[i]='g/m2'
          }else{
            if(sp.sla$OrigUnitStr[i]%in%c('m2/g')){
              sp.sla$OrigUnitStr[i]='m2/g'
            }else{
              if(sp.sla$OrigUnitStr[i]%in%c('mm2/g')){
                sp.sla$OrigUnitStr[i]='mm2/g'
              }
            }
          }
        }
      }
    }
  }
}

#Check the units after running the for loop
unique(sp.sla$OrigUnitStr)

#There are three null values--take out
sp.sla=subset(sp.sla,sp.sla$OrigValueStr!='NULL')

#Make SLA values numeric
sp.sla$OrigValueStr=as.numeric(sp.sla$OrigValueStr)

#Convert all SLA observations to the mm2/mg
for(i in 1:length(sp.sla$DatasetID)){
  #Indexes the unit column for each observation
  iunit=sp.sla$OrigUnitStr[i]
  #Indexes the value column for each observation
  isla=sp.sla$OrigValueStr[i]
  #If the units are in mm2/mg, then don't change the observation
  if(iunit=='mm2/mg'){
    sp.sla$OrigValueStr[i]=isla
    #If the units are in cm2/g, then multiply the value by 0.1
  }else{
    if(iunit=='cm2/g'){
      sp.sla$OrigValueStr[i]=isla*100/1000
      #If the units are in m2/kg, then don't change the observation
    }else{
      if(iunit=='m2/kg'){
        sp.sla$OrigValueStr[i]=isla*1000000/1000000
        #If the units are in g/cm2, then multiply the value by 10 and take the reciprocal
      }else{
        if(iunit=='g/cm2'){
          sp.sla$OrigValueStr[i]=(1/isla)*100/1000
          #If the units are in g/m2, then multiply the value by 0.001 and take the reciprocal
        }else{
          if(iunit=='g/m2'){
            sp.sla$OrigValueStr[i]=(1/isla)*1000000/1000
            #If the units are in m2/g, then multiply the value by 1000
          }else{
            if(iunit=='m2/g'){
              sp.sla$OrigValueStr[i]=isla*1000000/1000
              #If the units are in mm2/g, then multiply the value by 0.001
            }else{
              if(iunit=='mm2/g'){
                sp.sla$OrigValueStr[i]=isla/1000
                #If none of the above are true, then insert NA
              }else{
                sp.sla$OrigValueStr[i]=NA
              }
            }
          }
        }
      }
    }
  }
}

#Check to make sure there are no NA
subset(sp.sla,sp.sla$OrigValueStr=='NA')

#Use the same matrix as for the other observations
output=cbind(output,NA)

#Find the mean SLA for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.sla,AccSpeciesID==i)
  #Subset by only trait 3116 (= SLA)
  sp=subset(sp,TraitID==3116)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  sla=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,4]=sla
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)')
colnames(output)=colnames

df = as.data.frame(output)
df$`SLA (mm2/mg)` = as.numeric(as.character(df$`SLA (mm2/mg)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`SLA (mm2/mg)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`SLA (mm2/mg)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`SLA (mm2/mg)`))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`SLA (mm2/mg)`)) +
  coord_cartesian(ylim = c(0, 150))

#####################################
## Narea and species number (7767) ##
#####################################

#Read in the dataframe
sp.n=fread('7767.txt')
#Subset only necessary columns
sp.n=subset(sp.n,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.n=subset(sp.n,TraitID!=0)

#Look at the measurements for Narea
unique(sp.n$OriglName)
#All measurements appear consistent here

#Look at the units
unique(sp.n$OrigUnitStr)
#Make a table of the frequency of the units of Narea
unit=sort(table(sp.n$OrigUnitStr),decreasing=T)
#View the frequency table
View(unit)

#There are a bunch of observations that have no values (blank rows)--remove
sp.n=subset(sp.n,sp.n$OrigUnitStr!='')

#Make the units consistent
for(i in 1:length(sp.n$DatasetID)){
  #If units are mmol/m2, then make mmol/m2
  if(sp.n$OrigUnitStr[i]%in%c('mmol/m2')){
    sp.n$OrigUnitStr[i]='mmol/m2'
    #If units are g/m2, then make g/m2
  }else{
    if(sp.n$OrigUnitStr[i]%in%c('g / m2','g N m-2','g/m2','g m-2','g m2','gN m-2','gN.m-2')){
      sp.n$OrigUnitStr[i]='g/m2'
      #If units are mg/cm2, then make mg/cm2
    }else{
      if(sp.n$OrigUnitStr[i]%in%c('mg/cm2','mg cm-2')){
        sp.n$OrigUnitStr[i]='mg/cm2'
        #If none of the above, then make NA
      }else{
        sp.n$OrigUnitStr[i]=NA
      }
    }
  }
}

#Check the units after running the for loop
unique(sp.n$OrigUnitStr)

#Make Narea values numeric
sp.n$OrigValueStr=as.numeric(sp.n$OrigValueStr)

#Convert all Narea observations to g/m2
for(i in 1:length(sp.n$DatasetID)){
  #Indexes the unit for each observation
  iunit=sp.n$OrigUnitStr[i]
  #Indexes the value for each observation
  iN=sp.n$OrigValueStr[i]
  #If the units are mmol/m2, then multiply by 0.001 (mmol to mol) and 14 (mol to gN)
  if(iunit=='mmol/m2'){
    sp.n$OrigValueStr[i]=iN*0.001*14
    #If the units are g/m2, then keep the same value
  }else{
    if(iunit=='g/m2'){
      sp.n$OrigValueStr[i]=iN
      #If the units are mg/cm2, then multiply by 10
    }else{
      if(iunit=='mg/cm2'){
        sp.n$OrigValueStr[i]=iN*10
        #If none of the above applies, then make the value NA
      }else{
        sp.n$OrigValueStr[i]=NA
      }
    }
  }
}
##Note that the highest Narea is for crested wheat--probably from ag field

#Check to make sure there are no NA
subset(sp.n,sp.n$OrigValueStr=='NA')

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean Narea for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.n,AccSpeciesID==i)
  #Subset by only trait 50 (= Narea)
  sp=subset(sp,TraitID==50)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  nit=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,5]=nit
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`Narea (g/m2)` = as.numeric(as.character(df$`Narea (g/m2)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`Narea (g/m2)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`Narea (g/m2)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Narea (g/m2)`))

#####################################
## Aarea and species number (7768) ##
#####################################

#Read in the dataframe
sp.a=fread('7768.txt')
#Subset only necessary columns
sp.a=subset(sp.a,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.a=subset(sp.a,TraitID!=0)

#Look at the measurements of Aarea
unique(sp.a$OriglName)
#Make a table of the frequency of the measurements of Aarea
unit=sort(table(sp.a$OriglName),decreasing=T)
#View frequency table
View(unit)
#The measurements appear to be pretty much equivalent

#Look at the units
unique(sp.a$OrigUnitStr)
#Make a table of the frequency of the units of Aarea
unit=sort(table(sp.a$OrigUnitStr),decreasing=T)
#View the frequency table
View(unit)

#There are a bunch of observations that have no values (blank rows)--remove
sp.a=subset(sp.a,sp.a$OrigUnitStr!="")

#Make the units consistent
#I made the decisions to assume that \xb5 refers to micro, and that m-1 should have been m-2
for(i in 1:length(sp.a$DatasetID)){
  #If units are mmolCO2/m2/s then standardize to mmolCO2/m2/s
  if(sp.a$OrigUnitStr[i]%in%c('micro mol/m2/s','micromol/m2/s','micromoles/m2/s','umol CO2 m-2 s-1','umol/m2/s','micro mol m-2 s-1','micromol m-2 s-1','umol CO2/m2/s','umolCO2/m2-s','\xb5mol(CO2) m-2 s-1','micro mol CO2 m-2 s-1','umol CO2 / m2 / sec','umol CO2/m^2 s','micromolco2 m-2 s-1','micromol CO2 m-2 s-1','micromol CO2 m-1s-1','micromol. m-2. s-1','micro mol C m-2 s-1')){
    sp.a$OrigUnitStr[i]='umolCO2/m2/s'
    #If units are g/cm2/d, then keep the same
  }else{
    if(sp.a$OrigUnitStr[i]%in%c('g/cm2/d','g/cm2/day')){
      sp.a$OrigUnitStr[i]='g/cm2/d'
      #If units are g/m2/d, then keep the same
    }else{
      if(sp.a$OrigUnitStr[i]%in%c('g/m2/day','g/m2/d','g m-2 day-1')){
        sp.a$OrigUnitStr[i]='g/m2/d'
        #If none of the above apply, then make units NA
      }else{
        sp.a$OrigUnitStr[i]=NA
      }
    }
  }
}

#Check the units after running the for loop
unique(sp.a$OrigUnitStr)

#Make Aarea values numeric
sp.a$OrigValueStr=as.numeric(sp.a$OrigValueStr)

#Convert all Aarea observations to umol CO2/m2/s
for(i in 1:length(sp.a$DatasetID)){
  #Indexes the units of each observation
  iunit=sp.a$OrigUnitStr[i]
  #Indexes the value of each observation
  iA=sp.a$OrigValueStr[i]
  #If the units are micromol CO2/m2/s, then keep the same value
  if(iunit=='umolCO2/m2/s'){
    sp.a$OrigValueStr[i]=iA
    #If the units are g/cm2/d, then multiply by 1000000 (mol to micromol),divide by 44 (g to mol),divide by 0.0001 (cm2 to m2), and divide by 86400 (days to seconds)
  }else{
    if(iunit=='g/cm2/d'){
      sp.a$OrigValueStr[i]=iA*1000000/44/0.0001/86400
      #If the units are g/m2/d, then multiply by 1000000 (mol to micromol),divide by 44 (g to mol),divide by 86400 (days to s)
    }else{
      if(iunit=='g/m2/d'){
        sp.a$OrigValueStr[i]=iA*1000000/44/86400
        #If none of the above, assign NA to the value
      }else{
        sp.a$OrigValueStr[i]=NA
      }
    }
  }
}

#Check to make sure there are no NA
subset(sp.a,sp.a$OrigValueStr=='NA')

#Use the same matrix as for the other observations
output=cbind(output,NA)

test = 0
#Find the mean Aarea for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.a,AccSpeciesID==i)
  #Subset by only trait 53 (= Aarea)
  sp=subset(sp,TraitID==53)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  pho=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,6]=pho
  test = test + 1
  print(test / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`Aarea (mmol CO2/m2/s` = as.numeric(as.character(df$`Aarea (mmol CO2/m2/s`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`Aarea (mmol CO2/m2/s`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`Aarea (mmol CO2/m2/s`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Aarea (mmol CO2/m2/s`))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Aarea (mmol CO2/m2/s`)) +
  coord_cartesian(ylim = c(0, 50))

#####################################################
## Stem specific density and species number (7769) ##
#####################################################

#Read in the dataframe
sp.ssd=fread('7769.txt')
#Subset only necessary columns
sp.ssd=subset(sp.ssd,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.ssd=subset(sp.ssd,TraitID!=0)

#Look at the measurements
unique(sp.ssd$OriglName)
#Make a table of the frequency of the measurements of SSD
unit=sort(table(sp.ssd$OriglName),decreasing=T)
#View the frequency table
View(unit)
#Remove inconsistent measurements
sp.ssd=subset(sp.ssd,sp.ssd$OriglName!='WoodDensityMinGCm-3')
sp.ssd=subset(sp.ssd,sp.ssd$OriglName!='WoodDensityMaxGCm-3')

#Look at the units
unique(sp.ssd$OrigUnitStr)
#Make a table of the frequency of the units of SSD
unit=sort(table(sp.ssd$OrigUnitStr),decreasing=T)
#View the frequency table
View(unit)

#I have no idea what the unit t m3 is. If figured out, I will add back in
sp.ssd=subset(sp.ssd,sp.ssd$OrigUnitStr!='t m3')

#Make the units consistent
#Here i assume that g/(cm3*10) is actually g/cm3
for(i in 1:length(sp.ssd$DatasetID)){
  if(sp.ssd$OrigUnitStr[i]%in%c('kg/m3','kg m-3')){
    sp.ssd$OrigUnitStr[i]='kg/m3'
  }else{
    if(sp.ssd$OrigUnitStr[i]%in%c('g/dm3')){
      sp.ssd$OrigUnitStr[i]='g/dm3'
    }else{
      if(sp.ssd$OrigUnitStr[i]%in%c('g / cm3','g/cm^3','g/cm3','g_cm3-1','g cm-3','g*cm-3','g /(cm3*10)')){
        sp.ssd$OrigUnitStr[i]='g/cm3'
      }else{
        if(sp.ssd$OrigUnitStr[i]%in%c('mg/mm3','mg/mm^3','mg / mm3')){
          sp.ssd$OrigUnitStr[i]='mg/mm3'
        }else{
          if(sp.ssd$OrigUnitStr[i]%in%c('cm3/mg')){
            sp.ssd$OrigUnitStr[i]='cm3/mg'
          }else{
            sp.ssd$OrigUnitStr[i]=NA
          }
        }
      }
    }
  }
}

#Check the units after running the for loop
unique(sp.ssd$OrigUnitStr)

#Make SSD values numeric
sp.ssd$OrigValueStr=as.numeric(sp.ssd$OrigValueStr)

#Convert all SSD observations to g/cm3
for(i in 1:length(sp.ssd$DatasetID)){
  iunit=sp.ssd$OrigUnitStr[i]
  issd=sp.ssd$OrigValueStr[i]
  if(iunit=='kg/m3'){
    sp.ssd$OrigValueStr[i]=issd*1000/1000000
  }else{
    if(iunit=='g/dm3'){
      sp.ssd$OrigValueStr[i]=issd/1000
    }else{
      if(iunit=='g/cm3'){
        sp.ssd$OrigValueStr[i]=issd
      }else{
        if(iunit=='mg/mm3'){
          sp.ssd$OrigValueStr[i]=issd
        }else{
          if(iunit=='cm3/mg'){
            sp.ssd$OrigValueStr[i]=1/(issd/0.001)
          }else{
            sp.ssd$OrigValueStr[i]=NA
          }
        }
      }
    }
  }
}
##Interesting to note that the highest value is for Festuca?? Seems fake

#Check to make sure there are no NA
subset(sp.ssd,sp.ssd$OrigValueStr=='NA')

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean SSD for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.ssd,AccSpeciesID==i)
  #Subset by only trait 4 (= SSD)
  sp=subset(sp,TraitID==4)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  ssd=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,7]=ssd
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`SSD (g/cm3)` = as.numeric(as.character(df$`SSD (g/cm3)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`SSD (g/cm3)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`SSD (g/cm3)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`SSD (g/cm3)`))

#############################################################
## Sapwood specific conductivity and species number (7772) ##
#############################################################

#Read in the dataframe
sp.ssc=fread('7772.txt')
#Subset only necessary columns
sp.ssc=subset(sp.ssc,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.ssc=subset(sp.ssc,TraitID!=0)

#Look at the measurements of SSC
unique(sp.ssc$OriglName)
#All the measurements appear to be the same

#Look at the units
unique(sp.ssc$OrigUnitStr)
#Make a table of the frequency of the units of SSC
unit=sort(table(sp.ssc$OrigUnitStr),decreasing=T)
#View the frequency table
View(unit)

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean SSC for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.ssc,AccSpeciesID==i)
  #Subset by only trait 1096 (= SSC)
  sp=subset(sp,TraitID==1096)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  ssc=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,8]=ssc
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`SSC (kg/m/s/MPa)` = as.numeric(as.character(df$`SSC (kg/m/s/MPa)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`SSC (kg/m/s/MPa)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`SSC (kg/m/s/MPa)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`SSC (kg/m/s/MPa)`))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`SSC (kg/m/s/MPa)`)) +
  coord_cartesian(ylim = c(0, 50))

####################################################
## Specific root length and species number (7946) ##
####################################################

#Read in the dataframe
sp.srl=fread('7946.txt')
#Subset only necessary columns
sp.srl=subset(sp.srl,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.srl=subset(sp.srl,TraitID!=0)

#Look at the measurements of SRL
unique(sp.srl$OriglName)
#Make a table of the frequency of the measurements of SRL
unit=sort(table(sp.srl$OriglName),decreasing=T)
#View the frequency table
View(unit)
#Remove inconsistent measurements
sp.srl=subset(sp.srl,sp.srl$OriglName!='Min_Specific root length (SRL)')
sp.srl=subset(sp.srl,sp.srl$OriglName!='Max_Specific root length (SRL)')

#Look at the units
unique(sp.srl$OrigUnitStr)
#Make a table of the frequency of the units of SRL
unit=sort(table(sp.srl$OrigUnitStr),decreasing=T)
#View the frequency table
View(unit)

#Make the units consistent
for(i in 1:length(sp.srl$DatasetID)){
  #If the units are cm/g then make them say cm/g
  if(sp.srl$OrigUnitStr[i]=='cm/g'){
    sp.srl$OrigUnitStr[i]='cm/g'
    #If the units are m/g, then make them say m/g
  }else{
    if(sp.srl$OrigUnitStr[i]=='m/g'){
      sp.srl$OrigUnitStr[i]='m/g'
      #If the units are cm/mg, then keep them cm/mg
    }else{
      if(sp.srl$OrigUnitStr[i]=='g/cm'){
        sp.srl$OrigUnitStr[i]='g/cm'
        #If none of the above are true, make the units NA
      }else{
        sp.srl$OrigUnitStr[i]=NA
      }
    }
  }
}

#Check the units after running the for loop
unique(sp.srl$OrigUnitStr)

#Remove observation of 0 for computational reasons
sp.srl=subset(sp.srl,sp.srl$OrigValueStr!=0)

#Make SSD values numeric
sp.srl$OrigValueStr=as.numeric(sp.srl$OrigValueStr)

#Convert all SRL values to m/g
for(i in 1:length(sp.srl$DatasetID)){
  #Index the units for each observation
  iunit=sp.srl$OrigUnitStr[i]
  #Index the value for each observation
  isrl=sp.srl$OrigValueStr[i]
  #If the units are cm/g, then multiply the value by 0.01
  if(iunit=='cm/g'){
    sp.srl$OrigValueStr[i]=isrl*0.01
    #IF the units are m/g, keep the value the same
  }else{
    if(iunit=='m/g'){
      sp.srl$OrigValueStr[i]=isrl
      #If the units are cm/mg, multiply the value by 0.01 and divide by 0.001
    }else{
      if(iunit=='g/cm'){
        sp.srl$OrigValueStr[i]=(1/isrl)*0.01
        #If none of the above are true, make the value NA
      }else{
        sp.srl$OrigValueStr[i]=NA
      }
    }
  }
}

#Check to make sure there are no NA
subset(sp.srl,sp.srl$OrigValueStr=='NA')

#Use the same matrix as for the other observations
output=cbind(output,NA)

#Find the mean SRL for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.srl,AccSpeciesID==i)
  #Subset by only trait 614 (= SRL)
  sp=subset(sp,TraitID==614)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  srl=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,9]=srl
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`SRL (m/g)` = as.numeric(as.character(df$`SRL (m/g)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`SRL (m/g)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`SRL (m/g)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`SRL (m/g)`))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y =df$`SRL (m/g)`)) +
  coord_cartesian(ylim = c(0, 500))

##################################################
## Fine root diameter and species number (7775) ##
##################################################

#Read in the dataframe
sp.frd=fread('7775.txt')
#Subset only necessary columns
sp.frd=subset(sp.frd,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.frd=subset(sp.frd,TraitID!=0)

#Look at the measurements
unique(sp.frd$OriglName)
#Make a table of the frequency of the measurements of FRD
unit=sort(table(sp.frd$OriglName),decreasing=T)
#View frequency table
View(unit)
#Remove inconsistent measurements
sp.frd=subset(sp.frd,sp.frd$OriglName!='Min_Root diameter')
sp.frd=subset(sp.frd,sp.frd$OriglName!='Max_Root diameter')
sp.frd=subset(sp.frd,sp.frd$OriglName!='Upper quartile_Root diameter')
sp.frd=subset(sp.frd,sp.frd$OriglName!='Lower quartile_Root diameter')
sp.frd=subset(sp.frd,sp.frd$OriglName!='Modal_Root diameter')
#To keep the greatest number of observations, mean and median FRD were assumed to be ~the same

#Look at the units
unique(sp.frd$OrigUnitStr)
#Make a table of the frequency of the units of FRD
unit=sort(table(sp.frd$OrigUnitStr),decreasing=T)
#View the frequency table
View(unit)

#There are some non-numeric values--find them to be able to assign them values
value=sort(table(sp.frd$OrigValueStr),decreasing=T)
#View frequency table
View(value)

#Convert all observations into a single value. Ranges are condensed to median value
for(i in 1:length(sp.frd$DatasetID)){
  #If the range is 100-150, make 125
  if(sp.frd$OrigValueStr[i]=='100-150'){
    sp.frd$OrigValueStr[i]=125
    #If the range is less than 100, make 100
  }else{
    if(sp.frd$OrigValueStr[i]=='<100'){
      sp.frd$OrigValueStr[i]=100
      #If the range is 150-200, make 175
    }else{
      if(sp.frd$OrigValueStr[i]=='150-200'){
        sp.frd$OrigValueStr[i]=175
        #If the range is 200-300, make 250
      }else{
        if(sp.frd$OrigValueStr[i]=='200-300'){
          sp.frd$OrigValueStr[i]=250
          #If the range is 400-500, make 400
        }else{
          if(sp.frd$OrigValueStr[i]=='300-500'){
            sp.frd$OrigValueStr[i]=400
            #All other value should stay the same (assumed to already be a single value)
          }else{
            sp.frd$OrigValueStr[i]=sp.frd$OrigValueStr[i]
          }
        }
      }
    }
  }
}

#Make FRD values numeric
sp.frd$OrigValueStr=as.numeric(sp.frd$OrigValueStr)

#Convert all FRD values to mm
for(i in 1:length(sp.frd$DatasetID)){
  #Index the units for each observation
  iunit=sp.frd$OrigUnitStr[i]
  #Index the value for each observation
  ifrd=sp.frd$OrigValueStr[i]
  #If the units are micrometers, multiply the value by 0.001
  if(iunit=='microm'){
    sp.frd$OrigValueStr[i]=ifrd*0.001
    #If the units are in mm, keep the value the same
  }else{
    if(iunit=='mm'){
      sp.frd$OrigValueStr[i]=ifrd
      #If the units are anything else, make the value NA
    }else{
      sp.frd$OrigValueStr[i]=NA
    }
  }
}

#Check to make sure there are no NA
subset(sp.frd,sp.frd$OrigValueStr=='NA')

#Use the same matrix as for the other observations
output=cbind(output,NA)

#Find the mean FRD for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.frd,AccSpeciesID==i)
  #Subset by only trait 896 (= FRD)
  sp=subset(sp,TraitID==896)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  frd=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,10]=frd
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`FRD (mm)` = as.numeric(as.character(df$`FRD (mm)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`FRD (mm)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`FRD (mm)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`FRD (mm)`))

####################################################
## Coarse root diameter and species number (7776) ##
####################################################

#Read in the dataframe
sp.crd=fread('7776.txt')
#Subset only necessary columns
sp.crd=subset(sp.crd,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.crd=subset(sp.crd,TraitID!=0)

#Look at the units
unique(sp.crd$OrigUnitStr)

#Make CRD values numeric
sp.crd$OrigValueStr=as.numeric(sp.crd$OrigValueStr)

#Use the same matrix as for the other observations
output=cbind(output,NA)

#Find the mean CRD for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.crd,AccSpeciesID==i)
  #Subset by only trait 1526 (= CRD)
  sp=subset(sp,TraitID==1526)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  crd=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,11]=crd
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`CRD (mm)` = as.numeric(as.character(df$`CRD (mm)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`CRD (mm)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`CRD (mm)`)))

################################################################
## Mycorrhizal infection intensity  and species number (7953) ##
################################################################

#Read in the dataframe
sp.mc=fread('7953.txt')
#Subset only necessary columns
sp.mc=subset(sp.mc,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.mc=subset(sp.mc,TraitID!=0)

#Look at the units
unique(sp.mc$OrigUnitStr)

#Look at the observations
unique(sp.mc$OriglName)
#All measurements look consistent here

#Look at values
unique(sp.mc$OrigValueStr)

#Some values are categorical (e.g. 'normally mycorrhizal')--remove
sp.mc=sp.mc[!is.na(as.numeric(sp.mc$OrigValueStr))]

#Make density values numeric
sp.mc$OrigValueStr=as.numeric(sp.mc$OrigValueStr)

#Use the same matrix as for the other observations
output=cbind(output,NA)

#Find the mean MC for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.mc,AccSpeciesID==i)
  #Subset by only trait 1030 (= MC)
  sp=subset(sp,TraitID==1030)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  mc=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,12]=mc
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)', 'MC (%)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`MC (%)` = as.numeric(as.character(df$`MC (%)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`MC (%)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`MC (%)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, df$`MC (%)`))

###########################################################
## Root hydraulic conductivity and species number (7778) ##
###########################################################

#Read in the dataframe
sp.rhd=fread('7778.txt')
#Subset only necessary columns
sp.rhd=subset(sp.rhd,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.rhd=subset(sp.rhd,TraitID!=0)

#Look at the units
unique(sp.rhd$OrigUnitStr)

#Make RHD values numeric
sp.rhd$OrigValueStr=as.numeric(sp.rhd$OrigValueStr)

#Use the same matrix as for the other observations
output=cbind(output,NA)

#Find the mean RHD for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.rhd,AccSpeciesID==i)
  #Subset by only trait 1476 (= RHD)
  sp=subset(sp,TraitID==1476)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  rhd=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,13]=rhd
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)','MC (%)','RHC (cm3/cm/s/MPa*10^6)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`RHC (cm3/cm/s/MPa*10^6)` = as.numeric(as.character(df$`RHC (cm3/cm/s/MPa*10^6)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`RHC (cm3/cm/s/MPa*10^6)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`RHC (cm3/cm/s/MPa*10^6)`)))

#####################################################
## Stem conduit diameter and species number (7949) ##
#####################################################

#Read in the dataframe
sp.wvd=fread('7949.txt')
#Subset only necessary columns
sp.wvd=subset(sp.wvd,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.wvd=subset(sp.wvd,TraitID!=0)

#Look at the units
unique(sp.wvd$OrigUnitStr)

#Look at measurements
unique(sp.wvd$OriglName)
#Make a table of the frequency of the units of SLA
unit=sort(table(sp.wvd$OriglName),decreasing=T)
#View the frequency table
View(unit)

#There are relatively few observations using maximum and minimum diameters
#I am removing them to try to keep these observations as consistent as possible
#And hope this will not reduce the sample size too much for lianas
sp.wvd=subset(sp.wvd,sp.wvd$OriglName!='Stem vessel Dmax (\xb5m)')
sp.wvd=subset(sp.wvd,sp.wvd$OriglName!='Stem maximum vessel diameter')
sp.wvd=subset(sp.wvd,sp.wvd$OriglName!='Max. vessel diameter (micrometer)')
sp.wvd=subset(sp.wvd,sp.wvd$OriglName!='Min. vessel diameter (micrometer)')

#Make WVD values numeric
sp.wvd$OrigValueStr=as.numeric(sp.wvd$OrigValueStr)

#Use the same matrix as for the other observations
output=cbind(output,NA)

#Find the mean WVD for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.wvd,AccSpeciesID==i)
  #Subset by only trait 281 (= SCDb)
  sp=subset(sp,TraitID==281)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  wvd=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,14]=wvd
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)','MC (%)','RHC (cm3/cm/s/MPa*10^6)','WVD (um)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`WVD (um)` = as.numeric(as.character(df$`WVD (um)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`WVD (um)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`WVD (um)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`WVD (um)`))

####################################################
## Stem conduit density and species number (7948) ##
####################################################

#Read in the dataframe
sp.dens=fread('7948.txt')
#Subset only necessary columns
sp.dens=subset(sp.dens,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.dens=subset(sp.dens,TraitID!=0)

#Look at the units
unique(sp.dens$OrigUnitStr)

#Look at the observations
unique(sp.dens$OriglName)
#Make a table of the frequency of the observation methods
unit=sort(table(sp.dens$OriglName),decreasing=T)
#View frequency table
View(unit)

#I will assume that, unless otherwise stated, the observations correspond to average vessel density
#Therefore, I will remove minimum and maximum observations and keep everything else the same
sp.dens=subset(sp.dens,sp.dens$OriglName!='Max. vessel density (number/mm2)')
sp.dens=subset(sp.dens,sp.dens$OriglName!='Min. vessel density (number/mm2)')

#Make density values numeric
sp.dens$OrigValueStr=as.numeric(sp.dens$OrigValueStr)

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean vessel density for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.dens,AccSpeciesID==i)
  #Subset by only trait 169 (= conduit density)
  sp=subset(sp,TraitID==169)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  dens=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,15]=dens
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)','MC (%)','RHC (cm3/cm/s/MPa*10^6)','WVD (um)','Dens (1/mm^2)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`Dens (1/mm^2)` = as.numeric(as.character(df$`Dens (1/mm^2)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`Dens (1/mm^2)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`Dens (1/mm^2)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Dens (1/mm^2)`))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Dens (1/mm^2)`)) +
  coord_cartesian(ylim = c(0, 1000))

############################################
## Plant height and species number (7849) ##
############################################

#Read in the dataframe
sp.h=fread('7849.txt')
#Subset only necessary columns
sp.h=subset(sp.h,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.h=subset(sp.h,TraitID!=0)

#Look at the measurements
unique(sp.h$OriglName)
#Make a table of the frequency of the units of H
unit=sort(table(sp.h$OriglName),decreasing=T)
#View the frequency table
View(unit)
#Remove inconsistent measurements
sp.h=subset(sp.h,sp.h$OriglName!='Height: 3. typical maximum (cm)')
sp.h=subset(sp.h,sp.h$OriglName!='Maximum Height')
sp.h=subset(sp.h,sp.h$OriglName!='MaximumHeight')
sp.h=subset(sp.h,sp.h$OriglName!='Maximum plant height (m)')
sp.h=subset(sp.h,sp.h$OriglName!='Plant height [max]')
sp.h=subset(sp.h,sp.h$OriglName!='Plant height [min]')
sp.h=subset(sp.h,sp.h$OriglName!='MaximumHeightMaxM')
sp.h=subset(sp.h,sp.h$OriglName!='MaximumHeightMinM')
sp.h=subset(sp.h,sp.h$OriglName!='MaximumHeightM')
sp.h=subset(sp.h,sp.h$OriglName!='Hmax(m)')
sp.h=subset(sp.h,sp.h$OriglName!='Height max')
sp.h=subset(sp.h,sp.h$OriglName!='HEIGHT min')
sp.h=subset(sp.h,sp.h$OriglName!='max height observed at the site (cm)')
sp.h=subset(sp.h,sp.h$OriglName!='maximum PlantHeight observed at the site (m)')
sp.h=subset(sp.h,sp.h$OriglName!='MVH')
sp.h=subset(sp.h,sp.h$OriglName!='MaximumHeightExtremeM')
sp.h=subset(sp.h,sp.h$OriglName!='Maximum height (cm)')
sp.h=subset(sp.h,sp.h$OriglName!='Height max (m)')
sp.h=subset(sp.h,sp.h$OriglName!='maximum height')
sp.h=subset(sp.h,sp.h$OriglName!='Stature (max height)')
sp.h=subset(sp.h,sp.h$OriglName!='Plant_height_vegetative_max')
sp.h=subset(sp.h,sp.h$OriglName!='Plant_height_vegetative_min')
sp.h=subset(sp.h,sp.h$OriglName!='maxH')
sp.h=subset(sp.h,sp.h$OriglName!='Maximum plant height')
sp.h=subset(sp.h,sp.h$OriglName!='Hpot95_m')
sp.h=subset(sp.h,sp.h$OriglName!='Max. Height (cm)')
sp.h=subset(sp.h,sp.h$OriglName!='Max_Height')
sp.h=subset(sp.h,sp.h$OriglName!='Max height')
sp.h=subset(sp.h,sp.h$OriglName!='Max. Height')
sp.h=subset(sp.h,sp.h$OriglName!='Max adult height (m)')
sp.h=subset(sp.h,sp.h$OriglName!='max.height')
sp.h=subset(sp.h,sp.h$OriglName!='Flowering plant height, heighest leaf elongated')
sp.h=subset(sp.h,sp.h$OriglName!='Flowering plant height, heighest leaf not elongated')
sp.h=subset(sp.h,sp.h$OriglName!='max height cpg')
sp.h=subset(sp.h,sp.h$OriglName!='Max ht')
sp.h=subset(sp.h,sp.h$OriglName!='WholePlant HMAX')
sp.h=subset(sp.h,sp.h$OriglName!='Maximum canopy height')
sp.h=subset(sp.h,sp.h$OriglName!='WholePlant kHMax')
sp.h=subset(sp.h,sp.h$OriglName!='Height (seedling)')
sp.h=subset(sp.h,sp.h$OriglName!='MaxHeight')

#Look at the units
unique(sp.h$OrigUnitStr)
#Make a table of the frequency of the units of H
unit=sort(table(sp.h$OrigUnitStr),decreasing=T)
#View the frequency table
View(unit)

#Because there were so many different categorical values, I removed everything that isn't a number
sp.h=sp.h[!is.na(as.numeric(sp.h$OrigValueStr))]

#Make H values numeric
sp.h$OrigValueStr=as.numeric(sp.h$OrigValueStr)

#Convert all observations to m
for(i in 1:length(sp.h$DatasetID)){
  iunit=sp.h$OrigUnitStr[i]
  ih=sp.h$OrigValueStr[i]
  if(iunit=='m'){
    sp.h$OrigValueStr[i]=ih
  }else{
    if(iunit=='cm'){
      sp.h$OrigValueStr[i]=ih*0.01
    }else{
      if(iunit=='feet'){
        sp.h$OrigValueStr[i]=ih*0.3048
      }else{
        if(iunit=='mm'){
          sp.h$OrigValueStr[i]=ih*0.001
        }else{
          sp.h$OrigValueStr[i]=NA
        }
      }
    }
  }
}

#Check to make sure there are no NA
subset(sp.h,sp.h$OrigValueStr=='NA')

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean H for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.h,AccSpeciesID==i)
  #Subset by only trait 3106 (= H)
  sp=subset(sp,TraitID==3106)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  hh=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,16]=hh
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)','MC (%)','RHC (cm3/cm/s/MPa*10^6)','WVD (um)','Dens (1/mm)','Height (m)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`Height (m)` = as.numeric(as.character(df$`Height (m)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`Height (m)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`Height (m)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Height (m)`))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Height (m)`)) +
  coord_cartesian(ylim = c(0, 50))

#########################################
## Leaf area and species number (7850) ##
#########################################

#Read in the dataframe
sp.la=fread('7850.txt')
#Subset only necessary columns
sp.la=subset(sp.la,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.la=subset(sp.la,TraitID!=0)

#Look at measurements
unique(sp.la$OriglName)
#Look at frequency of measurements of LA
unit=sort(table(sp.la$OriglName),decreasing=T)
#View frequency table
View(unit)
#The descriptions are strange but everything seems consistent

#Look at values
value=sort(table(sp.la$OrigValueStr),decreasing=T)
#View frequency table for values
View(value)

#Take out nonsense values
sp.la=subset(sp.la,sp.la$OrigValueStr!='no')
sp.la=subset(sp.la,sp.la$OrigValueStr!='yes')

#Look at the units
unique(sp.la$OrigUnitStr)
#Make a table of the frequency of the units of SLA
unit=sort(table(sp.la$OrigUnitStr),decreasing=T)
#View the frequency table
View(unit)

#Make SLA values numeric
sp.la$OrigValueStr=as.numeric(sp.la$OrigValueStr)

#Convert all observations to cm2
for(i in 1:length(sp.la$DatasetID)){
  iunit=sp.la$OrigUnitStr[i]
  ila=sp.la$OrigValueStr[i]
  if(iunit=='mm2'){
    sp.la$OrigValueStr[i]=ila*0.01
  }else{
    if(iunit=='cm2'){
      sp.la$OrigValueStr[i]=ila
    }else{
      if(iunit=='m2'){
        sp.la$OrigValueStr[i]=ila*10000
      }else{
        sp.la$OrigValueStr[i]=NA
      }
    }
  }
}

#Check to make sure there are no NA
subset(sp.la,sp.la$OrigValueStr=='NA')

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean H for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.la,AccSpeciesID==i)
  #Subset by only trait 3113 (= LA)
  sp=subset(sp,TraitID==3113)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  la=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,17]=la
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)','MC (%)','RHC (cm3/cm/s/MPa*10^6)','WVD (um)','Dens (1/mm)','Height (m)','LA (cm^2)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`LA (cm^2)` = as.numeric(as.character(df$`LA (cm^2)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`LA (cm^2)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`LA (cm^2)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`LA (cm^2)`))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`LA (cm^2)`)) +
  coord_cartesian(ylim = c(0, 250))

#########################################
## Seed mass and species number (7868) ##
#########################################

#Read in the dataframe
sp.sm=fread('7868.txt')
#Subset only necessary columns
sp.sm=subset(sp.sm,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.sm=subset(sp.sm,TraitID!=0)

#Look at the measurements
unique(sp.sm$OriglName)
#Make a table of the frequency of the measurements of SM
unit=sort(table(sp.sm$OriglName),decreasing=T)
#View the frequency table
View(unit)
#Remove inconsistent measurements
#I assume that anything that doesn't say dry mass is original mass (keep orginal mass)
sp.sm=subset(sp.sm,sp.sm$OriglName!='dry seed mass (mg)')
sp.sm=subset(sp.sm,sp.sm$OriglName!='Seed dry weight MIT')
sp.sm=subset(sp.sm,sp.sm$OriglName!='DSPR_DRY')
sp.sm=subset(sp.sm,sp.sm$OriglName!='SWT, oven-dry seed weight (mg)')
sp.sm=subset(sp.sm,sp.sm$OriglName!='')
sp.sm=subset(sp.sm,sp.sm$OriglName!='mean dry wght (mg)')
sp.sm=subset(sp.sm,sp.sm$OriglName!='SeedMassMax')
sp.sm=subset(sp.sm,sp.sm$OriglName!='SeedMassMin')
sp.sm=subset(sp.sm,sp.sm$OriglName!='Mean dry mass per seed (g)')
sp.sm=subset(sp.sm,sp.sm$OriglName!='OriginalSeedMassMax')
sp.sm=subset(sp.sm,sp.sm$OriglName!='OriginalSeedMassMin')
sp.sm=subset(sp.sm,sp.sm$OriglName!='seed dry weight')
sp.sm=subset(sp.sm,sp.sm$OriglName!='Seed dry mass')
sp.sm=subset(sp.sm,sp.sm$OriglName!='WholePlant SSIZE')

#Look at the units
unique(sp.sm$OrigUnitStr)
#Make a table of the frequency of the units of SM
unit=sort(table(sp.sm$OrigUnitStr),decreasing=T)
#View the frequency table
View(unit)

#Look at values
value=sort(table(sp.sm$OrigValueStr),decreasing=T)
#View frequency table for values
View(value)

#Take out non-numeric observations (e.g. "very light")
sp.sm=sp.sm[!is.na(as.numeric(sp.sm$OrigValueStr))]

#Make units consistent
for(i in 1:length(sp.sm$DatasetID)){
  if(sp.sm$OrigUnitStr[i]%in%c('mg','mg per seed')){
    sp.sm$OrigUnitStr[i]='mg'
  }else{
    if(sp.sm$OrigUnitStr[i]%in%c('g / 1000 seeds','g/1000')){
      sp.sm$OrigUnitStr[i]='g/1000'
    }else{
      if(sp.sm$OrigUnitStr[i]=='micro gx10'){
        sp.sm$OrigUnitStr[i]='ug*10'
      }else{
        if(sp.sm$OrigUnitStr[i]%in%c('1/pound','1/lb','1/lb (lb~454g)')){
          sp.sm$OrigUnitStr[i]='1/lb'
        }else{
          if(sp.sm$OrigUnitStr[i]%in%c('g','gr')){
            sp.sm$OrigUnitStr[i]='g'
          }else{
            if(sp.sm$OrigUnitStr[i]=='mg/100 seeds'){
              sp.sm$OrigUnitStr[i]='mg/100'
            }else{
              if(sp.sm$OrigUnitStr[i]=='1/kg'){
                sp.sm$OrigUnitStr[i]='1/kg'
              }else{
                if(sp.sm$OrigUnitStr[i]=='1/g'){
                  sp.sm$OrigUnitStr[i]='1/g'
                }else{
                  if(sp.sm$OrigUnitStr[i]=='1/liter (liter~1,7kg)'){
                    sp.sm$OrigUnitStr[i]='1/L'
                  }else{
                    if(sp.sm$OrigUnitStr[i]=='g/1000'){
                      sp.sm$OrigUnitStr[i]='g/1000'
                    }else{
                      sp.sm$OrigUnitStr[i]=NA
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

#Check the units after running the for loop
unique(sp.sm$OrigUnitStr)

#Make SM values numeric
sp.sm$OrigValueStr=as.numeric(sp.sm$OrigValueStr)

#Convert all observations to mg
for(i in 1:length(sp.sm$DatasetID)){
  iunit=sp.sm$OrigUnitStr[i]
  ism=sp.sm$OrigValueStr[i]
  if(iunit=='mg'){
    sp.sm$OrigValueStr[i]=ism
  }else{
    if(iunit=='g/1000'){
      sp.sm$OrigValueStr[i]=ism/1000*1000
    }else{
      if(iunit=='ug*10'){
        sp.sm$OrigValueStr[i]=ism*0.001
      }else{
        if(iunit=='1/lb'){
          sp.sm$OrigValueStr[i]=(1/ism)*453592
        }else{
          if(iunit=='g'){
            sp.sm$OrigValueStr[i]=ism*1000
          }else{
            if(iunit=='mg/100'){
              sp.sm$OrigValueStr[i]=ism/100
            }else{
              if(iunit=='1/kg'){
                sp.sm$OrigValueStr[i]=(1/ism)*1000000
              }else{
                if(iunit=='1/g'){
                  sp.sm$OrigValueStr[i]=(1/ism)*1000
                }else{
                  if(iunit=='1/L'){
                    sp.sm$OrigValueStr[i]=(1/ism)*1.7*1000000
                  }else{
                    sp.sm$OrigValueStr[i]=NA
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

#Check to make sure there are no NA
subset(sp.sm,sp.sm$OrigValueStr=='NA')

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean SM for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.sm,AccSpeciesID==i)
  #Subset by only trait 26 (= SM)
  sp=subset(sp,TraitID==26)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  sm=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,18]=sm
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)','MC (%)','RHC (cm3/cm/s/MPa*10^6)','WVD (um)','Dens (1/mm^2)','Height (m)','LA (cm^2)','SM (mg)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`SM (mg)` = as.numeric(as.character(df$`SM (mg)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`SM (mg)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`SM (mg)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`SM (mg)`))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`SM (mg)`)) +
  coord_cartesian(ylim = c(0, 100))

#######################################
## Wood vessel element length (7950) ##
#######################################

#Read in the dataframe
sp.wvel=fread('7950.txt')
#Subset only necessary columns
sp.wvel=subset(sp.wvel,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.wvel=subset(sp.wvel,TraitID!=0)

#Look at the measurements
unique(sp.wvel$OriglName)
#Take out inconcsistent measurements
sp.wvel=subset(sp.wvel,sp.wvel$OriglName=='Av. vessel element length (micrometer)')

#Look at the units
unique(sp.wvel$OrigUnitStr)

#Make WVEL values numeric
sp.wvel$OrigValueStr=as.numeric(sp.wvel$OrigValueStr)

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean WVEL for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.wvel,AccSpeciesID==i)
  #Subset by only trait 282 (= WVEL)
  sp=subset(sp,TraitID==282)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  wvel=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,19]=wvel
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)','MC (%)','RHC (cm3/cm/s/MPa*10^6)','WVD (um)','Dens (1/mm^2)','Height (m)','LA (cm^2)','SM (mg)','WVEL (um)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`WVEL (um)` = as.numeric(as.character(df$`WVEL (um)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`WVEL (um)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`WVEL (um)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`WVEL (um)`))

###################################################
## Root rooting depth  and species number (7951) ##
###################################################

#Read in the dataframe
sp.rrd=fread('7951.txt')
#Subset only necessary columns
sp.rrd=subset(sp.rrd,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.rrd=subset(sp.rrd,TraitID!=0)

#Look at the measurements
unique(sp.rrd$OriglName)
#Make a table of the frequency of the measurements of RRD
unit=sort(table(sp.rrd$OriglName),decreasing=T)
#View the frequency table
View(unit)
#Remove inconsistent measurements
sp.rrd=subset(sp.rrd,sp.rrd$OriglName!='Root Depth, Minimum')
sp.rrd=subset(sp.rrd,sp.rrd$OriglName!='Max Rooting Depth (m)')
sp.rrd=subset(sp.rrd,sp.rrd$OriglName!='Rooting depth_Extrapolated 95 percent rooting depth')
sp.rrd=subset(sp.rrd,sp.rrd$OriglName!='Rooting depth_Extrapolated 50 percent rooting depth')
sp.rrd=subset(sp.rrd,sp.rrd$OriglName!='Rooting depth_Interpolated 95 percent rooting depth')
sp.rrd=subset(sp.rrd,sp.rrd$OriglName!='Rooting depth_min')
sp.rrd=subset(sp.rrd,sp.rrd$OriglName!='Rooting depth_max')
sp.rrd=subset(sp.rrd,sp.rrd$OriglName!='Rooting depth in the field (categorical)')

#Look at the units
unique(sp.rrd$OrigUnitStr)
#Make a table of the frequency of the units of RRD
unit=sort(table(sp.rrd$OrigUnitStr),decreasing=T)
#View the frequency table
View(unit)

#The unit 'text' literally means qualitative observation--remove
sp.rrd=subset(sp.rrd,sp.rrd$OrigUnitStr!='text')

#Make RRD values numeric
sp.rrd$OrigValueStr=as.numeric(sp.rrd$OrigValueStr)

#Convert all RRD observations to m
for(i in 1:length(sp.rrd$DatasetID)){
  iunit=sp.rrd$OrigUnitStr[i]
  irrd=sp.rrd$OrigValueStr[i]
  if(iunit=='m'){
    sp.rrd$OrigValueStr[i]=irrd
  }else{
    if(iunit=='cm'){
      sp.rrd$OrigValueStr[i]=irrd*0.01
    }else{
      if(iunit=='inches'){
        sp.rrd$OrigValueStr[i]=irrd*0.0254
      }else{
        sp.rrd$OrigValueStr[i]=NA
      }
    }
  }
}

#Check to make sure there are no NA
subset(sp.rrd,sp.rrd$OrigValueStr=='NA')

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean RRD for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.rrd,AccSpeciesID==i)
  #Subset by only trait 6 (= RRD)
  sp=subset(sp,TraitID==6)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  rrd=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,20]=rrd
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)','MC (%)','RHC (cm3/cm/s/MPa*10^6)','WVD (um)','Dens (1/mm^2)','Height (m)','LA (cm^2)','SM (mg)','WVEL (um)','RRD (m)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`RRD (m)` = as.numeric(as.character(df$`RRD (m)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`RRD (m)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`RRD (m)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`RRD (m)`))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`RRD (m)`)) +
  coord_cartesian(ylim = c(0, 20))

#######################################################
## Hydraulic vulnerability and species number (7975) ##
#######################################################

#Read in the dataframe
sp.hv=fread('7975.txt')
#Subset only necessary columns
sp.hv=subset(sp.hv,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName','Comment','DataName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.hv=subset(sp.hv,TraitID!=0)

#Look at the measurements
unique(sp.hv$OriglName)

#This trait is really a bunch of traits put together
#Separate P50==xylem pressure at which 50% of conductivity is lost
sp.p50=subset(sp.hv,sp.hv$OriglName==c('P50','P50 (MPa)','Xylem tension at 50% loss of hydraulic conductivity (MPa)','Choat et al. 2012 reported ?50 (MPa)','Mean ?50 (with all) (MPa)','Water potential at 50% loss of conductivity Psi_50 (MPa)'))
#Separate P88==xylem pressure at which 88% of conductivity is lost
sp.p88=subset(sp.hv,sp.hv$OriglName==c('P88','P88 (MPa)'))
#Separate Pmin==the most negative xylem pressure reported
sp.pmin=subset(sp.hv,sp.hv$OriglName==c('Pmin midday leaf','Pmin midday (Mpa)'))
#Separate P50 safety margin==difference between P50 and Pmin (lower number = higher risk strategy)
sp.p50sm=subset(sp.hv,sp.hv$OriglName==c('P50 safety margin'))
#Separate P88 safety margin==difference between P88 and Pmin
sp.p88sm=subset(sp.hv,sp.hv$OriglName==c('P88 safety margin'))

#Because P50 has more observations than the other traits, I'm continuing with that one

#Make P50 observations numeric
sp.p50$OrigValueStr=as.numeric(sp.p50$OrigValueStr)

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean P50 for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.p50,AccSpeciesID==i)
  #Subset by only trait 719 (= P50)
  sp=subset(sp,TraitID==719)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  p50=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,21]=p50
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)','MC (%)','RHC (cm3/cm/s/MPa*10^6)','WVD (um)','Dens (1/mm^2)','Height (m)','LA (cm^2)','SM (mg)','WVEL (um)','RRD (m)','P50 (MPa)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`P50 (MPa)` = as.numeric(as.character(df$`P50 (MPa)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`P50 (MPa)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`P50 (MPa)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`P50 (MPa)`))

############################################
## Leaf N (mass based) and species number ##
############################################

#Read in the dataframe
sp.nmass=fread('8096.txt')
#Subset only necessary columns
sp.nmass=subset(sp.nmass,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.nmass=subset(sp.nmass,TraitID!=0)

#Look at the measurements
unique(sp.nmass$OriglName)
#Make a frequency table of measurements
unit=sort(table(sp.nmass$OriglName),decreasing=T)
#View frequency table
View(unit)
#Remove inconsistent units
sp.nmass=subset(sp.nmass,sp.nmass$OriglName!='N_senesced_leaf')
sp.nmass=subset(sp.nmass,sp.nmass$OriglName!='Leaf nitrogen concentration predicted from NIRS')
sp.nmass=subset(sp.nmass,sp.nmass$OriglName!='LeafN Min (mg/mg *100)')
sp.nmass=subset(sp.nmass,sp.nmass$OriglName!='LeafN Max (mg/mg *100)')

#Look at the units
unique(sp.nmass$OrigUnitStr)
#Make a frequency table of units
unit=sort(table(sp.nmass$OrigUnitStr),decreasing=T)
#Look at frequency table
View(unit)
#Remove the rows with no observation
sp.nmass=subset(sp.nmass,sp.nmass$OrigUnitStr!='')
#Make units consistent
for(i in 1:length(sp.nmass$DatasetID)){
  if(sp.nmass$OrigUnitStr[i]%in%c('mg/g','mg g-1','mg_g-1','\xb5g mg-1','g/kg','mg/g dry mass','mg / g','mg N g-1','g kg-1')){
    sp.nmass$OrigUnitStr[i]='mg/g'
  }else{
    if(sp.nmass$OrigUnitStr[i]%in%c('mmol/g')){
      sp.nmass$OrigUnitStr[i]='mmol/g'
    }else{
      if(sp.nmass$OrigUnitStr[i]%in%c('g N g-1 DW','g/g','kg/kg')){
        sp.nmass$OrigUnitStr[i]='g/g'
      }else{
        if(sp.nmass$OrigUnitStr[i]%in%c('%','mg/mg *100','percent','% mass/mass')){
          sp.nmass$OrigUnitStr[i]='%'
        }else{
          sp.nmass$OrigUnitStr[i]=NA
        }
      }
    }
  }
}

#Check to make sure there are no NA
subset(sp.nmass,sp.nmass$OrigValueStr=='NA')

#There were some non-numeric values...remove (all were read as being less than zero)
sp.nmass=subset(sp.nmass,sp.nmass$OrigValueStr>0)

#Make Nmass values numeric
sp.nmass$OrigValueStr=as.numeric(sp.nmass$OrigValueStr)

#Convert all observations to mg/g
for(i in 1:length(sp.nmass$DatasetID)){
  iunit=sp.nmass$OrigUnitStr[i]
  inmass=sp.nmass$OrigValueStr[i]
  if(iunit=='mg/g'){
    sp.nmass$OrigValueStr[i]=inmass
  }else{
    if(iunit=='%'){
      sp.nmass$OrigValueStr[i]=inmass/100/0.001
    }else{
      if(iunit=='mmol/g'){
        sp.nmass$OrigValueStr[i]=inmass*0.001*14*1000
      }else{
        if(iunit=='g/g'){
          sp.nmass$OrigValueStr[i]=inmass*1000
        }else{
          sp.nmass$OrigValueStr[i]=NA
        }
      }
    }
  }
}

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean Nmass for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.nmass,AccSpeciesID==i)
  #Subset by only trait 14 (= Nmass)
  sp=subset(sp,TraitID==14)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  nmass=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,22]=nmass
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)','MC (%)','RHC (cm3/cm/s/MPa*10^6)','WVD (um)','Dens (1/mm^2)','Height (m)','LA (cm^2)','SM (mg)','WVEL (um)','RRD (m)','P50 (MPa)','Nmass (mg/g)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`Nmass (mg/g)` = as.numeric(as.character(df$`Nmass (mg/g)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`Nmass (mg/g)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`Nmass (mg/g)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Nmass (mg/g)`))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Nmass (mg/g)`)) +
  coord_cartesian(ylim = c(0, 100))

############################################
## Leaf P (mass based) and species number ##
############################################

#Read in the dataframe
sp.pmass=fread('8097.txt')
#Subset only necessary columns
sp.pmass=subset(sp.pmass,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond with leaf lifespan or growth type
sp.pmass=subset(sp.pmass,TraitID!=0)

#Look at the measurements
unique(sp.pmass$OriglName)
#Make a frequency table of measurements
unit=sort(table(sp.pmass$OriglName),decreasing=T)
#View frequency table
View(unit)
#Remove inconsistent units
sp.pmass=subset(sp.pmass,sp.pmass$OriglName!='P_senesced_leaf')

#Look at the units
unique(sp.pmass$OrigUnitStr)
#Make a frequency table of units
unit=sort(table(sp.pmass$OrigUnitStr),decreasing=T)
#Look at frequency table
View(unit)
#Remove the rows with no observation
sp.pmass=subset(sp.pmass,sp.pmass$OrigUnitStr!='')
#Remove ppm (only 3 obervations)
sp.pmass=subset(sp.pmass,sp.pmass$OrigUnitStr!='ppm')

#Make units consistent
for(i in 1:length(sp.pmass$DatasetID)){
  if(sp.pmass$OrigUnitStr[i]%in%c('mg/g','mg g-1','mg_g-1','g/kg')){
    sp.pmass$OrigUnitStr[i]='mg/g'
  }else{
    if(sp.pmass$OrigUnitStr[i]%in%c('%','percent','% mass/mass')){
      sp.pmass$OrigUnitStr[i]='%'
    }else{
      if(sp.pmass$OrigUnitStr[i]%in%c('g P g-1 DW','g/g')){
        sp.pmass$OrigUnitStr[i]='g/g'
      }else{
        if(sp.pmass$OrigUnitStr[i]%in%c('mg kg-1','mg/kg')){
          sp.pmass$OrigUnitStr[i]='mg/kg'
        }else{
          if(sp.pmass$OrigUnitStr[i]%in%c('mg/10g')){
            sp.pmass$OrigUnitStr[i]='mg/10g'
          }else{
            if(sp.pmass$OrigUnitStr[i]%in%c('mmol/kg')){
              sp.pmass$OrigUnitStr[i]='mmol/kg'
            }else{
              sp.pmass$OrigUnitStr[i]=NA
            }
          }
        }
      }
    }
  }
}

#Check to make sure there are no NA
subset(sp.pmass,sp.pmass$OrigValueStr=='NA')

#There were some non-numeric values (all were read as being less than zero)
sp.pmass=subset(sp.pmass,sp.pmass$OrigValueStr!='')
sp.pmass=subset(sp.pmass,sp.pmass$OrigValueStr>0)

#Make Pmass values numeric
sp.pmass$OrigValueStr=as.numeric(sp.pmass$OrigValueStr)

#Convert all observations to mg/g
for(i in 1:length(sp.pmass$DatasetID)){
  iunit=sp.pmass$OrigUnitStr[i]
  ipmass=sp.pmass$OrigValueStr[i]
  if(iunit=='mg/g'){
    sp.pmass$OrigValueStr[i]=ipmass
  }else{
    if(iunit=='%'){
      sp.pmass$OrigValueStr[i]=ipmass*1000/100
    }else{
      if(iunit=='g/g'){
        sp.pmass$OrigValueStr[i]=ipmass*1000
      }else{
        if(iunit=='mg/kg'){
          sp.pmass$OrigValueStr[i]=ipmass/1000
        }else{
          if(iunit=='mg/10g'){
            sp.pmass$OrigValueStr[i]=ipmass/10
          }else{
            if(iunit=='mmol/kg'){
              sp.pmass$OrigValueStr[i]=ipmass*0.001*15
            }else{
              sp.pmass$OrigValueStr[i]=NA
            }
          }
        }
      }
    }
  }
}

#Check for NA
subset(sp.pmass,sp.pmass$OrigValueStr=='NA')

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean Pmass for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.pmass,AccSpeciesID==i)
  #Subset by only trait 15 (= Pmass)
  sp=subset(sp,TraitID==15)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  pmass=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,23]=pmass
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)','MC (%)','RHC (cm3/cm/s/MPa*10^6)','WVD (um)','Dens (1/mm^2)','Height (m)','LA (cm^2)','SM (mg)','WVEL (um)','RRD (m)','P50 (MPa)','Nmass (mg/g)','Pmass (mg/g)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`Pmass (mg/g)` = as.numeric(as.character(df$`Pmass (mg/g)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`Pmass (mg/g)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`Pmass (mg/g)`)))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Pmass (mg/g)`))

ggplot() +
  geom_violin(aes(x = df$`Growth form`, y = df$`Pmass (mg/g)`)) +
  coord_cartesian(ylim = c(0, 10))

####################################
## Crown depth and species number ##
####################################

#Read in the dataframe
sp.ch=fread('8520.txt')
#Subset only necessary columns
sp.ch=subset(sp.ch,select=c('DatasetID','AccSpeciesID','ObservationID','TraitID','OrigValueStr','OrigUnitStr','OriglName'))
#Take out all the rows that do not correspond
sp.ch=subset(sp.ch,TraitID!=0)

#Look at the measurements
unique(sp.ch$OriglName)
#Make a frequency table of measurements
unit=sort(table(sp.ch$OriglName),decreasing=T)
#View frequency table
View(unit)

#Look at the units
unique(sp.ch$OrigUnitStr)
#Make a frequency table of units
unit=sort(table(sp.ch$OrigUnitStr),decreasing=T)
#Look at frequency table
View(unit)
#Remove the rows with no observation
sp.ch=subset(sp.ch,sp.ch$OrigUnitStr!='%')

#Make CH values numeric
sp.ch$OrigValueStr=as.numeric(sp.ch$OrigValueStr)

#Use the same matrix as for the other observations
output=cbind(output,NA)

count = 0
#Find the mean CH for each species and put into the output matrix
for(i in uniquespecies){
  #Subset by species ID number
  sp=subset(sp.ch,AccSpeciesID==i)
  #Subset by only trait 15 (= Pmass)
  sp=subset(sp,TraitID==773)
  #Make the values numeric
  sp$OrigValueStr=as.numeric(as.character(sp$OrigValueStr))
  #Find the mean for each species
  ch=mean(sp$OrigValueStr)
  #Put each observation in the output
  output[i,24]=ch
  count = count + 1
  print(count / length(uniquespecies) * 100)
}

#Add column names to output
colnames=c('Species ID','Growth form','Leaf lifespan (m)','SLA (mm2/mg)','Narea (g/m2)','Aarea (mmol CO2/m2/s','SSD (g/cm3)','SSC (kg/m/s/MPa)','SRL (m/g)','FRD (mm)','CRD (mm)','MC (%)','RHC (cm3/cm/s/MPa*10^6)','WVD (um)','Dens (1/mm^2)','Height (m)','LA (cm^2)','SM (mg)','WVEL (um)','RRD (m)','P50 (MPa)','Nmass (mg/g)','Pmass (mg/g)','Crown depth (m)')
colnames(output)=colnames

df = as.data.frame(output)
is.na(df) = is.na(df)
df$`Crown depth (m)` = as.numeric(as.character(df$`Crown depth (m)`))

length(which(df$`Growth form` == 'tree' & !is.na(df$`Crown depth (m)`)))
length(which(df$`Growth form` == 'liana' & !is.na(df$`Crown depth (m)`)))

write.csv(output, file = 'output_9oct.csv')
