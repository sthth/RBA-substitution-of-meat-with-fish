setwd("//dtu-storage/sthth/Documents/Case 1/Data")


#packages
library(data.table)
library(mc2d)
library(fitdistrplus)
library(goftest)


#settings
iters_var <- 100000
iters_unc <- 1000

ndvar(iters_var)
ndunc(iters_unc)


#population statistics
SEYLL_15_19 <- 74.54
SEYLL_20_24 <- 69.57
SEYLL_25_29 <- 64.6
SEYLL_30_34 <- 59.63
SEYLL_35_39 <- 54.67
SEYLL_40_44 <- 49.73
SEYLL_45_49 <- 44.81
SEYLL_50_54 <- 39.92
SEYLL_55_59 <- 35.07
SEYLL_60_64 <- 30.25
SEYLL_65_69 <- 25.49
SEYLL_70_74 <- 20.77
SEYLL_75_79 <- 16.43
SEYLL_80_84 <- 12.51
SEYLL_85 <- 7.6

LE_15_19m <- 79.07
LE_20_24m <- 79.16
LE_25_29m <- 79.28
LE_30_34m <- 79.41
LE_35_39m <- 79.56
LE_40_44m <- 79.77
LE_45_49m <- 80.09
LE_50_54m <- 80.54
LE_55_59m <- 81.23
LE_60_64m <- 82.17
LE_65_69m <- 83.40
LE_70_74m <- 84.87
LE_75_79m <- 86.66
LE_80_84m <- 88.90
LE_85m <- 95.43


LE_15_19w <- 82.94
LE_20_24w <- 82.99
LE_25_29w <- 83.04
LE_30_34w <- 83.11
LE_35_39w <- 83.19
LE_40_44w <- 83.33
LE_45_49w <- 83.52
LE_50_54w <- 83.84
LE_55_59w <- 84.33
LE_60_64w <- 85.01
LE_65_69w <- 85.92
LE_70_74w <- 87.02
LE_75_79w <- 88.43
LE_80_84w <- 90.35
LE_85w <- 96.10


pop_15_19m <- 181089 #Population no. for men by age
pop_20_24m <- 192234
pop_25_29m <- 176161
pop_30_34m <- 161526
pop_35_39m <- 180041 
pop_40_44m <- 195490
pop_45_49m <- 209022
pop_50_54m <- 199100
pop_55_59m <- 177815
pop_60_64m <- 166416
pop_65_69m <- 172808
pop_70_74m <- 130931
pop_75_79m <- 85272
pop_80_84m <- 51232
pop_85m <- 38575 #For 85+

pop_men <- pop_15_19m + pop_20_24m + pop_25_29m + pop_30_34m + pop_35_39m + pop_40_44m + pop_45_49m + pop_50_54m +
  pop_55_59m + pop_60_64m + pop_65_69m + pop_70_74m + pop_75_79m + pop_80_84m + pop_85m

pop_15_19w <- 171815 #for women
pop_20_24w <- 184806
pop_25_29w <- 170197
pop_30_34w <- 158343
pop_35_39w <- 178824 
pop_40_44w <- 194114
pop_45_49w <- 204582
pop_50_54w <- 195997
pop_55_59w <- 177966
pop_60_64w <- 170401
pop_65_69w <- 179072
pop_70_74w <- 142515
pop_75_79w <- 101122
pop_80_84w <- 70698
pop_85w <- 78904

pop_women <- pop_15_19w + pop_20_24w + pop_25_29w + pop_30_34w + pop_35_39w + pop_40_44w + pop_45_49w + pop_50_54w +
  pop_55_59w + pop_60_64w + pop_65_69w + pop_70_74w + pop_75_79w + pop_80_84w + pop_85w


#Stomach cancer incidence/mortality in Denmark, 2015

## Background risk of non-cardia stomach cancer for different age and sex groups - non-cardia assumed to comprise
## 1/3 of total stomach cancer cases (ref: Danish Cancer Society)

p15_19m <- 0
p20_24m <- 0
p25_29m <- 0
p30_34m <- 0.00001 * 0.3
p35_39m <- 0.00001 * 0.3
p40_44m <- 0.00004 * 0.3
p45_49m <- 0.00006 * 0.3
p50_54m <- 0.00007 * 0.3
p55_59m <- 0.00019 * 0.3
p60_64m <- 0.00025 * 0.3
p65_69m <- 0.00045 * 0.3
p70_74m <- 0.00052 * 0.3
p75_79m <- 0.00053 * 0.3
p80_84m <- 0.00115 * 0.3
p85m <- 0.00074 * 0.3


p15_19w <- 0
p20_24w <- 0.00001 * 0.3
p25_29w <- 0.00001 * 0.3
p30_34w <- 0
p35_39w <- 0
p40_44w <- 0.00002 * 0.3
p45_49w <- 0.00005 * 0.3
p50_54w <- 0.00004 * 0.3
p55_59w <- 0.00008 * 0.3
p60_64w <- 0.00014 * 0.3
p65_69w <- 0.00019 * 0.3
p70_74w <- 0.00021 * 0.3
p75_79w <- 0.00041 * 0.3
p80_84w <- 0.00043 * 0.3
p85w <- 0.00023 * 0.3

#Assume same mortality rate as other stomach cancers
fatal15_19m <- 0
fatal20_24m <- 0
fatal25_29m <- 0
fatal30_34m <- 0
fatal35_39m <- 1.1
fatal40_44m <- 0.5
fatal45_49m <- 0.166666667
fatal50_54m <- 0.857142857
fatal55_59m <- 0.557894737
fatal60_64m <- 0.696
fatal65_69m <- 0.52
fatal70_74m <- 0.538461538
fatal75_79m <- 0.756603774
fatal80_84m <- 0.58173913
fatal85m <- 1.006756757


fatal15_19w <- 0
fatal20_24w <- 0
fatal25_29w <- 0.6
fatal30_34w <- 0
fatal35_39w <- 0
fatal40_44w <- 0.75
fatal45_49w <- 0.6
fatal50_54w <- 0.5
fatal55_59w <- 0.625
fatal60_64w <- 0.507142857
fatal65_69w <- 0.473684211
fatal70_74w <- 0.519047619
fatal75_79w <- 0.497560976
fatal80_84w <- 0.844186047
fatal85w <- 0.939130435

#Stomach cancer natural history model
time_diagnosis <- 0.5
time_remission_cure <- 8
time_remission_death <- 0.267
time_disseminated <- 0.25
time_terminal <- 0.083
time_total_YLL <- time_diagnosis + time_remission_death + time_disseminated + time_terminal
time_total_YLD <- time_diagnosis + time_remission_cure


#disability weights (Salomon et al., 2015)
set.seed(1)
dw_diagnosis <- mcstoc(rpert, type = "U", min=0.193, mode=0.288, max=0.339)
dw_remission <- 0.2
dw_disseminated <- mcstoc(rpert, type = "U", min=0.307, mode=0.451, max=0.6)
dw_terminal <- mcstoc(rpert, type = "U", min=0.377, mode=0.54, max=0.687)


##### Ref scenario #####

FoodDaily <- read.csv("FoodDaily.csv")
MeatDaily <- FoodDaily[, c(1:4,34,45,58)]

agebreaks <- c(15,20,25,30,35,40,45,50,55,60,65,70,75,80)
agelabels= c("15-19","20-24","25-29","30-34","35-39", "40-44", "45-49", "50-54","55-59", "60-64","65-69", "70-74","75-79")

setDT(MeatDaily)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### DALY calc males #####

# Zero incidence for males 15-29 years of age

#Age 30-34, males

YLD <- p30_34m * dw_diagnosis * time_diagnosis + #Diagnosis
  p30_34m * (1-fatal30_34m) * dw_remission * time_remission_cure + #Remission, cure
  p30_34m * fatal30_34m * dw_remission * time_remission_death + #Remission, death
  p30_34m * fatal30_34m * dw_disseminated * time_disseminated + #Disseminated
  p30_34m * fatal30_34m * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p30_34m * fatal30_34m * (SEYLL_30_34 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY30_34totalm <- DALY * pop_30_34m


#Age 35-39, males

YLD <- p35_39m * dw_diagnosis * time_diagnosis + #Diagnosis
  p35_39m * (1-fatal35_39m) * dw_remission * time_remission_cure + #Remission, cure
  p35_39m * fatal35_39m * dw_remission * time_remission_death + #Remission, death
  p35_39m * fatal35_39m * dw_disseminated * time_disseminated + #Disseminated
  p35_39m * fatal35_39m * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p35_39m * fatal35_39m * (SEYLL_35_39 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY35_39totalm <- DALY * pop_35_39m


#Age 40-44, males

YLD <- p40_44m * dw_diagnosis * time_diagnosis + #Diagnosis
  p40_44m * (1-fatal40_44m) * dw_remission * time_remission_cure + #Remission, cure
  p40_44m * fatal40_44m * dw_remission * time_remission_death + #Remission, death
  p40_44m * fatal40_44m * dw_disseminated * time_disseminated + #Disseminated
  p40_44m * fatal40_44m * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p40_44m * fatal40_44m * (SEYLL_40_44 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY40_44totalm <- DALY * pop_40_44m



#Age 45-49, males

YLD <- p45_49m * dw_diagnosis * time_diagnosis + #Diagnosis
  p45_49m * (1-fatal45_49m) * dw_remission * time_remission_cure + #Remission, cure
  p45_49m * fatal45_49m * dw_remission * time_remission_death + #Remission, death
  p45_49m * fatal45_49m * dw_disseminated * time_disseminated + #Disseminated
  p45_49m * fatal45_49m * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p45_49m * fatal45_49m * (SEYLL_45_49 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY45_49totalm <- DALY * pop_45_49m


#Age 50-54, males

YLD <- p50_54m * dw_diagnosis * time_diagnosis + #Diagnosis
  p50_54m * (1-fatal50_54m) * dw_remission * time_remission_cure + #Remission, cure
  p50_54m * fatal50_54m * dw_remission * time_remission_death + #Remission, death
  p50_54m * fatal50_54m * dw_disseminated * time_disseminated + #Disseminated
  p50_54m * fatal50_54m * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p50_54m * fatal50_54m * (SEYLL_50_54 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY50_54totalm <- DALY * pop_50_54m


#Age 55-59, males

YLD <- p55_59m * dw_diagnosis * time_diagnosis + #Diagnosis
  p55_59m * (1-fatal55_59m) * dw_remission * time_remission_cure + #Remission, cure
  p55_59m * fatal55_59m * dw_remission * time_remission_death + #Remission, death
  p55_59m * fatal55_59m * dw_disseminated * time_disseminated + #Disseminated
  p55_59m * fatal55_59m * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p55_59m * fatal55_59m * (SEYLL_55_59 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY55_59totalm <- DALY * pop_55_59m


#Age 60-64, males

YLD <- p60_64m * dw_diagnosis * time_diagnosis + #Diagnosis
  p60_64m * (1-fatal60_64m) * dw_remission * time_remission_cure + #Remission, cure
  p60_64m * fatal60_64m * dw_remission * time_remission_death + #Remission, death
  p60_64m * fatal60_64m * dw_disseminated * time_disseminated + #Disseminated
  p60_64m * fatal60_64m * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p60_64m * fatal60_64m * (SEYLL_60_64 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY60_64totalm <- DALY * pop_60_64m


#Age 65-69, males

YLD <- p65_69m * dw_diagnosis * time_diagnosis + #Diagnosis
  p65_69m * (1-fatal65_69m) * dw_remission * time_remission_cure + #Remission, cure
  p65_69m * fatal65_69m * dw_remission * time_remission_death + #Remission, death
  p65_69m * fatal65_69m * dw_disseminated * time_disseminated + #Disseminated
  p65_69m * fatal65_69m * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p65_69m * fatal65_69m * (SEYLL_65_69 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY65_69totalm <- DALY * pop_65_69m



#Age 70-74, males

YLD <- p70_74m * dw_diagnosis * time_diagnosis + #Diagnosis
  p70_74m * (1-fatal70_74m) * dw_remission * time_remission_cure + #Remission, cure
  p70_74m * fatal70_74m * dw_remission * time_remission_death + #Remission, death
  p70_74m * fatal70_74m * dw_disseminated * time_disseminated + #Disseminated
  p70_74m * fatal70_74m * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p70_74m * fatal70_74m * (SEYLL_70_74 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY70_74totalm <- DALY * pop_70_74m


#Age 75-79, males

YLD <- p75_79m * dw_diagnosis * time_diagnosis + #Diagnosis
  p75_79m * (1-fatal75_79m) * dw_remission * time_remission_cure + #Remission, cure
  p75_79m * fatal75_79m * dw_remission * time_remission_death + #Remission, death
  p75_79m * fatal75_79m * dw_disseminated * time_disseminated + #Disseminated
  p75_79m * fatal75_79m * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p75_79m * fatal75_79m * (SEYLL_75_79 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY75_79totalm <- DALY * pop_75_79m


#Age 80-84, males

YLD <- p80_84m * dw_diagnosis * time_diagnosis + #Diagnosis
  p80_84m * (1-fatal80_84m) * dw_remission * time_remission_cure + #Remission, cure
  p80_84m * fatal80_84m * dw_remission * time_remission_death + #Remission, death
  p80_84m * fatal80_84m * dw_disseminated * time_disseminated + #Disseminated
  p80_84m * fatal80_84m * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p80_84m * fatal80_84m * (SEYLL_80_84 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY80_84totalm <- DALY * pop_80_84m


#Age 85+, males

YLD <- p85m * dw_diagnosis * time_diagnosis + #Diagnosis
  p85m * (1-fatal85m) * dw_remission * time_remission_cure + #Remission, cure
  p85m * fatal85m * dw_remission * time_remission_death + #Remission, death
  p85m * fatal85m * dw_disseminated * time_disseminated + #Disseminated
  p85m * fatal85m * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p85m * fatal85m * (SEYLL_85 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY85totalm <- DALY * pop_85m


##### DALY calc females #####

# Zero incidence for females 15-19 and 30-39 years of age

#Age 20-24, females

YLD <- p20_24w * dw_diagnosis * time_diagnosis + #Diagnosis
  p20_24w * (1-fatal20_24w) * dw_remission * time_remission_cure + #Remission, cure
  p20_24w * fatal20_24w * dw_remission * time_remission_death + #Remission, death
  p20_24w * fatal20_24w * dw_disseminated * time_disseminated + #Disseminated
  p20_24w * fatal20_24w * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p20_24w * fatal20_24w * (SEYLL_20_24 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY20_24totalw <- DALY * pop_20_24w


#Age 25-29, females

YLD <- p25_29w * dw_diagnosis * time_diagnosis + #Diagnosis
  p25_29w * (1-fatal25_29w) * dw_remission * time_remission_cure + #Remission, cure
  p25_29w * fatal25_29w * dw_remission * time_remission_death + #Remission, death
  p25_29w * fatal25_29w * dw_disseminated * time_disseminated + #Disseminated
  p25_29w * fatal25_29w * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p25_29w * fatal25_29w * (SEYLL_25_29 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY25_29totalw <- DALY * pop_25_29w


#Age 40-44, females

YLD <- p40_44w * dw_diagnosis * time_diagnosis + #Diagnosis
  p40_44w * (1-fatal40_44w) * dw_remission * time_remission_cure + #Remission, cure
  p40_44w * fatal40_44w * dw_remission * time_remission_death + #Remission, death
  p40_44w * fatal40_44w * dw_disseminated * time_disseminated + #Disseminated
  p40_44w * fatal40_44w * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p40_44w * fatal40_44w * (SEYLL_40_44 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY40_44totalw <- DALY * pop_40_44w


#Age 45-49, females

YLD <- p45_49w * dw_diagnosis * time_diagnosis + #Diagnosis
  p45_49w * (1-fatal45_49w) * dw_remission * time_remission_cure + #Remission, cure
  p45_49w * fatal45_49w * dw_remission * time_remission_death + #Remission, death
  p45_49w * fatal45_49w * dw_disseminated * time_disseminated + #Disseminated
  p45_49w * fatal45_49w * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p45_49w * fatal45_49w * (SEYLL_45_49 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY45_49totalw <- DALY * pop_45_49w



#Age 50-54, females

YLD <- p50_54w * dw_diagnosis * time_diagnosis + #Diagnosis
  p50_54w * (1-fatal50_54w) * dw_remission * time_remission_cure + #Remission, cure
  p50_54w * fatal50_54w * dw_remission * time_remission_death + #Remission, death
  p50_54w * fatal50_54w * dw_disseminated * time_disseminated + #Disseminated
  p50_54w * fatal50_54w * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p50_54w * fatal50_54w * (SEYLL_50_54 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY50_54totalw <- DALY * pop_50_54w


#Age 55-59, females

YLD <- p55_59w * dw_diagnosis * time_diagnosis + #Diagnosis
  p55_59w * (1-fatal55_59w) * dw_remission * time_remission_cure + #Remission, cure
  p55_59w * fatal55_59w * dw_remission * time_remission_death + #Remission, death
  p55_59w * fatal55_59w * dw_disseminated * time_disseminated + #Disseminated
  p55_59w * fatal55_59w * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p55_59w * fatal55_59w * (SEYLL_55_59 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY55_59totalw <- DALY * pop_55_59w



#Age 60-64, females

YLD <- p60_64w * dw_diagnosis * time_diagnosis + #Diagnosis
  p60_64w * (1-fatal60_64w) * dw_remission * time_remission_cure + #Remission, cure
  p60_64w * fatal60_64w * dw_remission * time_remission_death + #Remission, death
  p60_64w * fatal60_64w * dw_disseminated * time_disseminated + #Disseminated
  p60_64w * fatal60_64w * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p60_64w * fatal60_64w * (SEYLL_60_64 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY60_64totalw <- DALY * pop_60_64w


#Age 65-69, females

YLD <- p65_69w * dw_diagnosis * time_diagnosis + #Diagnosis
  p65_69w * (1-fatal65_69w) * dw_remission * time_remission_cure + #Remission, cure
  p65_69w * fatal65_69w * dw_remission * time_remission_death + #Remission, death
  p65_69w * fatal65_69w * dw_disseminated * time_disseminated + #Disseminated
  p65_69w * fatal65_69w * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p65_69w * fatal65_69w * (SEYLL_65_69 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY65_69totalw <- DALY * pop_65_69w


#Age 70-74, females

YLD <- p70_74w * dw_diagnosis * time_diagnosis + #Diagnosis
  p70_74w * (1-fatal70_74w) * dw_remission * time_remission_cure + #Remission, cure
  p70_74w * fatal70_74w * dw_remission * time_remission_death + #Remission, death
  p70_74w * fatal70_74w * dw_disseminated * time_disseminated + #Disseminated
  p70_74w * fatal70_74w * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p70_74w * fatal70_74w * (SEYLL_70_74 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY70_74totalw <- DALY * pop_70_74w


#Age 75-79, females

YLD <- p75_79w * dw_diagnosis * time_diagnosis + #Diagnosis
  p75_79w * (1-fatal75_79w) * dw_remission * time_remission_cure + #Remission, cure
  p75_79w * fatal75_79w * dw_remission * time_remission_death + #Remission, death
  p75_79w * fatal75_79w * dw_disseminated * time_disseminated + #Disseminated
  p75_79w * fatal75_79w * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p75_79w * fatal75_79w * (SEYLL_75_79 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY75_79totalw <- DALY * pop_75_79w


#Age 80-84, females

YLD <- p80_84w * dw_diagnosis * time_diagnosis + #Diagnosis
  p80_84w * (1-fatal80_84w) * dw_remission * time_remission_cure + #Remission, cure
  p80_84w * fatal80_84w * dw_remission * time_remission_death + #Remission, death
  p80_84w * fatal80_84w * dw_disseminated * time_disseminated + #Disseminated
  p80_84w * fatal80_84w * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p80_84w * fatal80_84w * (SEYLL_80_84 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY80_84totalw <- DALY * pop_80_84w


#Age 85+, females

YLD <- p85w * dw_diagnosis * time_diagnosis + #Diagnosis
  p85w * (1-fatal85w) * dw_remission * time_remission_cure + #Remission, cure
  p85w * fatal85w * dw_remission * time_remission_death + #Remission, death
  p85w * fatal85w * dw_disseminated * time_disseminated + #Disseminated
  p85w * fatal85w * dw_terminal * time_terminal #Terminal

#YLLs
YLL <- p85w * fatal85w * (SEYLL_85 - time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY85totalw <- DALY * pop_85w

###### Total DALYs ref #####

DALYtotalref <- DALY30_34totalm + DALY35_39totalm +  DALY40_44totalm + DALY45_49totalm + DALY50_54totalm +
  DALY55_59totalm + DALY60_64totalm + DALY65_69totalm +  DALY70_74totalm + DALY75_79totalm + DALY80_84totalm +
  DALY85totalm + 
  DALY20_24totalw + DALY25_29totalw + DALY40_44totalw + DALY45_49totalw + DALY50_54totalw + DALY55_59totalw +
  DALY60_64totalw + DALY65_69totalw + DALY70_74totalw + DALY75_79totalw + DALY80_84totalw + DALY85totalw

summary(DALYtotalref)



###### Total cases reference scenario #####

cases_ref <- p30_34m * pop_30_34m + p35_39m * pop_35_39m + p40_44m * pop_40_44m + p45_49m * pop_45_49m +
  p50_54m * pop_50_54m + p55_59m * pop_55_59m + p60_64m * pop_60_64m + p65_69m * pop_65_69m + p70_74m * pop_70_74m +
  p75_79m * pop_75_79m + p80_84m * pop_80_84m + p85m * pop_85m +
  p20_24w * pop_20_24w + p25_29w * pop_25_29w + p40_44w * pop_40_44w + p45_49w * pop_45_49w + p50_54w * pop_50_54w +
  p55_59w * pop_55_59w + p60_64w * pop_60_64w + p65_69w * pop_65_69w + p70_74w * pop_70_74w + p75_79w * pop_75_79w +
  p80_84w * pop_80_84w + p85w * pop_85w

summary(cases_ref)



##### Fitting exposures males #####



#Divide data into agegroups and fit exposure to distribution

#risk of SC = 0 for men 15-29 years

#Age 30-34, male
Proc30_34m <- subset(MeatDaily, agegroups=="30-34" & sex =="1", select = proc.meat)

#Probability of intake
probproc30_34m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc30_34m==0)/length(Proc30_34m$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc30_34m!=0)/length(Proc30_34m$proc.meat) #probability of consumption
probproc30_34m[1,] <- c(pr.no,pr.yes)

Proc30_34posm <- Proc30_34m[which(Proc30_34m!=0)]
Proc30_34posm <- as.vector(Proc30_34posm$proc.meat)

fit1_30_34m <- fitdist(Proc30_34posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc30_34posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_30_34m$estimate[1], sdlog=fit1_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34m, fitnames="lnorm") 

fit2_30_34m <- fitdist(Proc30_34posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc30_34posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Anderson-Darling test

# gamma has a bit better fit that lognormal

#Age 35-39, male
Proc35_39m <- subset(MeatDaily, agegroups=="35-39" & sex =="1", select = proc.meat)

#Probability of intake
probproc35_39m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc35_39m==0)/length(Proc35_39m$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc35_39m!=0)/length(Proc35_39m$proc.meat) #probability of consumption
probproc35_39m[1,] <- c(pr.no,pr.yes)

Proc35_39posm <- Proc35_39m[which(Proc35_39m!=0)]
Proc35_39posm <- as.vector(Proc35_39posm$proc.meat)

fit1_35_39m <- fitdist(Proc35_39posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc35_39posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_35_39m$estimate[1], sdlog=fit1_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39m, fitnames="lnorm") 

fit2_35_39m <- fitdist(Proc35_39posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc35_39posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Anderson-Darling test

# gamma has much better fit but lognormal still ok

#Age 40-44, male
Proc40_44m <- subset(MeatDaily, agegroups=="40-44" & sex =="1", select = proc.meat)

#Probability of intake
probproc40_44m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc40_44m==0)/length(Proc40_44m$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc40_44m!=0)/length(Proc40_44m$proc.meat) #probability of consumption
probproc40_44m[1,] <- c(pr.no,pr.yes)

Proc40_44posm <- Proc40_44m[which(Proc40_44m!=0)]
Proc40_44posm <- as.vector(Proc40_44posm$proc.meat)


fit1_40_44m <- fitdist(Proc40_44posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc40_44posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_40_44m$estimate[1], sdlog=fit1_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44m, fitnames="lnorm") 

fit2_40_44m <- fitdist(Proc40_44posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc40_44posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Anderson-Darling test

# gamma, lognormal below 0.05 but above 0.01

#Age 45-49, male
Proc45_49m <- subset(MeatDaily, agegroups=="45-49" & sex =="1", select = proc.meat)

#Probability of intake
probproc45_49m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc45_49m==0)/length(Proc45_49m$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc45_49m!=0)/length(Proc45_49m$proc.meat) #probability of consumption
probproc45_49m[1,] <- c(pr.no,pr.yes)

Proc45_49posm <- Proc45_49m[which(Proc45_49m!=0)]
Proc45_49posm <- as.vector(Proc45_49posm$proc.meat)


fit1_45_49m <- fitdist(Proc45_49posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc45_49posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_45_49m$estimate[1], sdlog=fit1_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49m, fitnames="lnorm") 

fit2_45_49m <- fitdist(Proc45_49posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc45_49posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Anderson-Darling test

# gamma has much better fit but lognormal above 0.05

#Age 50-54, male
Proc50_54m <- subset(MeatDaily, agegroups=="50-54" & sex =="1", select = proc.meat)

#Probability of intake
probproc50_54m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc50_54m==0)/length(Proc50_54m$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc50_54m!=0)/length(Proc50_54m$proc.meat) #probability of consumption
probproc50_54m[1,] <- c(pr.no,pr.yes)

Proc50_54posm <- Proc50_54m[which(Proc50_54m!=0)]
Proc50_54posm <- as.vector(Proc50_54posm$proc.meat)


fit1_50_54m <- fitdist(Proc50_54posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc50_54posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_50_54m$estimate[1], sdlog=fit1_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54m, fitnames="lnorm") 

fit2_50_54m <- fitdist(Proc50_54posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc50_54posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54m, fitnames="gamma") 


cvm.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Anderson-Darling test

# gamma has the best fit, lognormal below 0.05 but above 0.01

#Age 55-59, male
Proc55_59m <- subset(MeatDaily, agegroups=="55-59" & sex =="1", select = proc.meat)

#Probability of intake
probproc55_59m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc55_59m==0)/length(Proc55_59m$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc55_59m!=0)/length(Proc55_59m$proc.meat) #probability of consumption
probproc55_59m[1,] <- c(pr.no,pr.yes)

Proc55_59posm <- Proc55_59m[which(Proc55_59m!=0)]
Proc55_59posm <- as.vector(Proc55_59posm$proc.meat)


fit1_55_59m <- fitdist(Proc55_59posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc55_59posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_55_59m$estimate[1], sdlog=fit1_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59m, fitnames="lnorm") 

fit2_55_59m <- fitdist(Proc55_59posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc55_59posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59m, fitnames="gamma") 


cvm.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Anderson-Darling test

# gamma has the best fit, lognormal above 0.05

#Age 60-64, male
Proc60_64m <- subset(MeatDaily, agegroups=="60-64" & sex =="1", select = proc.meat)

#Probability of intake
probproc60_64m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc60_64m==0)/length(Proc60_64m$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc60_64m!=0)/length(Proc60_64m$proc.meat) #probability of consumption
probproc60_64m[1,] <- c(pr.no,pr.yes)

Proc60_64posm <- Proc60_64m[which(Proc60_64m!=0)]
Proc60_64posm <- as.vector(Proc60_64posm$proc.meat)


fit1_60_64m <- fitdist(Proc60_64posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc60_64posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_60_64m$estimate[1], sdlog=fit1_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64m, fitnames="lnorm") 

fit2_60_64m <- fitdist(Proc60_64posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc60_64posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Anderson-Darling test

# gamma has the best fit but lognormal above 0.05

#Age 65-69, male
Proc65_69m <- subset(MeatDaily, agegroups=="65-69" & sex =="1", select = proc.meat)

#Probability of intake
probproc65_69m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc65_69m==0)/length(Proc65_69m$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc65_69m!=0)/length(Proc65_69m$proc.meat) #probability of consumption
probproc65_69m[1,] <- c(pr.no,pr.yes)

Proc65_69posm <- Proc65_69m[which(Proc65_69m!=0)]
Proc65_69posm <- as.vector(Proc65_69posm$proc.meat)

fit1_65_69m <- fitdist(Proc65_69posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc65_69posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_65_69m$estimate[1], sdlog=fit1_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69m, fitnames="lnorm") 

fit2_65_69m <- fitdist(Proc65_69posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc65_69posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Anderson-Darling test

#  gamma has a much better fit but lognormal above 0.05

#Age 70-74, male
Proc70_74m <- subset(MeatDaily, agegroups=="70-74" & sex =="1", select = proc.meat)

#Probability of intake
probproc70_74m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc70_74m==0)/length(Proc70_74m$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc70_74m!=0)/length(Proc70_74m$proc.meat) #probability of consumption
probproc70_74m[1,] <- c(pr.no,pr.yes)

Proc70_74posm <- Proc70_74m[which(Proc70_74m!=0)]
Proc70_74posm <- as.vector(Proc70_74posm$proc.meat)

fit1_70_74m <- fitdist(Proc70_74posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc70_74posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_70_74m$estimate[1], sdlog=fit1_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74m, fitnames="lnorm") 

fit2_70_74m <- fitdist(Proc70_74posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc70_74posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Anderson-Darling test

# gamma has the best fit but lognormal above 0.05

#Age 75-79, male
Proc75_79m <- subset(MeatDaily, agegroups=="75-79" & sex =="1", select = proc.meat)

#Probability of intake
probproc75_79m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc75_79m==0)/length(Proc75_79m$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc75_79m!=0)/length(Proc75_79m$proc.meat) #probability of consumption
probproc75_79m[1,] <- c(pr.no,pr.yes)

Proc75_79posm <- Proc75_79m[which(Proc75_79m!=0)]
Proc75_79posm <- as.vector(Proc75_79posm$proc.meat)

fit1_75_79m <- fitdist(Proc75_79posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc75_79posm
plot(ecdf(t), lty=1)
x <- seq(0, 2000, length=1000)
lines(x, plnorm(x, meanlog=fit1_75_79m$estimate[1], sdlog=fit1_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79m, fitnames="lnorm") 

fit2_75_79m <- fitdist(Proc75_79posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc75_79posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Anderson-Darling test

# gamma and lognormal equally good


##### Fitting exposures females #####

#risk of SC = 0 for women 15-19 years


#Age 20-24, female
Proc20_24w <- subset(MeatDaily, agegroups=="20-24" & sex =="2", select = proc.meat)

#Probability of intake
probproc20_24w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc20_24w==0)/length(Proc20_24w$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc20_24w!=0)/length(Proc20_24w$proc.meat) #probability of consumption
probproc20_24w[1,] <- c(pr.no,pr.yes)

Proc20_24posw <- Proc20_24w[which(Proc20_24w!=0)]
Proc20_24posw <- as.vector(Proc20_24posw$proc.meat)

fit1_20_24w <- fitdist(Proc20_24posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc20_24posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_20_24w$estimate[1], sdlog=fit1_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24w, fitnames="lnorm") 

fit2_20_24w <- fitdist(Proc20_24posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc20_24posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Anderson-Darling test

# lognormal very good but gamma also above 0.05

#Age 25-29, female
Proc25_29w <- subset(MeatDaily, agegroups=="25-29" & sex =="2", select = proc.meat)

#Probability of intake
probproc25_29w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc25_29w==0)/length(Proc25_29w$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc25_29w!=0)/length(Proc25_29w$proc.meat) #probability of consumption
probproc25_29w[1,] <- c(pr.no,pr.yes)

Proc25_29posw <- Proc25_29w[which(Proc25_29w!=0)]
Proc25_29posw <- as.vector(Proc25_29posw$proc.meat)

fit1_25_29w <- fitdist(Proc25_29posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc25_29posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_25_29w$estimate[1], sdlog=fit1_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29w, fitnames="lnorm") 

fit2_25_29w <- fitdist(Proc25_29posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc25_29posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Anderson-Darling test

# gamma has the best fit but lognormal above 0.05



#risk of SC = 0 for women 30-39 years


#Age 40-44, female
Proc40_44w <- subset(MeatDaily, agegroups=="40-44" & sex =="2", select = proc.meat)

#Probability of intake
probproc40_44w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc40_44w==0)/length(Proc40_44w$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc40_44w!=0)/length(Proc40_44w$proc.meat) #probability of consumption
probproc40_44w[1,] <- c(pr.no,pr.yes)

Proc40_44posw <- Proc40_44w[which(Proc40_44w!=0)]
Proc40_44posw <- as.vector(Proc40_44posw$proc.meat)


fit1_40_44w <- fitdist(Proc40_44posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc40_44posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_40_44w$estimate[1], sdlog=fit1_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44w, fitnames="lnorm") 

fit2_40_44w <- fitdist(Proc40_44posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc40_44posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Anderson-Darling test

# gamma, lognormal above 0.05

#Age 45-49, female
Proc45_49w <- subset(MeatDaily, agegroups=="45-49" & sex =="2", select = proc.meat)

#Probability of intake
probproc45_49w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc45_49w==0)/length(Proc45_49w$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc45_49w!=0)/length(Proc45_49w$proc.meat) #probability of consumption
probproc45_49w[1,] <- c(pr.no,pr.yes)

Proc45_49posw <- Proc45_49w[which(Proc45_49w!=0)]
Proc45_49posw <- as.vector(Proc45_49posw$proc.meat)


fit1_45_49w <- fitdist(Proc45_49posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc45_49posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_45_49w$estimate[1], sdlog=fit1_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49w, fitnames="lnorm") 

fit2_45_49w <- fitdist(Proc45_49posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc45_49posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Anderson-Darling test

# gamma has much better fit but lognormal above 0.05

#Age 50-54, female
Proc50_54w <- subset(MeatDaily, agegroups=="50-54" & sex =="2", select = proc.meat)

#Probability of intake
probproc50_54w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc50_54w==0)/length(Proc50_54w$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc50_54w!=0)/length(Proc50_54w$proc.meat) #probability of consumption
probproc50_54w[1,] <- c(pr.no,pr.yes)

Proc50_54posw <- Proc50_54w[which(Proc50_54w!=0)]
Proc50_54posw <- as.vector(Proc50_54posw$proc.meat)


fit1_50_54w <- fitdist(Proc50_54posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc50_54posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_50_54w$estimate[1], sdlog=fit1_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54w, fitnames="lnorm") 

fit2_50_54w <- fitdist(Proc50_54posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc50_54posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54w, fitnames="gamma") 


cvm.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Anderson-Darling test

# gamma has the best fit, lognormal above 0.05

#Age 55-59, female
Proc55_59w <- subset(MeatDaily, agegroups=="55-59" & sex =="2", select = proc.meat)

#Probability of intake
probproc55_59w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc55_59w==0)/length(Proc55_59w$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc55_59w!=0)/length(Proc55_59w$proc.meat) #probability of consumption
probproc55_59w[1,] <- c(pr.no,pr.yes)

Proc55_59posw <- Proc55_59w[which(Proc55_59w!=0)]
Proc55_59posw <- as.vector(Proc55_59posw$proc.meat)


fit1_55_59w <- fitdist(Proc55_59posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc55_59posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_55_59w$estimate[1], sdlog=fit1_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59w, fitnames="lnorm") 

fit2_55_59w <- fitdist(Proc55_59posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc55_59posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59w, fitnames="gamma") 


cvm.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Anderson-Darling test

# gamma has the best fit, lognormal above 0.05

#Age 60-64, female
Proc60_64w <- subset(MeatDaily, agegroups=="60-64" & sex =="2", select = proc.meat)

#Probability of intake
probproc60_64w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc60_64w==0)/length(Proc60_64w$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc60_64w!=0)/length(Proc60_64w$proc.meat) #probability of consumption
probproc60_64w[1,] <- c(pr.no,pr.yes)

Proc60_64posw <- Proc60_64w[which(Proc60_64w!=0)]
Proc60_64posw <- as.vector(Proc60_64posw$proc.meat)


fit1_60_64w <- fitdist(Proc60_64posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc60_64posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_60_64w$estimate[1], sdlog=fit1_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64w, fitnames="lnorm") 

fit2_60_64w <- fitdist(Proc60_64posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc60_64posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Anderson-Darling test

# gamma has the best fit but lognormal above 0.05

#Age 65-69, female
Proc65_69w <- subset(MeatDaily, agegroups=="65-69" & sex =="2", select = proc.meat)

#Probability of intake
probproc65_69w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc65_69w==0)/length(Proc65_69w$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc65_69w!=0)/length(Proc65_69w$proc.meat) #probability of consumption
probproc65_69w[1,] <- c(pr.no,pr.yes)

Proc65_69posw <- Proc65_69w[which(Proc65_69w!=0)]
Proc65_69posw <- as.vector(Proc65_69posw$proc.meat)

fit1_65_69w <- fitdist(Proc65_69posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc65_69posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_65_69w$estimate[1], sdlog=fit1_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69w, fitnames="lnorm") 

fit2_65_69w <- fitdist(Proc65_69posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc65_69posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Anderson-Darling test

# gamma has a much better fit but lognormal above 0.05

#Age 70-74, female
Proc70_74w <- subset(MeatDaily, agegroups=="70-74" & sex =="2", select = proc.meat)

#Probability of intake
probproc70_74w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc70_74w==0)/length(Proc70_74w$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc70_74w!=0)/length(Proc70_74w$proc.meat) #probability of consumption
probproc70_74w[1,] <- c(pr.no,pr.yes)

Proc70_74posw <- Proc70_74w[which(Proc70_74w!=0)]
Proc70_74posw <- as.vector(Proc70_74posw$proc.meat)

fit1_70_74w <- fitdist(Proc70_74posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc70_74posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_70_74w$estimate[1], sdlog=fit1_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74w, fitnames="lnorm") 

fit2_70_74w <- fitdist(Proc70_74posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc70_74posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Anderson-Darling test

# gamma has the best fit but lognormal above 0.05

#Age 75-79, female
Proc75_79w <- subset(MeatDaily, agegroups=="75-79" & sex =="2", select = proc.meat)

#Probability of intake
probproc75_79w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc75_79w==0)/length(Proc75_79w$proc.meat) #probability of zero consumption
pr.yes <- sum(Proc75_79w!=0)/length(Proc75_79w$proc.meat) #probability of consumption
probproc75_79w[1,] <- c(pr.no,pr.yes)

Proc75_79posw <- Proc75_79w[which(Proc75_79w!=0)]
Proc75_79posw <- as.vector(Proc75_79posw$proc.meat)

fit1_75_79w <- fitdist(Proc75_79posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc75_79posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_75_79w$estimate[1], sdlog=fit1_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79w, fitnames="lnorm") 

fit2_75_79w <- fitdist(Proc75_79posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc75_79posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Anderson-Darling test

# gamma has the best fit, but lognormal above 0.05


##### Scenario 1 #####

Scenario1 <- read.csv("Scenario1.csv") # Since we're considering processed meat as a whole, the results for this endpoint will be the same
# for scenario 1,2,3,4

Alt_scenario <- Scenario1[,c(1:4,8:10)]
setDT(Alt_scenario)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Fitting exposures males #####

#risk of SC = 0 for men 15-29 years


#Age 30-34, male
Proc30_34m <- subset(Alt_scenario, agegroups=="30-34" & sex =="1", select = procmeat.new)

#Probability of intake
probproc_alt30_34m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc30_34m==0)/length(Proc30_34m$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc30_34m!=0)/length(Proc30_34m$procmeat.new) #probability of consumption
probproc_alt30_34m[1,] <- c(pr.no,pr.yes)

Proc30_34posm <- Proc30_34m[which(Proc30_34m!=0)]
Proc30_34posm <- as.vector(Proc30_34posm$procmeat.new)

fit1_alt_30_34m <- fitdist(Proc30_34posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc30_34posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_30_34m$estimate[1], sdlog=fit1_alt_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_30_34m, fitnames="lnorm") 

fit2_alt_30_34m <- fitdist(Proc30_34posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc30_34posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_alt_30_34m$estimate[1], fit2_alt_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_30_34m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_30_34m$estimate[1], fit1_alt_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_30_34m$estimate[1], fit1_alt_30_34m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_30_34m$estimate[1], fit2_alt_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_30_34m$estimate[1], fit2_alt_30_34m$estimate[2])  #Anderson-Darling test

# gamma has a bit better fit that lognormal

#Age 35-39, male
Proc35_39m <- subset(Alt_scenario, agegroups=="35-39" & sex =="1", select = procmeat.new)

#Probability of intake
probproc_alt35_39m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc35_39m==0)/length(Proc35_39m$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc35_39m!=0)/length(Proc35_39m$procmeat.new) #probability of consumption
probproc_alt35_39m[1,] <- c(pr.no,pr.yes)

Proc35_39posm <- Proc35_39m[which(Proc35_39m!=0)]
Proc35_39posm <- as.vector(Proc35_39posm$procmeat.new)

fit1_alt_35_39m <- fitdist(Proc35_39posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc35_39posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_35_39m$estimate[1], sdlog=fit1_alt_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_35_39m, fitnames="lnorm") 

fit2_alt_35_39m <- fitdist(Proc35_39posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc35_39posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_alt_35_39m$estimate[1], fit2_alt_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_35_39m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_35_39m$estimate[1], fit1_alt_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_35_39m$estimate[1], fit1_alt_35_39m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_35_39m$estimate[1], fit2_alt_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_35_39m$estimate[1], fit2_alt_35_39m$estimate[2])  #Anderson-Darling test

# gamma has much better fit but lognormal still above 0.05

#Age 40-44, male
Proc40_44m <- subset(Alt_scenario, agegroups=="40-44" & sex =="1", select = procmeat.new)

#Probability of intake
probproc_alt40_44m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc40_44m==0)/length(Proc40_44m$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc40_44m!=0)/length(Proc40_44m$procmeat.new) #probability of consumption
probproc_alt40_44m[1,] <- c(pr.no,pr.yes)

Proc40_44posm <- Proc40_44m[which(Proc40_44m!=0)]
Proc40_44posm <- as.vector(Proc40_44posm$procmeat.new)


fit1_alt_40_44m <- fitdist(Proc40_44posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc40_44posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_40_44m$estimate[1], sdlog=fit1_alt_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_40_44m, fitnames="lnorm") 

fit2_alt_40_44m <- fitdist(Proc40_44posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc40_44posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_alt_40_44m$estimate[1], fit2_alt_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_40_44m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_40_44m$estimate[1], fit1_alt_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_40_44m$estimate[1], fit1_alt_40_44m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_40_44m$estimate[1], fit2_alt_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_40_44m$estimate[1], fit2_alt_40_44m$estimate[2])  #Anderson-Darling test

# gamma, lognormal above 0.05

#Age 45-49, male
Proc45_49m <- subset(Alt_scenario, agegroups=="45-49" & sex =="1", select = procmeat.new)

#Probability of intake
probproc_alt45_49m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc45_49m==0)/length(Proc45_49m$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc45_49m!=0)/length(Proc45_49m$procmeat.new) #probability of consumption
probproc_alt45_49m[1,] <- c(pr.no,pr.yes)

Proc45_49posm <- Proc45_49m[which(Proc45_49m!=0)]
Proc45_49posm <- as.vector(Proc45_49posm$procmeat.new)


fit1_alt_45_49m <- fitdist(Proc45_49posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc45_49posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_45_49m$estimate[1], sdlog=fit1_alt_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_45_49m, fitnames="lnorm") 

fit2_alt_45_49m <- fitdist(Proc45_49posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc45_49posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_alt_45_49m$estimate[1], fit2_alt_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_45_49m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_45_49m$estimate[1], fit1_alt_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_45_49m$estimate[1], fit1_alt_45_49m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_45_49m$estimate[1], fit2_alt_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_45_49m$estimate[1], fit2_alt_45_49m$estimate[2])  #Anderson-Darling test

# gamma has much better fit but lognormal above 0.05

#Age 50-54, male
Proc50_54m <- subset(Alt_scenario, agegroups=="50-54" & sex =="1", select = procmeat.new)

#Probability of intake
probproc_alt50_54m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc50_54m==0)/length(Proc50_54m$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc50_54m!=0)/length(Proc50_54m$procmeat.new) #probability of consumption
probproc_alt50_54m[1,] <- c(pr.no,pr.yes)

Proc50_54posm <- Proc50_54m[which(Proc50_54m!=0)]
Proc50_54posm <- as.vector(Proc50_54posm$procmeat.new)


fit1_alt_50_54m <- fitdist(Proc50_54posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc50_54posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_50_54m$estimate[1], sdlog=fit1_alt_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_50_54m, fitnames="lnorm") 

fit2_alt_50_54m <- fitdist(Proc50_54posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc50_54posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_alt_50_54m$estimate[1], fit2_alt_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_50_54m, fitnames="gamma") 


cvm.test(t, plnorm, fit1_alt_50_54m$estimate[1], fit1_alt_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_50_54m$estimate[1], fit1_alt_50_54m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_50_54m$estimate[1], fit2_alt_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_50_54m$estimate[1], fit2_alt_50_54m$estimate[2])  #Anderson-Darling test

# gamma has the best fit, lognormal around 0.05 but above 0.01

#Age 55-59, male
Proc55_59m <- subset(Alt_scenario, agegroups=="55-59" & sex =="1", select = procmeat.new)

#Probability of intake
probproc_alt55_59m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc55_59m==0)/length(Proc55_59m$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc55_59m!=0)/length(Proc55_59m$procmeat.new) #probability of consumption
probproc_alt55_59m[1,] <- c(pr.no,pr.yes)

Proc55_59posm <- Proc55_59m[which(Proc55_59m!=0)]
Proc55_59posm <- as.vector(Proc55_59posm$procmeat.new)


fit1_alt_55_59m <- fitdist(Proc55_59posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc55_59posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_55_59m$estimate[1], sdlog=fit1_alt_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_55_59m, fitnames="lnorm") 

fit2_alt_55_59m <- fitdist(Proc55_59posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc55_59posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_alt_55_59m$estimate[1], fit2_alt_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_55_59m, fitnames="gamma") 


cvm.test(t, plnorm, fit1_alt_55_59m$estimate[1], fit1_alt_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_55_59m$estimate[1], fit1_alt_55_59m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_55_59m$estimate[1], fit2_alt_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_55_59m$estimate[1], fit2_alt_55_59m$estimate[2])  #Anderson-Darling test

# gamma has the best fit, lognormal above 0.05

#Age 60-64, male
Proc60_64m <- subset(Alt_scenario, agegroups=="60-64" & sex =="1", select = procmeat.new)

#Probability of intake
probproc_alt60_64m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc60_64m==0)/length(Proc60_64m$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc60_64m!=0)/length(Proc60_64m$procmeat.new) #probability of consumption
probproc_alt60_64m[1,] <- c(pr.no,pr.yes)

Proc60_64posm <- Proc60_64m[which(Proc60_64m!=0)]
Proc60_64posm <- as.vector(Proc60_64posm$procmeat.new)


fit1_alt_60_64m <- fitdist(Proc60_64posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc60_64posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_60_64m$estimate[1], sdlog=fit1_alt_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_60_64m, fitnames="lnorm") 

fit2_alt_60_64m <- fitdist(Proc60_64posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc60_64posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_alt_60_64m$estimate[1], fit2_alt_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_60_64m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_60_64m$estimate[1], fit1_alt_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_60_64m$estimate[1], fit1_alt_60_64m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_60_64m$estimate[1], fit2_alt_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_60_64m$estimate[1], fit2_alt_60_64m$estimate[2])  #Anderson-Darling test

# gamma has the best fit but lognormal above 0.05

#Age 65-69, male
Proc65_69m <- subset(Alt_scenario, agegroups=="65-69" & sex =="1", select = procmeat.new)

#Probability of intake
probproc_alt65_69m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc65_69m==0)/length(Proc65_69m$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc65_69m!=0)/length(Proc65_69m$procmeat.new) #probability of consumption
probproc_alt65_69m[1,] <- c(pr.no,pr.yes)

Proc65_69posm <- Proc65_69m[which(Proc65_69m!=0)]
Proc65_69posm <- as.vector(Proc65_69posm$procmeat.new)

fit1_alt_65_69m <- fitdist(Proc65_69posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc65_69posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_65_69m$estimate[1], sdlog=fit1_alt_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_65_69m, fitnames="lnorm") 

fit2_alt_65_69m <- fitdist(Proc65_69posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc65_69posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_alt_65_69m$estimate[1], fit2_alt_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_65_69m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_65_69m$estimate[1], fit1_alt_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_65_69m$estimate[1], fit1_alt_65_69m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_65_69m$estimate[1], fit2_alt_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_65_69m$estimate[1], fit2_alt_65_69m$estimate[2])  #Anderson-Darling test

# gamma has a much better fit but lognormal above 0.05

#Age 70-74, male
Proc70_74m <- subset(Alt_scenario, agegroups=="70-74" & sex =="1", select = procmeat.new)

#Probability of intake
probproc_alt70_74m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc70_74m==0)/length(Proc70_74m$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc70_74m!=0)/length(Proc70_74m$procmeat.new) #probability of consumption
probproc_alt70_74m[1,] <- c(pr.no,pr.yes)

Proc70_74posm <- Proc70_74m[which(Proc70_74m!=0)]
Proc70_74posm <- as.vector(Proc70_74posm$procmeat.new)

fit1_alt_70_74m <- fitdist(Proc70_74posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc70_74posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_70_74m$estimate[1], sdlog=fit1_alt_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_70_74m, fitnames="lnorm") 

fit2_alt_70_74m <- fitdist(Proc70_74posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc70_74posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_alt_70_74m$estimate[1], fit2_alt_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_70_74m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_70_74m$estimate[1], fit1_alt_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_70_74m$estimate[1], fit1_alt_70_74m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_70_74m$estimate[1], fit2_alt_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_70_74m$estimate[1], fit2_alt_70_74m$estimate[2])  #Anderson-Darling test

# gamma has the best fit, lognormal below 0.05 but above 0.01

#Age 75-79, male
Proc75_79m <- subset(Alt_scenario, agegroups=="75-79" & sex =="1", select = procmeat.new)

#Probability of intake
probproc_alt75_79m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc75_79m==0)/length(Proc75_79m$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc75_79m!=0)/length(Proc75_79m$procmeat.new) #probability of consumption
probproc_alt75_79m[1,] <- c(pr.no,pr.yes)

Proc75_79posm <- Proc75_79m[which(Proc75_79m!=0)]
Proc75_79posm <- as.vector(Proc75_79posm$procmeat.new)

fit1_alt_75_79m <- fitdist(Proc75_79posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc75_79posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_75_79m$estimate[1], sdlog=fit1_alt_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_75_79m, fitnames="lnorm") 

fit2_alt_75_79m <- fitdist(Proc75_79posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc75_79posm
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=40)
lines(x, pgamma(x, fit2_alt_75_79m$estimate[1], fit2_alt_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_75_79m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_75_79m$estimate[1], fit1_alt_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_75_79m$estimate[1], fit1_alt_75_79m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_75_79m$estimate[1], fit2_alt_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_75_79m$estimate[1], fit2_alt_75_79m$estimate[2])  #Anderson-Darling test

# gamma and lognormal equally good


##### Fitting exposures females ######


#risk of SC = 0 for women 15-19 years

#Age 20-24, female
Proc20_24w <- subset(Alt_scenario, agegroups=="20-24" & sex =="2", select = procmeat.new)

#Probability of intake
probproc_alt20_24w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc20_24w==0)/length(Proc20_24w$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc20_24w!=0)/length(Proc20_24w$procmeat.new) #probability of consumption
probproc_alt20_24w[1,] <- c(pr.no,pr.yes)

Proc20_24posw <- Proc20_24w[which(Proc20_24w!=0)]
Proc20_24posw <- as.vector(Proc20_24posw$procmeat.new)

fit1_alt_20_24w <- fitdist(Proc20_24posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc20_24posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_20_24w$estimate[1], sdlog=fit1_alt_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_20_24w, fitnames="lnorm") 

fit2_alt_20_24w <- fitdist(Proc20_24posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc20_24posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_alt_20_24w$estimate[1], fit2_alt_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_20_24w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_20_24w$estimate[1], fit1_alt_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_20_24w$estimate[1], fit1_alt_20_24w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_20_24w$estimate[1], fit2_alt_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_20_24w$estimate[1], fit2_alt_20_24w$estimate[2])  #Anderson-Darling test

#lognormal very good but gamma also above 0.05

#Age 25-29, female
Proc25_29w <- subset(Alt_scenario, agegroups=="25-29" & sex =="2", select = procmeat.new)

#Probability of intake
probproc_alt25_29w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc25_29w==0)/length(Proc25_29w$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc25_29w!=0)/length(Proc25_29w$procmeat.new) #probability of consumption
probproc_alt25_29w[1,] <- c(pr.no,pr.yes)

Proc25_29posw <- Proc25_29w[which(Proc25_29w!=0)]
Proc25_29posw <- as.vector(Proc25_29posw$procmeat.new)

fit1_alt_25_29w <- fitdist(Proc25_29posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc25_29posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_25_29w$estimate[1], sdlog=fit1_alt_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_25_29w, fitnames="lnorm")

fit2_alt_25_29w <- fitdist(Proc25_29posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc25_29posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_alt_25_29w$estimate[1], fit2_alt_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_25_29w, fitnames="gamma")

cvm.test(t, plnorm, fit1_alt_25_29w$estimate[1], fit1_alt_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_25_29w$estimate[1], fit1_alt_25_29w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_25_29w$estimate[1], fit2_alt_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_25_29w$estimate[1], fit2_alt_25_29w$estimate[2])  #Anderson-Darling test

# gamma has the best fit but lognormal above 0.05


#risk of SC = 0 for women 30-39 years


#Age 40-44, female
Proc40_44w <- subset(Alt_scenario, agegroups=="40-44" & sex =="2", select = procmeat.new)

#Probability of intake
probproc_alt40_44w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc40_44w==0)/length(Proc40_44w$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc40_44w!=0)/length(Proc40_44w$procmeat.new) #probability of consumption
probproc_alt40_44w[1,] <- c(pr.no,pr.yes)

Proc40_44posw <- Proc40_44w[which(Proc40_44w!=0)]
Proc40_44posw <- as.vector(Proc40_44posw$procmeat.new)


fit1_alt_40_44w <- fitdist(Proc40_44posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc40_44posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_40_44w$estimate[1], sdlog=fit1_alt_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_40_44w, fitnames="lnorm") 

fit2_alt_40_44w <- fitdist(Proc40_44posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc40_44posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_alt_40_44w$estimate[1], fit2_alt_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_40_44w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_40_44w$estimate[1], fit1_alt_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_40_44w$estimate[1], fit1_alt_40_44w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_40_44w$estimate[1], fit2_alt_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_40_44w$estimate[1], fit2_alt_40_44w$estimate[2])  #Anderson-Darling test

# gamma, lognormal above 0.05

#Age 45-49, female
Proc45_49w <- subset(Alt_scenario, agegroups=="45-49" & sex =="2", select = procmeat.new)

#Probability of intake
probproc_alt45_49w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc45_49w==0)/length(Proc45_49w$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc45_49w!=0)/length(Proc45_49w$procmeat.new) #probability of consumption
probproc_alt45_49w[1,] <- c(pr.no,pr.yes)

Proc45_49posw <- Proc45_49w[which(Proc45_49w!=0)]
Proc45_49posw <- as.vector(Proc45_49posw$procmeat.new)


fit1_alt_45_49w <- fitdist(Proc45_49posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc45_49posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_45_49w$estimate[1], sdlog=fit1_alt_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_45_49w, fitnames="lnorm") 

fit2_alt_45_49w <- fitdist(Proc45_49posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc45_49posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_alt_45_49w$estimate[1], fit2_alt_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_45_49w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_45_49w$estimate[1], fit1_alt_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_45_49w$estimate[1], fit1_alt_45_49w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_45_49w$estimate[1], fit2_alt_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_45_49w$estimate[1], fit2_alt_45_49w$estimate[2])  #Anderson-Darling test

# gamma has much better fit but lognormal just above 0.05

#Age 50-54, female
Proc50_54w <- subset(Alt_scenario, agegroups=="50-54" & sex =="2", select = procmeat.new)

#Probability of intake
probproc_alt50_54w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc50_54w==0)/length(Proc50_54w$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc50_54w!=0)/length(Proc50_54w$procmeat.new) #probability of consumption
probproc_alt50_54w[1,] <- c(pr.no,pr.yes)

Proc50_54posw <- Proc50_54w[which(Proc50_54w!=0)]
Proc50_54posw <- as.vector(Proc50_54posw$procmeat.new)


fit1_alt_50_54w <- fitdist(Proc50_54posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc50_54posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_50_54w$estimate[1], sdlog=fit1_alt_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_50_54w, fitnames="lnorm") 

fit2_alt_50_54w <- fitdist(Proc50_54posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc50_54posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_alt_50_54w$estimate[1], fit2_alt_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_50_54w, fitnames="gamma") 


cvm.test(t, plnorm, fit1_alt_50_54w$estimate[1], fit1_alt_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_50_54w$estimate[1], fit1_alt_50_54w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_50_54w$estimate[1], fit2_alt_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_50_54w$estimate[1], fit2_alt_50_54w$estimate[2])  #Anderson-Darling test

# gamma has the best fit, lognormal above 0.05

#Age 55-59, female
Proc55_59w <- subset(Alt_scenario, agegroups=="55-59" & sex =="2", select = procmeat.new)

#Probability of intake
probproc_alt55_59w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc55_59w==0)/length(Proc55_59w$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc55_59w!=0)/length(Proc55_59w$procmeat.new) #probability of consumption
probproc_alt55_59w[1,] <- c(pr.no,pr.yes)

Proc55_59posw <- Proc55_59w[which(Proc55_59w!=0)]
Proc55_59posw <- as.vector(Proc55_59posw$procmeat.new)


fit1_alt_55_59w <- fitdist(Proc55_59posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc55_59posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_55_59w$estimate[1], sdlog=fit1_alt_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_55_59w, fitnames="lnorm") 

fit2_alt_55_59w <- fitdist(Proc55_59posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc55_59posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_alt_55_59w$estimate[1], fit2_alt_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_55_59w, fitnames="gamma") 


cvm.test(t, plnorm, fit1_alt_55_59w$estimate[1], fit1_alt_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_55_59w$estimate[1], fit1_alt_55_59w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_55_59w$estimate[1], fit2_alt_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_55_59w$estimate[1], fit2_alt_55_59w$estimate[2])  #Anderson-Darling test

# gamma has the best fit, lognormal above 0.05

#Age 60-64, female
Proc60_64w <- subset(Alt_scenario, agegroups=="60-64" & sex =="2", select = procmeat.new)

#Probability of intake
probproc_alt60_64w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc60_64w==0)/length(Proc60_64w$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc60_64w!=0)/length(Proc60_64w$procmeat.new) #probability of consumption
probproc_alt60_64w[1,] <- c(pr.no,pr.yes)

Proc60_64posw <- Proc60_64w[which(Proc60_64w!=0)]
Proc60_64posw <- as.vector(Proc60_64posw$procmeat.new)


fit1_alt_60_64w <- fitdist(Proc60_64posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc60_64posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_60_64w$estimate[1], sdlog=fit1_alt_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_60_64w, fitnames="lnorm") 

fit2_alt_60_64w <- fitdist(Proc60_64posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc60_64posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_alt_60_64w$estimate[1], fit2_alt_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_60_64w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_60_64w$estimate[1], fit1_alt_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_60_64w$estimate[1], fit1_alt_60_64w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_60_64w$estimate[1], fit2_alt_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_60_64w$estimate[1], fit2_alt_60_64w$estimate[2])  #Anderson-Darling test

# gamma has the best fit but lognormal  just above 0.05

#Age 65-69, female
Proc65_69w <- subset(Alt_scenario, agegroups=="65-69" & sex =="2", select = procmeat.new)

#Probability of intake
probproc_alt65_69w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc65_69w==0)/length(Proc65_69w$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc65_69w!=0)/length(Proc65_69w$procmeat.new) #probability of consumption
probproc_alt65_69w[1,] <- c(pr.no,pr.yes)

Proc65_69posw <- Proc65_69w[which(Proc65_69w!=0)]
Proc65_69posw <- as.vector(Proc65_69posw$procmeat.new)

fit1_alt_65_69w <- fitdist(Proc65_69posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc65_69posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_65_69w$estimate[1], sdlog=fit1_alt_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_65_69w, fitnames="lnorm") 

fit2_alt_65_69w <- fitdist(Proc65_69posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc65_69posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_alt_65_69w$estimate[1], fit2_alt_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_65_69w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_65_69w$estimate[1], fit1_alt_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_65_69w$estimate[1], fit1_alt_65_69w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_65_69w$estimate[1], fit2_alt_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_65_69w$estimate[1], fit2_alt_65_69w$estimate[2])  #Anderson-Darling test

# gamma and lognormal are equally good

#Age 70-74, female
Proc70_74w <- subset(Alt_scenario, agegroups=="70-74" & sex =="2", select = procmeat.new)

#Probability of intake
probproc_alt70_74w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc70_74w==0)/length(Proc70_74w$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc70_74w!=0)/length(Proc70_74w$procmeat.new) #probability of consumption
probproc_alt70_74w[1,] <- c(pr.no,pr.yes)

Proc70_74posw <- Proc70_74w[which(Proc70_74w!=0)]
Proc70_74posw <- as.vector(Proc70_74posw$procmeat.new)

fit1_alt_70_74w <- fitdist(Proc70_74posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc70_74posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_70_74w$estimate[1], sdlog=fit1_alt_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_70_74w, fitnames="lnorm") 

fit2_alt_70_74w <- fitdist(Proc70_74posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc70_74posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_alt_70_74w$estimate[1], fit2_alt_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_70_74w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_70_74w$estimate[1], fit1_alt_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_70_74w$estimate[1], fit1_alt_70_74w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_70_74w$estimate[1], fit2_alt_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_70_74w$estimate[1], fit2_alt_70_74w$estimate[2])  #Anderson-Darling test

# gamma has the best fit but lognormal above 0.05

#Age 75-79, female
Proc75_79w <- subset(Alt_scenario, agegroups=="75-79" & sex =="2", select = procmeat.new)

#Probability of intake
probproc_alt75_79w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Proc75_79w==0)/length(Proc75_79w$procmeat.new) #probability of zero consumption
pr.yes <- sum(Proc75_79w!=0)/length(Proc75_79w$procmeat.new) #probability of consumption
probproc_alt75_79w[1,] <- c(pr.no,pr.yes)

Proc75_79posw <- Proc75_79w[which(Proc75_79w!=0)]
Proc75_79posw <- as.vector(Proc75_79posw$procmeat.new)

fit1_alt_75_79w <- fitdist(Proc75_79posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Proc75_79posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_75_79w$estimate[1], sdlog=fit1_alt_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_75_79w, fitnames="lnorm") 

fit2_alt_75_79w <- fitdist(Proc75_79posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Proc75_79posw
plot(ecdf(t2), lty=1)
x <- seq(0, 2000, length=90)
lines(x, pgamma(x, fit2_alt_75_79w$estimate[1], fit2_alt_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_75_79w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_75_79w$estimate[1], fit1_alt_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_75_79w$estimate[1], fit1_alt_75_79w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_75_79w$estimate[1], fit2_alt_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_75_79w$estimate[1], fit2_alt_75_79w$estimate[2])  #Anderson-Darling test

# gamma and lognormal are equally good


##### DALY calc males #####

max(MeatDaily$proc.meat)
# [1] 467.1429

# truncate at 500 g processed/day. Truncation did not change mean.


# risk of SC = 0 for 15-29 year old males


# Age 30-34, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc30_34m[,1],probproc30_34m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_30_34m$estimate[1], fit2_30_34m$estimate[2], rtrunc = T, linf = 0, lsup = 500)
r <- mcstoc(rpert, type = "U", min = 0.0002, mode = 0.0036, max = 0.0076) #1.18, 95% CI: 1.01-1.38 converted to RR/g/day
RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt30_34m[,1],probproc_alt30_34m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_30_34m$estimate[1], fit2_alt_30_34m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p30_34m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_30_34m <- pEffect * pop_30_34m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal30_34m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal30_34m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal30_34m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal30_34m * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal30_34m * (SEYLL_30_34 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY30_34totalm_alt <- DALY * pop_30_34m
summary(DALY30_34totalm_alt)



summary(DALY30_34totalm)



# Age 35-39, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc35_39m[,1],probproc35_39m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_35_39m$estimate[1], fit2_35_39m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt35_39m[,1],probproc_alt35_39m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_35_39m$estimate[1], fit2_alt_35_39m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p35_39m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_35_39m <- pEffect * pop_35_39m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal35_39m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal35_39m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal35_39m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal35_39m * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal35_39m * (SEYLL_35_39 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY35_39totalm_alt <- DALY * pop_35_39m
summary(DALY35_39totalm_alt)


summary(DALY35_39totalm)



# Age 40-44, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc40_44m[,1],probproc40_44m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_40_44m$estimate[1], fit2_40_44m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt40_44m[,1],probproc_alt40_44m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_40_44m$estimate[1], fit2_alt_40_44m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p40_44m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_40_44m <- pEffect * pop_40_44m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal40_44m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal40_44m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal40_44m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal40_44m * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal40_44m * (SEYLL_40_44 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY40_44totalm_alt <- DALY * pop_40_44m
summary(DALY40_44totalm_alt)


summary(DALY40_44totalm)



# Age 45-49, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc45_49m[,1],probproc45_49m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_45_49m$estimate[1], fit2_45_49m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt45_49m[,1],probproc_alt45_49m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_45_49m$estimate[1], fit2_alt_45_49m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p45_49m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_45_49m <- pEffect * pop_45_49m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal45_49m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal45_49m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal45_49m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal45_49m * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal45_49m * (SEYLL_45_49 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY45_49totalm_alt <- DALY * pop_45_49m
summary(DALY45_49totalm_alt)


summary(DALY45_49totalm)



# Age 50-54, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc50_54m[,1],probproc50_54m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_50_54m$estimate[1], fit2_50_54m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt50_54m[,1],probproc_alt50_54m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_50_54m$estimate[1], fit2_alt_50_54m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p50_54m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_50_54m <- pEffect * pop_50_54m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal50_54m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal50_54m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal50_54m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal50_54m * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal50_54m * (SEYLL_50_54 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY50_54totalm_alt <- DALY * pop_50_54m
summary(DALY50_54totalm_alt)


summary(DALY50_54totalm)


# Age 55-59, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc55_59m[,1],probproc55_59m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_55_59m$estimate[1], fit2_55_59m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt55_59m[,1],probproc_alt55_59m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_55_59m$estimate[1], fit2_alt_55_59m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p55_59m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_55_59m <- pEffect * pop_55_59m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal55_59m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal55_59m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal55_59m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal55_59m * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal55_59m * (SEYLL_55_59 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY55_59totalm_alt <- DALY * pop_55_59m
summary(DALY55_59totalm_alt)


summary(DALY55_59totalm)


# Age 60-64, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc60_64m[,1],probproc60_64m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_60_64m$estimate[1], fit2_60_64m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt60_64m[,1],probproc_alt60_64m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_60_64m$estimate[1], fit2_alt_60_64m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p60_64m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_60_64m <- pEffect * pop_60_64m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal60_64m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal60_64m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal60_64m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal60_64m * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal60_64m * (SEYLL_60_64 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY60_64totalm_alt <- DALY * pop_60_64m
summary(DALY60_64totalm_alt)


summary(DALY60_64totalm)


# Age 65-69, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc65_69m[,1],probproc65_69m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_65_69m$estimate[1], fit2_65_69m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt65_69m[,1],probproc_alt65_69m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_65_69m$estimate[1], fit2_alt_65_69m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p65_69m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_65_69m <- pEffect * pop_65_69m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal65_69m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal65_69m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal65_69m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal65_69m * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal65_69m * (SEYLL_65_69 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY65_69totalm_alt <- DALY * pop_65_69m
summary(DALY65_69totalm_alt)


summary(DALY65_69totalm)



# Age 70-74, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc70_74m[,1],probproc70_74m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_70_74m$estimate[1], fit2_70_74m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt70_74m[,1],probproc_alt70_74m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_70_74m$estimate[1], fit2_alt_70_74m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p70_74m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_70_74m <- pEffect * pop_70_74m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal70_74m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal70_74m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal70_74m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal70_74m * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal70_74m * (SEYLL_70_74 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY70_74totalm_alt <- DALY * pop_70_74m
summary(DALY70_74totalm_alt)


summary(DALY70_74totalm)


# Age 75-79, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc75_79m[,1],probproc75_79m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_75_79m$estimate[1], fit2_75_79m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt75_79m[,1],probproc_alt75_79m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_75_79m$estimate[1], fit2_alt_75_79m$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p75_79m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_75_79m <- pEffect * pop_75_79m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal75_79m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal75_79m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal75_79m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal75_79m * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal75_79m * (SEYLL_75_79 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY75_79totalm_alt <- DALY * pop_75_79m
summary(DALY75_79totalm_alt)


summary(DALY75_79totalm)



# Age 80-84, males
# Extrapolate RR from 75 year olds to >80 years

pEffect <- p80_84m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_80_84m <- pEffect * pop_80_84m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal80_84m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal80_84m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal80_84m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal80_84m * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal80_84m * (SEYLL_80_84 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY80_84totalm_alt <- DALY * pop_80_84m
summary(DALY80_84totalm_alt)


summary(DALY80_84totalm)




# Age 85+, males
# Extrapolate RR from 75 year olds to >80 years

pEffect <- p85m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_85m <- pEffect * pop_85m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal85m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal85m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal85m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal85m * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal85m * (SEYLL_85 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY85totalm_alt <- DALY * pop_85m
summary(DALY85totalm_alt)


summary(DALY85totalm)


##### DALY calc females #####

# risk of SC = 0 for 15-19 and 30-39 yearold females


# Age 20-24, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc20_24w[,1],probproc20_24w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_20_24w$estimate[1], fit2_20_24w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt20_24w[,1],probproc_alt20_24w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_20_24w$estimate[1], fit2_alt_20_24w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)



pEffect <- p20_24w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_20_24w <- pEffect * pop_20_24w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal20_24w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal20_24w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal20_24w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal20_24w * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal20_24w * (SEYLL_20_24 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY20_24totalw_alt <- DALY * pop_20_24w
summary(DALY20_24totalw_alt)


summary(DALY20_24totalw)


# Age 25-29, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc25_29w[,1],probproc25_29w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_25_29w$estimate[1], fit2_25_29w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt25_29w[,1],probproc_alt25_29w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_25_29w$estimate[1], fit2_alt_25_29w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p25_29w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_25_29w <- pEffect * pop_25_29w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal25_29w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal25_29w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal25_29w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal25_29w * dw_terminal * time_terminal #terminal
  
  
#YLLs
YLL <- pEffect * fatal25_29w * (SEYLL_25_29 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY25_29totalw_alt <- DALY * pop_25_29w
summary(DALY25_29totalw_alt)


summary(DALY25_29totalw)



# Age 30-34, females, risk of SC = 0

# Age 35-39, females, risk of SC = 0


# Age 40-44, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc40_44w[,1],probproc40_44w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_40_44w$estimate[1], fit2_40_44w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt40_44w[,1],probproc_alt40_44w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_40_44w$estimate[1], fit2_alt_40_44w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p40_44w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_40_44w <- pEffect * pop_40_44w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal40_44w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal40_44w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal40_44w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal40_44w * dw_terminal * time_terminal #terminal
 

#YLLs
YLL <- pEffect * fatal40_44w * (SEYLL_40_44 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY40_44totalw_alt <- DALY * pop_40_44w
summary(DALY40_44totalw_alt)


summary(DALY40_44totalw)



# Age 45-49, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc45_49w[,1],probproc45_49w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_45_49w$estimate[1], fit2_45_49w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt45_49w[,1],probproc_alt45_49w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_45_49w$estimate[1], fit2_alt_45_49w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p45_49w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_45_49w <- pEffect * pop_45_49w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal45_49w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal45_49w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal45_49w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal45_49w * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal45_49w * (SEYLL_45_49 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY45_49totalw_alt <- DALY * pop_45_49w
summary(DALY45_49totalw_alt)



summary(DALY45_49totalw)




# Age 50-54, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc50_54w[,1],probproc50_54w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_50_54w$estimate[1], fit2_50_54w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt50_54w[,1],probproc_alt50_54w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_50_54w$estimate[1], fit2_alt_50_54w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p50_54w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_50_54w <- pEffect * pop_50_54w


YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal50_54w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal50_54w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal50_54w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal50_54w * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal50_54w * (SEYLL_50_54 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY50_54totalw_alt <- DALY * pop_50_54w
summary(DALY50_54totalw_alt)


summary(DALY50_54totalw)



# Age 55-59, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc55_59w[,1],probproc55_59w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_55_59w$estimate[1], fit2_55_59w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt55_59w[,1],probproc_alt55_59w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_55_59w$estimate[1], fit2_alt_55_59w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p55_59w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_55_59w <- pEffect * pop_55_59w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal55_59w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal55_59w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal55_59w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal55_59w * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal55_59w * (SEYLL_55_59 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY55_59totalw_alt <- DALY * pop_55_59w
summary(DALY55_59totalw_alt)


summary(DALY55_59totalw)



# Age 60-64, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc60_64w[,1],probproc60_64w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_60_64w$estimate[1], fit2_60_64w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt60_64w[,1],probproc_alt60_64w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_60_64w$estimate[1], fit2_alt_60_64w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p60_64w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_60_64w <- pEffect * pop_60_64w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal60_64w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal60_64w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal60_64w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal60_64w * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal60_64w * (SEYLL_60_64 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY60_64totalw_alt <- DALY * pop_60_64w
summary(DALY60_64totalw_alt)


summary(DALY60_64totalw)


# Age 65-69, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc65_69w[,1],probproc65_69w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_65_69w$estimate[1], fit2_65_69w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt65_69w[,1],probproc_alt65_69w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_65_69w$estimate[1], fit2_alt_65_69w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p65_69w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_65_69w <- pEffect * pop_65_69w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal65_69w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal65_69w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal65_69w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal65_69w * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal65_69w * (SEYLL_65_69 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY65_69totalw_alt <- DALY * pop_65_69w
summary(DALY65_69totalw_alt)


summary(DALY65_69totalw)




# Age 70-74, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc70_74w[,1],probproc70_74w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_70_74w$estimate[1], fit2_70_74w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt70_74w[,1],probproc_alt70_74w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_70_74w$estimate[1], fit2_alt_70_74w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p70_74w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_70_74w <- pEffect * pop_70_74w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal70_74w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal70_74w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal70_74w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal70_74w * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal70_74w * (SEYLL_70_74 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY70_74totalw_alt <- DALY * pop_70_74w
summary(DALY70_74totalw_alt)


summary(DALY70_74totalw)



# Age 75-79, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc75_79w[,1],probproc75_79w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_75_79w$estimate[1], fit2_75_79w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probproc_alt75_79w[,1],probproc_alt75_79w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_75_79w$estimate[1], fit2_alt_75_79w$estimate[2], rtrunc = T, linf = 0, lsup = 500)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p75_79w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_75_79w <- pEffect * pop_75_79w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal75_79w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal75_79w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal75_79w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal75_79w * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal75_79w * (SEYLL_75_79 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY75_79totalw_alt <- DALY * pop_75_79w
summary(DALY75_79totalw_alt)


summary(DALY75_79totalw)



# Age 80-84, females

#Extrapolate intakes for 75 year-olds to >80

pEffect <- p80_84w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_80_84w <- pEffect * pop_80_84w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal80_84w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal80_84w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal80_84w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal80_84w * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal80_84w * (SEYLL_80_84 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY80_84totalw_alt <- DALY * pop_80_84w
summary(DALY80_84totalw_alt)


summary(DALY80_84totalw)



# Age 85+, females

#Extrapolate intakes for 75 year-olds to >80

pEffect <- p85w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_85w <- pEffect * pop_85w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal85w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal85w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal85w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal85w * dw_terminal * time_terminal #terminal


#YLLs
YLL <- pEffect * fatal85w * (SEYLL_85 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY85totalw_alt <- DALY * pop_85w
summary(DALY85totalw_alt)


summary(DALY85totalw)


##### Total DALYs Alt #####

DALYtotal_alt <- DALY30_34totalm_alt + DALY35_39totalm_alt +  DALY40_44totalm_alt + DALY45_49totalm_alt +
  DALY50_54totalm_alt + DALY55_59totalm_alt + DALY60_64totalm_alt + DALY65_69totalm_alt +  DALY70_74totalm_alt +
  DALY75_79totalm_alt + DALY80_84totalm_alt + DALY85totalm_alt + 
  DALY20_24totalw_alt + DALY25_29totalw_alt +  DALY40_44totalw_alt + DALY45_49totalw_alt + DALY50_54totalw_alt +
  DALY55_59totalw_alt + DALY60_64totalw_alt + DALY65_69totalw_alt +  DALY70_74totalw_alt + DALY75_79totalw_alt +
  DALY80_84totalw_alt + DALY85totalw_alt

summary(DALYtotal_alt)

##### Total cases alternative scenario #####

cases_alt <- cases_30_34m + cases_35_39m + cases_40_44m + cases_45_49m + cases_50_54m + cases_55_59m + cases_60_64m +
  cases_65_69m + cases_70_74m + cases_75_79m + cases_80_84m + cases_85m +
  cases_20_24w + cases_25_29w + cases_40_44w + cases_45_49w + cases_50_54w + cases_55_59w + cases_60_64w +
  cases_65_69w + cases_70_74w + cases_75_79w + cases_80_84w + cases_85w

summary(cases_alt)


##### Extra number of cases #####

cases_diff <- cases_alt - cases_ref

summary(cases_diff)



##### Per 100,000 #####

pop15_85_plus <- pop_men + pop_women


tDALY_ref_proc_SC <- DALYtotalref
summary(tDALY_ref_proc_SC)


tDALY_alt_proc_SC <-DALYtotal_alt
summary(tDALY_alt_proc_SC)


tDALY_ref_proc_SC_100000 <- tDALY_ref_proc_SC/pop15_85_plus*1e+05
summary(tDALY_ref_proc_SC_100000)


tDALY_alt_proc_SC_100000 <- tDALY_alt_proc_SC/pop15_85_plus*1e+05
summary(tDALY_alt_proc_SC_100000)



dDALY_alt_proc_SC <- tDALY_alt_proc_SC - tDALY_ref_proc_SC
summary(dDALY_alt_proc_SC)


dDALY_alt_proc_SC_100000 <- tDALY_alt_proc_SC_100000 - tDALY_ref_proc_SC_100000
summary(dDALY_alt_proc_SC_100000)


##### Save #####

tDALY_ref_proc_SC_unc <- apply(tDALY_ref_proc_SC, 2, function(x) x*1) #get data out of mc2d
write.csv(tDALY_ref_proc_SC_unc, "tDALY_ref_proc_SC_unc.csv")


tDALY_alt_proc_SC_matrix <- apply(tDALY_alt_proc_SC, 2, function(x) x*1) #get data out of mc2d
tDALY_alt_proc_SC_unc <- apply(tDALY_alt_proc_SC_matrix, 2, function(x) mean(x)) #only uncertainty dimension
write.csv(tDALY_alt_proc_SC_unc, "tDALY_alt_proc_SC_unc.csv")
