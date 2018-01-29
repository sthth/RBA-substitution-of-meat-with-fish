memory.limit(200000)

#packages
library(data.table)
library(mc2d)
library(fitdistrplus)
library(goftest)

#settings
iters_var <- 1e+05
iters_unc <- 1e+03

ndvar(iters_var)
ndunc(iters_unc)

#Population statistics
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


#CRC cancer incidence/mortality in Denmark, 2015

## Background risk of CRC for different age and sex groups

p15_19m <- 0.00003
p20_24m <- 0.00004
p25_29m <- 0.00002
p30_34m <- 0.00008
p35_39m <- 0.00008
p40_44m <- 0.00019
p45_49m <- 0.00034
p50_54m <- 0.00088
p55_59m <- 0.0013
p60_64m <- 0.00239
p65_69m <- 0.00409
p70_74m <- 0.0056
p75_79m <- 0.00642
p80_84m <- 0.00731
p85m <- 0.00683

p15_19w <- 0.00002
p20_24w <- 0.00002
p25_29w <- 0.00001
p30_34w <- 0.00006
p35_39w <- 0.00008
p40_44w <- 0.00014
p45_49w <- 0.00026
p50_54w <- 0.00066
p55_59w <- 0.00103
p60_64w <- 0.00154
p65_69w <- 0.00233
p70_74w <- 0.0032
p75_79w <- 0.00381
p80_84w <- 0.00435
p85w <- 0.00398

## Case fatality for CRC for different age and sex groups

fatal15_19m <- 0
fatal20_24m <- 0
fatal25_29m <- 0
fatal30_34m <- 0.075
fatal35_39m <- 0.075
fatal40_44m <- 0.052631579
fatal45_49m <- 0.285294118
fatal50_54m <- 0.152272727
fatal55_59m <- 0.19
fatal60_64m <- 0.238912134
fatal65_69m <- 0.164303178
fatal70_74m <- 0.241964286
fatal75_79m <- 0.269314642
fatal80_84m <- 0.379069767
fatal85m <- 0.59795022

fatal15_19w <- 0
fatal20_24w <- 0
fatal25_29w <- 0
fatal30_34w <- 0.516666667
fatal35_39w <- 0.2125
fatal40_44w <- 0.142857143
fatal45_49w <- 0.323076923
fatal50_54w <- 0.166666667
fatal55_59w <- 0.196116505
fatal60_64w <- 0.267532468
fatal65_69w <- 0.195708155
fatal70_74w <- 0.25875
fatal75_79w <- 0.300787402
fatal80_84w <- 0.449655172
fatal85w <- 0.749497487

p_stoma <- 0.13 #probability of sequelae: stoma


#Natural history model CRC (Soerjomataram et al. 2012)
time_diagnosis <- 1.08
time_remission_cure <- 7
time_remission_death <- 1.267
time_disseminated <- 0.25
time_terminal <- 0.083
time_total_YLL <- time_diagnosis + time_remission_death + time_disseminated + time_terminal
time_total_YLD <- time_diagnosis + time_remission_cure


#disability weights (Salomon et al., 2015)
set.seed(1)
dw_diagnosis <- mcstoc(rpert, type = "U", min=0.193, mode=0.288, max=0.339)
dw_remission <- 0.2
dw_stoma <- mcstoc(rpert, type = "U", min=0.063, mode=0.095, max=0.131)
dw_disseminated <- mcstoc(rpert, type = "U", min=0.307, mode=0.451, max=0.6)
dw_terminal <- mcstoc(rpert, type = "U", min=0.377, mode=0.54, max=0.687)


#RR dose-response
r <- mcstoc(rpert, type = "U", min = 0.0005, mode = 0.0017, max = 0.0031) #1.17, 95% CI: 1.05-1.31 converted to RR/g/day (WCRF/AICR CUP SLR 2011)




##### Ref scenario #####

FoodDaily <- read.csv("FoodDaily.csv")
MeatDaily <- FoodDaily[, c(1:4,34,45,58)]


agebreaks <- c(15,20,25,30,35,40,45,50,55,60,65,70,75,80)
agelabels= c("15-19","20-24","25-29","30-34","35-39", "40-44", "45-49", "50-54","55-59", "60-64","65-69", "70-74","75-79")

setDT(MeatDaily)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### DALY calc males #####

#Age 15-19, males

#Diagnosis
YLD_diagnosis <- p15_19m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p15_19m * (1-fatal15_19m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p15_19m * fatal15_19m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p15_19m * fatal15_19m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p15_19m * fatal15_19m * dw_terminal * time_terminal

#Sequela: stoma - assume that individuals get the stoma after the diagnosis and remission phase
YLD_stoma <- p15_19m * (1-fatal15_19m) * p_stoma * dw_stoma * (LE_15_19m - (17 + time_total_YLD))

#Total YLD
YLD_15_19m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_15_19m <- p15_19m * fatal15_19m * (SEYLL_15_19 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY15_19m <- YLL_15_19m + YLD_15_19m

DALY15_19totalm <- DALY15_19m * pop_15_19m


#Age 20-24, males

#Diagnosis
YLD_diagnosis <- p20_24m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p20_24m * (1-fatal20_24m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p20_24m * fatal20_24m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p20_24m * fatal20_24m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p20_24m * fatal20_24m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p20_24m * (1-fatal20_24m) * p_stoma * dw_stoma * (LE_20_24m - (22 + time_total_YLD))

#Total YLD
YLD_20_24m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_20_24m <- p20_24m * fatal20_24m * (SEYLL_20_24 - time_total_YLL)

#DALYs
DALY20_24m <- YLL_20_24m + YLD_20_24m

DALY20_24totalm <- DALY20_24m * pop_20_24m

#Age 25-29, males

#Diagnosis
YLD_diagnosis <- p25_29m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p25_29m * (1-fatal25_29m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p25_29m * fatal25_29m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p25_29m * fatal25_29m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p25_29m * fatal25_29m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p25_29m * (1-fatal25_29m) * p_stoma * dw_stoma * (LE_25_29m - (27 + time_total_YLD))

#Total YLD
YLD_25_29m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_25_29m <- p25_29m * fatal25_29m * (SEYLL_25_29 - time_total_YLL)

#DALYs
DALY25_29m <- YLL_25_29m + YLD_25_29m

DALY25_29totalm <- DALY25_29m * pop_25_29m

#Age 30-34, males

#Diagnosis
YLD_diagnosis <- p30_34m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p30_34m * (1-fatal30_34m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p30_34m * fatal30_34m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p30_34m * fatal30_34m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p30_34m * fatal30_34m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p30_34m * (1-fatal30_34m) * p_stoma * dw_stoma * (LE_30_34m - (32 + time_total_YLD))

#Total YLD
YLD_30_34m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_30_34m <- p30_34m * fatal30_34m * (SEYLL_30_34 - time_total_YLL)

#DALYs
DALY30_34m <- YLL_30_34m + YLD_30_34m

DALY30_34totalm <- DALY30_34m * pop_30_34m


#Age 35-39, males

#Diagnosis
YLD_diagnosis <- p35_39m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p35_39m * (1-fatal35_39m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p35_39m * fatal35_39m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p35_39m * fatal35_39m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p35_39m * fatal35_39m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p35_39m * (1-fatal35_39m) * p_stoma * dw_stoma * (LE_35_39m - (37 + time_total_YLD))

#Total YLD
YLD_35_39m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_35_39m <- p35_39m * fatal35_39m * (SEYLL_35_39 - time_total_YLL)

#DALYs
DALY35_39m <- YLL_35_39m + YLD_35_39m

DALY35_39totalm <- DALY35_39m * pop_35_39m


#Age 40-44, males

#Diagnosis
YLD_diagnosis <- p40_44m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p40_44m * (1-fatal40_44m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p40_44m * fatal40_44m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p40_44m * fatal40_44m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p40_44m * fatal40_44m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p40_44m * (1-fatal40_44m) * p_stoma * dw_stoma * (LE_40_44m - (42 + time_total_YLD))

#Total YLD
YLD_40_44m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_40_44m <- p40_44m * fatal40_44m * (SEYLL_40_44 - time_total_YLL)

#DALYs
DALY40_44m <- YLL_40_44m + YLD_40_44m

DALY40_44totalm <- DALY40_44m * pop_40_44m


#Age 45-49, males

#Diagnosis
YLD_diagnosis <- p45_49m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p45_49m * (1-fatal45_49m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p45_49m * fatal45_49m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p45_49m * fatal45_49m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p45_49m * fatal45_49m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p45_49m * (1-fatal45_49m) * p_stoma * dw_stoma * (LE_45_49m - (47 + time_total_YLD))

#Total YLD
YLD_45_49m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_45_49m <- p45_49m * fatal45_49m * (SEYLL_45_49 - time_total_YLL)

#DALYs
DALY45_49m <- YLL_45_49m + YLD_45_49m

DALY45_49totalm <- DALY45_49m * pop_45_49m


#Age 50-54, males

#Diagnosis
YLD_diagnosis <- p50_54m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p50_54m * (1-fatal50_54m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p50_54m * fatal50_54m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p50_54m * fatal50_54m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p50_54m * fatal50_54m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p50_54m * (1-fatal50_54m) * p_stoma * dw_stoma * (LE_50_54m - (52 + time_total_YLD))

#Total YLD
YLD_50_54m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_50_54m <- p50_54m * fatal50_54m * (SEYLL_50_54 - time_total_YLL)

#DALYs
DALY50_54m <- YLL_50_54m + YLD_50_54m

DALY50_54totalm <- DALY50_54m * pop_50_54m


#Age 55-59, males

#Diagnosis
YLD_diagnosis <- p55_59m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p55_59m * (1-fatal55_59m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p55_59m * fatal55_59m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p55_59m * fatal55_59m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p55_59m * fatal55_59m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p55_59m * (1-fatal55_59m) * p_stoma * dw_stoma * (LE_55_59m - (57 + time_total_YLD))

#Total YLD
YLD_55_59m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_55_59m <- p55_59m * fatal55_59m * (SEYLL_55_59 - time_total_YLL)

#DALYs
DALY55_59m <- YLL_55_59m + YLD_55_59m

DALY55_59totalm <- DALY55_59m * pop_55_59m


#Age 60-64, males

#Diagnosis
YLD_diagnosis <- p60_64m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p60_64m * (1-fatal60_64m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p60_64m * fatal60_64m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p60_64m * fatal60_64m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p60_64m * fatal60_64m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p60_64m * (1-fatal60_64m) * p_stoma * dw_stoma * (LE_60_64m - (62 + time_total_YLD))

#Total YLD
YLD_60_64m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_60_64m <- p60_64m * fatal60_64m * (SEYLL_60_64 - time_total_YLL)

#DALYs
DALY60_64m <- YLL_60_64m + YLD_60_64m

DALY60_64totalm <- DALY60_64m * pop_60_64m


#Age 65-69, males

#Diagnosis
YLD_diagnosis <- p65_69m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p65_69m * (1-fatal65_69m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p65_69m * fatal65_69m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p65_69m * fatal65_69m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p65_69m * fatal65_69m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p65_69m * (1-fatal65_69m) * p_stoma * dw_stoma * (LE_65_69m - (67 + time_total_YLD))

#Total YLD
YLD_65_69m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_65_69m <- p65_69m * fatal65_69m * (SEYLL_65_69 - time_total_YLL)

#DALYs
DALY65_69m <- YLL_65_69m + YLD_65_69m

DALY65_69totalm <- DALY65_69m * pop_65_69m


#Age 70-74, males

#Diagnosis
YLD_diagnosis <- p70_74m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p70_74m * (1-fatal70_74m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p70_74m * fatal70_74m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p70_74m * fatal70_74m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p70_74m * fatal70_74m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p70_74m * (1-fatal70_74m) * p_stoma * dw_stoma * (LE_70_74m - (72 + time_total_YLD))

#Total YLD
YLD_70_74m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_70_74m <- p70_74m * fatal70_74m * (SEYLL_70_74 - time_total_YLL)

#DALYs
DALY70_74m <- YLL_70_74m + YLD_70_74m

DALY70_74totalm <- DALY70_74m * pop_70_74m


#Age 75-79, males

#Diagnosis
YLD_diagnosis <- p75_79m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p75_79m * (1-fatal75_79m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p75_79m * fatal75_79m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p75_79m * fatal75_79m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p75_79m * fatal75_79m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p75_79m * (1-fatal75_79m) * p_stoma * dw_stoma * (LE_75_79m - (77 + time_total_YLD))

#Total YLD
YLD_75_79m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_75_79m <- p75_79m * fatal75_79m * (SEYLL_75_79 - time_total_YLL)

#DALYs
DALY75_79m <- YLL_75_79m + YLD_75_79m

DALY75_79totalm <- DALY75_79m * pop_75_79m


#Age 80-84, males

#Diagnosis
YLD_diagnosis <- p80_84m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p80_84m * (1-fatal80_84m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p80_84m * fatal80_84m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p80_84m * fatal80_84m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p80_84m * fatal80_84m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p80_84m * (1-fatal80_84m) * p_stoma * dw_stoma * (LE_80_84m - (82 + time_total_YLD))

#Total YLD
YLD_80_84m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_80_84m <- p80_84m * fatal80_84m * (SEYLL_80_84 - time_total_YLL)

#DALYs
DALY80_84m <- YLL_80_84m + YLD_80_84m

DALY80_84totalm <- DALY80_84m * pop_80_84m


#Age 85+, males

#Diagnosis
YLD_diagnosis <- p85m * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p85m * (1-fatal85m) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p85m * fatal85m * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p85m * fatal85m * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p85m * fatal85m * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p85m * (1-fatal85m) * p_stoma * dw_stoma * (LE_85m - (87 + time_total_YLD))

#Total YLD
YLD_85m <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_85m <- p85m * fatal85m * (SEYLL_85 - time_total_YLL) 

#DALYs
DALY85m <- YLL_85m + YLD_85m

DALY85totalm <- DALY85m * pop_85m


##### DALY calc females #####

#Age 15-19, females

#Diagnosis
YLD_diagnosis <- p15_19w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p15_19w * (1-fatal15_19w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p15_19w * fatal15_19w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p15_19w * fatal15_19w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p15_19w * fatal15_19w * dw_terminal * time_terminal

#Sequela: stoma - assume that individuals get the stoma after the diagnosis and remission phase
YLD_stoma <- p15_19w * (1-fatal15_19w) * p_stoma * dw_stoma * (LE_15_19w - (17 + time_total_YLD))

#Total YLD
YLD_15_19w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_15_19w <- p15_19w * fatal15_19w * (SEYLL_15_19 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY15_19w <- YLL_15_19w + YLD_15_19w

DALY15_19totalw <- DALY15_19w * pop_15_19w


#Age 20-24, females

#Diagnosis
YLD_diagnosis <- p20_24w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p20_24w * (1-fatal20_24w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p20_24w * fatal20_24w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p20_24w * fatal20_24w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p20_24w * fatal20_24w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p20_24w * (1-fatal20_24w) * p_stoma * dw_stoma * (LE_20_24w - (22 + time_total_YLD))

#Total YLD
YLD_20_24w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_20_24w <- p20_24w * fatal20_24w * (SEYLL_20_24 - time_total_YLL)

#DALYs
DALY20_24w <- YLL_20_24w + YLD_20_24w

DALY20_24totalw <- DALY20_24w * pop_20_24w

#Age 25-29, females

#Diagnosis
YLD_diagnosis <- p25_29w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p25_29w * (1-fatal25_29w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p25_29w * fatal25_29w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p25_29w * fatal25_29w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p25_29w * fatal25_29w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p25_29w * (1-fatal25_29w) * p_stoma * dw_stoma * (LE_25_29w - (27 + time_total_YLD))

#Total YLD
YLD_25_29w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_25_29w <- p25_29w * fatal25_29w * (SEYLL_25_29 - time_total_YLL)

#DALYs
DALY25_29w <- YLL_25_29w + YLD_25_29w

DALY25_29totalw <- DALY25_29w * pop_25_29w

#Age 30-34, females

#Diagnosis
YLD_diagnosis <- p30_34w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p30_34w * (1-fatal30_34w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p30_34w * fatal30_34w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p30_34w * fatal30_34w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p30_34w * fatal30_34w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p30_34w * (1-fatal30_34w) * p_stoma * dw_stoma * (LE_30_34w - (32 + time_total_YLD))

#Total YLD
YLD_30_34w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_30_34w <- p30_34w * fatal30_34w * (SEYLL_30_34 - time_total_YLL)

#DALYs
DALY30_34w <- YLL_30_34w + YLD_30_34w

DALY30_34totalw <- DALY30_34w * pop_30_34w


#Age 35-39, females

#Diagnosis
YLD_diagnosis <- p35_39w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p35_39w * (1-fatal35_39w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p35_39w * fatal35_39w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p35_39w * fatal35_39w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p35_39w * fatal35_39w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p35_39w * (1-fatal35_39w) * p_stoma * dw_stoma * (LE_35_39w - (37 + time_total_YLD))

#Total YLD
YLD_35_39w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_35_39w <- p35_39w * fatal35_39w * (SEYLL_35_39 - time_total_YLL)

#DALYs
DALY35_39w <- YLL_35_39w + YLD_35_39w

DALY35_39totalw <- DALY35_39w * pop_35_39w


#Age 40-44, females

#Diagnosis
YLD_diagnosis <- p40_44w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p40_44w * (1-fatal40_44w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p40_44w * fatal40_44w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p40_44w * fatal40_44w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p40_44w * fatal40_44w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p40_44w * (1-fatal40_44w) * p_stoma * dw_stoma * (LE_40_44w - (42 + time_total_YLD))

#Total YLD
YLD_40_44w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_40_44w <- p40_44w * fatal40_44w * (SEYLL_40_44 - time_total_YLL)

#DALYs
DALY40_44w <- YLL_40_44w + YLD_40_44w

DALY40_44totalw <- DALY40_44w * pop_40_44w


#Age 45-49, females

#Diagnosis
YLD_diagnosis <- p45_49w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p45_49w * (1-fatal45_49w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p45_49w * fatal45_49w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p45_49w * fatal45_49w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p45_49w * fatal45_49w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p45_49w * (1-fatal45_49w) * p_stoma * dw_stoma * (LE_45_49w - (47 + time_total_YLD))

#Total YLD
YLD_45_49w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_45_49w <- p45_49w * fatal45_49w * (SEYLL_45_49 - time_total_YLL)

#DALYs
DALY45_49w <- YLL_45_49w + YLD_45_49w

DALY45_49totalw <- DALY45_49w * pop_45_49w


#Age 50-54, females

#Diagnosis
YLD_diagnosis <- p50_54w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p50_54w * (1-fatal50_54w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p50_54w * fatal50_54w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p50_54w * fatal50_54w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p50_54w * fatal50_54w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p50_54w * (1-fatal50_54w) * p_stoma * dw_stoma * (LE_50_54w - (52 + time_total_YLD))

#Total YLD
YLD_50_54w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_50_54w <- p50_54w * fatal50_54w * (SEYLL_50_54 - time_total_YLL)

#DALYs
DALY50_54w <- YLL_50_54w + YLD_50_54w

DALY50_54totalw <- DALY50_54w * pop_50_54w


#Age 55-59, females

#Diagnosis
YLD_diagnosis <- p55_59w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p55_59w * (1-fatal55_59w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p55_59w * fatal55_59w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p55_59w * fatal55_59w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p55_59w * fatal55_59w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p55_59w * (1-fatal55_59w) * p_stoma * dw_stoma * (LE_55_59w - (57 + time_total_YLD))

#Total YLD
YLD_55_59w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_55_59w <- p55_59w * fatal55_59w * (SEYLL_55_59 - time_total_YLL)

#DALYs
DALY55_59w <- YLL_55_59w + YLD_55_59w

DALY55_59totalw <- DALY55_59w * pop_55_59w


#Age 60-64, females

#Diagnosis
YLD_diagnosis <- p60_64w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p60_64w * (1-fatal60_64w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p60_64w * fatal60_64w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p60_64w * fatal60_64w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p60_64w * fatal60_64w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p60_64w * (1-fatal60_64w) * p_stoma * dw_stoma * (LE_60_64w - (62 + time_total_YLD))

#Total YLD
YLD_60_64w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_60_64w <- p60_64w * fatal60_64w * (SEYLL_60_64 - time_total_YLL)

#DALYs
DALY60_64w <- YLL_60_64w + YLD_60_64w

DALY60_64totalw <- DALY60_64w * pop_60_64w


#Age 65-69, females

#Diagnosis
YLD_diagnosis <- p65_69w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p65_69w * (1-fatal65_69w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p65_69w * fatal65_69w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p65_69w * fatal65_69w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p65_69w * fatal65_69w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p65_69w * (1-fatal65_69w) * p_stoma * dw_stoma * (LE_65_69w - (67 + time_total_YLD))

#Total YLD
YLD_65_69w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_65_69w <- p65_69w * fatal65_69w * (SEYLL_65_69 - time_total_YLL)

#DALYs
DALY65_69w <- YLL_65_69w + YLD_65_69w

DALY65_69totalw <- DALY65_69w * pop_65_69w


#Age 70-74, females

#Diagnosis
YLD_diagnosis <- p70_74w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p70_74w * (1-fatal70_74w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p70_74w * fatal70_74w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p70_74w * fatal70_74w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p70_74w * fatal70_74w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p70_74w * (1-fatal70_74w) * p_stoma * dw_stoma * (LE_70_74w - (72 + time_total_YLD))

#Total YLD
YLD_70_74w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_70_74w <- p70_74w * fatal70_74w * (SEYLL_70_74 - time_total_YLL)

#DALYs
DALY70_74w <- YLL_70_74w + YLD_70_74w

DALY70_74totalw <- DALY70_74w * pop_70_74w


#Age 75-79, females

#Diagnosis
YLD_diagnosis <- p75_79w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p75_79w * (1-fatal75_79w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p75_79w * fatal75_79w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p75_79w * fatal75_79w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p75_79w * fatal75_79w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p75_79w * (1-fatal75_79w) * p_stoma * dw_stoma * (LE_75_79w - (77 + time_total_YLD))

#Total YLD
YLD_75_79w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_75_79w <- p75_79w * fatal75_79w * (SEYLL_75_79 - time_total_YLL)

#DALYs
DALY75_79w <- YLL_75_79w + YLD_75_79w

DALY75_79totalw <- DALY75_79w * pop_75_79w


#Age 80-84, females

#Diagnosis
YLD_diagnosis <- p80_84w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p80_84w * (1-fatal80_84w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p80_84w * fatal80_84w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p80_84w * fatal80_84w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p80_84w * fatal80_84w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p80_84w * (1-fatal80_84w) * p_stoma * dw_stoma * (LE_80_84w - (82 + time_total_YLD))

#Total YLD
YLD_80_84w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_80_84w <- p80_84w * fatal80_84w * (SEYLL_80_84 - time_total_YLL)

#DALYs
DALY80_84w <- YLL_80_84w + YLD_80_84w

DALY80_84totalw <- DALY80_84w * pop_80_84w


#Age 85+, females

#Diagnosis
YLD_diagnosis <- p85w * dw_diagnosis * time_diagnosis

#Remission, cure
YLD_remission_cure <- p85w * (1-fatal85w) * dw_remission * time_remission_cure

#Remission, death
YLD_remission_death <- p85w * fatal85w * dw_remission * time_remission_death

#Disseminated
YLD_disseminated <- p85w * fatal85w * dw_disseminated * time_disseminated

#Terminal
YLD_terminal <- p85w * fatal85w * dw_terminal * time_terminal

#Sequela: stoma
YLD_stoma <- p85w * (1-fatal85w) * p_stoma * dw_stoma * (LE_85w - (87 + time_total_YLD))

#Total YLD
YLD_85w <- YLD_diagnosis + YLD_remission_cure + YLD_remission_death + YLD_disseminated + YLD_terminal + YLD_stoma

#YLLs
YLL_85w <- p85w * fatal85w * (SEYLL_85 - time_total_YLL) 

#DALYs
DALY85w <- YLL_85w + YLD_85w

DALY85totalw <- DALY85w * pop_85w


##### Total DALYs ref #####

DALYtotalref <- DALY15_19totalm + DALY20_24totalm + DALY25_29totalm + DALY30_34totalm + DALY35_39totalm +
  DALY40_44totalm + DALY45_49totalm + DALY50_54totalm + DALY55_59totalm + DALY60_64totalm + DALY65_69totalm +
  DALY70_74totalm + DALY75_79totalm + DALY80_84totalm + DALY85totalm + 
  DALY15_19totalw + DALY20_24totalw + DALY25_29totalw + DALY30_34totalw + DALY35_39totalw +
  DALY40_44totalw + DALY45_49totalw + DALY50_54totalw + DALY55_59totalw + DALY60_64totalw + DALY65_69totalw +
  DALY70_74totalw + DALY75_79totalw + DALY80_84totalw + DALY85totalw

summary(DALYtotalref)


##### Total cases reference scenario #####

cases_ref <- p15_19m * pop_15_19m + p20_24m * pop_20_24m + p25_29m * pop_25_29m + p30_34m * pop_30_34m + p35_39m * pop_35_39m +
  p40_44m * pop_40_44m + p45_49m * pop_45_49m + p50_54m * pop_50_54m + p55_59m * pop_55_59m + p60_64m * pop_60_64m + 
  p65_69m * pop_65_69m + p70_74m * pop_70_74m + p75_79m * pop_75_79m + p80_84m * pop_80_84m + p85m * pop_85m +
  p15_19w * pop_15_19w + p20_24w * pop_20_24w + p25_29w * pop_25_29w + p30_34w * pop_30_34w + p35_39w * pop_35_39w +
  p40_44w * pop_40_44w + p45_49w * pop_45_49w + p50_54w * pop_50_54w + p55_59w * pop_55_59w + p60_64w * pop_60_64w + 
  p65_69w * pop_65_69w + p70_74w * pop_70_74w + p75_79w * pop_75_79w + p80_84w * pop_80_84w + p85w * pop_85w

summary(cases_ref)


##### Fitting exposures males #####

#Divide data into agegroups and fit exposure to distribution

#Age 15-19, male
red15_19m <- subset(MeatDaily, agegroups=="15-19" & sex =="1", select = red.meat)

#Probability of intake
probred15_19m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red15_19m==0)/length(red15_19m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red15_19m!=0)/length(red15_19m$red.meat) #probability of DHA.EPA exp
probred15_19m[1,] <- c(pr.no,pr.yes)

red15_19posm <- red15_19m[which(red15_19m!=0)]
red15_19posm <- as.vector(red15_19posm$red.meat)


fit1_15_19m <- fitdist(red15_19posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red15_19posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_15_19m$estimate[1], sdlog=fit1_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19m, fitnames="lnorm") 

fit2_15_19m <- fitdist(red15_19posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red15_19posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2])  #Anderson-Darling test


#Age 20-24, male
red20_24m <- subset(MeatDaily, agegroups=="20-24" & sex =="1", select = red.meat)

#Probability of intake
probred20_24m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red20_24m==0)/length(red20_24m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red20_24m!=0)/length(red20_24m$red.meat) #probability of DHA.EPA exp
probred20_24m[1,] <- c(pr.no,pr.yes)

red20_24posm <- red20_24m[which(red20_24m!=0)]
red20_24posm <- as.vector(red20_24posm$red.meat)

fit1_20_24m <- fitdist(red20_24posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red20_24posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_20_24m$estimate[1], sdlog=fit1_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24m, fitnames="lnorm") 

fit2_20_24m <- fitdist(red20_24posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red20_24posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Anderson-Darling test


#Age 25-29, male
red25_29m <- subset(MeatDaily, agegroups=="25-29" & sex =="1", select = red.meat)

#Probability of intake
probred25_29m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red25_29m==0)/length(red25_29m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red25_29m!=0)/length(red25_29m$red.meat) #probability of DHA.EPA exp
probred25_29m[1,] <- c(pr.no,pr.yes)

red25_29posm <- red25_29m[which(red25_29m!=0)]
red25_29posm <- as.vector(red25_29posm$red.meat)

fit1_25_29m <- fitdist(red25_29posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red25_29posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_25_29m$estimate[1], sdlog=fit1_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29m, fitnames="lnorm") 

fit2_25_29m <- fitdist(red25_29posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red25_29posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2])  #Anderson-Darling test


#Age 30-34, male
red30_34m <- subset(MeatDaily, agegroups=="30-34" & sex =="1", select = red.meat)

#Probability of intake
probred30_34m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red30_34m==0)/length(red30_34m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red30_34m!=0)/length(red30_34m$red.meat) #probability of DHA.EPA exp
probred30_34m[1,] <- c(pr.no,pr.yes)

red30_34posm <- red30_34m[which(red30_34m!=0)]
red30_34posm <- as.vector(red30_34posm$red.meat)

fit1_30_34m <- fitdist(red30_34posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red30_34posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_30_34m$estimate[1], sdlog=fit1_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34m, fitnames="lnorm") 

fit2_30_34m <- fitdist(red30_34posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red30_34posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Anderson-Darling test


#Age 35-39, male
red35_39m <- subset(MeatDaily, agegroups=="35-39" & sex =="1", select = red.meat)

#Probability of intake
probred35_39m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red35_39m==0)/length(red35_39m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red35_39m!=0)/length(red35_39m$red.meat) #probability of DHA.EPA exp
probred35_39m[1,] <- c(pr.no,pr.yes)

red35_39posm <- red35_39m[which(red35_39m!=0)]
red35_39posm <- as.vector(red35_39posm$red.meat)

fit1_35_39m <- fitdist(red35_39posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red35_39posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_35_39m$estimate[1], sdlog=fit1_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39m, fitnames="lnorm") 

fit2_35_39m <- fitdist(red35_39posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red35_39posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Anderson-Darling test


#Age 40-44, male
red40_44m <- subset(MeatDaily, agegroups=="40-44" & sex =="1", select = red.meat)

#Probability of intake
probred40_44m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red40_44m==0)/length(red40_44m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red40_44m!=0)/length(red40_44m$red.meat) #probability of DHA.EPA exp
probred40_44m[1,] <- c(pr.no,pr.yes)

red40_44posm <- red40_44m[which(red40_44m!=0)]
red40_44posm <- as.vector(red40_44posm$red.meat)


fit1_40_44m <- fitdist(red40_44posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red40_44posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_40_44m$estimate[1], sdlog=fit1_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44m, fitnames="lnorm") 

fit2_40_44m <- fitdist(red40_44posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red40_44posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Anderson-Darling test


#Age 45-49, male
red45_49m <- subset(MeatDaily, agegroups=="45-49" & sex =="1", select = red.meat)

#Probability of intake
probred45_49m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red45_49m==0)/length(red45_49m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red45_49m!=0)/length(red45_49m$red.meat) #probability of DHA.EPA exp
probred45_49m[1,] <- c(pr.no,pr.yes)

red45_49posm <- red45_49m[which(red45_49m!=0)]
red45_49posm <- as.vector(red45_49posm$red.meat)


fit1_45_49m <- fitdist(red45_49posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red45_49posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_45_49m$estimate[1], sdlog=fit1_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49m, fitnames="lnorm") 

fit2_45_49m <- fitdist(red45_49posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red45_49posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Anderson-Darling test


#Age 50-54, male
red50_54m <- subset(MeatDaily, agegroups=="50-54" & sex =="1", select = red.meat)

#Probability of intake
probred50_54m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red50_54m==0)/length(red50_54m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red50_54m!=0)/length(red50_54m$red.meat) #probability of DHA.EPA exp
probred50_54m[1,] <- c(pr.no,pr.yes)

red50_54posm <- red50_54m[which(red50_54m!=0)]
red50_54posm <- as.vector(red50_54posm$red.meat)


fit1_50_54m <- fitdist(red50_54posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red50_54posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_50_54m$estimate[1], sdlog=fit1_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54m, fitnames="lnorm") 

fit2_50_54m <- fitdist(red50_54posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red50_54posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54m, fitnames="gamma") 


cvm.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Anderson-Darling test


#Age 55-59, male
red55_59m <- subset(MeatDaily, agegroups=="55-59" & sex =="1", select = red.meat)

#Probability of intake
probred55_59m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red55_59m==0)/length(red55_59m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red55_59m!=0)/length(red55_59m$red.meat) #probability of DHA.EPA exp
probred55_59m[1,] <- c(pr.no,pr.yes)

red55_59posm <- red55_59m[which(red55_59m!=0)]
red55_59posm <- as.vector(red55_59posm$red.meat)


fit1_55_59m <- fitdist(red55_59posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red55_59posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_55_59m$estimate[1], sdlog=fit1_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59m, fitnames="lnorm") 

fit2_55_59m <- fitdist(red55_59posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red55_59posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59m, fitnames="gamma") 


cvm.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Anderson-Darling test


#Age 60-64, male
red60_64m <- subset(MeatDaily, agegroups=="60-64" & sex =="1", select = red.meat)

#Probability of intake
probred60_64m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red60_64m==0)/length(red60_64m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red60_64m!=0)/length(red60_64m$red.meat) #probability of DHA.EPA exp
probred60_64m[1,] <- c(pr.no,pr.yes)

red60_64posm <- red60_64m[which(red60_64m!=0)]
red60_64posm <- as.vector(red60_64posm$red.meat)


fit1_60_64m <- fitdist(red60_64posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red60_64posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_60_64m$estimate[1], sdlog=fit1_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64m, fitnames="lnorm") 

fit2_60_64m <- fitdist(red60_64posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red60_64posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Anderson-Darling test


#Age 65-69, male
red65_69m <- subset(MeatDaily, agegroups=="65-69" & sex =="1", select = red.meat)

#Probability of intake
probred65_69m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red65_69m==0)/length(red65_69m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red65_69m!=0)/length(red65_69m$red.meat) #probability of DHA.EPA exp
probred65_69m[1,] <- c(pr.no,pr.yes)

red65_69posm <- red65_69m[which(red65_69m!=0)]
red65_69posm <- as.vector(red65_69posm$red.meat)

fit1_65_69m <- fitdist(red65_69posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red65_69posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_65_69m$estimate[1], sdlog=fit1_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69m, fitnames="lnorm") 

fit2_65_69m <- fitdist(red65_69posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red65_69posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Anderson-Darling test


#Age 70-74, male
red70_74m <- subset(MeatDaily, agegroups=="70-74" & sex =="1", select = red.meat)

#Probability of intake
probred70_74m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red70_74m==0)/length(red70_74m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red70_74m!=0)/length(red70_74m$red.meat) #probability of DHA.EPA exp
probred70_74m[1,] <- c(pr.no,pr.yes)

red70_74posm <- red70_74m[which(red70_74m!=0)]
red70_74posm <- as.vector(red70_74posm$red.meat)

fit1_70_74m <- fitdist(red70_74posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red70_74posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_70_74m$estimate[1], sdlog=fit1_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74m, fitnames="lnorm") 

fit2_70_74m <- fitdist(red70_74posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red70_74posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Anderson-Darling test


#Age 75-79, male
red75_79m <- subset(MeatDaily, agegroups=="75-79" & sex =="1", select = red.meat)

#Probability of intake
probred75_79m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red75_79m==0)/length(red75_79m$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red75_79m!=0)/length(red75_79m$red.meat) #probability of DHA.EPA exp
probred75_79m[1,] <- c(pr.no,pr.yes)

red75_79posm <- red75_79m[which(red75_79m!=0)]
red75_79posm <- as.vector(red75_79posm$red.meat)

fit1_75_79m <- fitdist(red75_79posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red75_79posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_75_79m$estimate[1], sdlog=fit1_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79m, fitnames="lnorm") 

fit2_75_79m <- fitdist(red75_79posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red75_79posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Anderson-Darling test


##### Fitting exposures females #####

#Divide data into agegroups and fit exposure to distribution

#Age 15-19, female
red15_19w <- subset(MeatDaily, agegroups=="15-19" & sex =="2", select = red.meat)

#Probability of intake
probred15_19w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red15_19w==0)/length(red15_19w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red15_19w!=0)/length(red15_19w$red.meat) #probability of DHA.EPA exp
probred15_19w[1,] <- c(pr.no,pr.yes)

red15_19posw <- red15_19w[which(red15_19w!=0)]
red15_19posw <- as.vector(red15_19posw$red.meat)


fit1_15_19w <- fitdist(red15_19posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red15_19posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_15_19w$estimate[1], sdlog=fit1_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19w, fitnames="lnorm") 

fit2_15_19w <- fitdist(red15_19posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red15_19posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19w, fitnames="gamma") 


cvm.test(t, plnorm, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2])  #Anderson-Darling test


#Age 20-24, female
red20_24w <- subset(MeatDaily, agegroups=="20-24" & sex =="2", select = red.meat)

#Probability of intake
probred20_24w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red20_24w==0)/length(red20_24w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red20_24w!=0)/length(red20_24w$red.meat) #probability of DHA.EPA exp
probred20_24w[1,] <- c(pr.no,pr.yes)

red20_24posw <- red20_24w[which(red20_24w!=0)]
red20_24posw <- as.vector(red20_24posw$red.meat)

fit1_20_24w <- fitdist(red20_24posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red20_24posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_20_24w$estimate[1], sdlog=fit1_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24w, fitnames="lnorm") 

fit2_20_24w <- fitdist(red20_24posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red20_24posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Anderson-Darling test


#Age 25-29, female
red25_29w <- subset(MeatDaily, agegroups=="25-29" & sex =="2", select = red.meat)

#Probability of intake
probred25_29w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red25_29w==0)/length(red25_29w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red25_29w!=0)/length(red25_29w$red.meat) #probability of DHA.EPA exp
probred25_29w[1,] <- c(pr.no,pr.yes)

red25_29posw <- red25_29w[which(red25_29w!=0)]
red25_29posw <- as.vector(red25_29posw$red.meat)

fit1_25_29w <- fitdist(red25_29posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red25_29posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_25_29w$estimate[1], sdlog=fit1_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29w, fitnames="lnorm") 

fit2_25_29w <- fitdist(red25_29posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red25_29posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Anderson-Darling test


#Age 30-34, female
red30_34w <- subset(MeatDaily, agegroups=="30-34" & sex =="2", select = red.meat)

#Probability of intake
probred30_34w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red30_34w==0)/length(red30_34w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red30_34w!=0)/length(red30_34w$red.meat) #probability of DHA.EPA exp
probred30_34w[1,] <- c(pr.no,pr.yes)

red30_34posw <- red30_34w[which(red30_34w!=0)]
red30_34posw <- as.vector(red30_34posw$red.meat)

fit1_30_34w <- fitdist(red30_34posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red30_34posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_30_34w$estimate[1], sdlog=fit1_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34w, fitnames="lnorm") 

fit2_30_34w <- fitdist(red30_34posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red30_34posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2])  #Anderson-Darling test


#Age 35-39, female
red35_39w <- subset(MeatDaily, agegroups=="35-39" & sex =="2", select = red.meat)

#Probability of intake
probred35_39w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red35_39w==0)/length(red35_39w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red35_39w!=0)/length(red35_39w$red.meat) #probability of DHA.EPA exp
probred35_39w[1,] <- c(pr.no,pr.yes)

red35_39posw <- red35_39w[which(red35_39w!=0)]
red35_39posw <- as.vector(red35_39posw$red.meat)

fit1_35_39w <- fitdist(red35_39posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red35_39posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_35_39w$estimate[1], sdlog=fit1_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39w, fitnames="lnorm") 

fit2_35_39w <- fitdist(red35_39posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red35_39posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Anderson-Darling test


#Age 40-44, female
red40_44w <- subset(MeatDaily, agegroups=="40-44" & sex =="2", select = red.meat)

#Probability of intake
probred40_44w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red40_44w==0)/length(red40_44w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red40_44w!=0)/length(red40_44w$red.meat) #probability of DHA.EPA exp
probred40_44w[1,] <- c(pr.no,pr.yes)

red40_44posw <- red40_44w[which(red40_44w!=0)]
red40_44posw <- as.vector(red40_44posw$red.meat)


fit1_40_44w <- fitdist(red40_44posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red40_44posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_40_44w$estimate[1], sdlog=fit1_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44w, fitnames="lnorm") 

fit2_40_44w <- fitdist(red40_44posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red40_44posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Anderson-Darling test


#Age 45-49, female
red45_49w <- subset(MeatDaily, agegroups=="45-49" & sex =="2", select = red.meat)

#Probability of intake
probred45_49w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red45_49w==0)/length(red45_49w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red45_49w!=0)/length(red45_49w$red.meat) #probability of DHA.EPA exp
probred45_49w[1,] <- c(pr.no,pr.yes)

red45_49posw <- red45_49w[which(red45_49w!=0)]
red45_49posw <- as.vector(red45_49posw$red.meat)


fit1_45_49w <- fitdist(red45_49posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red45_49posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_45_49w$estimate[1], sdlog=fit1_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49w, fitnames="lnorm") 

fit2_45_49w <- fitdist(red45_49posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red45_49posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Anderson-Darling test


#Age 50-54, female
red50_54w <- subset(MeatDaily, agegroups=="50-54" & sex =="2", select = red.meat)

#Probability of intake
probred50_54w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red50_54w==0)/length(red50_54w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red50_54w!=0)/length(red50_54w$red.meat) #probability of DHA.EPA exp
probred50_54w[1,] <- c(pr.no,pr.yes)

red50_54posw <- red50_54w[which(red50_54w!=0)]
red50_54posw <- as.vector(red50_54posw$red.meat)


fit1_50_54w <- fitdist(red50_54posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red50_54posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_50_54w$estimate[1], sdlog=fit1_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54w, fitnames="lnorm") 

fit2_50_54w <- fitdist(red50_54posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red50_54posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54w, fitnames="gamma") 


cvm.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Anderson-Darling test


#Age 55-59, female
red55_59w <- subset(MeatDaily, agegroups=="55-59" & sex =="2", select = red.meat)

#Probability of intake
probred55_59w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red55_59w==0)/length(red55_59w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red55_59w!=0)/length(red55_59w$red.meat) #probability of DHA.EPA exp
probred55_59w[1,] <- c(pr.no,pr.yes)

red55_59posw <- red55_59w[which(red55_59w!=0)]
red55_59posw <- as.vector(red55_59posw$red.meat)


fit1_55_59w <- fitdist(red55_59posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red55_59posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_55_59w$estimate[1], sdlog=fit1_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59w, fitnames="lnorm") 

fit2_55_59w <- fitdist(red55_59posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red55_59posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59w, fitnames="gamma") 


cvm.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Anderson-Darling test


#Age 60-64, female
red60_64w <- subset(MeatDaily, agegroups=="60-64" & sex =="2", select = red.meat)

#Probability of intake
probred60_64w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red60_64w==0)/length(red60_64w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red60_64w!=0)/length(red60_64w$red.meat) #probability of DHA.EPA exp
probred60_64w[1,] <- c(pr.no,pr.yes)

red60_64posw <- red60_64w[which(red60_64w!=0)]
red60_64posw <- as.vector(red60_64posw$red.meat)


fit1_60_64w <- fitdist(red60_64posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red60_64posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_60_64w$estimate[1], sdlog=fit1_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64w, fitnames="lnorm") 

fit2_60_64w <- fitdist(red60_64posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red60_64posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Anderson-Darling test


#Age 65-69, female
red65_69w <- subset(MeatDaily, agegroups=="65-69" & sex =="2", select = red.meat)

#Probability of intake
probred65_69w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red65_69w==0)/length(red65_69w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red65_69w!=0)/length(red65_69w$red.meat) #probability of DHA.EPA exp
probred65_69w[1,] <- c(pr.no,pr.yes)

red65_69posw <- red65_69w[which(red65_69w!=0)]
red65_69posw <- as.vector(red65_69posw$red.meat)

fit1_65_69w <- fitdist(red65_69posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red65_69posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_65_69w$estimate[1], sdlog=fit1_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69w, fitnames="lnorm") 

fit2_65_69w <- fitdist(red65_69posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red65_69posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Anderson-Darling test


#Age 70-74, female
red70_74w <- subset(MeatDaily, agegroups=="70-74" & sex =="2", select = red.meat)

#Probability of intake
probred70_74w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red70_74w==0)/length(red70_74w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red70_74w!=0)/length(red70_74w$red.meat) #probability of DHA.EPA exp
probred70_74w[1,] <- c(pr.no,pr.yes)

red70_74posw <- red70_74w[which(red70_74w!=0)]
red70_74posw <- as.vector(red70_74posw$red.meat)

fit1_70_74w <- fitdist(red70_74posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red70_74posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_70_74w$estimate[1], sdlog=fit1_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74w, fitnames="lnorm") 

fit2_70_74w <- fitdist(red70_74posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red70_74posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Anderson-Darling test


#Age 75-79, female
red75_79w <- subset(MeatDaily, agegroups=="75-79" & sex =="2", select = red.meat)

#Probability of intake
probred75_79w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red75_79w==0)/length(red75_79w$red.meat) #probability of zero DHA.EPA exp
pr.yes <- sum(red75_79w!=0)/length(red75_79w$red.meat) #probability of DHA.EPA exp
probred75_79w[1,] <- c(pr.no,pr.yes)

red75_79posw <- red75_79w[which(red75_79w!=0)]
red75_79posw <- as.vector(red75_79posw$red.meat)

fit1_75_79w <- fitdist(red75_79posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red75_79posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_75_79w$estimate[1], sdlog=fit1_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79w, fitnames="lnorm") 

fit2_75_79w <- fitdist(red75_79posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red75_79posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Anderson-Darling test


##### Scenario 1 #####

Scenario1 <- read.csv("Scenario1.csv") # Since we're considering redessed meat as a whole, the results for this endpoint will be the same
# for scenario 1,2,3,4

Alt_scenario <- Scenario1[,c(1:4,8:10)]
setDT(Alt_scenario)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Fitting exposures males #####

#Age 15-19, male
red15_19m <- subset(Alt_scenario, agegroups=="15-19" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt15_19m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red15_19m==0)/length(red15_19m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red15_19m!=0)/length(red15_19m$redmeat.new) #probability of DHA.EPA exp
probred_alt15_19m[1,] <- c(pr.no,pr.yes)

red15_19posm <- red15_19m[which(red15_19m!=0)]
red15_19posm <- as.vector(red15_19posm$redmeat.new)

fit1_alt_15_19m <- fitdist(red15_19posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red15_19posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_15_19m$estimate[1], sdlog=fit1_alt_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_15_19m, fitnames="lnorm") 

fit2_alt_15_19m <- fitdist(red15_19posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red15_19posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_15_19m$estimate[1], fit2_alt_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_15_19m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_15_19m$estimate[1], fit1_alt_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_15_19m$estimate[1], fit1_alt_15_19m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_15_19m$estimate[1], fit2_alt_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_15_19m$estimate[1], fit2_alt_15_19m$estimate[2])  #Anderson-Darling test


#Age 20-24, male
red20_24m <- subset(Alt_scenario, agegroups=="20-24" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt20_24m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red20_24m==0)/length(red20_24m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red20_24m!=0)/length(red20_24m$redmeat.new) #probability of DHA.EPA exp
probred_alt20_24m[1,] <- c(pr.no,pr.yes)

red20_24posm <- red20_24m[which(red20_24m!=0)]
red20_24posm <- as.vector(red20_24posm$redmeat.new)

fit1_alt_20_24m <- fitdist(red20_24posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red20_24posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_20_24m$estimate[1], sdlog=fit1_alt_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_20_24m, fitnames="lnorm") 

fit2_alt_20_24m <- fitdist(red20_24posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red20_24posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_20_24m$estimate[1], fit2_alt_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_20_24m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_20_24m$estimate[1], fit1_alt_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_20_24m$estimate[1], fit1_alt_20_24m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_20_24m$estimate[1], fit2_alt_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_20_24m$estimate[1], fit2_alt_20_24m$estimate[2])  #Anderson-Darling test


#Age 25-29, male
red25_29m <- subset(Alt_scenario, agegroups=="25-29" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt25_29m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red25_29m==0)/length(red25_29m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red25_29m!=0)/length(red25_29m$redmeat.new) #probability of DHA.EPA exp
probred_alt25_29m[1,] <- c(pr.no,pr.yes)

red25_29posm <- red25_29m[which(red25_29m!=0)]
red25_29posm <- as.vector(red25_29posm$redmeat.new)

fit1_alt_25_29m <- fitdist(red25_29posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red25_29posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_25_29m$estimate[1], sdlog=fit1_alt_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_25_29m, fitnames="lnorm") 

fit2_alt_25_29m <- fitdist(red25_29posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red25_29posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_25_29m$estimate[1], fit2_alt_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_25_29m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_25_29m$estimate[1], fit1_alt_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_25_29m$estimate[1], fit1_alt_25_29m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_25_29m$estimate[1], fit2_alt_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_25_29m$estimate[1], fit2_alt_25_29m$estimate[2])  #Anderson-Darling test


#Age 30-34, male
red30_34m <- subset(Alt_scenario, agegroups=="30-34" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt30_34m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red30_34m==0)/length(red30_34m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red30_34m!=0)/length(red30_34m$redmeat.new) #probability of DHA.EPA exp
probred_alt30_34m[1,] <- c(pr.no,pr.yes)

red30_34posm <- red30_34m[which(red30_34m!=0)]
red30_34posm <- as.vector(red30_34posm$redmeat.new)

fit1_alt_30_34m <- fitdist(red30_34posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red30_34posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_30_34m$estimate[1], sdlog=fit1_alt_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_30_34m, fitnames="lnorm") 

fit2_alt_30_34m <- fitdist(red30_34posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red30_34posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_30_34m$estimate[1], fit2_alt_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_30_34m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_30_34m$estimate[1], fit1_alt_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_30_34m$estimate[1], fit1_alt_30_34m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_30_34m$estimate[1], fit2_alt_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_30_34m$estimate[1], fit2_alt_30_34m$estimate[2])  #Anderson-Darling test


#Age 35-39, male
red35_39m <- subset(Alt_scenario, agegroups=="35-39" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt35_39m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red35_39m==0)/length(red35_39m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red35_39m!=0)/length(red35_39m$redmeat.new) #probability of DHA.EPA exp
probred_alt35_39m[1,] <- c(pr.no,pr.yes)

red35_39posm <- red35_39m[which(red35_39m!=0)]
red35_39posm <- as.vector(red35_39posm$redmeat.new)

fit1_alt_35_39m <- fitdist(red35_39posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red35_39posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_35_39m$estimate[1], sdlog=fit1_alt_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_35_39m, fitnames="lnorm") 

fit2_alt_35_39m <- fitdist(red35_39posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red35_39posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_35_39m$estimate[1], fit2_alt_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_35_39m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_35_39m$estimate[1], fit1_alt_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_35_39m$estimate[1], fit1_alt_35_39m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_35_39m$estimate[1], fit2_alt_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_35_39m$estimate[1], fit2_alt_35_39m$estimate[2])  #Anderson-Darling test


#Age 40-44, male
red40_44m <- subset(Alt_scenario, agegroups=="40-44" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt40_44m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red40_44m==0)/length(red40_44m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red40_44m!=0)/length(red40_44m$redmeat.new) #probability of DHA.EPA exp
probred_alt40_44m[1,] <- c(pr.no,pr.yes)

red40_44posm <- red40_44m[which(red40_44m!=0)]
red40_44posm <- as.vector(red40_44posm$redmeat.new)


fit1_alt_40_44m <- fitdist(red40_44posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red40_44posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_40_44m$estimate[1], sdlog=fit1_alt_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_40_44m, fitnames="lnorm") 

fit2_alt_40_44m <- fitdist(red40_44posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red40_44posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_40_44m$estimate[1], fit2_alt_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_40_44m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_40_44m$estimate[1], fit1_alt_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_40_44m$estimate[1], fit1_alt_40_44m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_40_44m$estimate[1], fit2_alt_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_40_44m$estimate[1], fit2_alt_40_44m$estimate[2])  #Anderson-Darling test


#Age 45-49, male
red45_49m <- subset(Alt_scenario, agegroups=="45-49" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt45_49m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red45_49m==0)/length(red45_49m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red45_49m!=0)/length(red45_49m$redmeat.new) #probability of DHA.EPA exp
probred_alt45_49m[1,] <- c(pr.no,pr.yes)

red45_49posm <- red45_49m[which(red45_49m!=0)]
red45_49posm <- as.vector(red45_49posm$redmeat.new)


fit1_alt_45_49m <- fitdist(red45_49posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red45_49posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_45_49m$estimate[1], sdlog=fit1_alt_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_45_49m, fitnames="lnorm") 

fit2_alt_45_49m <- fitdist(red45_49posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red45_49posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_45_49m$estimate[1], fit2_alt_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_45_49m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_45_49m$estimate[1], fit1_alt_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_45_49m$estimate[1], fit1_alt_45_49m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_45_49m$estimate[1], fit2_alt_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_45_49m$estimate[1], fit2_alt_45_49m$estimate[2])  #Anderson-Darling test


#Age 50-54, male
red50_54m <- subset(Alt_scenario, agegroups=="50-54" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt50_54m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red50_54m==0)/length(red50_54m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red50_54m!=0)/length(red50_54m$redmeat.new) #probability of DHA.EPA exp
probred_alt50_54m[1,] <- c(pr.no,pr.yes)

red50_54posm <- red50_54m[which(red50_54m!=0)]
red50_54posm <- as.vector(red50_54posm$redmeat.new)


fit1_alt_50_54m <- fitdist(red50_54posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red50_54posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_50_54m$estimate[1], sdlog=fit1_alt_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_50_54m, fitnames="lnorm") 

fit2_alt_50_54m <- fitdist(red50_54posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red50_54posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_50_54m$estimate[1], fit2_alt_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_50_54m, fitnames="gamma") 


cvm.test(t, plnorm, fit1_alt_50_54m$estimate[1], fit1_alt_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_50_54m$estimate[1], fit1_alt_50_54m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_50_54m$estimate[1], fit2_alt_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_50_54m$estimate[1], fit2_alt_50_54m$estimate[2])  #Anderson-Darling test


#Age 55-59, male
red55_59m <- subset(Alt_scenario, agegroups=="55-59" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt55_59m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red55_59m==0)/length(red55_59m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red55_59m!=0)/length(red55_59m$redmeat.new) #probability of DHA.EPA exp
probred_alt55_59m[1,] <- c(pr.no,pr.yes)

red55_59posm <- red55_59m[which(red55_59m!=0)]
red55_59posm <- as.vector(red55_59posm$redmeat.new)


fit1_alt_55_59m <- fitdist(red55_59posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red55_59posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_55_59m$estimate[1], sdlog=fit1_alt_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_55_59m, fitnames="lnorm") 

fit2_alt_55_59m <- fitdist(red55_59posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red55_59posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_55_59m$estimate[1], fit2_alt_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_55_59m, fitnames="gamma") 


cvm.test(t, plnorm, fit1_alt_55_59m$estimate[1], fit1_alt_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_55_59m$estimate[1], fit1_alt_55_59m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_55_59m$estimate[1], fit2_alt_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_55_59m$estimate[1], fit2_alt_55_59m$estimate[2])  #Anderson-Darling test


#Age 60-64, male
red60_64m <- subset(Alt_scenario, agegroups=="60-64" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt60_64m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red60_64m==0)/length(red60_64m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red60_64m!=0)/length(red60_64m$redmeat.new) #probability of DHA.EPA exp
probred_alt60_64m[1,] <- c(pr.no,pr.yes)

red60_64posm <- red60_64m[which(red60_64m!=0)]
red60_64posm <- as.vector(red60_64posm$redmeat.new)


fit1_alt_60_64m <- fitdist(red60_64posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red60_64posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_60_64m$estimate[1], sdlog=fit1_alt_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_60_64m, fitnames="lnorm") 

fit2_alt_60_64m <- fitdist(red60_64posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red60_64posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_60_64m$estimate[1], fit2_alt_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_60_64m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_60_64m$estimate[1], fit1_alt_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_60_64m$estimate[1], fit1_alt_60_64m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_60_64m$estimate[1], fit2_alt_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_60_64m$estimate[1], fit2_alt_60_64m$estimate[2])  #Anderson-Darling test


#Age 65-69, male
red65_69m <- subset(Alt_scenario, agegroups=="65-69" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt65_69m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red65_69m==0)/length(red65_69m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red65_69m!=0)/length(red65_69m$redmeat.new) #probability of DHA.EPA exp
probred_alt65_69m[1,] <- c(pr.no,pr.yes)

red65_69posm <- red65_69m[which(red65_69m!=0)]
red65_69posm <- as.vector(red65_69posm$redmeat.new)

fit1_alt_65_69m <- fitdist(red65_69posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red65_69posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_65_69m$estimate[1], sdlog=fit1_alt_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_65_69m, fitnames="lnorm") 

fit2_alt_65_69m <- fitdist(red65_69posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red65_69posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_65_69m$estimate[1], fit2_alt_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_65_69m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_65_69m$estimate[1], fit1_alt_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_65_69m$estimate[1], fit1_alt_65_69m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_65_69m$estimate[1], fit2_alt_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_65_69m$estimate[1], fit2_alt_65_69m$estimate[2])  #Anderson-Darling test


#Age 70-74, male
red70_74m <- subset(Alt_scenario, agegroups=="70-74" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt70_74m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red70_74m==0)/length(red70_74m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red70_74m!=0)/length(red70_74m$redmeat.new) #probability of DHA.EPA exp
probred_alt70_74m[1,] <- c(pr.no,pr.yes)

red70_74posm <- red70_74m[which(red70_74m!=0)]
red70_74posm <- as.vector(red70_74posm$redmeat.new)

fit1_alt_70_74m <- fitdist(red70_74posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red70_74posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_70_74m$estimate[1], sdlog=fit1_alt_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_70_74m, fitnames="lnorm") 

fit2_alt_70_74m <- fitdist(red70_74posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red70_74posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_70_74m$estimate[1], fit2_alt_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_70_74m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_70_74m$estimate[1], fit1_alt_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_70_74m$estimate[1], fit1_alt_70_74m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_70_74m$estimate[1], fit2_alt_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_70_74m$estimate[1], fit2_alt_70_74m$estimate[2])  #Anderson-Darling test


#Age 75-79, male
red75_79m <- subset(Alt_scenario, agegroups=="75-79" & sex =="1", select = redmeat.new)

#Probability of intake
probred_alt75_79m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red75_79m==0)/length(red75_79m$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red75_79m!=0)/length(red75_79m$redmeat.new) #probability of DHA.EPA exp
probred_alt75_79m[1,] <- c(pr.no,pr.yes)

red75_79posm <- red75_79m[which(red75_79m!=0)]
red75_79posm <- as.vector(red75_79posm$redmeat.new)

fit1_alt_75_79m <- fitdist(red75_79posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red75_79posm
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_75_79m$estimate[1], sdlog=fit1_alt_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_75_79m, fitnames="lnorm") 

fit2_alt_75_79m <- fitdist(red75_79posm, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red75_79posm
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=40)
lines(x, pgamma(x, fit2_alt_75_79m$estimate[1], fit2_alt_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_75_79m, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_75_79m$estimate[1], fit1_alt_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_75_79m$estimate[1], fit1_alt_75_79m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_75_79m$estimate[1], fit2_alt_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_75_79m$estimate[1], fit2_alt_75_79m$estimate[2])  #Anderson-Darling test


##### Fitting exposures females #####


#Age 15-19, female
red15_19w <- subset(Alt_scenario, agegroups=="15-19" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt15_19w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red15_19w==0)/length(red15_19w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red15_19w!=0)/length(red15_19w$redmeat.new) #probability of DHA.EPA exp
probred_alt15_19w[1,] <- c(pr.no,pr.yes)

red15_19posw <- red15_19w[which(red15_19w!=0)]
red15_19posw <- as.vector(red15_19posw$redmeat.new)

fit1_alt_15_19w <- fitdist(red15_19posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red15_19posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_15_19w$estimate[1], sdlog=fit1_alt_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_15_19w, fitnames="lnorm") 

fit2_alt_15_19w <- fitdist(red15_19posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red15_19posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_15_19w$estimate[1], fit2_alt_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_15_19w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_15_19w$estimate[1], fit1_alt_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_15_19w$estimate[1], fit1_alt_15_19w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_15_19w$estimate[1], fit2_alt_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_15_19w$estimate[1], fit2_alt_15_19w$estimate[2])  #Anderson-Darling test


#Age 20-24, female
red20_24w <- subset(Alt_scenario, agegroups=="20-24" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt20_24w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red20_24w==0)/length(red20_24w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red20_24w!=0)/length(red20_24w$redmeat.new) #probability of DHA.EPA exp
probred_alt20_24w[1,] <- c(pr.no,pr.yes)

red20_24posw <- red20_24w[which(red20_24w!=0)]
red20_24posw <- as.vector(red20_24posw$redmeat.new)

fit1_alt_20_24w <- fitdist(red20_24posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red20_24posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_20_24w$estimate[1], sdlog=fit1_alt_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_20_24w, fitnames="lnorm") 

fit2_alt_20_24w <- fitdist(red20_24posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red20_24posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_20_24w$estimate[1], fit2_alt_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_20_24w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_20_24w$estimate[1], fit1_alt_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_20_24w$estimate[1], fit1_alt_20_24w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_20_24w$estimate[1], fit2_alt_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_20_24w$estimate[1], fit2_alt_20_24w$estimate[2])  #Anderson-Darling test


#Age 25-29, female
red25_29w <- subset(Alt_scenario, agegroups=="25-29" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt25_29w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red25_29w==0)/length(red25_29w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red25_29w!=0)/length(red25_29w$redmeat.new) #probability of DHA.EPA exp
probred_alt25_29w[1,] <- c(pr.no,pr.yes)

red25_29posw <- red25_29w[which(red25_29w!=0)]
red25_29posw <- as.vector(red25_29posw$redmeat.new)

fit1_alt_25_29w <- fitdist(red25_29posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red25_29posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_25_29w$estimate[1], sdlog=fit1_alt_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_25_29w, fitnames="lnorm") 

fit2_alt_25_29w <- fitdist(red25_29posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red25_29posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_25_29w$estimate[1], fit2_alt_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_25_29w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_25_29w$estimate[1], fit1_alt_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_25_29w$estimate[1], fit1_alt_25_29w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_25_29w$estimate[1], fit2_alt_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_25_29w$estimate[1], fit2_alt_25_29w$estimate[2])  #Anderson-Darling test


#Age 30-34, female
red30_34w <- subset(Alt_scenario, agegroups=="30-34" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt30_34w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red30_34w==0)/length(red30_34w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red30_34w!=0)/length(red30_34w$redmeat.new) #probability of DHA.EPA exp
probred_alt30_34w[1,] <- c(pr.no,pr.yes)

red30_34posw <- red30_34w[which(red30_34w!=0)]
red30_34posw <- as.vector(red30_34posw$redmeat.new)

fit1_alt_30_34w <- fitdist(red30_34posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red30_34posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_30_34w$estimate[1], sdlog=fit1_alt_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_30_34w, fitnames="lnorm") 

fit2_alt_30_34w <- fitdist(red30_34posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red30_34posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_30_34w$estimate[1], fit2_alt_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_30_34w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_30_34w$estimate[1], fit1_alt_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_30_34w$estimate[1], fit1_alt_30_34w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_30_34w$estimate[1], fit2_alt_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_30_34w$estimate[1], fit2_alt_30_34w$estimate[2])  #Anderson-Darling test


#Age 35-39, female
red35_39w <- subset(Alt_scenario, agegroups=="35-39" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt35_39w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red35_39w==0)/length(red35_39w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red35_39w!=0)/length(red35_39w$redmeat.new) #probability of DHA.EPA exp
probred_alt35_39w[1,] <- c(pr.no,pr.yes)

red35_39posw <- red35_39w[which(red35_39w!=0)]
red35_39posw <- as.vector(red35_39posw$redmeat.new)

fit1_alt_35_39w <- fitdist(red35_39posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red35_39posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_35_39w$estimate[1], sdlog=fit1_alt_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_35_39w, fitnames="lnorm") 

fit2_alt_35_39w <- fitdist(red35_39posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red35_39posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_35_39w$estimate[1], fit2_alt_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_35_39w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_35_39w$estimate[1], fit1_alt_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_35_39w$estimate[1], fit1_alt_35_39w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_35_39w$estimate[1], fit2_alt_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_35_39w$estimate[1], fit2_alt_35_39w$estimate[2])  #Anderson-Darling test


#Age 40-44, female
red40_44w <- subset(Alt_scenario, agegroups=="40-44" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt40_44w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red40_44w==0)/length(red40_44w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red40_44w!=0)/length(red40_44w$redmeat.new) #probability of DHA.EPA exp
probred_alt40_44w[1,] <- c(pr.no,pr.yes)

red40_44posw <- red40_44w[which(red40_44w!=0)]
red40_44posw <- as.vector(red40_44posw$redmeat.new)


fit1_alt_40_44w <- fitdist(red40_44posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red40_44posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_40_44w$estimate[1], sdlog=fit1_alt_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_40_44w, fitnames="lnorm") 

fit2_alt_40_44w <- fitdist(red40_44posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red40_44posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_40_44w$estimate[1], fit2_alt_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_40_44w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_40_44w$estimate[1], fit1_alt_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_40_44w$estimate[1], fit1_alt_40_44w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_40_44w$estimate[1], fit2_alt_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_40_44w$estimate[1], fit2_alt_40_44w$estimate[2])  #Anderson-Darling test


#Age 45-49, female
red45_49w <- subset(Alt_scenario, agegroups=="45-49" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt45_49w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red45_49w==0)/length(red45_49w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red45_49w!=0)/length(red45_49w$redmeat.new) #probability of DHA.EPA exp
probred_alt45_49w[1,] <- c(pr.no,pr.yes)

red45_49posw <- red45_49w[which(red45_49w!=0)]
red45_49posw <- as.vector(red45_49posw$redmeat.new)


fit1_alt_45_49w <- fitdist(red45_49posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red45_49posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_45_49w$estimate[1], sdlog=fit1_alt_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_45_49w, fitnames="lnorm") 

fit2_alt_45_49w <- fitdist(red45_49posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red45_49posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_45_49w$estimate[1], fit2_alt_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_45_49w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_45_49w$estimate[1], fit1_alt_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_45_49w$estimate[1], fit1_alt_45_49w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_45_49w$estimate[1], fit2_alt_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_45_49w$estimate[1], fit2_alt_45_49w$estimate[2])  #Anderson-Darling test


#Age 50-54, female
red50_54w <- subset(Alt_scenario, agegroups=="50-54" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt50_54w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red50_54w==0)/length(red50_54w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red50_54w!=0)/length(red50_54w$redmeat.new) #probability of DHA.EPA exp
probred_alt50_54w[1,] <- c(pr.no,pr.yes)

red50_54posw <- red50_54w[which(red50_54w!=0)]
red50_54posw <- as.vector(red50_54posw$redmeat.new)


fit1_alt_50_54w <- fitdist(red50_54posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red50_54posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_50_54w$estimate[1], sdlog=fit1_alt_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_50_54w, fitnames="lnorm") 

fit2_alt_50_54w <- fitdist(red50_54posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red50_54posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_50_54w$estimate[1], fit2_alt_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_50_54w, fitnames="gamma") 


cvm.test(t, plnorm, fit1_alt_50_54w$estimate[1], fit1_alt_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_50_54w$estimate[1], fit1_alt_50_54w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_50_54w$estimate[1], fit2_alt_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_50_54w$estimate[1], fit2_alt_50_54w$estimate[2])  #Anderson-Darling test


#Age 55-59, female
red55_59w <- subset(Alt_scenario, agegroups=="55-59" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt55_59w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red55_59w==0)/length(red55_59w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red55_59w!=0)/length(red55_59w$redmeat.new) #probability of DHA.EPA exp
probred_alt55_59w[1,] <- c(pr.no,pr.yes)

red55_59posw <- red55_59w[which(red55_59w!=0)]
red55_59posw <- as.vector(red55_59posw$redmeat.new)


fit1_alt_55_59w <- fitdist(red55_59posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red55_59posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_55_59w$estimate[1], sdlog=fit1_alt_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_55_59w, fitnames="lnorm") 

fit2_alt_55_59w <- fitdist(red55_59posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red55_59posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_55_59w$estimate[1], fit2_alt_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_55_59w, fitnames="gamma") 


cvm.test(t, plnorm, fit1_alt_55_59w$estimate[1], fit1_alt_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_55_59w$estimate[1], fit1_alt_55_59w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_55_59w$estimate[1], fit2_alt_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_55_59w$estimate[1], fit2_alt_55_59w$estimate[2])  #Anderson-Darling test


#Age 60-64, female
red60_64w <- subset(Alt_scenario, agegroups=="60-64" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt60_64w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red60_64w==0)/length(red60_64w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red60_64w!=0)/length(red60_64w$redmeat.new) #probability of DHA.EPA exp
probred_alt60_64w[1,] <- c(pr.no,pr.yes)

red60_64posw <- red60_64w[which(red60_64w!=0)]
red60_64posw <- as.vector(red60_64posw$redmeat.new)


fit1_alt_60_64w <- fitdist(red60_64posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red60_64posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_60_64w$estimate[1], sdlog=fit1_alt_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_60_64w, fitnames="lnorm") 

fit2_alt_60_64w <- fitdist(red60_64posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red60_64posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_60_64w$estimate[1], fit2_alt_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_60_64w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_60_64w$estimate[1], fit1_alt_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_60_64w$estimate[1], fit1_alt_60_64w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_60_64w$estimate[1], fit2_alt_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_60_64w$estimate[1], fit2_alt_60_64w$estimate[2])  #Anderson-Darling test


#Age 65-69, female
red65_69w <- subset(Alt_scenario, agegroups=="65-69" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt65_69w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red65_69w==0)/length(red65_69w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red65_69w!=0)/length(red65_69w$redmeat.new) #probability of DHA.EPA exp
probred_alt65_69w[1,] <- c(pr.no,pr.yes)

red65_69posw <- red65_69w[which(red65_69w!=0)]
red65_69posw <- as.vector(red65_69posw$redmeat.new)

fit1_alt_65_69w <- fitdist(red65_69posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red65_69posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_65_69w$estimate[1], sdlog=fit1_alt_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_65_69w, fitnames="lnorm") 

fit2_alt_65_69w <- fitdist(red65_69posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red65_69posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_65_69w$estimate[1], fit2_alt_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_65_69w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_65_69w$estimate[1], fit1_alt_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_65_69w$estimate[1], fit1_alt_65_69w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_65_69w$estimate[1], fit2_alt_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_65_69w$estimate[1], fit2_alt_65_69w$estimate[2])  #Anderson-Darling test


#Age 70-74, female
red70_74w <- subset(Alt_scenario, agegroups=="70-74" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt70_74w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red70_74w==0)/length(red70_74w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red70_74w!=0)/length(red70_74w$redmeat.new) #probability of DHA.EPA exp
probred_alt70_74w[1,] <- c(pr.no,pr.yes)

red70_74posw <- red70_74w[which(red70_74w!=0)]
red70_74posw <- as.vector(red70_74posw$redmeat.new)

fit1_alt_70_74w <- fitdist(red70_74posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red70_74posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_70_74w$estimate[1], sdlog=fit1_alt_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_70_74w, fitnames="lnorm") 

fit2_alt_70_74w <- fitdist(red70_74posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red70_74posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_70_74w$estimate[1], fit2_alt_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_70_74w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_70_74w$estimate[1], fit1_alt_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_70_74w$estimate[1], fit1_alt_70_74w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_70_74w$estimate[1], fit2_alt_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_70_74w$estimate[1], fit2_alt_70_74w$estimate[2])  #Anderson-Darling test


#Age 75-79, female
red75_79w <- subset(Alt_scenario, agegroups=="75-79" & sex =="2", select = redmeat.new)

#Probability of intake
probred_alt75_79w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(red75_79w==0)/length(red75_79w$redmeat.new) #probability of zero DHA.EPA exp
pr.yes <- sum(red75_79w!=0)/length(red75_79w$redmeat.new) #probability of DHA.EPA exp
probred_alt75_79w[1,] <- c(pr.no,pr.yes)

red75_79posw <- red75_79w[which(red75_79w!=0)]
red75_79posw <- as.vector(red75_79posw$redmeat.new)

fit1_alt_75_79w <- fitdist(red75_79posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- red75_79posw
plot(ecdf(t), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, plnorm(x, meanlog=fit1_alt_75_79w$estimate[1], sdlog=fit1_alt_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit1_alt_75_79w, fitnames="lnorm") 

fit2_alt_75_79w <- fitdist(red75_79posw, 'gamma') #Fit gamma distribution to intake amounts
t2 <- red75_79posw
plot(ecdf(t2), lty=1)
x <- seq(0, 4000, length=1000)
lines(x, pgamma(x, fit2_alt_75_79w$estimate[1], fit2_alt_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit2_alt_75_79w, fitnames="gamma") 

cvm.test(t, plnorm, fit1_alt_75_79w$estimate[1], fit1_alt_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_alt_75_79w$estimate[1], fit1_alt_75_79w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_alt_75_79w$estimate[1], fit2_alt_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_alt_75_79w$estimate[1], fit2_alt_75_79w$estimate[2])  #Anderson-Darling test



##### DALY calc males #####

max(MeatDaily$red.meat)
# [1] 332.7429

# Gamma distribution has the best fit for both the reference and alternative scenario
# Truncate at 400 g red meat/day to avoid extreme values. Truncation did not change mean.


# Age 15-19, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred15_19m[,1],probred15_19m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_15_19m$estimate[1], fit2_15_19m$estimate[2], rtrunc = T, linf = 0, lsup = 400)
RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt15_19m[,1],probred_alt15_19m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_15_19m$estimate[1], fit2_alt_15_19m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p15_19m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_15_19m <- pEffect * pop_15_19m


YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal15_19m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal15_19m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal15_19m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal15_19m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal15_19m) * p_stoma * dw_stoma * (LE_15_19m - (17 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal15_19m * (SEYLL_15_19 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY15_19totalm_alt <- DALY * pop_15_19m
summary(DALY15_19totalm_alt)

summary(DALY15_19totalm)

# Age 20-24, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred20_24m[,1],probred20_24m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_20_24m$estimate[1], fit2_20_24m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt20_24m[,1],probred_alt20_24m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_20_24m$estimate[1], fit2_alt_20_24m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p20_24m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_20_24m <- pEffect * pop_20_24m


YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal20_24m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal20_24m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal20_24m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal20_24m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal20_24m) * p_stoma * dw_stoma * (LE_20_24m - (22 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal20_24m * (SEYLL_20_24 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY20_24totalm_alt <- DALY * pop_20_24m
summary(DALY20_24totalm_alt)

summary(DALY20_24totalm)


# Age 25-29, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred25_29m[,1],probred25_29m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_25_29m$estimate[1], fit2_25_29m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt25_29m[,1],probred_alt25_29m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_25_29m$estimate[1], fit2_alt_25_29m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p25_29m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_25_29m <- pEffect * pop_25_29m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal25_29m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal25_29m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal25_29m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal25_29m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal25_29m) * p_stoma * dw_stoma * (LE_25_29m - (27 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal25_29m * (SEYLL_25_29 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY25_29totalm_alt <- DALY * pop_25_29m
summary(DALY25_29totalm_alt)

summary(DALY25_29totalm)


# Age 30-34, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred30_34m[,1],probred30_34m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_30_34m$estimate[1], fit2_30_34m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt30_34m[,1],probred_alt30_34m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_30_34m$estimate[1], fit2_alt_30_34m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p30_34m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_30_34m <- pEffect * pop_30_34m


YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal30_34m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal30_34m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal30_34m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal30_34m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal30_34m) * p_stoma * dw_stoma * (LE_30_34m - (32 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal30_34m * (SEYLL_30_34 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY30_34totalm_alt <- DALY * pop_30_34m
summary(DALY30_34totalm_alt)


summary(DALY30_34totalm)


# Age 35-39, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred35_39m[,1],probred35_39m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_35_39m$estimate[1], fit2_35_39m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt35_39m[,1],probred_alt35_39m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_35_39m$estimate[1], fit2_alt_35_39m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p35_39m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_35_39m <- pEffect * pop_35_39m


YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal35_39m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal35_39m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal35_39m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal35_39m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal35_39m) * p_stoma * dw_stoma * (LE_35_39m - (37 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal35_39m * (SEYLL_35_39 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY35_39totalm_alt <- DALY * pop_35_39m
summary(DALY35_39totalm_alt)

summary(DALY35_39totalm)


# Age 40-44, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred40_44m[,1],probred40_44m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_40_44m$estimate[1], fit2_40_44m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt40_44m[,1],probred_alt40_44m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_40_44m$estimate[1], fit2_alt_40_44m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p40_44m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_40_44m <- pEffect * pop_40_44m


YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal40_44m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal40_44m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal40_44m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal40_44m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal40_44m) * p_stoma * dw_stoma * (LE_40_44m - (42 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal40_44m * (SEYLL_40_44 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY40_44totalm_alt <- DALY * pop_40_44m
summary(DALY40_44totalm_alt)

summary(DALY40_44totalm)


# Age 45-49, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred45_49m[,1],probred45_49m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_45_49m$estimate[1], fit2_45_49m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt45_49m[,1],probred_alt45_49m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_45_49m$estimate[1], fit2_alt_45_49m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p45_49m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_45_49m <- pEffect * pop_45_49m


YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal45_49m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal45_49m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal45_49m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal45_49m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal45_49m) * p_stoma * dw_stoma * (LE_45_49m - (47 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal45_49m * (SEYLL_45_49 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY45_49totalm_alt <- DALY * pop_45_49m
summary(DALY45_49totalm_alt)

summary(DALY45_49totalm)



# Age 50-54, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred50_54m[,1],probred50_54m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_50_54m$estimate[1], fit2_50_54m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt50_54m[,1],probred_alt50_54m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_50_54m$estimate[1], fit2_alt_50_54m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p50_54m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_50_54m <- pEffect * pop_50_54m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal50_54m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal50_54m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal50_54m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal50_54m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal50_54m) * p_stoma * dw_stoma * (LE_50_54m - (52 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal50_54m * (SEYLL_50_54 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY50_54totalm_alt <- DALY * pop_50_54m
summary(DALY50_54totalm_alt)

summary(DALY50_54totalm)

# Age 55-59, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred55_59m[,1],probred55_59m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_55_59m$estimate[1], fit2_55_59m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt55_59m[,1],probred_alt55_59m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_55_59m$estimate[1], fit2_alt_55_59m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p55_59m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_55_59m <- pEffect * pop_55_59m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal55_59m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal55_59m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal55_59m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal55_59m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal55_59m) * p_stoma * dw_stoma * (LE_55_59m - (57 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal55_59m * (SEYLL_55_59 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY55_59totalm_alt <- DALY * pop_55_59m
summary(DALY55_59totalm_alt)

summary(DALY55_59totalm)

# Age 60-64, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred60_64m[,1],probred60_64m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_60_64m$estimate[1], fit2_60_64m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt60_64m[,1],probred_alt60_64m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_60_64m$estimate[1], fit2_alt_60_64m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p60_64m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_60_64m <- pEffect * pop_60_64m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal60_64m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal60_64m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal60_64m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal60_64m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal60_64m) * p_stoma * dw_stoma * (LE_60_64m - (62 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal60_64m * (SEYLL_60_64 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY60_64totalm_alt <- DALY * pop_60_64m
summary(DALY60_64totalm_alt)

summary(DALY60_64totalm)


# Age 65-69, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred65_69m[,1],probred65_69m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_65_69m$estimate[1], fit2_65_69m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt65_69m[,1],probred_alt65_69m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_65_69m$estimate[1], fit2_alt_65_69m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p65_69m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_65_69m <- pEffect * pop_65_69m


YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal65_69m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal65_69m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal65_69m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal65_69m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal65_69m) * p_stoma * dw_stoma * (LE_65_69m - (67 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal65_69m * (SEYLL_65_69 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY65_69totalm_alt <- DALY * pop_65_69m
summary(DALY65_69totalm_alt)

summary(DALY65_69totalm)


# Age 70-74, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred70_74m[,1],probred70_74m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_70_74m$estimate[1], fit2_70_74m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt70_74m[,1],probred_alt70_74m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_70_74m$estimate[1], fit2_alt_70_74m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p70_74m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_70_74m <- pEffect * pop_70_74m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal70_74m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal70_74m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal70_74m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal70_74m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal70_74m) * p_stoma * dw_stoma * (LE_70_74m - (72 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal70_74m * (SEYLL_70_74 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY70_74totalm_alt <- DALY * pop_70_74m
summary(DALY70_74totalm_alt)

summary(DALY70_74totalm)

# Age 75-79, males

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred75_79m[,1],probred75_79m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_75_79m$estimate[1], fit2_75_79m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt75_79m[,1],probred_alt75_79m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_75_79m$estimate[1], fit2_alt_75_79m$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p75_79m/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_75_79m <- pEffect * pop_75_79m

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal75_79m) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal75_79m * dw_remission * time_remission_death + #remission_death
  pEffect * fatal75_79m * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal75_79m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal75_79m) * p_stoma * dw_stoma * (LE_75_79m - (77 + time_total_YLD)) #stoma

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
  pEffect * fatal80_84m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal80_84m) * p_stoma * dw_stoma * (LE_80_84m - (82 + time_total_YLD)) #stoma

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
  pEffect * fatal85m * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal85m) * p_stoma * dw_stoma * (LE_85m - (87 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal85m * (SEYLL_85 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY85totalm_alt <- DALY * pop_85m
summary(DALY85totalm_alt)

summary(DALY85totalm)

################################################ DALY calc females #####################################################

# Age 15-19, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred15_19w[,1],probred15_19w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_15_19w$estimate[1], fit2_15_19w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt15_19w[,1],probred_alt15_19w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_15_19w$estimate[1], fit2_alt_15_19w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p15_19w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_15_19w <- pEffect * pop_15_19w


YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal15_19w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal15_19w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal15_19w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal15_19w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal15_19w) * p_stoma * dw_stoma * (LE_15_19w - (17 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal15_19w * (SEYLL_15_19 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY15_19totalw_alt <- DALY * pop_15_19w
summary(DALY15_19totalw_alt)

summary(DALY15_19totalw)

# Age 20-24, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred20_24w[,1],probred20_24w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_20_24w$estimate[1], fit2_20_24w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt20_24w[,1],probred_alt20_24w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_20_24w$estimate[1], fit2_alt_20_24w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p20_24w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_20_24w <- pEffect * pop_20_24w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal20_24w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal20_24w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal20_24w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal20_24w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal20_24w) * p_stoma * dw_stoma * (LE_20_24w - (22 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal20_24w * (SEYLL_20_24 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY20_24totalw_alt <- DALY * pop_20_24w
summary(DALY20_24totalw_alt)

summary(DALY20_24totalw)

# Age 25-29, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred25_29w[,1],probred25_29w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_25_29w$estimate[1], fit2_25_29w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt25_29w[,1],probred_alt25_29w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_25_29w$estimate[1], fit2_alt_25_29w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p25_29w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_25_29w <- pEffect * pop_25_29w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal25_29w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal25_29w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal25_29w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal25_29w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal25_29w) * p_stoma * dw_stoma * (LE_25_29w - (27 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal25_29w * (SEYLL_25_29 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY25_29totalw_alt <- DALY * pop_25_29w
summary(DALY25_29totalw_alt)

summary(DALY25_29totalw)


# Age 30-34, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred30_34w[,1],probred30_34w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_30_34w$estimate[1], fit2_30_34w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt30_34w[,1],probred_alt30_34w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_30_34w$estimate[1], fit2_alt_30_34w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p30_34w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_30_34w <- pEffect * pop_30_34w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal30_34w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal30_34w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal30_34w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal30_34w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal30_34w) * p_stoma * dw_stoma * (LE_30_34w - (32 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal30_34w * (SEYLL_30_34 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY30_34totalw_alt <- DALY * pop_30_34w
summary(DALY30_34totalw_alt)

summary(DALY30_34totalw)


# Age 35-39, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred35_39w[,1],probred35_39w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_35_39w$estimate[1], fit2_35_39w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt35_39w[,1],probred_alt35_39w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_35_39w$estimate[1], fit2_alt_35_39w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p35_39w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_35_39w <- pEffect * pop_35_39w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal35_39w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal35_39w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal35_39w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal35_39w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal35_39w) * p_stoma * dw_stoma * (LE_35_39w - (37 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal35_39w * (SEYLL_35_39 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY35_39totalw_alt <- DALY * pop_35_39w
summary(DALY35_39totalw_alt)

summary(DALY35_39totalw)


# Age 40-44, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred40_44w[,1],probred40_44w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_40_44w$estimate[1], fit2_40_44w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt40_44w[,1],probred_alt40_44w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_40_44w$estimate[1], fit2_alt_40_44w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p40_44w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_40_44w <- pEffect * pop_40_44w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal40_44w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal40_44w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal40_44w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal40_44w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal40_44w) * p_stoma * dw_stoma * (LE_40_44w - (42 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal40_44w * (SEYLL_40_44 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY40_44totalw_alt <- DALY * pop_40_44w
summary(DALY40_44totalw_alt)

summary(DALY40_44totalw)


# Age 45-49, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred45_49w[,1],probred45_49w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_45_49w$estimate[1], fit2_45_49w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt45_49w[,1],probred_alt45_49w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_45_49w$estimate[1], fit2_alt_45_49w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)

pEffect <- p45_49w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_45_49w <- pEffect * pop_45_49w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal45_49w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal45_49w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal45_49w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal45_49w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal45_49w) * p_stoma * dw_stoma * (LE_45_49w - (47 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal45_49w * (SEYLL_45_49 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY45_49totalw_alt <- DALY * pop_45_49w
summary(DALY45_49totalw_alt)


summary(DALY45_49totalw)



# Age 50-54, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred50_54w[,1],probred50_54w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_50_54w$estimate[1], fit2_50_54w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt50_54w[,1],probred_alt50_54w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_50_54w$estimate[1], fit2_alt_50_54w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p50_54w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_50_54w <- pEffect * pop_50_54w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal50_54w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal50_54w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal50_54w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal50_54w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal50_54w) * p_stoma * dw_stoma * (LE_50_54w - (52 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal50_54w * (SEYLL_50_54 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY50_54totalw_alt <- DALY * pop_50_54w
summary(DALY50_54totalw_alt)


summary(DALY50_54totalw)



# Age 55-59, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred55_59w[,1],probred55_59w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_55_59w$estimate[1], fit2_55_59w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt55_59w[,1],probred_alt55_59w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_55_59w$estimate[1], fit2_alt_55_59w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p55_59w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_55_59w <- pEffect * pop_55_59w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal55_59w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal55_59w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal55_59w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal55_59w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal55_59w) * p_stoma * dw_stoma * (LE_55_59w - (57 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal55_59w * (SEYLL_55_59 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY55_59totalw_alt <- DALY * pop_55_59w
summary(DALY55_59totalw_alt)


summary(DALY55_59totalw)



# Age 60-64, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred60_64w[,1],probred60_64w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_60_64w$estimate[1], fit2_60_64w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt60_64w[,1],probred_alt60_64w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_60_64w$estimate[1], fit2_alt_60_64w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p60_64w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_60_64w <- pEffect * pop_60_64w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal60_64w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal60_64w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal60_64w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal60_64w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal60_64w) * p_stoma * dw_stoma * (LE_60_64w - (62 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal60_64w * (SEYLL_60_64 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY60_64totalw_alt <- DALY * pop_60_64w
summary(DALY60_64totalw_alt)


summary(DALY60_64totalw)


# Age 65-69, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred65_69w[,1],probred65_69w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_65_69w$estimate[1], fit2_65_69w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt65_69w[,1],probred_alt65_69w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_65_69w$estimate[1], fit2_alt_65_69w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p65_69w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_65_69w <- pEffect * pop_65_69w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal65_69w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal65_69w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal65_69w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal65_69w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal65_69w) * p_stoma * dw_stoma * (LE_65_69w - (67 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal65_69w * (SEYLL_65_69 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY65_69totalw_alt <- DALY * pop_65_69w
summary(DALY65_69totalw_alt)


summary(DALY65_69totalw)




# Age 70-74, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred70_74w[,1],probred70_74w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_70_74w$estimate[1], fit2_70_74w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt70_74w[,1],probred_alt70_74w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_70_74w$estimate[1], fit2_alt_70_74w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p70_74w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_70_74w <- pEffect * pop_70_74w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal70_74w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal70_74w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal70_74w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal70_74w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal70_74w) * p_stoma * dw_stoma * (LE_70_74w - (72 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal70_74w * (SEYLL_70_74 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY70_74totalw_alt <- DALY * pop_70_74w
summary(DALY70_74totalw_alt)


summary(DALY70_74totalw)



# Age 75-79, females

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred75_79w[,1],probred75_79w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_75_79w$estimate[1], fit2_75_79w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_ref <- ifexp * exp * r + 1

summary(RR_ref)


set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probred_alt75_79w[,1],probred_alt75_79w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_alt_75_79w$estimate[1], fit2_alt_75_79w$estimate[2], rtrunc = T, linf = 0, lsup = 400)

RR_alt <- ifexp * exp * r + 1

summary(RR_alt)


pEffect <- p75_79w/RR_ref * RR_alt  # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)

summary(pEffect)

cases_75_79w <- pEffect * pop_75_79w

YLD <- pEffect * dw_diagnosis * time_diagnosis + #diagnosis
  pEffect * (1-fatal75_79w) * dw_remission * time_remission_cure + #remission_cure
  pEffect * fatal75_79w * dw_remission * time_remission_death + #remission_death
  pEffect * fatal75_79w * dw_disseminated * time_disseminated + #disseminated
  pEffect * fatal75_79w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal75_79w) * p_stoma * dw_stoma * (LE_75_79w - (72 + time_total_YLD)) #stoma

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
  pEffect * fatal80_84w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal80_84w) * p_stoma * dw_stoma * (LE_80_84w - (82 + time_total_YLD)) #stoma

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
  pEffect * fatal85w * dw_terminal * time_terminal + #terminal
  pEffect * (1-fatal85w) * p_stoma * dw_stoma * (LE_85w - (87 + time_total_YLD)) #stoma

#YLLs
YLL <- pEffect * fatal85w * (SEYLL_85 -  time_total_YLL) #SEYLL for age of onset

#DALYs
DALY <- YLL + YLD

DALY85totalw_alt <- DALY * pop_85w
summary(DALY85totalw_alt)


summary(DALY85totalw)


###### Total DALYs Alt #####

DALYtotal_alt <- DALY15_19totalm_alt + DALY20_24totalm_alt + DALY25_29totalm_alt + DALY30_34totalm_alt + DALY35_39totalm_alt +
  DALY40_44totalm_alt + DALY45_49totalm_alt + DALY50_54totalm_alt + DALY55_59totalm_alt + DALY60_64totalm_alt + DALY65_69totalm_alt +
  DALY70_74totalm_alt + DALY75_79totalm_alt + DALY80_84totalm_alt + DALY85totalm_alt + 
  DALY15_19totalw_alt + DALY20_24totalw_alt + DALY25_29totalw_alt + DALY30_34totalw_alt + DALY35_39totalw_alt +
  DALY40_44totalw_alt + DALY45_49totalw_alt + DALY50_54totalw_alt + DALY55_59totalw_alt + DALY60_64totalw_alt + DALY65_69totalw_alt +
  DALY70_74totalw_alt + DALY75_79totalw_alt + DALY80_84totalw_alt + DALY85totalw_alt

summary(DALYtotal_alt)


##### Total cases alternative scenario #####

cases_alt <- cases_15_19m + cases_20_24m + cases_25_29m + cases_30_34m + cases_35_39m + cases_40_44m + cases_45_49m + 
  cases_50_54m + cases_55_59m + cases_60_64m + cases_65_69m + cases_70_74m + cases_75_79m + cases_80_84m + cases_85m +
  cases_15_19w + cases_20_24w + cases_25_29w + cases_30_34w + cases_35_39w + cases_40_44w + cases_45_49w + 
  cases_50_54w + cases_55_59w + cases_60_64w + cases_65_69w + cases_70_74w + cases_75_79w + cases_80_84w + cases_85w

summary(cases_alt)


##### Extra number of cases #####

cases_diff <- cases_alt - cases_ref

summary(cases_diff)


##### Per 100,000 #####

pop15_85_plus <- pop_men + pop_women

tDALY_ref_red_CRC <- DALYtotalref
summary(tDALY_ref_red_CRC)


tDALY_alt_red_CRC <-DALYtotal_alt
summary(tDALY_alt_red_CRC)


tDALY_ref_red_CRC_100000 <- tDALY_ref_red_CRC/pop15_85_plus*1e+05
summary(tDALY_ref_red_CRC_100000)


tDALY_alt_red_CRC_100000 <- tDALY_alt_red_CRC/pop15_85_plus*1e+05
summary(tDALY_alt_red_CRC_100000)


dDALY_alt_red_CRC <- tDALY_alt_red_CRC - tDALY_ref_red_CRC
summary(dDALY_alt_red_CRC)

dDALY_alt_red_CRC_100000 <- dDALY_alt_red_CRC/pop15_85_plus*1e+05
summary(dDALY_alt_red_CRC_100000)


##### Save #####

tDALY_ref_red_CRC_unc <- apply(tDALY_ref_red_CRC, 2, function(x) x*1) #get data out of mc2d
write.csv(tDALY_ref_red_CRC_unc, "tDALY_ref_red_CRC_unc.csv")


tDALY_alt_red_CRC_matrix <- apply(tDALY_alt_red_CRC, 2, function(x) x*1) #get data out of mc2d
tDALY_alt_red_CRC_unc <- apply(tDALY_alt_red_CRC_matrix, 2, function(x) mean(x)) #only uncertainty dimension
write.csv(tDALY_alt_red_CRC_unc, "tDALY_alt_red_CRC_unc.csv")

