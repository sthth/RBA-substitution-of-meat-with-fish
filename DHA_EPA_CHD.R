memory.limit(56000)

#packages
library(dtplyr)
library(dplyr)
library(tidyr)
library(data.table)
library(mc2d)
library(fitdistrplus)
library(goftest)

#settings
iters_var <- 1e+05
iters_unc <- 1e+03

ndvar(iters_var)
ndunc(iters_unc)

#Prob. of fatal IHD, Danish register of causes of death 2015
pIHD_15_19m <- 0     #Men, 15-19y
pIHD_20_24m <- 5e-06 #Men, 20-24y
pIHD_25_29m <- 0     #Men, 25-29y
pIHD_30_34m <- 6e-06 #Etc.
pIHD_35_39m <- 2.2e-05
pIHD_40_44m <- 6.1e-05
pIHD_45_49m <- 0.00016 
pIHD_50_54m <- 0.000288 
pIHD_55_59m <- 0.000493 
pIHD_60_64m <- 0.000884
pIHD_65_69m <- 0.001332
pIHD_70_74m <- 0.002063 
pIHD_75_79m <- 0.003766
pIHD_80_84m <- 0.006765
pIHD_85m <- 0.01641

pIHD_15_19w <- 0 #women, 15-19y
pIHD_20_24w <- 0 #women, 20-24y
pIHD_25_29w <- 0 #women, 25-29y
pIHD_30_34w <- 0 #Etc.
pIHD_35_39w <- 6e-06
pIHD_40_44w <- 2e-05
pIHD_45_49w <- 6.4e-05
pIHD_50_54w <- 0.000101
pIHD_55_59w <- 0.00019 
pIHD_60_64w <- 0.000241
pIHD_65_69w <- 4e-04
pIHD_70_74w <- 0.000808 
pIHD_75_79w <- 0.001622
pIHD_80_84w <- 0.00341
pIHD_85w <- 0.011312

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



# Mozaffarian & Rimm 2006: 36% total reduction in risk of CHD mortality at an intake of 250mg/d EPA+DHA compared to no intake of EPA+DHA.
# Little additional risk reduction was observed at intakes higher than 250mg/d. Based on both primary and secondary prevention
# but relationship confirmed by Harris et al. 2009 in their meta-analysis on only primary prevention studies: Similar risk reduction of
# approximately 35% reduction in cardiac mortality by consumption of 250-500mg/d EPA+DHA was observed in this meta-analysis.
# Linear regression: y = -0,0014x + 1 i.e. a decrease in RR of 0.0014 per mg DHA+EPA per day

# CHD.gradient <- -0.0014 #Decrease in RR per mg DHA+EPA per day - up to 250 mg/day

# From dose-response model of DHA+EPA and risk of fatal CHD
# RR <- 0.0014*x+1

r <- mcstoc(rpert, type = "U", min = -0.0008, mode = -0.0014, max = -0.002)


##### Reference scenario #####

FishDaily <- read.csv("FishDaily.csv")

#Water loss has NOT been taken into account in FoodDaily/FishDaily dataset so amounts are prepared amounts.
#Water loss HAS been taken into account in the substitution scenarios

summary(FishDaily$total.fish) #Total fish consumption for the whole population

#Convert fish consumption amounts into raw weights - assume that all fish is prepared and have a water loss of 20%
wloss.fish <- 1/0.8
FishDaily[,8:23] <- FishDaily[,8:23]*wloss.fish
FishDaily$Totalfish.raw <- rowSums(FishDaily[,8:23])

summary(FishDaily$Totalfish.raw) #For the whole population


DHA_EPA_fish <- read.csv2("DHA_EPA_Fish.csv")
DHA_EPA_fish <- DHA_EPA_fish[,c(1:3)]


DHA_EPA_fish <- tbl_df(DHA_EPA_fish)

DHA_EPA_fish <-  DHA_EPA_fish %>%
  mutate(DHA.EPA = (EPA.g.g + DHAg.g)*1000) %>% #sum columns with concentrations of EPA and DHA + convert to mg/g
  dplyr::select(-c(EPA.g.g, DHAg.g)) #Delete columsn with concentrations of EPA and DHA

DHA_EPA_fish <- as.data.frame(DHA_EPA_fish)
rownames(DHA_EPA_fish) <- DHA_EPA_fish$X #Define row names
DHA_EPA_fish$X <- NULL #Delete first column (now row names)


#Calculate DHA exposure from fish (mg/day)
FishintakeDHA_EPA <- t(FishDaily[,c(8:23)]) #Create data set only with intakes of the different fish species
FishintakeDHA_EPA <- t(FishintakeDHA_EPA * DHA_EPA_fish[1:16,]) ##Multiply intakes of fish species with MeHg conc of the different species to obtain dataset with MeHg exposures from the different species

DHA_EPAcontribution <- colMeans(FishintakeDHA_EPA) #Mean daily exposure in mg/day

barplot(t(DHA_EPA_fish), main = "DHA+EPA concentration", las=2, ylim = c(0,30), ylab = 'mg DHA+EPA/g fish')
barplot(DHA_EPAcontribution, main = "DHA+EPA contribution", las=2, ylim = c(0,150),  ylab = 'mg DHA+EPA/day') #Display mean daily DHA exposure in mg/day - take both consumption and concentration into account

FishDaily$DHA.EPA<- rowSums(FishintakeDHA_EPA) #Add sum of DHA exposures (mg/day) for each individual to the FishDaily dataset

hist(FishDaily$Totalfish.raw, freq = F, breaks = 'fd', main = 'Fish consumption', xlab = 'g/day', ylim = c(0,0.05))

summary(FishDaily$DHA.EPA)

DHA.EPA_ref <- FishDaily[,c(1:4,24:25)]


#Divide individuals according to age groups and sex
agebreaks <- c(15,20,25,30,35,40,45,50,55,60,65,70,75,80)
agelabels= c("15-19","20-24","25-29","30-34","35-39", "40-44", "45-49", "50-54","55-59", "60-64","65-69", "70-74","75-79")

setDT(DHA.EPA_ref)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]



#####    mc2d DALY calc. males  #####

#We assume that current intake of DHA+EPA is reflected in the incidence of fatal CHD in 2015


### Age 15-19, males ###

DALY_15_19refm <- pIHD_15_19m * SEYLL_15_19
summary(DALY_15_19refm)

tDALY_15_19refm <- DALY_15_19refm * pop_15_19m 
summary(tDALY_15_19refm)


### Age 20-24, males ###

DALY_20_24refm <- pIHD_20_24m * SEYLL_20_24
summary(DALY_20_24refm)

tDALY_20_24refm <- DALY_20_24refm * pop_20_24m 
summary(tDALY_20_24refm)


### Age 25-29, males ###

DALY_25_29refm <- pIHD_25_29m * SEYLL_25_29
summary(DALY_25_29refm)

tDALY_25_29refm <- DALY_25_29refm * pop_25_29m 
summary(tDALY_25_29refm)

### Age 30-34, males ###

DALY_30_34refm <- pIHD_30_34m * SEYLL_30_34
summary(DALY_30_34refm)

tDALY_30_34refm <- DALY_30_34refm * pop_30_34m 
summary(tDALY_30_34refm)


### Age 35-39, males ###

DALY_35_39refm <- pIHD_35_39m * SEYLL_35_39
summary(DALY_35_39refm)

tDALY_35_39refm <- DALY_35_39refm * pop_35_39m 
summary(tDALY_35_39refm)


### Age 40-44, males ###

DALY_40_44refm <- pIHD_40_44m * SEYLL_40_44
summary(DALY_40_44refm)

tDALY_40_44refm <- DALY_40_44refm * pop_40_44m 
summary(tDALY_40_44refm)


### Age 45-49, males ###

DALY_45_49refm <- pIHD_45_49m * SEYLL_45_49
summary(DALY_45_49refm)

tDALY_45_49refm <- DALY_45_49refm * pop_45_49m 
summary(tDALY_45_49refm)


### Age 50-54, males ###

DALY_50_54refm <- pIHD_50_54m * SEYLL_50_54
summary(DALY_50_54refm)

tDALY_50_54refm <- DALY_50_54refm * pop_50_54m 
summary(tDALY_50_54refm)


### Age 55-59, males ###

DALY_55_59refm <- pIHD_55_59m * SEYLL_55_59
summary(DALY_55_59refm)

tDALY_55_59refm <- DALY_55_59refm * pop_55_59m 
summary(tDALY_55_59refm)


### Age 60-64, males ###

DALY_60_64refm <- pIHD_60_64m * SEYLL_60_64
summary(DALY_60_64refm)

tDALY_60_64refm <- DALY_60_64refm * pop_60_64m 
summary(tDALY_60_64refm)


### Age 65-69, males ###

DALY_65_69refm <- pIHD_65_69m * SEYLL_65_69
summary(DALY_65_69refm)

tDALY_65_69refm <- DALY_65_69refm * pop_65_69m 
summary(tDALY_65_69refm)


### Age 70-74, males ###

DALY_70_74refm <- pIHD_70_74m * SEYLL_70_74
summary(DALY_70_74refm)

tDALY_70_74refm <- DALY_70_74refm * pop_70_74m 
summary(tDALY_70_74refm)


### Age 75-79, males ###

DALY_75_79refm <- pIHD_75_79m * SEYLL_75_79
summary(DALY_75_79refm)

tDALY_75_79refm <- DALY_75_79refm * pop_75_79m 
summary(tDALY_75_79refm)


### Age 80+, males ###

DALY_80_84refm <- pIHD_80_84m * SEYLL_80_84
summary(DALY_80_84refm)

tDALY_80_84refm <- DALY_80_84refm * pop_80_84m 
summary(tDALY_80_84refm)


### Age 85+, males ###

DALY_85refm <- pIHD_85m * SEYLL_85
summary(DALY_85refm)

tDALY_85refm <- DALY_85refm * pop_85m 
summary(tDALY_85refm)


###############################################    mc2d DALY calc. females  ############################################################

### Age 15-19, females ###

DALY_15_19refw <- pIHD_15_19w * SEYLL_15_19
summary(DALY_15_19refw)

tDALY_15_19refw <- DALY_15_19refw * pop_15_19w 
summary(tDALY_15_19refw)

### Age 20-24, females ###

DALY_20_24refw <- pIHD_20_24w * SEYLL_20_24
summary(DALY_20_24refw)

tDALY_20_24refw <- DALY_20_24refw * pop_20_24w 
summary(tDALY_20_24refw)

### Age 25-29, females ###

DALY_25_29refw <- pIHD_25_29w * SEYLL_25_29
summary(DALY_25_29refw)

tDALY_25_29refw <- DALY_25_29refw * pop_25_29m 
summary(tDALY_25_29refw)

### Age 30-34, females ###

DALY_30_34refw <- pIHD_30_34w * SEYLL_30_34
summary(DALY_30_34refw)

tDALY_30_34refw <- DALY_30_34refw * pop_30_34w 
summary(tDALY_30_34refw)

### Age 35-39, females ###

DALY_35_39refw <- pIHD_35_39w * SEYLL_35_39
summary(DALY_35_39refw)

tDALY_35_39refw <- DALY_35_39refw * pop_35_39w 
summary(tDALY_35_39refw)


### Age 40-44, females ###

DALY_40_44refw <- pIHD_40_44w * SEYLL_40_44
summary(DALY_40_44refw)

tDALY_40_44refw <- DALY_40_44refw * pop_40_44w 
summary(tDALY_40_44refw)


### Age 45-49, females ###

DALY_45_49refw <- pIHD_45_49w * SEYLL_45_49
summary(DALY_45_49refw)

tDALY_45_49refw <- DALY_45_49refw * pop_45_49w 
summary(tDALY_45_49refw)


### Age 50-54, females ###

DALY_50_54refw <- pIHD_50_54w * SEYLL_50_54
summary(DALY_50_54refw)

tDALY_50_54refw <- DALY_50_54refw * pop_50_54w 
summary(tDALY_50_54refw)


### Age 55-59, females ###

DALY_55_59refw <- pIHD_55_59w * SEYLL_55_59
summary(DALY_55_59refw)

tDALY_55_59refw <- DALY_55_59refw * pop_55_59w 
summary(tDALY_55_59refw)


### Age 60-64, females ###

DALY_60_64refw <- pIHD_60_64w * SEYLL_60_64
summary(DALY_60_64refw)

tDALY_60_64refw <- DALY_60_64refw * pop_60_64w 
summary(tDALY_60_64refw)


### Age 65-69, females ###

DALY_65_69refw <- pIHD_65_69w * SEYLL_65_69
summary(DALY_65_69refw)

tDALY_65_69refw <- DALY_65_69refw * pop_65_69w 
summary(tDALY_65_69refw)


### Age 70-74, females ###

DALY_70_74refw <- pIHD_70_74w * SEYLL_70_74
summary(DALY_70_74refw)

tDALY_70_74refw <- DALY_70_74refw * pop_70_74w 
summary(tDALY_70_74refw)



### Age 75-79, females ###

DALY_75_79refw <- pIHD_75_79w * SEYLL_75_79
summary(DALY_75_79refw)

tDALY_75_79refw <- DALY_75_79refw * pop_75_79w 
summary(tDALY_75_79refw)


### Age 80-84, females ###

DALY_80_84refw <- pIHD_80_84w * SEYLL_80_84
summary(DALY_80_84refw)

tDALY_80_84refw <- DALY_80_84refw * pop_80_84w 
summary(tDALY_80_84refw)


### Age 85+, females ###

DALY_85refw <- pIHD_85w * SEYLL_85
summary(DALY_85refw)

tDALY_85refw <- DALY_85refw * pop_85w 
summary(tDALY_85refw)


##### Total DALY #####

DALYtotal.CHD <- tDALY_15_19refm + tDALY_20_24refm + tDALY_25_29refm + tDALY_30_34refm + tDALY_35_39refm + tDALY_40_44refm +
  tDALY_45_49refm + tDALY_45_49refm + tDALY_50_54refm + tDALY_55_59refm + tDALY_60_64refm + tDALY_65_69refm + 
  tDALY_70_74refm + tDALY_75_79refm + tDALY_80_84refm + tDALY_85refm +
  tDALY_15_19refw + tDALY_20_24refw + tDALY_25_29refw + tDALY_30_34refw + tDALY_35_39refw + tDALY_40_44refw + tDALY_45_49refw +
  tDALY_45_49refw + tDALY_50_54refw + tDALY_55_59refw + tDALY_60_64refw + tDALY_65_69refw + tDALY_70_74refw +
  tDALY_75_79refw + tDALY_80_84refw + tDALY_85refw

summary(DALYtotal.CHD)


##### Fitting exposures males #####

#Divide data into agegroups and fit exposure to distribution

#Probability = 0 for males of 15-19 years

#Age 20-24, male
DHA.EPA_20_24m <- subset(DHA.EPA_ref, agegroups=="20-24" & sex=="1", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_20_24m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_20_24m==0)/length(DHA.EPA_20_24m$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_20_24m!=0)/length(DHA.EPA_20_24m$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_20_24m[1,] <- c(pr.no,pr.yes)

DHA.EPA_20_24_posm <- DHA.EPA_20_24m[which(DHA.EPA_20_24m!=0)] #Only intakes that are greater than zero
DHA.EPA_20_24_posm <- as.vector(DHA.EPA_20_24_posm$DHA.EPA)

fit1_20_24m <- fitdist(DHA.EPA_20_24_posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_20_24_posm
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_20_24m$estimate[1], sdlog=fit1_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24m, fitnames="lnorm") 

fit2_20_24m <- fitdist(DHA.EPA_20_24_posm, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_20_24_posm
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_20_24m$estimate[1],fit2_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24m, fitnames="gamma")


cvm.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Anderson-Darling test


#Probability = 0 for males of 25-29 years

#Age 30-34, male
DHA.EPA_30_34m <- subset(DHA.EPA_ref, agegroups=="30-34" & sex=="1", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_30_34m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_30_34m==0)/length(DHA.EPA_30_34m$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_30_34m!=0)/length(DHA.EPA_30_34m$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_30_34m[1,] <- c(pr.no,pr.yes)


DHA.EPA_30_34_posm <- DHA.EPA_30_34m[which(DHA.EPA_30_34m!=0)] #Only intakes that are greater than zero
DHA.EPA_30_34_posm <- as.vector(DHA.EPA_30_34_posm$DHA.EPA)

fit1_30_34m <- fitdist(DHA.EPA_30_34_posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_30_34_posm
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_30_34m$estimate[1], sdlog=fit1_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34m, fitnames="lnorm") 

fit2_30_34m <- fitdist(DHA.EPA_30_34_posm, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_30_34_posm
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_30_34m$estimate[1],fit2_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34m, fitnames="gamma")

cvm.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Anderson-Darling test


#Age 35-39, male
DHA.EPA_35_39m <- subset(DHA.EPA_ref, agegroups=="35-39" & sex=="1", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_35_39m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_35_39m==0)/length(DHA.EPA_35_39m$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_35_39m!=0)/length(DHA.EPA_35_39m$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_35_39m[1,] <- c(pr.no,pr.yes)


DHA.EPA_35_39_posm <- DHA.EPA_35_39m[which(DHA.EPA_35_39m!=0)] #Only intakes that are greater than zero
DHA.EPA_35_39_posm <- as.vector(DHA.EPA_35_39_posm$DHA.EPA)

fit1_35_39m <- fitdist(DHA.EPA_35_39_posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_35_39_posm
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_35_39m$estimate[1], sdlog=fit1_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39m, fitnames="lnorm") 

fit2_35_39m <- fitdist(DHA.EPA_35_39_posm, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_35_39_posm
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_35_39m$estimate[1],fit2_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39m, fitnames="gamma")

cvm.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Anderson-Darling test


#Age 40-44, male
DHA.EPA_40_44m <- subset(DHA.EPA_ref, agegroups=="40-44" & sex=="1", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_40_44m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_40_44m==0)/length(DHA.EPA_40_44m$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_40_44m!=0)/length(DHA.EPA_40_44m$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_40_44m[1,] <- c(pr.no,pr.yes)

DHA.EPA_40_44_posm <- DHA.EPA_40_44m[which(DHA.EPA_40_44m!=0)] #Only intakes that are greater than zero
DHA.EPA_40_44_posm <- as.vector(DHA.EPA_40_44_posm$DHA.EPA)

fit1_40_44m <- fitdist(DHA.EPA_40_44_posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_40_44_posm
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_40_44m$estimate[1], sdlog=fit1_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44m, fitnames="lnorm") 

fit2_40_44m <- fitdist(DHA.EPA_40_44_posm, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_40_44_posm
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_40_44m$estimate[1],fit2_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44m, fitnames="gamma")

cvm.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Anderson-Darling test


#Age 45-49, male
DHA.EPA_45_49m <- subset(DHA.EPA_ref, agegroups=="45-49" & sex=="1", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_45_49m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_45_49m==0)/length(DHA.EPA_45_49m$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_45_49m!=0)/length(DHA.EPA_45_49m$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_45_49m[1,] <- c(pr.no,pr.yes)

DHA.EPA_45_49_posm <- DHA.EPA_45_49m[which(DHA.EPA_45_49m!=0)] #Only intakes that are greater than zero
DHA.EPA_45_49_posm <- as.vector(DHA.EPA_45_49_posm$DHA.EPA)

fit1_45_49m <- fitdist(DHA.EPA_45_49_posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_45_49_posm
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_45_49m$estimate[1], sdlog=fit1_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49m, fitnames="lnorm") 

fit2_45_49m <- fitdist(DHA.EPA_45_49_posm, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_45_49_posm
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_45_49m$estimate[1],fit2_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49m, fitnames="gamma")

cvm.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Anderson-Darling test


#Age 50-54, male
DHA.EPA_50_54m <- subset(DHA.EPA_ref, agegroups=="50-54" & sex=="1", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_50_54m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_50_54m==0)/length(DHA.EPA_50_54m$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_50_54m!=0)/length(DHA.EPA_50_54m$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_50_54m[1,] <- c(pr.no,pr.yes)

DHA.EPA_50_54_posm <- DHA.EPA_50_54m[which(DHA.EPA_50_54m!=0)] #Only intakes that are greater than zero
DHA.EPA_50_54_posm <- as.vector(DHA.EPA_50_54_posm$DHA.EPA)

fit1_50_54m <- fitdist(DHA.EPA_50_54_posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_50_54_posm
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_50_54m$estimate[1], sdlog=fit1_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54m, fitnames="lnorm") 

fit2_50_54m <- fitdist(DHA.EPA_50_54_posm, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_50_54_posm
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_50_54m$estimate[1],fit2_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54m, fitnames="gamma")

cvm.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Anderson-Darling test


#Age 55-59, male
DHA.EPA_55_59m <- subset(DHA.EPA_ref, agegroups=="55-59" & sex=="1", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_55_59m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_55_59m==0)/length(DHA.EPA_55_59m$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_55_59m!=0)/length(DHA.EPA_55_59m$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_55_59m[1,] <- c(pr.no,pr.yes)

DHA.EPA_55_59_posm <- DHA.EPA_55_59m[which(DHA.EPA_55_59m!=0)] #Only intakes that are greater than zero
DHA.EPA_55_59_posm <- as.vector(DHA.EPA_55_59_posm$DHA.EPA)

fit1_55_59m <- fitdist(DHA.EPA_55_59_posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_55_59_posm
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_55_59m$estimate[1], sdlog=fit1_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59m, fitnames="lnorm") 

fit2_55_59m <- fitdist(DHA.EPA_55_59_posm, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_55_59_posm
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_55_59m$estimate[1],fit2_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59m, fitnames="gamma")

cvm.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Anderson-Darling test


#Age 60-64, male
DHA.EPA_60_64m <- subset(DHA.EPA_ref, agegroups=="60-64" & sex=="1", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_60_64m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_60_64m==0)/length(DHA.EPA_60_64m$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_60_64m!=0)/length(DHA.EPA_60_64m$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_60_64m[1,] <- c(pr.no,pr.yes)

DHA.EPA_60_64_posm <- DHA.EPA_60_64m[which(DHA.EPA_60_64m!=0)] #Only intakes that are greater than zero
DHA.EPA_60_64_posm <- as.vector(DHA.EPA_60_64_posm$DHA.EPA)

fit1_60_64m <- fitdist(DHA.EPA_60_64_posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_60_64_posm
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_60_64m$estimate[1], sdlog=fit1_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64m, fitnames="lnorm") 

fit2_60_64m <- fitdist(DHA.EPA_60_64_posm, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_60_64_posm
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_60_64m$estimate[1],fit2_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64m, fitnames="gamma")

cvm.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Anderson-Darling test


#Age 65-69, male
DHA.EPA_65_69m <- subset(DHA.EPA_ref, agegroups=="65-69" & sex=="1", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_65_69m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_65_69m==0)/length(DHA.EPA_65_69m$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_65_69m!=0)/length(DHA.EPA_65_69m$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_65_69m[1,] <- c(pr.no,pr.yes)

DHA.EPA_65_69_posm <- DHA.EPA_65_69m[which(DHA.EPA_65_69m!=0)] #Only intakes that are greater than zero
DHA.EPA_65_69_posm <- as.vector(DHA.EPA_65_69_posm$DHA.EPA)

fit1_65_69m <- fitdist(DHA.EPA_65_69_posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_65_69_posm
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_65_69m$estimate[1], sdlog=fit1_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69m, fitnames="lnorm") 

fit2_65_69m <- fitdist(DHA.EPA_65_69_posm, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_65_69_posm
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_65_69m$estimate[1],fit2_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69m, fitnames="gamma")

cvm.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Anderson-Darling test


#Age 70-74, male
DHA.EPA_70_74m <- subset(DHA.EPA_ref, agegroups=="70-74" & sex=="1", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_70_74m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_70_74m==0)/length(DHA.EPA_70_74m$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_70_74m!=0)/length(DHA.EPA_70_74m$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_70_74m[1,] <- c(pr.no,pr.yes)

DHA.EPA_70_74_posm <- DHA.EPA_70_74m[which(DHA.EPA_70_74m!=0)] #Only intakes that are greater than zero
DHA.EPA_70_74_posm <- as.vector(DHA.EPA_70_74_posm$DHA.EPA)

fit1_70_74m <- fitdist(DHA.EPA_70_74_posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_70_74_posm
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_70_74m$estimate[1], sdlog=fit1_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74m, fitnames="lnorm") 

fit2_70_74m <- fitdist(DHA.EPA_70_74_posm, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_70_74_posm
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_70_74m$estimate[1],fit2_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74m, fitnames="gamma")

cvm.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Anderson-Darling test


#Age 75-84, male
DHA.EPA_75_84m <- subset(DHA.EPA_ref, agegroups=="75-79" & sex=="1", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_75_84m <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_75_84m==0)/length(DHA.EPA_75_84m$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_75_84m!=0)/length(DHA.EPA_75_84m$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_75_84m[1,] <- c(pr.no,pr.yes)

DHA.EPA_75_84_posm <- DHA.EPA_75_84m[which(DHA.EPA_75_84m!=0)] #Only intakes that are greater than zero
DHA.EPA_75_84_posm <- as.vector(DHA.EPA_75_84_posm$DHA.EPA)

fit1_75_84m <- fitdist(DHA.EPA_75_84_posm, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_75_84_posm
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_75_84m$estimate[1], sdlog=fit1_75_84m$estimate[2]), col="red", lty=1)
gofstat(fit1_75_84m, fitnames="lnorm") 

fit2_75_84m <- fitdist(DHA.EPA_75_84_posm, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_75_84_posm
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_75_84m$estimate[1],fit2_75_84m$estimate[2]), col="red", lty=1)
gofstat(fit2_75_84m, fitnames="gamma")

cvm.test(t, plnorm, fit1_75_84m$estimate[1], fit1_75_84m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_84m$estimate[1], fit1_75_84m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_84m$estimate[1], fit2_75_84m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_84m$estimate[1], fit2_75_84m$estimate[2])  #Anderson-Darling test


##### Fitting exposures females #####

#Divide data into agegroups and fit exposure to distribution

#Probability = 0 for males of 15-34 years

#Age 35-39, female
DHA.EPA_35_39w <- subset(DHA.EPA_ref, agegroups=="35-39" & sex=="2", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_35_39w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_35_39w==0)/length(DHA.EPA_35_39w$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_35_39w!=0)/length(DHA.EPA_35_39w$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_35_39w[1,] <- c(pr.no,pr.yes)

DHA.EPA_35_39_posw <- DHA.EPA_35_39w[which(DHA.EPA_35_39w!=0)] #Only intakes that are greater than zero
DHA.EPA_35_39_posw <- as.vector(DHA.EPA_35_39_posw$DHA.EPA)

fit1_35_39w <- fitdist(DHA.EPA_35_39_posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_35_39_posw
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_35_39w$estimate[1], sdlog=fit1_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39w, fitnames="lnorm") 

fit2_35_39w <- fitdist(DHA.EPA_35_39_posw, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_35_39_posw
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_35_39w$estimate[1],fit2_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39w, fitnames="gamma")

cvm.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Anderson-Darling test


#Age 40-44, female
DHA.EPA_40_44w <- subset(DHA.EPA_ref, agegroups=="40-44" & sex=="2", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_40_44w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_40_44w==0)/length(DHA.EPA_40_44w$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_40_44w!=0)/length(DHA.EPA_40_44w$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_40_44w[1,] <- c(pr.no,pr.yes)

DHA.EPA_40_44_posw <- DHA.EPA_40_44w[which(DHA.EPA_40_44w!=0)] #Only intakes that are greater than zero
DHA.EPA_40_44_posw <- as.vector(DHA.EPA_40_44_posw$DHA.EPA)

fit1_40_44w <- fitdist(DHA.EPA_40_44_posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_40_44_posw
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_40_44w$estimate[1], sdlog=fit1_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44w, fitnames="lnorm") 

fit2_40_44w <- fitdist(DHA.EPA_40_44_posw, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_40_44_posw
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_40_44w$estimate[1],fit2_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44w, fitnames="gamma")

cvm.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Anderson-Darling test


#Age 45-49, female
DHA.EPA_45_49w <- subset(DHA.EPA_ref, agegroups=="45-49" & sex=="2", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_45_49w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_45_49w==0)/length(DHA.EPA_45_49w$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_45_49w!=0)/length(DHA.EPA_45_49w$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_45_49w[1,] <- c(pr.no,pr.yes)

DHA.EPA_45_49_posw <- DHA.EPA_45_49w[which(DHA.EPA_45_49w!=0)] #Only intakes that are greater than zero
DHA.EPA_45_49_posw <- as.vector(DHA.EPA_45_49_posw$DHA.EPA)

fit1_45_49w <- fitdist(DHA.EPA_45_49_posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_45_49_posw
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_45_49w$estimate[1], sdlog=fit1_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49w, fitnames="lnorm") 

fit2_45_49w <- fitdist(DHA.EPA_45_49_posw, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_45_49_posw
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_45_49w$estimate[1],fit2_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49w, fitnames="gamma")

cvm.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Anderson-Darling test


#Age 50-54, female
DHA.EPA_50_54w <- subset(DHA.EPA_ref, agegroups=="50-54" & sex=="2", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_50_54w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_50_54w==0)/length(DHA.EPA_50_54w$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_50_54w!=0)/length(DHA.EPA_50_54w$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_50_54w[1,] <- c(pr.no,pr.yes)

DHA.EPA_50_54_posw <- DHA.EPA_50_54w[which(DHA.EPA_50_54w!=0)] #Only intakes that are greater than zero
DHA.EPA_50_54_posw <- as.vector(DHA.EPA_50_54_posw$DHA.EPA)

fit1_50_54w <- fitdist(DHA.EPA_50_54_posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_50_54_posw
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_50_54w$estimate[1], sdlog=fit1_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54w, fitnames="lnorm") 

fit2_50_54w <- fitdist(DHA.EPA_50_54_posw, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_50_54_posw
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_50_54w$estimate[1],fit2_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54w, fitnames="gamma")

cvm.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Anderson-Darling test


#Age 55-59, female
DHA.EPA_55_59w <- subset(DHA.EPA_ref, agegroups=="55-59" & sex=="2", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_55_59w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_55_59w==0)/length(DHA.EPA_55_59w$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_55_59w!=0)/length(DHA.EPA_55_59w$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_55_59w[1,] <- c(pr.no,pr.yes)

DHA.EPA_55_59_posw <- DHA.EPA_55_59w[which(DHA.EPA_55_59w!=0)] #Only intakes that are greater than zero
DHA.EPA_55_59_posw <- as.vector(DHA.EPA_55_59_posw$DHA.EPA)

fit1_55_59w <- fitdist(DHA.EPA_55_59_posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_55_59_posw
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_55_59w$estimate[1], sdlog=fit1_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59w, fitnames="lnorm") 

fit2_55_59w <- fitdist(DHA.EPA_55_59_posw, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_55_59_posw
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_55_59w$estimate[1],fit2_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59w, fitnames="gamma")

cvm.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Anderson-Darling test


#Age 60-64, female
DHA.EPA_60_64w <- subset(DHA.EPA_ref, agegroups=="60-64" & sex=="2", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_60_64w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_60_64w==0)/length(DHA.EPA_60_64w$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_60_64w!=0)/length(DHA.EPA_60_64w$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_60_64w[1,] <- c(pr.no,pr.yes)

DHA.EPA_60_64_posw <- DHA.EPA_60_64w[which(DHA.EPA_60_64w!=0)] #Only intakes that are greater than zero
DHA.EPA_60_64_posw <- as.vector(DHA.EPA_60_64_posw$DHA.EPA)

fit1_60_64w <- fitdist(DHA.EPA_60_64_posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_60_64_posw
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_60_64w$estimate[1], sdlog=fit1_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64w, fitnames="lnorm") 

fit2_60_64w <- fitdist(DHA.EPA_60_64_posw, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_60_64_posw
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_60_64w$estimate[1],fit2_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64w, fitnames="gamma")

cvm.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Anderson-Darling test


#Age 65-69, female
DHA.EPA_65_69w <- subset(DHA.EPA_ref, agegroups=="65-69" & sex=="2", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_65_69w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_65_69w==0)/length(DHA.EPA_65_69w$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_65_69w!=0)/length(DHA.EPA_65_69w$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_65_69w[1,] <- c(pr.no,pr.yes)

DHA.EPA_65_69_posw <- DHA.EPA_65_69w[which(DHA.EPA_65_69w!=0)] #Only intakes that are greater than zero
DHA.EPA_65_69_posw <- as.vector(DHA.EPA_65_69_posw$DHA.EPA)

fit1_65_69w <- fitdist(DHA.EPA_65_69_posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_65_69_posw
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_65_69w$estimate[1], sdlog=fit1_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69w, fitnames="lnorm") 

fit2_65_69w <- fitdist(DHA.EPA_65_69_posw, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_65_69_posw
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_65_69w$estimate[1],fit2_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69w, fitnames="gamma")

cvm.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Anderson-Darling test


#Age 70-74, female
DHA.EPA_70_74w <- subset(DHA.EPA_ref, agegroups=="70-74" & sex=="2", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_70_74w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_70_74w==0)/length(DHA.EPA_70_74w$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_70_74w!=0)/length(DHA.EPA_70_74w$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_70_74w[1,] <- c(pr.no,pr.yes)

DHA.EPA_70_74_posw <- DHA.EPA_70_74w[which(DHA.EPA_70_74w!=0)] #Only intakes that are greater than zero
DHA.EPA_70_74_posw <- as.vector(DHA.EPA_70_74_posw$DHA.EPA)

fit1_70_74w <- fitdist(DHA.EPA_70_74_posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_70_74_posw
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_70_74w$estimate[1], sdlog=fit1_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74w, fitnames="lnorm") 

fit2_70_74w <- fitdist(DHA.EPA_70_74_posw, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_70_74_posw
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_70_74w$estimate[1],fit2_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74w, fitnames="gamma")

cvm.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Anderson-Darling test


#Age 75-84, female
DHA.EPA_75_84w <- subset(DHA.EPA_ref, agegroups=="75-79" & sex=="2", select = DHA.EPA)

#Probability of fish intake
probDHA.EPA_75_84w <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(DHA.EPA_75_84w==0)/length(DHA.EPA_75_84w$DHA.EPA) #probability of zero DHA.EPA exp
pr.yes <- sum(DHA.EPA_75_84w!=0)/length(DHA.EPA_75_84w$DHA.EPA) #probability of DHA.EPA exp
probDHA.EPA_75_84w[1,] <- c(pr.no,pr.yes)

DHA.EPA_75_84_posw <- DHA.EPA_75_84w[which(DHA.EPA_75_84w!=0)] #Only intakes that are greater than zero
DHA.EPA_75_84_posw <- as.vector(DHA.EPA_75_84_posw$DHA.EPA)

fit1_75_84w <- fitdist(DHA.EPA_75_84_posw, 'lnorm') #Fit lognormal distribution to intake amounts
t <- DHA.EPA_75_84_posw
plot(ecdf(t), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, plnorm(x, meanlog=fit1_75_84w$estimate[1], sdlog=fit1_75_84w$estimate[2]), col="red", lty=1)
gofstat(fit1_75_84w, fitnames="lnorm") 

fit2_75_84w <- fitdist(DHA.EPA_75_84_posw, 'gamma', start=list(shape = 1, rate = 1), lower =0.001) #Fit gamma distribution to intake amounts
t2 <- DHA.EPA_75_84_posw
plot(ecdf(t2), lty=1)
x <- seq(0, 100000, length=1000)
lines(x, pgamma(x, fit2_75_84w$estimate[1],fit2_75_84w$estimate[2]), col="red", lty=1)
gofstat(fit2_75_84w, fitnames="gamma")

cvm.test(t, plnorm, fit1_75_84w$estimate[1], fit1_75_84w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_84w$estimate[1], fit1_75_84w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_84w$estimate[1], fit2_75_84w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_84w$estimate[1], fit2_75_84w$estimate[2])  #Anderson-Darling test


###### Peffect ref males #####

# calculate RR(ref) and Peffect(ref) for Peffect(alt) calc in alt scen

#Prob = 0 for men, 15-19y

### Age 20-24, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_20_24m[,1],probDHA.EPA_20_24m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect20_24m <- pIHD_20_24m/RR
summary(pEffect20_24m)


#Prob = 0 for men, 25-29

### Age 30-34, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_30_34m[,1],probDHA.EPA_30_34m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect30_34m <- pIHD_30_34m/RR
summary(pEffect30_34m)


### Age 35-39, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_35_39m[,1],probDHA.EPA_35_39m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect35_39m <- pIHD_35_39m/RR
summary(pEffect35_39m)


### Age 40-44, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_40_44m[,1],probDHA.EPA_40_44m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect40_44m <- pIHD_40_44m/RR
summary(pEffect40_44m)


### Age 45-49, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_45_49m[,1],probDHA.EPA_45_49m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect45_49m <- pIHD_45_49m/RR
summary(pEffect45_49m)



### Age 50-54, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_50_54m[,1],probDHA.EPA_50_54m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect50_54m <- pIHD_50_54m/RR
summary(pEffect50_54m)



### Age 55-59, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_55_59m[,1],probDHA.EPA_55_59m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect55_59m <- pIHD_55_59m/RR
summary(pEffect55_59m)



### Age 60-64, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_60_64m[,1],probDHA.EPA_60_64m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect60_64m <- pIHD_60_64m/RR
summary(pEffect60_64m)



### Age 65-69, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_65_69m[,1],probDHA.EPA_65_69m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)

pEffect65_69m <- pIHD_65_69m/RR
summary(pEffect65_69m)


### Age 70-74, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_70_74m[,1],probDHA.EPA_70_74m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect70_74m <- pIHD_70_74m/RR
summary(pEffect70_74m)



### Age 75-79, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_75_84m[,1],probDHA.EPA_75_84m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_75_84m$estimate[1], fit2_75_84m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect75_79m <- pIHD_75_79m/RR
summary(pEffect75_79m)


### Age 80-84, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_75_84m[,1],probDHA.EPA_75_84m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_75_84m$estimate[1], fit2_75_84m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect80_84m <- pIHD_80_84m/RR
summary(pEffect80_84m)



### Age 85+, males ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_75_84m[,1],probDHA.EPA_75_84m[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_75_84m$estimate[1], fit2_75_84m$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect85m <- pIHD_85m/RR
summary(pEffect85m)


##### Peffect ref females #####

#Prob = 0 for women < 35

### Age 35-39, females ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_35_39w[,1],probDHA.EPA_35_39w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect35_39w <- pIHD_35_39w/RR
summary(pEffect35_39w)




### Age 40-44, females ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_40_44w[,1],probDHA.EPA_40_44w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect40_44w <- pIHD_40_44w/RR
summary(pEffect40_44w)


### Age 45-49, females ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_45_49w[,1],probDHA.EPA_45_49w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect45_49w <- pIHD_45_49w/RR
summary(pEffect45_49w)



### Age 50-54, females ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_50_54w[,1],probDHA.EPA_50_54w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect50_54w <- pIHD_50_54w/RR
summary(pEffect50_54w)



### Age 55-59, females ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_55_59w[,1],probDHA.EPA_55_59w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect55_59w <- pIHD_55_59w/RR
summary(pEffect55_59m)


### Age 60-64, females ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_60_64w[,1],probDHA.EPA_60_64w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect60_64w <- pIHD_60_64w/RR
summary(pEffect60_64w)



### Age 65-69, females ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_65_69w[,1],probDHA.EPA_65_69w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect65_69w <- pIHD_65_69w/RR
summary(pEffect65_69w)



### Age 70-74, females ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_70_74w[,1],probDHA.EPA_70_74w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect70_74w <- pIHD_70_74w/RR
summary(pEffect70_74w)


### Age 75-79, females ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_75_84w[,1],probDHA.EPA_75_84w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_75_84w$estimate[1], fit2_75_84w$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect75_79w <- pIHD_75_79w/RR
summary(pEffect75_79w)



### Age 80-84, females ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_75_84w[,1],probDHA.EPA_75_84w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_75_84w$estimate[1], fit2_75_84w$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect80_84w <- pIHD_80_84w/RR
summary(pEffect80_84w)


### Age 85+, females ###

set.seed(1)
ifexp <- mcstoc(rempiricalD, values = c(0,1), prob = c(probDHA.EPA_75_84w[,1],probDHA.EPA_75_84w[,2]))
exp <- mcstoc(rgamma, type = "V", fit2_75_84w$estimate[1], fit2_75_84w$estimate[2])
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- ifexp * exp * r + 1
summary(RR)


pEffect85w <- pIHD_85w/RR
summary(pEffect85w)


##### Cases ref #####

cases_ref <- pIHD_20_24m * pop_20_24m + pIHD_30_34m * pop_30_34m + pIHD_35_39m * pop_35_39m + pIHD_40_44m * pop_40_44m +
  pIHD_45_49m * pop_45_49m + pIHD_50_54m * pop_50_54m + pIHD_55_59m * pop_55_59m + pIHD_60_64m * pop_60_64m +
  pIHD_65_69m * pop_65_69m + pIHD_70_74m * pop_70_74m + pIHD_75_79m * pop_75_79m + pIHD_80_84m * pop_80_84m +
  pIHD_85m * pop_85m +
  pIHD_35_39w * pop_35_39w + pIHD_40_44w * pop_40_44w + pIHD_45_49w * pop_45_49w + pIHD_50_54w * pop_50_54w +
  pIHD_55_59w * pop_55_59w + pIHD_60_64w * pop_60_64w + pIHD_65_69w * pop_65_69w + pIHD_70_74w * pop_70_74w +
  pIHD_75_79w * pop_75_79w + pIHD_80_84w * pop_80_84w + pIHD_85w * pop_85w


##### Scenario 1 #####

Scenario1 <- read.csv("Scenario1.csv")
# fish intakes are raw weights, ready for exposure assessment.

#Calculate DHA exposure from fish (mg/day)
Scen1_DHA_EPA <- t(Scenario1[,c(55:70)]) #Create data set only with intakes of the different fish species
Scen1_DHA_EPA <- t(Scen1_DHA_EPA * DHA_EPA_fish[1:16,]) ##Multiply intakes of fish species with MeHg conc of the different species to obtain dataset with MeHg exposures from the different species


Scenario1$DHA.EPA <- rowSums(Scen1_DHA_EPA) #Add sum of DHA exposures for each individual to the Scenario1.DHA dataset
Scenario1_CHD <- Scenario1[,c(1:4,54,71)]

#Divide individuals according to age groups and sex

setDT(Scenario1_CHD)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]



###### Exposure males ######

#Divide data into gender- and agegroups

#Age 20-24, male
DHA.EPA_20_24m <- subset(Scenario1_CHD, agegroups=="20-24" & sex=="1", select = DHA.EPA)
DHA.EPA_20_24m <- as.vector(DHA.EPA_20_24m$DHA.EPA)

#Age 30-34, male
DHA.EPA_30_34m <- subset(Scenario1_CHD, agegroups=="30-34" & sex=="1", select = DHA.EPA)
DHA.EPA_30_34m <- as.vector(DHA.EPA_30_34m$DHA.EPA)

#Age 35-39, male
DHA.EPA_35_39m <- subset(Scenario1_CHD, agegroups=="35-39" & sex=="1", select = DHA.EPA)
DHA.EPA_35_39m <- as.vector(DHA.EPA_35_39m$DHA.EPA)

#Age 40-44, male
DHA.EPA_40_44m <- subset(Scenario1_CHD, agegroups=="40-44" & sex=="1", select = DHA.EPA)
DHA.EPA_40_44m <- as.vector(DHA.EPA_40_44m$DHA.EPA)

#Age 45-49, male
DHA.EPA_45_49m <- subset(Scenario1_CHD, agegroups=="45-49" & sex=="1", select = DHA.EPA)
DHA.EPA_45_49m <- as.vector(DHA.EPA_45_49m$DHA.EPA)

#Age 50-54, male
DHA.EPA_50_54m <- subset(Scenario1_CHD, agegroups=="50-54" & sex=="1", select = DHA.EPA)
DHA.EPA_50_54m <- as.vector(DHA.EPA_50_54m$DHA.EPA)

#Age 55-59, male
DHA.EPA_55_59m <- subset(Scenario1_CHD, agegroups=="55-59" & sex=="1", select = DHA.EPA)
DHA.EPA_55_59m <- as.vector(DHA.EPA_55_59m$DHA.EPA)

#Age 60-64, male
DHA.EPA_60_64m <- subset(Scenario1_CHD, agegroups=="60-64" & sex=="1", select = DHA.EPA)
DHA.EPA_60_64m <- as.vector(DHA.EPA_60_64m$DHA.EPA)

#Age 65-69, male
DHA.EPA_65_69m <- subset(Scenario1_CHD, agegroups=="65-69" & sex=="1", select = DHA.EPA)
DHA.EPA_65_69m <- as.vector(DHA.EPA_65_69m$DHA.EPA)

#Age 70-74, male
DHA.EPA_70_74m <- subset(Scenario1_CHD, agegroups=="70-74" & sex=="1", select = DHA.EPA)
DHA.EPA_70_74m <- as.vector(DHA.EPA_70_74m$DHA.EPA)

#Age 75-79, male
DHA.EPA_75_79m <- subset(Scenario1_CHD, agegroups=="75-79" & sex=="1", select = DHA.EPA)
DHA.EPA_75_79m <- as.vector(DHA.EPA_75_79m$DHA.EPA)


##### Exposure females #####

#Age 35-39, female
DHA.EPA_35_39w <- subset(Scenario1_CHD, agegroups=="35-39" & sex=="2", select = DHA.EPA)
DHA.EPA_35_39w <- as.vector(DHA.EPA_35_39w$DHA.EPA)

#Age 40-44, female
DHA.EPA_40_44w <- subset(Scenario1_CHD, agegroups=="40-44" & sex=="2", select = DHA.EPA)
DHA.EPA_40_44w <- as.vector(DHA.EPA_40_44w$DHA.EPA)

#Age 45-49, female
DHA.EPA_45_49w <- subset(Scenario1_CHD, agegroups=="45-49" & sex=="2", select = DHA.EPA)
DHA.EPA_45_49w <- as.vector(DHA.EPA_45_49w$DHA.EPA)

#Age 50-54, female
DHA.EPA_50_54w <- subset(Scenario1_CHD, agegroups=="50-54" & sex=="2", select = DHA.EPA)
DHA.EPA_50_54w <- as.vector(DHA.EPA_50_54w$DHA.EPA)

#Age 55-59, female
DHA.EPA_55_59w <- subset(Scenario1_CHD, agegroups=="55-59" & sex=="2", select = DHA.EPA)
DHA.EPA_55_59w <- as.vector(DHA.EPA_55_59w$DHA.EPA)

#Age 60-64, female
DHA.EPA_60_64w <- subset(Scenario1_CHD, agegroups=="60-64" & sex=="2", select = DHA.EPA)
DHA.EPA_60_64w <- as.vector(DHA.EPA_60_64w$DHA.EPA)

#Age 65-69, female
DHA.EPA_65_69w <- subset(Scenario1_CHD, agegroups=="65-69" & sex=="2", select = DHA.EPA)
DHA.EPA_65_69w <- as.vector(DHA.EPA_65_69w$DHA.EPA)

#Age 70-74, female
DHA.EPA_70_74w <- subset(Scenario1_CHD, agegroups=="70-74" & sex=="2", select = DHA.EPA)
DHA.EPA_70_74w <- as.vector(DHA.EPA_70_74w$DHA.EPA)

#Age 75-79, female
DHA.EPA_75_79w <- subset(Scenario1_CHD, agegroups=="75-79" & sex=="2", select = DHA.EPA)
DHA.EPA_75_79w <- as.vector(DHA.EPA_75_79w$DHA.EPA)


##### DALY scen 1 male #####

#Prob = 0 for men, 15-19y

### Age 20-24, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_20_24m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect20_24m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_20_24m <- pEffect * pop_20_24m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_20_24scen1m <- pEffect*SEYLL_20_24
summary(DALY_20_24scen1m)


tDALY_20_24scen1m <- DALY_20_24scen1m * pop_20_24m 
summary(tDALY_20_24scen1m)


#Prob = 0 for men, 25-29y

### Age 30-34, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_30_34m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect30_34m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_30_34m <- pEffect * pop_30_34m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_30_34scen1m <- pEffect*SEYLL_30_34
summary(DALY_30_34scen1m)


tDALY_30_34scen1m <- DALY_30_34scen1m * pop_30_34m 
summary(tDALY_30_34scen1m)


### Age 35-39, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_35_39m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect35_39m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_35_39m <- pEffect * pop_35_39m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_35_39scen1m <- pEffect*SEYLL_35_39
summary(DALY_35_39scen1m)


tDALY_35_39scen1m <- DALY_35_39scen1m * pop_35_39m 
summary(tDALY_35_39scen1m)


### Age 40-44, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_40_44m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect40_44m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_40_44m <- pEffect * pop_40_44m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_40_44scen1m <- pEffect*SEYLL_40_44
summary(DALY_40_44scen1m)


tDALY_40_44scen1m <- DALY_40_44scen1m * pop_40_44m 
summary(tDALY_40_44scen1m)


### Age 45-49, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_45_49m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect45_49m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_45_49m <- pEffect * pop_45_49m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_45_49scen1m <- pEffect*SEYLL_45_49
summary(DALY_45_49scen1m)


tDALY_45_49scen1m <- DALY_45_49scen1m * pop_45_49m 
summary(tDALY_45_49scen1m)


### Age 50-54, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_50_54m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect50_54m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_50_54m <- pEffect * pop_50_54m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_50_54scen1m <- pEffect*SEYLL_50_54
summary(DALY_50_54scen1m)


tDALY_50_54scen1m <- DALY_50_54scen1m * pop_50_54m 
summary(tDALY_50_54scen1m)


### Age 55-59, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_55_59m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect55_59m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_55_59m <- pEffect * pop_55_59m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_55_59scen1m <- pEffect*SEYLL_55_59
summary(DALY_55_59scen1m)


tDALY_55_59scen1m <- DALY_55_59scen1m * pop_55_59m 
summary(tDALY_55_59scen1m)


### Age 60-64, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_60_64m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect60_64m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_60_64m <- pEffect * pop_60_64m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_60_64scen1m <- pEffect*SEYLL_60_64
summary(DALY_60_64scen1m)


tDALY_60_64scen1m <- DALY_60_64scen1m * pop_60_64m 
summary(tDALY_60_64scen1m)


### Age 65-69, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_65_69m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect65_69m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_65_69m <- pEffect * pop_65_69m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_65_69scen1m <- pEffect*SEYLL_65_69
summary(DALY_65_69scen1m)


tDALY_65_69scen1m <- DALY_65_69scen1m * pop_65_69m 
summary(tDALY_65_69scen1m)



### Age 70-74, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_70_74m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect70_74m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_70_74m <- pEffect * pop_70_74m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_70_74scen1m <- pEffect*SEYLL_70_74
summary(DALY_70_74scen1m)


tDALY_70_74scen1m <- DALY_70_74scen1m * pop_70_74m 
summary(tDALY_70_74scen1m)



### Age 75-79, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect75_79m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_75_79m <- pEffect * pop_75_79m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_75_79scen1m <- pEffect*SEYLL_75_79
summary(DALY_75_79scen1m)


tDALY_75_79scen1m <- DALY_75_79scen1m * pop_75_79m 
summary(tDALY_75_79scen1m)



### Age 80-84, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79m) #Extrapolate consumption data to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect80_84m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_80_84m <- pEffect * pop_80_84m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_80_84scen1m <- pEffect*SEYLL_80_84
summary(DALY_80_84scen1m)


tDALY_80_84scen1m <- DALY_80_84scen1m * pop_80_84m 
summary(tDALY_80_84scen1m)


### Age 85, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79m) #Extrapolate consumption data to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)

pEffect <- pEffect85m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_85m <- pEffect * pop_85m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_85scen1m <- pEffect*SEYLL_85
summary(DALY_85scen1m)


tDALY_85scen1m <- DALY_85scen1m * pop_85m 
summary(tDALY_85scen1m)


##### DALY scen 1 females #####


### Age 35-39, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_35_39w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect35_39w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_35_39w <- pEffect * pop_35_39w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_35_39scen1w <- pEffect*SEYLL_35_39
summary(DALY_35_39scen1w)


tDALY_35_39scen1w <- DALY_35_39scen1w * pop_35_39w 
summary(tDALY_35_39scen1w)


### Age 40-44, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_40_44w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect40_44w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_40_44w <- pEffect * pop_40_44w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_40_44scen1w <- pEffect*SEYLL_40_44
summary(DALY_40_44scen1w)


tDALY_40_44scen1w <- DALY_40_44scen1w * pop_40_44w 
summary(tDALY_40_44scen1w)


### Age 45-49, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_45_49w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect45_49w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_45_49w <- pEffect * pop_45_49w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_45_49scen1w <- pEffect*SEYLL_45_49
summary(DALY_45_49scen1w)


tDALY_45_49scen1w <- DALY_45_49scen1w * pop_45_49w 
summary(tDALY_45_49scen1w)


### Age 50-54, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_50_54w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect50_54w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_50_54w <- pEffect * pop_50_54w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_50_54scen1w <- pEffect*SEYLL_50_54
summary(DALY_50_54scen1w)


tDALY_50_54scen1w <- DALY_50_54scen1w * pop_50_54w 
summary(tDALY_50_54scen1w)


### Age 55-59, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_55_59w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect55_59w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_55_59w <- pEffect * pop_55_59w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_55_59scen1w <- pEffect*SEYLL_55_59
summary(DALY_55_59scen1w)


tDALY_55_59scen1w <- DALY_55_59scen1w * pop_55_59w 
summary(tDALY_55_59scen1w)



### Age 60-64, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_60_64w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect60_64w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_60_64w <- pEffect * pop_60_64w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_60_64scen1w <- pEffect*SEYLL_60_64
summary(DALY_60_64scen1w)


tDALY_60_64scen1w <- DALY_60_64scen1w * pop_60_64w 
summary(tDALY_60_64scen1w)


### Age 65-69, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_65_69w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect65_69w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_65_69w <- pEffect * pop_65_69w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_65_69scen1w <- pEffect*SEYLL_65_69
summary(DALY_65_69scen1w)


tDALY_65_69scen1w <- DALY_65_69scen1w * pop_65_69w 
summary(tDALY_65_69scen1w)



### Age 70-74, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_70_74w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect70_74w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_70_74w <- pEffect * pop_70_74w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_70_74scen1w <- pEffect*SEYLL_70_74
summary(DALY_70_74scen1w)


tDALY_70_74scen1w <- DALY_70_74scen1w * pop_70_74w 
summary(tDALY_70_74scen1w)


### Age 75-79, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect75_79w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_75_79w <- pEffect * pop_75_79w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_75_79scen1w <- pEffect*SEYLL_75_79
summary(DALY_75_79scen1w)


tDALY_75_79scen1w <- DALY_75_79scen1w * pop_75_79w 
summary(tDALY_75_79scen1w)


### Age 80-84, females ### 

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79w) #Extrapolate consumption data to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect80_84w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_80_84w <- pEffect * pop_80_84w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_80_84scen1w <- pEffect*SEYLL_80_84
summary(DALY_80_84scen1w)


tDALY_80_84scen1w <- DALY_80_84scen1w * pop_80_84w 
summary(tDALY_80_84scen1w)


### Age 85, females ### 

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79w) #Extrapolate consumption data to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect85w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_85w <- pEffect * pop_85w



#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_85scen1w <- pEffect*SEYLL_85
summary(DALY_85scen1w)


tDALY_85scen1w <- DALY_85scen1w * pop_85w 
summary(tDALY_85scen1w)


##### Total DALY #####

DALYtotal.CHD.scen1 <- tDALY_20_24scen1m + tDALY_30_34scen1m + tDALY_35_39scen1m + tDALY_40_44scen1m + tDALY_45_49scen1m + tDALY_45_49scen1m + tDALY_50_54scen1m + tDALY_55_59scen1m + 
  tDALY_60_64scen1m + tDALY_65_69scen1m + tDALY_70_74scen1m + tDALY_75_79scen1m + tDALY_80_84scen1m + tDALY_85scen1m +
  tDALY_35_39scen1w + tDALY_40_44scen1w + tDALY_45_49scen1w + tDALY_45_49scen1w + tDALY_50_54scen1w + tDALY_55_59scen1w + 
  tDALY_60_64scen1w + tDALY_65_69scen1w + tDALY_70_74scen1w + tDALY_75_79scen1w + tDALY_80_84scen1w + tDALY_85scen1w

summary(DALYtotal.CHD.scen1)


##### Cases scenario 1 ######

cases_scen1 <- cases_20_24m + cases_30_34m + cases_35_39m + cases_40_44m + cases_45_49m + cases_50_54m + cases_55_59m +
  cases_60_64m + cases_65_69m + cases_70_74m + cases_75_79m + cases_80_84m + cases_85m + 
  cases_35_39w + cases_40_44w + cases_45_49w + cases_50_54w + cases_55_59w + cases_60_64w + cases_65_69w +
  cases_70_74w + cases_75_79w + cases_80_84w + cases_85w

cases_diff_scen1 <- cases_scen1 - cases_ref
summary(cases_diff_scen1)


###### Scenario 2 ######

Scenario2 <- read.csv("Scenario2.csv")
# fish intakes are raw weights, ready for exposure assessment.

#Calculate DHA exposure from fish (mg/day)
Scen2_DHA_EPA <- t(Scenario2[,c(55:60)]) #Create data set only with intakes of the different fish species
Scen2_DHA_EPA <- t(Scen2_DHA_EPA * DHA_EPA_fish[c(3,9,10,11,12,16),]) ##Multiply intakes of fish species with MeHg conc of the different species to obtain dataset with MeHg exposures from the different species


Scenario2$DHA.EPA <- rowSums(Scen2_DHA_EPA) #Add sum of DHA exposures for each individual to the Scenario2.DHA dataset
Scenario2_CHD <- Scenario2[,c(1:4,54,61)]

#Divide individuals according to age groups and sex

setDT(Scenario2_CHD)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Exposure males #####

#Divide data into gender- and agegroups

#Prob = 0  for males 15-19y and 25-29y

#Age 20-24, male
DHA.EPA_20_24m <- subset(Scenario2_CHD, agegroups=="20-24" & sex=="1", select = DHA.EPA)
DHA.EPA_20_24m <- as.vector(DHA.EPA_20_24m$DHA.EPA)

#Age 30-34, male
DHA.EPA_30_34m <- subset(Scenario2_CHD, agegroups=="30-34" & sex=="1", select = DHA.EPA)
DHA.EPA_30_34m <- as.vector(DHA.EPA_30_34m$DHA.EPA)

#Age 35-39, male
DHA.EPA_35_39m <- subset(Scenario2_CHD, agegroups=="35-39" & sex=="1", select = DHA.EPA)
DHA.EPA_35_39m <- as.vector(DHA.EPA_35_39m$DHA.EPA)

#Age 40-44, male
DHA.EPA_40_44m <- subset(Scenario2_CHD, agegroups=="40-44" & sex=="1", select = DHA.EPA)
DHA.EPA_40_44m <- as.vector(DHA.EPA_40_44m$DHA.EPA)

#Age 45-49, male
DHA.EPA_45_49m <- subset(Scenario2_CHD, agegroups=="45-49" & sex=="1", select = DHA.EPA)
DHA.EPA_45_49m <- as.vector(DHA.EPA_45_49m$DHA.EPA)

#Age 50-54, male
DHA.EPA_50_54m <- subset(Scenario2_CHD, agegroups=="50-54" & sex=="1", select = DHA.EPA)
DHA.EPA_50_54m <- as.vector(DHA.EPA_50_54m$DHA.EPA)

#Age 55-59, male
DHA.EPA_55_59m <- subset(Scenario2_CHD, agegroups=="55-59" & sex=="1", select = DHA.EPA)
DHA.EPA_55_59m <- as.vector(DHA.EPA_55_59m$DHA.EPA)

#Age 60-64, male
DHA.EPA_60_64m <- subset(Scenario2_CHD, agegroups=="60-64" & sex=="1", select = DHA.EPA)
DHA.EPA_60_64m <- as.vector(DHA.EPA_60_64m$DHA.EPA)

#Age 65-69, male
DHA.EPA_65_69m <- subset(Scenario2_CHD, agegroups=="65-69" & sex=="1", select = DHA.EPA)
DHA.EPA_65_69m <- as.vector(DHA.EPA_65_69m$DHA.EPA)

#Age 70-74, male
DHA.EPA_70_74m <- subset(Scenario2_CHD, agegroups=="70-74" & sex=="1", select = DHA.EPA)
DHA.EPA_70_74m <- as.vector(DHA.EPA_70_74m$DHA.EPA)

#Age 75-79, male
DHA.EPA_75_79m <- subset(Scenario2_CHD, agegroups=="75-79" & sex=="1", select = DHA.EPA)
DHA.EPA_75_79m <- as.vector(DHA.EPA_75_79m$DHA.EPA)


##### Exposure females #####

#Prob = 0 for females > 35y

#Age 35-39, female
DHA.EPA_35_39w <- subset(Scenario2_CHD, agegroups=="35-39" & sex=="2", select = DHA.EPA)
DHA.EPA_35_39w <- as.vector(DHA.EPA_35_39w$DHA.EPA)

#Age 40-44, female
DHA.EPA_40_44w <- subset(Scenario2_CHD, agegroups=="40-44" & sex=="2", select = DHA.EPA)
DHA.EPA_40_44w <- as.vector(DHA.EPA_40_44w$DHA.EPA)

#Age 45-49, female
DHA.EPA_45_49w <- subset(Scenario2_CHD, agegroups=="45-49" & sex=="2", select = DHA.EPA)
DHA.EPA_45_49w <- as.vector(DHA.EPA_45_49w$DHA.EPA)

#Age 50-54, female
DHA.EPA_50_54w <- subset(Scenario2_CHD, agegroups=="50-54" & sex=="2", select = DHA.EPA)
DHA.EPA_50_54w <- as.vector(DHA.EPA_50_54w$DHA.EPA)

#Age 55-59, female
DHA.EPA_55_59w <- subset(Scenario2_CHD, agegroups=="55-59" & sex=="2", select = DHA.EPA)
DHA.EPA_55_59w <- as.vector(DHA.EPA_55_59w$DHA.EPA)

#Age 60-64, female
DHA.EPA_60_64w <- subset(Scenario2_CHD, agegroups=="60-64" & sex=="2", select = DHA.EPA)
DHA.EPA_60_64w <- as.vector(DHA.EPA_60_64w$DHA.EPA)

#Age 65-69, female
DHA.EPA_65_69w <- subset(Scenario2_CHD, agegroups=="65-69" & sex=="2", select = DHA.EPA)
DHA.EPA_65_69w <- as.vector(DHA.EPA_65_69w$DHA.EPA)

#Age 70-74, female
DHA.EPA_70_74w <- subset(Scenario2_CHD, agegroups=="70-74" & sex=="2", select = DHA.EPA)
DHA.EPA_70_74w <- as.vector(DHA.EPA_70_74w$DHA.EPA)

#Age 75-79, female
DHA.EPA_75_79w <- subset(Scenario2_CHD, agegroups=="75-79" & sex=="2", select = DHA.EPA)
DHA.EPA_75_79w <- as.vector(DHA.EPA_75_79w$DHA.EPA)


##### DALY scen 2 male #####

#Prob = 0 for men, 15-19y

### Age 20-24, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_20_24m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect20_24m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_20_24m <- pEffect * pop_20_24m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_20_24scen2m <- pEffect*SEYLL_20_24
summary(DALY_20_24scen2m)


tDALY_20_24scen2m <- DALY_20_24scen2m * pop_20_24m 
summary(tDALY_20_24scen2m)


#Prob = 0 for men, 25-29y

### Age 30-34, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_30_34m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect30_34m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_30_34m <- pEffect * pop_30_34m



#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_30_34scen2m <- pEffect*SEYLL_30_34
summary(DALY_30_34scen2m)


tDALY_30_34scen2m <- DALY_30_34scen2m * pop_30_34m 
summary(tDALY_30_34scen2m)


### Age 35-39, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_35_39m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect35_39m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_35_39m <- pEffect * pop_35_39m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_35_39scen2m <- pEffect*SEYLL_35_39
summary(DALY_35_39scen2m)


tDALY_35_39scen2m <- DALY_35_39scen2m * pop_35_39m 
summary(tDALY_35_39scen2m)


### Age 40-44, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_40_44m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect40_44m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_40_44m <- pEffect * pop_40_44m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_40_44scen2m <- pEffect*SEYLL_40_44
summary(DALY_40_44scen2m)


tDALY_40_44scen2m <- DALY_40_44scen2m * pop_40_44m 
summary(tDALY_40_44scen2m)


### Age 45-49, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_45_49m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect45_49m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_45_49m <- pEffect * pop_45_49m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_45_49scen2m <- pEffect*SEYLL_45_49
summary(DALY_45_49scen2m)


tDALY_45_49scen2m <- DALY_45_49scen2m * pop_45_49m 
summary(tDALY_45_49scen2m)


### Age 50-54, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_50_54m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect50_54m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_50_54m <- pEffect * pop_50_54m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_50_54scen2m <- pEffect*SEYLL_50_54
summary(DALY_50_54scen2m)


tDALY_50_54scen2m <- DALY_50_54scen2m * pop_50_54m 
summary(tDALY_50_54scen2m)


### Age 55-59, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_55_59m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect55_59m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_55_59m <- pEffect * pop_55_59m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_55_59scen2m <- pEffect*SEYLL_55_59
summary(DALY_55_59scen2m)


tDALY_55_59scen2m <- DALY_55_59scen2m * pop_55_59m 
summary(tDALY_55_59scen2m)


### Age 60-64, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_60_64m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect60_64m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_60_64m <- pEffect * pop_60_64m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_60_64scen2m <- pEffect*SEYLL_60_64
summary(DALY_60_64scen2m)


tDALY_60_64scen2m <- DALY_60_64scen2m * pop_60_64m 
summary(tDALY_60_64scen2m)


### Age 65-69, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_65_69m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect65_69m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_65_69m <- pEffect * pop_65_69m



#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_65_69scen2m <- pEffect*SEYLL_65_69
summary(DALY_65_69scen2m)


tDALY_65_69scen2m <- DALY_65_69scen2m * pop_65_69m 
summary(tDALY_65_69scen2m)



### Age 70-74, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_70_74m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect70_74m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_70_74m <- pEffect * pop_70_74m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_70_74scen2m <- pEffect*SEYLL_70_74
summary(DALY_70_74scen2m)


tDALY_70_74scen2m <- DALY_70_74scen2m * pop_70_74m 
summary(tDALY_70_74scen2m)


### Age 75-79, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)

pEffect <- pEffect75_79m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_75_79m <- pEffect * pop_75_79m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_75_79scen2m <- pEffect*SEYLL_75_79
summary(DALY_75_79scen2m)


tDALY_75_79scen2m <- DALY_75_79scen2m * pop_75_79m 
summary(tDALY_75_79scen2m)


### Age 80-84, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79m) #extrapolate to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect80_84m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_80_84m <- pEffect * pop_80_84m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_80_84scen2m <- pEffect*SEYLL_80_84
summary(DALY_80_84scen2m)


tDALY_80_84scen2m <- DALY_80_84scen2m * pop_80_84m 
summary(tDALY_80_84scen2m)


### Age 85, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79m) #extrapolate to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect85m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_85m <- pEffect * pop_85m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_85scen2m <- pEffect*SEYLL_85
summary(DALY_85scen2m)


tDALY_85scen2m <- DALY_85scen2m * pop_85m 
summary(tDALY_85scen2m)



################################################ DALY scen 2 females ########################################################

#Prob = 0 for women <35y

### Age 35-39, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_35_39w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect35_39w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_35_39w <- pEffect * pop_35_39w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_35_39scen2w <- pEffect*SEYLL_35_39
summary(DALY_35_39scen2w)


tDALY_35_39scen2w <- DALY_35_39scen2w * pop_35_39w 
summary(tDALY_35_39scen2w)


### Age 40-44, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_40_44w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect40_44w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_40_44w <- pEffect * pop_40_44w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_40_44scen2w <- pEffect*SEYLL_40_44
summary(DALY_40_44scen2w)


tDALY_40_44scen2w <- DALY_40_44scen2w * pop_40_44w 
summary(tDALY_40_44scen2w)


### Age 45-49, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_45_49w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect45_49w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_45_49w <- pEffect * pop_45_49w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_45_49scen2w <- pEffect*SEYLL_45_49
summary(DALY_45_49scen2w)


tDALY_45_49scen2w <- DALY_45_49scen2w * pop_45_49w 
summary(tDALY_45_49scen2w)


### Age 50-54, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_50_54w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect50_54w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_50_54w <- pEffect * pop_50_54w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_50_54scen2w <- pEffect*SEYLL_50_54
summary(DALY_50_54scen2w)


tDALY_50_54scen2w <- DALY_50_54scen2w * pop_50_54w 
summary(tDALY_50_54scen2w)


### Age 55-59, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_55_59w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect55_59w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_55_59w <- pEffect * pop_55_59w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_55_59scen2w <- pEffect*SEYLL_55_59
summary(DALY_55_59scen2w)


tDALY_55_59scen2w <- DALY_55_59scen2w * pop_55_59w 
summary(tDALY_55_59scen2w)


### Age 60-64, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_60_64w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect60_64w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_60_64w <- pEffect * pop_60_64w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_60_64scen2w <- pEffect*SEYLL_60_64
summary(DALY_60_64scen2w)


tDALY_60_64scen2w <- DALY_60_64scen2w * pop_60_64w 
summary(tDALY_60_64scen2w)


### Age 65-69, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_65_69w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect65_69w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_65_69w <- pEffect * pop_65_69w

#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_65_69scen2w <- pEffect*SEYLL_65_69
summary(DALY_65_69scen2w)


tDALY_65_69scen2w <- DALY_65_69scen2w * pop_65_69w 
summary(tDALY_65_69scen2w)



### Age 70-74, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_70_74w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)

pEffect <- pEffect70_74w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_70_74w <- pEffect * pop_70_74w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_70_74scen2w <- pEffect*SEYLL_70_74
summary(DALY_70_74scen2w)


tDALY_70_74scen2w <- DALY_70_74scen2w * pop_70_74w 
summary(tDALY_70_74scen2w)


### Age 75-79, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect75_79w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_75_79w <- pEffect * pop_75_79w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_75_79scen2w <- pEffect*SEYLL_75_79
summary(DALY_75_79scen2w)


tDALY_75_79scen2w <- DALY_75_79scen2w * pop_75_79w 
summary(tDALY_75_79scen2w)



### Age 80-84, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79w) #Extrapolate intakes to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect80_84w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_80_84w <- pEffect * pop_80_84w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_80_84scen2w <- pEffect*SEYLL_80_84
summary(DALY_80_84scen2w)


tDALY_80_84scen2w <- DALY_80_84scen2w * pop_80_84w 
summary(tDALY_80_84scen2w)


### Age 85+, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79w) #Extrapolate intakes to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect85w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_85w <- pEffect * pop_85w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_85scen2w <- pEffect*SEYLL_85
summary(DALY_85scen2w)


tDALY_85scen2w <- DALY_85scen2w * pop_85w 
summary(tDALY_85scen2w)


##### Total DALY #####

DALYtotal.CHD.scen2 <- tDALY_20_24scen2m + tDALY_30_34scen2m + tDALY_35_39scen2m + tDALY_40_44scen2m + tDALY_45_49scen2m +
  tDALY_45_49scen2m + tDALY_50_54scen2m + tDALY_55_59scen2m + tDALY_60_64scen2m + tDALY_65_69scen2m + tDALY_70_74scen2m +
  tDALY_75_79scen2m + tDALY_80_84scen2m + tDALY_85scen2m + tDALY_35_39scen2w + tDALY_40_44scen2w + tDALY_45_49scen2w +
  tDALY_45_49scen2w + tDALY_50_54scen2w + tDALY_55_59scen2w + tDALY_60_64scen2w + tDALY_65_69scen2w + tDALY_70_74scen2w + 
  tDALY_75_79scen2w + tDALY_80_84scen2w + tDALY_85scen2w

summary(DALYtotal.CHD.scen2)



##### Cases scenario 2 #####

cases_scen2 <- cases_20_24m + cases_30_34m + cases_35_39m + cases_40_44m + cases_45_49m + cases_50_54m + cases_55_59m +
  cases_60_64m + cases_65_69m + cases_70_74m + cases_75_79m + cases_80_84m + cases_85m + 
  cases_35_39w + cases_40_44w + cases_45_49w + cases_50_54w + cases_55_59w + cases_60_64w + cases_65_69w +
  cases_70_74w + cases_75_79w + cases_80_84w + cases_85w

cases_diff_scen2 <- cases_scen2 - cases_ref
summary(cases_diff_scen2)



##### Scenario 3 #####

Scenario3 <- read.csv("Scenario3.csv")
# fish intakes are raw weights, ready for exposure assessment.

#Calculate DHA exposure from fish (mg/day)
Scen3_DHA_EPA <- t(Scenario3[,c(55:64)]) #Create data set only with intakes of the different fish species
Scen3_DHA_EPA <- t(Scen3_DHA_EPA * DHA_EPA_fish[c(1,2,4,5,6,7,8,13,14,15),]) ##Multiply intakes of fish species with MeHg conc of the different species to obtain dataset with MeHg exposures from the different species


Scenario3$DHA.EPA <- rowSums(Scen3_DHA_EPA) #Add sum of DHA exposures for each individual to the Scenario3.DHA dataset
Scenario3_CHD <- Scenario3[,c(1:4,54,65)]

#Divide individuals according to age groups and sex

setDT(Scenario3_CHD)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Exposure males ######

#Divide data into gender- and agegroups

#Prob=0 for men, 15-19y

#Age 20-24, male
DHA.EPA_20_24m <- subset(Scenario3_CHD, agegroups=="20-24" & sex=="1", select = DHA.EPA)
DHA.EPA_20_24m <- as.vector(DHA.EPA_20_24m$DHA.EPA)

#Prob=0 for men, 25-29y

#Age 30-34, male
DHA.EPA_30_34m <- subset(Scenario3_CHD, agegroups=="30-34" & sex=="1", select = DHA.EPA)
DHA.EPA_30_34m <- as.vector(DHA.EPA_30_34m$DHA.EPA)

#Age 35-39, male
DHA.EPA_35_39m <- subset(Scenario3_CHD, agegroups=="35-39" & sex=="1", select = DHA.EPA)
DHA.EPA_35_39m <- as.vector(DHA.EPA_35_39m$DHA.EPA)

#Age 40-44, male
DHA.EPA_40_44m <- subset(Scenario3_CHD, agegroups=="40-44" & sex=="1", select = DHA.EPA)
DHA.EPA_40_44m <- as.vector(DHA.EPA_40_44m$DHA.EPA)

#Age 45-49, male
DHA.EPA_45_49m <- subset(Scenario3_CHD, agegroups=="45-49" & sex=="1", select = DHA.EPA)
DHA.EPA_45_49m <- as.vector(DHA.EPA_45_49m$DHA.EPA)

#Age 50-54, male
DHA.EPA_50_54m <- subset(Scenario3_CHD, agegroups=="50-54" & sex=="1", select = DHA.EPA)
DHA.EPA_50_54m <- as.vector(DHA.EPA_50_54m$DHA.EPA)

#Age 55-59, male
DHA.EPA_55_59m <- subset(Scenario3_CHD, agegroups=="55-59" & sex=="1", select = DHA.EPA)
DHA.EPA_55_59m <- as.vector(DHA.EPA_55_59m$DHA.EPA)

#Age 60-64, male
DHA.EPA_60_64m <- subset(Scenario3_CHD, agegroups=="60-64" & sex=="1", select = DHA.EPA)
DHA.EPA_60_64m <- as.vector(DHA.EPA_60_64m$DHA.EPA)

#Age 65-69, male
DHA.EPA_65_69m <- subset(Scenario3_CHD, agegroups=="65-69" & sex=="1", select = DHA.EPA)
DHA.EPA_65_69m <- as.vector(DHA.EPA_65_69m$DHA.EPA)

#Age 70-74, male
DHA.EPA_70_74m <- subset(Scenario3_CHD, agegroups=="70-74" & sex=="1", select = DHA.EPA)
DHA.EPA_70_74m <- as.vector(DHA.EPA_70_74m$DHA.EPA)

#Age 75-79, male
DHA.EPA_75_79m <- subset(Scenario3_CHD, agegroups=="75-79" & sex=="1", select = DHA.EPA)
DHA.EPA_75_79m <- as.vector(DHA.EPA_75_79m$DHA.EPA)

###################################################### Exposure females #########################################################

#Prob=0 for women <35

#Age 35-39, female
DHA.EPA_35_39w <- subset(Scenario3_CHD, agegroups=="35-39" & sex=="2", select = DHA.EPA)
DHA.EPA_35_39w <- as.vector(DHA.EPA_35_39w$DHA.EPA)

#Age 40-44, female
DHA.EPA_40_44w <- subset(Scenario3_CHD, agegroups=="40-44" & sex=="2", select = DHA.EPA)
DHA.EPA_40_44w <- as.vector(DHA.EPA_40_44w$DHA.EPA)

#Age 45-49, female
DHA.EPA_45_49w <- subset(Scenario3_CHD, agegroups=="45-49" & sex=="2", select = DHA.EPA)
DHA.EPA_45_49w <- as.vector(DHA.EPA_45_49w$DHA.EPA)

#Age 50-54, female
DHA.EPA_50_54w <- subset(Scenario3_CHD, agegroups=="50-54" & sex=="2", select = DHA.EPA)
DHA.EPA_50_54w <- as.vector(DHA.EPA_50_54w$DHA.EPA)

#Age 55-59, female
DHA.EPA_55_59w <- subset(Scenario3_CHD, agegroups=="55-59" & sex=="2", select = DHA.EPA)
DHA.EPA_55_59w <- as.vector(DHA.EPA_55_59w$DHA.EPA)

#Age 60-64, female
DHA.EPA_60_64w <- subset(Scenario3_CHD, agegroups=="60-64" & sex=="2", select = DHA.EPA)
DHA.EPA_60_64w <- as.vector(DHA.EPA_60_64w$DHA.EPA)

#Age 65-69, female
DHA.EPA_65_69w <- subset(Scenario3_CHD, agegroups=="65-69" & sex=="2", select = DHA.EPA)
DHA.EPA_65_69w <- as.vector(DHA.EPA_65_69w$DHA.EPA)

#Age 70-74, female
DHA.EPA_70_74w <- subset(Scenario3_CHD, agegroups=="70-74" & sex=="2", select = DHA.EPA)
DHA.EPA_70_74w <- as.vector(DHA.EPA_70_74w$DHA.EPA)

#Age 75-79, female
DHA.EPA_75_79w <- subset(Scenario3_CHD, agegroups=="75-79" & sex=="2", select = DHA.EPA)
DHA.EPA_75_79w <- as.vector(DHA.EPA_75_79w$DHA.EPA)

##################################################### DALY scen 3 male #########################################################

#Prob = 0 for men, 15-19y

### Age 20-24, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_20_24m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect20_24m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_20_24m <- pEffect * pop_20_24m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_20_24scen3m <- pEffect*SEYLL_20_24
summary(DALY_20_24scen3m)


tDALY_20_24scen3m <- DALY_20_24scen3m * pop_20_24m 
summary(tDALY_20_24scen3m)


#Prob = 0 for men, 25-29y

### Age 30-34, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_30_34m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect30_34m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_30_34m <- pEffect * pop_30_34m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_30_34scen3m <- pEffect*SEYLL_30_34
summary(DALY_30_34scen3m)


tDALY_30_34scen3m <- DALY_30_34scen3m * pop_30_34m 
summary(tDALY_30_34scen3m)


### Age 35-39, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_35_39m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect35_39m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_35_39m <- pEffect * pop_35_39m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_35_39scen3m <- pEffect*SEYLL_35_39
summary(DALY_35_39scen3m)


tDALY_35_39scen3m <- DALY_35_39scen3m * pop_35_39m 
summary(tDALY_35_39scen3m)


### Age 40-44, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_40_44m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect40_44m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_40_44m <- pEffect * pop_40_44m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_40_44scen3m <- pEffect*SEYLL_40_44
summary(DALY_40_44scen3m)


tDALY_40_44scen3m <- DALY_40_44scen3m * pop_40_44m 
summary(tDALY_40_44scen3m)


### Age 45-49, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_45_49m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect45_49m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_45_49m <- pEffect * pop_45_49m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_45_49scen3m <- pEffect*SEYLL_45_49
summary(DALY_45_49scen3m)


tDALY_45_49scen3m <- DALY_45_49scen3m * pop_45_49m 
summary(tDALY_45_49scen3m)


### Age 50-54, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_50_54m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect50_54m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_50_54m <- pEffect * pop_50_54m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_50_54scen3m <- pEffect*SEYLL_50_54
summary(DALY_50_54scen3m)


tDALY_50_54scen3m <- DALY_50_54scen3m * pop_50_54m 
summary(tDALY_50_54scen3m)


### Age 55-59, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_55_59m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect55_59m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_55_59m <- pEffect * pop_55_59m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_55_59scen3m <- pEffect*SEYLL_55_59
summary(DALY_55_59scen3m)


tDALY_55_59scen3m <- DALY_55_59scen3m * pop_55_59m 
summary(tDALY_55_59scen3m)


### Age 60-64, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_60_64m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect60_64m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_60_64m <- pEffect * pop_60_64m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_60_64scen3m <- pEffect*SEYLL_60_64
summary(DALY_60_64scen3m)


tDALY_60_64scen3m <- DALY_60_64scen3m * pop_60_64m 
summary(tDALY_60_64scen3m)


### Age 65-69, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_65_69m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect65_69m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_65_69m <- pEffect * pop_65_69m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_65_69scen3m <- pEffect*SEYLL_65_69
summary(DALY_65_69scen3m)


tDALY_65_69scen3m <- DALY_65_69scen3m * pop_65_69m 
summary(tDALY_65_69scen3m)


### Age 70-74, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_70_74m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect70_74m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_70_74m <- pEffect * pop_70_74m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_70_74scen3m <- pEffect*SEYLL_70_74
summary(DALY_70_74scen3m)


tDALY_70_74scen3m <- DALY_70_74scen3m * pop_70_74m 
summary(tDALY_70_74scen3m)


### Age 75-79, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect75_79m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_75_79m <- pEffect * pop_75_79m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_75_79scen3m <- pEffect*SEYLL_75_79
summary(DALY_75_79scen3m)


tDALY_75_79scen3m <- DALY_75_79scen3m * pop_75_79m 
summary(tDALY_75_79scen3m)



### Age 80-84, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79m) #Extrapolate intake to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect80_84m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_80_84m <- pEffect * pop_80_84m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_80_84scen3m <- pEffect*SEYLL_80_84
summary(DALY_80_84scen3m)


tDALY_80_84scen3m <- DALY_80_84scen3m * pop_80_84m 
summary(tDALY_80_84scen3m)



### Age 85+, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79m) #Extrapolate intake to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect85m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_85m <- pEffect * pop_85m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_85scen3m <- pEffect*SEYLL_85
summary(DALY_85scen3m)


tDALY_85scen3m <- DALY_85scen3m * pop_85m 
summary(tDALY_85scen3m)


################################################ DALY scen 3 females ########################################################

#Prob = 0 for women < 35

### Age 35-39, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_35_39w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect35_39w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_35_39w <- pEffect * pop_35_39w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_35_39scen3w <- pEffect*SEYLL_35_39
summary(DALY_35_39scen3w)


tDALY_35_39scen3w <- DALY_35_39scen3w * pop_35_39w
summary(tDALY_35_39scen3w)


### Age 40-44, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_40_44w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect40_44w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_40_44w <- pEffect * pop_40_44w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_40_44scen3w <- pEffect*SEYLL_40_44
summary(DALY_40_44scen3w)


tDALY_40_44scen3w <- DALY_40_44scen3w * pop_40_44w
summary(tDALY_40_44scen3w)


### Age 45-49, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_45_49w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect45_49w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_45_49w <- pEffect * pop_45_49w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_45_49scen3w <- pEffect*SEYLL_45_49
summary(DALY_45_49scen3w)


tDALY_45_49scen3w <- DALY_45_49scen3w * pop_45_49w
summary(tDALY_45_49scen3w)


### Age 50-54, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_50_54w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect50_54w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_50_54w <- pEffect * pop_50_54w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_50_54scen3w <- pEffect*SEYLL_50_54
summary(DALY_50_54scen3w)


tDALY_50_54scen3w <- DALY_50_54scen3w * pop_50_54w
summary(tDALY_50_54scen3w)


### Age 55-59, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_55_59w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect55_59w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_55_59w <- pEffect * pop_55_59w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_55_59scen3w <- pEffect*SEYLL_55_59
summary(DALY_55_59scen3w)


tDALY_55_59scen3w <- DALY_55_59scen3w * pop_55_59w
summary(tDALY_55_59scen3w)


### Age 60-64, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_60_64w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect60_64w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_60_64w <- pEffect * pop_60_64w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_60_64scen3w <- pEffect*SEYLL_60_64
summary(DALY_60_64scen3w)


tDALY_60_64scen3w <- DALY_60_64scen3w * pop_60_64w
summary(tDALY_60_64scen3w)


### Age 65-69, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_65_69w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect65_69w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_65_69w <- pEffect * pop_65_69w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_65_69scen3w <- pEffect*SEYLL_65_69
summary(DALY_65_69scen3w)


tDALY_65_69scen3w <- DALY_65_69scen3w * pop_65_69w
summary(tDALY_65_69scen3w)



### Age 70-74, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_70_74w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect70_74w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_70_74w <- pEffect * pop_70_74w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_70_74scen3w <- pEffect*SEYLL_70_74
summary(DALY_70_74scen3w)


tDALY_70_74scen3w <- DALY_70_74scen3w * pop_70_74w
summary(tDALY_70_74scen3w)


### Age 75-79, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect75_79w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_75_79w <- pEffect * pop_75_79w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_75_79scen3w <- pEffect*SEYLL_75_79
summary(DALY_75_79scen3w)

tDALY_75_79scen3w <- DALY_75_79scen3w * pop_75_79w
summary(tDALY_75_79scen3w)

### Age 80-84, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79w) #Extrapolate intakes to > 75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect80_84w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_80_84w <- pEffect * pop_80_84w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_80_84scen3w <- pEffect*SEYLL_80_84
summary(DALY_80_84scen3w)

tDALY_80_84scen3w <- DALY_80_84scen3w * pop_80_84w
summary(tDALY_80_84scen3w)


### Age 85+, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79w) #Extrapolate intakes to > 75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)

pEffect <- pEffect85w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_85w <- pEffect * pop_85w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_85scen3w <- pEffect*SEYLL_85
summary(DALY_85scen3w)


tDALY_85scen3w <- DALY_85scen3w * pop_85w
summary(tDALY_85scen3w)


##### Total DALY #####

DALYtotal.CHD.scen3 <- tDALY_20_24scen3m + tDALY_30_34scen3m +tDALY_35_39scen3m + tDALY_40_44scen3m + tDALY_45_49scen3m + tDALY_45_49scen3m + tDALY_50_54scen3m + tDALY_55_59scen3m +
  tDALY_60_64scen3m + tDALY_65_69scen3m + tDALY_70_74scen3m + tDALY_75_79scen3m + tDALY_80_84scen3m + tDALY_85scen3m +
  tDALY_35_39scen3w + tDALY_40_44scen3w + tDALY_45_49scen3w + tDALY_45_49scen3w + tDALY_50_54scen3w + tDALY_55_59scen3w +
  tDALY_60_64scen3w + tDALY_65_69scen3w + tDALY_70_74scen3w + tDALY_75_79scen3w + tDALY_80_84scen3w + tDALY_85scen3w

summary(DALYtotal.CHD.scen3)


##### Cases scenario 3 #####

cases_scen3 <- cases_20_24m + cases_30_34m + cases_35_39m + cases_40_44m + cases_45_49m + cases_50_54m + cases_55_59m +
  cases_60_64m + cases_65_69m + cases_70_74m + cases_75_79m + cases_80_84m + cases_85m + 
  cases_35_39w + cases_40_44w + cases_45_49w + cases_50_54w + cases_55_59w + cases_60_64w + cases_65_69w +
  cases_70_74w + cases_75_79w + cases_80_84w + cases_85w

cases_diff_scen3 <- cases_scen3 - cases_ref
summary(cases_diff_scen3)


##### Scenario 4 #####

Scenario4 <- read.csv("Scenario4.csv")
# fish intakes are raw weights, ready for exposure assessment.

#Calculate DHA exposure from fish (mg/day)
Scen4_DHA_EPA <- t(Scenario4[,55]) #Create data set only with intakes of the different fish species
Scen4_DHA_EPA <- t(Scen4_DHA_EPA * DHA_EPA_fish[8,]) ##Multiply intakes of fish species with MeHg conc of the different species to obtain dataset with MeHg exposures from the different species

Scenario4$DHA.EPA <- rowSums(Scen4_DHA_EPA) #Add sum of DHA exposures for each individual to the Scenario4.DHA dataset
Scenario4_CHD <- Scenario4[,c(1:4,54,56)]

#Divide individuals according to age groups and sex

setDT(Scenario4_CHD)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]



##### Exposure males #####

#Divide data into gender- and agegroups

#Prob=0 for men, 15-19y

#Age 20-24, male
DHA.EPA_20_24m <- subset(Scenario4_CHD, agegroups=="20-24" & sex=="1", select = DHA.EPA)
DHA.EPA_20_24m <- as.vector(DHA.EPA_20_24m$DHA.EPA)

#Prob=0 for men, 25-29y

#Age 30-34, male
DHA.EPA_30_34m <- subset(Scenario4_CHD, agegroups=="30-34" & sex=="1", select = DHA.EPA)
DHA.EPA_30_34m <- as.vector(DHA.EPA_30_34m$DHA.EPA)

#Age 35-39, male
DHA.EPA_35_39m <- subset(Scenario4_CHD, agegroups=="35-39" & sex=="1", select = DHA.EPA)
DHA.EPA_35_39m <- as.vector(DHA.EPA_35_39m$DHA.EPA)

#Age 40-44, male
DHA.EPA_40_44m <- subset(Scenario4_CHD, agegroups=="40-44" & sex=="1", select = DHA.EPA)
DHA.EPA_40_44m <- as.vector(DHA.EPA_40_44m$DHA.EPA)

#Age 45-49, male
DHA.EPA_45_49m <- subset(Scenario4_CHD, agegroups=="45-49" & sex=="1", select = DHA.EPA)
DHA.EPA_45_49m <- as.vector(DHA.EPA_45_49m$DHA.EPA)

#Age 50-54, male
DHA.EPA_50_54m <- subset(Scenario4_CHD, agegroups=="50-54" & sex=="1", select = DHA.EPA)
DHA.EPA_50_54m <- as.vector(DHA.EPA_50_54m$DHA.EPA)

#Age 55-59, male
DHA.EPA_55_59m <- subset(Scenario4_CHD, agegroups=="55-59" & sex=="1", select = DHA.EPA)
DHA.EPA_55_59m <- as.vector(DHA.EPA_55_59m$DHA.EPA)

#Age 60-64, male
DHA.EPA_60_64m <- subset(Scenario4_CHD, agegroups=="60-64" & sex=="1", select = DHA.EPA)
DHA.EPA_60_64m <- as.vector(DHA.EPA_60_64m$DHA.EPA)

#Age 65-69, male
DHA.EPA_65_69m <- subset(Scenario4_CHD, agegroups=="65-69" & sex=="1", select = DHA.EPA)
DHA.EPA_65_69m <- as.vector(DHA.EPA_65_69m$DHA.EPA)

#Age 70-74, male
DHA.EPA_70_74m <- subset(Scenario4_CHD, agegroups=="70-74" & sex=="1", select = DHA.EPA)
DHA.EPA_70_74m <- as.vector(DHA.EPA_70_74m$DHA.EPA)

#Age 75-79, male
DHA.EPA_75_79m <- subset(Scenario4_CHD, agegroups=="75-79" & sex=="1", select = DHA.EPA)
DHA.EPA_75_79m <- as.vector(DHA.EPA_75_79m$DHA.EPA)


###### Exposure females #####

#Prob=0 for women<35

#Age 35-39, female
DHA.EPA_35_39w <- subset(Scenario4_CHD, agegroups=="35-39" & sex=="2", select = DHA.EPA)
DHA.EPA_35_39w <- as.vector(DHA.EPA_35_39w$DHA.EPA)

#Age 40-44, female
DHA.EPA_40_44w <- subset(Scenario4_CHD, agegroups=="40-44" & sex=="2", select = DHA.EPA)
DHA.EPA_40_44w <- as.vector(DHA.EPA_40_44w$DHA.EPA)

#Age 45-49, female
DHA.EPA_45_49w <- subset(Scenario4_CHD, agegroups=="45-49" & sex=="2", select = DHA.EPA)
DHA.EPA_45_49w <- as.vector(DHA.EPA_45_49w$DHA.EPA)

#Age 50-54, female
DHA.EPA_50_54w <- subset(Scenario4_CHD, agegroups=="50-54" & sex=="2", select = DHA.EPA)
DHA.EPA_50_54w <- as.vector(DHA.EPA_50_54w$DHA.EPA)

#Age 55-59, female
DHA.EPA_55_59w <- subset(Scenario4_CHD, agegroups=="55-59" & sex=="2", select = DHA.EPA)
DHA.EPA_55_59w <- as.vector(DHA.EPA_55_59w$DHA.EPA)

#Age 60-64, female
DHA.EPA_60_64w <- subset(Scenario4_CHD, agegroups=="60-64" & sex=="2", select = DHA.EPA)
DHA.EPA_60_64w <- as.vector(DHA.EPA_60_64w$DHA.EPA)

#Age 65-69, female
DHA.EPA_65_69w <- subset(Scenario4_CHD, agegroups=="65-69" & sex=="2", select = DHA.EPA)
DHA.EPA_65_69w <- as.vector(DHA.EPA_65_69w$DHA.EPA)

#Age 70-74, female
DHA.EPA_70_74w <- subset(Scenario4_CHD, agegroups=="70-74" & sex=="2", select = DHA.EPA)
DHA.EPA_70_74w <- as.vector(DHA.EPA_70_74w$DHA.EPA)

#Age 75-79, female
DHA.EPA_75_79w <- subset(Scenario4_CHD, agegroups=="75-79" & sex=="2", select = DHA.EPA)
DHA.EPA_75_79w <- as.vector(DHA.EPA_75_79w$DHA.EPA)


###### DALY scen 4 male ######

#Prob = 0 for men, 15-19y

### Age 20-24, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_20_24m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect20_24m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_20_24m <- pEffect * pop_20_24m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_20_24scen4m <- pEffect*SEYLL_20_24
summary(DALY_20_24scen4m)

tDALY_20_24scen4m <- DALY_20_24scen4m * pop_20_24m 
summary(tDALY_20_24scen4m)

#Prob = 0 for men, 25-29y

### Age 30-34, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_30_34m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect30_34m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_30_34m <- pEffect * pop_30_34m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_30_34scen4m <- pEffect*SEYLL_30_34
summary(DALY_30_34scen4m)


tDALY_30_34scen4m <- DALY_30_34scen4m * pop_30_34m 
summary(tDALY_30_34scen4m)


### Age 35-39, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_35_39m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect35_39m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_35_39m <- pEffect * pop_35_39m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_35_39scen4m <- pEffect*SEYLL_35_39
summary(DALY_35_39scen4m)


tDALY_35_39scen4m <- DALY_35_39scen4m * pop_35_39m 
summary(tDALY_35_39scen4m)


### Age 40-44, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_40_44m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect40_44m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_40_44m <- pEffect * pop_40_44m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_40_44scen4m <- pEffect*SEYLL_40_44
summary(DALY_40_44scen4m)


tDALY_40_44scen4m <- DALY_40_44scen4m * pop_40_44m 
summary(tDALY_40_44scen4m)


### Age 45-49, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_45_49m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect45_49m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_45_49m <- pEffect * pop_45_49m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_45_49scen4m <- pEffect*SEYLL_45_49
summary(DALY_45_49scen4m)


tDALY_45_49scen4m <- DALY_45_49scen4m * pop_45_49m 
summary(tDALY_45_49scen4m)


### Age 50-54, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_50_54m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect50_54m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_50_54m <- pEffect * pop_50_54m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_50_54scen4m <- pEffect*SEYLL_50_54
summary(DALY_50_54scen4m)


tDALY_50_54scen4m <- DALY_50_54scen4m * pop_50_54m 
summary(tDALY_50_54scen4m)


### Age 55-59, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_55_59m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect55_59m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_55_59m <- pEffect * pop_55_59m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_55_59scen4m <- pEffect*SEYLL_55_59
summary(DALY_55_59scen4m)


tDALY_55_59scen4m <- DALY_55_59scen4m * pop_55_59m 
summary(tDALY_55_59scen4m)


### Age 60-64, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_60_64m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect60_64m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_60_64m <- pEffect * pop_60_64m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_60_64scen4m <- pEffect*SEYLL_60_64
summary(DALY_60_64scen4m)


tDALY_60_64scen4m <- DALY_60_64scen4m * pop_60_64m 
summary(tDALY_60_64scen4m)


### Age 65-69, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_65_69m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect65_69m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_65_69m <- pEffect * pop_65_69m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_65_69scen4m <- pEffect*SEYLL_65_69
summary(DALY_65_69scen4m)


tDALY_65_69scen4m <- DALY_65_69scen4m * pop_65_69m 
summary(tDALY_65_69scen4m)


### Age 70-74, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_70_74m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)

pEffect <- pEffect70_74m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_70_74m <- pEffect * pop_70_74m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_70_74scen4m <- pEffect*SEYLL_70_74
summary(DALY_70_74scen4m)


tDALY_70_74scen4m <- DALY_70_74scen4m * pop_70_74m 
summary(tDALY_70_74scen4m)


### Age 75-79, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79m)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect75_79m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_75_79m <- pEffect * pop_75_79m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_75_79scen4m <- pEffect*SEYLL_75_79
summary(DALY_75_79scen4m)


tDALY_75_79scen4m <- DALY_75_79scen4m * pop_75_79m 
summary(tDALY_75_79scen4m)


### Age 80-84, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79m) #Extrapolate intakes to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect80_84m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_80_84m <- pEffect * pop_80_84m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_80_84scen4m <- pEffect*SEYLL_80_84
summary(DALY_80_84scen4m)


tDALY_80_84scen4m <- DALY_80_84scen4m * pop_80_84m 
summary(tDALY_80_84scen4m)




### Age 85+, males ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79m) #Extrapolate intakes to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect85m * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_85m <- pEffect * pop_85m


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_85scen4m <- pEffect*SEYLL_85
summary(DALY_85scen4m)


tDALY_85scen4m <- DALY_85scen4m * pop_85m 
summary(tDALY_85scen4m)



##### DALY scen 4 females #####

#Prob = 0 for women < 35

### Age 35-39, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_35_39w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect35_39w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_35_39w <- pEffect * pop_35_39w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_35_39scen4w <- pEffect*SEYLL_35_39
summary(DALY_35_39scen4w)


tDALY_35_39scen4w <- DALY_35_39scen4w * pop_35_39w
summary(tDALY_35_39scen4w)


### Age 40-44, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_40_44w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect40_44w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_40_44w <- pEffect * pop_40_44w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_40_44scen4w <- pEffect*SEYLL_40_44
summary(DALY_40_44scen4w)


tDALY_40_44scen4w <- DALY_40_44scen4w * pop_40_44w
summary(tDALY_40_44scen4w)


### Age 45-49, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_45_49w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)



pEffect <- pEffect45_49w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_45_49w <- pEffect * pop_45_49w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_45_49scen4w <- pEffect*SEYLL_45_49
summary(DALY_45_49scen4w)


tDALY_45_49scen4w <- DALY_45_49scen4w * pop_45_49w
summary(tDALY_45_49scen4w)


### Age 50-54, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_50_54w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect50_54w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_50_54w <- pEffect * pop_50_54w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_50_54scen4w <- pEffect*SEYLL_50_54
summary(DALY_50_54scen4w)


tDALY_50_54scen4w <- DALY_50_54scen4w * pop_50_54w
summary(tDALY_50_54scen4w)


### Age 55-59, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_55_59w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect55_59w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_55_59w <- pEffect * pop_55_59w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_55_59scen4w <- pEffect*SEYLL_55_59
summary(DALY_55_59scen4w)


tDALY_55_59scen4w <- DALY_55_59scen4w * pop_55_59w
summary(tDALY_55_59scen4w)


### Age 60-64, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_60_64w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect60_64w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_60_64w <- pEffect * pop_60_64w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_60_64scen4w <- pEffect*SEYLL_60_64
summary(DALY_60_64scen4w)


tDALY_60_64scen4w <- DALY_60_64scen4w * pop_60_64w
summary(tDALY_60_64scen4w)


### Age 65-69, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_65_69w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect65_69w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_65_69w <- pEffect * pop_65_69w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_65_69scen4w <- pEffect*SEYLL_65_69
summary(DALY_65_69scen4w)


tDALY_65_69scen4w <- DALY_65_69scen4w * pop_65_69w
summary(tDALY_65_69scen4w)



### Age 70-74, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_70_74w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)



pEffect <- pEffect70_74w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_70_74w <- pEffect * pop_70_74w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_70_74scen4w <- pEffect*SEYLL_70_74
summary(DALY_70_74scen4w)


tDALY_70_74scen4w <- DALY_70_74scen4w * pop_70_74w
summary(tDALY_70_74scen4w)


### Age 75-79, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79w)
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect75_79w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_75_79w <- pEffect * pop_75_79w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_75_79scen4w <- pEffect*SEYLL_75_79
summary(DALY_75_79scen4w)


tDALY_75_79scen4w <- DALY_75_79scen4w * pop_75_79w
summary(tDALY_75_79scen4w)




### Age 80-84, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79w) #Extrapolate intakes to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect80_84w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_80_84w <- pEffect * pop_80_84w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_80_84scen4w <- pEffect*SEYLL_80_84
summary(DALY_80_84scen4w)


tDALY_80_84scen4w <- DALY_80_84scen4w * pop_80_84w
summary(tDALY_80_84scen4w)



### Age 85+, females ###

set.seed(1)
exp <- mcstoc(rempiricalD, type = "V", values = DHA.EPA_75_79w) #Extrapolate intakes to >75y
exp <- unmc(exp)
exp_trunc <- numeric(iters_var)
for(i in 1:iters_var){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<250, x, 250)
  exp_trunc[i] <- y
}
exp <- mcstoc(rempiricalD, type = "V", values = exp_trunc)

RR <- exp * r + 1
summary(RR)


pEffect <- pEffect85w * RR # new absolute risk: Peffect(ref)/RR(ref) * RR(alt)
summary(pEffect)

cases_85w <- pEffect * pop_85w


#Prec = 0,  Pdie = 1, YLDdie = 0 therefore YLDdie*w = 0
#Final DALY equation: DALY = pEffect*(LE(a)-CA)
#where LE(a) is the life expectancy at age a, and CA is the age of onset

DALY_85scen4w <- pEffect*SEYLL_85
summary(DALY_85scen4w)


tDALY_85scen4w <- DALY_85scen4w * pop_85w
summary(tDALY_85scen4w)


##### Total DALY #####

DALYtotal.CHD.scen4 <- tDALY_20_24scen4m + tDALY_30_34scen4m + tDALY_35_39scen4m + tDALY_40_44scen4m + tDALY_45_49scen4m + tDALY_45_49scen4m + tDALY_50_54scen4m + tDALY_55_59scen4m +
  tDALY_60_64scen4m + tDALY_65_69scen4m + tDALY_70_74scen4m + tDALY_75_79scen4m + tDALY_80_84scen4m + tDALY_85scen4m +
  tDALY_35_39scen4w + tDALY_40_44scen4w + tDALY_45_49scen4w + tDALY_45_49scen4w + tDALY_50_54scen4w + tDALY_55_59scen4w +
  tDALY_60_64scen4w + tDALY_65_69scen4w + tDALY_70_74scen4w + tDALY_75_79scen4w + tDALY_80_84scen4w + tDALY_85scen4w

summary(DALYtotal.CHD.scen4)


##### Cases scenario 4 #####

cases_scen4 <- cases_20_24m + cases_30_34m + cases_35_39m + cases_40_44m + cases_45_49m + cases_50_54m + cases_55_59m +
  cases_60_64m + cases_65_69m + cases_70_74m + cases_75_79m + cases_80_84m + cases_85m + 
  cases_35_39w + cases_40_44w + cases_45_49w + cases_50_54w + cases_55_59w + cases_60_64w + cases_65_69w +
  cases_70_74w + cases_75_79w + cases_80_84w + cases_85w

cases_diff_scen4 <- cases_scen4 - cases_ref
summary(cases_diff_scen4)



###### Scenario comparison #####

##Total DALYs

tDALY_ref_fish_CHD <- DALYtotal.CHD
summary(tDALY_ref_fish_CHD)


tDALY_scen1_fish_CHD <- DALYtotal.CHD.scen1
summary(tDALY_scen1_fish_CHD)


tDALY_scen2_fish_CHD <- DALYtotal.CHD.scen2
summary(tDALY_scen2_fish_CHD)


tDALY_scen3_fish_CHD <- DALYtotal.CHD.scen3
summary(tDALY_scen3_fish_CHD)


tDALY_scen4_fish_CHD <- DALYtotal.CHD.scen4
summary(tDALY_scen4_fish_CHD)


##Per 100,000

pop15_85_plus <- 4697068 #Danish population >= 15 y

tDALY_ref_fish_CHD_100000 <- tDALY_ref_fish_CHD/pop15_85_plus*1e+05
summary(tDALY_ref_fish_CHD_100000)


tDALY_scen1_fish_CHD_100000 <- tDALY_scen1_fish_CHD/pop15_85_plus*1e+05
summary(tDALY_scen1_fish_CHD_100000)


tDALY_scen2_fish_CHD_100000 <- tDALY_scen2_fish_CHD/pop15_85_plus*1e+05
summary(tDALY_scen2_fish_CHD_100000)


tDALY_scen3_fish_CHD_100000 <- tDALY_scen3_fish_CHD/pop15_85_plus*1e+05
summary(tDALY_scen3_fish_CHD_100000)


tDALY_scen4_fish_CHD_100000 <- tDALY_scen4_fish_CHD/pop15_85_plus*1e+05
summary(tDALY_scen4_fish_CHD_100000)


##Total DALY difference

dDALY_scen1_fish_CHD <- tDALY_scen1_fish_CHD - tDALY_ref_fish_CHD
summary(dDALY_scen1_fish_CHD)


dDALY_scen2_fish_CHD <- tDALY_scen2_fish_CHD - tDALY_ref_fish_CHD
summary(dDALY_scen2_fish_CHD)


dDALY_scen3_fish_CHD <- tDALY_scen3_fish_CHD - tDALY_ref_fish_CHD
summary(dDALY_scen3_fish_CHD)


dDALY_scen4_fish_CHD <- tDALY_scen4_fish_CHD - tDALY_ref_fish_CHD
summary(dDALY_scen4_fish_CHD)


##DALY difference per 100,000

dDALY_scen1_fish_CHD_100000 <- tDALY_scen1_fish_CHD_100000 - tDALY_ref_fish_CHD_100000
summary(dDALY_scen1_fish_CHD_100000)


dDALY_scen2_fish_CHD_100000 <- tDALY_scen2_fish_CHD_100000 - tDALY_ref_fish_CHD_100000
summary(dDALY_scen2_fish_CHD_100000)


dDALY_scen3_fish_CHD_100000 <- tDALY_scen3_fish_CHD_100000 - tDALY_ref_fish_CHD_100000
summary(dDALY_scen3_fish_CHD_100000)


dDALY_scen4_fish_CHD_100000 <- tDALY_scen4_fish_CHD_100000 - tDALY_ref_fish_CHD_100000
summary(dDALY_scen4_fish_CHD_100000)


##### Save #####

tDALY_ref_fish_CHD_matrix <- tDALY_ref_fish_CHD
write.csv(tDALY_ref_fish_CHD_matrix, "tDALY_ref_fish_CHD.csv")

tDALY_scen1_fish_CHD_matrix <- apply(tDALY_scen1_fish_CHD, 2, function(x) x*1) #get data out of mc2d
tDALY_scen1_fish_CHD_unc <- apply(tDALY_scen1_fish_CHD_matrix, 2, function(x) mean(x)) #only uncertainty dimension
write.csv(tDALY_scen1_fish_CHD_unc, "tDALY_scen1_fish_CHD_unc.csv")


tDALY_scen2_fish_CHD_matrix <- apply(tDALY_scen2_fish_CHD, 2, function(x) x*1)
tDALY_scen2_fish_CHD_unc <- apply(tDALY_scen2_fish_CHD_matrix, 2, function(x) mean(x))
write.csv(tDALY_scen2_fish_CHD_unc, "tDALY_scen2_fish_CHD_unc.csv")


tDALY_scen3_fish_CHD_matrix <- apply(tDALY_scen3_fish_CHD, 2, function(x) x*1)
tDALY_scen3_fish_CHD_unc <- apply(tDALY_scen3_fish_CHD_matrix, 2, function(x) mean(x))
write.csv(tDALY_scen3_fish_CHD_unc, "tDALY_scen3_fish_CHD_unc.csv")


tDALY_scen4_fish_CHD_matrix <- apply(tDALY_scen4_fish_CHD, 2, function(x) x*1)
tDALY_scen4_fish_CHD_unc <- apply(tDALY_scen4_fish_CHD_matrix, 2, function(x) mean(x))
write.csv(tDALY_scen4_fish_CHD_unc, "tDALY_scen4_fish_CHD_unc.csv")
