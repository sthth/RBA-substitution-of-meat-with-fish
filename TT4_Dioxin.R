memory.limit(56000)

##packages
library(mc2d)
library(data.table)
library(fitdistrplus)
library(goftest)


## settings
nvar <- 1e5
nunc <- 1e3
set.seed(1)


## helpers
mean_median_ci <-
  function(x) {
    c(mean = mean(x),
      median = median(x)
      quantile(x, probs = c(0.025, 0.975)))
  }



## Population statistics

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

pop15_85_plus <- 4697068



## IPRA 


#Assume that 20% decrease in TT4 cause hypothyroidism
#CED for 20% decrease in TT4 due to dioxin exposure


# critical effect dose - animals

#CES 0.20 Hill
CEDM_animal <- 15.91
CEDL_animal <- 2.908
CEDU_animal <- 49.86


set.seed(1)
CED_animal <- rpert(nunc, CEDL_animal, CEDM_animal, CEDU_animal)
hist(CED_animal)

# critical effect dose - humans
CED_human <- CED_animal / 163*1000 ##US EPA CF for 90 day exposure to convert animal CED into human CED and to convert ng to pg
hist(CED_human)

# extrapolation factor - intraspecies
EF_intra <- rlnorm(nvar, log(1), log(1.98))
hist(EF_intra)

str(CED_human) 
str(t(CED_human)) #1 row, 1000 columns (uncertainty)
head(CED_human)

## individual critical effect dose : ICED = CED(human)/EF(intra)
## .. each col = uncertainty simulation
## .. each row = variability simulation
## .. ICED[var, unc]
ICED <- apply(t(CED_human), 2, function(x) x / EF_intra)
str(ICED)


## DALY input parameters

# YLL = 0
# YLD = I * D * DW

# Assume treatment within one year D = 1

D <- 1

# Hypothyroidism DW: 0.019 (0.010-0.032), Salomon 2015

dw <- rpert(nunc, min=0.010, mode=0.019, max=0.032)



##### Reference scenario #####

Dioxin_exp <- read.csv("Dioxin_exposure.csv")


agebreaks <- c(15,20,25,30,35,40,45,50,55,60,65,70,75,80)
agelabels= c("15-19","20-24","25-29","30-34","35-39", "40-44", "45-49", "50-54","55-59", "60-64","65-69", "70-74","75-79")

setDT(Dioxin_exp)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Fitting exposures male #####


#Age 15-19, males
Dioxin_15_19m <- subset(Dioxin_exp, agegroups=="15-19" & sex=="1", select = totaldioxin_bw)
Dioxin_15_19m <- as.vector(Dioxin_15_19m$totaldioxin_bw)


fit1_15_19m <- fitdist(Dioxin_15_19m, "lnorm")
t <- Dioxin_15_19m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_15_19m$estimate[1], sdlog = fit1_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19m, fitnames="lnorm") 


fit2_15_19m <- fitdist(Dioxin_15_19m, 'gamma') 
t2 <- Dioxin_15_19m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_15_19m$estimate[1],fit2_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19m, fitnames="gamma")


cvm.test(t, plnorm, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2])  #Anderson-Darling test


#Age 20-24, males
Dioxin_20_24m <- subset(Dioxin_exp, agegroups=="20-24" & sex=="1", select = totaldioxin_bw)
Dioxin_20_24m <- as.vector(Dioxin_20_24m$totaldioxin_bw)

fit1_20_24m <- fitdist(Dioxin_20_24m, "lnorm")
t <- Dioxin_20_24m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_20_24m$estimate[1], sdlog = fit1_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24m, fitnames="lnorm") 


fit2_20_24m <- fitdist(Dioxin_20_24m, 'gamma') 
t2 <- Dioxin_20_24m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_20_24m$estimate[1],fit2_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24m, fitnames="gamma")

cvm.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Anderson-Darling test


#Age 25-29, males
Dioxin_25_29m <- subset(Dioxin_exp, agegroups=="25-29" & sex=="1", select = totaldioxin_bw)
Dioxin_25_29m <- as.vector(Dioxin_25_29m$totaldioxin_bw)


fit1_25_29m <- fitdist(Dioxin_25_29m, "lnorm")
t <- Dioxin_25_29m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_25_29m$estimate[1], sdlog = fit1_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29m, fitnames="lnorm") 


fit2_25_29m <- fitdist(Dioxin_25_29m, 'gamma') 
t2 <- Dioxin_25_29m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_25_29m$estimate[1],fit2_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29m, fitnames="gamma")


cvm.test(t, plnorm, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2])  #Anderson-Darling test


#Age 30-34, males
Dioxin_30_34m <- subset(Dioxin_exp, agegroups=="30-34" & sex=="1", select = totaldioxin_bw)
Dioxin_30_34m <- as.vector(Dioxin_30_34m$totaldioxin_bw)


fit1_30_34m <- fitdist(Dioxin_30_34m, "lnorm")
t <- Dioxin_30_34m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_30_34m$estimate[1], sdlog = fit1_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34m, fitnames="lnorm") 


fit2_30_34m <- fitdist(Dioxin_30_34m, 'gamma') 
t2 <- Dioxin_30_34m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_30_34m$estimate[1],fit2_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34m, fitnames="gamma")


cvm.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Anderson-Darling test


#Age 35-39, males
Dioxin_35_39m <- subset(Dioxin_exp, agegroups=="35-39" & sex=="1", select = totaldioxin_bw)
Dioxin_35_39m <- as.vector(Dioxin_35_39m$totaldioxin_bw)


fit1_35_39m <- fitdist(Dioxin_35_39m, "lnorm")
t <- Dioxin_35_39m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_35_39m$estimate[1], sdlog = fit1_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39m, fitnames="lnorm") 


fit2_35_39m <- fitdist(Dioxin_35_39m, 'gamma') 
t2 <- Dioxin_35_39m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_35_39m$estimate[1],fit2_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39m, fitnames="gamma")


cvm.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Anderson-Darling test


#Age 40-44, males
Dioxin_40_44m <- subset(Dioxin_exp, agegroups=="40-44" & sex=="1", select = totaldioxin_bw)
Dioxin_40_44m <- as.vector(Dioxin_40_44m$totaldioxin_bw)


fit1_40_44m <- fitdist(Dioxin_40_44m, "lnorm")
t <- Dioxin_40_44m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_40_44m$estimate[1], sdlog = fit1_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44m, fitnames="lnorm") 


fit2_40_44m <- fitdist(Dioxin_40_44m, 'gamma') 
t2 <- Dioxin_40_44m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_40_44m$estimate[1],fit2_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44m, fitnames="gamma")


cvm.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Anderson-Darling test


#Age 45-49, males
Dioxin_45_49m <- subset(Dioxin_exp, agegroups=="45-49" & sex=="1", select = totaldioxin_bw)
Dioxin_45_49m <- as.vector(Dioxin_45_49m$totaldioxin_bw)


fit1_45_49m <- fitdist(Dioxin_45_49m, "lnorm")
t <- Dioxin_45_49m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_45_49m$estimate[1], sdlog = fit1_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49m, fitnames="lnorm") 


fit2_45_49m <- fitdist(Dioxin_45_49m, 'gamma') 
t2 <- Dioxin_45_49m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_45_49m$estimate[1],fit2_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49m, fitnames="gamma")


cvm.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Anderson-Darling test


#Age 50-54, males
Dioxin_50_54m <- subset(Dioxin_exp, agegroups=="50-54" & sex=="1", select = totaldioxin_bw)
Dioxin_50_54m <- as.vector(Dioxin_50_54m$totaldioxin_bw)


fit1_50_54m <- fitdist(Dioxin_50_54m, "lnorm")
t <- Dioxin_50_54m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_50_54m$estimate[1], sdlog = fit1_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54m, fitnames="lnorm") 


fit2_50_54m <- fitdist(Dioxin_50_54m, 'gamma') 
t2 <- Dioxin_50_54m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_50_54m$estimate[1],fit2_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54m, fitnames="gamma")


cvm.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Anderson-Darling test


#Age 55-59, males
Dioxin_55_59m <- subset(Dioxin_exp, agegroups=="55-59" & sex=="1", select = totaldioxin_bw)
Dioxin_55_59m <- as.vector(Dioxin_55_59m$totaldioxin_bw)


fit1_55_59m <- fitdist(Dioxin_55_59m, "lnorm")
t <- Dioxin_55_59m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_55_59m$estimate[1], sdlog = fit1_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59m, fitnames="lnorm") 


fit2_55_59m <- fitdist(Dioxin_55_59m, 'gamma') 
t2 <- Dioxin_55_59m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_55_59m$estimate[1],fit2_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59m, fitnames="gamma")


cvm.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Anderson-Darling test



#Age 60-64, males
Dioxin_60_64m <- subset(Dioxin_exp, agegroups=="60-64" & sex=="1", select = totaldioxin_bw)
Dioxin_60_64m <- as.vector(Dioxin_60_64m$totaldioxin_bw)


fit1_60_64m <- fitdist(Dioxin_60_64m, "lnorm")
t <- Dioxin_60_64m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_60_64m$estimate[1], sdlog = fit1_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64m, fitnames="lnorm") 


fit2_60_64m <- fitdist(Dioxin_60_64m, 'gamma') 
t2 <- Dioxin_60_64m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_60_64m$estimate[1],fit2_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64m, fitnames="gamma")


cvm.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Anderson-Darling test


#Age 65-69, males
Dioxin_65_69m <- subset(Dioxin_exp, agegroups=="65-69" & sex=="1", select = totaldioxin_bw)
Dioxin_65_69m <- as.vector(Dioxin_65_69m$totaldioxin_bw)


fit1_65_69m <- fitdist(Dioxin_65_69m, "lnorm")
t <- Dioxin_65_69m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_65_69m$estimate[1], sdlog = fit1_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69m, fitnames="lnorm") 


fit2_65_69m <- fitdist(Dioxin_65_69m, 'gamma') 
t2 <- Dioxin_65_69m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_65_69m$estimate[1],fit2_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69m, fitnames="gamma")


cvm.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Anderson-Darling test



#Age 70-74, males
Dioxin_70_74m <- subset(Dioxin_exp, agegroups=="70-74" & sex=="1", select = totaldioxin_bw)
Dioxin_70_74m <- as.vector(Dioxin_70_74m$totaldioxin_bw)


fit1_70_74m <- fitdist(Dioxin_70_74m, "lnorm")
t <- Dioxin_70_74m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_70_74m$estimate[1], sdlog = fit1_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74m, fitnames="lnorm") 


fit2_70_74m <- fitdist(Dioxin_70_74m, 'gamma') 
t2 <- Dioxin_70_74m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_70_74m$estimate[1],fit2_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74m, fitnames="gamma")


cvm.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Anderson-Darling test



#Age 75-79, males
Dioxin_75_79m <- subset(Dioxin_exp, agegroups=="75-79" & sex=="1", select = totaldioxin_bw)
Dioxin_75_79m <- as.vector(Dioxin_75_79m$totaldioxin_bw)


fit1_75_79m <- fitdist(Dioxin_75_79m, "lnorm")
t <- Dioxin_75_79m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_75_79m$estimate[1], sdlog = fit1_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79m, fitnames="lnorm") 


fit2_75_79m <- fitdist(Dioxin_75_79m, 'gamma') 
t2 <- Dioxin_75_79m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_75_79m$estimate[1],fit2_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79m, fitnames="gamma")


cvm.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Anderson-Darling test



##### Fitting exposures female #####

#Age 15-19, female
Dioxin_15_19w <- subset(Dioxin_exp, agegroups=="15-19" & sex=="2", select = totaldioxin_bw)
Dioxin_15_19w <- as.vector(Dioxin_15_19w$totaldioxin_bw)

fit1_15_19w <- fitdist(Dioxin_15_19w, "lnorm")
t <- Dioxin_15_19w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_15_19w$estimate[1], sdlog = fit1_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19w, fitnames="lnorm") 


fit2_15_19w <- fitdist(Dioxin_15_19w, 'gamma') 
t2 <- Dioxin_15_19w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_15_19w$estimate[1],fit2_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19w, fitnames="gamma")

cvm.test(t, plnorm, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2])  #Anderson-Darling test


#Age 20-24, female
Dioxin_20_24w <- subset(Dioxin_exp, agegroups=="20-24" & sex=="2", select = totaldioxin_bw)
Dioxin_20_24w <- as.vector(Dioxin_20_24w$totaldioxin_bw)

fit1_20_24w <- fitdist(Dioxin_20_24w, "lnorm")
t <- Dioxin_20_24w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_20_24w$estimate[1], sdlog = fit1_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24w, fitnames="lnorm") 

fit2_20_24w <- fitdist(Dioxin_20_24w, 'gamma') 
t2 <- Dioxin_20_24w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_20_24w$estimate[1],fit2_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24w, fitnames="gamma")

cvm.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Anderson-Darling test


#Age 25-29, female
Dioxin_25_29w <- subset(Dioxin_exp, agegroups=="25-29" & sex=="2", select = totaldioxin_bw)
Dioxin_25_29w <- as.vector(Dioxin_25_29w$totaldioxin_bw)

fit1_25_29w <- fitdist(Dioxin_25_29w, "lnorm")
t <- Dioxin_25_29w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_25_29w$estimate[1], sdlog = fit1_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29w, fitnames="lnorm") 

fit2_25_29w <- fitdist(Dioxin_25_29w, 'gamma')
t2 <- Dioxin_25_29w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_25_29w$estimate[1],fit2_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29w, fitnames="gamma")

cvm.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Anderson-Darling test


#Age 30-34, female
Dioxin_30_34w <- subset(Dioxin_exp, agegroups=="30-34" & sex=="2", select = totaldioxin_bw)
Dioxin_30_34w <- as.vector(Dioxin_30_34w$totaldioxin_bw)

fit1_30_34w <- fitdist(Dioxin_30_34w, "lnorm")
t <- Dioxin_30_34w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_30_34w$estimate[1], sdlog = fit1_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34w, fitnames="lnorm") 

fit2_30_34w <- fitdist(Dioxin_30_34w, 'gamma') 
t2 <- Dioxin_30_34w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_30_34w$estimate[1],fit2_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34w, fitnames="gamma")

cvm.test(t, plnorm, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2])  #Anderson-Darling test


#Age 35-39, female
Dioxin_35_39w <- subset(Dioxin_exp, agegroups=="35-39" & sex=="2", select = totaldioxin_bw)
Dioxin_35_39w <- as.vector(Dioxin_35_39w$totaldioxin_bw)

fit1_35_39w <- fitdist(Dioxin_35_39w, "lnorm")
t <- Dioxin_35_39w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_35_39w$estimate[1], sdlog = fit1_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39w, fitnames="lnorm") 

fit2_35_39w <- fitdist(Dioxin_35_39w, 'gamma') 
t2 <- Dioxin_35_39w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_35_39w$estimate[1],fit2_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39w, fitnames="gamma")

cvm.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Anderson-Darling test


#Age 40-44, female
Dioxin_40_44w <- subset(Dioxin_exp, agegroups=="40-44" & sex=="2", select = totaldioxin_bw)
Dioxin_40_44w <- as.vector(Dioxin_40_44w$totaldioxin_bw)

fit1_40_44w <- fitdist(Dioxin_40_44w, "lnorm")
t <- Dioxin_40_44w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_40_44w$estimate[1], sdlog = fit1_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44w, fitnames="lnorm") 

fit2_40_44w <- fitdist(Dioxin_40_44w, 'gamma')
t2 <- Dioxin_40_44w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_40_44w$estimate[1],fit2_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44w, fitnames="gamma")

cvm.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Anderson-Darling test


#Age 45-49, female
Dioxin_45_49w <- subset(Dioxin_exp, agegroups=="45-49" & sex=="2", select = totaldioxin_bw)
Dioxin_45_49w <- as.vector(Dioxin_45_49w$totaldioxin_bw)

fit1_45_49w <- fitdist(Dioxin_45_49w, "lnorm")
t <- Dioxin_45_49w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_45_49w$estimate[1], sdlog = fit1_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49w, fitnames="lnorm") 

fit2_45_49w <- fitdist(Dioxin_45_49w, 'gamma')
t2 <- Dioxin_45_49w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_45_49w$estimate[1],fit2_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49w, fitnames="gamma")

cvm.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Anderson-Darling test


#Age 50-54, female
Dioxin_50_54w <- subset(Dioxin_exp, agegroups=="50-54" & sex=="2", select = totaldioxin_bw)
Dioxin_50_54w <- as.vector(Dioxin_50_54w$totaldioxin_bw)

fit1_50_54w <- fitdist(Dioxin_50_54w, "lnorm")
t <- Dioxin_50_54w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_50_54w$estimate[1], sdlog = fit1_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54w, fitnames="lnorm") 

fit2_50_54w <- fitdist(Dioxin_50_54w, 'gamma')
t2 <- Dioxin_50_54w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_50_54w$estimate[1],fit2_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54w, fitnames="gamma")

cvm.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Anderson-Darling test


#Age 55-59, female
Dioxin_55_59w <- subset(Dioxin_exp, agegroups=="55-59" & sex=="2", select = totaldioxin_bw)
Dioxin_55_59w <- as.vector(Dioxin_55_59w$totaldioxin_bw)

fit1_55_59w <- fitdist(Dioxin_55_59w, "lnorm")
t <- Dioxin_55_59w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_55_59w$estimate[1], sdlog = fit1_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59w, fitnames="lnorm") 

fit2_55_59w <- fitdist(Dioxin_55_59w, 'gamma')
t2 <- Dioxin_55_59w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_55_59w$estimate[1],fit2_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59w, fitnames="gamma")

cvm.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Anderson-Darling test



#Age 60-64, female
Dioxin_60_64w <- subset(Dioxin_exp, agegroups=="60-64" & sex=="2", select = totaldioxin_bw)
Dioxin_60_64w <- as.vector(Dioxin_60_64w$totaldioxin_bw)

fit1_60_64w <- fitdist(Dioxin_60_64w, "lnorm")
t <- Dioxin_60_64w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_60_64w$estimate[1], sdlog = fit1_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64w, fitnames="lnorm") 

fit2_60_64w <- fitdist(Dioxin_60_64w, 'gamma')
t2 <- Dioxin_60_64w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_60_64w$estimate[1],fit2_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64w, fitnames="gamma")

cvm.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Anderson-Darling test


#Age 65-69, female
Dioxin_65_69w <- subset(Dioxin_exp, agegroups=="65-69" & sex=="2", select = totaldioxin_bw)
Dioxin_65_69w <- as.vector(Dioxin_65_69w$totaldioxin_bw)

fit1_65_69w <- fitdist(Dioxin_65_69w, "lnorm")
t <- Dioxin_65_69w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_65_69w$estimate[1], sdlog = fit1_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69w, fitnames="lnorm") 

fit2_65_69w <- fitdist(Dioxin_65_69w, 'gamma')
t2 <- Dioxin_65_69w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_65_69w$estimate[1],fit2_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69w, fitnames="gamma")

cvm.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Anderson-Darling test


#Age 70-74, female
Dioxin_70_74w <- subset(Dioxin_exp, agegroups=="70-74" & sex=="2", select = totaldioxin_bw)
Dioxin_70_74w <- as.vector(Dioxin_70_74w$totaldioxin_bw)

fit1_70_74w <- fitdist(Dioxin_70_74w, "lnorm")
t <- Dioxin_70_74w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_70_74w$estimate[1], sdlog = fit1_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74w, fitnames="lnorm") 

fit2_70_74w <- fitdist(Dioxin_70_74w, 'gamma')
t2 <- Dioxin_70_74w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_70_74w$estimate[1],fit2_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74w, fitnames="gamma")

cvm.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Anderson-Darling test


#Age 75-79, female
Dioxin_75_79w <- subset(Dioxin_exp, agegroups=="75-79" & sex=="2", select = totaldioxin_bw)
Dioxin_75_79w <- as.vector(Dioxin_75_79w$totaldioxin_bw)

fit1_75_79w <- fitdist(Dioxin_75_79w, "lnorm")
t <- Dioxin_75_79w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_75_79w$estimate[1], sdlog = fit1_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79w, fitnames="lnorm") 

fit2_75_79w <- fitdist(Dioxin_75_79w, 'gamma')
t2 <- Dioxin_75_79w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_75_79w$estimate[1],fit2_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79w, fitnames="gamma")

cvm.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Anderson-Darling test


##### Prob of hypothyroidism males #####


# Assume exposure is variable - Lognormal has the best fit


# Age 15-19, males
set.seed(1)
exp <- rlnorm(nvar, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2]) 

pTT4_15_19m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 20-24, males
set.seed(1)
exp <- rlnorm(nvar, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2]) 

pTT4_20_24m <- apply(ICED, 2, function(x) mean(x < exp))


# Age 25-29, males
set.seed(1)
exp <- rlnorm(nvar, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2]) 

pTT4_25_29m <- apply(ICED, 2, function(x) mean(x < exp))


# Age 30-34, males
set.seed(1)
exp <- rlnorm(nvar, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2]) 

pTT4_30_34m <- apply(ICED, 2, function(x) mean(x < exp))


# Age 35-39, males
set.seed(1)
exp <- rlnorm(nvar, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2]) 

pTT4_35_39m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 40-44, males
set.seed(1)
exp <- rlnorm(nvar, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2]) 

pTT4_40_44m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 45-49, males
set.seed(1)
exp <- rlnorm(nvar, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2]) 

pTT4_45_49m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 50-54, males
set.seed(1)
exp <- rlnorm(nvar, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2]) 

pTT4_50_54m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 55-59, males
set.seed(1)
exp <- rlnorm(nvar, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2]) 

pTT4_55_59m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 60-64, males
set.seed(1)
exp <- rlnorm(nvar, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2]) 

pTT4_60_64m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 65-69, males
set.seed(1)
exp <- rlnorm(nvar, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2]) 

pTT4_65_69m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 70-74, males
set.seed(1)
exp <- rlnorm(nvar, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2]) 

pTT4_70_74m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 75-79, males
set.seed(1)
exp <- rlnorm(nvar, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2]) 

pTT4_75_79m <- apply(ICED, 2, function(x) mean(x < exp))


#Assume age 80-84 and 85+ have the same prob as 75-79



##### Prob of hypothyroidism females #####

# Assume exposure is variable

# Age 15-19, females
set.seed(1)
exp <- rlnorm(nvar, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2]) 

pTT4_15_19w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 20-24, females
set.seed(1)
exp <- rlnorm(nvar, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2]) 

pTT4_20_24w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 25-29, females
set.seed(1)
exp <- rlnorm(nvar, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2]) 

pTT4_25_29w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 30-34, females
set.seed(1)
exp <- rlnorm(nvar, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2]) 

pTT4_30_34w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 35-39, females
set.seed(1)
exp <- rlnorm(nvar, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2]) 

pTT4_35_39w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 40-44, females
set.seed(1)
exp <- rlnorm(nvar, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2]) 

pTT4_40_44w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 45-49, females
set.seed(1)
exp <- rlnorm(nvar, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2]) 

pTT4_45_49w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 50-54, females
set.seed(1)
exp <- rlnorm(nvar, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2]) 

pTT4_50_54w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 55-59, females
set.seed(1)
exp <- rlnorm(nvar, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2]) 

pTT4_55_59w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 60-64, females
set.seed(1)
exp <- rlnorm(nvar, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2]) 

pTT4_60_64w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 65-69, females
set.seed(1)
exp <- rlnorm(nvar, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2]) 

pTT4_65_69w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 70-74, females
set.seed(1)
exp <- rlnorm(nvar, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2]) 

pTT4_70_74w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 75-79, females
set.seed(1)
exp <- rlnorm(nvar, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2]) 

pTT4_75_79w <- apply(ICED, 2, function(x) mean(x < exp))



#Assume age 80-84 and 85+ have the same prob as 75-79


#####  DALY calculation  #####

##### Males #####

#Age 15-19, male

DALY15_19m <- pTT4_15_19m * pop_15_19m * D * dw
mean_median_ci(DALY15_19m)


#Age 20-24y, male

DALY20_24m <- pTT4_20_24m * pop_20_24m * D * dw
mean_median_ci(DALY20_24m)


#Age 25-29y, male

DALY25_29m <- pTT4_25_29m * pop_25_29m * D * dw
mean_median_ci(DALY25_29m)


#Age 30-34y, male

DALY30_34m <- pTT4_30_34m * pop_30_34m * D * dw
mean_median_ci(DALY30_34m)


#Age 35-39y, male

DALY35_39m <- pTT4_35_39m * pop_35_39m * D * dw
mean_median_ci(DALY35_39m)



#Age 40-44y, male

DALY40_44m <- pTT4_40_44m * pop_40_44m * D * dw
mean_median_ci(DALY40_44m)


#Age 45-49y, male

DALY45_49m <- pTT4_45_49m * pop_45_49m * D * dw
mean_median_ci(DALY45_49m)


#Age 50-54y, male

DALY50_54m <- pTT4_50_54m * pop_50_54m * D * dw
mean_median_ci(DALY50_54m)


#Age 55-59y, male

DALY55_59m <- pTT4_55_59m * pop_55_59m * D * dw
mean_median_ci(DALY55_59m)


#Age 60-64y, male

DALY60_64m <- pTT4_60_64m * pop_60_64m * D * dw
mean_median_ci(DALY60_64m)


#Age 65-69y, male

DALY65_69m <- pTT4_65_69m * pop_65_69m * D * dw
mean_median_ci(DALY65_69m)


#Age 70-74y, male

DALY70_74m <- pTT4_70_74m * pop_70_74m * D * dw
mean_median_ci(DALY70_74m)


#Age 75-79y, male

DALY75_79m <- pTT4_75_79m * pop_75_79m * D * dw
mean_median_ci(DALY75_79m)


#Age 80-84y, male

DALY80_84m <- pTT4_75_79m * pop_80_84m * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY80_84m)


#Age 85+y, male

DALY85m <- pTT4_75_79m * pop_85m * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY85m)


################################################## Females ###########################################################

#Age 15-19, female

DALY15_19w <- pTT4_15_19w * pop_15_19w * D * dw
mean_median_ci(DALY15_19w)


#Age 20-24y, female

DALY20_24w <- pTT4_20_24w * pop_20_24w * D * dw
mean_median_ci(DALY20_24w)


#Age 25-29y, female

DALY25_29w <- pTT4_25_29w * pop_25_29w * D * dw
mean_median_ci(DALY25_29w)


#Age 30-34y, female

DALY30_34w <- pTT4_30_34w * pop_30_34w * D * dw
mean_median_ci(DALY30_34w)


#Age 35-39y, female

DALY35_39w <- pTT4_35_39w * pop_35_39w * D * dw
mean_median_ci(DALY35_39w)



#Age 40-44y, female

DALY40_44w <- pTT4_40_44w * pop_40_44w * D * dw
mean_median_ci(DALY40_44w)


#Age 45-49y, female

DALY45_49w <- pTT4_45_49w * pop_45_49w * D * dw
mean_median_ci(DALY45_49w)


#Age 50-54y, female

DALY50_54w <- pTT4_50_54w * pop_50_54w * D * dw
mean_median_ci(DALY50_54w)


#Age 55-59y, female

DALY55_59w <- pTT4_55_59w * pop_55_59w * D * dw
mean_median_ci(DALY55_59w)


#Age 60-64y, female

DALY60_64w <- pTT4_60_64w * pop_60_64w * D * dw
mean_median_ci(DALY60_64w)


#Age 65-69y, female

DALY65_69w <- pTT4_65_69w * pop_65_69w * D * dw
mean_median_ci(DALY65_69w)


#Age 70-74y, female

DALY70_74w <- pTT4_70_74w * pop_70_74w * D * dw
mean_median_ci(DALY70_74w)


#Age 75-79y, female

DALY75_79w <- pTT4_75_79w * pop_75_79w * D * dw
mean_median_ci(DALY75_79w)


#Age 80-84y, female

DALY80_84w <- pTT4_75_79w * pop_80_84w * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY80_84w)


#Age 85+y, female

DALY85w <- pTT4_75_79w * pop_85w * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY85w)


##### Total DALYs ref #####

DALYtotal <- DALY15_19m + DALY20_24m + DALY25_29m + DALY30_34m + DALY35_39m + DALY40_44m + DALY45_49m + DALY50_54m +
  DALY55_59m + DALY60_64m + DALY65_69m + DALY70_74m + DALY75_79m + DALY80_84m + DALY85m +
  DALY15_19w + DALY20_24w + DALY25_29w + DALY30_34w + DALY35_39w + DALY40_44w + DALY45_49w + DALY50_54w + DALY55_59w +
  DALY60_64w + DALY65_69w + DALY70_74w + DALY75_79w + DALY80_84w + DALY85w


mean_median_ci(DALYtotal)


##### Cases reference scenario #####

cases_ref <- pTT4_15_19m * pop_15_19m + pTT4_20_24m * pop_20_24m + pTT4_25_29m * pop_25_29m + pTT4_30_34m * pop_30_34m +
  pTT4_35_39m * pop_35_39m + pTT4_40_44m * pop_40_44m + pTT4_45_49m * pop_45_49m + pTT4_50_54m * pop_50_54m + pTT4_55_59m * pop_55_59m +
  pTT4_60_64m * pop_60_64m + pTT4_65_69m * pop_65_69m + pTT4_70_74m * pop_70_74m + pTT4_75_79m * pop_75_79m + pTT4_75_79m * pop_80_84m +
  pTT4_75_79m * pop_85m +
  pTT4_15_19w * pop_15_19w + pTT4_20_24w * pop_20_24w + pTT4_25_29w * pop_25_29w + pTT4_30_34w * pop_30_34w +
  pTT4_35_39w * pop_35_39w + pTT4_40_44w * pop_40_44w + pTT4_45_49w * pop_45_49w + pTT4_50_54w * pop_50_54w + pTT4_55_59w * pop_55_59w +
  pTT4_60_64w * pop_60_64w + pTT4_65_69w * pop_65_69w + pTT4_70_74w * pop_70_74w + pTT4_75_79w * pop_75_79w + pTT4_75_79w * pop_80_84w +
  pTT4_75_79w * pop_85w
  
  
mean_median_ci(cases_ref)


##### Scenario 1 #####

Dioxin_exp1 <- read.csv("Dioxin_exp_scen1.csv")

setDT(Dioxin_exp1)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Fitting exposures male#####

#Age 15-19, males
Dioxin_15_19m <- subset(Dioxin_exp1, agegroups=="15-19" & sex=="1", select = totaldioxin_bw)
Dioxin_15_19m <- as.vector(Dioxin_15_19m$totaldioxin_bw)


fit1_15_19m <- fitdist(Dioxin_15_19m, "lnorm")
t <- Dioxin_15_19m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_15_19m$estimate[1], sdlog = fit1_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19m, fitnames="lnorm") 


fit2_15_19m <- fitdist(Dioxin_15_19m, 'gamma') 
t2 <- Dioxin_15_19m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_15_19m$estimate[1],fit2_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19m, fitnames="gamma")


cvm.test(t, plnorm, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2])  #Anderson-Darling test


#Age 20-24, males
Dioxin_20_24m <- subset(Dioxin_exp1, agegroups=="20-24" & sex=="1", select = totaldioxin_bw)
Dioxin_20_24m <- as.vector(Dioxin_20_24m$totaldioxin_bw)



fit1_20_24m <- fitdist(Dioxin_20_24m, "lnorm")
t <- Dioxin_20_24m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_20_24m$estimate[1], sdlog = fit1_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24m, fitnames="lnorm") 


fit2_20_24m <- fitdist(Dioxin_20_24m, 'gamma') 
t2 <- Dioxin_20_24m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_20_24m$estimate[1],fit2_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24m, fitnames="gamma")


cvm.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Anderson-Darling test


#Age 25-29, males
Dioxin_25_29m <- subset(Dioxin_exp1, agegroups=="25-29" & sex=="1", select = totaldioxin_bw)
Dioxin_25_29m <- as.vector(Dioxin_25_29m$totaldioxin_bw)



fit1_25_29m <- fitdist(Dioxin_25_29m, "lnorm")
t <- Dioxin_25_29m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_25_29m$estimate[1], sdlog = fit1_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29m, fitnames="lnorm") 


fit2_25_29m <- fitdist(Dioxin_25_29m, 'gamma') 
t2 <- Dioxin_25_29m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_25_29m$estimate[1],fit2_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29m, fitnames="gamma")


cvm.test(t, plnorm, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2])  #Anderson-Darling test


#Age 30-34, males
Dioxin_30_34m <- subset(Dioxin_exp1, agegroups=="30-34" & sex=="1", select = totaldioxin_bw)
Dioxin_30_34m <- as.vector(Dioxin_30_34m$totaldioxin_bw)



fit1_30_34m <- fitdist(Dioxin_30_34m, "lnorm")
t <- Dioxin_30_34m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_30_34m$estimate[1], sdlog = fit1_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34m, fitnames="lnorm") 


fit2_30_34m <- fitdist(Dioxin_30_34m, 'gamma') 
t2 <- Dioxin_30_34m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_30_34m$estimate[1],fit2_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34m, fitnames="gamma")


cvm.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Anderson-Darling test


#Age 35-39, males
Dioxin_35_39m <- subset(Dioxin_exp1, agegroups=="35-39" & sex=="1", select = totaldioxin_bw)
Dioxin_35_39m <- as.vector(Dioxin_35_39m$totaldioxin_bw)



fit1_35_39m <- fitdist(Dioxin_35_39m, "lnorm")
t <- Dioxin_35_39m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_35_39m$estimate[1], sdlog = fit1_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39m, fitnames="lnorm") 


fit2_35_39m <- fitdist(Dioxin_35_39m, 'gamma') 
t2 <- Dioxin_35_39m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_35_39m$estimate[1],fit2_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39m, fitnames="gamma")


cvm.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Anderson-Darling test


#Age 40-44, males
Dioxin_40_44m <- subset(Dioxin_exp1, agegroups=="40-44" & sex=="1", select = totaldioxin_bw)
Dioxin_40_44m <- as.vector(Dioxin_40_44m$totaldioxin_bw)



fit1_40_44m <- fitdist(Dioxin_40_44m, "lnorm")
t <- Dioxin_40_44m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_40_44m$estimate[1], sdlog = fit1_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44m, fitnames="lnorm") 


fit2_40_44m <- fitdist(Dioxin_40_44m, 'gamma') 
t2 <- Dioxin_40_44m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_40_44m$estimate[1],fit2_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44m, fitnames="gamma")


cvm.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Anderson-Darling test


#Age 45-49, males
Dioxin_45_49m <- subset(Dioxin_exp1, agegroups=="45-49" & sex=="1", select = totaldioxin_bw)
Dioxin_45_49m <- as.vector(Dioxin_45_49m$totaldioxin_bw)



fit1_45_49m <- fitdist(Dioxin_45_49m, "lnorm")
t <- Dioxin_45_49m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_45_49m$estimate[1], sdlog = fit1_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49m, fitnames="lnorm") 


fit2_45_49m <- fitdist(Dioxin_45_49m, 'gamma') 
t2 <- Dioxin_45_49m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_45_49m$estimate[1],fit2_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49m, fitnames="gamma")


cvm.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Anderson-Darling test


#Age 50-54, males
Dioxin_50_54m <- subset(Dioxin_exp1, agegroups=="50-54" & sex=="1", select = totaldioxin_bw)
Dioxin_50_54m <- as.vector(Dioxin_50_54m$totaldioxin_bw)



fit1_50_54m <- fitdist(Dioxin_50_54m, "lnorm")
t <- Dioxin_50_54m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_50_54m$estimate[1], sdlog = fit1_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54m, fitnames="lnorm") 


fit2_50_54m <- fitdist(Dioxin_50_54m, 'gamma') 
t2 <- Dioxin_50_54m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_50_54m$estimate[1],fit2_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54m, fitnames="gamma")


cvm.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Anderson-Darling test


#Age 55-59, males
Dioxin_55_59m <- subset(Dioxin_exp1, agegroups=="55-59" & sex=="1", select = totaldioxin_bw)
Dioxin_55_59m <- as.vector(Dioxin_55_59m$totaldioxin_bw)



fit1_55_59m <- fitdist(Dioxin_55_59m, "lnorm")
t <- Dioxin_55_59m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_55_59m$estimate[1], sdlog = fit1_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59m, fitnames="lnorm") 


fit2_55_59m <- fitdist(Dioxin_55_59m, 'gamma') 
t2 <- Dioxin_55_59m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_55_59m$estimate[1],fit2_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59m, fitnames="gamma")


cvm.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Anderson-Darling test



#Age 60-64, males
Dioxin_60_64m <- subset(Dioxin_exp1, agegroups=="60-64" & sex=="1", select = totaldioxin_bw)
Dioxin_60_64m <- as.vector(Dioxin_60_64m$totaldioxin_bw)



fit1_60_64m <- fitdist(Dioxin_60_64m, "lnorm")
t <- Dioxin_60_64m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_60_64m$estimate[1], sdlog = fit1_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64m, fitnames="lnorm") 


fit2_60_64m <- fitdist(Dioxin_60_64m, 'gamma') 
t2 <- Dioxin_60_64m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_60_64m$estimate[1],fit2_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64m, fitnames="gamma")


cvm.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Anderson-Darling test


#Age 65-69, males
Dioxin_65_69m <- subset(Dioxin_exp1, agegroups=="65-69" & sex=="1", select = totaldioxin_bw)
Dioxin_65_69m <- as.vector(Dioxin_65_69m$totaldioxin_bw)



fit1_65_69m <- fitdist(Dioxin_65_69m, "lnorm")
t <- Dioxin_65_69m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_65_69m$estimate[1], sdlog = fit1_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69m, fitnames="lnorm") 


fit2_65_69m <- fitdist(Dioxin_65_69m, 'gamma') 
t2 <- Dioxin_65_69m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_65_69m$estimate[1],fit2_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69m, fitnames="gamma")


cvm.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Anderson-Darling test



#Age 70-74, males
Dioxin_70_74m <- subset(Dioxin_exp1, agegroups=="70-74" & sex=="1", select = totaldioxin_bw)
Dioxin_70_74m <- as.vector(Dioxin_70_74m$totaldioxin_bw)



fit1_70_74m <- fitdist(Dioxin_70_74m, "lnorm")
t <- Dioxin_70_74m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_70_74m$estimate[1], sdlog = fit1_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74m, fitnames="lnorm") 


fit2_70_74m <- fitdist(Dioxin_70_74m, 'gamma') 
t2 <- Dioxin_70_74m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_70_74m$estimate[1],fit2_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74m, fitnames="gamma")


cvm.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Anderson-Darling test



#Age 75-79, males
Dioxin_75_79m <- subset(Dioxin_exp1, agegroups=="75-79" & sex=="1", select = totaldioxin_bw)
Dioxin_75_79m <- as.vector(Dioxin_75_79m$totaldioxin_bw)



fit1_75_79m <- fitdist(Dioxin_75_79m, "lnorm")
t <- Dioxin_75_79m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_75_79m$estimate[1], sdlog = fit1_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79m, fitnames="lnorm") 


fit2_75_79m <- fitdist(Dioxin_75_79m, 'gamma') 
t2 <- Dioxin_75_79m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_75_79m$estimate[1],fit2_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79m, fitnames="gamma")


cvm.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Anderson-Darling test



##### Fitting exposures female #####


#Age 15-19, female
Dioxin_15_19w <- subset(Dioxin_exp1, agegroups=="15-19" & sex=="2", select = totaldioxin_bw)
Dioxin_15_19w <- as.vector(Dioxin_15_19w$totaldioxin_bw)

fit1_15_19w <- fitdist(Dioxin_15_19w, "lnorm")
t <- Dioxin_15_19w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_15_19w$estimate[1], sdlog = fit1_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19w, fitnames="lnorm") 


fit2_15_19w <- fitdist(Dioxin_15_19w, 'gamma') 
t2 <- Dioxin_15_19w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_15_19w$estimate[1],fit2_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19w, fitnames="gamma")

cvm.test(t, plnorm, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2])  #Anderson-Darling test


#Age 20-24, female
Dioxin_20_24w <- subset(Dioxin_exp1, agegroups=="20-24" & sex=="2", select = totaldioxin_bw)
Dioxin_20_24w <- as.vector(Dioxin_20_24w$totaldioxin_bw)

fit1_20_24w <- fitdist(Dioxin_20_24w, "lnorm")
t <- Dioxin_20_24w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_20_24w$estimate[1], sdlog = fit1_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24w, fitnames="lnorm") 

fit2_20_24w <- fitdist(Dioxin_20_24w, 'gamma') 
t2 <- Dioxin_20_24w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_20_24w$estimate[1],fit2_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24w, fitnames="gamma")

cvm.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Anderson-Darling test


#Age 25-29, female
Dioxin_25_29w <- subset(Dioxin_exp1, agegroups=="25-29" & sex=="2", select = totaldioxin_bw)
Dioxin_25_29w <- as.vector(Dioxin_25_29w$totaldioxin_bw)

fit1_25_29w <- fitdist(Dioxin_25_29w, "lnorm")
t <- Dioxin_25_29w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_25_29w$estimate[1], sdlog = fit1_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29w, fitnames="lnorm") 

fit2_25_29w <- fitdist(Dioxin_25_29w, 'gamma')
t2 <- Dioxin_25_29w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_25_29w$estimate[1],fit2_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29w, fitnames="gamma")

cvm.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Anderson-Darling test


#Age 30-34, female
Dioxin_30_34w <- subset(Dioxin_exp1, agegroups=="30-34" & sex=="2", select = totaldioxin_bw)
Dioxin_30_34w <- as.vector(Dioxin_30_34w$totaldioxin_bw)

fit1_30_34w <- fitdist(Dioxin_30_34w, "lnorm")
t <- Dioxin_30_34w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_30_34w$estimate[1], sdlog = fit1_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34w, fitnames="lnorm") 

fit2_30_34w <- fitdist(Dioxin_30_34w, 'gamma') 
t2 <- Dioxin_30_34w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_30_34w$estimate[1],fit2_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34w, fitnames="gamma")

cvm.test(t, plnorm, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2])  #Anderson-Darling test


#Age 35-39, female
Dioxin_35_39w <- subset(Dioxin_exp1, agegroups=="35-39" & sex=="2", select = totaldioxin_bw)
Dioxin_35_39w <- as.vector(Dioxin_35_39w$totaldioxin_bw)

fit1_35_39w <- fitdist(Dioxin_35_39w, "lnorm")
t <- Dioxin_35_39w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_35_39w$estimate[1], sdlog = fit1_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39w, fitnames="lnorm") 

fit2_35_39w <- fitdist(Dioxin_35_39w, 'gamma') 
t2 <- Dioxin_35_39w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_35_39w$estimate[1],fit2_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39w, fitnames="gamma")

cvm.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Anderson-Darling test


#Age 40-44, female
Dioxin_40_44w <- subset(Dioxin_exp1, agegroups=="40-44" & sex=="2", select = totaldioxin_bw)
Dioxin_40_44w <- as.vector(Dioxin_40_44w$totaldioxin_bw)

fit1_40_44w <- fitdist(Dioxin_40_44w, "lnorm")
t <- Dioxin_40_44w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_40_44w$estimate[1], sdlog = fit1_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44w, fitnames="lnorm") 

fit2_40_44w <- fitdist(Dioxin_40_44w, 'gamma')
t2 <- Dioxin_40_44w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_40_44w$estimate[1],fit2_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44w, fitnames="gamma")

cvm.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Anderson-Darling test


#Age 45-49, female
Dioxin_45_49w <- subset(Dioxin_exp1, agegroups=="45-49" & sex=="2", select = totaldioxin_bw)
Dioxin_45_49w <- as.vector(Dioxin_45_49w$totaldioxin_bw)

fit1_45_49w <- fitdist(Dioxin_45_49w, "lnorm")
t <- Dioxin_45_49w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_45_49w$estimate[1], sdlog = fit1_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49w, fitnames="lnorm") 

fit2_45_49w <- fitdist(Dioxin_45_49w, 'gamma')
t2 <- Dioxin_45_49w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_45_49w$estimate[1],fit2_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49w, fitnames="gamma")

cvm.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Anderson-Darling test


#Age 50-54, female
Dioxin_50_54w <- subset(Dioxin_exp1, agegroups=="50-54" & sex=="2", select = totaldioxin_bw)
Dioxin_50_54w <- as.vector(Dioxin_50_54w$totaldioxin_bw)

fit1_50_54w <- fitdist(Dioxin_50_54w, "lnorm")
t <- Dioxin_50_54w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_50_54w$estimate[1], sdlog = fit1_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54w, fitnames="lnorm") 

fit2_50_54w <- fitdist(Dioxin_50_54w, 'gamma')
t2 <- Dioxin_50_54w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_50_54w$estimate[1],fit2_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54w, fitnames="gamma")

cvm.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Anderson-Darling test


#Age 55-59, female
Dioxin_55_59w <- subset(Dioxin_exp1, agegroups=="55-59" & sex=="2", select = totaldioxin_bw)
Dioxin_55_59w <- as.vector(Dioxin_55_59w$totaldioxin_bw)

fit1_55_59w <- fitdist(Dioxin_55_59w, "lnorm")
t <- Dioxin_55_59w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_55_59w$estimate[1], sdlog = fit1_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59w, fitnames="lnorm") 

fit2_55_59w <- fitdist(Dioxin_55_59w, 'gamma')
t2 <- Dioxin_55_59w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_55_59w$estimate[1],fit2_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59w, fitnames="gamma")

cvm.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Anderson-Darling test



#Age 60-64, female
Dioxin_60_64w <- subset(Dioxin_exp1, agegroups=="60-64" & sex=="2", select = totaldioxin_bw)
Dioxin_60_64w <- as.vector(Dioxin_60_64w$totaldioxin_bw)

fit1_60_64w <- fitdist(Dioxin_60_64w, "lnorm")
t <- Dioxin_60_64w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_60_64w$estimate[1], sdlog = fit1_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64w, fitnames="lnorm") 

fit2_60_64w <- fitdist(Dioxin_60_64w, 'gamma')
t2 <- Dioxin_60_64w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_60_64w$estimate[1],fit2_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64w, fitnames="gamma")

cvm.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Anderson-Darling test


#Age 65-69, female
Dioxin_65_69w <- subset(Dioxin_exp1, agegroups=="65-69" & sex=="2", select = totaldioxin_bw)
Dioxin_65_69w <- as.vector(Dioxin_65_69w$totaldioxin_bw)

fit1_65_69w <- fitdist(Dioxin_65_69w, "lnorm")
t <- Dioxin_65_69w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_65_69w$estimate[1], sdlog = fit1_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69w, fitnames="lnorm") 

fit2_65_69w <- fitdist(Dioxin_65_69w, 'gamma')
t2 <- Dioxin_65_69w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_65_69w$estimate[1],fit2_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69w, fitnames="gamma")

cvm.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Anderson-Darling test


#Age 70-74, female
Dioxin_70_74w <- subset(Dioxin_exp1, agegroups=="70-74" & sex=="2", select = totaldioxin_bw)
Dioxin_70_74w <- as.vector(Dioxin_70_74w$totaldioxin_bw)

fit1_70_74w <- fitdist(Dioxin_70_74w, "lnorm")
t <- Dioxin_70_74w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_70_74w$estimate[1], sdlog = fit1_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74w, fitnames="lnorm") 

fit2_70_74w <- fitdist(Dioxin_70_74w, 'gamma')
t2 <- Dioxin_70_74w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_70_74w$estimate[1],fit2_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74w, fitnames="gamma")

cvm.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Anderson-Darling test


#Age 75-79, female
Dioxin_75_79w <- subset(Dioxin_exp1, agegroups=="75-79" & sex=="2", select = totaldioxin_bw)
Dioxin_75_79w <- as.vector(Dioxin_75_79w$totaldioxin_bw)

fit1_75_79w <- fitdist(Dioxin_75_79w, "lnorm")
t <- Dioxin_75_79w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_75_79w$estimate[1], sdlog = fit1_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79w, fitnames="lnorm") 

fit2_75_79w <- fitdist(Dioxin_75_79w, 'gamma')
t2 <- Dioxin_75_79w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_75_79w$estimate[1],fit2_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79w, fitnames="gamma")

cvm.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Anderson-Darling test



##### Prob of hypothyroidism males #####

# Assume exposure is variable - Lognormal has the best fit


# Age 15-19, males
set.seed(1)
exp <- rlnorm(nvar, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2]) 

pTT4_15_19m <- apply(ICED, 2, function(x) mean(x < exp))


# Age 20-24, males
set.seed(1)
exp <- rlnorm(nvar, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2]) 

pTT4_20_24m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 25-29, males
set.seed(1)
exp <- rlnorm(nvar, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2]) 

pTT4_25_29m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 30-34, males
set.seed(1)
exp <- rlnorm(nvar, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2]) 

pTT4_30_34m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 35-39, males
set.seed(1)
exp <- rlnorm(nvar, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2]) 

pTT4_35_39m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 40-44, males
set.seed(1)
exp <- rlnorm(nvar, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2]) 

pTT4_40_44m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 45-49, males
set.seed(1)
exp <- rlnorm(nvar, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2]) 

pTT4_45_49m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 50-54, males
set.seed(1)
exp <- rlnorm(nvar, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2]) 

pTT4_50_54m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 55-59, males
set.seed(1)
exp <- rlnorm(nvar, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2]) 

pTT4_55_59m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 60-64, males
set.seed(1)
exp <- rlnorm(nvar, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2]) 

pTT4_60_64m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 65-69, males
set.seed(1)
exp <- rlnorm(nvar, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2]) 

pTT4_65_69m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 70-74, males
set.seed(1)
exp <- rlnorm(nvar, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2]) 

pTT4_70_74m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 75-79, males
set.seed(1)
exp <- rlnorm(nvar, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2]) 

pTT4_75_79m <- apply(ICED, 2, function(x) mean(x < exp))


#Assume age 80-84 and 85+ has the same prob as 75-79


##### Prob of hypothyroidism females #####


# Age 15-19, females
set.seed(1)
exp <- rlnorm(nvar, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2]) 

pTT4_15_19w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 20-24, females
set.seed(1)
exp <- rlnorm(nvar, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2]) 

pTT4_20_24w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 25-29, females
set.seed(1)
exp <- rlnorm(nvar, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2]) 

pTT4_25_29w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 30-34, females
set.seed(1)
exp <- rlnorm(nvar, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2]) 

pTT4_30_34w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 35-39, females
set.seed(1)
exp <- rlnorm(nvar, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2]) 

pTT4_35_39w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 40-44, females
set.seed(1)
exp <- rlnorm(nvar, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2]) 

pTT4_40_44w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 45-49, females
set.seed(1)
exp <- rlnorm(nvar, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2]) 

pTT4_45_49w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 50-54, females
set.seed(1)
exp <- rlnorm(nvar, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2]) 

pTT4_50_54w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 55-59, females
set.seed(1)
exp <- rlnorm(nvar, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2]) 

pTT4_55_59w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 60-64, females
set.seed(1)
exp <- rlnorm(nvar, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2]) 

pTT4_60_64w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 65-69, females
set.seed(1)
exp <- rlnorm(nvar, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2]) 

pTT4_65_69w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 70-74, females
set.seed(1)
exp <- rlnorm(nvar, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2]) 

pTT4_70_74w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 75-79, females
set.seed(1)
exp <- rlnorm(nvar, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2]) 

pTT4_75_79w <- apply(ICED, 2, function(x) mean(x < exp))


#Assume age 80-84 and 85+ has the same prob as 75-79


##### DALY calculation #####

##### Males #####

#Age 15-19, male

DALY15_19m_scen1 <- pTT4_15_19m * pop_15_19m * D * dw
mean_median_ci(DALY15_19m_scen1)


#Age 20-24y, male

DALY20_24m_scen1 <- pTT4_20_24m * pop_20_24m * D * dw
mean_median_ci(DALY20_24m_scen1)


#Age 25-29y, male

DALY25_29m_scen1 <- pTT4_25_29m * pop_25_29m * D * dw
mean_median_ci(DALY25_29m_scen1)


#Age 30-34y, male

DALY30_34m_scen1 <- pTT4_30_34m * pop_30_34m * D * dw
mean_median_ci(DALY30_34m_scen1)


#Age 35-39y, male

DALY35_39m_scen1 <- pTT4_35_39m * pop_35_39m * D * dw
mean_median_ci(DALY35_39m_scen1)



#Age 40-44y, male

DALY40_44m_scen1 <- pTT4_40_44m * pop_40_44m * D * dw
mean_median_ci(DALY40_44m_scen1)


#Age 45-49y, male

DALY45_49m_scen1 <- pTT4_45_49m * pop_45_49m * D * dw
mean_median_ci(DALY45_49m_scen1)


#Age 50-54y, male

DALY50_54m_scen1 <- pTT4_50_54m * pop_50_54m * D * dw
mean_median_ci(DALY50_54m_scen1)


#Age 55-59y, male

DALY55_59m_scen1 <- pTT4_55_59m * pop_55_59m * D * dw
mean_median_ci(DALY55_59m_scen1)


#Age 60-64y, male

DALY60_64m_scen1 <- pTT4_60_64m * pop_60_64m * D * dw
mean_median_ci(DALY60_64m_scen1)


#Age 65-69y, male

DALY65_69m_scen1 <- pTT4_65_69m * pop_65_69m * D * dw
mean_median_ci(DALY65_69m_scen1)


#Age 70-74y, male

DALY70_74m_scen1 <- pTT4_70_74m * pop_70_74m * D * dw
mean_median_ci(DALY70_74m_scen1)


#Age 75-79y, male

DALY75_79m_scen1 <- pTT4_75_79m * pop_75_79m * D * dw
mean_median_ci(DALY75_79m_scen1)


#Age 80-84y, male

DALY80_84m_scen1 <- pTT4_75_79m * pop_80_84m * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY80_84m_scen1)


#Age 85+y, male

DALY85m_scen1 <- pTT4_75_79m * pop_85m * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY85m_scen1)


##### Females #####

#Age 15-19, female

DALY15_19w_scen1 <- pTT4_15_19w * pop_15_19w * D * dw
mean_median_ci(DALY15_19w_scen1)


#Age 20-24y, female

DALY20_24w_scen1 <- pTT4_20_24w * pop_20_24w * D * dw
mean_median_ci(DALY20_24w_scen1)


#Age 25-29y, female

DALY25_29w_scen1 <- pTT4_25_29w * pop_25_29w * D * dw
mean_median_ci(DALY25_29w_scen1)


#Age 30-34y, female

DALY30_34w_scen1 <- pTT4_30_34w * pop_30_34w * D * dw
mean_median_ci(DALY30_34w_scen1)


#Age 35-39y, female

DALY35_39w_scen1 <- pTT4_35_39w * pop_35_39w * D * dw
mean_median_ci(DALY35_39w_scen1)



#Age 40-44y, female

DALY40_44w_scen1 <- pTT4_40_44w * pop_40_44w * D * dw
mean_median_ci(DALY40_44w_scen1)


#Age 45-49y, female

DALY45_49w_scen1 <- pTT4_45_49w * pop_45_49w * D * dw
mean_median_ci(DALY45_49w_scen1)


#Age 50-54y, female

DALY50_54w_scen1 <- pTT4_50_54w * pop_50_54w * D * dw
mean_median_ci(DALY50_54w_scen1)


#Age 55-59y, female

DALY55_59w_scen1 <- pTT4_55_59w * pop_55_59w * D * dw
mean_median_ci(DALY55_59w_scen1)


#Age 60-64y, female

DALY60_64w_scen1 <- pTT4_60_64w * pop_60_64w * D * dw
mean_median_ci(DALY60_64w_scen1)


#Age 65-69y, female

DALY65_69w_scen1 <- pTT4_65_69w * pop_65_69w * D * dw
mean_median_ci(DALY65_69w_scen1)


#Age 70-74y, female

DALY70_74w_scen1 <- pTT4_70_74w * pop_70_74w * D * dw
mean_median_ci(DALY70_74w_scen1)


#Age 75-79y, female

DALY75_79w_scen1 <- pTT4_75_79w * pop_75_79w * D * dw
mean_median_ci(DALY75_79w_scen1)


#Age 80-84y, female

DALY80_84w_scen1 <- pTT4_75_79w * pop_80_84w * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY80_84w_scen1)


#Age 85+y, female

DALY85w_scen1 <- pTT4_75_79w * pop_85w * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY85w_scen1)


##### Total DALYs scen 1 #####

DALYtotal_scen1 <- DALY15_19m_scen1 + DALY20_24m_scen1 + DALY25_29m_scen1 + DALY30_34m_scen1 + DALY35_39m_scen1 +
  DALY40_44m_scen1 + DALY45_49m_scen1 + DALY50_54m_scen1 +  DALY55_59m_scen1 + DALY60_64m_scen1 + DALY65_69m_scen1 +
  DALY70_74m_scen1 + DALY75_79m_scen1 + DALY80_84m_scen1 + DALY85m_scen1 +
  DALY15_19w_scen1 + DALY20_24w_scen1 + DALY25_29w_scen1 + DALY30_34w_scen1 + DALY35_39w_scen1 + DALY40_44w_scen1 +
  DALY45_49w_scen1 + DALY50_54w_scen1 + DALY55_59w_scen1 + DALY60_64w_scen1 + DALY65_69w_scen1 + DALY70_74w_scen1 +
  DALY75_79w_scen1 + DALY80_84w_scen1 + DALY85w_scen1


mean_median_ci(DALYtotal_scen1)


##### Cases scenario 1 #####

cases_scen1 <- pTT4_15_19m * pop_15_19m + pTT4_20_24m * pop_20_24m + pTT4_25_29m * pop_25_29m + pTT4_30_34m * pop_30_34m +
  pTT4_35_39m * pop_35_39m + pTT4_40_44m * pop_40_44m + pTT4_45_49m * pop_45_49m + pTT4_50_54m * pop_50_54m + pTT4_55_59m * pop_55_59m +
  pTT4_60_64m * pop_60_64m + pTT4_65_69m * pop_65_69m + pTT4_70_74m * pop_70_74m + pTT4_75_79m * pop_75_79m + pTT4_75_79m * pop_80_84m +
  pTT4_75_79m * pop_85m +
  pTT4_15_19w * pop_15_19w + pTT4_20_24w * pop_20_24w + pTT4_25_29w * pop_25_29w + pTT4_30_34w * pop_30_34w +
  pTT4_35_39w * pop_35_39w + pTT4_40_44w * pop_40_44w + pTT4_45_49w * pop_45_49w + pTT4_50_54w * pop_50_54w + pTT4_55_59w * pop_55_59w +
  pTT4_60_64w * pop_60_64w + pTT4_65_69w * pop_65_69w + pTT4_70_74w * pop_70_74w + pTT4_75_79w * pop_75_79w + pTT4_75_79w * pop_80_84w +
  pTT4_75_79w * pop_85w


mean_median_ci(cases_scen1)


##### Scenario 2 #####

Dioxin_exp2 <- read.csv("Dioxin_exp_scen2.csv")

setDT(Dioxin_exp2)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Fitting exposures male#####

#Age 15-19, males
Dioxin_15_19m <- subset(Dioxin_exp2, agegroups=="15-19" & sex=="1", select = totaldioxin_bw)
Dioxin_15_19m <- as.vector(Dioxin_15_19m$totaldioxin_bw)


fit1_15_19m <- fitdist(Dioxin_15_19m, "lnorm")
t <- Dioxin_15_19m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_15_19m$estimate[1], sdlog = fit1_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19m, fitnames="lnorm") 


fit2_15_19m <- fitdist(Dioxin_15_19m, 'gamma') 
t2 <- Dioxin_15_19m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_15_19m$estimate[1],fit2_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19m, fitnames="gamma")


cvm.test(t, plnorm, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2])  #Anderson-Darling test


#Age 20-24, males
Dioxin_20_24m <- subset(Dioxin_exp2, agegroups=="20-24" & sex=="1", select = totaldioxin_bw)
Dioxin_20_24m <- as.vector(Dioxin_20_24m$totaldioxin_bw)


fit1_20_24m <- fitdist(Dioxin_20_24m, "lnorm")
t <- Dioxin_20_24m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_20_24m$estimate[1], sdlog = fit1_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24m, fitnames="lnorm") 


fit2_20_24m <- fitdist(Dioxin_20_24m, 'gamma') 
t2 <- Dioxin_20_24m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_20_24m$estimate[1],fit2_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24m, fitnames="gamma")


cvm.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Anderson-Darling test


#Age 25-29, males
Dioxin_25_29m <- subset(Dioxin_exp2, agegroups=="25-29" & sex=="1", select = totaldioxin_bw)
Dioxin_25_29m <- as.vector(Dioxin_25_29m$totaldioxin_bw)


fit1_25_29m <- fitdist(Dioxin_25_29m, "lnorm")
t <- Dioxin_25_29m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_25_29m$estimate[1], sdlog = fit1_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29m, fitnames="lnorm") 


fit2_25_29m <- fitdist(Dioxin_25_29m, 'gamma') 
t2 <- Dioxin_25_29m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_25_29m$estimate[1],fit2_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29m, fitnames="gamma")


cvm.test(t, plnorm, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2])  #Anderson-Darling test


#Age 30-34, males
Dioxin_30_34m <- subset(Dioxin_exp2, agegroups=="30-34" & sex=="1", select = totaldioxin_bw)
Dioxin_30_34m <- as.vector(Dioxin_30_34m$totaldioxin_bw)


fit1_30_34m <- fitdist(Dioxin_30_34m, "lnorm")
t <- Dioxin_30_34m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_30_34m$estimate[1], sdlog = fit1_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34m, fitnames="lnorm") 


fit2_30_34m <- fitdist(Dioxin_30_34m, 'gamma') 
t2 <- Dioxin_30_34m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_30_34m$estimate[1],fit2_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34m, fitnames="gamma")


cvm.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Anderson-Darling test


#Age 35-39, males
Dioxin_35_39m <- subset(Dioxin_exp2, agegroups=="35-39" & sex=="1", select = totaldioxin_bw)
Dioxin_35_39m <- as.vector(Dioxin_35_39m$totaldioxin_bw)


fit1_35_39m <- fitdist(Dioxin_35_39m, "lnorm")
t <- Dioxin_35_39m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_35_39m$estimate[1], sdlog = fit1_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39m, fitnames="lnorm") 


fit2_35_39m <- fitdist(Dioxin_35_39m, 'gamma') 
t2 <- Dioxin_35_39m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_35_39m$estimate[1],fit2_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39m, fitnames="gamma")


cvm.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Anderson-Darling test


#Age 40-44, males
Dioxin_40_44m <- subset(Dioxin_exp2, agegroups=="40-44" & sex=="1", select = totaldioxin_bw)
Dioxin_40_44m <- as.vector(Dioxin_40_44m$totaldioxin_bw)


fit1_40_44m <- fitdist(Dioxin_40_44m, "lnorm")
t <- Dioxin_40_44m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_40_44m$estimate[1], sdlog = fit1_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44m, fitnames="lnorm") 


fit2_40_44m <- fitdist(Dioxin_40_44m, 'gamma') 
t2 <- Dioxin_40_44m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_40_44m$estimate[1],fit2_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44m, fitnames="gamma")


cvm.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Anderson-Darling test


#Age 45-49, males
Dioxin_45_49m <- subset(Dioxin_exp2, agegroups=="45-49" & sex=="1", select = totaldioxin_bw)
Dioxin_45_49m <- as.vector(Dioxin_45_49m$totaldioxin_bw)


fit1_45_49m <- fitdist(Dioxin_45_49m, "lnorm")
t <- Dioxin_45_49m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_45_49m$estimate[1], sdlog = fit1_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49m, fitnames="lnorm") 


fit2_45_49m <- fitdist(Dioxin_45_49m, 'gamma') 
t2 <- Dioxin_45_49m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_45_49m$estimate[1],fit2_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49m, fitnames="gamma")


cvm.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Anderson-Darling test


#Age 50-54, males
Dioxin_50_54m <- subset(Dioxin_exp2, agegroups=="50-54" & sex=="1", select = totaldioxin_bw)
Dioxin_50_54m <- as.vector(Dioxin_50_54m$totaldioxin_bw)


fit1_50_54m <- fitdist(Dioxin_50_54m, "lnorm")
t <- Dioxin_50_54m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_50_54m$estimate[1], sdlog = fit1_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54m, fitnames="lnorm") 


fit2_50_54m <- fitdist(Dioxin_50_54m, 'gamma') 
t2 <- Dioxin_50_54m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_50_54m$estimate[1],fit2_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54m, fitnames="gamma")


cvm.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Anderson-Darling test


#Age 55-59, males
Dioxin_55_59m <- subset(Dioxin_exp2, agegroups=="55-59" & sex=="1", select = totaldioxin_bw)
Dioxin_55_59m <- as.vector(Dioxin_55_59m$totaldioxin_bw)


fit1_55_59m <- fitdist(Dioxin_55_59m, "lnorm")
t <- Dioxin_55_59m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_55_59m$estimate[1], sdlog = fit1_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59m, fitnames="lnorm") 


fit2_55_59m <- fitdist(Dioxin_55_59m, 'gamma') 
t2 <- Dioxin_55_59m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_55_59m$estimate[1],fit2_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59m, fitnames="gamma")


cvm.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Anderson-Darling test



#Age 60-64, males
Dioxin_60_64m <- subset(Dioxin_exp2, agegroups=="60-64" & sex=="1", select = totaldioxin_bw)
Dioxin_60_64m <- as.vector(Dioxin_60_64m$totaldioxin_bw)


fit1_60_64m <- fitdist(Dioxin_60_64m, "lnorm")
t <- Dioxin_60_64m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_60_64m$estimate[1], sdlog = fit1_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64m, fitnames="lnorm") 


fit2_60_64m <- fitdist(Dioxin_60_64m, 'gamma') 
t2 <- Dioxin_60_64m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_60_64m$estimate[1],fit2_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64m, fitnames="gamma")


cvm.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Anderson-Darling test


#Age 65-69, males
Dioxin_65_69m <- subset(Dioxin_exp2, agegroups=="65-69" & sex=="1", select = totaldioxin_bw)
Dioxin_65_69m <- as.vector(Dioxin_65_69m$totaldioxin_bw)


fit1_65_69m <- fitdist(Dioxin_65_69m, "lnorm")
t <- Dioxin_65_69m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_65_69m$estimate[1], sdlog = fit1_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69m, fitnames="lnorm") 


fit2_65_69m <- fitdist(Dioxin_65_69m, 'gamma') 
t2 <- Dioxin_65_69m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_65_69m$estimate[1],fit2_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69m, fitnames="gamma")


cvm.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Anderson-Darling test



#Age 70-74, males
Dioxin_70_74m <- subset(Dioxin_exp2, agegroups=="70-74" & sex=="1", select = totaldioxin_bw)
Dioxin_70_74m <- as.vector(Dioxin_70_74m$totaldioxin_bw)


fit1_70_74m <- fitdist(Dioxin_70_74m, "lnorm")
t <- Dioxin_70_74m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_70_74m$estimate[1], sdlog = fit1_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74m, fitnames="lnorm") 


fit2_70_74m <- fitdist(Dioxin_70_74m, 'gamma') 
t2 <- Dioxin_70_74m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_70_74m$estimate[1],fit2_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74m, fitnames="gamma")


cvm.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Anderson-Darling test



#Age 75-79, males
Dioxin_75_79m <- subset(Dioxin_exp2, agegroups=="75-79" & sex=="1", select = totaldioxin_bw)
Dioxin_75_79m <- as.vector(Dioxin_75_79m$totaldioxin_bw)


fit1_75_79m <- fitdist(Dioxin_75_79m, "lnorm")
t <- Dioxin_75_79m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_75_79m$estimate[1], sdlog = fit1_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79m, fitnames="lnorm") 


fit2_75_79m <- fitdist(Dioxin_75_79m, 'gamma') 
t2 <- Dioxin_75_79m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_75_79m$estimate[1],fit2_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79m, fitnames="gamma")


cvm.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Anderson-Darling test



##### Fitting exposures female #####


#Age 15-19, female
Dioxin_15_19w <- subset(Dioxin_exp2, agegroups=="15-19" & sex=="2", select = totaldioxin_bw)
Dioxin_15_19w <- as.vector(Dioxin_15_19w$totaldioxin_bw)

fit1_15_19w <- fitdist(Dioxin_15_19w, "lnorm")
t <- Dioxin_15_19w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_15_19w$estimate[1], sdlog = fit1_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19w, fitnames="lnorm") 


fit2_15_19w <- fitdist(Dioxin_15_19w, 'gamma') 
t2 <- Dioxin_15_19w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_15_19w$estimate[1],fit2_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19w, fitnames="gamma")

cvm.test(t, plnorm, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2])  #Anderson-Darling test


#Age 20-24, female
Dioxin_20_24w <- subset(Dioxin_exp2, agegroups=="20-24" & sex=="2", select = totaldioxin_bw)
Dioxin_20_24w <- as.vector(Dioxin_20_24w$totaldioxin_bw)

fit1_20_24w <- fitdist(Dioxin_20_24w, "lnorm")
t <- Dioxin_20_24w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_20_24w$estimate[1], sdlog = fit1_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24w, fitnames="lnorm") 

fit2_20_24w <- fitdist(Dioxin_20_24w, 'gamma') 
t2 <- Dioxin_20_24w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_20_24w$estimate[1],fit2_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24w, fitnames="gamma")

cvm.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Anderson-Darling test


#Age 25-29, female
Dioxin_25_29w <- subset(Dioxin_exp2, agegroups=="25-29" & sex=="2", select = totaldioxin_bw)
Dioxin_25_29w <- as.vector(Dioxin_25_29w$totaldioxin_bw)

fit1_25_29w <- fitdist(Dioxin_25_29w, "lnorm")
t <- Dioxin_25_29w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_25_29w$estimate[1], sdlog = fit1_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29w, fitnames="lnorm") 

fit2_25_29w <- fitdist(Dioxin_25_29w, 'gamma')
t2 <- Dioxin_25_29w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_25_29w$estimate[1],fit2_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29w, fitnames="gamma")

cvm.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Anderson-Darling test


#Age 30-34, female
Dioxin_30_34w <- subset(Dioxin_exp2, agegroups=="30-34" & sex=="2", select = totaldioxin_bw)
Dioxin_30_34w <- as.vector(Dioxin_30_34w$totaldioxin_bw)

fit1_30_34w <- fitdist(Dioxin_30_34w, "lnorm")
t <- Dioxin_30_34w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_30_34w$estimate[1], sdlog = fit1_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34w, fitnames="lnorm") 

fit2_30_34w <- fitdist(Dioxin_30_34w, 'gamma') 
t2 <- Dioxin_30_34w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_30_34w$estimate[1],fit2_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34w, fitnames="gamma")

cvm.test(t, plnorm, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2])  #Anderson-Darling test


#Age 35-39, female
Dioxin_35_39w <- subset(Dioxin_exp2, agegroups=="35-39" & sex=="2", select = totaldioxin_bw)
Dioxin_35_39w <- as.vector(Dioxin_35_39w$totaldioxin_bw)

fit1_35_39w <- fitdist(Dioxin_35_39w, "lnorm")
t <- Dioxin_35_39w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_35_39w$estimate[1], sdlog = fit1_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39w, fitnames="lnorm") 

fit2_35_39w <- fitdist(Dioxin_35_39w, 'gamma') 
t2 <- Dioxin_35_39w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_35_39w$estimate[1],fit2_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39w, fitnames="gamma")

cvm.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Anderson-Darling test


#Age 40-44, female
Dioxin_40_44w <- subset(Dioxin_exp2, agegroups=="40-44" & sex=="2", select = totaldioxin_bw)
Dioxin_40_44w <- as.vector(Dioxin_40_44w$totaldioxin_bw)

fit1_40_44w <- fitdist(Dioxin_40_44w, "lnorm")
t <- Dioxin_40_44w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_40_44w$estimate[1], sdlog = fit1_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44w, fitnames="lnorm") 

fit2_40_44w <- fitdist(Dioxin_40_44w, 'gamma')
t2 <- Dioxin_40_44w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_40_44w$estimate[1],fit2_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44w, fitnames="gamma")

cvm.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Anderson-Darling test


#Age 45-49, female
Dioxin_45_49w <- subset(Dioxin_exp2, agegroups=="45-49" & sex=="2", select = totaldioxin_bw)
Dioxin_45_49w <- as.vector(Dioxin_45_49w$totaldioxin_bw)

fit1_45_49w <- fitdist(Dioxin_45_49w, "lnorm")
t <- Dioxin_45_49w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_45_49w$estimate[1], sdlog = fit1_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49w, fitnames="lnorm") 

fit2_45_49w <- fitdist(Dioxin_45_49w, 'gamma')
t2 <- Dioxin_45_49w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_45_49w$estimate[1],fit2_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49w, fitnames="gamma")

cvm.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Anderson-Darling test


#Age 50-54, female
Dioxin_50_54w <- subset(Dioxin_exp2, agegroups=="50-54" & sex=="2", select = totaldioxin_bw)
Dioxin_50_54w <- as.vector(Dioxin_50_54w$totaldioxin_bw)

fit1_50_54w <- fitdist(Dioxin_50_54w, "lnorm")
t <- Dioxin_50_54w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_50_54w$estimate[1], sdlog = fit1_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54w, fitnames="lnorm") 

fit2_50_54w <- fitdist(Dioxin_50_54w, 'gamma')
t2 <- Dioxin_50_54w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_50_54w$estimate[1],fit2_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54w, fitnames="gamma")

cvm.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Anderson-Darling test


#Age 55-59, female
Dioxin_55_59w <- subset(Dioxin_exp2, agegroups=="55-59" & sex=="2", select = totaldioxin_bw)
Dioxin_55_59w <- as.vector(Dioxin_55_59w$totaldioxin_bw)

fit1_55_59w <- fitdist(Dioxin_55_59w, "lnorm")
t <- Dioxin_55_59w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_55_59w$estimate[1], sdlog = fit1_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59w, fitnames="lnorm") 

fit2_55_59w <- fitdist(Dioxin_55_59w, 'gamma')
t2 <- Dioxin_55_59w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_55_59w$estimate[1],fit2_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59w, fitnames="gamma")

cvm.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Anderson-Darling test



#Age 60-64, female
Dioxin_60_64w <- subset(Dioxin_exp2, agegroups=="60-64" & sex=="2", select = totaldioxin_bw)
Dioxin_60_64w <- as.vector(Dioxin_60_64w$totaldioxin_bw)

fit1_60_64w <- fitdist(Dioxin_60_64w, "lnorm")
t <- Dioxin_60_64w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_60_64w$estimate[1], sdlog = fit1_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64w, fitnames="lnorm") 

fit2_60_64w <- fitdist(Dioxin_60_64w, 'gamma')
t2 <- Dioxin_60_64w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_60_64w$estimate[1],fit2_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64w, fitnames="gamma")

cvm.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Anderson-Darling test


#Age 65-69, female
Dioxin_65_69w <- subset(Dioxin_exp2, agegroups=="65-69" & sex=="2", select = totaldioxin_bw)
Dioxin_65_69w <- as.vector(Dioxin_65_69w$totaldioxin_bw)

fit1_65_69w <- fitdist(Dioxin_65_69w, "lnorm")
t <- Dioxin_65_69w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_65_69w$estimate[1], sdlog = fit1_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69w, fitnames="lnorm") 

fit2_65_69w <- fitdist(Dioxin_65_69w, 'gamma')
t2 <- Dioxin_65_69w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_65_69w$estimate[1],fit2_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69w, fitnames="gamma")

cvm.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Anderson-Darling test


#Age 70-74, female
Dioxin_70_74w <- subset(Dioxin_exp2, agegroups=="70-74" & sex=="2", select = totaldioxin_bw)
Dioxin_70_74w <- as.vector(Dioxin_70_74w$totaldioxin_bw)

fit1_70_74w <- fitdist(Dioxin_70_74w, "lnorm")
t <- Dioxin_70_74w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_70_74w$estimate[1], sdlog = fit1_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74w, fitnames="lnorm") 

fit2_70_74w <- fitdist(Dioxin_70_74w, 'gamma')
t2 <- Dioxin_70_74w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_70_74w$estimate[1],fit2_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74w, fitnames="gamma")

cvm.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Anderson-Darling test


#Age 75-79, female
Dioxin_75_79w <- subset(Dioxin_exp2, agegroups=="75-79" & sex=="2", select = totaldioxin_bw)
Dioxin_75_79w <- as.vector(Dioxin_75_79w$totaldioxin_bw)

fit1_75_79w <- fitdist(Dioxin_75_79w, "lnorm")
t <- Dioxin_75_79w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_75_79w$estimate[1], sdlog = fit1_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79w, fitnames="lnorm") 

fit2_75_79w <- fitdist(Dioxin_75_79w, 'gamma')
t2 <- Dioxin_75_79w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_75_79w$estimate[1],fit2_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79w, fitnames="gamma")

cvm.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Anderson-Darling test


##### Prob of hypothyroidism males #####

# Assume exposure is variable - Lognormal has the best fit


# Age 15-19, males
set.seed(1)
exp <- rlnorm(nvar, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2]) 

pTT4_15_19m <- apply(ICED, 2, function(x) mean(x < exp))


# Age 20-24, males
set.seed(1)
exp <- rlnorm(nvar, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2]) 

pTT4_20_24m <- apply(ICED, 2, function(x) mean(x < exp))


# Age 25-29, males
set.seed(1)
exp <- rlnorm(nvar, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2]) 

pTT4_25_29m <- apply(ICED, 2, function(x) mean(x < exp))


# Age 30-34, males
set.seed(1)
exp <- rlnorm(nvar, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2]) 

pTT4_30_34m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 35-39, males
set.seed(1)
exp <- rlnorm(nvar, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2]) 

pTT4_35_39m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 40-44, males
set.seed(1)
exp <- rlnorm(nvar, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2]) 

pTT4_40_44m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 45-49, males
set.seed(1)
exp <- rlnorm(nvar, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2]) 

pTT4_45_49m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 50-54, males
set.seed(1)
exp <- rlnorm(nvar, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2]) 

pTT4_50_54m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 55-59, males
set.seed(1)
exp <- rlnorm(nvar, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2]) 

pTT4_55_59m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 60-64, males
set.seed(1)
exp <- rlnorm(nvar, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2]) 

pTT4_60_64m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 65-69, males
set.seed(1)
exp <- rlnorm(nvar, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2]) 

pTT4_65_69m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 70-74, males
set.seed(1)
exp <- rlnorm(nvar, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2]) 

pTT4_70_74m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 75-79, males
set.seed(1)
exp <- rlnorm(nvar, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2]) 

pTT4_75_79m <- apply(ICED, 2, function(x) mean(x < exp))


#Assume age 80-84 and 85+ has the same prob as 75-79


##### Prob of hypothyroidism females #####


# Assume exposure is variable

# Age 15-19, females
set.seed(1)
exp <- rlnorm(nvar, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2]) 

pTT4_15_19w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 20-24, females
set.seed(1)
exp <- rlnorm(nvar, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2]) 

pTT4_20_24w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 25-29, females
set.seed(1)
exp <- rlnorm(nvar, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2]) 

pTT4_25_29w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 30-34, females
set.seed(1)
exp <- rlnorm(nvar, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2]) 

pTT4_30_34w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 35-39, females
set.seed(1)
exp <- rlnorm(nvar, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2]) 

pTT4_35_39w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 40-44, females
set.seed(1)
exp <- rlnorm(nvar, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2]) 

pTT4_40_44w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 45-49, females
set.seed(1)
exp <- rlnorm(nvar, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2]) 

pTT4_45_49w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 50-54, females
set.seed(1)
exp <- rlnorm(nvar, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2]) 

pTT4_50_54w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 55-59, females
set.seed(1)
exp <- rlnorm(nvar, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2]) 

pTT4_55_59w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 60-64, females
set.seed(1)
exp <- rlnorm(nvar, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2]) 

pTT4_60_64w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 65-69, females
set.seed(1)
exp <- rlnorm(nvar, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2]) 

pTT4_65_69w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 70-74, females
set.seed(1)
exp <- rlnorm(nvar, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2]) 

pTT4_70_74w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 75-79, females
set.seed(1)
exp <- rlnorm(nvar, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2]) 

pTT4_75_79w <- apply(ICED, 2, function(x) mean(x < exp))


#Assume age 80-84 and 85+ has the same prob as 75-79


##### DALY calculation #####

##### Males #####

#Age 15-19, male

DALY15_19m_scen2 <- pTT4_15_19m * pop_15_19m * D * dw
mean_median_ci(DALY15_19m_scen2)


#Age 20-24y, male

DALY20_24m_scen2 <- pTT4_20_24m * pop_20_24m * D * dw
mean_median_ci(DALY20_24m_scen2)


#Age 25-29y, male

DALY25_29m_scen2 <- pTT4_25_29m * pop_25_29m * D * dw
mean_median_ci(DALY25_29m_scen2)


#Age 30-34y, male

DALY30_34m_scen2 <- pTT4_30_34m * pop_30_34m * D * dw
mean_median_ci(DALY30_34m_scen2)


#Age 35-39y, male

DALY35_39m_scen2 <- pTT4_35_39m * pop_35_39m * D * dw
mean_median_ci(DALY35_39m_scen2)



#Age 40-44y, male

DALY40_44m_scen2 <- pTT4_40_44m * pop_40_44m * D * dw
mean_median_ci(DALY40_44m_scen2)


#Age 45-49y, male

DALY45_49m_scen2 <- pTT4_45_49m * pop_45_49m * D * dw
mean_median_ci(DALY45_49m_scen2)


#Age 50-54y, male

DALY50_54m_scen2 <- pTT4_50_54m * pop_50_54m * D * dw
mean_median_ci(DALY50_54m_scen2)


#Age 55-59y, male

DALY55_59m_scen2 <- pTT4_55_59m * pop_55_59m * D * dw
mean_median_ci(DALY55_59m_scen2)


#Age 60-64y, male

DALY60_64m_scen2 <- pTT4_60_64m * pop_60_64m * D * dw
mean_median_ci(DALY60_64m_scen2)


#Age 65-69y, male

DALY65_69m_scen2 <- pTT4_65_69m * pop_65_69m * D * dw
mean_median_ci(DALY65_69m_scen2)


#Age 70-74y, male

DALY70_74m_scen2 <- pTT4_70_74m * pop_70_74m * D * dw
mean_median_ci(DALY70_74m_scen2)


#Age 75-79y, male

DALY75_79m_scen2 <- pTT4_75_79m * pop_75_79m * D * dw
mean_median_ci(DALY75_79m_scen2)


#Age 80-84y, male

DALY80_84m_scen2 <- pTT4_75_79m * pop_80_84m * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY80_84m_scen2)


#Age 85+y, male

DALY85m_scen2 <- pTT4_75_79m * pop_85m * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY85m_scen2)


##### Females #####

#Age 15-19, female

DALY15_19w_scen2 <- pTT4_15_19w * pop_15_19w * D * dw
mean_median_ci(DALY15_19w_scen2)


#Age 20-24y, female

DALY20_24w_scen2 <- pTT4_20_24w * pop_20_24w * D * dw
mean_median_ci(DALY20_24w_scen2)


#Age 25-29y, female

DALY25_29w_scen2 <- pTT4_25_29w * pop_25_29w * D * dw
mean_median_ci(DALY25_29w_scen2)


#Age 30-34y, female

DALY30_34w_scen2 <- pTT4_30_34w * pop_30_34w * D * dw
mean_median_ci(DALY30_34w_scen2)


#Age 35-39y, female

DALY35_39w_scen2 <- pTT4_35_39w * pop_35_39w * D * dw
mean_median_ci(DALY35_39w_scen2)



#Age 40-44y, female

DALY40_44w_scen2 <- pTT4_40_44w * pop_40_44w * D * dw
mean_median_ci(DALY40_44w_scen2)


#Age 45-49y, female

DALY45_49w_scen2 <- pTT4_45_49w * pop_45_49w * D * dw
mean_median_ci(DALY45_49w_scen2)


#Age 50-54y, female

DALY50_54w_scen2 <- pTT4_50_54w * pop_50_54w * D * dw
mean_median_ci(DALY50_54w_scen2)


#Age 55-59y, female

DALY55_59w_scen2 <- pTT4_55_59w * pop_55_59w * D * dw
mean_median_ci(DALY55_59w_scen2)


#Age 60-64y, female

DALY60_64w_scen2 <- pTT4_60_64w * pop_60_64w * D * dw
mean_median_ci(DALY60_64w_scen2)


#Age 65-69y, female

DALY65_69w_scen2 <- pTT4_65_69w * pop_65_69w * D * dw
mean_median_ci(DALY65_69w_scen2)


#Age 70-74y, female

DALY70_74w_scen2 <- pTT4_70_74w * pop_70_74w * D * dw
mean_median_ci(DALY70_74w_scen2)


#Age 75-79y, female

DALY75_79w_scen2 <- pTT4_75_79w * pop_75_79w * D * dw
mean_median_ci(DALY75_79w_scen2)


#Age 80-84y, female

DALY80_84w_scen2 <- pTT4_75_79w * pop_80_84w * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY80_84w_scen2)


#Age 85+y, female

DALY85w_scen2 <- pTT4_75_79w * pop_85w * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY85w_scen2)


##### Total DALYs scen 2 #####

DALYtotal_scen2 <- DALY15_19m_scen2 + DALY20_24m_scen2 + DALY25_29m_scen2 + DALY30_34m_scen2 + DALY35_39m_scen2 +
  DALY40_44m_scen2 + DALY45_49m_scen2 + DALY50_54m_scen2 +  DALY55_59m_scen2 + DALY60_64m_scen2 + DALY65_69m_scen2 +
  DALY70_74m_scen2 + DALY75_79m_scen2 + DALY80_84m_scen2 + DALY85m_scen2 +
  DALY15_19w_scen2 + DALY20_24w_scen2 + DALY25_29w_scen2 + DALY30_34w_scen2 + DALY35_39w_scen2 + DALY40_44w_scen2 +
  DALY45_49w_scen2 + DALY50_54w_scen2 + DALY55_59w_scen2 + DALY60_64w_scen2 + DALY65_69w_scen2 + DALY70_74w_scen2 +
  DALY75_79w_scen2 + DALY80_84w_scen2 + DALY85w_scen2


mean_median_ci(DALYtotal_scen2)


##### Cases scenario 2 #####

cases_scen2 <- pTT4_15_19m * pop_15_19m + pTT4_20_24m * pop_20_24m + pTT4_25_29m * pop_25_29m + pTT4_30_34m * pop_30_34m +
  pTT4_35_39m * pop_35_39m + pTT4_40_44m * pop_40_44m + pTT4_45_49m * pop_45_49m + pTT4_50_54m * pop_50_54m + pTT4_55_59m * pop_55_59m +
  pTT4_60_64m * pop_60_64m + pTT4_65_69m * pop_65_69m + pTT4_70_74m * pop_70_74m + pTT4_75_79m * pop_75_79m + pTT4_75_79m * pop_80_84m +
  pTT4_75_79m * pop_85m +
  pTT4_15_19w * pop_15_19w + pTT4_20_24w * pop_20_24w + pTT4_25_29w * pop_25_29w + pTT4_30_34w * pop_30_34w +
  pTT4_35_39w * pop_35_39w + pTT4_40_44w * pop_40_44w + pTT4_45_49w * pop_45_49w + pTT4_50_54w * pop_50_54w + pTT4_55_59w * pop_55_59w +
  pTT4_60_64w * pop_60_64w + pTT4_65_69w * pop_65_69w + pTT4_70_74w * pop_70_74w + pTT4_75_79w * pop_75_79w + pTT4_75_79w * pop_80_84w +
  pTT4_75_79w * pop_85w


mean_median_ci(cases_scen2)


##### Scenario 3 #####

Dioxin_exp3 <- read.csv("Dioxin_exp_scen3.csv")

setDT(Dioxin_exp3)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Fitting exposures male #####

#Age 15-19, males
Dioxin_15_19m <- subset(Dioxin_exp3, agegroups=="15-19" & sex=="1", select = totaldioxin_bw)
Dioxin_15_19m <- as.vector(Dioxin_15_19m$totaldioxin_bw)


fit1_15_19m <- fitdist(Dioxin_15_19m, "lnorm")
t <- Dioxin_15_19m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_15_19m$estimate[1], sdlog = fit1_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19m, fitnames="lnorm") 


fit2_15_19m <- fitdist(Dioxin_15_19m, 'gamma') 
t2 <- Dioxin_15_19m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_15_19m$estimate[1],fit2_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19m, fitnames="gamma")


cvm.test(t, plnorm, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2])  #Anderson-Darling test


#Age 20-24, males
Dioxin_20_24m <- subset(Dioxin_exp3, agegroups=="20-24" & sex=="1", select = totaldioxin_bw)
Dioxin_20_24m <- as.vector(Dioxin_20_24m$totaldioxin_bw)


fit1_20_24m <- fitdist(Dioxin_20_24m, "lnorm")
t <- Dioxin_20_24m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_20_24m$estimate[1], sdlog = fit1_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24m, fitnames="lnorm") 


fit2_20_24m <- fitdist(Dioxin_20_24m, 'gamma') 
t2 <- Dioxin_20_24m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_20_24m$estimate[1],fit2_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24m, fitnames="gamma")


cvm.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Anderson-Darling test


#Age 25-29, males
Dioxin_25_29m <- subset(Dioxin_exp3, agegroups=="25-29" & sex=="1", select = totaldioxin_bw)
Dioxin_25_29m <- as.vector(Dioxin_25_29m$totaldioxin_bw)


fit1_25_29m <- fitdist(Dioxin_25_29m, "lnorm")
t <- Dioxin_25_29m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_25_29m$estimate[1], sdlog = fit1_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29m, fitnames="lnorm") 


fit2_25_29m <- fitdist(Dioxin_25_29m, 'gamma') 
t2 <- Dioxin_25_29m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_25_29m$estimate[1],fit2_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29m, fitnames="gamma")


cvm.test(t, plnorm, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2])  #Anderson-Darling test


#Age 30-34, males
Dioxin_30_34m <- subset(Dioxin_exp3, agegroups=="30-34" & sex=="1", select = totaldioxin_bw)
Dioxin_30_34m <- as.vector(Dioxin_30_34m$totaldioxin_bw)


fit1_30_34m <- fitdist(Dioxin_30_34m, "lnorm")
t <- Dioxin_30_34m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_30_34m$estimate[1], sdlog = fit1_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34m, fitnames="lnorm") 


fit2_30_34m <- fitdist(Dioxin_30_34m, 'gamma') 
t2 <- Dioxin_30_34m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_30_34m$estimate[1],fit2_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34m, fitnames="gamma")


cvm.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Anderson-Darling test


#Age 35-39, males
Dioxin_35_39m <- subset(Dioxin_exp3, agegroups=="35-39" & sex=="1", select = totaldioxin_bw)
Dioxin_35_39m <- as.vector(Dioxin_35_39m$totaldioxin_bw)


fit1_35_39m <- fitdist(Dioxin_35_39m, "lnorm")
t <- Dioxin_35_39m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_35_39m$estimate[1], sdlog = fit1_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39m, fitnames="lnorm") 


fit2_35_39m <- fitdist(Dioxin_35_39m, 'gamma') 
t2 <- Dioxin_35_39m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_35_39m$estimate[1],fit2_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39m, fitnames="gamma")


cvm.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Anderson-Darling test


#Age 40-44, males
Dioxin_40_44m <- subset(Dioxin_exp3, agegroups=="40-44" & sex=="1", select = totaldioxin_bw)
Dioxin_40_44m <- as.vector(Dioxin_40_44m$totaldioxin_bw)


fit1_40_44m <- fitdist(Dioxin_40_44m, "lnorm")
t <- Dioxin_40_44m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_40_44m$estimate[1], sdlog = fit1_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44m, fitnames="lnorm") 


fit2_40_44m <- fitdist(Dioxin_40_44m, 'gamma') 
t2 <- Dioxin_40_44m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_40_44m$estimate[1],fit2_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44m, fitnames="gamma")


cvm.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Anderson-Darling test


#Age 45-49, males
Dioxin_45_49m <- subset(Dioxin_exp3, agegroups=="45-49" & sex=="1", select = totaldioxin_bw)
Dioxin_45_49m <- as.vector(Dioxin_45_49m$totaldioxin_bw)


fit1_45_49m <- fitdist(Dioxin_45_49m, "lnorm")
t <- Dioxin_45_49m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_45_49m$estimate[1], sdlog = fit1_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49m, fitnames="lnorm") 


fit2_45_49m <- fitdist(Dioxin_45_49m, 'gamma') 
t2 <- Dioxin_45_49m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_45_49m$estimate[1],fit2_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49m, fitnames="gamma")


cvm.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Anderson-Darling test


#Age 50-54, males
Dioxin_50_54m <- subset(Dioxin_exp3, agegroups=="50-54" & sex=="1", select = totaldioxin_bw)
Dioxin_50_54m <- as.vector(Dioxin_50_54m$totaldioxin_bw)


fit1_50_54m <- fitdist(Dioxin_50_54m, "lnorm")
t <- Dioxin_50_54m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_50_54m$estimate[1], sdlog = fit1_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54m, fitnames="lnorm") 


fit2_50_54m <- fitdist(Dioxin_50_54m, 'gamma') 
t2 <- Dioxin_50_54m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_50_54m$estimate[1],fit2_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54m, fitnames="gamma")


cvm.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Anderson-Darling test


#Age 55-59, males
Dioxin_55_59m <- subset(Dioxin_exp3, agegroups=="55-59" & sex=="1", select = totaldioxin_bw)
Dioxin_55_59m <- as.vector(Dioxin_55_59m$totaldioxin_bw)


fit1_55_59m <- fitdist(Dioxin_55_59m, "lnorm")
t <- Dioxin_55_59m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_55_59m$estimate[1], sdlog = fit1_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59m, fitnames="lnorm") 


fit2_55_59m <- fitdist(Dioxin_55_59m, 'gamma') 
t2 <- Dioxin_55_59m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_55_59m$estimate[1],fit2_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59m, fitnames="gamma")


cvm.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Anderson-Darling test



#Age 60-64, males
Dioxin_60_64m <- subset(Dioxin_exp3, agegroups=="60-64" & sex=="1", select = totaldioxin_bw)
Dioxin_60_64m <- as.vector(Dioxin_60_64m$totaldioxin_bw)


fit1_60_64m <- fitdist(Dioxin_60_64m, "lnorm")
t <- Dioxin_60_64m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_60_64m$estimate[1], sdlog = fit1_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64m, fitnames="lnorm") 


fit2_60_64m <- fitdist(Dioxin_60_64m, 'gamma') 
t2 <- Dioxin_60_64m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_60_64m$estimate[1],fit2_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64m, fitnames="gamma")


cvm.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Anderson-Darling test


#Age 65-69, males
Dioxin_65_69m <- subset(Dioxin_exp3, agegroups=="65-69" & sex=="1", select = totaldioxin_bw)
Dioxin_65_69m <- as.vector(Dioxin_65_69m$totaldioxin_bw)


fit1_65_69m <- fitdist(Dioxin_65_69m, "lnorm")
t <- Dioxin_65_69m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_65_69m$estimate[1], sdlog = fit1_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69m, fitnames="lnorm") 


fit2_65_69m <- fitdist(Dioxin_65_69m, 'gamma') 
t2 <- Dioxin_65_69m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_65_69m$estimate[1],fit2_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69m, fitnames="gamma")


cvm.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Anderson-Darling test



#Age 70-74, males
Dioxin_70_74m <- subset(Dioxin_exp3, agegroups=="70-74" & sex=="1", select = totaldioxin_bw)
Dioxin_70_74m <- as.vector(Dioxin_70_74m$totaldioxin_bw)


fit1_70_74m <- fitdist(Dioxin_70_74m, "lnorm")
t <- Dioxin_70_74m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_70_74m$estimate[1], sdlog = fit1_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74m, fitnames="lnorm") 


fit2_70_74m <- fitdist(Dioxin_70_74m, 'gamma') 
t2 <- Dioxin_70_74m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_70_74m$estimate[1],fit2_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74m, fitnames="gamma")


cvm.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Anderson-Darling test



#Age 75-79, males
Dioxin_75_79m <- subset(Dioxin_exp3, agegroups=="75-79" & sex=="1", select = totaldioxin_bw)
Dioxin_75_79m <- as.vector(Dioxin_75_79m$totaldioxin_bw)


fit1_75_79m <- fitdist(Dioxin_75_79m, "lnorm")
t <- Dioxin_75_79m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_75_79m$estimate[1], sdlog = fit1_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79m, fitnames="lnorm") 


fit2_75_79m <- fitdist(Dioxin_75_79m, 'gamma') 
t2 <- Dioxin_75_79m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_75_79m$estimate[1],fit2_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79m, fitnames="gamma")


cvm.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Anderson-Darling test



##### Fitting exposures female #####


#Age 15-19, female
Dioxin_15_19w <- subset(Dioxin_exp3, agegroups=="15-19" & sex=="2", select = totaldioxin_bw)
Dioxin_15_19w <- as.vector(Dioxin_15_19w$totaldioxin_bw)

fit1_15_19w <- fitdist(Dioxin_15_19w, "lnorm")
t <- Dioxin_15_19w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_15_19w$estimate[1], sdlog = fit1_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19w, fitnames="lnorm") 


fit2_15_19w <- fitdist(Dioxin_15_19w, 'gamma') 
t2 <- Dioxin_15_19w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_15_19w$estimate[1],fit2_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19w, fitnames="gamma")

cvm.test(t, plnorm, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2])  #Anderson-Darling test


#Age 20-24, female
Dioxin_20_24w <- subset(Dioxin_exp3, agegroups=="20-24" & sex=="2", select = totaldioxin_bw)
Dioxin_20_24w <- as.vector(Dioxin_20_24w$totaldioxin_bw)

fit1_20_24w <- fitdist(Dioxin_20_24w, "lnorm")
t <- Dioxin_20_24w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_20_24w$estimate[1], sdlog = fit1_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24w, fitnames="lnorm") 

fit2_20_24w <- fitdist(Dioxin_20_24w, 'gamma') 
t2 <- Dioxin_20_24w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_20_24w$estimate[1],fit2_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24w, fitnames="gamma")

cvm.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Anderson-Darling test


#Age 25-29, female
Dioxin_25_29w <- subset(Dioxin_exp3, agegroups=="25-29" & sex=="2", select = totaldioxin_bw)
Dioxin_25_29w <- as.vector(Dioxin_25_29w$totaldioxin_bw)

fit1_25_29w <- fitdist(Dioxin_25_29w, "lnorm")
t <- Dioxin_25_29w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_25_29w$estimate[1], sdlog = fit1_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29w, fitnames="lnorm") 

fit2_25_29w <- fitdist(Dioxin_25_29w, 'gamma')
t2 <- Dioxin_25_29w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_25_29w$estimate[1],fit2_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29w, fitnames="gamma")

cvm.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Anderson-Darling test


#Age 30-34, female
Dioxin_30_34w <- subset(Dioxin_exp3, agegroups=="30-34" & sex=="2", select = totaldioxin_bw)
Dioxin_30_34w <- as.vector(Dioxin_30_34w$totaldioxin_bw)

fit1_30_34w <- fitdist(Dioxin_30_34w, "lnorm")
t <- Dioxin_30_34w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_30_34w$estimate[1], sdlog = fit1_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34w, fitnames="lnorm") 

fit2_30_34w <- fitdist(Dioxin_30_34w, 'gamma') 
t2 <- Dioxin_30_34w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_30_34w$estimate[1],fit2_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34w, fitnames="gamma")

cvm.test(t, plnorm, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2])  #Anderson-Darling test


#Age 35-39, female
Dioxin_35_39w <- subset(Dioxin_exp3, agegroups=="35-39" & sex=="2", select = totaldioxin_bw)
Dioxin_35_39w <- as.vector(Dioxin_35_39w$totaldioxin_bw)

fit1_35_39w <- fitdist(Dioxin_35_39w, "lnorm")
t <- Dioxin_35_39w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_35_39w$estimate[1], sdlog = fit1_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39w, fitnames="lnorm") 

fit2_35_39w <- fitdist(Dioxin_35_39w, 'gamma') 
t2 <- Dioxin_35_39w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_35_39w$estimate[1],fit2_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39w, fitnames="gamma")

cvm.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Anderson-Darling test


#Age 40-44, female
Dioxin_40_44w <- subset(Dioxin_exp3, agegroups=="40-44" & sex=="2", select = totaldioxin_bw)
Dioxin_40_44w <- as.vector(Dioxin_40_44w$totaldioxin_bw)

fit1_40_44w <- fitdist(Dioxin_40_44w, "lnorm")
t <- Dioxin_40_44w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_40_44w$estimate[1], sdlog = fit1_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44w, fitnames="lnorm") 

fit2_40_44w <- fitdist(Dioxin_40_44w, 'gamma')
t2 <- Dioxin_40_44w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_40_44w$estimate[1],fit2_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44w, fitnames="gamma")

cvm.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Anderson-Darling test


#Age 45-49, female
Dioxin_45_49w <- subset(Dioxin_exp3, agegroups=="45-49" & sex=="2", select = totaldioxin_bw)
Dioxin_45_49w <- as.vector(Dioxin_45_49w$totaldioxin_bw)

fit1_45_49w <- fitdist(Dioxin_45_49w, "lnorm")
t <- Dioxin_45_49w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_45_49w$estimate[1], sdlog = fit1_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49w, fitnames="lnorm") 

fit2_45_49w <- fitdist(Dioxin_45_49w, 'gamma')
t2 <- Dioxin_45_49w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_45_49w$estimate[1],fit2_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49w, fitnames="gamma")

cvm.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Anderson-Darling test


#Age 50-54, female
Dioxin_50_54w <- subset(Dioxin_exp3, agegroups=="50-54" & sex=="2", select = totaldioxin_bw)
Dioxin_50_54w <- as.vector(Dioxin_50_54w$totaldioxin_bw)

fit1_50_54w <- fitdist(Dioxin_50_54w, "lnorm")
t <- Dioxin_50_54w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_50_54w$estimate[1], sdlog = fit1_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54w, fitnames="lnorm") 

fit2_50_54w <- fitdist(Dioxin_50_54w, 'gamma')
t2 <- Dioxin_50_54w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_50_54w$estimate[1],fit2_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54w, fitnames="gamma")

cvm.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Anderson-Darling test


#Age 55-59, female
Dioxin_55_59w <- subset(Dioxin_exp3, agegroups=="55-59" & sex=="2", select = totaldioxin_bw)
Dioxin_55_59w <- as.vector(Dioxin_55_59w$totaldioxin_bw)

fit1_55_59w <- fitdist(Dioxin_55_59w, "lnorm")
t <- Dioxin_55_59w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_55_59w$estimate[1], sdlog = fit1_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59w, fitnames="lnorm") 

fit2_55_59w <- fitdist(Dioxin_55_59w, 'gamma')
t2 <- Dioxin_55_59w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_55_59w$estimate[1],fit2_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59w, fitnames="gamma")

cvm.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Anderson-Darling test



#Age 60-64, female
Dioxin_60_64w <- subset(Dioxin_exp3, agegroups=="60-64" & sex=="2", select = totaldioxin_bw)
Dioxin_60_64w <- as.vector(Dioxin_60_64w$totaldioxin_bw)

fit1_60_64w <- fitdist(Dioxin_60_64w, "lnorm")
t <- Dioxin_60_64w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_60_64w$estimate[1], sdlog = fit1_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64w, fitnames="lnorm") 

fit2_60_64w <- fitdist(Dioxin_60_64w, 'gamma')
t2 <- Dioxin_60_64w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_60_64w$estimate[1],fit2_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64w, fitnames="gamma")

cvm.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Anderson-Darling test


#Age 65-69, female
Dioxin_65_69w <- subset(Dioxin_exp3, agegroups=="65-69" & sex=="2", select = totaldioxin_bw)
Dioxin_65_69w <- as.vector(Dioxin_65_69w$totaldioxin_bw)

fit1_65_69w <- fitdist(Dioxin_65_69w, "lnorm")
t <- Dioxin_65_69w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_65_69w$estimate[1], sdlog = fit1_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69w, fitnames="lnorm") 

fit2_65_69w <- fitdist(Dioxin_65_69w, 'gamma')
t2 <- Dioxin_65_69w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_65_69w$estimate[1],fit2_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69w, fitnames="gamma")

cvm.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Anderson-Darling test


#Age 70-74, female
Dioxin_70_74w <- subset(Dioxin_exp3, agegroups=="70-74" & sex=="2", select = totaldioxin_bw)
Dioxin_70_74w <- as.vector(Dioxin_70_74w$totaldioxin_bw)

fit1_70_74w <- fitdist(Dioxin_70_74w, "lnorm")
t <- Dioxin_70_74w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_70_74w$estimate[1], sdlog = fit1_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74w, fitnames="lnorm") 

fit2_70_74w <- fitdist(Dioxin_70_74w, 'gamma')
t2 <- Dioxin_70_74w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_70_74w$estimate[1],fit2_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74w, fitnames="gamma")

cvm.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Anderson-Darling test


#Age 75-79, female
Dioxin_75_79w <- subset(Dioxin_exp3, agegroups=="75-79" & sex=="2", select = totaldioxin_bw)
Dioxin_75_79w <- as.vector(Dioxin_75_79w$totaldioxin_bw)

fit1_75_79w <- fitdist(Dioxin_75_79w, "lnorm")
t <- Dioxin_75_79w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_75_79w$estimate[1], sdlog = fit1_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79w, fitnames="lnorm") 

fit2_75_79w <- fitdist(Dioxin_75_79w, 'gamma')
t2 <- Dioxin_75_79w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_75_79w$estimate[1],fit2_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79w, fitnames="gamma")

cvm.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Anderson-Darling test


##### Prob of hypothyroidism males #####


# Assume exposure is variable - Lognormal has the best fit


# Age 15-19, males
set.seed(1)
exp <- rlnorm(nvar, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2]) 

pTT4_15_19m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 20-24, males
set.seed(1)
exp <- rlnorm(nvar, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2]) 

pTT4_20_24m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 25-29, males
set.seed(1)
exp <- rlnorm(nvar, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2]) 

pTT4_25_29m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 30-34, males
set.seed(1)
exp <- rlnorm(nvar, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2]) 

pTT4_30_34m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 35-39, males
set.seed(1)
exp <- rlnorm(nvar, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2]) 

pTT4_35_39m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 40-44, males
set.seed(1)
exp <- rlnorm(nvar, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2]) 

pTT4_40_44m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 45-49, males
set.seed(1)
exp <- rlnorm(nvar, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2]) 

pTT4_45_49m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 50-54, males
set.seed(1)
exp <- rlnorm(nvar, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2]) 

pTT4_50_54m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 55-59, males
set.seed(1)
exp <- rlnorm(nvar, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2]) 

pTT4_55_59m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 60-64, males
set.seed(1)
exp <- rlnorm(nvar, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2]) 

pTT4_60_64m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 65-69, males
set.seed(1)
exp <- rlnorm(nvar, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2]) 

pTT4_65_69m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 70-74, males
set.seed(1)
exp <- rlnorm(nvar, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2]) 

pTT4_70_74m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 75-79, males
set.seed(1)
exp <- rlnorm(nvar, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2]) 

pTT4_75_79m <- apply(ICED, 2, function(x) mean(x < exp))


#Assume age 80-84 and 85+ has the same prob as 75-79


##### Prob of hypothyroidism females #####


# Assume exposure is variable

# Age 15-19, females
set.seed(1)
exp <- rlnorm(nvar, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2]) 

pTT4_15_19w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 20-24, females
set.seed(1)
exp <- rlnorm(nvar, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2]) 

pTT4_20_24w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 25-29, females
set.seed(1)
exp <- rlnorm(nvar, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2]) 

pTT4_25_29w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 30-34, females
set.seed(1)
exp <- rlnorm(nvar, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2]) 

pTT4_30_34w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 35-39, females
set.seed(1)
exp <- rlnorm(nvar, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2]) 

pTT4_35_39w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 40-44, females
set.seed(1)
exp <- rlnorm(nvar, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2]) 

pTT4_40_44w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 45-49, females
set.seed(1)
exp <- rlnorm(nvar, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2]) 

pTT4_45_49w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 50-54, females
set.seed(1)
exp <- rlnorm(nvar, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2]) 

pTT4_50_54w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 55-59, females
set.seed(1)
exp <- rlnorm(nvar, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2]) 

pTT4_55_59w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 60-64, females
set.seed(1)
exp <- rlnorm(nvar, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2]) 

pTT4_60_64w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 65-69, females
set.seed(1)
exp <- rlnorm(nvar, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2]) 

pTT4_65_69w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 70-74, females
set.seed(1)
exp <- rlnorm(nvar, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2]) 

pTT4_70_74w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 75-79, females
set.seed(1)
exp <- rlnorm(nvar, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2]) 

pTT4_75_79w <- apply(ICED, 2, function(x) mean(x < exp))


#Assume age 80-84 and 85+ has the same prob as 75-79


##### DALY calculation #####

##### Males #####

#Age 15-19, male

DALY15_19m_scen3 <- pTT4_15_19m * pop_15_19m * D * dw
mean_median_ci(DALY15_19m_scen3)


#Age 20-24y, male

DALY20_24m_scen3 <- pTT4_20_24m * pop_20_24m * D * dw
mean_median_ci(DALY20_24m_scen3)


#Age 25-29y, male

DALY25_29m_scen3 <- pTT4_25_29m * pop_25_29m * D * dw
mean_median_ci(DALY25_29m_scen3)


#Age 30-34y, male

DALY30_34m_scen3 <- pTT4_30_34m * pop_30_34m * D * dw
mean_median_ci(DALY30_34m_scen3)


#Age 35-39y, male

DALY35_39m_scen3 <- pTT4_35_39m * pop_35_39m * D * dw
mean_median_ci(DALY35_39m_scen3)



#Age 40-44y, male

DALY40_44m_scen3 <- pTT4_40_44m * pop_40_44m * D * dw
mean_median_ci(DALY40_44m_scen3)


#Age 45-49y, male

DALY45_49m_scen3 <- pTT4_45_49m * pop_45_49m * D * dw
mean_median_ci(DALY45_49m_scen3)


#Age 50-54y, male

DALY50_54m_scen3 <- pTT4_50_54m * pop_50_54m * D * dw
mean_median_ci(DALY50_54m_scen3)


#Age 55-59y, male

DALY55_59m_scen3 <- pTT4_55_59m * pop_55_59m * D * dw
mean_median_ci(DALY55_59m_scen3)


#Age 60-64y, male

DALY60_64m_scen3 <- pTT4_60_64m * pop_60_64m * D * dw
mean_median_ci(DALY60_64m_scen3)


#Age 65-69y, male

DALY65_69m_scen3 <- pTT4_65_69m * pop_65_69m * D * dw
mean_median_ci(DALY65_69m_scen3)


#Age 70-74y, male

DALY70_74m_scen3 <- pTT4_70_74m * pop_70_74m * D * dw
mean_median_ci(DALY70_74m_scen3)


#Age 75-79y, male

DALY75_79m_scen3 <- pTT4_75_79m * pop_75_79m * D * dw
mean_median_ci(DALY75_79m_scen3)


#Age 80-84y, male

DALY80_84m_scen3 <- pTT4_75_79m * pop_80_84m * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY80_84m_scen3)


#Age 85+y, male

DALY85m_scen3 <- pTT4_75_79m * pop_85m * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY85m_scen3)


##### Females #####

#Age 15-19, female

DALY15_19w_scen3 <- pTT4_15_19w * pop_15_19w * D * dw
mean_median_ci(DALY15_19w_scen3)


#Age 20-24y, female

DALY20_24w_scen3 <- pTT4_20_24w * pop_20_24w * D * dw
mean_median_ci(DALY20_24w_scen3)


#Age 25-29y, female

DALY25_29w_scen3 <- pTT4_25_29w * pop_25_29w * D * dw
mean_median_ci(DALY25_29w_scen3)


#Age 30-34y, female

DALY30_34w_scen3 <- pTT4_30_34w * pop_30_34w * D * dw
mean_median_ci(DALY30_34w_scen3)


#Age 35-39y, female

DALY35_39w_scen3 <- pTT4_35_39w * pop_35_39w * D * dw
mean_median_ci(DALY35_39w_scen3)



#Age 40-44y, female

DALY40_44w_scen3 <- pTT4_40_44w * pop_40_44w * D * dw
mean_median_ci(DALY40_44w_scen3)


#Age 45-49y, female

DALY45_49w_scen3 <- pTT4_45_49w * pop_45_49w * D * dw
mean_median_ci(DALY45_49w_scen3)


#Age 50-54y, female

DALY50_54w_scen3 <- pTT4_50_54w * pop_50_54w * D * dw
mean_median_ci(DALY50_54w_scen3)


#Age 55-59y, female

DALY55_59w_scen3 <- pTT4_55_59w * pop_55_59w * D * dw
mean_median_ci(DALY55_59w_scen3)


#Age 60-64y, female

DALY60_64w_scen3 <- pTT4_60_64w * pop_60_64w * D * dw
mean_median_ci(DALY60_64w_scen3)


#Age 65-69y, female

DALY65_69w_scen3 <- pTT4_65_69w * pop_65_69w * D * dw
mean_median_ci(DALY65_69w_scen3)


#Age 70-74y, female

DALY70_74w_scen3 <- pTT4_70_74w * pop_70_74w * D * dw
mean_median_ci(DALY70_74w_scen3)


#Age 75-79y, female

DALY75_79w_scen3 <- pTT4_75_79w * pop_75_79w * D * dw
mean_median_ci(DALY75_79w_scen3)


#Age 80-84y, female

DALY80_84w_scen3 <- pTT4_75_79w * pop_80_84w * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY80_84w_scen3)


#Age 85+y, female

DALY85w_scen3 <- pTT4_75_79w * pop_85w * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY85w_scen3)


##### Total DALYs scen 3 #####

DALYtotal_scen3 <- DALY15_19m_scen3 + DALY20_24m_scen3 + DALY25_29m_scen3 + DALY30_34m_scen3 + DALY35_39m_scen3 +
  DALY40_44m_scen3 + DALY45_49m_scen3 + DALY50_54m_scen3 +  DALY55_59m_scen3 + DALY60_64m_scen3 + DALY65_69m_scen3 +
  DALY70_74m_scen3 + DALY75_79m_scen3 + DALY80_84m_scen3 + DALY85m_scen3 +
  DALY15_19w_scen3 + DALY20_24w_scen3 + DALY25_29w_scen3 + DALY30_34w_scen3 + DALY35_39w_scen3 + DALY40_44w_scen3 +
  DALY45_49w_scen3 + DALY50_54w_scen3 + DALY55_59w_scen3 + DALY60_64w_scen3 + DALY65_69w_scen3 + DALY70_74w_scen3 +
  DALY75_79w_scen3 + DALY80_84w_scen3 + DALY85w_scen3


mean_median_ci(DALYtotal_scen3)


##### Cases scenario 3 #####

cases_scen3 <- pTT4_15_19m * pop_15_19m + pTT4_20_24m * pop_20_24m + pTT4_25_29m * pop_25_29m + pTT4_30_34m * pop_30_34m +
  pTT4_35_39m * pop_35_39m + pTT4_40_44m * pop_40_44m + pTT4_45_49m * pop_45_49m + pTT4_50_54m * pop_50_54m + pTT4_55_59m * pop_55_59m +
  pTT4_60_64m * pop_60_64m + pTT4_65_69m * pop_65_69m + pTT4_70_74m * pop_70_74m + pTT4_75_79m * pop_75_79m + pTT4_75_79m * pop_80_84m +
  pTT4_75_79m * pop_85m +
  pTT4_15_19w * pop_15_19w + pTT4_20_24w * pop_20_24w + pTT4_25_29w * pop_25_29w + pTT4_30_34w * pop_30_34w +
  pTT4_35_39w * pop_35_39w + pTT4_40_44w * pop_40_44w + pTT4_45_49w * pop_45_49w + pTT4_50_54w * pop_50_54w + pTT4_55_59w * pop_55_59w +
  pTT4_60_64w * pop_60_64w + pTT4_65_69w * pop_65_69w + pTT4_70_74w * pop_70_74w + pTT4_75_79w * pop_75_79w + pTT4_75_79w * pop_80_84w +
  pTT4_75_79w * pop_85w


mean_median_ci(cases_scen3)


##### Scenario 4 #####

Dioxin_exp4 <- read.csv("Dioxin_exp_scen4.csv")

setDT(Dioxin_exp4)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Fitting exposures male #####

#Age 15-19, males
Dioxin_15_19m <- subset(Dioxin_exp4, agegroups=="15-19" & sex=="1", select = totaldioxin_bw)
Dioxin_15_19m <- as.vector(Dioxin_15_19m$totaldioxin_bw)


fit1_15_19m <- fitdist(Dioxin_15_19m, "lnorm")
t <- Dioxin_15_19m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_15_19m$estimate[1], sdlog = fit1_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19m, fitnames="lnorm") 


fit2_15_19m <- fitdist(Dioxin_15_19m, 'gamma') 
t2 <- Dioxin_15_19m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_15_19m$estimate[1],fit2_15_19m$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19m, fitnames="gamma")


cvm.test(t, plnorm, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19m$estimate[1], fit2_15_19m$estimate[2])  #Anderson-Darling test


#Age 20-24, males
Dioxin_20_24m <- subset(Dioxin_exp4, agegroups=="20-24" & sex=="1", select = totaldioxin_bw)
Dioxin_20_24m <- as.vector(Dioxin_20_24m$totaldioxin_bw)


fit1_20_24m <- fitdist(Dioxin_20_24m, "lnorm")
t <- Dioxin_20_24m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_20_24m$estimate[1], sdlog = fit1_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24m, fitnames="lnorm") 


fit2_20_24m <- fitdist(Dioxin_20_24m, 'gamma') 
t2 <- Dioxin_20_24m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_20_24m$estimate[1],fit2_20_24m$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24m, fitnames="gamma")


cvm.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24m$estimate[1], fit2_20_24m$estimate[2])  #Anderson-Darling test


#Age 25-29, males
Dioxin_25_29m <- subset(Dioxin_exp4, agegroups=="25-29" & sex=="1", select = totaldioxin_bw)
Dioxin_25_29m <- as.vector(Dioxin_25_29m$totaldioxin_bw)


fit1_25_29m <- fitdist(Dioxin_25_29m, "lnorm")
t <- Dioxin_25_29m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_25_29m$estimate[1], sdlog = fit1_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29m, fitnames="lnorm") 


fit2_25_29m <- fitdist(Dioxin_25_29m, 'gamma') 
t2 <- Dioxin_25_29m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_25_29m$estimate[1],fit2_25_29m$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29m, fitnames="gamma")


cvm.test(t, plnorm, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29m$estimate[1], fit2_25_29m$estimate[2])  #Anderson-Darling test


#Age 30-34, males
Dioxin_30_34m <- subset(Dioxin_exp4, agegroups=="30-34" & sex=="1", select = totaldioxin_bw)
Dioxin_30_34m <- as.vector(Dioxin_30_34m$totaldioxin_bw)


fit1_30_34m <- fitdist(Dioxin_30_34m, "lnorm")
t <- Dioxin_30_34m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_30_34m$estimate[1], sdlog = fit1_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34m, fitnames="lnorm") 


fit2_30_34m <- fitdist(Dioxin_30_34m, 'gamma') 
t2 <- Dioxin_30_34m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_30_34m$estimate[1],fit2_30_34m$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34m, fitnames="gamma")


cvm.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34m$estimate[1], fit2_30_34m$estimate[2])  #Anderson-Darling test


#Age 35-39, males
Dioxin_35_39m <- subset(Dioxin_exp4, agegroups=="35-39" & sex=="1", select = totaldioxin_bw)
Dioxin_35_39m <- as.vector(Dioxin_35_39m$totaldioxin_bw)


fit1_35_39m <- fitdist(Dioxin_35_39m, "lnorm")
t <- Dioxin_35_39m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_35_39m$estimate[1], sdlog = fit1_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39m, fitnames="lnorm") 


fit2_35_39m <- fitdist(Dioxin_35_39m, 'gamma') 
t2 <- Dioxin_35_39m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_35_39m$estimate[1],fit2_35_39m$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39m, fitnames="gamma")


cvm.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39m$estimate[1], fit2_35_39m$estimate[2])  #Anderson-Darling test


#Age 40-44, males
Dioxin_40_44m <- subset(Dioxin_exp4, agegroups=="40-44" & sex=="1", select = totaldioxin_bw)
Dioxin_40_44m <- as.vector(Dioxin_40_44m$totaldioxin_bw)


fit1_40_44m <- fitdist(Dioxin_40_44m, "lnorm")
t <- Dioxin_40_44m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_40_44m$estimate[1], sdlog = fit1_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44m, fitnames="lnorm") 


fit2_40_44m <- fitdist(Dioxin_40_44m, 'gamma') 
t2 <- Dioxin_40_44m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_40_44m$estimate[1],fit2_40_44m$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44m, fitnames="gamma")


cvm.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44m$estimate[1], fit2_40_44m$estimate[2])  #Anderson-Darling test


#Age 45-49, males
Dioxin_45_49m <- subset(Dioxin_exp4, agegroups=="45-49" & sex=="1", select = totaldioxin_bw)
Dioxin_45_49m <- as.vector(Dioxin_45_49m$totaldioxin_bw)


fit1_45_49m <- fitdist(Dioxin_45_49m, "lnorm")
t <- Dioxin_45_49m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_45_49m$estimate[1], sdlog = fit1_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49m, fitnames="lnorm") 


fit2_45_49m <- fitdist(Dioxin_45_49m, 'gamma') 
t2 <- Dioxin_45_49m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_45_49m$estimate[1],fit2_45_49m$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49m, fitnames="gamma")


cvm.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49m$estimate[1], fit2_45_49m$estimate[2])  #Anderson-Darling test


#Age 50-54, males
Dioxin_50_54m <- subset(Dioxin_exp4, agegroups=="50-54" & sex=="1", select = totaldioxin_bw)
Dioxin_50_54m <- as.vector(Dioxin_50_54m$totaldioxin_bw)


fit1_50_54m <- fitdist(Dioxin_50_54m, "lnorm")
t <- Dioxin_50_54m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_50_54m$estimate[1], sdlog = fit1_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54m, fitnames="lnorm") 


fit2_50_54m <- fitdist(Dioxin_50_54m, 'gamma') 
t2 <- Dioxin_50_54m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_50_54m$estimate[1],fit2_50_54m$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54m, fitnames="gamma")


cvm.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54m$estimate[1], fit2_50_54m$estimate[2])  #Anderson-Darling test


#Age 55-59, males
Dioxin_55_59m <- subset(Dioxin_exp4, agegroups=="55-59" & sex=="1", select = totaldioxin_bw)
Dioxin_55_59m <- as.vector(Dioxin_55_59m$totaldioxin_bw)


fit1_55_59m <- fitdist(Dioxin_55_59m, "lnorm")
t <- Dioxin_55_59m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_55_59m$estimate[1], sdlog = fit1_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59m, fitnames="lnorm") 


fit2_55_59m <- fitdist(Dioxin_55_59m, 'gamma') 
t2 <- Dioxin_55_59m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_55_59m$estimate[1],fit2_55_59m$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59m, fitnames="gamma")


cvm.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59m$estimate[1], fit2_55_59m$estimate[2])  #Anderson-Darling test



#Age 60-64, males
Dioxin_60_64m <- subset(Dioxin_exp4, agegroups=="60-64" & sex=="1", select = totaldioxin_bw)
Dioxin_60_64m <- as.vector(Dioxin_60_64m$totaldioxin_bw)


fit1_60_64m <- fitdist(Dioxin_60_64m, "lnorm")
t <- Dioxin_60_64m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_60_64m$estimate[1], sdlog = fit1_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64m, fitnames="lnorm") 


fit2_60_64m <- fitdist(Dioxin_60_64m, 'gamma') 
t2 <- Dioxin_60_64m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_60_64m$estimate[1],fit2_60_64m$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64m, fitnames="gamma")


cvm.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64m$estimate[1], fit2_60_64m$estimate[2])  #Anderson-Darling test


#Age 65-69, males
Dioxin_65_69m <- subset(Dioxin_exp4, agegroups=="65-69" & sex=="1", select = totaldioxin_bw)
Dioxin_65_69m <- as.vector(Dioxin_65_69m$totaldioxin_bw)


fit1_65_69m <- fitdist(Dioxin_65_69m, "lnorm")
t <- Dioxin_65_69m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_65_69m$estimate[1], sdlog = fit1_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69m, fitnames="lnorm") 


fit2_65_69m <- fitdist(Dioxin_65_69m, 'gamma') 
t2 <- Dioxin_65_69m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_65_69m$estimate[1],fit2_65_69m$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69m, fitnames="gamma")


cvm.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69m$estimate[1], fit2_65_69m$estimate[2])  #Anderson-Darling test



#Age 70-74, males
Dioxin_70_74m <- subset(Dioxin_exp4, agegroups=="70-74" & sex=="1", select = totaldioxin_bw)
Dioxin_70_74m <- as.vector(Dioxin_70_74m$totaldioxin_bw)


fit1_70_74m <- fitdist(Dioxin_70_74m, "lnorm")
t <- Dioxin_70_74m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_70_74m$estimate[1], sdlog = fit1_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74m, fitnames="lnorm") 


fit2_70_74m <- fitdist(Dioxin_70_74m, 'gamma') 
t2 <- Dioxin_70_74m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_70_74m$estimate[1],fit2_70_74m$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74m, fitnames="gamma")


cvm.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74m$estimate[1], fit2_70_74m$estimate[2])  #Anderson-Darling test



#Age 75-79, males
Dioxin_75_79m <- subset(Dioxin_exp4, agegroups=="75-79" & sex=="1", select = totaldioxin_bw)
Dioxin_75_79m <- as.vector(Dioxin_75_79m$totaldioxin_bw)


fit1_75_79m <- fitdist(Dioxin_75_79m, "lnorm")
t <- Dioxin_75_79m
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_75_79m$estimate[1], sdlog = fit1_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79m, fitnames="lnorm") 


fit2_75_79m <- fitdist(Dioxin_75_79m, 'gamma') 
t2 <- Dioxin_75_79m
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_75_79m$estimate[1],fit2_75_79m$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79m, fitnames="gamma")


cvm.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79m$estimate[1], fit2_75_79m$estimate[2])  #Anderson-Darling test



##### Fitting exposures female #####


#Age 15-19, female
Dioxin_15_19w <- subset(Dioxin_exp4, agegroups=="15-19" & sex=="2", select = totaldioxin_bw)
Dioxin_15_19w <- as.vector(Dioxin_15_19w$totaldioxin_bw)

fit1_15_19w <- fitdist(Dioxin_15_19w, "lnorm")
t <- Dioxin_15_19w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_15_19w$estimate[1], sdlog = fit1_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19w, fitnames="lnorm") 


fit2_15_19w <- fitdist(Dioxin_15_19w, 'gamma') 
t2 <- Dioxin_15_19w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_15_19w$estimate[1],fit2_15_19w$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19w, fitnames="gamma")

cvm.test(t, plnorm, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19w$estimate[1], fit2_15_19w$estimate[2])  #Anderson-Darling test


#Age 20-24, female
Dioxin_20_24w <- subset(Dioxin_exp4, agegroups=="20-24" & sex=="2", select = totaldioxin_bw)
Dioxin_20_24w <- as.vector(Dioxin_20_24w$totaldioxin_bw)

fit1_20_24w <- fitdist(Dioxin_20_24w, "lnorm")
t <- Dioxin_20_24w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_20_24w$estimate[1], sdlog = fit1_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24w, fitnames="lnorm") 

fit2_20_24w <- fitdist(Dioxin_20_24w, 'gamma') 
t2 <- Dioxin_20_24w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_20_24w$estimate[1],fit2_20_24w$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24w, fitnames="gamma")

cvm.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24w$estimate[1], fit2_20_24w$estimate[2])  #Anderson-Darling test


#Age 25-29, female
Dioxin_25_29w <- subset(Dioxin_exp4, agegroups=="25-29" & sex=="2", select = totaldioxin_bw)
Dioxin_25_29w <- as.vector(Dioxin_25_29w$totaldioxin_bw)

fit1_25_29w <- fitdist(Dioxin_25_29w, "lnorm")
t <- Dioxin_25_29w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_25_29w$estimate[1], sdlog = fit1_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29w, fitnames="lnorm") 

fit2_25_29w <- fitdist(Dioxin_25_29w, 'gamma')
t2 <- Dioxin_25_29w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_25_29w$estimate[1],fit2_25_29w$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29w, fitnames="gamma")

cvm.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29w$estimate[1], fit2_25_29w$estimate[2])  #Anderson-Darling test


#Age 30-34, female
Dioxin_30_34w <- subset(Dioxin_exp4, agegroups=="30-34" & sex=="2", select = totaldioxin_bw)
Dioxin_30_34w <- as.vector(Dioxin_30_34w$totaldioxin_bw)

fit1_30_34w <- fitdist(Dioxin_30_34w, "lnorm")
t <- Dioxin_30_34w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_30_34w$estimate[1], sdlog = fit1_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34w, fitnames="lnorm") 

fit2_30_34w <- fitdist(Dioxin_30_34w, 'gamma') 
t2 <- Dioxin_30_34w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_30_34w$estimate[1],fit2_30_34w$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34w, fitnames="gamma")

cvm.test(t, plnorm, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34w$estimate[1], fit2_30_34w$estimate[2])  #Anderson-Darling test


#Age 35-39, female
Dioxin_35_39w <- subset(Dioxin_exp4, agegroups=="35-39" & sex=="2", select = totaldioxin_bw)
Dioxin_35_39w <- as.vector(Dioxin_35_39w$totaldioxin_bw)

fit1_35_39w <- fitdist(Dioxin_35_39w, "lnorm")
t <- Dioxin_35_39w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_35_39w$estimate[1], sdlog = fit1_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39w, fitnames="lnorm") 

fit2_35_39w <- fitdist(Dioxin_35_39w, 'gamma') 
t2 <- Dioxin_35_39w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_35_39w$estimate[1],fit2_35_39w$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39w, fitnames="gamma")

cvm.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39w$estimate[1], fit2_35_39w$estimate[2])  #Anderson-Darling test


#Age 40-44, female
Dioxin_40_44w <- subset(Dioxin_exp4, agegroups=="40-44" & sex=="2", select = totaldioxin_bw)
Dioxin_40_44w <- as.vector(Dioxin_40_44w$totaldioxin_bw)

fit1_40_44w <- fitdist(Dioxin_40_44w, "lnorm")
t <- Dioxin_40_44w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_40_44w$estimate[1], sdlog = fit1_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44w, fitnames="lnorm") 

fit2_40_44w <- fitdist(Dioxin_40_44w, 'gamma')
t2 <- Dioxin_40_44w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_40_44w$estimate[1],fit2_40_44w$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44w, fitnames="gamma")

cvm.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44w$estimate[1], fit2_40_44w$estimate[2])  #Anderson-Darling test


#Age 45-49, female
Dioxin_45_49w <- subset(Dioxin_exp4, agegroups=="45-49" & sex=="2", select = totaldioxin_bw)
Dioxin_45_49w <- as.vector(Dioxin_45_49w$totaldioxin_bw)

fit1_45_49w <- fitdist(Dioxin_45_49w, "lnorm")
t <- Dioxin_45_49w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_45_49w$estimate[1], sdlog = fit1_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49w, fitnames="lnorm") 

fit2_45_49w <- fitdist(Dioxin_45_49w, 'gamma')
t2 <- Dioxin_45_49w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_45_49w$estimate[1],fit2_45_49w$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49w, fitnames="gamma")

cvm.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49w$estimate[1], fit2_45_49w$estimate[2])  #Anderson-Darling test


#Age 50-54, female
Dioxin_50_54w <- subset(Dioxin_exp4, agegroups=="50-54" & sex=="2", select = totaldioxin_bw)
Dioxin_50_54w <- as.vector(Dioxin_50_54w$totaldioxin_bw)

fit1_50_54w <- fitdist(Dioxin_50_54w, "lnorm")
t <- Dioxin_50_54w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_50_54w$estimate[1], sdlog = fit1_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit1_50_54w, fitnames="lnorm") 

fit2_50_54w <- fitdist(Dioxin_50_54w, 'gamma')
t2 <- Dioxin_50_54w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_50_54w$estimate[1],fit2_50_54w$estimate[2]), col="red", lty=1)
gofstat(fit2_50_54w, fitnames="gamma")

cvm.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_50_54w$estimate[1], fit2_50_54w$estimate[2])  #Anderson-Darling test


#Age 55-59, female
Dioxin_55_59w <- subset(Dioxin_exp4, agegroups=="55-59" & sex=="2", select = totaldioxin_bw)
Dioxin_55_59w <- as.vector(Dioxin_55_59w$totaldioxin_bw)

fit1_55_59w <- fitdist(Dioxin_55_59w, "lnorm")
t <- Dioxin_55_59w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_55_59w$estimate[1], sdlog = fit1_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit1_55_59w, fitnames="lnorm") 

fit2_55_59w <- fitdist(Dioxin_55_59w, 'gamma')
t2 <- Dioxin_55_59w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_55_59w$estimate[1],fit2_55_59w$estimate[2]), col="red", lty=1)
gofstat(fit2_55_59w, fitnames="gamma")

cvm.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_55_59w$estimate[1], fit2_55_59w$estimate[2])  #Anderson-Darling test



#Age 60-64, female
Dioxin_60_64w <- subset(Dioxin_exp4, agegroups=="60-64" & sex=="2", select = totaldioxin_bw)
Dioxin_60_64w <- as.vector(Dioxin_60_64w$totaldioxin_bw)

fit1_60_64w <- fitdist(Dioxin_60_64w, "lnorm")
t <- Dioxin_60_64w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_60_64w$estimate[1], sdlog = fit1_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit1_60_64w, fitnames="lnorm") 

fit2_60_64w <- fitdist(Dioxin_60_64w, 'gamma')
t2 <- Dioxin_60_64w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_60_64w$estimate[1],fit2_60_64w$estimate[2]), col="red", lty=1)
gofstat(fit2_60_64w, fitnames="gamma")

cvm.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_60_64w$estimate[1], fit2_60_64w$estimate[2])  #Anderson-Darling test


#Age 65-69, female
Dioxin_65_69w <- subset(Dioxin_exp4, agegroups=="65-69" & sex=="2", select = totaldioxin_bw)
Dioxin_65_69w <- as.vector(Dioxin_65_69w$totaldioxin_bw)

fit1_65_69w <- fitdist(Dioxin_65_69w, "lnorm")
t <- Dioxin_65_69w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_65_69w$estimate[1], sdlog = fit1_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit1_65_69w, fitnames="lnorm") 

fit2_65_69w <- fitdist(Dioxin_65_69w, 'gamma')
t2 <- Dioxin_65_69w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_65_69w$estimate[1],fit2_65_69w$estimate[2]), col="red", lty=1)
gofstat(fit2_65_69w, fitnames="gamma")

cvm.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_65_69w$estimate[1], fit2_65_69w$estimate[2])  #Anderson-Darling test


#Age 70-74, female
Dioxin_70_74w <- subset(Dioxin_exp4, agegroups=="70-74" & sex=="2", select = totaldioxin_bw)
Dioxin_70_74w <- as.vector(Dioxin_70_74w$totaldioxin_bw)

fit1_70_74w <- fitdist(Dioxin_70_74w, "lnorm")
t <- Dioxin_70_74w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_70_74w$estimate[1], sdlog = fit1_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit1_70_74w, fitnames="lnorm") 

fit2_70_74w <- fitdist(Dioxin_70_74w, 'gamma')
t2 <- Dioxin_70_74w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_70_74w$estimate[1],fit2_70_74w$estimate[2]), col="red", lty=1)
gofstat(fit2_70_74w, fitnames="gamma")

cvm.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_70_74w$estimate[1], fit2_70_74w$estimate[2])  #Anderson-Darling test


#Age 75-79, female
Dioxin_75_79w <- subset(Dioxin_exp4, agegroups=="75-79" & sex=="2", select = totaldioxin_bw)
Dioxin_75_79w <- as.vector(Dioxin_75_79w$totaldioxin_bw)

fit1_75_79w <- fitdist(Dioxin_75_79w, "lnorm")
t <- Dioxin_75_79w
plot(ecdf(t), lty=1)
x <- seq(0,4, length=1000)
lines(x, plnorm(x, meanlog = fit1_75_79w$estimate[1], sdlog = fit1_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit1_75_79w, fitnames="lnorm") 

fit2_75_79w <- fitdist(Dioxin_75_79w, 'gamma')
t2 <- Dioxin_75_79w
plot(ecdf(t2), lty=1)
x <- seq(0, 4, length=1000)
lines(x, pgamma(x, fit2_75_79w$estimate[1],fit2_75_79w$estimate[2]), col="red", lty=1)
gofstat(fit2_75_79w, fitnames="gamma")

cvm.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_75_79w$estimate[1], fit2_75_79w$estimate[2])  #Anderson-Darling test



##### Prob of hypothyroidism males #####


# Assume exposure is variable

# Age 15-19, males
set.seed(1)
exp <- rlnorm(nvar, fit1_15_19m$estimate[1], fit1_15_19m$estimate[2]) 

pTT4_15_19m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 20-24, males
set.seed(1)
exp <- rlnorm(nvar, fit1_20_24m$estimate[1], fit1_20_24m$estimate[2]) 

pTT4_20_24m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 25-29, males
set.seed(1)
exp <- rlnorm(nvar, fit1_25_29m$estimate[1], fit1_25_29m$estimate[2]) 

pTT4_25_29m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 30-34, males
set.seed(1)
exp <- rlnorm(nvar, fit1_30_34m$estimate[1], fit1_30_34m$estimate[2]) 

pTT4_30_34m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 35-39, males
set.seed(1)
exp <- rlnorm(nvar, fit1_35_39m$estimate[1], fit1_35_39m$estimate[2]) 

pTT4_35_39m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 40-44, males
set.seed(1)
exp <- rlnorm(nvar, fit1_40_44m$estimate[1], fit1_40_44m$estimate[2]) 

pTT4_40_44m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 45-49, males
set.seed(1)
exp <- rlnorm(nvar, fit1_45_49m$estimate[1], fit1_45_49m$estimate[2]) 

pTT4_45_49m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 50-54, males
set.seed(1)
exp <- rlnorm(nvar, fit1_50_54m$estimate[1], fit1_50_54m$estimate[2]) 

pTT4_50_54m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 55-59, males
set.seed(1)
exp <- rlnorm(nvar, fit1_55_59m$estimate[1], fit1_55_59m$estimate[2]) 

pTT4_55_59m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 60-64, males
set.seed(1)
exp <- rlnorm(nvar, fit1_60_64m$estimate[1], fit1_60_64m$estimate[2]) 

pTT4_60_64m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 65-69, males
set.seed(1)
exp <- rlnorm(nvar, fit1_65_69m$estimate[1], fit1_65_69m$estimate[2]) 

pTT4_65_69m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 70-74, males
set.seed(1)
exp <- rlnorm(nvar, fit1_70_74m$estimate[1], fit1_70_74m$estimate[2]) 

pTT4_70_74m <- apply(ICED, 2, function(x) mean(x < exp))



# Age 75-79, males
set.seed(1)
exp <- rlnorm(nvar, fit1_75_79m$estimate[1], fit1_75_79m$estimate[2]) 

pTT4_75_79m <- apply(ICED, 2, function(x) mean(x < exp))


#Assume age 80-84 and 85+ has the same prob as 75-79



##### Prob of hypothyroidism females #####


# Assume exposure is variable

# Age 15-19, females
set.seed(1)
exp <- rlnorm(nvar, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2]) 

pTT4_15_19w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 20-24, females
set.seed(1)
exp <- rlnorm(nvar, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2]) 

pTT4_20_24w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 25-29, females
set.seed(1)
exp <- rlnorm(nvar, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2]) 

pTT4_25_29w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 30-34, females
set.seed(1)
exp <- rlnorm(nvar, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2]) 

pTT4_30_34w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 35-39, females
set.seed(1)
exp <- rlnorm(nvar, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2]) 

pTT4_35_39w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 40-44, females
set.seed(1)
exp <- rlnorm(nvar, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2]) 

pTT4_40_44w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 45-49, females
set.seed(1)
exp <- rlnorm(nvar, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2]) 

pTT4_45_49w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 50-54, females
set.seed(1)
exp <- rlnorm(nvar, fit1_50_54w$estimate[1], fit1_50_54w$estimate[2]) 

pTT4_50_54w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 55-59, females
set.seed(1)
exp <- rlnorm(nvar, fit1_55_59w$estimate[1], fit1_55_59w$estimate[2]) 

pTT4_55_59w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 60-64, females
set.seed(1)
exp <- rlnorm(nvar, fit1_60_64w$estimate[1], fit1_60_64w$estimate[2]) 

pTT4_60_64w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 65-69, females
set.seed(1)
exp <- rlnorm(nvar, fit1_65_69w$estimate[1], fit1_65_69w$estimate[2]) 

pTT4_65_69w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 70-74, females
set.seed(1)
exp <- rlnorm(nvar, fit1_70_74w$estimate[1], fit1_70_74w$estimate[2]) 

pTT4_70_74w <- apply(ICED, 2, function(x) mean(x < exp))



# Age 75-79, females
set.seed(1)
exp <- rlnorm(nvar, fit1_75_79w$estimate[1], fit1_75_79w$estimate[2]) 

pTT4_75_79w <- apply(ICED, 2, function(x) mean(x < exp))



#Assume age 80-84 and 85+ has the same prob as 75-79



##### DALY calculation #####

##### Males #####

#Age 15-19, male

DALY15_19m_scen4 <- pTT4_15_19m * pop_15_19m * D * dw
mean_median_ci(DALY15_19m_scen4)


#Age 20-24y, male

DALY20_24m_scen4 <- pTT4_20_24m * pop_20_24m * D * dw
mean_median_ci(DALY20_24m_scen4)


#Age 25-29y, male

DALY25_29m_scen4 <- pTT4_25_29m * pop_25_29m * D * dw
mean_median_ci(DALY25_29m_scen4)


#Age 30-34y, male

DALY30_34m_scen4 <- pTT4_30_34m * pop_30_34m * D * dw
mean_median_ci(DALY30_34m_scen4)


#Age 35-39y, male

DALY35_39m_scen4 <- pTT4_35_39m * pop_35_39m * D * dw
mean_median_ci(DALY35_39m_scen4)



#Age 40-44y, male

DALY40_44m_scen4 <- pTT4_40_44m * pop_40_44m * D * dw
mean_median_ci(DALY40_44m_scen4)


#Age 45-49y, male

DALY45_49m_scen4 <- pTT4_45_49m * pop_45_49m * D * dw
mean_median_ci(DALY45_49m_scen4)


#Age 50-54y, male

DALY50_54m_scen4 <- pTT4_50_54m * pop_50_54m * D * dw
mean_median_ci(DALY50_54m_scen4)


#Age 55-59y, male

DALY55_59m_scen4 <- pTT4_55_59m * pop_55_59m * D * dw
mean_median_ci(DALY55_59m_scen4)


#Age 60-64y, male

DALY60_64m_scen4 <- pTT4_60_64m * pop_60_64m * D * dw
mean_median_ci(DALY60_64m_scen4)


#Age 65-69y, male

DALY65_69m_scen4 <- pTT4_65_69m * pop_65_69m * D * dw
mean_median_ci(DALY65_69m_scen4)


#Age 70-74y, male

DALY70_74m_scen4 <- pTT4_70_74m * pop_70_74m * D * dw
mean_median_ci(DALY70_74m_scen4)


#Age 75-79y, male

DALY75_79m_scen4 <- pTT4_75_79m * pop_75_79m * D * dw
mean_median_ci(DALY75_79m_scen4)


#Age 80-84y, male

DALY80_84m_scen4 <- pTT4_75_79m * pop_80_84m * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY80_84m_scen4)


#Age 85+y, male

DALY85m_scen4 <- pTT4_75_79m * pop_85m * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY85m_scen4)


##### Females #####

#Age 15-19, female

DALY15_19w_scen4 <- pTT4_15_19w * pop_15_19w * D * dw
mean_median_ci(DALY15_19w_scen4)


#Age 20-24y, female

DALY20_24w_scen4 <- pTT4_20_24w * pop_20_24w * D * dw
mean_median_ci(DALY20_24w_scen4)


#Age 25-29y, female

DALY25_29w_scen4 <- pTT4_25_29w * pop_25_29w * D * dw
mean_median_ci(DALY25_29w_scen4)


#Age 30-34y, female

DALY30_34w_scen4 <- pTT4_30_34w * pop_30_34w * D * dw
mean_median_ci(DALY30_34w_scen4)


#Age 35-39y, female

DALY35_39w_scen4 <- pTT4_35_39w * pop_35_39w * D * dw
mean_median_ci(DALY35_39w_scen4)



#Age 40-44y, female

DALY40_44w_scen4 <- pTT4_40_44w * pop_40_44w * D * dw
mean_median_ci(DALY40_44w_scen4)


#Age 45-49y, female

DALY45_49w_scen4 <- pTT4_45_49w * pop_45_49w * D * dw
mean_median_ci(DALY45_49w_scen4)


#Age 50-54y, female

DALY50_54w_scen4 <- pTT4_50_54w * pop_50_54w * D * dw
mean_median_ci(DALY50_54w_scen4)


#Age 55-59y, female

DALY55_59w_scen4 <- pTT4_55_59w * pop_55_59w * D * dw
mean_median_ci(DALY55_59w_scen4)


#Age 60-64y, female

DALY60_64w_scen4 <- pTT4_60_64w * pop_60_64w * D * dw
mean_median_ci(DALY60_64w_scen4)


#Age 65-69y, female

DALY65_69w_scen4 <- pTT4_65_69w * pop_65_69w * D * dw
mean_median_ci(DALY65_69w_scen4)


#Age 70-74y, female

DALY70_74w_scen4 <- pTT4_70_74w * pop_70_74w * D * dw
mean_median_ci(DALY70_74w_scen4)


#Age 75-79y, female

DALY75_79w_scen4 <- pTT4_75_79w * pop_75_79w * D * dw
mean_median_ci(DALY75_79w_scen4)


#Age 80-84y, female

DALY80_84w_scen4 <- pTT4_75_79w * pop_80_84w * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY80_84w_scen4)


#Age 85+y, female

DALY85w_scen4 <- pTT4_75_79w * pop_85w * D * dw #Extrapolate probability of hypothyroidism for 75 yearolds to > 75
mean_median_ci(DALY85w_scen4)


##### Total DALYs scen 4 #####

DALYtotal_scen4 <- DALY15_19m_scen4 + DALY20_24m_scen4 + DALY25_29m_scen4 + DALY30_34m_scen4 + DALY35_39m_scen4 +
  DALY40_44m_scen4 + DALY45_49m_scen4 + DALY50_54m_scen4 +  DALY55_59m_scen4 + DALY60_64m_scen4 + DALY65_69m_scen4 +
  DALY70_74m_scen4 + DALY75_79m_scen4 + DALY80_84m_scen4 + DALY85m_scen4 +
  DALY15_19w_scen4 + DALY20_24w_scen4 + DALY25_29w_scen4 + DALY30_34w_scen4 + DALY35_39w_scen4 + DALY40_44w_scen4 +
  DALY45_49w_scen4 + DALY50_54w_scen4 + DALY55_59w_scen4 + DALY60_64w_scen4 + DALY65_69w_scen4 + DALY70_74w_scen4 +
  DALY75_79w_scen4 + DALY80_84w_scen4 + DALY85w_scen4


mean_median_ci(DALYtotal_scen4)


##### Cases scenario 4 #####

cases_scen4 <- pTT4_15_19m * pop_15_19m + pTT4_20_24m * pop_20_24m + pTT4_25_29m * pop_25_29m + pTT4_30_34m * pop_30_34m +
  pTT4_35_39m * pop_35_39m + pTT4_40_44m * pop_40_44m + pTT4_45_49m * pop_45_49m + pTT4_50_54m * pop_50_54m + pTT4_55_59m * pop_55_59m +
  pTT4_60_64m * pop_60_64m + pTT4_65_69m * pop_65_69m + pTT4_70_74m * pop_70_74m + pTT4_75_79m * pop_75_79m + pTT4_75_79m * pop_80_84m +
  pTT4_75_79m * pop_85m +
  pTT4_15_19w * pop_15_19w + pTT4_20_24w * pop_20_24w + pTT4_25_29w * pop_25_29w + pTT4_30_34w * pop_30_34w +
  pTT4_35_39w * pop_35_39w + pTT4_40_44w * pop_40_44w + pTT4_45_49w * pop_45_49w + pTT4_50_54w * pop_50_54w + pTT4_55_59w * pop_55_59w +
  pTT4_60_64w * pop_60_64w + pTT4_65_69w * pop_65_69w + pTT4_70_74w * pop_70_74w + pTT4_75_79w * pop_75_79w + pTT4_75_79w * pop_80_84w +
  pTT4_75_79w * pop_85w


mean_median_ci(cases_scen4)


#### Total DALYs all scenarios #####

tDALY_ref_diox_TT4 <- DALYtotal
mean_median_ci(tDALY_ref_diox_TT4)


tDALY_scen1_diox_TT4 <- DALYtotal_scen1
mean_median_ci(tDALY_scen1_diox_TT4)


tDALY_scen2_diox_TT4 <- DALYtotal_scen2
mean_median_ci(tDALY_scen2_diox_TT4)


tDALY_scen3_diox_TT4 <- DALYtotal_scen3
mean_median_ci(tDALY_scen3_diox_TT4)


tDALY_scen4_diox_TT4 <- DALYtotal_scen4
mean_median_ci(tDALY_scen4_diox_TT4)



##### Delta DALYs #####

dDALY_scen1_diox_TT4 <- tDALY_scen1_diox_TT4 - tDALY_ref_diox_TT4
mean_median_ci(dDALY_scen1_diox_TT4)


dDALY_scen2_diox_TT4 <- tDALY_scen2_diox_TT4 - tDALY_ref_diox_TT4
mean_median_ci(dDALY_scen2_diox_TT4)


dDALY_scen3_diox_TT4 <- tDALY_scen3_diox_TT4 - tDALY_ref_diox_TT4
mean_median_ci(dDALY_scen3_diox_TT4)


dDALY_scen4_diox_TT4 <- tDALY_scen4_diox_TT4 - tDALY_ref_diox_TT4
mean_median_ci(dDALY_scen4_diox_TT4)



##### Delta DALYs per 100,000 #####

dDALY_scen1_diox_TT4_100000 <- dDALY_scen1_diox_TT4/pop15_85_plus*1e+05
mean_median_ci(dDALY_scen1_diox_TT4_100000)


dDALY_scen2_diox_TT4_100000 <- dDALY_scen2_diox_TT4/pop15_85_plus*1e+05
mean_median_ci(dDALY_scen2_diox_TT4_100000)


dDALY_scen3_diox_TT4_100000 <- dDALY_scen3_diox_TT4/pop15_85_plus*1e+05
mean_median_ci(dDALY_scen3_diox_TT4_100000)


dDALY_scen4_diox_TT4_100000 <- dDALY_scen4_diox_TT4/pop15_85_plus*1e+05
mean_median_ci(dDALY_scen4_diox_TT4_100000)



##### Extra number of cases #####

cases_diff_scen1 <- cases_scen1 - cases_ref
mean_median_ci(cases_diff_scen1)


cases_diff_scen2 <- cases_scen2 - cases_ref
mean_median_ci(cases_diff_scen2)


cases_diff_scen3 <- cases_scen3 - cases_ref
mean_median_ci(cases_diff_scen3)


cases_diff_scen4 <- cases_scen4 - cases_ref
mean_median_ci(cases_diff_scen4)




###### Extra cases per 100,000 #####

cases_diff_scen1_100000 <- cases_diff_scen1/pop15_85_plus*1e+05
mean_median_ci(cases_diff_scen1_100000)


cases_diff_scen2_100000 <- cases_diff_scen2/pop15_85_plus*1e+05
mean_median_ci(cases_diff_scen2_100000)


cases_diff_scen3_100000 <- cases_diff_scen3/pop15_85_plus*1e+05
mean_median_ci(cases_diff_scen3_100000)


cases_diff_scen4_100000 <- cases_diff_scen4/pop15_85_plus*1e+05
mean_median_ci(cases_diff_scen4_100000)



##### Save #####


write.csv(tDALY_ref_diox_TT4, "tDALY_ref_diox_TT4.csv")
write.csv(tDALY_scen1_diox_TT4, "tDALY_scen1_diox_TT4.csv")
write.csv(tDALY_scen2_diox_TT4, "tDALY_scen2_diox_TT4.csv")
write.csv(tDALY_scen3_diox_TT4, "tDALY_scen3_diox_TT4.csv")
write.csv(tDALY_scen4_diox_TT4, "tDALY_scen4_diox_TT4.csv")


