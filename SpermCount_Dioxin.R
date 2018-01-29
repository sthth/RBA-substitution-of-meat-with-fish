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


## Population numbers

women15_19 <- 171815 
women20_24 <- 184806
women25_29 <- 170197
women30_34 <- 158343
women35_39 <- 178824
women40_44 <- 194114
women45_49 <- 204582
women15_49 <- 1262681

pop15_85_plus <- 4697068


## IPRA

#Assume that 10% decrease in sperm count cause infertility
#CED for 10% decrease in sperm count due to prenatal exposure


# critical effect dose - animals

#CED 20% Hill
CEDM_animal <- 743.09
CEDL_animal <- 108.3
CEDU_animal <- 27710

set.seed(1)
CED_animal <- rpert(nunc, CEDL_animal, CEDM_animal, CEDU_animal)
hist(CED_animal)

## critical effect dose - humans
CED_human <- CED_animal / 3082*1000
hist(CED_human)

## extrapolation factor - intraspecies
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

Pbaby15_19 <- 0.0023 #probability of giving birth in 2015
Pbaby20_24 <- 0.029
Pbaby25_29 <- 0.1027
Pbaby30_34 <- 0.13
Pbaby35_39 <- 0.064
Pbaby40_44 <- 0.0139
Pbaby45_49 <- 0.0008


pmale <- 0.51 #probability of pregnant woman giving birth to a boy

#YLL=0

#YLD = Pbaby * Pmale * PPinfetile * D * DW

#Duration=29 (49-20)

D=29

#DW for primary infertility (Salomon 2015)
set.seed(1)
dw <- rpert(nunc, min=0.003, mode=0.008, max=0.015)


##### Reference scenario #####

Dioxin_exp <- read.csv("Dioxin_exposure.csv")

agebreaks <- c(15,20,25,30,35,40,45,50,55,60,65,70,75,80)
agelabels= c("15-19","20-24","25-29","30-34","35-39", "40-44", "45-49", "50-54","55-59", "60-64","65-69", "70-74","75-79")

setDT(Dioxin_exp)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]



##### Fitting exposures #####

#Endpoint only relevant for women in the fertile age (15-49y)

#Age 15-19
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


#Age 20-24
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


#Age 25-29
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


#Age 30-34
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


#Age 35-39
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


#Age 40-44
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


#Age 45-49
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



##### Probability of infertility #####

# Assume exposure is variable, lognormal has the best fit


### Age 15-19 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_15_19 <- apply(ICED, 2, function(x) mean(x < exp)) #probability in each column of ICED>exp (the probability in each column takes account of the variability in that column)



### Age 20-24 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_20_24 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 25-29 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_25_29 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 30-34 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_30_34 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 35-39 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_35_39 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 40-44 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_40_44 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 45-49 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_45_49 <- apply(ICED, 2, function(x) mean(x < exp))



###### DALY calculation  #####


#Age 15-19y

DALY15_19w <- Pbaby15_19 * pmale * pinfert_15_19 * D * dw * women15_19
mean_median_ci(DALY15_19w)



#Age 20-24y

DALY20_24w <- Pbaby20_24 * pmale * pinfert_20_24 * D * dw * women20_24
mean_median_ci(DALY20_24w)



#Age 25-29y

DALY25_29w <- Pbaby25_29 * pmale * pinfert_25_29 * D * dw * women25_29
mean_median_ci(DALY25_29w)



#Age 30-34y

DALY30_34w <- Pbaby30_34 * pmale * pinfert_30_34 * D * dw * women30_34
mean_median_ci(DALY30_34w)



#Age 35-39y

DALY35_39w <- Pbaby35_39 * pmale * pinfert_35_39 * D * dw * women35_39
mean_median_ci(DALY35_39w)



#Age 40-44y

DALY40_44w <- Pbaby40_44 * pmale * pinfert_40_44 * D * dw * women40_44
mean_median_ci(DALY40_44w)



#Age 45-49y

DALY45_49w <- Pbaby45_49 * pmale * pinfert_45_49 * D * dw * women45_49
mean_median_ci(DALY45_49w)



##### Total DALYs ref #####

DALYtotal <- DALY15_19w + DALY20_24w + DALY25_29w + DALY30_34w + DALY35_39w + DALY40_44w + DALY45_49w
mean_median_ci(DALYtotal)


##### Cases reference scenario #####

cases_ref <- Pbaby15_19 * pmale * pinfert_15_19 * women15_19 + Pbaby20_24 * pmale * pinfert_20_24 * women20_24 +
  Pbaby25_29 * pmale * pinfert_25_29 * women25_29 + Pbaby30_34 * pmale * pinfert_30_34 * women30_34 +
  Pbaby35_39 * pmale * pinfert_35_39 * women35_39 + Pbaby40_44 * pmale * pinfert_40_44 * women40_44 +
  Pbaby45_49 * pmale * pinfert_45_49 * women45_49

mean_median_ci(cases_ref)



##### Scenario 1 #####

Dioxin_exp1 <- read.csv("Dioxin_exp_scen1.csv")

setDT(Dioxin_exp1)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Fitting exposures #####

#Endpoint only relevant for women in the fertile age (15-49y)

#Age 15-19
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


#Age 20-24
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


#Age 25-29
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


#Age 30-34
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


#Age 35-39
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


#Age 40-44
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


#Age 45-49
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



##### Probability of infertility #####


# Assume exposure is variable, lognormal has the best fit

### Age 15-19 ###


## Assume exposure is variable

set.seed(1)
exp <- rlnorm(nvar, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_15_19 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 20-24 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_20_24 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 25-29 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_25_29 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 30-34 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_30_34 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 35-39 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_35_39 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 40-44 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_40_44 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 45-49 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_45_49 <- apply(ICED, 2, function(x) mean(x < exp))



##### DALY calculation #####

#Age 15-19y

DALY15_19w_scen1 <- Pbaby15_19 * pmale * pinfert_15_19 * D * dw * women15_19
mean_median_ci(DALY15_19w_scen1)


#Age 20-24y

DALY20_24w_scen1 <- Pbaby20_24 * pmale * pinfert_20_24 * D * dw * women20_24
mean_median_ci(DALY20_24w_scen1)


#Age 25-29y

DALY25_29w_scen1 <- Pbaby25_29 * pmale * pinfert_25_29 * D * dw * women25_29
mean_median_ci(DALY25_29w_scen1)


#Age 30-34y

DALY30_34w_scen1 <- Pbaby30_34 * pmale * pinfert_30_34 * D * dw * women30_34
mean_median_ci(DALY30_34w_scen1)


#Age 35-39y

DALY35_39w_scen1 <- Pbaby35_39 * pmale * pinfert_35_39 * D * dw * women35_39
mean_median_ci(DALY35_39w_scen1)


#Age 40-44y

DALY40_44w_scen1 <- Pbaby40_44 * pmale * pinfert_40_44 * D * dw * women40_44
mean_median_ci(DALY40_44w_scen1)


#Age 45-49y

DALY45_49w_scen1 <- Pbaby45_49 * pmale * pinfert_45_49 * D * dw * women45_49
mean_median_ci(DALY45_49w_scen1)


#### Total DALYs scenario 1 #####

DALYtotal_scen1 <- DALY15_19w_scen1 + DALY20_24w_scen1 + DALY25_29w_scen1 + DALY30_34w_scen1 + DALY35_39w_scen1 +
  DALY40_44w_scen1 + DALY45_49w_scen1

mean_median_ci(DALYtotal_scen1)



##### Cases scenario 1 #####

cases_scen1 <- Pbaby15_19 * pmale * pinfert_15_19 * women15_19 + Pbaby20_24 * pmale * pinfert_20_24 * women20_24 +
  Pbaby25_29 * pmale * pinfert_25_29 * women25_29 + Pbaby30_34 * pmale * pinfert_30_34 * women30_34 +
  Pbaby35_39 * pmale * pinfert_35_39 * women35_39 + Pbaby40_44 * pmale * pinfert_40_44 * women40_44 +
  Pbaby45_49 * pmale * pinfert_45_49 * women45_49

mean_median_ci(cases_scen1)


##### Scenario 2 #####

Dioxin_exp2 <- read.csv("Dioxin_exp_scen2.csv")


setDT(Dioxin_exp2)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]



##### Fitting exposures #####

#Endpoint only relevant for women in the fertile age (15-49y)

#Age 15-19
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


#Age 20-24
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


#Age 25-29
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


#Age 30-34
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


#Age 35-39
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


#Age 40-44
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


#Age 45-49
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




##### Probability of infertility #####

# Assume exposure is variable, lognormal has the best fit



### Age 15-19 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_15_19 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 20-24 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_20_24 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 25-29 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_25_29 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 30-34 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_30_34 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 35-39 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_35_39 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 40-44 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_40_44 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 45-49 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_45_49 <- apply(ICED, 2, function(x) mean(x < exp))



##### DALY calculation #####

#Age 15-19y

DALY15_19w_scen2 <- Pbaby15_19 * pmale * pinfert_15_19 * D * dw * women15_19
mean_median_ci(DALY15_19w_scen2)


#Age 20-24y

DALY20_24w_scen2 <- Pbaby20_24 * pmale * pinfert_20_24 * D * dw * women20_24
mean_median_ci(DALY20_24w_scen2)


#Age 25-29y

DALY25_29w_scen2 <- Pbaby25_29 * pmale * pinfert_25_29 * D * dw * women25_29
mean_median_ci(DALY25_29w_scen2)


#Age 30-34y

DALY30_34w_scen2 <- Pbaby30_34 * pmale * pinfert_30_34 * D * dw * women30_34
mean_median_ci(DALY30_34w_scen2)


#Age 35-39y

DALY35_39w_scen2 <- Pbaby35_39 * pmale * pinfert_35_39 * D * dw * women35_39
mean_median_ci(DALY35_39w_scen2)


#Age 40-44y

DALY40_44w_scen2 <- Pbaby40_44 * pmale * pinfert_40_44 * D * dw * women40_44
mean_median_ci(DALY40_44w_scen2)


#Age 45-49y

DALY45_49w_scen2 <- Pbaby45_49 * pmale * pinfert_45_49 * D * dw * women45_49
mean_median_ci(DALY45_49w_scen2)



##### Total DALYs scenario 2 #####

DALYtotal_scen2 <- DALY15_19w_scen2 + DALY20_24w_scen2 + DALY25_29w_scen2 + DALY30_34w_scen2 + DALY35_39w_scen2 +
  DALY40_44w_scen2 + DALY45_49w_scen2

mean_median_ci(DALYtotal_scen2)




##### Cases scenario 2 #####

cases_scen2 <- Pbaby15_19 * pmale * pinfert_15_19 * women15_19 + Pbaby20_24 * pmale * pinfert_20_24 * women20_24 +
  Pbaby25_29 * pmale * pinfert_25_29 * women25_29 + Pbaby30_34 * pmale * pinfert_30_34 * women30_34 +
  Pbaby35_39 * pmale * pinfert_35_39 * women35_39 + Pbaby40_44 * pmale * pinfert_40_44 * women40_44 +
  Pbaby45_49 * pmale * pinfert_45_49 * women45_49

mean_median_ci(cases_scen2)



##### Scenario 3 #####

Dioxin_exp3 <- read.csv("Dioxin_exp_scen3.csv")

setDT(Dioxin_exp3)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]



##### Fitting exposures #####

#Endpoint only relevant for women in the fertile age (15-49y)

#Age 15-19
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


#Age 20-24
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


#Age 25-29
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


#Age 30-34
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


#Age 35-39
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


#Age 40-44
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


#Age 45-49
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



##### Probability of infertility #####

# Assume exposure is variable, lognormal has the best fit


### Age 15-19 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_15_19 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 20-24 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_20_24 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 25-29 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_25_29 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 30-34 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_30_34 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 35-39 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_35_39 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 40-44 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_40_44 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 45-49 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_45_49 <- apply(ICED, 2, function(x) mean(x < exp))



##### DALY calculation #####

#Age 15-19y

DALY15_19w_scen3 <- Pbaby15_19 * pmale * pinfert_15_19 * D * dw * women15_19
mean_median_ci(DALY15_19w_scen3)


#Age 20-24y

DALY20_24w_scen3 <- Pbaby20_24 * pmale * pinfert_20_24 * D * dw * women20_24
mean_median_ci(DALY20_24w_scen3)


#Age 25-29y

DALY25_29w_scen3 <- Pbaby25_29 * pmale * pinfert_25_29 * D * dw * women25_29
mean_median_ci(DALY25_29w_scen3)


#Age 30-34y

DALY30_34w_scen3 <- Pbaby30_34 * pmale * pinfert_30_34 * D * dw * women30_34
mean_median_ci(DALY30_34w_scen3)


#Age 35-39y

DALY35_39w_scen3 <- Pbaby35_39 * pmale * pinfert_35_39 * D * dw * women35_39
mean_median_ci(DALY35_39w_scen3)


#Age 40-44y

DALY40_44w_scen3 <- Pbaby40_44 * pmale * pinfert_40_44 * D * dw * women40_44
mean_median_ci(DALY40_44w_scen3)


#Age 45-49y

DALY45_49w_scen3 <- Pbaby45_49 * pmale * pinfert_45_49 * D * dw * women45_49
mean_median_ci(DALY45_49w_scen3)



##### Total DALYs scenario 3 #####

DALYtotal_scen3 <- DALY15_19w_scen3 + DALY20_24w_scen3 + DALY25_29w_scen3 + DALY30_34w_scen3 + DALY35_39w_scen3 +
  DALY40_44w_scen3 + DALY45_49w_scen3

mean_median_ci(DALYtotal_scen3)



##### Cases scenario 3 #####

cases_scen3 <- Pbaby15_19 * pmale * pinfert_15_19 * women15_19 + Pbaby20_24 * pmale * pinfert_20_24 * women20_24 +
  Pbaby25_29 * pmale * pinfert_25_29 * women25_29 + Pbaby30_34 * pmale * pinfert_30_34 * women30_34 +
  Pbaby35_39 * pmale * pinfert_35_39 * women35_39 + Pbaby40_44 * pmale * pinfert_40_44 * women40_44 +
  Pbaby45_49 * pmale * pinfert_45_49 * women45_49

mean_median_ci(cases_scen3)



##### Scenario 4 #####

Dioxin_exp4 <- read.csv("Dioxin_exp_scen4.csv")


setDT(Dioxin_exp4)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]



##### Fitting exposures #####

#Endpoint only relevant for women in the fertile age (15-49y)

#Age 15-19
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


#Age 20-24
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


#Age 25-29
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


#Age 30-34
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


#Age 35-39
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


#Age 40-44
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


#Age 45-49
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



##### Probability of infertility #####

# Assume exposure is variable, lognormal has the best fit


### Age 15-19 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_15_19w$estimate[1], fit1_15_19w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_15_19 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 20-24 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_20_24w$estimate[1], fit1_20_24w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_20_24 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 25-29 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_25_29w$estimate[1], fit1_25_29w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_25_29 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 30-34 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_30_34w$estimate[1], fit1_30_34w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_30_34 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 35-39 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_35_39w$estimate[1], fit1_35_39w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_35_39 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 40-44 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_40_44w$estimate[1], fit1_40_44w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_40_44 <- apply(ICED, 2, function(x) mean(x < exp))



### Age 45-49 ###

set.seed(1)
exp <- rlnorm(nvar, fit1_45_49w$estimate[1], fit1_45_49w$estimate[2]) #potential truncation at 6 pg/kg bw/day - not necessary?

pinfert_45_49 <- apply(ICED, 2, function(x) mean(x < exp))



##### DALY calculation #####

#Age 15-19y

DALY15_19w_scen4 <- Pbaby15_19 * pmale * pinfert_15_19 * D * dw * women15_19
mean_median_ci(DALY15_19w_scen4)


#Age 20-24y

DALY20_24w_scen4 <- Pbaby20_24 * pmale * pinfert_20_24 * D * dw * women20_24
mean_median_ci(DALY20_24w_scen4)


#Age 25-29y

DALY25_29w_scen4 <- Pbaby25_29 * pmale * pinfert_25_29 * D * dw * women25_29
mean_median_ci(DALY25_29w_scen4)


#Age 30-34y

DALY30_34w_scen4 <- Pbaby30_34 * pmale * pinfert_30_34 * D * dw * women30_34
mean_median_ci(DALY30_34w_scen4)


#Age 35-39y

DALY35_39w_scen4 <- Pbaby35_39 * pmale * pinfert_35_39 * D * dw * women35_39
mean_median_ci(DALY35_39w_scen4)


#Age 40-44y

DALY40_44w_scen4 <- Pbaby40_44 * pmale * pinfert_40_44 * D * dw * women40_44
mean_median_ci(DALY40_44w_scen4)


#Age 45-49y

DALY45_49w_scen4 <- Pbaby45_49 * pmale * pinfert_45_49 * D * dw * women45_49
mean_median_ci(DALY45_49w_scen4)


##### Total DALYs scenario 4 #####

DALYtotal_scen4 <- DALY15_19w_scen4 + DALY20_24w_scen4 + DALY25_29w_scen4 + DALY30_34w_scen4 + DALY35_39w_scen4 +
  DALY40_44w_scen4 + DALY45_49w_scen4

mean_median_ci(DALYtotal_scen4)



##### Cases scenario 4 #####

cases_scen4 <- Pbaby15_19 * pmale * pinfert_15_19 * women15_19 + Pbaby20_24 * pmale * pinfert_20_24 * women20_24 +
  Pbaby25_29 * pmale * pinfert_25_29 * women25_29 + Pbaby30_34 * pmale * pinfert_30_34 * women30_34 +
  Pbaby35_39 * pmale * pinfert_35_39 * women35_39 + Pbaby40_44 * pmale * pinfert_40_44 * women40_44 +
  Pbaby45_49 * pmale * pinfert_45_49 * women45_49

mean_median_ci(cases_scen4)



##### Total DALYs all scenarios #####


tDALY_ref_diox_fert <- DALYtotal
mean_median_ci(tDALY_ref_diox_fert)



tDALY_scen1_diox_fert <- DALYtotal_scen1
mean_median_ci(tDALY_scen1_diox_fert)



tDALY_scen2_diox_fert <- DALYtotal_scen2
mean_median_ci(tDALY_scen2_diox_fert)



tDALY_scen3_diox_fert <- DALYtotal_scen3
mean_median_ci(tDALY_scen3_diox_fert)



tDALY_scen4_diox_fert <- DALYtotal_scen4
mean_median_ci(tDALY_scen4_diox_fert)



##### Delta DALYs #####

dDALY_scen1_diox_fert <- tDALY_scen1_diox_fert - tDALY_ref_diox_fert
mean_median_ci(dDALY_scen1_diox_fert)


dDALY_scen2_diox_fert <- tDALY_scen2_diox_fert - tDALY_ref_diox_fert
mean_median_ci(dDALY_scen2_diox_fert)


dDALY_scen3_diox_fert <- tDALY_scen3_diox_fert - tDALY_ref_diox_fert
mean_median_ci(dDALY_scen3_diox_fert)


dDALY_scen4_diox_fert <- tDALY_scen4_diox_fert - tDALY_ref_diox_fert
mean_median_ci(dDALY_scen4_diox_fert)



##### Delta DALYs per 100,000 #####

dDALY_scen1_diox_fert_100000 <- dDALY_scen1_diox_fert/pop15_85_plus*1e+05
mean_median_ci(dDALY_scen1_diox_fert_100000)


dDALY_scen2_diox_fert_100000 <- dDALY_scen2_diox_fert/pop15_85_plus*1e+05
mean_median_ci(dDALY_scen2_diox_fert_100000)


dDALY_scen3_diox_fert_100000 <- dDALY_scen3_diox_fert/pop15_85_plus*1e+05
mean_median_ci(dDALY_scen3_diox_fert_100000)


dDALY_scen4_diox_fert_100000 <- dDALY_scen4_diox_fert/pop15_85_plus*1e+05
mean_median_ci(dDALY_scen4_diox_fert_100000)



##### Extra number of cases #####

cases_diff_scen1 <- cases_scen1 - cases_ref
mean_median_ci(cases_diff_scen1)


cases_diff_scen2 <- cases_scen2 - cases_ref
mean_median_ci(cases_diff_scen2)


cases_diff_scen3 <- cases_scen3 - cases_ref
mean_median_ci(cases_diff_scen3)


cases_diff_scen4 <- cases_scen4 - cases_ref
mean_median_ci(cases_diff_scen4)



##### Extra cases per 100,000 #####

cases_diff_scen1_100000 <- cases_diff_scen1/pop15_85_plus*1e+05
mean_median_ci(cases_diff_scen1_100000)


cases_diff_scen2_100000 <- cases_diff_scen2/pop15_85_plus*1e+05
mean_median_ci(cases_diff_scen2_100000)


cases_diff_scen3_100000 <- cases_diff_scen3/pop15_85_plus*1e+05
mean_median_ci(cases_diff_scen3_100000)


cases_diff_scen4_100000 <- cases_diff_scen4/pop15_85_plus*1e+05
mean_median_ci(cases_diff_scen4_100000)




##### Save #####


write.csv(tDALY_ref_diox_fert, "tDALY_ref_diox_fert.csv")
write.csv(tDALY_scen1_diox_fert, "tDALY_scen1_diox_fert.csv")
write.csv(tDALY_scen2_diox_fert, "tDALY_scen2_diox_fert.csv")
write.csv(tDALY_scen3_diox_fert, "tDALY_scen3_diox_fert.csv")
write.csv(tDALY_scen4_diox_fert, "tDALY_scen4_diox_fert.csv")