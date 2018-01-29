memory.limit(56000)


#Packages
library(data.table)
library(fitdistrplus)
library(goftest)
library(mc2d)

## helpers
mean_median_ci <-
  function(x) {
    c(mean = mean(x),
      quantile(x, probs = c(0.025, 0.975)))
  }



## settings
nvar <- 1e5
nunc <- 1e3
ndvar(nvar)
ndunc(nunc)


#Population numbers
women15_19 <- 171815
women20_24 <- 184806
women25_29 <- 170197
women30_34 <- 158343
women35_39 <- 178824
women40_44 <- 194114
women45_49 <- 204582
women15_49 <- 1262681


#Probability of giving birth
Pbaby15_19 <- 0.0023
Pbaby20_24 <- 0.029
Pbaby25_29 <- 0.1027
Pbaby30_34 <- 0.13
Pbaby35_39 <- 0.064
Pbaby40_44 <- 0.0139
Pbaby45_49 <- 0.0008


LE <- 80.6

births2015 <- 58205 #number of babies born in 2015 (statistikbanken)

ndvar(100000)
ndunc(1000)

#Danish life expectancy at birth
LE <- 80.6


#Disability weights intellectual disabilities
wIQ85 <- 0 #IQ > 85
wIQ70.85sim <- rpert(nunc, min=0.005, mode=0.011, max=0.020) #IQ 70-85
wIQ50.69sim <- rpert(nunc, min=0.026, mode=0.043, max=0.064) #IQ 50-69
wIQ35.49sim <- rpert(nunc, min=0.066, mode=0.100, max=0.142) #IQ 35-49
wIQ20.34sim <- rpert(nunc, min=0.107, mode=0.160, max=0.226) #IQ 20-34
wIQ20sim <- rpert(nunc, min=0.133, mode=0.200, max=0.283) #IQ < 20


#Dose-response
set.seed(1)
r <- rpert(nunc, min = -19.5, mode = -8.5, max = -1.5) #Dose-response Zeilmaker et al., 2013
IQnorm <- rnorm(nvar, mean = 100, sd = 15) #Definition of IQ



##### Reference scenario #####

FishDaily <- read.csv("FishDaily.csv")

#Water loss has NOT been taken into account in FoodDaily/FishDaily dataset so amounts are prepared amounts.
#Water loss HAS been taken into account in the substitution scenarios


#Convert fish consumption amounts into raw weights - assume that all fish is prepared and have a water loss of 20%
wloss.fish <- 1/0.8
FishDaily[,8:23] <- FishDaily[,8:23]*wloss.fish
FishDaily$Totalfish <- rowSums(FishDaily[,8:23]) #Test if it gives the same as total.fish * 1.25
names(FishDaily)[24] <- 'Totalfish.raw' #Rename to specify that it's raw


MeHg_table <- read.csv2("Mercury_fish.csv") #Mercury concentrations in fish and shellfish (ug/g)


#EFSA: We can assume that 100% of Hg in fish is MeHg - 80% for crustaceans and molluscs. Data for Hg in shrimps, crabs and mussels
#are only Hg - multiply mean Hg concentration in these with 0.8 
MeHg_table[,4:6] <- MeHg_table[,4:6]*0.8 #Convert 80% of Hg in shrimp, crab and mussel to MeHg

MeHg_table <-  as.data.frame(t(MeHg_table))

#Calculate MeHg exposure from fish (ug/day)
FishMeHg <- t(FishDaily[,c(8:23)]) #Create data set only with intakes of the different fish species
FishMeHg <- t(FishMeHg * MeHg_table[1:16,]) ##Multiply intakes of fish species with MeHg conc of the different species to obtain dataset with MeHg exposures from the different species


FishDaily$MeHg <- rowSums(FishMeHg) #Add sum of MeHg exposures for each individual to the FishDaily dataset

#MeHg/kg bw/week 
FishDaily$MeHg_bw <- FishDaily$MeHg / FishDaily$weight


#Exclude males, and women above 49y of age (due to data on female fertility) - women below 15 have already been excluded

FishDaily <- subset(FishDaily, age < 50 & sex == '2')


agebreaks <- c(15,20,25,30,35,40,45,50)
agelabels= c("15-19","20-24","25-29","30-34","35-39","40-44","45-49")

setDT(FishDaily)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


#################################################### Fitting exposures ###########################################################

#Divide data into agegroups and fit exposure to distribution

#Age 15-19
MeHg15_19 <- subset(FishDaily$MeHg_bw, FishDaily$agegroups=="15-19")

#Probability of fish intake
probMeHg15_19 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(MeHg15_19==0)/length(MeHg15_19) #probability of zero MeHg exp
pr.yes <- sum(MeHg15_19!=0)/length(MeHg15_19) #probability of MeHg exp
probMeHg15_19[1,] <- c(pr.no,pr.yes)


MeHg15_19_pos <- MeHg15_19[which(MeHg15_19!=0)] #Only intakes that are greater than zero

fit1_15_19 <- fitdist(MeHg15_19_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- MeHg15_19_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_15_19$estimate["meanlog"], sdlog=fit1_15_19$estimate["sdlog"]), col="red", lty=1)
gofstat(fit1_15_19, fitnames="lnorm") 

fit2_15_19 <- fitdist(MeHg15_19_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- MeHg15_19_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_15_19$estimate[1],fit2_15_19$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19, fitnames="gamma")

library(goftest)
cvm.test(t, plnorm, fit1_15_19$estimate[1], fit1_15_19$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19$estimate[1], fit1_15_19$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19$estimate[1], fit2_15_19$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19$estimate[1], fit2_15_19$estimate[2])  #Anderson-Darling test


#Age 20-24
MeHg20_24 <- subset(FishDaily$MeHg_bw, FishDaily$agegroups=="20-24")

probMeHg20_24 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(MeHg20_24==0)/length(MeHg20_24) #probability of zero MeHg exp
pr.yes <- sum(MeHg20_24!=0)/length(MeHg20_24) #probability of MeHg exp
probMeHg20_24[1,] <- c(pr.no,pr.yes)

MeHg20_24_pos <- MeHg20_24[which(MeHg20_24!=0)] #Only intakes that are greater than zero

fit1_20_24 <- fitdist(MeHg20_24_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- MeHg20_24_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_20_24$estimate["meanlog"], sdlog=fit1_20_24$estimate["sdlog"]), col="red", lty=1)
gofstat(fit1_20_24, fitnames="lnorm") 

fit2_20_24 <- fitdist(MeHg20_24_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- MeHg20_24_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_20_24$estimate[1],fit2_20_24$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24, fitnames="gamma")

cvm.test(t, plnorm, fit1_20_24$estimate[1], fit1_20_24$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24$estimate[1], fit1_20_24$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24$estimate[1], fit2_20_24$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24$estimate[1], fit2_20_24$estimate[2])  #Anderson-Darling test



#Age 25-29
MeHg25_29 <- subset(FishDaily$MeHg_bw, FishDaily$agegroups=="25-29")

probMeHg25_29 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(MeHg25_29==0)/length(MeHg25_29) #probability of zero MeHg exp
pr.yes <- sum(MeHg25_29!=0)/length(MeHg25_29) #probability of MeHg exp
probMeHg25_29[1,] <- c(pr.no,pr.yes)

MeHg25_29_pos <- MeHg25_29[which(MeHg25_29!=0)] #Only intakes that are greater than zero

fit1_25_29 <- fitdist(MeHg25_29_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- MeHg25_29_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_25_29$estimate["meanlog"], sdlog=fit1_25_29$estimate["sdlog"]), col="red", lty=1)
gofstat(fit1_25_29, fitnames="lnorm")

fit2_25_29 <- fitdist(MeHg25_29_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- MeHg25_29_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_25_29$estimate[1],fit2_25_29$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29, fitnames="gamma")

cvm.test(t, plnorm, fit1_25_29$estimate[1], fit1_25_29$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29$estimate[1], fit1_25_29$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29$estimate[1], fit2_25_29$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29$estimate[1], fit2_25_29$estimate[2])  #Anderson-Darling test



#Age 30-34
MeHg30_34 <- subset(FishDaily$MeHg_bw, FishDaily$agegroups=="30-34")

probMeHg30_34 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(MeHg30_34==0)/length(MeHg30_34) #probability of zero MeHg exp
pr.yes <- sum(MeHg30_34!=0)/length(MeHg30_34) #probability of MeHg exp
probMeHg30_34[1,] <- c(pr.no,pr.yes)

MeHg30_34_pos <- MeHg30_34[which(MeHg30_34!=0)] #Only intakes that are greater than zero

fit1_30_34 <- fitdist(MeHg30_34_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- MeHg30_34_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_30_34$estimate["meanlog"], sdlog=fit1_30_34$estimate["sdlog"]), col="red", lty=1)
gofstat(fit1_30_34, fitnames="lnorm") 

fit2_30_34 <- fitdist(MeHg30_34_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- MeHg30_34_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_30_34$estimate[1],fit2_30_34$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34, fitnames="gamma")

cvm.test(t, plnorm, fit1_30_34$estimate[1], fit1_30_34$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34$estimate[1], fit1_30_34$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34$estimate[1], fit2_30_34$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34$estimate[1], fit2_30_34$estimate[2])  #Anderson-Darling test



#Age 35-39
MeHg35_39 <- subset(FishDaily$MeHg_bw, FishDaily$agegroups=="35-39")

probMeHg35_39 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(MeHg35_39==0)/length(MeHg35_39) #probability of zero MeHg exp
pr.yes <- sum(MeHg35_39!=0)/length(MeHg35_39) #probability of MeHg exp
probMeHg35_39[1,] <- c(pr.no,pr.yes)

MeHg35_39_pos <- MeHg35_39[which(MeHg35_39!=0)] #Only intakes that are greater than zero

fit1_35_39 <- fitdist(MeHg35_39_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- MeHg35_39_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_35_39$estimate["meanlog"], sdlog=fit1_35_39$estimate["sdlog"]), col="red", lty=1)
gofstat(fit1_35_39, fitnames="lnorm") 

fit2_35_39 <- fitdist(MeHg35_39_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- MeHg35_39_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_35_39$estimate[1],fit2_35_39$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39, fitnames="gamma")

cvm.test(t, plnorm, fit1_35_39$estimate[1], fit1_35_39$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39$estimate[1], fit1_35_39$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39$estimate[1], fit2_35_39$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39$estimate[1], fit2_35_39$estimate[2])  #Anderson-Darling test



#Age 40-44
MeHg40_44 <- subset(FishDaily$MeHg_bw, FishDaily$agegroups=="40-44")

probMeHg40_44 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(MeHg40_44==0)/length(MeHg40_44) #probability of zero MeHg exp
pr.yes <- sum(MeHg40_44!=0)/length(MeHg40_44) #probability of MeHg exp
probMeHg40_44[1,] <- c(pr.no,pr.yes)

MeHg40_44_pos <- MeHg40_44[which(MeHg40_44!=0)] #Only intakes that are greater than zero

fit1_40_44 <- fitdist(MeHg40_44_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- MeHg40_44_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_40_44$estimate["meanlog"], sdlog=fit1_40_44$estimate["sdlog"]), col="red", lty=1)
gofstat(fit1_40_44, fitnames="lnorm") 

fit2_40_44 <- fitdist(MeHg40_44_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- MeHg40_44_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_40_44$estimate[1],fit2_40_44$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44, fitnames="gamma")

cvm.test(t, plnorm, fit1_40_44$estimate[1], fit1_40_44$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44$estimate[1], fit1_40_44$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44$estimate[1], fit2_40_44$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44$estimate[1], fit2_40_44$estimate[2])  #Anderson-Darling test



#Age 45-49
MeHg45_49 <- subset(FishDaily$MeHg_bw, FishDaily$agegroups=="45-49")

probMeHg45_49 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(MeHg45_49==0)/length(MeHg45_49) #probability of zero MeHg exp
pr.yes <- sum(MeHg45_49!=0)/length(MeHg45_49) #probability of MeHg exp
probMeHg45_49[1,] <- c(pr.no,pr.yes)

MeHg45_49_pos <- MeHg45_49[which(MeHg45_49!=0)] #Only intakes that are greater than zero

fit1_45_49 <- fitdist(MeHg45_49_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- MeHg45_49_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_45_49$estimate["meanlog"], sdlog=fit1_45_49$estimate["sdlog"]), col="red", lty=1)
gofstat(fit1_45_49, fitnames="lnorm") 

fit2_45_49 <- fitdist(MeHg45_49_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- MeHg45_49_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_45_49$estimate[1],fit2_45_49$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49, fitnames="gamma")

cvm.test(t, plnorm, fit1_45_49$estimate[1], fit1_45_49$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49$estimate[1], fit1_45_49$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49$estimate[1], fit2_45_49$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49$estimate[1], fit2_45_49$estimate[2])  #Anderson-Darling test



##### DALY calculation #####

# Lognormal and Gamma almost equally good fits, but gamma dist does not have as extreme values as lognormal 
# thus we chose to use gamma dist


ifexp <- sample(c(0,1), nvar, prob = probMeHg15_19, replace = T)
exp <- rgamma(nvar, fit2_15_19$estimate[1], fit2_15_19$estimate[2])


str(ifexp)
str(exp)
str(IQnorm)
str(r)


str(t(r)) #1 row, 1000 columns - columns: uncertainty

## IQ change: IQdiff = ifexp * exp * r
## .. each col = uncertainty simulation
## .. each row = variability simulation
## .. IQchange[var, unc]

IQdiff15_19 <- apply(t(r), 2, function(x) x * ifexp * exp)
str(IQdiff15_19)
IQdiff15_19[c(1:10), c(1:10)]

IQnew <- apply(t(r), 2, function(x) x * ifexp * exp + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal15_19 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline15_19 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild15_19 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate15_19 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere15_19 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound15_19 <- apply(IQnew, 2, function(x) mean(x < 20))

sum(mean(IQnormal15_19)+mean(IQborderline15_19)+mean(IQmild15_19)+mean(IQmoderate15_19)+mean(IQsevere15_19)+mean(IQprofound15_19))


DALY15_19mc <- Pbaby15_19*(IQnormal15_19*wIQ85 + IQborderline15_19*wIQ70.85sim + IQmild15_19*wIQ50.69sim +
                             IQmoderate15_19*wIQ35.49sim + IQsevere15_19*wIQ20.34sim +
                             IQprofound15_19*wIQ20sim) * LE * women15_19 #per total women in DK in this age group

mean_median_ci(DALY15_19mc)


### Age 20-24 ###

ifexp <- sample(c(0,1), nvar, prob = probMeHg20_24, replace = T)
exp <- rgamma(nvar, fit2_20_24$estimate[1], fit2_20_24$estimate[2])

IQdiff20_24 <- apply(t(r), 2, function(x) x * ifexp * exp)

IQnew <- apply(IQdiff20_24, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal20_24 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline20_24 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild20_24 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate20_24 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere20_24 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound20_24 <- apply(IQnew, 2, function(x) mean(x < 20))

set.seed(1)
DALY20_24mc <- Pbaby20_24*(IQnormal20_24*wIQ85 + IQborderline20_24*wIQ70.85sim + IQmild20_24*wIQ50.69sim +
                             IQmoderate20_24*wIQ35.49sim + IQsevere20_24*wIQ20.34sim +
                             IQprofound20_24*wIQ20sim) * LE * women20_24 #per total women in DK in this age group

mean_median_ci(DALY20_24mc)


### Age 25-29 ###

set.seed(1)
ifexp <- sample(c(0,1), nvar, prob = probMeHg25_29, replace = T)
exp <- rgamma(nvar, fit2_25_29$estimate[1], fit2_25_29$estimate[2])

IQdiff25_29 <- apply(t(r), 2, function(x) x * ifexp * exp)

IQnew <- apply(IQdiff25_29, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal25_29 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline25_29 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild25_29 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate25_29 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere25_29 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound25_29 <- apply(IQnew, 2, function(x) mean(x < 20))

set.seed(1)
DALY25_29mc <- Pbaby25_29*(IQnormal25_29*wIQ85 + IQborderline25_29*wIQ70.85sim + IQmild25_29*wIQ50.69sim +
                             IQmoderate25_29*wIQ35.49sim + IQsevere25_29*wIQ20.34sim +
                             IQprofound25_29*wIQ20sim) * LE * women25_29 #per total women in DK in this age group

mean_median_ci(DALY25_29mc)

### Age 30-34 ###

set.seed(1)
ifexp <- sample(c(0,1), nvar, prob = probMeHg30_34, replace = T)
exp <- rgamma(nvar, fit2_30_34$estimate[1], fit2_30_34$estimate[2])

IQdiff30_34 <- apply(t(r), 2, function(x) x * ifexp * exp)

IQnew <- apply(IQdiff30_34, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal30_34 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline30_34 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild30_34 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate30_34 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere30_34 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound30_34 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY30_34mc <- Pbaby30_34*(IQnormal30_34*wIQ85 + IQborderline30_34*wIQ70.85sim + IQmild30_34*wIQ50.69sim +
                             IQmoderate30_34*wIQ35.49sim + IQsevere30_34*wIQ20.34sim +
                             IQprofound30_34*wIQ20sim) * LE * women30_34 #per total women in DK in this age group

mean_median_ci(DALY30_34mc)


### Age 35-39 ###

set.seed(1)
ifexp <- sample(c(0,1), nvar, prob = probMeHg35_39, replace = T)
exp <- rgamma(nvar, fit2_35_39$estimate[1], fit2_35_39$estimate[2])

IQdiff35_39 <- apply(t(r), 2, function(x) x * ifexp * exp)

IQnew <- apply(IQdiff35_39, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 


IQnormal35_39 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline35_39 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild35_39 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate35_39 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere35_39 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound35_39 <- apply(IQnew, 2, function(x) mean(x < 20))

set.seed(1)
DALY35_39mc <- Pbaby35_39*(IQnormal35_39*wIQ85 + IQborderline35_39*wIQ70.85sim + IQmild35_39*wIQ50.69sim +
                             IQmoderate35_39*wIQ35.49sim + IQsevere35_39*wIQ20.34sim +
                             IQprofound35_39*wIQ20sim) * LE * women35_39 #per total women in DK in this age group

mean_median_ci(DALY35_39mc)


### Age 40-44 ###

set.seed(1)
ifexp <- sample(c(0,1), nvar, prob = probMeHg40_44, replace = T)
exp <- rgamma(nvar, fit2_40_44$estimate[1], fit2_40_44$estimate[2])

IQdiff40_44 <- apply(t(r), 2, function(x) x * ifexp * exp)

IQnew <- apply(IQdiff40_44, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal40_44 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline40_44 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild40_44 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate40_44 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere40_44 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound40_44 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY40_44mc <- Pbaby40_44*(IQnormal40_44*wIQ85 + IQborderline40_44*wIQ70.85sim + IQmild40_44*wIQ50.69sim +
                             IQmoderate40_44*wIQ35.49sim + IQsevere40_44*wIQ20.34sim +
                             IQprofound40_44*wIQ20sim) * LE * women40_44 #per total women in DK in this age group

mean_median_ci(DALY40_44mc)


### Age 45-49 ###

set.seed(1)
ifexp <- sample(c(0,1), nvar, prob = probMeHg45_49, replace = T)
exp <- rgamma(nvar, fit2_45_49$estimate[1], fit2_45_49$estimate[2])

IQdiff45_49 <- apply(t(r), 2, function(x) x * ifexp * exp)

IQnew <- apply(IQdiff45_49, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal45_49 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline45_49 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild45_49 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate45_49 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere45_49 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound45_49 <- apply(IQnew, 2, function(x) mean(x < 20))

set.seed(1)
DALY45_49mc <- Pbaby45_49*(IQnormal45_49*wIQ85 + IQborderline45_49*wIQ70.85sim + IQmild45_49*wIQ50.69sim +
                             IQmoderate45_49*wIQ35.49sim + IQsevere45_49*wIQ20.34sim +
                             IQprofound45_49*wIQ20sim) * LE * women45_49 #per total women in DK in this age group

mean_median_ci(DALY45_49mc)


##### Total DALYs ref #####

tDALY_ref_MeHg_IQ <- DALY15_19mc + DALY20_24mc + DALY25_29mc + DALY30_34mc + DALY35_39mc + DALY40_44mc + DALY45_49mc
mean_median_ci(tDALY_ref_MeHg_IQ)


IQchange <- (IQdiff15_19 + IQdiff20_24 + IQdiff25_29 + IQdiff30_34 + IQdiff35_39 + IQdiff40_44 + IQdiff45_49)/7


##### Number of cases #####

#Number of extra cases within each IQ class and total due to MeHg exposure
cases_borderline <- Pbaby15_19*IQborderline15_19*women15_19 + Pbaby20_24*IQborderline20_24*women20_24 + Pbaby25_29*IQborderline25_29*women25_29 + 
  Pbaby30_34*IQborderline30_34*women30_34 + Pbaby35_39*IQborderline35_39*women35_39 + Pbaby40_44*IQborderline40_44*women40_44 +
  Pbaby45_49*IQborderline45_49*women45_49
mean_median_ci(cases_borderline)


cases_mild <- Pbaby15_19*IQmild15_19*women15_19 + Pbaby20_24*IQmild20_24*women20_24 + Pbaby25_29*IQmild25_29*women25_29 + 
  Pbaby30_34*IQmild30_34*women30_34 + Pbaby35_39*IQmild35_39*women35_39 + Pbaby40_44*IQmild40_44*women40_44 +
  Pbaby45_49*IQmild45_49*women45_49
mean_median_ci(cases_mild)


cases_moderate <- Pbaby15_19*IQmoderate15_19*women15_19 + Pbaby20_24*IQmoderate20_24*women20_24 + Pbaby25_29*IQmoderate25_29*women25_29 + 
  Pbaby30_34*IQmoderate30_34*women30_34 + Pbaby35_39*IQmoderate35_39*women35_39 + Pbaby40_44*IQmoderate40_44*women40_44 +
  Pbaby45_49*IQmoderate45_49*women45_49
mean_median_ci(cases_moderate)


cases_severe <- Pbaby15_19*IQsevere15_19*women15_19 + Pbaby20_24*IQsevere20_24*women20_24 + Pbaby25_29*IQsevere25_29*women25_29 + 
  Pbaby30_34*IQsevere30_34*women30_34 + Pbaby35_39*IQsevere35_39*women35_39 + Pbaby40_44*IQsevere40_44*women40_44 +
  Pbaby45_49*IQsevere45_49*women45_49
mean_median_ci(cases_severe)


cases_profound <- Pbaby15_19*IQprofound15_19*women15_19 + Pbaby20_24*IQprofound20_24*women20_24 + Pbaby25_29*IQprofound25_29*women25_29 + 
  Pbaby30_34*IQprofound30_34*women30_34 + Pbaby35_39*IQprofound35_39*women35_39 + Pbaby40_44*IQprofound40_44*women40_44 +
  Pbaby45_49*IQprofound45_49*women45_49
mean_median_ci(cases_profound)



#Summing to total number of extra cases with intellectual disability weight due to MeHg exposure per 58,205 newborns
totalcases_ref <- cases_borderline + cases_mild + cases_moderate + cases_severe + cases_profound
mean_median_ci(totalcases_ref)



##### Scenario 1 #####

Scenario1 <- read.csv("Scenario1.csv")
# fish intakes are raw weights, ready for exposure assessment.

#Exclude males, and women above 49y of age (due to data on female fertility; women below 15y have already been excluded in FishDaily dataset)

Scenario1 <- subset(Scenario1, age < 50 & sex == '2')


Scen1.MeHg <- t(Scenario1[,c(55:70)]) #Create data set only with intakes of the different fish species
Scen1.MeHg <- t(Scen1.MeHg * MeHg_table[1:16,]) ##Multiply intakes of fish species with MeHg conc of the different species to obtain dataset with MeHg exposures from the different species


Scenario1$MeHg <- rowSums(Scen1.MeHg) #Add sum of MeHg exposures for each individual to the Scenario1 dataset


#MeHg/kg bw/week for the subpopulation
Scenario1$MeHg_bw <- Scenario1$MeHg / Scenario1$weight


Scenario1_IQ <- Scenario1[,c(1:4,52,71:72)]

setDT(Scenario1_IQ)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Exposures #####

#Divide data into agegroups. Use empirical distribution

#Age 15-19
MeHg15_19 <- subset(Scenario1_IQ$MeHg_bw, Scenario1_IQ$agegroups=="15-19")


#Age 20-24
MeHg20_24 <- subset(Scenario1_IQ$MeHg_bw, Scenario1_IQ$agegroups=="20-24")


#Age 25-29
MeHg25_29 <- subset(Scenario1_IQ$MeHg_bw, Scenario1_IQ$agegroups=="25-29")


#Age 30-34
MeHg30_34 <- subset(Scenario1_IQ$MeHg_bw, Scenario1_IQ$agegroups=="30-34")


#Age 35-39
MeHg35_39 <- subset(Scenario1_IQ$MeHg_bw, Scenario1_IQ$agegroups=="35-39")


#Age 40-44
MeHg40_44 <- subset(Scenario1_IQ$MeHg_bw, Scenario1_IQ$agegroups=="40-44")


#Age 45-49
MeHg45_49 <- subset(Scenario1_IQ$MeHg_bw, Scenario1_IQ$agegroups=="45-49")


##### DALY calculation #####

## Age 15-19 ##

set.seed(1)
exp <- sample(MeHg15_19, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff15_19 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff15_19, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal15_19 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline15_19 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild15_19 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate15_19 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere15_19 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound15_19 <- apply(IQnew, 2, function(x) mean(x < 20))

sum(mean(IQnormal15_19)+mean(IQborderline15_19)+mean(IQmild15_19)+mean(IQmoderate15_19)+mean(IQsevere15_19)+mean(IQprofound15_19))
# [1] 1


set.seed(1)
DALY15_19mc <- Pbaby15_19*(IQnormal15_19*wIQ85 + IQborderline15_19*wIQ70.85sim + IQmild15_19*wIQ50.69sim +
                             IQmoderate15_19*wIQ35.49sim + IQsevere15_19*wIQ20.34sim +
                             IQprofound15_19*wIQ20sim) * LE * women15_19 #per total women in DK in this age group

mean_median_ci(DALY15_19mc)



### Age 20-24 ###

set.seed(1)
exp <- sample(MeHg20_24, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff20_24 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff20_24, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal20_24 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline20_24 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild20_24 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate20_24 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere20_24 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound20_24 <- apply(IQnew, 2, function(x) mean(x < 20))

set.seed(1)
DALY20_24mc <- Pbaby20_24*(IQnormal20_24*wIQ85 + IQborderline20_24*wIQ70.85sim + IQmild20_24*wIQ50.69sim +
                             IQmoderate20_24*wIQ35.49sim + IQsevere20_24*wIQ20.34sim +
                             IQprofound20_24*wIQ20sim) * LE * women20_24 #per total women in DK in this age group

mean_median_ci(DALY20_24mc)


### Age 25-29 ###

set.seed(1)
exp <- sample(MeHg25_29, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff25_29 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff25_29, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal25_29 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline25_29 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild25_29 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate25_29 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere25_29 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound25_29 <- apply(IQnew, 2, function(x) mean(x < 20))

set.seed(1)
DALY25_29mc <- Pbaby25_29*(IQnormal25_29*wIQ85 + IQborderline25_29*wIQ70.85sim + IQmild25_29*wIQ50.69sim +
                             IQmoderate25_29*wIQ35.49sim + IQsevere25_29*wIQ20.34sim +
                             IQprofound25_29*wIQ20sim) * LE * women25_29 #per total women in DK in this age group

mean_median_ci(DALY25_29mc)



### Age 30-34 ###

set.seed(1)
exp <- sample(MeHg30_34, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff30_34 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff30_34, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal30_34 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline30_34 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild30_34 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate30_34 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere30_34 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound30_34 <- apply(IQnew, 2, function(x) mean(x < 20))

set.seed(1)
DALY30_34mc <- Pbaby30_34*(IQnormal30_34*wIQ85 + IQborderline30_34*wIQ70.85sim + IQmild30_34*wIQ50.69sim +
                             IQmoderate30_34*wIQ35.49sim + IQsevere30_34*wIQ20.34sim +
                             IQprofound30_34*wIQ20sim) * LE * women30_34 #per total women in DK in this age group

mean_median_ci(DALY30_34mc)


### Age 35-39 ###

set.seed(1)
exp <- sample(MeHg35_39, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff35_39 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff35_39, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal35_39 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline35_39 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild35_39 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate35_39 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere35_39 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound35_39 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY35_39mc <- Pbaby35_39*(IQnormal35_39*wIQ85 + IQborderline35_39*wIQ70.85sim + IQmild35_39*wIQ50.69sim +
                             IQmoderate35_39*wIQ35.49sim + IQsevere35_39*wIQ20.34sim +
                             IQprofound35_39*wIQ20sim) * LE * women35_39 #per total women in DK in this age group

mean_median_ci(DALY35_39mc)


### Age 40-44 ###

set.seed(1)
exp <- sample(MeHg40_44, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff40_44 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff40_44, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal40_44 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline40_44 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild40_44 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate40_44 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere40_44 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound40_44 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY40_44mc <- Pbaby40_44*(IQnormal40_44*wIQ85 + IQborderline40_44*wIQ70.85sim + IQmild40_44*wIQ50.69sim +
                             IQmoderate40_44*wIQ35.49sim + IQsevere40_44*wIQ20.34sim +
                             IQprofound40_44*wIQ20sim) * LE * women40_44 #per total women in DK in this age group

mean_median_ci(DALY40_44mc)


### Age 45-49 ###

set.seed(1)
exp <- sample(MeHg45_49, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff45_49 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff45_49, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal45_49 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline45_49 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild45_49 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate45_49 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere45_49 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound45_49 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY45_49mc <- Pbaby45_49*(IQnormal45_49*wIQ85 + IQborderline45_49*wIQ70.85sim + IQmild45_49*wIQ50.69sim +
                             IQmoderate45_49*wIQ35.49sim + IQsevere45_49*wIQ20.34sim +
                             IQprofound45_49*wIQ20sim) * LE * women45_49 #per total women in DK in this age group

mean_median_ci(DALY45_49mc)


##### Total DALYs #####

tDALY_scen1_MeHg_IQ <- DALY15_19mc + DALY20_24mc + DALY25_29mc + DALY30_34mc + DALY35_39mc + DALY40_44mc + DALY45_49mc
mean_median_ci(tDALY_scen1_MeHg_IQ)

median(tDALY_scen1_MeHg_IQ)


IQchange.scen1 <- (IQdiff15_19 + IQdiff20_24 + IQdiff25_29 + IQdiff30_34 + IQdiff35_39 + IQdiff40_44 + IQdiff45_49)/7



##### Number of cases #####

#Number of extra cases within each IQ class and total due to MeHg exposure
cases_borderline <- Pbaby15_19*IQborderline15_19*women15_19 + Pbaby20_24*IQborderline20_24*women20_24 + Pbaby25_29*IQborderline25_29*women25_29 + 
  Pbaby30_34*IQborderline30_34*women30_34 + Pbaby35_39*IQborderline35_39*women35_39 + Pbaby40_44*IQborderline40_44*women40_44 +
  Pbaby45_49*IQborderline45_49*women45_49
mean_median_ci(cases_borderline)


cases_mild <- Pbaby15_19*IQmild15_19*women15_19 + Pbaby20_24*IQmild20_24*women20_24 + Pbaby25_29*IQmild25_29*women25_29 + 
  Pbaby30_34*IQmild30_34*women30_34 + Pbaby35_39*IQmild35_39*women35_39 + Pbaby40_44*IQmild40_44*women40_44 +
  Pbaby45_49*IQmild45_49*women45_49
mean_median_ci(cases_mild)


cases_moderate <- Pbaby15_19*IQmoderate15_19*women15_19 + Pbaby20_24*IQmoderate20_24*women20_24 + Pbaby25_29*IQmoderate25_29*women25_29 + 
  Pbaby30_34*IQmoderate30_34*women30_34 + Pbaby35_39*IQmoderate35_39*women35_39 + Pbaby40_44*IQmoderate40_44*women40_44 +
  Pbaby45_49*IQmoderate45_49*women45_49
mean_median_ci(cases_moderate)


cases_severe <- Pbaby15_19*IQsevere15_19*women15_19 + Pbaby20_24*IQsevere20_24*women20_24 + Pbaby25_29*IQsevere25_29*women25_29 + 
  Pbaby30_34*IQsevere30_34*women30_34 + Pbaby35_39*IQsevere35_39*women35_39 + Pbaby40_44*IQsevere40_44*women40_44 +
  Pbaby45_49*IQsevere45_49*women45_49
mean_median_ci(cases_severe)


cases_profound <- Pbaby15_19*IQprofound15_19*women15_19 + Pbaby20_24*IQprofound20_24*women20_24 + Pbaby25_29*IQprofound25_29*women25_29 + 
  Pbaby30_34*IQprofound30_34*women30_34 + Pbaby35_39*IQprofound35_39*women35_39 + Pbaby40_44*IQprofound40_44*women40_44 +
  Pbaby45_49*IQprofound45_49*women45_49
mean_median_ci(cases_profound)


#Summing to total number of extra cases with intellectual disability weight due to MeHg exposure per 58,205 newborns
totalcases_scen1 <- cases_borderline + cases_mild + cases_moderate + cases_severe + cases_profound
mean_median_ci(totalcases_scen1)



##### Scenario 2 #####

Scenario2 <- read.csv("Scenario2.csv")
# fish intakes are raw weights, ready for exposure assessment.

#Exclude males, and women above 49y of age (due to data on female fertility; women below 15y have already been excluded in FishDaily dataset)

Scenario2 <- subset(Scenario2, age < 50 & sex == '2')


Scen2.MeHg <- t(Scenario2[,c(55:60)]) #Create data set only with intakes of the different fish species
Scen2.MeHg <- t(Scen2.MeHg * MeHg_table[c(3,9,10,11,12,16),]) ##Multiply intakes of fish species with MeHg conc of the different species to obtain dataset with MeHg exposures from the different species


Scenario2$MeHg <- rowSums(Scen2.MeHg) #Add sum of MeHg exposures for each individual to the Scenario2 dataset

#MeHg/kg bw/week 
Scenario2$MeHg_bw <- Scenario2$MeHg / Scenario2$weight


Scenario2_IQ <- Scenario2[,c(1:4,54,61:62)]

setDT(Scenario2_IQ)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]



##### Exposures #####

#Divide data into agegroups. Use empirical distribution

#Age 15-19

MeHg15_19 <- subset(Scenario2_IQ$MeHg_bw, Scenario2_IQ$agegroups=="15-19")

#We don't need to model the probability of MeHg exposure since all individuals now consume fish

#Age 20-24
MeHg20_24 <- subset(Scenario2_IQ$MeHg_bw, Scenario2_IQ$agegroups=="20-24")


#Age 25-29
MeHg25_29 <- subset(Scenario2_IQ$MeHg_bw, Scenario2_IQ$agegroups=="25-29")


#Age 30-34
MeHg30_34 <- subset(Scenario2_IQ$MeHg_bw, Scenario2_IQ$agegroups=="30-34")


#Age 35-39
MeHg35_39 <- subset(Scenario2_IQ$MeHg_bw, Scenario2_IQ$agegroups=="35-39")


#Age 40-44
MeHg40_44 <- subset(Scenario2_IQ$MeHg_bw, Scenario2_IQ$agegroups=="40-44")


#Age 45-49
MeHg45_49 <- subset(Scenario2_IQ$MeHg_bw, Scenario2_IQ$agegroups=="45-49")


###### DALY calculation #####

#Age 15-49

set.seed(1)
exp <- sample(MeHg15_19, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff15_19 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff15_19, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal15_19 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline15_19 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild15_19 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate15_19 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere15_19 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound15_19 <- apply(IQnew, 2, function(x) mean(x < 20))

sum(mean(IQnormal15_19)+mean(IQborderline15_19)+mean(IQmild15_19)+mean(IQmoderate15_19)+mean(IQsevere15_19)+mean(IQprofound15_19))
# [1] 1


set.seed(1)
DALY15_19mc <- Pbaby15_19*(IQnormal15_19*wIQ85 + IQborderline15_19*wIQ70.85sim + IQmild15_19*wIQ50.69sim +
                             IQmoderate15_19*wIQ35.49sim + IQsevere15_19*wIQ20.34sim +
                             IQprofound15_19*wIQ20sim) * LE * women15_19 #per total women in DK in this age group

mean_median_ci(DALY15_19mc)



### Age 20-24 ###

set.seed(1)
exp <- sample(MeHg20_24, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff20_24 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff20_24, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal20_24 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline20_24 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild20_24 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate20_24 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere20_24 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound20_24 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY20_24mc <- Pbaby20_24*(IQnormal20_24*wIQ85 + IQborderline20_24*wIQ70.85sim + IQmild20_24*wIQ50.69sim +
                             IQmoderate20_24*wIQ35.49sim + IQsevere20_24*wIQ20.34sim +
                             IQprofound20_24*wIQ20sim) * LE * women20_24 #per total women in DK in this age group

mean_median_ci(DALY20_24mc)


### Age 25-29 ###

set.seed(1)
exp <- sample(MeHg25_29, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff25_29 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff25_29, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal25_29 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline25_29 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild25_29 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate25_29 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere25_29 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound25_29 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY25_29mc <- Pbaby25_29*(IQnormal25_29*wIQ85 + IQborderline25_29*wIQ70.85sim + IQmild25_29*wIQ50.69sim +
                             IQmoderate25_29*wIQ35.49sim + IQsevere25_29*wIQ20.34sim +
                             IQprofound25_29*wIQ20sim) * LE * women25_29 #per total women in DK in this age group

mean_median_ci(DALY25_29mc)


### Age 30-34 ###

set.seed(1)
exp <- sample(MeHg30_34, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff30_34 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff30_34, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal30_34 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline30_34 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild30_34 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate30_34 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere30_34 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound30_34 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY30_34mc <- Pbaby30_34*(IQnormal30_34*wIQ85 + IQborderline30_34*wIQ70.85sim + IQmild30_34*wIQ50.69sim +
                             IQmoderate30_34*wIQ35.49sim + IQsevere30_34*wIQ20.34sim +
                             IQprofound30_34*wIQ20sim) * LE * women30_34 #per total women in DK in this age group

mean_median_ci(DALY30_34mc)


### Age 35-39 ###

set.seed(1)
exp <- sample(MeHg35_39, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff35_39 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff35_39, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal35_39 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline35_39 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild35_39 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate35_39 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere35_39 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound35_39 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY35_39mc <- Pbaby35_39*(IQnormal35_39*wIQ85 + IQborderline35_39*wIQ70.85sim + IQmild35_39*wIQ50.69sim +
                             IQmoderate35_39*wIQ35.49sim + IQsevere35_39*wIQ20.34sim +
                             IQprofound35_39*wIQ20sim) * LE * women35_39 #per total women in DK in this age group

mean_median_ci(DALY35_39mc)


### Age 40-44 ###

set.seed(1)
exp <- sample(MeHg40_44, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff40_44 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff40_44, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal40_44 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline40_44 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild40_44 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate40_44 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere40_44 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound40_44 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY40_44mc <- Pbaby40_44*(IQnormal40_44*wIQ85 + IQborderline40_44*wIQ70.85sim + IQmild40_44*wIQ50.69sim +
                             IQmoderate40_44*wIQ35.49sim + IQsevere40_44*wIQ20.34sim +
                             IQprofound40_44*wIQ20sim) * LE * women40_44 #per total women in DK in this age group

mean_median_ci(DALY40_44mc)


### Age 45-49 ###

set.seed(1)
exp <- sample(MeHg45_49, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff45_49 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff45_49, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal45_49 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline45_49 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild45_49 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate45_49 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere45_49 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound45_49 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY45_49mc <- Pbaby45_49*(IQnormal45_49*wIQ85 + IQborderline45_49*wIQ70.85sim + IQmild45_49*wIQ50.69sim +
                             IQmoderate45_49*wIQ35.49sim + IQsevere45_49*wIQ20.34sim +
                             IQprofound45_49*wIQ20sim) * LE * women45_49 #per total women in DK in this age group

mean_median_ci(DALY45_49mc)


##### Total DALYs #####

tDALY_scen2_MeHg_IQ <- DALY15_19mc + DALY20_24mc + DALY25_29mc + DALY30_34mc + DALY35_39mc + DALY40_44mc + DALY45_49mc
mean_median_ci(tDALY_scen2_MeHg_IQ)


IQchange.scen2 <- (IQdiff15_19 + IQdiff20_24 + IQdiff25_29 + IQdiff30_34 + IQdiff35_39 + IQdiff40_44 + IQdiff45_49)/7



##### Number of cases #####

#Number of extra cases within each IQ class and total due to MeHg exposure
cases_borderline <- Pbaby15_19*IQborderline15_19*women15_19 + Pbaby20_24*IQborderline20_24*women20_24 + Pbaby25_29*IQborderline25_29*women25_29 + 
  Pbaby30_34*IQborderline30_34*women30_34 + Pbaby35_39*IQborderline35_39*women35_39 + Pbaby40_44*IQborderline40_44*women40_44 +
  Pbaby45_49*IQborderline45_49*women45_49
mean_median_ci(cases_borderline)


cases_mild <- Pbaby15_19*IQmild15_19*women15_19 + Pbaby20_24*IQmild20_24*women20_24 + Pbaby25_29*IQmild25_29*women25_29 + 
  Pbaby30_34*IQmild30_34*women30_34 + Pbaby35_39*IQmild35_39*women35_39 + Pbaby40_44*IQmild40_44*women40_44 +
  Pbaby45_49*IQmild45_49*women45_49
mean_median_ci(cases_mild)


cases_moderate <- Pbaby15_19*IQmoderate15_19*women15_19 + Pbaby20_24*IQmoderate20_24*women20_24 + Pbaby25_29*IQmoderate25_29*women25_29 + 
  Pbaby30_34*IQmoderate30_34*women30_34 + Pbaby35_39*IQmoderate35_39*women35_39 + Pbaby40_44*IQmoderate40_44*women40_44 +
  Pbaby45_49*IQmoderate45_49*women45_49
mean_median_ci(cases_moderate)


cases_severe <- Pbaby15_19*IQsevere15_19*women15_19 + Pbaby20_24*IQsevere20_24*women20_24 + Pbaby25_29*IQsevere25_29*women25_29 + 
  Pbaby30_34*IQsevere30_34*women30_34 + Pbaby35_39*IQsevere35_39*women35_39 + Pbaby40_44*IQsevere40_44*women40_44 +
  Pbaby45_49*IQsevere45_49*women45_49
mean_median_ci(cases_severe)


cases_profound <- Pbaby15_19*IQprofound15_19*women15_19 + Pbaby20_24*IQprofound20_24*women20_24 + Pbaby25_29*IQprofound25_29*women25_29 + 
  Pbaby30_34*IQprofound30_34*women30_34 + Pbaby35_39*IQprofound35_39*women35_39 + Pbaby40_44*IQprofound40_44*women40_44 +
  Pbaby45_49*IQprofound45_49*women45_49
mean_median_ci(cases_profound)


#Summing to total number of extra cases with intellectual disability weight due to MeHg exposure per 58,205 newborns
totalcases_scen2 <- cases_borderline + cases_mild + cases_moderate + cases_severe + cases_profound
mean_median_ci(totalcases_scen2)



###### Scenario 3 #####

Scenario3 <- read.csv("Scenario3.csv")
# fish intakes are raw weights, ready for exposure assessment.

#Exclude males, and women above 49y of age (due to data on female fertility; women below 15y have already been excluded in FishDaily dataset)

Scenario3 <- subset(Scenario3, age < 50 & sex == '2')

Scen3.MeHg <- t(Scenario3[,c(55:64)]) #Create data set only with intakes of the different fish species
Scen3.MeHg <- t(Scen3.MeHg * MeHg_table[c(1,2,4,5,6,7,8,13,14,15),]) ##Multiply intakes of fish species with MeHg conc of the different species to obtain dataset with MeHg exposures from the different species

Scenario3$MeHg <- rowSums(Scen3.MeHg) #Add sum of MeHg exposures for each individual to the Scenario3 dataset


#MeHg/kg bw/week 
Scenario3$MeHg_bw <- Scenario3$MeHg / Scenario3$weight


Scenario3_IQ <- Scenario3[,c(1:4,54,65:66)]

setDT(Scenario3_IQ)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Exposures #####

#Divide data into agegroups. Use empirical distribution

#Age 15-19
MeHg15_19 <- subset(Scenario3_IQ$MeHg_bw, Scenario3_IQ$agegroups=="15-19")

#We don't need to model the probability of MeHg exposure since all individuals now consume fish


#Age 20-24
MeHg20_24 <- subset(Scenario3_IQ$MeHg_bw, Scenario3_IQ$agegroups=="20-24")


#Age 25-29
MeHg25_29 <- subset(Scenario3_IQ$MeHg_bw, Scenario3_IQ$agegroups=="25-29")


#Age 30-34
MeHg30_34 <- subset(Scenario3_IQ$MeHg_bw, Scenario3_IQ$agegroups=="30-34")


#Age 35-39
MeHg35_39 <- subset(Scenario3_IQ$MeHg_bw, Scenario3_IQ$agegroups=="35-39")


#Age 40-44
MeHg40_44 <- subset(Scenario3_IQ$MeHg_bw, Scenario3_IQ$agegroups=="40-44")


#Age 45-49
MeHg45_49 <- subset(Scenario3_IQ$MeHg_bw, Scenario3_IQ$agegroups=="45-49")



##### DALY calculation #####

## Age 15-19 ##

set.seed(1)
exp <- sample(MeHg15_19, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff15_19 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff15_19, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal15_19 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline15_19 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild15_19 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate15_19 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere15_19 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound15_19 <- apply(IQnew, 2, function(x) mean(x < 20))

sum(mean(IQnormal15_19)+mean(IQborderline15_19)+mean(IQmild15_19)+mean(IQmoderate15_19)+mean(IQsevere15_19)+mean(IQprofound15_19))
# [1] 1


set.seed(1)
DALY15_19mc <- Pbaby15_19*(IQnormal15_19*wIQ85 + IQborderline15_19*wIQ70.85sim + IQmild15_19*wIQ50.69sim +
                             IQmoderate15_19*wIQ35.49sim + IQsevere15_19*wIQ20.34sim +
                             IQprofound15_19*wIQ20sim) * LE * women15_19 #per total women in DK in this age group

mean_median_ci(DALY15_19mc)



### Age 20-24 ###

set.seed(1)
exp <- sample(MeHg20_24, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff20_24 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff20_24, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal20_24 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline20_24 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild20_24 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate20_24 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere20_24 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound20_24 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY20_24mc <- Pbaby20_24*(IQnormal20_24*wIQ85 + IQborderline20_24*wIQ70.85sim + IQmild20_24*wIQ50.69sim +
                             IQmoderate20_24*wIQ35.49sim + IQsevere20_24*wIQ20.34sim +
                             IQprofound20_24*wIQ20sim) * LE * women20_24 #per total women in DK in this age group

mean_median_ci(DALY20_24mc)


### Age 25-29 ###

set.seed(1)
exp <- sample(MeHg25_29, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff25_29 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff25_29, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal25_29 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline25_29 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild25_29 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate25_29 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere25_29 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound25_29 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY25_29mc <- Pbaby25_29*(IQnormal25_29*wIQ85 + IQborderline25_29*wIQ70.85sim + IQmild25_29*wIQ50.69sim +
                             IQmoderate25_29*wIQ35.49sim + IQsevere25_29*wIQ20.34sim +
                             IQprofound25_29*wIQ20sim) * LE * women25_29 #per total women in DK in this age group

mean_median_ci(DALY25_29mc)



### Age 30-34 ###

set.seed(1)
exp <- sample(MeHg30_34, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff30_34 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff30_34, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal30_34 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline30_34 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild30_34 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate30_34 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere30_34 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound30_34 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY30_34mc <- Pbaby30_34*(IQnormal30_34*wIQ85 + IQborderline30_34*wIQ70.85sim + IQmild30_34*wIQ50.69sim +
                             IQmoderate30_34*wIQ35.49sim + IQsevere30_34*wIQ20.34sim +
                             IQprofound30_34*wIQ20sim) * LE * women30_34 #per total women in DK in this age group

mean_median_ci(DALY30_34mc)


### Age 35-39 ###

set.seed(1)
exp <- sample(MeHg35_39, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff35_39 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff35_39, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal35_39 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline35_39 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild35_39 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate35_39 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere35_39 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound35_39 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY35_39mc <- Pbaby35_39*(IQnormal35_39*wIQ85 + IQborderline35_39*wIQ70.85sim + IQmild35_39*wIQ50.69sim +
                             IQmoderate35_39*wIQ35.49sim + IQsevere35_39*wIQ20.34sim +
                             IQprofound35_39*wIQ20sim) * LE * women35_39 #per total women in DK in this age group

mean_median_ci(DALY35_39mc)


### Age 40-44 ###

set.seed(1)
exp <- sample(MeHg40_44, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff40_44 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff40_44, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal40_44 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline40_44 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild40_44 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate40_44 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere40_44 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound40_44 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY40_44mc <- Pbaby40_44*(IQnormal40_44*wIQ85 + IQborderline40_44*wIQ70.85sim + IQmild40_44*wIQ50.69sim +
                             IQmoderate40_44*wIQ35.49sim + IQsevere40_44*wIQ20.34sim +
                             IQprofound40_44*wIQ20sim) * LE * women40_44 #per total women in DK in this age group

mean_median_ci(DALY40_44mc)



### Age 45-49 ###

set.seed(1)
exp <- sample(MeHg45_49, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff45_49 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff45_49, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal45_49 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline45_49 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild45_49 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate45_49 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere45_49 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound45_49 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY45_49mc <- Pbaby45_49*(IQnormal45_49*wIQ85 + IQborderline45_49*wIQ70.85sim + IQmild45_49*wIQ50.69sim +
                             IQmoderate45_49*wIQ35.49sim + IQsevere45_49*wIQ20.34sim +
                             IQprofound45_49*wIQ20sim) * LE * women45_49 #per total women in DK in this age group

mean_median_ci(DALY45_49mc)



##### Total DALY #####

tDALY_scen3_MeHg_IQ <- DALY15_19mc + DALY20_24mc + DALY25_29mc + DALY30_34mc + DALY35_39mc + DALY40_44mc + DALY45_49mc
mean_median_ci(tDALY_scen3_MeHg_IQ)


IQchange.scen3 <- (IQdiff15_19 + IQdiff20_24 + IQdiff25_29 + IQdiff30_34 + IQdiff35_39 + IQdiff40_44 + IQdiff45_49)/7


##### Number of cases #####

#Number of extra cases within each IQ class and total due to MeHg exposure
cases_borderline <- Pbaby15_19*IQborderline15_19*women15_19 + Pbaby20_24*IQborderline20_24*women20_24 + Pbaby25_29*IQborderline25_29*women25_29 + 
  Pbaby30_34*IQborderline30_34*women30_34 + Pbaby35_39*IQborderline35_39*women35_39 + Pbaby40_44*IQborderline40_44*women40_44 +
  Pbaby45_49*IQborderline45_49*women45_49
mean_median_ci(cases_borderline)


cases_mild <- Pbaby15_19*IQmild15_19*women15_19 + Pbaby20_24*IQmild20_24*women20_24 + Pbaby25_29*IQmild25_29*women25_29 + 
  Pbaby30_34*IQmild30_34*women30_34 + Pbaby35_39*IQmild35_39*women35_39 + Pbaby40_44*IQmild40_44*women40_44 +
  Pbaby45_49*IQmild45_49*women45_49
mean_median_ci(cases_mild)


cases_moderate <- Pbaby15_19*IQmoderate15_19*women15_19 + Pbaby20_24*IQmoderate20_24*women20_24 + Pbaby25_29*IQmoderate25_29*women25_29 + 
  Pbaby30_34*IQmoderate30_34*women30_34 + Pbaby35_39*IQmoderate35_39*women35_39 + Pbaby40_44*IQmoderate40_44*women40_44 +
  Pbaby45_49*IQmoderate45_49*women45_49
mean_median_ci(cases_moderate)


cases_severe <- Pbaby15_19*IQsevere15_19*women15_19 + Pbaby20_24*IQsevere20_24*women20_24 + Pbaby25_29*IQsevere25_29*women25_29 + 
  Pbaby30_34*IQsevere30_34*women30_34 + Pbaby35_39*IQsevere35_39*women35_39 + Pbaby40_44*IQsevere40_44*women40_44 +
  Pbaby45_49*IQsevere45_49*women45_49
mean_median_ci(cases_severe)


cases_profound <- Pbaby15_19*IQprofound15_19*women15_19 + Pbaby20_24*IQprofound20_24*women20_24 + Pbaby25_29*IQprofound25_29*women25_29 + 
  Pbaby30_34*IQprofound30_34*women30_34 + Pbaby35_39*IQprofound35_39*women35_39 + Pbaby40_44*IQprofound40_44*women40_44 +
  Pbaby45_49*IQprofound45_49*women45_49
mean_median_ci(cases_profound)


#Summing to total number of extra cases with intellectual disability weight due to MeHg exposure per 58,205 newborns
totalcases_scen3 <- cases_borderline + cases_mild + cases_moderate + cases_severe + cases_profound
mean_median_ci(totalcases_scen3)



##### Scenario 4 #####

Scenario4 <- read.csv("Scenario4.csv")
# fish intakes are raw weights, ready for exposure assessment.

#Exclude males, and women above 49y of age (due to data on female fertility; women below 15y have already been excluded in FishDaily dataset)

Scenario4 <- subset(Scenario4, age < 50 & sex == '2')


Scen4.MeHg <- t(Scenario4[,c(55)]) #Create data set only with intakes of the different fish species
Scen4.MeHg <- t(Scen4.MeHg * MeHg_table[8,]) ##Multiply intakes of fish species with MeHg conc of the different species to obtain dataset with MeHg exposures from the different species

Scenario4$MeHg <- rowSums(Scen4.MeHg) #Add sum of MeHg exposures for each individual to the Scenario3 dataset

#MeHg/kg bw/week 
Scenario4$MeHg_bw <- Scenario4$MeHg / Scenario4$weight


Scenario4_IQ <- Scenario4[,c(1:4,54,56:57)]

setDT(Scenario4_IQ)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Exposures #####

#Divide data into agegroups. Use empirical distribution

#Age 15-19
MeHg15_19 <- subset(Scenario4_IQ$MeHg_bw, Scenario4_IQ$agegroups=="15-19")
#We don't need to model the probability of MeHg exposure since all individuals now consume fish


#Age 20-24
MeHg20_24 <- subset(Scenario4_IQ$MeHg_bw, Scenario4_IQ$agegroups=="20-24")


#Age 25-29
MeHg25_29 <- subset(Scenario4_IQ$MeHg_bw, Scenario4_IQ$agegroups=="25-29")


#Age 30-34
MeHg30_34 <- subset(Scenario4_IQ$MeHg_bw, Scenario4_IQ$agegroups=="30-34")


#Age 35-39
MeHg35_39 <- subset(Scenario4_IQ$MeHg_bw, Scenario4_IQ$agegroups=="35-39")


#Age 40-44
MeHg40_44 <- subset(Scenario4_IQ$MeHg_bw, Scenario4_IQ$agegroups=="40-44")


#Age 45-49
MeHg45_49 <- subset(Scenario4_IQ$MeHg_bw, Scenario4_IQ$agegroups=="45-49")



##### DALY calculation #####

## Age 15-19 ##

set.seed(1)
exp <- sample(MeHg15_19, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff15_19 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff15_19, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal15_19 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline15_19 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild15_19 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate15_19 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere15_19 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound15_19 <- apply(IQnew, 2, function(x) mean(x < 20))

sum(mean(IQnormal15_19)+mean(IQborderline15_19)+mean(IQmild15_19)+mean(IQmoderate15_19)+mean(IQsevere15_19)+mean(IQprofound15_19))
# [1] 1


set.seed(1)
DALY15_19mc <- Pbaby15_19*(IQnormal15_19*wIQ85 + IQborderline15_19*wIQ70.85sim + IQmild15_19*wIQ50.69sim +
                             IQmoderate15_19*wIQ35.49sim + IQsevere15_19*wIQ20.34sim +
                             IQprofound15_19*wIQ20sim) * LE * women15_19 #per total women in DK in this age group

mean_median_ci(DALY15_19mc)


### Age 20-24 ###

set.seed(1)
exp <- sample(MeHg20_24, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff20_24 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff20_24, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal20_24 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline20_24 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild20_24 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate20_24 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere20_24 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound20_24 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY20_24mc <- Pbaby20_24*(IQnormal20_24*wIQ85 + IQborderline20_24*wIQ70.85sim + IQmild20_24*wIQ50.69sim +
                             IQmoderate20_24*wIQ35.49sim + IQsevere20_24*wIQ20.34sim +
                             IQprofound20_24*wIQ20sim) * LE * women20_24 #per total women in DK in this age group

mean_median_ci(DALY20_24mc)


### Age 25-29 ###

set.seed(1)
exp <- sample(MeHg25_29, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff25_29 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff25_29, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal25_29 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline25_29 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild25_29 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate25_29 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere25_29 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound25_29 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY25_29mc <- Pbaby25_29*(IQnormal25_29*wIQ85 + IQborderline25_29*wIQ70.85sim + IQmild25_29*wIQ50.69sim +
                             IQmoderate25_29*wIQ35.49sim + IQsevere25_29*wIQ20.34sim +
                             IQprofound25_29*wIQ20sim) * LE * women25_29 #per total women in DK in this age group

mean_median_ci(DALY25_29mc)



### Age 30-34 ###

set.seed(1)
exp <- sample(MeHg30_34, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff30_34 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff30_34, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal30_34 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline30_34 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild30_34 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate30_34 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere30_34 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound30_34 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY30_34mc <- Pbaby30_34*(IQnormal30_34*wIQ85 + IQborderline30_34*wIQ70.85sim + IQmild30_34*wIQ50.69sim +
                             IQmoderate30_34*wIQ35.49sim + IQsevere30_34*wIQ20.34sim +
                             IQprofound30_34*wIQ20sim) * LE * women30_34 #per total women in DK in this age group

mean_median_ci(DALY30_34mc)


### Age 35-39 ###

set.seed(1)
exp <- sample(MeHg35_39, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff35_39 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff35_39, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal35_39 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline35_39 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild35_39 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate35_39 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere35_39 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound35_39 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY35_39mc <- Pbaby35_39*(IQnormal35_39*wIQ85 + IQborderline35_39*wIQ70.85sim + IQmild35_39*wIQ50.69sim +
                             IQmoderate35_39*wIQ35.49sim + IQsevere35_39*wIQ20.34sim +
                             IQprofound35_39*wIQ20sim) * LE * women35_39 #per total women in DK in this age group

mean_median_ci(DALY35_39mc)


### Age 40-44 ###

set.seed(1)
exp <- sample(MeHg40_44, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff40_44 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff40_44, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal40_44 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline40_44 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild40_44 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate40_44 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere40_44 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound40_44 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY40_44mc <- Pbaby40_44*(IQnormal40_44*wIQ85 + IQborderline40_44*wIQ70.85sim + IQmild40_44*wIQ50.69sim +
                             IQmoderate40_44*wIQ35.49sim + IQsevere40_44*wIQ20.34sim +
                             IQprofound40_44*wIQ20sim) * LE * women40_44 #per total women in DK in this age group

mean_median_ci(DALY40_44mc)



### Age 45-49 ###

set.seed(1)
exp <- sample(MeHg45_49, nvar, replace = T) #Everyone are exposed so ifexp = 1

IQdiff45_49 <- apply(t(r), 2, function(x) x *  exp)

IQnew <- apply(IQdiff45_49, 2, function(x) x + IQnorm) #Add IQ dist (var) each col (uncertainty) 

IQnormal45_49 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline45_49 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild45_49 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate45_49 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere45_49 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound45_49 <- apply(IQnew, 2, function(x) mean(x < 20))


set.seed(1)
DALY45_49mc <- Pbaby45_49*(IQnormal45_49*wIQ85 + IQborderline45_49*wIQ70.85sim + IQmild45_49*wIQ50.69sim +
                             IQmoderate45_49*wIQ35.49sim + IQsevere45_49*wIQ20.34sim +
                             IQprofound45_49*wIQ20sim) * LE * women45_49 #per total women in DK in this age group

mean_median_ci(DALY45_49mc)



##### Total DALYs #####

tDALY_scen4_MeHg_IQ <- DALY15_19mc + DALY20_24mc + DALY25_29mc + DALY30_34mc + DALY35_39mc + DALY40_44mc + DALY45_49mc
mean_median_ci(tDALY_scen4_MeHg_IQ)



IQchange.scen4 <- (IQdiff15_19 + IQdiff20_24 + IQdiff25_29 + IQdiff30_34 + IQdiff35_39 + IQdiff40_44 + IQdiff45_49)/7



##### Number of cases #####

#Number of extra cases within each IQ class and total due to MeHg exposure
cases_borderline <- Pbaby15_19*IQborderline15_19*women15_19 + Pbaby20_24*IQborderline20_24*women20_24 + Pbaby25_29*IQborderline25_29*women25_29 + 
  Pbaby30_34*IQborderline30_34*women30_34 + Pbaby35_39*IQborderline35_39*women35_39 + Pbaby40_44*IQborderline40_44*women40_44 +
  Pbaby45_49*IQborderline45_49*women45_49
mean_median_ci(cases_borderline)


cases_mild <- Pbaby15_19*IQmild15_19*women15_19 + Pbaby20_24*IQmild20_24*women20_24 + Pbaby25_29*IQmild25_29*women25_29 + 
  Pbaby30_34*IQmild30_34*women30_34 + Pbaby35_39*IQmild35_39*women35_39 + Pbaby40_44*IQmild40_44*women40_44 +
  Pbaby45_49*IQmild45_49*women45_49
mean_median_ci(cases_mild)


cases_moderate <- Pbaby15_19*IQmoderate15_19*women15_19 + Pbaby20_24*IQmoderate20_24*women20_24 + Pbaby25_29*IQmoderate25_29*women25_29 + 
  Pbaby30_34*IQmoderate30_34*women30_34 + Pbaby35_39*IQmoderate35_39*women35_39 + Pbaby40_44*IQmoderate40_44*women40_44 +
  Pbaby45_49*IQmoderate45_49*women45_49
mean_median_ci(cases_moderate)


cases_severe <- Pbaby15_19*IQsevere15_19*women15_19 + Pbaby20_24*IQsevere20_24*women20_24 + Pbaby25_29*IQsevere25_29*women25_29 + 
  Pbaby30_34*IQsevere30_34*women30_34 + Pbaby35_39*IQsevere35_39*women35_39 + Pbaby40_44*IQsevere40_44*women40_44 +
  Pbaby45_49*IQsevere45_49*women45_49
mean_median_ci(cases_severe)


cases_profound <- Pbaby15_19*IQprofound15_19*women15_19 + Pbaby20_24*IQprofound20_24*women20_24 + Pbaby25_29*IQprofound25_29*women25_29 + 
  Pbaby30_34*IQprofound30_34*women30_34 + Pbaby35_39*IQprofound35_39*women35_39 + Pbaby40_44*IQprofound40_44*women40_44 +
  Pbaby45_49*IQprofound45_49*women45_49
mean_median_ci(cases_profound)


#Summing to total number of extra cases with intellectual disability weight due to MeHg exposure per 58,205 newborns
totalcases_scen4 <- cases_borderline + cases_mild + cases_moderate + cases_severe + cases_profound
mean_median_ci(totalcases_scen4)



##### Extra number of cases #####

cases_diff_scen1 <- totalcases_scen1 - totalcases_ref
mean_median_ci(cases_diff_scen1)


cases_diff_scen2 <- totalcases_scen2 - totalcases_ref
mean_median_ci(cases_diff_scen2)


cases_diff_scen3 <- totalcases_scen3 - totalcases_ref
mean_median_ci(cases_diff_scen3)


cases_diff_scen4 <- totalcases_scen4 - totalcases_ref
mean_median_ci(cases_diff_scen4)



##### Delta DALY #####

dDALY_MeHg_IQ_scen1 <- tDALY_scen1_MeHg_IQ - tDALY_ref_MeHg_IQ
mean_median_ci(dDALY_MeHg_IQ_scen1)


dDALY_MeHg_IQ_scen2 <- tDALY_scen2_MeHg_IQ - tDALY_ref_MeHg_IQ
mean_median_ci(dDALY_MeHg_IQ_scen2)


dDALY_MeHg_IQ_scen3 <- tDALY_scen3_MeHg_IQ - tDALY_ref_MeHg_IQ
mean_median_ci(dDALY_MeHg_IQ_scen3)


dDALY_MeHg_IQ_scen4 <- tDALY_scen4_MeHg_IQ - tDALY_ref_MeHg_IQ
mean_median_ci(dDALY_MeHg_IQ_scen4)



##### Save #####

write.csv(tDALY_ref_MeHg_IQ, "tDALY_ref_MeHg_IQ.csv")

write.csv(tDALY_scen1_MeHg_IQ, "tDALY_scen1_MeHg_IQ.csv")

write.csv(tDALY_scen2_MeHg_IQ, "tDALY_scen2_MeHg_IQ.csv")

write.csv(tDALY_scen3_MeHg_IQ, "tDALY_scen3_MeHg_IQ.csv")

write.csv(tDALY_scen4_MeHg_IQ, "tDALY_scen4_MeHg_IQ.csv")


IQchange.fishIQ.scen1_unc <- apply(IQchange.scen1, 2, function(x) x*1) #uncertainty dimension
write.csv(IQchange.fishIQ.scen1_unc, "IQchange.fishIQ.scen1.csv")

IQchange.fishIQ.scen2_unc <- apply(IQchange.scen2, 2, function(x) x*1) #uncertainty dimension
write.csv(IQchange.fishIQ.scen2_unc, "IQchange.fishIQ.scen2.csv")

IQchange.fishIQ.scen3_unc <- apply(IQchange.scen3, 2, function(x) x*1) #uncertainty dimension
write.csv(IQchange.fishIQ.scen3_unc, "IQchange.fishIQ.scen3.csv")

IQchange.fishIQ.scen4_unc <- apply(IQchange.scen4, 2, function(x) x*1) #uncertainty dimension
write.csv(IQchange.fishIQ.scen4_unc, "IQchange.fishIQ.scen4.csv")
