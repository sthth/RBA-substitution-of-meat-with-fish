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
r <- rpert(nunc, min = 0.104, mode = 0.152, max = 0.212) #Dose-response Hibbeln et al., 2007
IQnorm <- rnorm(nvar, mean = 100, sd = 15) #Definition of IQ



##### Reference scenario #####

FishDaily <- read.csv("FishDaily.csv")


#FoodDaily/FishDaily fish amounts are PREPARED AMOUNTS for the ref scenario - prepared amounts also available in alt scenarios.
#No need to convert back into prepared weights
#Webpanel 1 for the Hibbeln paper shows that the dose-response is based on prepared fish intakes


#Exclude males, and women above 49y of age (due to data on female fertility) - women below 15 have already been excluded

FishDaily <- subset(FishDaily, sex == 2 & age < 50)


agebreaks <- c(15,20,25,30,35,40,45,50)
agelabels= c("15-19","20-24","25-29","30-34","35-39","40-44","45-49")

setDT(FishDaily)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Fitting exposures #####

#Divide data into agegroups and fit exposure to distribution

#Age 15-19
Fish15_19 <- subset(FishDaily, agegroups == "15-19", select = total.fish)

#Probability of fish intake
probfish15_19 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Fish15_19==0)/length(Fish15_19$total.fish) #probability of zero DHA exp
pr.yes <- sum(Fish15_19!=0)/length(Fish15_19$total.fish) #probability of DHA exp
probfish15_19[1,] <- c(pr.no,pr.yes)



Fish15_19_pos <- Fish15_19[which(Fish15_19!=0)] #Only intakes that are greater than zero
Fish15_19_pos <- as.vector(Fish15_19_pos$total.fish)

fit1_15_19 <- fitdist(Fish15_19_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Fish15_19_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_15_19$estimate[1], sdlog=fit1_15_19$estimate[2]), col="red", lty=1)
gofstat(fit1_15_19, fitnames="lnorm") 

fit2_15_19 <- fitdist(Fish15_19_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Fish15_19_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_15_19$estimate[1],fit2_15_19$estimate[2]), col="red", lty=1)
gofstat(fit2_15_19, fitnames="gamma")


cvm.test(t, plnorm, fit1_15_19$estimate[1], fit1_15_19$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_15_19$estimate[1], fit1_15_19$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_15_19$estimate[1], fit2_15_19$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_15_19$estimate[1], fit2_15_19$estimate[2])  #Anderson-Darling test


#Age 20-24
Fish20_24 <- subset(FishDaily, agegroups == "20-24", select = total.fish)

#Probability of fish intake
probfish20_24 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Fish20_24==0)/length(Fish20_24$total.fish) #probability of zero DHA exp
pr.yes <- sum(Fish20_24!=0)/length(Fish20_24$total.fish) #probability of DHA exp
probfish20_24[1,] <- c(pr.no,pr.yes)

Fish20_24_pos <- Fish20_24[which(Fish20_24!=0)] #Only intakes that are greater than zero
Fish20_24_pos <- as.vector(Fish20_24_pos$total.fish)

fit1_20_24 <- fitdist(Fish20_24_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Fish20_24_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_20_24$estimate[1], sdlog=fit1_20_24$estimate[2]), col="red", lty=1)
gofstat(fit1_20_24, fitnames="lnorm") 

fit2_20_24 <- fitdist(Fish20_24_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Fish20_24_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_20_24$estimate[1],fit2_20_24$estimate[2]), col="red", lty=1)
gofstat(fit2_20_24, fitnames="gamma")

cvm.test(t, plnorm, fit1_20_24$estimate[1], fit1_20_24$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_20_24$estimate[1], fit1_20_24$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_20_24$estimate[1], fit2_20_24$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_20_24$estimate[1], fit2_20_24$estimate[2])  #Anderson-Darling test


#Age 25-29
Fish25_29 <- subset(FishDaily, agegroups == "25-29", select = total.fish)

#Probability of fish intake
probfish25_29 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Fish25_29==0)/length(Fish25_29$total.fish) #probability of zero DHA exp
pr.yes <- sum(Fish25_29!=0)/length(Fish25_29$total.fish) #probability of DHA exp
probfish25_29[1,] <- c(pr.no,pr.yes)

Fish25_29_pos <- Fish25_29[which(Fish25_29!=0)] #Only intakes that are greater than zero
Fish25_29_pos <- as.vector(Fish25_29_pos$total.fish)

fit1_25_29 <- fitdist(Fish25_29_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Fish25_29_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_25_29$estimate[1], sdlog=fit1_25_29$estimate[2]), col="red", lty=1)
gofstat(fit1_25_29, fitnames="lnorm") 

fit2_25_29 <- fitdist(Fish25_29_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Fish25_29_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_25_29$estimate[1],fit2_25_29$estimate[2]), col="red", lty=1)
gofstat(fit2_25_29, fitnames="gamma")

cvm.test(t, plnorm, fit1_25_29$estimate[1], fit1_25_29$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_25_29$estimate[1], fit1_25_29$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_25_29$estimate[1], fit2_25_29$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_25_29$estimate[1], fit2_25_29$estimate[2])  #Anderson-Darling test



#Age 30-34
Fish30_34 <- subset(FishDaily, agegroups == "30-34", select = total.fish)

#Probability of fish intake
probfish30_34 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Fish30_34==0)/length(Fish30_34$total.fish) #probability of zero DHA exp
pr.yes <- sum(Fish30_34!=0)/length(Fish30_34$total.fish) #probability of DHA exp
probfish30_34[1,] <- c(pr.no,pr.yes)

Fish30_34_pos <- Fish30_34[which(Fish30_34!=0)] #Only intakes that are greater than zero
Fish30_34_pos <- as.vector(Fish30_34_pos$total.fish)

fit1_30_34 <- fitdist(Fish30_34_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Fish30_34_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_30_34$estimate[1], sdlog=fit1_30_34$estimate[2]), col="red", lty=1)
gofstat(fit1_30_34, fitnames="lnorm") 

fit2_30_34 <- fitdist(Fish30_34_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Fish30_34_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_30_34$estimate[1],fit2_30_34$estimate[2]), col="red", lty=1)
gofstat(fit2_30_34, fitnames="gamma")

cvm.test(t, plnorm, fit1_30_34$estimate[1], fit1_30_34$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_30_34$estimate[1], fit1_30_34$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_30_34$estimate[1], fit2_30_34$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_30_34$estimate[1], fit2_30_34$estimate[2])  #Anderson-Darling test


#Age 35-39
Fish35_39 <- subset(FishDaily, agegroups == "35-39", select = total.fish)

#Probability of fish intake
probfish35_39 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Fish35_39==0)/length(Fish35_39$total.fish) #probability of zero DHA exp
pr.yes <- sum(Fish35_39!=0)/length(Fish35_39$total.fish) #probability of DHA exp
probfish35_39[1,] <- c(pr.no,pr.yes)

Fish35_39_pos <- Fish35_39[which(Fish35_39!=0)] #Only intakes that are greater than zero
Fish35_39_pos <- as.vector(Fish35_39_pos$total.fish)

fit1_35_39 <- fitdist(Fish35_39_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Fish35_39_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_35_39$estimate[1], sdlog=fit1_35_39$estimate[2]), col="red", lty=1)
gofstat(fit1_35_39, fitnames="lnorm") 

fit2_35_39 <- fitdist(Fish35_39_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Fish35_39_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_35_39$estimate[1],fit2_35_39$estimate[2]), col="red", lty=1)
gofstat(fit2_35_39, fitnames="gamma")

cvm.test(t, plnorm, fit1_35_39$estimate[1], fit1_35_39$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_35_39$estimate[1], fit1_35_39$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_35_39$estimate[1], fit2_35_39$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_35_39$estimate[1], fit2_35_39$estimate[2])  #Anderson-Darling test



#Age 40-44
Fish40_44 <- subset(FishDaily, agegroups == "40-44", select = total.fish)

#Probability of fish intake
probfish40_44 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Fish40_44==0)/length(Fish40_44$total.fish) #probability of zero DHA exp
pr.yes <- sum(Fish40_44!=0)/length(Fish40_44$total.fish) #probability of DHA exp
probfish40_44[1,] <- c(pr.no,pr.yes)

Fish40_44_pos <- Fish40_44[which(Fish40_44!=0)] #Only intakes that are greater than zero
Fish40_44_pos <- as.vector(Fish40_44_pos$total.fish)

fit1_40_44 <- fitdist(Fish40_44_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Fish40_44_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_40_44$estimate[1], sdlog=fit1_40_44$estimate[2]), col="red", lty=1)
gofstat(fit1_40_44, fitnames="lnorm") 

fit2_40_44 <- fitdist(Fish40_44_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Fish40_44_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_40_44$estimate[1],fit2_40_44$estimate[2]), col="red", lty=1)
gofstat(fit2_40_44, fitnames="gamma")

cvm.test(t, plnorm, fit1_40_44$estimate[1], fit1_40_44$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_40_44$estimate[1], fit1_40_44$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_40_44$estimate[1], fit2_40_44$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_40_44$estimate[1], fit2_40_44$estimate[2])  #Anderson-Darling test


#Age 45-49
Fish45_49 <- subset(FishDaily, agegroups == "45-49", select = total.fish)

#Probability of fish intake
probfish45_49 <- matrix(nrow=1,ncol=2) #create empty matrix
pr.no <- sum(Fish45_49==0)/length(Fish45_49$total.fish) #probability of zero DHA exp
pr.yes <- sum(Fish45_49!=0)/length(Fish45_49$total.fish) #probability of DHA exp
probfish45_49[1,] <- c(pr.no,pr.yes)

Fish45_49_pos <- Fish45_49[which(Fish45_49!=0)] #Only intakes that are greater than zero
Fish45_49_pos <- as.vector(Fish45_49_pos$total.fish)

fit1_45_49 <- fitdist(Fish45_49_pos, 'lnorm') #Fit lognormal distribution to intake amounts
t <- Fish45_49_pos
plot(ecdf(t), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, plnorm(x, meanlog=fit1_45_49$estimate[1], sdlog=fit1_45_49$estimate[2]), col="red", lty=1)
gofstat(fit1_45_49, fitnames="lnorm") 

fit2_45_49 <- fitdist(Fish45_49_pos, 'gamma') #Fit gamma distribution to intake amounts
t2 <- Fish45_49_pos
plot(ecdf(t2), lty=1)
x <- seq(0, 0.2, length=1000)
lines(x, pgamma(x, fit2_45_49$estimate[1],fit2_45_49$estimate[2]), col="red", lty=1)
gofstat(fit2_45_49, fitnames="gamma")

cvm.test(t, plnorm, fit1_45_49$estimate[1], fit1_45_49$estimate[2])  #Cram?r-von Mises test
ad.test(t, plnorm, fit1_45_49$estimate[1], fit1_45_49$estimate[2])  #Anderson-Darling test

cvm.test(t2, pgamma, fit2_45_49$estimate[1], fit2_45_49$estimate[2])  #Cram?r-von Mises test
ad.test(t2, pgamma, fit2_45_49$estimate[1], fit2_45_49$estimate[2])  #Anderson-Darling test


##### DALY calcalculation #####

#Gamma distribution has the best fit for exposure
#Truncate exposure distribution at 30.5 g fish/day due to upper limit of dose-response

### Age 15-19 ###
set.seed(1)
ifexp <- sample(c(0,1), nvar, prob = probfish15_19, replace = T)

exp <- rgamma(nvar, fit2_15_19$estimate[1], fit2_15_19$estimate[2])
exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}


str(ifexp)
str(exp_trunc)
str(IQnorm)
str(r)


str(t(r)) #1 row, 1000 columns - columns: uncertainty

## IQ change: IQdiff = ifexp * exp * r
## .. each col = uncertainty simulation
## .. each row = variability simulation
## .. IQchange[var, unc]

IQdiff15_19 <- apply(t(r), 2, function(x) x * ifexp * exp_trunc)
str(IQdiff15_19)

IQnew <- apply(t(r), 2, function(x) x * ifexp * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal15_19 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline15_19 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild15_19 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate15_19 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere15_19 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound15_19 <- apply(IQnew, 2, function(x) mean(x < 20))

sum(mean(IQnormal15_19)+mean(IQborderline15_19)+mean(IQmild15_19)+mean(IQmoderate15_19)+mean(IQsevere15_19)+mean(IQprofound15_19))


DALY15_19 <- Pbaby15_19*(IQnormal15_19*wIQ85 + IQborderline15_19*wIQ70.85sim + IQmild15_19*wIQ50.69sim +
                             IQmoderate15_19*wIQ35.49sim + IQsevere15_19*wIQ20.34sim +
                             IQprofound15_19*wIQ20sim) * LE * women15_19 #per total women in DK in this age group

mean_median_ci(DALY15_19)



### Age 20-24 ###
set.seed(1)
ifexp <- sample(c(0,1), nvar, prob = probfish20_24, replace = T)

exp <- rgamma(nvar, fit2_20_24$estimate[1], fit2_20_24$estimate[2])
exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}


IQdiff20_24 <- apply(t(r), 2, function(x) x * ifexp * exp_trunc)
str(IQdiff20_24)

IQnew <- apply(t(r), 2, function(x) x * ifexp * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal20_24 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline20_24 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild20_24 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate20_24 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere20_24 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound20_24 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY20_24 <- Pbaby20_24*(IQnormal20_24*wIQ85 + IQborderline20_24*wIQ70.85sim + IQmild20_24*wIQ50.69sim +
                             IQmoderate20_24*wIQ35.49sim + IQsevere20_24*wIQ20.34sim +
                             IQprofound20_24*wIQ20sim) * LE * women20_24 #per total women in DK in this age group

mean_median_ci(DALY20_24)


### Age 25-29 ###
set.seed(1)
ifexp <- sample(c(0,1), nvar, prob = probfish25_29, replace = T)

exp <- rgamma(nvar, fit2_25_29$estimate[1], fit2_25_29$estimate[2])
exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}



IQdiff25_29 <- apply(t(r), 2, function(x) x * ifexp * exp_trunc)
str(IQdiff25_29)

IQnew <- apply(t(r), 2, function(x) x * ifexp * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal25_29 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline25_29 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild25_29 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate25_29 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere25_29 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound25_29 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY25_29 <- Pbaby25_29*(IQnormal25_29*wIQ85 + IQborderline25_29*wIQ70.85sim + IQmild25_29*wIQ50.69sim +
                             IQmoderate25_29*wIQ35.49sim + IQsevere25_29*wIQ20.34sim +
                             IQprofound25_29*wIQ20sim) * LE * women25_29 #per total women in DK in this age group

mean_median_ci(DALY25_29)


### Age 30-34 ###
set.seed(1)
ifexp <- sample(c(0,1), nvar, prob = probfish30_34, replace = T)

exp <- rgamma(nvar, fit2_30_34$estimate[1], fit2_30_34$estimate[2])
exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}


IQdiff30_34 <- apply(t(r), 2, function(x) x * ifexp * exp_trunc)
str(IQdiff30_34)

IQnew <- apply(t(r), 2, function(x) x * ifexp * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal30_34 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline30_34 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild30_34 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate30_34 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere30_34 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound30_34 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY30_34 <- Pbaby30_34*(IQnormal30_34*wIQ85 + IQborderline30_34*wIQ70.85sim + IQmild30_34*wIQ50.69sim +
                             IQmoderate30_34*wIQ35.49sim + IQsevere30_34*wIQ20.34sim +
                             IQprofound30_34*wIQ20sim) * LE * women30_34 #per total women in DK in this age group

mean_median_ci(DALY30_34)



### Age 35-39 ###
set.seed(1)
ifexp <- sample(c(0,1), nvar, prob = probfish35_39, replace = T)

exp <- rgamma(nvar, fit2_35_39$estimate[1], fit2_35_39$estimate[2])
exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}


IQdiff35_39 <- apply(t(r), 2, function(x) x * ifexp * exp_trunc)
str(IQdiff35_39)

IQnew <- apply(t(r), 2, function(x) x * ifexp * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal35_39 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline35_39 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild35_39 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate35_39 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere35_39 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound35_39 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY35_39 <- Pbaby35_39*(IQnormal35_39*wIQ85 + IQborderline35_39*wIQ70.85sim + IQmild35_39*wIQ50.69sim +
                             IQmoderate35_39*wIQ35.49sim + IQsevere35_39*wIQ20.34sim +
                             IQprofound35_39*wIQ20sim) * LE * women35_39 #per total women in DK in this age group

mean_median_ci(DALY35_39)


### Age 40-44 ###
set.seed(1)
ifexp <- sample(c(0,1), nvar, prob = probfish40_44, replace = T)

exp <- rgamma(nvar, fit2_40_44$estimate[1], fit2_40_44$estimate[2])
exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}


IQdiff40_44 <- apply(t(r), 2, function(x) x * ifexp * exp_trunc)
str(IQdiff40_44)

IQnew <- apply(t(r), 2, function(x) x * ifexp * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal40_44 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline40_44 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild40_44 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate40_44 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere40_44 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound40_44 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY40_44 <- Pbaby40_44*(IQnormal40_44*wIQ85 + IQborderline40_44*wIQ70.85sim + IQmild40_44*wIQ50.69sim +
                             IQmoderate40_44*wIQ35.49sim + IQsevere40_44*wIQ20.34sim +
                             IQprofound40_44*wIQ20sim) * LE * women40_44 #per total women in DK in this age group

mean_median_ci(DALY40_44)


### Age 45-49 ###
set.seed(1)
ifexp <- sample(c(0,1), nvar, prob = probfish45_49, replace = T)

exp <- rgamma(nvar, fit2_45_49$estimate[1], fit2_45_49$estimate[2])
exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}


IQdiff45_49 <- apply(t(r), 2, function(x) x * ifexp * exp_trunc)
str(IQdiff45_49)

IQnew <- apply(t(r), 2, function(x) x * ifexp * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal45_49 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline45_49 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild45_49 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate45_49 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere45_49 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound45_49 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY45_49 <- Pbaby45_49*(IQnormal45_49*wIQ85 + IQborderline45_49*wIQ70.85sim + IQmild45_49*wIQ50.69sim +
                             IQmoderate45_49*wIQ35.49sim + IQsevere45_49*wIQ20.34sim +
                             IQprofound45_49*wIQ20sim) * LE * women45_49 #per total women in DK in this age group

mean_median_ci(DALY45_49)


##### Total DALYs #####

tDALY_ref_fish_IQ <- DALY15_19 + DALY20_24 + DALY25_29 + DALY30_34 + DALY35_39 + DALY40_44 + DALY45_49
mean_median_ci(tDALY_ref_fish_IQ)

IQchange.fishIQ <- (IQdiff15_19 + IQdiff20_24 + IQdiff25_29 + IQdiff30_34 + IQdiff35_39 + IQdiff40_44 + IQdiff45_49)/7


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


#Summing to total number of cases
totalcases_ref <- cases_borderline + cases_mild + cases_moderate + cases_severe + cases_profound
mean_median_ci(totalcases_ref)


##### Scenario 1 #####

Alt_scenario <- read.csv("Scenario1.csv") # Since we're considering fish as a whole, the results for this endpoint will be the same
# for scenario 1,2,3,4

# fish intakes are raw weights and need to be converted to prepared weights
#Exclude males, and women above 49y of age (due to data on female fertility; women below 15y have already been excluded in FishDaily dataset)

Alt_scenario <- subset(Alt_scenario, sex == 2 & age < 50)

Alt_scenario <- Alt_scenario[,c(1:5)] #We choose the prepared fish weights since dose-response is for prepared weight

setDT(Alt_scenario)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]


##### Fitting exposures scen1 #####

#Divide data into agegroups. Use empirical distribution

#Age 15-19
Fish15_19 <- subset(Alt_scenario, agegroups=="15-19", select = total.fish.new)
Fish15_19 <- as.vector(Fish15_19$total.fish.new)

#We don't need to model the probability of MeHg exposure since all individuals now consume fish

#Age 20-24
Fish20_24 <- subset(Alt_scenario, agegroups=="20-24", select = total.fish.new)
Fish20_24 <- as.vector(Fish20_24$total.fish.new)

#Age 25-29
Fish25_29 <- subset(Alt_scenario, agegroups=="25-29", select = total.fish.new)
Fish25_29 <- as.vector(Fish25_29$total.fish.new)

#Age 30-34
Fish30_34 <- subset(Alt_scenario, agegroups=="30-34", select = total.fish.new)
Fish30_34 <- as.vector(Fish30_34$total.fish.new)

#Age 35-39
Fish35_39 <- subset(Alt_scenario, agegroups=="35-39", select = total.fish.new)
Fish35_39 <- as.vector(Fish35_39$total.fish.new)

#Age 40-44
Fish40_44 <- subset(Alt_scenario, agegroups=="40-44", select = total.fish.new)
Fish40_44 <- as.vector(Fish40_44$total.fish.new)

#Age 45-49
Fish45_49 <- subset(Alt_scenario, agegroups=="45-49", select = total.fish.new)
Fish45_49 <- as.vector(Fish45_49$total.fish.new)



##### DALY calculation #####

#Age 15-19
set.seed(1)

exp <- sample(Fish15_19, nvar, replace = T) #Everyone are exposed so ifexp = 1

exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}


IQdiff15_19 <- apply(t(r), 2, function(x) x * exp_trunc)
str(IQdiff15_19)

IQnew <- apply(t(r), 2, function(x) x * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal15_19 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline15_19 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild15_19 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate15_19 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere15_19 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound15_19 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY15_19 <- Pbaby15_19*(IQnormal15_19*wIQ85 + IQborderline15_19*wIQ70.85sim + IQmild15_19*wIQ50.69sim +
                             IQmoderate15_19*wIQ35.49sim + IQsevere15_19*wIQ20.34sim +
                             IQprofound15_19*wIQ20sim) * LE * women15_19 #per total women in DK in this age group

mean_median_ci(DALY15_19)


#Age 20-24

set.seed(1)

exp <- sample(Fish20_24, nvar, replace = T) #Everyone are exposed so ifexp = 1

exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}



IQdiff20_24 <- apply(t(r), 2, function(x) x * exp_trunc)
str(IQdiff20_24)

IQnew <- apply(t(r), 2, function(x) x * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal20_24 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline20_24 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild20_24 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate20_24 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere20_24 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound20_24 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY20_24 <- Pbaby20_24*(IQnormal20_24*wIQ85 + IQborderline20_24*wIQ70.85sim + IQmild20_24*wIQ50.69sim +
                             IQmoderate20_24*wIQ35.49sim + IQsevere20_24*wIQ20.34sim +
                             IQprofound20_24*wIQ20sim) * LE * women20_24 #per total women in DK in this age group

mean_median_ci(DALY20_24)


#Age 25-29

set.seed(1)

exp <- sample(Fish25_29, nvar, replace = T) #Everyone are exposed so ifexp = 1

exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}


IQdiff25_29 <- apply(t(r), 2, function(x) x * exp_trunc)
str(IQdiff25_29)

IQnew <- apply(t(r), 2, function(x) x * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal25_29 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline25_29 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild25_29 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate25_29 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere25_29 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound25_29 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY25_29 <- Pbaby25_29*(IQnormal25_29*wIQ85 + IQborderline25_29*wIQ70.85sim + IQmild25_29*wIQ50.69sim +
                             IQmoderate25_29*wIQ35.49sim + IQsevere25_29*wIQ20.34sim +
                             IQprofound25_29*wIQ20sim) * LE * women25_29 #per total women in DK in this age group

mean_median_ci(DALY25_29)


#Age 30-34

set.seed(1)

exp <- sample(Fish30_34, nvar, replace = T) #Everyone are exposed so ifexp = 1

exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}


IQdiff30_34 <- apply(t(r), 2, function(x) x * exp_trunc)
str(IQdiff30_34)

IQnew <- apply(t(r), 2, function(x) x * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal30_34 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline30_34 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild30_34 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate30_34 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere30_34 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound30_34 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY30_34 <- Pbaby30_34*(IQnormal30_34*wIQ85 + IQborderline30_34*wIQ70.85sim + IQmild30_34*wIQ50.69sim +
                             IQmoderate30_34*wIQ35.49sim + IQsevere30_34*wIQ20.34sim +
                             IQprofound30_34*wIQ20sim) * LE * women30_34 #per total women in DK in this age group

mean_median_ci(DALY30_34)



#Age 35-39

set.seed(1)

exp <- sample(Fish35_39, nvar, replace = T) #Everyone are exposed so ifexp = 1

exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}



IQdiff35_39 <- apply(t(r), 2, function(x) x * exp_trunc)
str(IQdiff35_39)

IQnew <- apply(t(r), 2, function(x) x * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal35_39 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline35_39 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild35_39 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate35_39 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere35_39 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound35_39 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY35_39 <- Pbaby35_39*(IQnormal35_39*wIQ85 + IQborderline35_39*wIQ70.85sim + IQmild35_39*wIQ50.69sim +
                           IQmoderate35_39*wIQ35.49sim + IQsevere35_39*wIQ20.34sim +
                           IQprofound35_39*wIQ20sim) * LE * women35_39 #per total women in DK in this age group

mean_median_ci(DALY35_39)


#Age 40-44

set.seed(1)

exp <- sample(Fish40_44, nvar, replace = T) #Everyone are exposed so ifexp = 1

exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}


IQdiff40_44 <- apply(t(r), 2, function(x) x * exp_trunc)
str(IQdiff40_44)

IQnew <- apply(t(r), 2, function(x) x * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal40_44 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline40_44 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild40_44 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate40_44 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere40_44 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound40_44 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY40_44 <- Pbaby40_44*(IQnormal40_44*wIQ85 + IQborderline40_44*wIQ70.85sim + IQmild40_44*wIQ50.69sim +
                           IQmoderate40_44*wIQ35.49sim + IQsevere40_44*wIQ20.34sim +
                           IQprofound40_44*wIQ20sim) * LE * women40_44 #per total women in DK in this age group

mean_median_ci(DALY40_44)



#Age 45-49

set.seed(1)

exp <- sample(Fish45_49, nvar, replace = T) #Everyone are exposed so ifexp = 1

exp_trunc <- numeric(nvar)
for(i in 1:nvar){
  x <- sample(exp, 1, replace = F)
  y <- ifelse(x<30.5, x,30.5)
  exp_trunc[i] <- y
}


IQdiff45_49 <- apply(t(r), 2, function(x) x * exp_trunc)
str(IQdiff45_49)

IQnew <- apply(t(r), 2, function(x) x * exp_trunc + IQnorm) #multiply each col (r, uncertainty) by prob. of exp, pos. exp and add IQ dist
str(IQnew)


IQnormal45_49 <- apply(IQnew, 2, function(x) mean(x >= 85))
IQborderline45_49 <- apply(IQnew, 2, function(x) mean(x >= 70 & x < 85))
IQmild45_49 <- apply(IQnew, 2, function(x) mean(x >= 50 & x < 70))
IQmoderate45_49 <- apply(IQnew, 2, function(x) mean(x >= 35 & x < 50))
IQsevere45_49 <- apply(IQnew, 2, function(x) mean(x >= 20 & x < 35))
IQprofound45_49 <- apply(IQnew, 2, function(x) mean(x < 20))


DALY45_49 <- Pbaby45_49*(IQnormal45_49*wIQ85 + IQborderline45_49*wIQ70.85sim + IQmild45_49*wIQ50.69sim +
                           IQmoderate45_49*wIQ35.49sim + IQsevere45_49*wIQ20.34sim +
                           IQprofound45_49*wIQ20sim) * LE * women45_49 #per total women in DK in this age group

mean_median_ci(DALY45_49)


##### Total DALYs alt scen #####

tDALY_alt_fish_IQ <- DALY15_19 + DALY20_24 + DALY25_29 + DALY30_34 + DALY35_39 + DALY40_44 + DALY45_49

mean_median_ci(tDALY_alt_fish_IQ)


IQchange.fishIQ.scen1 <- (IQdiff15_19 + IQdiff20_24 + IQdiff25_29 + IQdiff30_34 + IQdiff35_39 + IQdiff40_44 + IQdiff45_49)/7


##### Delta DALY #####

dDALY_fish_IQ <- tDALY_alt_fish_IQ - tDALY_ref_fish_IQ

mean_median_ci(dDALY_fish_IQ)



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
totalcases_alt <- cases_borderline + cases_mild + cases_moderate + cases_severe + cases_profound
mean_median_ci(totalcases_alt)


##### Extra number of cases #####

cases_diff_alt <- totalcases_alt - totalcases_ref
mean_median_ci(cases_diff_alt)



##### Save #####

write.csv(tDALY_ref_fish_IQ, "tDALY_ref_fish_IQ.csv")

write.csv(tDALY_alt_fish_IQ, "tDALY_alt_fish_IQ.csv")

IQchange.fishIQ.scen1_unc <- apply(IQchange.fishIQ.scen1, 2, function(x) x*1) #get uncertainty dimension
write.csv(IQchange.fishIQ.scen1_unc, "IQchange.fishIQ.scen1.csv")
