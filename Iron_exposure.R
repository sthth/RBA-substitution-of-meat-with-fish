
### We use a simple fraction calculation to estimate the fraction of individuals below AR


#Packages
library(data.table)
library(plyr)
library(fitdistrplus)
library(goftest)


##### Reference scenario #####

##### Exposure from fish #####

FishDaily <- read.csv("FishDaily.csv")
#NB! Water loss has NOT been taken into account in FoodDaily/FishDaily dataset so amounts are prepared amounts.
#Water loss HAS been taken into account in the substitution Scenarios


#Convert fish consumption amounts into raw weights - assume that all fish is prepared and have a water loss of 20%
wloss.fish <- 1/0.8
FishDaily[,8:23] <- FishDaily[,8:23]*wloss.fish
FishDaily$Totalfish.raw <- rowSums(FishDaily[,8:23]) #Test if it gives the same as total.fish * 1.25


Nutrients_fish <- read.csv2("NutrientData_fish.csv") 

rownames(Nutrients_fish) <- Nutrients_fish[,1]
Nutrients_fish[,1] <- NULL

Iron_fish <- Nutrients_fish[,c(1)]

Iron_exp_fish <- t(FishDaily[,c(8:23)]) #Create data set only with intakes of the different fish species
Iron_exp_fish <- Iron_exp_fish * Iron_fish


Iron_exp <- FishDaily[,1:4]

agebreaks <- c(15,20,25,30,35,40,45,50,55,60,65,70,75,80)
agelabels= c("15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79")

setDT(Iron_exp)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

Iron_exp$Iron_fish <- colSums(Iron_exp_fish) #Add sum of exposures for each individual to the FishDaily dataset



##### Exposure from red and processed meat #####

FoodDaily <- read.csv("FoodDaily.csv")

#NB! Water loss has NOT been taken into account in FoodDaily/FishDaily dataset so amounts are prepared amounts.
#Water loss HAS been taken into account in the substitution Scenarios

#Convert red meat amounts into raw weights - 25% water loss. Processed meat do not need to be converted because we have
#fat contents for individual processed meats (and also concentration data)
wloss.redmeat <- 1/0.75
FoodDaily[,46:57] <- FoodDaily[,46:57]*wloss.redmeat


Nutrients_meat <- read.csv2("NutrientData_meat.csv") 

rownames(Nutrients_meat) <- Nutrients_meat[,1]
Nutrients_meat[,1] <- NULL

Iron_meat <- Nutrients_meat[,c(1)]

Iron_exp_meat <- t(FoodDaily[,c(35:44,46:57)]) #Create data set only with intakes of the different fish species
Iron_exp_meat <- Iron_exp_meat * Iron_meat


Iron_exp$Iron_meat <- colSums(Iron_exp_meat) #Add sum of exposures for each individual to the FishDaily dataset


##### Total exp from fish and meat #####

Iron_exp$Iron_total <- Iron_exp$Iron_fish + Iron_exp$Iron_meat #Only total iron from fish and meat


##### Scenario 1 #####

Scenario1 <- read.csv("Scenario1.csv")
# fish intakes are raw weights, ready for exposure assessment.



##### Exposure from fish #####

scen1.Iron_fish <- t(Scenario1[,c(55:70)]) #Create data set only with intakes of the different fish species

scen1.Iron_fish <- scen1.Iron_fish * Iron_fish ##Multiply intakes of fish species with concentrations

Iron_exp_scen1 <- Scenario1[,c(1:4, 54)]

setDT(Iron_exp_scen1)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

Iron_exp_scen1$Iron_fish <- colSums(scen1.Iron_fish) #Add sum of exposures for each individual to the Scenario1 dataset



##### Exposure from red and processed meat #####

#Meat intake is the same in all 4 alternative scenarios

Alt.Iron_meat <- t(Scenario1[,c(31:40,42:53)])

Alt.Iron_meat <- Alt.Iron_meat * Iron_meat

Iron_exp_scen1$Iron_meat <- colSums(Alt.Iron_meat) #Add sum of exposures for each individual to the Scenario1 dataset



##### Total exp from fish and meat #####

Iron_exp_scen1$Iron_total <- Iron_exp_scen1$Iron_fish + Iron_exp_scen1$Iron_meat #Only total iron from fish and meat



##### Change in iron intake #####

Iron_exp_scen1$Iron_diff <- Iron_exp_scen1$Iron_total - Iron_exp$Iron_total

summary(Iron_exp_scen1$Iron_diff)



Nutrient_status <- read.csv2("NutrientDataTotal_IP.csv") #Current status of micronutrients
Nutrient_status[is.na(Nutrient_status)] <- 0


Mean_status <- aggregate(Nutrient_status, by=list(Nutrient_status$PROJEKTID), 
                         FUN=mean)


FoodDaily <- read.csv("FoodDaily.csv")
FoodDaily_IP <- as.data.frame(FoodDaily[,1])
names(FoodDaily_IP)[1] <- 'PROJEKTID'

Status_IP_mean <- merge(Mean_status,FoodDaily_IP,by="PROJEKTID") #only the 2811 individuals we are considering (>15y, reported all 7 days)

Status_Iron <- cbind(Status_IP_mean[,c(1,8)])#X8 is iron

Iron_exp_scen1 <- merge(Iron_exp_scen1, Status_Iron, by="PROJEKTID")

names(Iron_exp_scen1)[11] <- "Status"

summary(Iron_exp_scen1$Status) #Current iron intake


#New iron status
Iron_exp_scen1$Status_new <- Iron_exp_scen1$Status + Iron_exp_scen1$Iron_diff

summary(Iron_exp_scen1$Status_new)


Iron_male_scen1 <- subset(Iron_exp_scen1, sex=="1", select = c(age,agegroups,Status,Status_new))


Iron_female_scen1 <- subset(Iron_exp_scen1, sex=="2", select = c(age,agegroups,Status,Status_new))



##### Fraction < RI reference #####

#NNR 2012: Recommended intake (RI) of iron is 7 mg/day for men and 6 mg/day for post-menupausal women. RI is 10 mg/day for other women.
#Divide data into agegroups and fit exposure to distribution

#Males
mean(Iron_male_scen1$Status)


sd(Iron_male_scen1$Status)


median(Iron_male_scen1$Status)


mean(Iron_male_scen1$Status < 7)



#Females > 49 (postmenupausal)
Iron_pmfemale_scen1 <- subset(Iron_female_scen1, age > 49)

mean(Iron_pmfemale_scen1$Status)


sd(Iron_pmfemale_scen1$Status)


median(Iron_pmfemale_scen1$Status)


mean(Iron_pmfemale_scen1$Status < 6)



#Females <= 49 (fertile, non-postmenupausal)

Iron_ffemale_scen1 <- subset(Iron_female_scen1, age < 50)

mean(Iron_ffemale_scen1$Status)


sd(Iron_ffemale_scen1$Status)


median(Iron_ffemale_scen1$Status)


mean(Iron_ffemale_scen1$Status < 10)




##### Fraction < RI scenario 1 #####

#Males

mean(Iron_male_scen1$Status_new)


sd(Iron_male_scen1$Status_new)


median(Iron_male_scen1$Status_new)


mean(Iron_male_scen1$Status_new < 7)



#Females > 49 (postmenupausal)

mean(Iron_pmfemale_scen1$Status_new)


sd(Iron_pmfemale_scen1$Status_new)


median(Iron_pmfemale_scen1$Status_new)


mean(Iron_pmfemale_scen1$Status_new < 6)



#Females <= 49 (fertile, non-postmenupausal)

mean(Iron_ffemale_scen1$Status_new)


sd(Iron_ffemale_scen1$Status_new)


median(Iron_ffemale_scen1$Status_new)


mean(Iron_ffemale_scen1$Status_new < 10)



##### Scenario 2 #####

Scenario2 <- read.csv("Scenario2.csv")
# fish intakes are raw weights, ready for exposure assessment.


##### Exposure from fish #####

scen2.Iron_fish <- t(Scenario2[,c(55:60)]) #Create data set only with intakes of the different fish species

scen2.Iron_fish <- scen2.Iron_fish * Iron_fish[c(3,9,10,11,12,16)] ##Multiply intakes of fish species with concentrations

Iron_exp_scen2 <- Scenario2[,c(1:4, 54)]

setDT(Iron_exp_scen2)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

Iron_exp_scen2$Iron_fish <- colSums(scen2.Iron_fish) #Add sum of exposures for each individual to the Scenario2 dataset


##### Exposure from red and processed meat #####

#Meat intake is the same in all 4 alternative scenarios

Iron_exp_scen2$Iron_meat <- colSums(Alt.Iron_meat) #Add sum of exposures for each individual to the Scenario2 dataset


##### Total exp from fish and meat #####


Iron_exp_scen2$Iron_total <- Iron_exp_scen2$Iron_fish + Iron_exp_scen2$Iron_meat #Only total iron from fish and meat


##### Change in iron intake #####

Iron_exp_scen2$Iron_diff <- Iron_exp_scen2$Iron_total - Iron_exp$Iron_total

summary(Iron_exp_scen2$Iron_diff)


Iron_exp_scen2 <- merge(Iron_exp_scen2, Status_Iron, by="PROJEKTID")

names(Iron_exp_scen2)[11] <- "Status"


#New iron status
Iron_exp_scen2$Status_new <- Iron_exp_scen2$Status + Iron_exp_scen2$Iron_diff

summary(Iron_exp_scen2$Status_new)



Iron_male_scen2 <- subset(Iron_exp_scen2, sex=="1", select = c(age,agegroups,Status,Status_new))


Iron_female_scen2 <- subset(Iron_exp_scen2, sex=="2", select = c(age,agegroups,Status,Status_new))



##### Fraction < RI scenario 2 #####

#Males

mean(Iron_male_scen2$Status_new)

sd(Iron_male_scen2$Status_new)

median(Iron_male_scen2$Status_new)

mean(Iron_male_scen2$Status_new < 7)


#Females > 49 (postmenupausal)

Iron_pmfemale_scen2 <- subset(Iron_female_scen2, age > 49)

mean(Iron_pmfemale_scen2$Status_new)

sd(Iron_pmfemale_scen2$Status_new)

median(Iron_pmfemale_scen2$Status_new)

mean(Iron_pmfemale_scen2$Status_new < 6)


#Females <= 49 (fertile, non-postmenupausal)

Iron_ffemale_scen2 <- subset(Iron_female_scen2, age < 50)

mean(Iron_ffemale_scen2$Status_new)

sd(Iron_ffemale_scen2$Status_new)

median(Iron_ffemale_scen2$Status_new)

mean(Iron_ffemale_scen2$Status_new < 10)



##### Scenario 3 #####

Scenario3 <- read.csv("Scenario3.csv")
# fish intakes are raw weights, ready for exposure assessment.



##### Exposure from fish #####

scen3.Iron_fish <- t(Scenario3[,c(55:64)]) #Create data set only with intakes of the different fish species

scen3.Iron_fish <- scen3.Iron_fish * Iron_fish[c(1,2,4,5,6,7,8,13,14,15)] ##Multiply intakes of fish species with concentrations

Iron_exp_scen3 <- Scenario3[,c(1:4, 54)]

setDT(Iron_exp_scen3)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

Iron_exp_scen3$Iron_fish <- colSums(scen3.Iron_fish) #Add sum of exposures for each individual to the Scenario3 dataset



##### Exposure from red and processed meat #####

#Meat intake is the same in all 4 alternative scenarios

Iron_exp_scen3$Iron_meat <- colSums(Alt.Iron_meat) #Add sum of exposures for each individual to the Scenario3 dataset



##### Total exp from fish and meat #####

Iron_exp_scen3$Iron_total <- Iron_exp_scen3$Iron_fish + Iron_exp_scen3$Iron_meat #Only total iron from fish and meat



##### Change in iron intake #####

Iron_exp_scen3$Iron_diff <- Iron_exp_scen3$Iron_total - Iron_exp$Iron_total

summary(Iron_exp_scen3$Iron_diff)


Iron_exp_scen3 <- merge(Iron_exp_scen3, Status_Iron, by="PROJEKTID")

names(Iron_exp_scen3)[11] <- "Status"


#New iron status
Iron_exp_scen3$Status_new <- Iron_exp_scen3$Status + Iron_exp_scen3$Iron_diff

summary(Iron_exp_scen3$Status_new) 


Iron_male_scen3 <- subset(Iron_exp_scen3, sex=="1", select = c(age,agegroups,Status,Status_new))


Iron_female_scen3 <- subset(Iron_exp_scen3, sex=="2", select = c(age,agegroups,Status,Status_new))



##### Fraction < RI scenario 3 #####

#Males

mean(Iron_male_scen3$Status_new)

sd(Iron_male_scen3$Status_new)

median(Iron_male_scen3$Status_new)

mean(Iron_male_scen3$Status_new < 7)


#Females > 49 (postmenupausal)
Iron_pmfemale_scen3 <- subset(Iron_female_scen3, age > 49)

mean(Iron_pmfemale_scen3$Status_new)

sd(Iron_pmfemale_scen3$Status_new)

median(Iron_pmfemale_scen3$Status_new)

mean(Iron_pmfemale_scen3$Status_new < 6)


#Females <= 49 (fertile, non-postmenupausal)

Iron_ffemale_scen3 <- subset(Iron_female_scen3, age < 50)

mean(Iron_ffemale_scen3$Status_new)

sd(Iron_ffemale_scen3$Status_new)

median(Iron_ffemale_scen3$Status_new)

mean(Iron_ffemale_scen3$Status_new < 10)




##### Scenario 4 #####

Scenario4 <- read.csv("Scenario4.csv")
# fish intakes are raw weights, ready for exposure assessment.


##### Exposure from fish #####

scen4.Iron_fish <- t(Scenario4[,55]) #Create data set only with intakes of the different fish species

scen4.Iron_fish <- scen4.Iron_fish * Iron_fish[8] ##Multiply intakes of fish species with concentrations

Iron_exp_scen4 <- Scenario4[,c(1:4, 54)]

setDT(Iron_exp_scen4)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

Iron_exp_scen4$Iron_fish <- colSums(scen4.Iron_fish) #Add sum of exposures for each individual to the Scenario4 dataset



##### Exposure from red and processed meat #####

#Meat intake is the same in all 4 alternative scenarios

Iron_exp_scen4$Iron_meat <- colSums(Alt.Iron_meat) #Add sum of exposures for each individual to the Scenario4 dataset



##### Total exp from fish and meat #####


Iron_exp_scen4$Iron_total <- Iron_exp_scen4$Iron_fish + Iron_exp_scen4$Iron_meat #Only total iron from fish and meat



##### Change in iron intake #####

Iron_exp_scen4$Iron_diff <- Iron_exp_scen4$Iron_total - Iron_exp$Iron_total

summary(Iron_exp_scen4$Iron_diff)


Iron_exp_scen4 <- merge(Iron_exp_scen4, Status_Iron, by="PROJEKTID")

names(Iron_exp_scen4)[11] <- "Status"


#New iron status
Iron_exp_scen4$Status_new <- Iron_exp_scen4$Status + Iron_exp_scen4$Iron_diff

summary(Iron_exp_scen4$Status_new) 


Iron_male_scen4 <- subset(Iron_exp_scen4, sex=="1", select = c(age,agegroups,Status,Status_new))



Iron_female_scen4 <- subset(Iron_exp_scen4, sex=="2", select = c(age,agegroups,Status,Status_new))



##### Fraction < RI scenario 4 #####

#Males

mean(Iron_male_scen4$Status_new)

sd(Iron_male_scen4$Status_new)

median(Iron_male_scen4$Status_new)

mean(Iron_male_scen4$Status_new < 7)


#Females > 49 (postmenupausal)

Iron_pmfemale_scen4 <- subset(Iron_female_scen4, age > 49)

mean(Iron_pmfemale_scen4$Status_new)

sd(Iron_pmfemale_scen4$Status_new)

median(Iron_pmfemale_scen4$Status_new)

mean(Iron_pmfemale_scen4$Status_new < 6)


#Females <= 49 (fertile, non-postmenupausal)

Iron_ffemale_scen4 <- subset(Iron_female_scen4, age < 50)

mean(Iron_ffemale_scen4$Status_new)

sd(Iron_ffemale_scen4$Status_new)

median(Iron_ffemale_scen4$Status_new)

mean(Iron_ffemale_scen4$Status_new < 10)

