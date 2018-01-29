
### We use a simple fraction calculation to estimate the fraction of individuals below AR


#Packages
library(fitdistrplus)
library(data.table)
library(plyr)

##### Reference scenario #####

#### Exposure from fish ####

FishDaily <- read.csv("FishDaily.csv")
#NB! Water loss has NOT been taken into account in FoodDaily/FishDaily dataset so amounts are prepared amounts.
#Water loss HAS been taken into account in the substitution Scenarios


#Convert fish consumption amounts into raw weights - assume that all fish is prepared and have a water loss of 20%
wloss.fish <- 1/0.8
FishDaily[,8:23] <- FishDaily[,8:23]*wloss.fish
FishDaily$Totalfish.raw <- rowSums(FishDaily[,8:23])


Nutrients_fish <- read.csv2("NutrientData_fish.csv") 

rownames(Nutrients_fish) <- Nutrients_fish[,1]
Nutrients_fish[,1] <- NULL

VitD_fish <- Nutrients_fish[,c(2)]

VitD_exp_fish <- t(FishDaily[,c(8:23)]) #Create data set only with intakes of the different fish species
VitD_exp_fish <- VitD_exp_fish * VitD_fish

VitD_exp <- FishDaily[,1:4]

agebreaks <- c(15,20,25,30,35,40,45,50,55,60,65,70,75,80)
agelabels= c("15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79")

setDT(VitD_exp)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

VitD_exp$VitD_fish <- colSums(VitD_exp_fish) #Add sum of exposures for each individual to the FishDaily dataset


#### Exposure from red and processed meat ####

FoodDaily <- read.csv("FoodDaily.csv")

#NB! Water loss has NOT been taken into account in FoodDaily/FishDaily dataset so amounts are prepared amounts.
#Water loss HAS been taken into account in the substitution Scenarios

#Convert red meat amounts into raw weights - 25% water loss. Processed meat do not need to be converted because we have
#fat contents for individual processed meats (and also concentration data)
wloss.redmeat <- 1/0.75
FoodDaily[,46:57] <- FoodDaily[,46:57]*wloss.redmeat

mean(rowSums(FoodDaily[,46:57])) 
mean(FoodDaily$red.meat.raw)


Nutrients_meat <- read.csv2("NutrientData_meat.csv") 

rownames(Nutrients_meat) <- Nutrients_meat[,1]
Nutrients_meat[,1] <- NULL

VitD_meat <- Nutrients_meat[,c(2)]

VitD_exp_meat <- t(FoodDaily[,c(35:44,46:57)]) #Create data set only with intakes of the different fish species
VitD_exp_meat <- VitD_exp_meat * VitD_meat


VitD_exp$VitD_meat <- colSums(VitD_exp_meat) #Add sum of exposures for each individual to the FishDaily dataset


#### Total exp from fish and meat ####

VitD_exp$VitD_total <- VitD_exp$VitD_fish + VitD_exp$VitD_meat #Only total vit D from fish and meat



##### Scenario 1 #####

Scenario1 <- read.csv("Scenario1.csv")
# fish intakes are raw weights, ready for exposure assessment.


#### Exposure from fish ####

scen1.vitD_fish <- t(Scenario1[,c(55:70)]) #Create data set only with intakes of the different fish species

scen1.vitD_fish <- scen1.vitD_fish * VitD_fish ##Multiply intakes of fish species with concentrations

VitD_exp_scen1 <- Scenario1[,c(1:4, 54)]

setDT(VitD_exp_scen1)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

VitD_exp_scen1$VitD_fish <- colSums(scen1.vitD_fish) #Add sum of exposures for each individual to the Scenario1 dataset



#### Exposure from red and processed meat ####

#Meat intake is the same in all 4 alternative scenarios

Alt.vitD_meat <- t(Scenario1[,c(31:40,42:53)])

Alt.vitD_meat <- Alt.vitD_meat * VitD_meat

VitD_exp_scen1$VitD_meat <- colSums(Alt.vitD_meat) #Add sum of exposures for each individual to the Scenario1 dataset


#### Total exp from fish and meat ####

VitD_exp_scen1$VitD_total <- VitD_exp_scen1$VitD_fish + VitD_exp_scen1$VitD_meat #Only total vit D from fish and meat


#### Scen 1 change in vit D intake ####

VitD_exp_scen1$VitD_diff <- VitD_exp_scen1$VitD_total - VitD_exp$VitD_total

summary(VitD_exp_scen1$VitD_diff)


Nutrient_status <- read.csv2("NutrientDataTotal_IP.csv") #Current status of micronutrients
Nutrient_status[is.na(Nutrient_status)] <- 0


#Test no. of reporting days - create new variable "count"
setDT(Nutrient_status)[, count:=.N, by = .(PROJEKTID, PROJEKTID)]

Mean_status <- aggregate(Nutrient_status, by=list(Nutrient_status$PROJEKTID), 
                         FUN=mean)


FoodDaily <- read.csv("FoodDaily.csv")
FoodDaily_IP <- as.data.frame(FoodDaily[,1])
names(FoodDaily_IP)[1] <- 'PROJEKTID'

Status_IP_mean <- merge(Mean_status,FoodDaily_IP,by="PROJEKTID") #only the 2811 individuals we are considering (>15y, reported all 7 days)

Status_vitD <- cbind(Status_IP_mean[,c(1,9)])#vitamin D

VitD_exp_scen1 <- merge(VitD_exp_scen1, Status_vitD, by="PROJEKTID")

names(VitD_exp_scen1)[11] <- "Status"

summary(VitD_exp_scen1$Status) #Current vit D intake

tapply(VitD_exp_scen1$Status, VitD_exp_scen1$sex, mean)

VitD_exp_scen1$Status_new <- VitD_exp_scen1$Status + VitD_exp_scen1$VitD_diff

summary(VitD_exp_scen1$Status_new) 

tapply(VitD_exp_scen1$Status_new, VitD_exp_scen1$sex, mean)



##### Fraction < RI reference #####


#NNR 2012: Recommended intake (RI) of vitamin D i 10 ug/day; AR 7.5 ug/day

#Males
vitD_male_ref <- subset(VitD_exp_scen1, sex == 1, select = Status)
mean(vitD_male_ref$Status)


sd(vitD_male_ref$Status)


median(vitD_male_ref$Status)


mean(vitD_male_ref$Status < 7.5)



#Females

vitD_female_ref <- subset(VitD_exp_scen1, sex == 2, select = Status)
mean(vitD_female_ref$Status)


sd(vitD_female_ref$Status)


median(vitD_female_ref$Status)


mean(vitD_female_ref$Status < 7.5)




##### Fraction < RI scenario 1 #####

#Males
vitD_male_scen1 <- subset(VitD_exp_scen1, sex == 1, select = Status_new)
mean(vitD_male_scen1$Status_new)


sd(vitD_male_scen1$Status_new)


median(vitD_male_scen1$Status_new)


mean(vitD_male_scen1$Status_new < 7.5)



#Females

vitD_female_scen1 <- subset(VitD_exp_scen1, sex == 2, select = Status_new)
mean(vitD_female_scen1$Status_new)


sd(vitD_female_scen1$Status_new)


median(vitD_female_scen1$Status_new)


mean(vitD_female_scen1$Status_new < 7.5)




##### Scenario 2 #####

Scenario2 <- read.csv("Scenario2.csv")
# fish intakes are raw weights, ready for exposure assessment.


#### Exposure from fish ####

scen2.vitD_fish <- t(Scenario2[,c(55:60)]) #Create data set only with intakes of the different fish species

scen2.vitD_fish <- scen2.vitD_fish * VitD_fish[c(3,9,10,11,12,16)] ##Multiply intakes of fish species with concentrations

VitD_exp_scen2 <- Scenario2[,c(1:4, 54)]

setDT(VitD_exp_scen2)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

VitD_exp_scen2$VitD_fish <- colSums(scen2.vitD_fish) #Add sum of exposures for each individual to the Scenario2 dataset


#### Exposure from red and processed meat ####

#Meat intake is the same in all 4 alternative scenarios

VitD_exp_scen2$VitD_meat <- colSums(Alt.vitD_meat) #Add sum of exposures for each individual to the Scenario2 dataset


#### Total exp from fish and meat ####


VitD_exp_scen2$VitD_total <- VitD_exp_scen2$VitD_fish + VitD_exp_scen2$VitD_meat #Only total vit D from fish and meat


#### Scen 2 change in vit D intake ####

VitD_exp_scen2$VitD_diff <- VitD_exp_scen2$VitD_total - VitD_exp$VitD_total

summary(VitD_exp_scen2$VitD_diff)


VitD_exp_scen2 <- merge(VitD_exp_scen2, Status_vitD, by="PROJEKTID")

names(VitD_exp_scen2)[11] <- "Status"

summary(VitD_exp_scen2$Status) #Current vit D intake


VitD_exp_scen2$Status_new <- VitD_exp_scen2$Status + VitD_exp_scen2$VitD_diff

summary(VitD_exp_scen2$Status_new)


##### Fraction < RI scenario 2 #####

#Males
vitD_male_scen2 <- subset(VitD_exp_scen2, sex == 1, select = Status_new)
mean(vitD_male_scen2$Status_new)


sd(vitD_male_scen2$Status_new)


median(vitD_male_scen2$Status_new)


mean(vitD_male_scen2$Status_new < 7.5)



#Females

vitD_female_scen2 <- subset(VitD_exp_scen2, sex == 2, select = Status_new)
mean(vitD_female_scen2$Status_new)


sd(vitD_female_scen2$Status_new)


median(vitD_female_scen2$Status_new)


mean(vitD_female_scen2$Status < 7.5)



##### Scenario 3 #####

Scenario3 <- read.csv("Scenario3.csv")
# fish intakes are raw weights, ready for exposure assessment.


#### Exposure from fish ####

scen3.vitD_fish <- t(Scenario3[,c(55:64)]) #Create data set only with intakes of the different fish species

scen3.vitD_fish <- scen3.vitD_fish * VitD_fish[c(1,2,4,5,6,7,8,13,14,15)] ##Multiply intakes of fish species with concentrations

VitD_exp_scen3 <- Scenario3[,c(1:4, 54)]

setDT(VitD_exp_scen3)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

VitD_exp_scen3$VitD_fish <- colSums(scen3.vitD_fish) #Add sum of exposures for each individual to the Scenario3 dataset



#### Exposure from red and processed meat ####

#Meat intake is the same in all 4 alternative scenarios

VitD_exp_scen3$VitD_meat <- colSums(Alt.vitD_meat) #Add sum of exposures for each individual to the Scenario3 dataset


#### Total exp from fish and meat ####


VitD_exp_scen3$VitD_total <- VitD_exp_scen3$VitD_fish + VitD_exp_scen3$VitD_meat #Only total vit D from fish and meat


#### Change in vit D intake ####

VitD_exp_scen3$VitD_diff <- VitD_exp_scen3$VitD_total - VitD_exp$VitD_total

summary(VitD_exp_scen3$VitD_diff)


VitD_exp_scen3 <- merge(VitD_exp_scen3, Status_vitD, by="PROJEKTID")

names(VitD_exp_scen3)[11] <- "Status"

summary(VitD_exp_scen3$Status) #Current vit D intake


VitD_exp_scen3$Status_new <- VitD_exp_scen3$Status + VitD_exp_scen3$VitD_diff

summary(VitD_exp_scen3$Status_new) 


##### Fraction < RI scenario 3 #####

#Males
vitD_male_scen3 <- subset(VitD_exp_scen3, sex == 1, select = Status_new)
mean(vitD_male_scen3$Status_new)


sd(vitD_male_scen3$Status_new)


median(vitD_male_scen3$Status_new)


mean(vitD_male_scen3$Status_new < 7.5)



#Females

vitD_female_scen3 <- subset(VitD_exp_scen3, sex == 2, select = Status_new)
mean(vitD_female_scen3$Status_new)


sd(vitD_female_scen3$Status_new)


median(vitD_female_scen3$Status_new)


mean(vitD_female_scen3$Status < 7.5)



##### Scenario 4 #####

Scenario4 <- read.csv("Scenario4.csv")
# fish intakes are raw weights, ready for exposure assessment.


##### Exposure from fish #####

scen4.vitD_fish <- t(Scenario4[,55]) #Create data set only with intakes of the different fish species

scen4.vitD_fish <- scen4.vitD_fish * VitD_fish[8] ##Multiply intakes of fish species with concentrations

VitD_exp_scen4 <- Scenario4[,c(1:4, 54)]

setDT(VitD_exp_scen4)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

VitD_exp_scen4$VitD_fish <- colSums(scen4.vitD_fish) #Add sum of exposures for each individual to the Scenario4 dataset


##### Exposure from meat #####

#Meat intake is the same in all 4 alternative scenarios

VitD_exp_scen4$VitD_meat <- colSums(Alt.vitD_meat) #Add sum of exposures for each individual to the Scenario4 dataset


##### Total exp from fish and meat #####


VitD_exp_scen4$VitD_total <- VitD_exp_scen4$VitD_fish + VitD_exp_scen4$VitD_meat #Only total vit D from fish and meat


##### Scen 2 change in vit D intake #####

VitD_exp_scen4$VitD_diff <- VitD_exp_scen4$VitD_total - VitD_exp$VitD_total

summary(VitD_exp_scen4$VitD_diff)


VitD_exp_scen4 <- merge(VitD_exp_scen4, Status_vitD, by="PROJEKTID")

names(VitD_exp_scen4)[11] <- "Status"

summary(VitD_exp_scen4$Status) #Current vit D intake


tapply(VitD_exp_scen4$Status, VitD_exp_scen4$agegroups, mean)


VitD_exp_scen4$Status_new <- VitD_exp_scen4$Status + VitD_exp_scen4$VitD_diff

summary(VitD_exp_scen4$Status_new) #Note: negative values



##### Fraction < RI scenario 4 #####

#Males
vitD_male_scen4 <- subset(VitD_exp_scen4, sex == 1, select = Status_new)
mean(vitD_male_scen4$Status_new)


sd(vitD_male_scen4$Status_new)


median(vitD_male_scen4$Status_new)


mean(vitD_male_scen4$Status_new < 7.5)



#Females

vitD_female_scen4 <- subset(VitD_exp_scen4, sex == 2, select = Status_new)
mean(vitD_female_scen4$Status_new)


sd(vitD_female_scen4$Status_new)


median(vitD_female_scen4$Status_new)


mean(vitD_female_scen4$Status < 7.5)

