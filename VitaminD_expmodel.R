setwd("//dtu-storage/sthth/Documents/Case 1/Data")


#Packages
library(fitdistrplus)
library(data.table)
library(plyr)

######################################################### Reference scenario ################################################################

############################################################# Fish ##########################################################################

FishDaily <- read.csv("FishDaily.csv")
#NB! Water loss has NOT been taken into account in FoodDaily/FishDaily dataset so amounts are prepared amounts.
#Water loss HAS been taken into account in the substitution Scenarios

summary(FishDaily$total.fish)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   7.464  23.249  31.499  45.966 225.586

summary(FoodDaily$total.fish) #For the whole population
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   7.464  23.249  31.499  45.966 225.586  

#Convert fish consumption amounts into raw weights - assume that all fish is prepared and have a water loss of 20%
wloss.fish <- 1/0.8
FishDaily[,8:23] <- FishDaily[,8:23]*wloss.fish
FishDaily$Totalfish.raw <- rowSums(FishDaily[,8:23]) #Test if it gives the same as total.fish * 1.25


summary(FishDaily$Totalfish.raw) #For the whole population
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    9.33   29.06   39.37   57.46  281.98  

Nutrients_fish <- read.csv2("NutrientData_fish.csv") 

rownames(Nutrients_fish) <- Nutrients_fish[,1]
Nutrients_fish[,1] <- NULL

VitD_fish <- Nutrients_fish[,c(2)]

VitD_exp_fish <- t(FishDaily[,c(8:23)]) #Create data set only with intakes of the different fish species
VitD_exp_fish <- VitD_exp_fish * VitD_fish

VitD.contr.fish <- rowMeans(VitD_exp_fish) #Mean daily exposure in ?g/day

barplot(VitD.contr.fish, main = "Vitamin D contribution", las=2, ylim = c(0,1),  ylab = '?g/day') #Contribution to mean daily exposure ?g VitD/day

VitD_exp <- FishDaily[,1:4]

agebreaks <- c(15,20,25,30,35,40,45,50,55,60,65,70,75,80)
agelabels= c("15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79")

setDT(VitD_exp)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

VitD_exp$VitD_fish <- colSums(VitD_exp_fish) #Add sum of exposures for each individual to the FishDaily dataset

summary(VitD_exp$VitD_fish)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1645  1.2503  2.3262  3.4262 24.9037 

tapply(VitD_exp$VitD_fish, VitD_exp$agegroups, mean)
# 15-19     20-24     25-29     30-34     35-39     40-44     45-49     50-54     55-59     60-64     65-69     70-74     75-79 
# 0.9172004 1.4487999 1.7439746 1.6783105 1.7970700 1.8898510 2.0346494 2.6026212 2.6136370 3.3228925 3.4805401 3.6041749 4.2014302 


############################################################# Meat ##########################################################################

FoodDaily <- read.csv("FoodDaily.csv")

#NB! Water loss has NOT been taken into account in FoodDaily/FishDaily dataset so amounts are prepared amounts.
#Water loss HAS been taken into account in the substitution Scenarios

#Convert red meat amounts into raw weights - 25% water loss. Processed meat do not need to be converted because we have
#fat contents for individual processed meats (and also concentration data)
wloss.redmeat <- 1/0.75
FoodDaily[,46:57] <- FoodDaily[,46:57]*wloss.redmeat

mean(rowSums(FoodDaily[,46:57])) #Test if it gives the same as red.meat.raw (approx. - not exactly equal due to assumptions on proportions of meat)
mean(FoodDaily$red.meat.raw)

summary(FoodDaily$red.meat.raw) #For the whole population
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   58.35   90.77  101.50  132.39  443.66 

summary(FoodDaily$red.meat)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   43.76   68.08   76.13   99.29  332.74

summary(FoodDaily$proc.meat)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   12.50   28.29   39.04   54.21  467.14  


Nutrients_meat <- read.csv2("NutrientData_meat.csv") 

rownames(Nutrients_meat) <- Nutrients_meat[,1]
Nutrients_meat[,1] <- NULL

VitD_meat <- Nutrients_meat[,c(2)]

VitD_exp_meat <- t(FoodDaily[,c(35:44,46:57)]) #Create data set only with intakes of the different fish species
VitD_exp_meat <- VitD_exp_meat * VitD_meat

VitD.contr.meat <- rowMeans(VitD_exp_meat) #Mean daily exposure in ?g/day

barplot(VitD.contr.meat, main = "Vitamin D contribution", las=2, ylim = c(0,0.25),  ylab = '?g/day') 

VitD_exp$VitD_meat <- colSums(VitD_exp_meat) #Add sum of exposures for each individual to the FishDaily dataset

summary(VitD_exp$VitD_meat)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.4584  0.6906  0.7829  1.0158  3.4800 

tapply(VitD_exp$VitD_meat, VitD_exp$agegroups, mean)
# 15-19     20-24     25-29     30-34     35-39     40-44     45-49     50-54     55-59     60-64     65-69     70-74     75-79 
# 0.8060094 0.7450743 0.8485123 0.8571913 0.8503350 0.8281512 0.7897611 0.8176342 0.7815916 0.7611351 0.6910015 0.6448728 0.7527620 

############################################# Total exp from fish and meat ############################################

VitD_exp$VitD_total <- VitD_exp$VitD_fish + VitD_exp$VitD_meat #Only total vit D from fish and meat


################################################ Scenario 1 ###########################################################

Scenario1 <- read.csv("Scenario1.csv")
# fish intakes are raw weights, ready for exposure assessment.

summary(Scenario1$total.fish.new)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 50.00   50.00   50.00   56.52   50.00  225.59 

summary(Scenario1$fish.total.raw)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 62.50   62.50   62.50   70.65   62.50  281.98

################################################## Fish ###############################################################

scen1.vitD_fish <- t(Scenario1[,c(55:70)]) #Create data set only with intakes of the different fish species

scen1.vitD_fish <- scen1.vitD_fish * VitD_fish ##Multiply intakes of fish species with concentrations

VitD_exp_scen1 <- Scenario1[,c(1:4, 54)]

setDT(VitD_exp_scen1)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

VitD_exp_scen1$VitD_fish <- colSums(scen1.vitD_fish) #Add sum of exposures for each individual to the Scenario1 dataset

summary(VitD_exp_scen1$VitD_fish)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.575   3.429   3.429   3.992   3.742  18.038 

tapply(VitD_exp_scen1$VitD_fish, VitD_exp_scen1$agegroups, mean)
# 15-19    20-24    25-29    30-34    35-39    40-44    45-49    50-54    55-59    60-64    65-69    70-74    75-79 
# 3.561494 3.844000 3.719097 3.763062 3.872813 3.890911 3.842565 4.063538 4.116313 4.286590 4.286872 4.444481 4.314000

################################################## Meat ###############################################################

#Meat intake is the same in all 4 alternative scenarios

Alt.vitD_meat <- t(Scenario1[,c(31:40,42:53)])

Alt.vitD_meat <- Alt.vitD_meat * VitD_meat

VitD_exp_scen1$VitD_meat <- colSums(Alt.vitD_meat) #Add sum of exposures for each individual to the Scenario1 dataset

summary(VitD_exp_scen1$VitD_meat)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.3336  0.5747  0.6544  0.8851  3.2055 

tapply(VitD_exp_scen1$VitD_meat, VitD_exp_scen1$agegroups, mean)
# 15-19     20-24     25-29     30-34     35-39     40-44     45-49     50-54     55-59     60-64     65-69     70-74     75-79 
# 0.6051155 0.5971415 0.7031919 0.6949869 0.6770383 0.6794014 0.6459724 0.7074587 0.6693204 0.6808557 0.6151622 0.5754677 0.6696064 

########################################### Total exp from fish and meat ##############################################


VitD_exp_scen1$VitD_total <- VitD_exp_scen1$VitD_fish + VitD_exp_scen1$VitD_meat #Only total vit D from fish and meat

############################################ Change in vit D intake ###################################################

VitD_exp_scen1$VitD_diff <- VitD_exp_scen1$VitD_total - VitD_exp$VitD_total

summary(VitD_exp_scen1$VitD_diff)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -12.0877   0.4523   2.2292   1.5374   3.0595  10.9263 

tapply(VitD_exp_scen1$VitD_diff, VitD_exp_scen1$agegroups, mean)
# 15-19      20-24      25-29      30-34      35-39      40-44      45-49      50-54      55-59      60-64      65-69      70-74      75-79 
# 2.44340017 2.24726754 1.82980205 1.92254728 1.90244670 1.85230977 1.66412641 1.35074137 1.39040516 0.88341782 0.73049263 0.77090091 0.02941426 

Nutrient_status <- read.csv2("NutrientDataTotal_IP.csv") #Current status of micronutrients
Nutrient_status[is.na(Nutrient_status)] <- 0

#Test no. of reporting days - create new variable "count"
setDT(Nutrient_status)[, count:=.N, by = .(PROJEKTID, PROJEKTID)]

count(Nutrient_status,c('count'))
# count  freq
# 1     7 26642

Mean_status <- aggregate(Nutrient_status, by=list(Nutrient_status$PROJEKTID), 
                         FUN=mean)
#3952 individuals

FoodDaily <- read.csv("FoodDaily.csv")
FoodDaily_IP <- as.data.frame(FoodDaily[,1])
names(FoodDaily_IP)[1] <- 'PROJEKTID'

Status_IP_mean <- merge(Mean_status,FoodDaily_IP,by="PROJEKTID") #only the 2811 individuals we are considering (>15y, reported all 7 days)

Status_vitD <- cbind(Status_IP_mean[,c(1,9)])#X23 is vitamin D

VitD_exp_scen1 <- merge(VitD_exp_scen1, Status_vitD, by="PROJEKTID")

names(VitD_exp_scen1)[11] <- "Status"

summary(VitD_exp_scen1$Status) #Current vit D intake
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.360   2.122   3.281   4.726   5.838  49.343 

tapply(VitD_exp_scen1$Status, VitD_exp_scen1$agegroups, mean)
# 15-19    20-24    25-29    30-34    35-39    40-44    45-49    50-54    55-59    60-64    65-69    70-74    75-79 
# 2.973812 3.771233 4.229046 4.312944 4.035335 4.286943 4.478281 5.197517 5.457295 5.713219 5.848397 5.612916 6.035042 

tapply(VitD_exp_scen1$Status, VitD_exp_scen1$sex, mean)
# 1        2 
# 5.296205 4.175419 

VitD_exp_scen1$Status_new <- VitD_exp_scen1$Status + VitD_exp_scen1$VitD_diff

tapply(VitD_exp_scen1$Status_new, VitD_exp_scen1$agegroups, mean)
# 15-19    20-24    25-29    30-34    35-39    40-44    45-49    50-54    55-59    60-64    65-69    70-74    75-79 
# 5.417212 6.018500 6.058849 6.235491 5.937782 6.139253 6.142407 6.548258 6.847701 6.596804 6.578890 6.383817 6.064456 

summary(VitD_exp_scen1$Status_new) #Note: negative values
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.04357  4.53044  5.32732  6.26303  6.78635 44.77446

tapply(VitD_exp_scen1$Status_new, VitD_exp_scen1$sex, mean)
# 1        2 
# 6.784484 5.760194 

########################################### Fraction < RI reference ###################################################

#NNR 2012: Recommended intake (RI) of vitamin D i 10 ug/day; AR 7.5 ug/day

#Males
vitD_male_ref <- subset(VitD_exp_scen1, sex == 1, select = Status)
mean(vitD_male_ref$Status)
# [1] 5.313006

sd(vitD_male_ref$Status)
# [1] 4.692671

median(vitD_male_ref$Status)
# [1] 3.669864

mean(vitD_male_ref$Status < 7.5)
# [1] 0.8021739


#Females

vitD_female_ref <- subset(VitD_exp_scen1, sex == 2, select = Status)
mean(vitD_female_ref$Status)
# [1] 4.198355

sd(vitD_female_ref$Status)
# [1] 3.704536

median(vitD_female_ref$Status)
# [1] 2.944767

mean(vitD_female_ref$Status < 7.5)
# [1] 0.85884



########################################### Fraction < RI scenario 1 ##################################################

#Males
vitD_male_scen1 <- subset(VitD_exp_scen1, sex == 1, select = Status_new)
mean(vitD_male_scen1$Status_new)
# [1] 6.801253

sd(vitD_male_scen1$Status_new)
# [1] 3.491093

median(vitD_male_scen1$Status_new)
# [1] 5.765313

mean(vitD_male_scen1$Status_new < 7.5)
# [1] 0.7666667


#Females

vitD_female_scen1 <- subset(VitD_exp_scen1, sex == 2, select = Status_new)
mean(vitD_female_scen1$Status_new)
# [1] 5.78313

sd(vitD_female_scen1$Status_new)
# [1] 2.695201

median(vitD_female_scen1$Status_new)
# [1] 4.957429

mean(vitD_female_scen1$Status_new < 7.5)
# [1] 0.8469602




################################################ Scenario 2 ###########################################################

Scenario2 <- read.csv("Scenario2.csv")
# fish intakes are raw weights, ready for exposure assessment.

summary(Scenario2$total.fish.new)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 50.00   50.00   50.00   56.52   50.00  225.59 

summary(Scenario2$fish.total.raw)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 62.50   62.50   62.50   70.65   62.50  281.98

################################################## Fish ###############################################################

scen2.vitD_fish <- t(Scenario2[,c(55:60)]) #Create data set only with intakes of the different fish species

scen2.vitD_fish <- scen2.vitD_fish * VitD_fish[c(3,9,10,11,12,16)] ##Multiply intakes of fish species with concentrations

VitD_exp_scen2 <- Scenario2[,c(1:4, 54)]

setDT(VitD_exp_scen2)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

VitD_exp_scen2$VitD_fish <- colSums(scen2.vitD_fish) #Add sum of exposures for each individual to the Scenario2 dataset

summary(VitD_exp_scen2$VitD_fish)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.322   5.989   5.989   6.861   6.317  27.169 

tapply(VitD_exp_scen2$VitD_fish, VitD_exp_scen2$agegroups, mean)
# 15-19    20-24    25-29    30-34    35-39    40-44    45-49    50-54    55-59    60-64    65-69    70-74    75-79 
# 6.199263 6.640462 6.421656 6.524476 6.657330 6.696252 6.606387 6.959915 7.102290 7.337851 7.352182 7.523351 7.190634 

################################################## Meat ###############################################################

#Meat intake is the same in all 4 alternative scenarios

VitD_exp_scen2$VitD_meat <- colSums(Alt.vitD_meat) #Add sum of exposures for each individual to the Scenario2 dataset


########################################### Total exp from fish and meat ##############################################


VitD_exp_scen2$VitD_total <- VitD_exp_scen2$VitD_fish + VitD_exp_scen2$VitD_meat #Only total vit D from fish and meat


############################################ Change in vit D intake ###################################################

VitD_exp_scen2$VitD_diff <- VitD_exp_scen2$VitD_total - VitD_exp$VitD_total

summary(VitD_exp_scen2$VitD_diff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.3732  3.2934  4.9386  4.4068  5.6475 20.3780  

tapply(VitD_exp_scen2$VitD_diff, VitD_exp_scen2$agegroups, mean)
# 15-19    20-24    25-29    30-34    35-39    40-44    45-49    50-54    55-59    60-64    65-69    70-74    75-79 
# 5.081169 5.043729 4.532361 4.683961 4.686963 4.657652 4.427949 4.247118 4.376382 3.934679 3.795803 3.849771 2.906048

VitD_exp_scen2 <- merge(VitD_exp_scen2, Status_vitD, by="PROJEKTID")

names(VitD_exp_scen2)[11] <- "Status"

summary(VitD_exp_scen2$Status) #Current vit D intake
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.360   2.122   3.281   4.726   5.838  49.343 

tapply(VitD_exp_scen2$Status, VitD_exp_scen2$agegroups, mean)
# 15-19    20-24    25-29    30-34    35-39    40-44    45-49    50-54    55-59    60-64    65-69    70-74    75-79 
# 2.973812 3.771233 4.229046 4.312944 4.035335 4.286943 4.478281 5.197517 5.457295 5.713219 5.848397 5.612916 6.035042 

VitD_exp_scen2$Status_new <- VitD_exp_scen2$Status + VitD_exp_scen2$VitD_diff

summary(VitD_exp_scen2$Status_new) #Note: negative values
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.473   7.100   7.912   9.132   9.572  50.677 


########################################### Fraction < RI scenario 2 ##################################################

#Males
vitD_male_scen2 <- subset(VitD_exp_scen2, sex == 1, select = Status_new)
mean(vitD_male_scen2$Status_new)
# [1] 9.779809

sd(vitD_male_scen2$Status_new)
# [1] 4.267333

median(vitD_male_scen2$Status_new)
# [1] 8.376606

mean(vitD_male_scen2$Status_new < 7.5)
# [1] 0.2528986


#Females

vitD_female_scen2 <- subset(VitD_exp_scen2, sex == 2, select = Status_new)
mean(vitD_female_scen2$Status_new)
# [1] 8.547292

sd(vitD_female_scen2$Status_new)
# [1] 3.218545

median(vitD_female_scen2$Status_new)
# [1] 7.542961

mean(vitD_female_scen2$Status < 7.5)
# [1] 0.4884696

################################################ Scenario 3 ###########################################################

Scenario3 <- read.csv("Scenario3.csv")
# fish intakes are raw weights, ready for exposure assessment.

summary(Scenario3$total.fish.new)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 50.00   50.00   50.00   56.52   50.00  225.59 

summary(Scenario3$fish.total.raw)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 62.50   62.50   62.50   70.65   62.50  281.98

################################################## Fish ###############################################################

scen3.vitD_fish <- t(Scenario3[,c(55:64)]) #Create data set only with intakes of the different fish species

scen3.vitD_fish <- scen3.vitD_fish * VitD_fish[c(1,2,4,5,6,7,8,13,14,15)] ##Multiply intakes of fish species with concentrations

VitD_exp_scen3 <- Scenario3[,c(1:4, 54)]

setDT(VitD_exp_scen3)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

VitD_exp_scen3$VitD_fish <- colSums(scen3.vitD_fish) #Add sum of exposures for each individual to the Scenario3 dataset

summary(VitD_exp_scen3$VitD_fish)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4461  1.0818  1.0818  1.3093  1.2368  7.6197 

tapply(VitD_exp_scen3$VitD_fish, VitD_exp_scen3$agegroups, mean)
# 15-19    20-24    25-29    30-34    35-39    40-44    45-49    50-54    55-59    60-64    65-69    70-74    75-79 
# 1.133499 1.245871 1.206738 1.209047 1.269833 1.272294 1.259446 1.343598 1.337947 1.419171 1.413111 1.509121 1.514520 

################################################## Meat ###############################################################

#Meat intake is the same in all 4 alternative scenarios

VitD_exp_scen3$VitD_meat <- colSums(Alt.vitD_meat) #Add sum of exposures for each individual to the Scenario3 dataset


########################################### Total exp from fish and meat ##############################################


VitD_exp_scen3$VitD_total <- VitD_exp_scen3$VitD_fish + VitD_exp_scen3$VitD_meat #Only total vit D from fish and meat


############################################ Change in vit D intake ###################################################

VitD_exp_scen3$VitD_diff <- VitD_exp_scen3$VitD_total - VitD_exp$VitD_total

summary(VitD_exp_scen3$VitD_diff)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -21.9814  -2.1849  -0.2789  -1.1454   0.6927   3.9819  

tapply(VitD_exp_scen3$VitD_diff, VitD_exp_scen3$agegroups, mean)
# 15-19       20-24       25-29       30-34       35-39       40-44       45-49       50-54       55-59       60-64       65-69 
# 0.01540467 -0.35086190 -0.68255670 -0.63146809 -0.70053389 -0.76630655 -0.91899197 -1.36919880 -1.38796159 -1.98400064 -2.14326838 
# 70-74       75-79 
# -2.16445926 -2.77006550 

VitD_exp_scen3 <- merge(VitD_exp_scen3, Status_vitD, by="PROJEKTID")

names(VitD_exp_scen3)[11] <- "Status"

summary(VitD_exp_scen3$Status) #Current vit D intake
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.360   2.122   3.281   4.726   5.838  49.343 

tapply(VitD_exp_scen3$Status, VitD_exp_scen3$agegroups, mean)
# 15-19    20-24    25-29    30-34    35-39    40-44    45-49    50-54    55-59    60-64    65-69    70-74    75-79 
# 2.973812 3.771233 4.229046 4.312944 4.035335 4.286943 4.478281 5.197517 5.457295 5.713219 5.848397 5.612916 6.035042 

VitD_exp_scen3$Status_new <- VitD_exp_scen3$Status + VitD_exp_scen3$VitD_diff

summary(VitD_exp_scen3$Status_new) #Note: negative values
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -4.757   2.090   2.850   3.580   4.105  39.347 


########################################### Fraction < RI scenario 3 ##################################################

#Males
vitD_male_scen3 <- subset(VitD_exp_scen3, sex == 1, select = Status_new)
mean(vitD_male_scen3$Status_new)
# [1] 4.008995

sd(vitD_male_scen3$Status_new)
# [1] 3.017512

median(vitD_male_scen3$Status_new)
# [1] 3.230044

mean(vitD_male_scen3$Status_new < 7.5)
# [1] 0.9137681


#Females

vitD_female_scen3 <- subset(VitD_exp_scen3, sex == 2, select = Status_new)
mean(vitD_female_scen3$Status_new)
# [1] 3.20596

sd(vitD_female_scen3$Status_new)
# [1] 2.356175

median(vitD_female_scen3$Status_new)
# [1] 2.537636

mean(vitD_female_scen3$Status < 7.5)
# [1] 0.944095

################################################ Scenario 4 ###########################################################

Scenario4 <- read.csv("Scenario4.csv")
# fish intakes are raw weights, ready for exposure assessment.

summary(Scenario4$total.fish.new)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 50.00   50.00   50.00   56.52   50.00  225.59 

summary(Scenario4$fish.total.raw)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 62.50   62.50   62.50   70.65   62.50  281.98

################################################## Fish ###############################################################

scen4.vitD_fish <- t(Scenario4[,55]) #Create data set only with intakes of the different fish species

scen4.vitD_fish <- scen4.vitD_fish * VitD_fish[8] ##Multiply intakes of fish species with concentrations

VitD_exp_scen4 <- Scenario4[,c(1:4, 54)]

setDT(VitD_exp_scen4)[ , agegroups:= cut(age, breaks= agebreaks, right= FALSE, labels= agelabels)]

VitD_exp_scen4$VitD_fish <- colSums(scen4.vitD_fish) #Add sum of exposures for each individual to the Scenario4 dataset

summary(VitD_exp_scen4$VitD_fish)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.645   1.645   1.645   1.860   1.645   7.423 

tapply(VitD_exp_scen4$VitD_fish, VitD_exp_scen4$agegroups, mean)
# 15-19    20-24    25-29    30-34    35-39    40-44    45-49    50-54    55-59    60-64    65-69    70-74    75-79 
# 1.698033 1.807596 1.747352 1.781386 1.804784 1.817102 1.791212 1.881101 1.931378 1.982318 1.989333 2.013181 1.898468  

################################################## Meat ###############################################################

#Meat intake is the same in all 4 alternative scenarios

VitD_exp_scen4$VitD_meat <- colSums(Alt.vitD_meat) #Add sum of exposures for each individual to the Scenario4 dataset


########################################### Total exp from fish and meat ##############################################


VitD_exp_scen4$VitD_total <- VitD_exp_scen4$VitD_fish + VitD_exp_scen4$VitD_meat #Only total vit D from fish and meat


############################################ Change in vit D intake ###################################################

VitD_exp_scen4$VitD_diff <- VitD_exp_scen4$VitD_total - VitD_exp$VitD_total

summary(VitD_exp_scen4$VitD_diff)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -17.4145  -1.6512   0.2890  -0.5948   1.2591   4.1116 

tapply(VitD_exp_scen4$VitD_diff, VitD_exp_scen4$agegroups, mean)
# 15-19       20-24       25-29       30-34       35-39       40-44       45-49       50-54       55-59       60-64       65-69 
# 0.57993879  0.21086327 -0.14194328 -0.05912938 -0.16558274 -0.22149829 -0.38722625 -0.83169627 -0.79453068 -1.42085432 -1.56704624 
# 70-74       75-79 
# -1.66039893 -2.38611751 

VitD_exp_scen4 <- merge(VitD_exp_scen4, Status_vitD, by="PROJEKTID")

names(VitD_exp_scen4)[11] <- "Status"

summary(VitD_exp_scen4$Status) #Current vit D intake
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.360   2.122   3.281   4.726   5.838  49.343 

tapply(VitD_exp_scen4$Status, VitD_exp_scen4$agegroups, mean)
# 15-19    20-24    25-29    30-34    35-39    40-44    45-49    50-54    55-59    60-64    65-69    70-74    75-79 
# 2.973812 3.771233 4.229046 4.312944 4.035335 4.286943 4.478281 5.197517 5.457295 5.713219 5.848397 5.612916 6.035042 

VitD_exp_scen4$Status_new <- VitD_exp_scen4$Status + VitD_exp_scen4$VitD_diff

summary(VitD_exp_scen4$Status_new) #Note: negative values
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -3.267   2.607   3.381   4.131   4.647  40.620 


########################################### Fraction < RI scenario 4 ##################################################

#Males
vitD_male_scen4 <- subset(VitD_exp_scen4, sex == 1, select = Status_new)
mean(vitD_male_scen4$Status_new)
# [1] 4.56917

sd(vitD_male_scen4$Status_new)
# [1] 3.222737

median(vitD_male_scen4$Status_new)
# [1] 3.771658

mean(vitD_male_scen4$Status_new < 7.5)
# [1] 0.8833333


#Females

vitD_female_scen4 <- subset(VitD_exp_scen4, sex == 2, select = Status_new)
mean(vitD_female_scen4$Status_new)
# [1] 3.747368

sd(vitD_female_scen4$Status_new)
# [1] 2.503932

median(vitD_female_scen4$Status_new)
# [1] 3.079548

mean(vitD_female_scen4$Status < 7.5)
# [1] 0.9287212