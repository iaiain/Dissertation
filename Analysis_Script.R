# Using Multilevel Logistic Regression with lme4 to analyze the ICEMR PRISM Cohort data
# data on clinepidb.org
# Libraries                       ----                      
library(tidyverse)
library(forcats)
library(lubridate)
library(lme4)
library(performance)  # for calculating pseudo-R^2 using Tjur, with r2()
library(DHARMa) # for model diagnostics

# Explanation of Data after running all intro sections: ----
#   Households: one row is one household, original data from file with factored categorical vars
#   H.Measures: one row is one household measure, original data from file with factored categorical vars (UNUSED)
#   Participants: one row is one participant, original data from file with factored categorical vars
#   P.Measures: one row is one participant measure, original data from file with factored categorical vars (- one row for missing diagnosis)
#   Samples : one row is one sample, original data from file with factored categorical vars (UNUSED)
#   Kihihi, Walukuba, and Nagongera: subsets of the original Households table, filtered on subcounty

#   master: the main table for the mlms - one row for each observation and has data on the observations, the 
#           participant, and on their household. The primary outcome variable is Malaria

# Intro                           ----
Households    <- read.delim("PRISM_cohort_Households.txt")
H.Measures    <- read.delim("PRISM_cohort_Household_repeated_measures.txt")
Participants  <- read.delim("PRISM_cohort_Participants.txt")
P.Measures    <- read.delim("PRISM_cohort_Participant_repeated_measures.txt")
Samples       <- read.delim("PRISM_cohort_Samples.txt")

# Cleaning and Factoring          ----
# summary(Households)
Households <- Households %>% 
  mutate(
    Air.bricks..EUPATH_0000018. = factor(Air.bricks..EUPATH_0000018.),
    Animal.drawn.cart..EUPATH_0000166. = factor(Animal.drawn.cart..EUPATH_0000166.),
    Bank.account..EUPATH_0000167. = factor(Bank.account..EUPATH_0000167.),
    Bed..ENVO_00000501. = factor(Bed..ENVO_00000501.),
    Bicycle..ENVO_01000614. = factor(Bicycle..ENVO_01000614.),
    Boat.with.a.motor..EUPATH_0000179. = factor(Boat.with.a.motor..EUPATH_0000179.),
    Boat.without.a.motor..EUPATH_0000170. = factor(Boat.without.a.motor..EUPATH_0000170.),
    Car.or.truck..EUPATH_0000171. = factor(Car.or.truck..EUPATH_0000171.),
    Cassette.player..ENVO_01000578. = factor(Cassette.player..ENVO_01000578.),
    Chair..ENVO_01000586. = factor(Chair..ENVO_01000586.),
    Clock..ENVO_01000596. = factor(Clock..ENVO_01000596.),
    Cooking.fuel..EUPATH_0000023. = factor(Cooking.fuel..EUPATH_0000023.),
    Country..ENVO_00000009. = factor(Country..ENVO_00000009.),
    Cupboard..ENVO_01000595. = factor(Cupboard..ENVO_01000595.),
    Drinking.water.source..ENVO_00003064. = factor(Drinking.water.source..ENVO_00003064.),
    Dwelling.type..ENVO_01000744. = factor(Dwelling.type..ENVO_01000744.),
    Dwelling.type..ENVO_01000744. = fct_relevel(Dwelling.type..ENVO_01000744.,"Traditional"),   # want traditional as reference
    Eaves..ENVO_01000825.= factor(Eaves..ENVO_01000825.),
    Eaves..ENVO_01000825. = fct_relevel(Eaves..ENVO_01000825.,"Open"),
    Electricity..EUPATH_0021084. = factor(Electricity..EUPATH_0021084.),
    Floor.material..EUPATH_0000006. = factor(Floor.material..EUPATH_0000006.),
    Earth.floor = case_when(
      Floor.material..EUPATH_0000006. == "Earth and dung" ~ "Yes",
      Floor.material..EUPATH_0000006. == "Earth or sand" ~ "Yes",
      Floor.material..EUPATH_0000006. == "Bricks" ~ "No",
      Floor.material..EUPATH_0000006. == "Cement/concrete" ~ "No",
      Floor.material..EUPATH_0000006. == "Stones" ~ "No",
      Floor.material..EUPATH_0000006. == "Parquet or polished wood" ~ "No",
      ),
    Earth.floor = factor(Earth.floor),
    Food.problems.per.week..EUPATH_0000029. = factor(Food.problems.per.week..EUPATH_0000029.),
    Household.wealth.index..categorical..EUPATH_0000143.= factor(Household.wealth.index..categorical..EUPATH_0000143.),
    Human.waste.facilities..EUPATH_0000335. = factor(Human.waste.facilities..EUPATH_0000335.),
    Landline.phone..ENVO_01000582. = factor(Landline.phone..ENVO_01000582.),
    Lighting.source..OBI_0400065. = factor(Lighting.source..OBI_0400065.),
    Mobile.phone..ENVO_01000581. = factor(Mobile.phone..ENVO_01000581.),
    Motorcycle.or.scooter..ENVO_01000615. = factor(Motorcycle.or.scooter..ENVO_01000615.),
    Radio..ENVO_01000577. = factor(Radio..ENVO_01000577.),
    Refrigerator..ENVO_01000583. = factor(Refrigerator..ENVO_01000583.),
    Roof.material..EUPATH_0000003. = factor(Roof.material..EUPATH_0000003.),
    Roof.material..EUPATH_0000003. = fct_relevel(Roof.material..EUPATH_0000003.,"Thatch"),
    Sofa..ENVO_01000588. = factor(Sofa..ENVO_01000588.),
    Sub.county.in.Uganda..EUPATH_0000054. = factor(Sub.county.in.Uganda..EUPATH_0000054.),
    Table..ENVO_01000584. = factor(Table..ENVO_01000584.),
    Television..ENVO_01000579. = factor(Television..ENVO_01000579.),
    Wall.material..EUPATH_0000009. = factor(Wall.material..EUPATH_0000009.),
    Wall.material..EUPATH_0000009. = fct_relevel(Wall.material..EUPATH_0000009.,"Mud"),
    Watch..EUPATH_0000186. = factor(Watch..EUPATH_0000186.)
  )

# Cleaning Participants
# summary(Participants)
Participants <- Participants %>% 
  mutate(
    Alpha.thalassemia.genotype..EUPATH_0000034. = factor(Alpha.thalassemia.genotype..EUPATH_0000034.),
    CD36.genotype..EUPATH_0000737. = factor(CD36.genotype..EUPATH_0000737.),
    Cause.of.death..EUPATH_0020001. = factor(Cause.of.death..EUPATH_0020001.),
    Death.location..EUPATH_0010067. = factor(Death.location..EUPATH_0010067.),
    G6PD.genotype..EUPATH_0000033. = factor(G6PD.genotype..EUPATH_0000033.),
    HbS.genotype..EUPATH_0000035. = factor(HbS.genotype..EUPATH_0000035.),
    Reason.for.withdrawal..EUPATH_0000208. = factor(Reason.for.withdrawal..EUPATH_0000208.),
    Sex..PATO_0000047. = factor(Sex..PATO_0000047.),
    Timing.of.enrollment..EUPATH_0000219. = factor(Timing.of.enrollment..EUPATH_0000219.),
    Death.date..EUPATH_0000150. = as.Date(Death.date..EUPATH_0000150.),
    Enrollment.date..EUPATH_0000151. = as.Date(Enrollment.date..EUPATH_0000151.),
    Last.date.observed..EUPATH_0000152. = as.Date(Last.date.observed..EUPATH_0000152.)
  )

# Cleaning H.Measures
# summary(H.Measures)
H.Measures <- H.Measures %>% 
  mutate(
    Collection.date..EUPATH_0020003. = as.Date(Collection.date..EUPATH_0020003.)
  )

# Cleaning P.Measures
# summary(P.Measures)
P.Measures <- P.Measures %>% filter(Malaria.diagnosis..EUPATH_0000090. != "")
P.Measures <- P.Measures %>% 
  mutate(
    Abdominal.pain..HP_0002027. = factor(Abdominal.pain..HP_0002027.),
    Admitting.hospital..EUPATH_0000318. = factor(Admitting.hospital..EUPATH_0000318.),
    Anorexia..SYMP_0000523. = factor(Anorexia..SYMP_0000523.),
    Antimalarial.medication..EUPATH_0000058. = factor(Antimalarial.medication..EUPATH_0000058.),
    Complicated.malaria..EUPATH_0000040. = factor(Complicated.malaria..EUPATH_0000040.),
    Complicated.malaria.signs..EUPATH_0020254. = factor(Complicated.malaria.signs..EUPATH_0020254.),
    Cough..SYMP_0000614. = factor(Cough..SYMP_0000614.),
    Diagnosis.at.hospitalization..EUPATH_0000638. = factor(Diagnosis.at.hospitalization..EUPATH_0000638.),
    Diarrhea..HP_0002014. = factor(Diarrhea..HP_0002014.),
    Fatigue..SYMP_0019177. = factor(Fatigue..SYMP_0019177.),
    Febrile..EUPATH_0000097. = factor(Febrile..EUPATH_0000097.),
    Headache..HP_0002315. = factor(Headache..HP_0002315.),
    ITN.last.night..EUPATH_0000216. = case_when(
      ITN.last.night..EUPATH_0000216. == "Yes" ~ "Yes",
      ITN.last.night..EUPATH_0000216. == "No" ~ "No",
      ITN.last.night..EUPATH_0000216. == "" ~ "NA"
    ),
    ITN.last.night..EUPATH_0000216. = factor(ITN.last.night..EUPATH_0000216.),
    Jaundice..HP_0000952. = factor(Jaundice..HP_0000952.),
    Joint.pains..HP_0002829. = factor(Joint.pains..HP_0002829.),
    Malaria.diagnosis..EUPATH_0000090. = factor(Malaria.diagnosis..EUPATH_0000090.),
    Malaria.diagnosis.and.parasite.status..EUPATH_0000338. = 
      factor(Malaria.diagnosis.and.parasite.status..EUPATH_0000338.),
    Muscle.aches..HP_0003326. = factor(Muscle.aches..HP_0003326.),
    Non.malaria.medication..EUPATH_0000059. = factor(Non.malaria.medication..EUPATH_0000059.),
    Observation.type..BFO_0000015. = factor(Observation.type..BFO_0000015.),
    Other.diagnosis..EUPATH_0000317. = factor(Other.diagnosis..EUPATH_0000317.),
    Other.medical.complaint..EUPATH_0020002. = factor(Other.medical.complaint..EUPATH_0020002.),
    Seizures..SYMP_0000124. = factor(Seizures..SYMP_0000124.),
    Severe.malaria.signs..EUPATH_0020184. = factor(Severe.malaria.signs..EUPATH_0020184.),
    Subjective.fever..EUPATH_0000100. = factor(Subjective.fever..EUPATH_0000100.),
    Vomiting..HP_0002013. = factor(Vomiting..HP_0002013.),
    Observation.date..EUPATH_0004991. = as.Date(Observation.date..EUPATH_0004991.),
    Hospital.discharge.date..EUPATH_0000320. = as.Date(Hospital.discharge.date..EUPATH_0000320.)
  )

# Cleaning Samples
# summary(Samples)
Samples <- Samples %>% 
  mutate(
    Buffy.coat.sample..OBIB_0000036. = factor(Buffy.coat.sample..OBIB_0000036.),
    Erythrocyte.sample..OBI_2000016. = factor(Erythrocyte.sample..OBI_2000016.),
    Filter.paper.sample..EUPATH_0000127. = factor(Filter.paper.sample..EUPATH_0000127.),
    Peripheral.blood.mononuclear.cell.sample..EUPATH_0000128. = 
      factor(Peripheral.blood.mononuclear.cell.sample..EUPATH_0000128.),
    Plasma.sample..OBI_0100016. = factor(Plasma.sample..OBI_0100016.),
    Plasmodium.asexual.stages..by.microscopy..EUPATH_0000048. = 
      factor(Plasmodium.asexual.stages..by.microscopy..EUPATH_0000048.),
    Plasmodium.gametocytes..by.microscopy..EUPATH_0000207. = 
      factor(Plasmodium.gametocytes..by.microscopy..EUPATH_0000207.),
    Plasmodium..by.LAMP..EUPATH_0000487. = factor(Plasmodium..by.LAMP..EUPATH_0000487.)
  )

# Derived Variables               ----

# PARTICIPANTS
Participants$Years.at.risk <- 
  difftime(Participants$Last.date.observed..EUPATH_0000152.,
           Participants$Enrollment.date..EUPATH_0000151., units = "days") # will convert to yrs

# clean up new column
Participants <- Participants %>% mutate(Years.at.risk = as.numeric(Years.at.risk)/365)

# P.MEASURES
P.Measures$Observation.month <- format(P.Measures$Observation.date..EUPATH_0004991., "%m")
P.Measures$Observation.month <- as.integer(P.Measures$Observation.month)

P.Measures <- P.Measures %>% 
  mutate(Peak.season = case_when(
      Observation.month > 3 & Observation.month < 7 ~ "Yes",
      Observation.month > 9 ~ "Yes",
      Observation.month < 4 ~ "No",
      Observation.month > 6 & Observation.month < 10 ~ "No"),
    Peak.season = factor(Peak.season))

# loop over P.Measures, for each measure grab the Part_Id.
# if this measure represents a positive malaria case, increment the Participants$Malaria.cases variable
Participants$Observations       <- rep(c(0),each=nrow(Participants)) # STILL NEED TO ADD THIS TO DERIVED VARIABLES LIST
Participants$Malaria.cases      <- rep(c(0),each=nrow(Participants))
Participants$Malaria.peak.cases <- rep(c(0),each=nrow(Participants))
for (i in 1:nrow(P.Measures)) {
  PID <- P.Measures$Participant_Id[i]   # get participant ID
  Prow <- which(Participants$Participant_Id == PID)   # get the row in participants table for this participant
  Participants$Observations[Prow] <- Participants$Observations[Prow] + 1    # increment observations count
  if (P.Measures$Malaria.diagnosis..EUPATH_0000090.[i] == "Yes"){
    Participants$Malaria.cases[Prow] <- Participants$Malaria.cases[Prow] + 1    # if malaria + yes increment case count
    if (P.Measures$Peak.season[i] == "Yes"){
      Participants$Malaria.peak.cases[Prow] <- Participants$Malaria.peak.cases[Prow] + 1
    }
  }
}
remove(i,PID,Prow) # need to refresh these values before rerunning the loop unless I want to mess things up
# now i can calculate Had.malaria!
Participants <- Participants %>% mutate(Had.malaria = case_when(
  Malaria.cases > 0 ~ "Yes",Malaria.cases == 0 ~ "No"))

# clean up new column
Participants <- Participants %>% mutate(Had.malaria = factor(Had.malaria))

# HOUSEHOLDS
# add empty columns in Households
Households$Participants <- rep(c(0),each=nrow(Households))        # total number of participants in study
Households$Observations <- rep(c(0),each=nrow(Households))        # total observations
Households$Person.years <- rep(c(0),each=nrow(Households))        # total person years of participants
Households$Malaria.cases <- rep(c(0),each=nrow(Households))       # number of malaria cases
Households$Malaria.peak.cases <- rep(c(0),each=nrow(Households))  # number of malaria cases in peak season
Households$Age.yrs.tot <- rep(c(0),each=nrow(Households))         # cumulative age of participants (to calc avg age)
# CANNOT RUN THIS LOOP TWICE or it will multiply the true values of the new columns x number of loop runs
# to reset: rerun the lines above that set the columns to 0s, then run loop once
for (i in 1:nrow(Participants)) {
  HID <- Participants$Household_Id[i]
  Hrow <- which(Households$Household_Id == HID)
  Households$Participants[Hrow]       <- Households$Participants[Hrow] + 1     # increment num of Participants
  Households$Observations[Hrow]       <- Households$Observations[Hrow] + Participants$Observations[i]
  Households$Person.years[Hrow]       <- Households$Person.years[Hrow] + Participants$Years.at.risk[i]
  Households$Malaria.cases[Hrow]      <- Households$Malaria.cases[Hrow] + Participants$Malaria.cases[i]
  Households$Malaria.peak.cases[Hrow] <- Households$Malaria.peak.cases[Hrow] + Participants$Malaria.peak.cases[i]
  Households$Age.yrs.tot[Hrow]        <- Households$Age.yrs.tot[Hrow] + Participants$Age.at.enrollment..years...EUPATH_0000120.[i]
}
remove(i,HID,Hrow) # need to refresh these values before rerunning the loop unless I want to mess things up

Households$Age.mean <- Households$Age.yrs.tot / Households$Participants
Households <- Households %>% select(-Age.yrs.tot)

# Functions                       ----
patDemo2 <- function(location){
  subset <- master %>% filter(Subcounty == location)
  nParts <- subset %>% select(PID) %>% unique() %>% nrow()
  nObs <-   subset %>% nrow()
  nObsP <-  subset %>% filter(Peak.season == "Yes") %>% nrow()
# loop to get avg age, number of male participants, pys sum, malaria cases, malaria peak cases
  pAge <- sexRatio <- nPys <- nMal <- nMalP <- 0
  for (i in 1:nObs){
    pAge <- pAge + subset$Age[i] # add age of part to counter of total age, to be averaged later
    if (subset$Sex[i] == "Male"){
      sexRatio <- sexRatio + 1
    } # increment sexRatio which will be a count of male parts for now
    nPys <- nPys + subset$Years.at.risk[i] # add pys of part to counter of total pys
    if (subset$Malaria[i] == "Yes"){
      nMal <- nMal + 1 # increment number of malaria cases
      if (subset$Peak.season[i] == "Yes"){
      nMalP <- nMalP + 1 # increment number of malaria peak cases
      }
    }
  }
  pAge <- pAge/nObs
  sexRatio <- sexRatio/(nObs-sexRatio)
  col <- c(nParts,round(pAge,digits = 3),round(sexRatio,digits = 3),round(nPys/1000,digits = 3),nObs,nObsP,nMal,nMalP)
}
houseDemographics <- function(subset){
  nhouses <- subset %>% nrow()
  # loop over subset for the rest of the variables
  modern <- rich <- middle <- poor <- totmeals <- earth <- 0
  for (i in 1:nhouses){
    if (subset$Dwelling.type..ENVO_01000744.[i] == "Modern"){modern<-modern+1}
    if (subset$Household.wealth.index..categorical..EUPATH_0000143.[i] == "Least poor"){rich<-rich+1}
    if (subset$Household.wealth.index..categorical..EUPATH_0000143.[i] == "Middle"){middle<-middle+1}
    if (subset$Household.wealth.index..categorical..EUPATH_0000143.[i] == "Poorest"){poor<-poor+1}
    if (subset$Earth.floor[i] == "Yes"){earth<-earth+1}
    totmeals <- totmeals + subset$Meals.per.day..EUPATH_0000027.[i]
  }
  meals <- totmeals/nhouses %>% round(.,digits = 3)
  col <- c(nhouses,modern,rich,middle,poor,earth,round(meals,digits = 3))
}
calcVPCs <- function(model){
  # assumes model is the output of glmer with three levels, level 2 = PID and level 3 = HID
  sig.h <- VarCorr(model)$HID[1]
  sig.p <- VarCorr(model)$PID[1]
  constant <- (pi^2)/3
  VPC.h <- round((sig.h / (sig.h + sig.p + constant)),digits = 3)
  VPC.p <- round((sig.h + sig.p) / (sig.h + sig.p + constant),digits = 3)
  VPC.i <- 1-VPC.p
  # output is a row in a future table where the columns are VPC_i, VPC_p, and VPC_h
  c(VPC.i,(VPC.p-VPC.h),VPC.h)
}
modOutputs <- function(model){
  print(summary(model))
  writeLines("\nFixed Effects Odds Ratios")
  print(fixef(model) %>% exp() %>% round(.,digits = 3))
  writeLines("\n95% Confidence Intervals for Odds Ratios")
  print(confint(model,method="Wald") %>% exp() %>% round(.,digits = 3))
  writeLines("\nTjur's Pseudo-R2")
  print(r2(model))
  writeLines("\nEstimates of VPCs")
  print(calcVPCs(model))
}
# Table Joins                     ----
# only taking a few columns from each of the Households, Participants, and P.Measures tables to keep this manageable
houseSub <- Households %>% 
  select(Household_Id,Sub.county.in.Uganda..EUPATH_0000054.,Dwelling.type..ENVO_01000744.,
         Household.wealth.index..categorical..EUPATH_0000143.,Earth.floor,Meals.per.day..EUPATH_0000027.)
partSub <- Participants %>% 
  select(Participant_Id,Household_Id,Sex..PATO_0000047.,Years.at.risk,Malaria.cases,
         Malaria.peak.cases,Had.malaria)
PartPlus <- full_join(partSub, houseSub, by = "Household_Id")
pmSub <- P.Measures %>% select(Participant_repeated_measure_Id,Participant_Id,
                               Observation.date..EUPATH_0004991.,Age..years...OBI_0001169.,
                               Malaria.diagnosis..EUPATH_0000090.,Time.since.enrollment..days...EUPATH_0000191.,
                               Peak.season,ITN.last.night..EUPATH_0000216.)
master <- full_join(pmSub, PartPlus, by = "Participant_Id")
#rename some columns
colnames(master) <- c("MeasureID","PID","Observation.date","Age","Malaria","Days.since.enrollment",
                      "Peak.season","ITN","HID","Sex","Years.at.risk","Cases","Peak.cases","Had.malaria",
                      "Subcounty","Type","Wealth","Earth.floor","Meals.per.day")
master <- master %>% select(MeasureID,PID,HID,Subcounty,Malaria,Type,Peak.season,Age,Sex,Wealth,
                            Meals.per.day,Earth.floor,Cases,Peak.cases,Years.at.risk,Had.malaria,
                            Observation.date,Days.since.enrollment,ITN)

remove(houseSub,partSub,PartPlus,pmSub)
master$Year <- format(master$Observation.date,"%Y")
master <- master %>% mutate(Year = factor(Year))

# End of Intro                    ----

# -------------------------------------------------------------------------
        # Stage I - Demographics #
# -------------------------------------------------------------------------
# Patient Demographics Table      ----
kihCol <- patDemo2("Kihihi")
walCol <- patDemo2("Walukuba")
nagCol <- patDemo2("Nagongera")
patDemLabels <- c("Participants", "Average Age", "Male:Female Ratio","Person Years (/1000)",
                  "Observations","Observations: Peak Season","Malaria Cases","Malaria Cases: Peak Season")
patDemTable <- data.frame(patDemLabels,kihCol,walCol,nagCol)
colnames(patDemTable) <- c("Patient Demographics", "Kihihi", "Walukuba", "Nagongera")
remove(patDemLabels,kihCol,walCol,nagCol)

# Household Demographics Table    ----

# Create location-specific subsets of Households
Nagongera <- Households %>% filter(Sub.county.in.Uganda..EUPATH_0000054. == "Nagongera")
Kihihi    <- Households %>% filter(Sub.county.in.Uganda..EUPATH_0000054. == "Kihihi")
Walukuba  <- Households %>% filter(Sub.county.in.Uganda..EUPATH_0000054. == "Walukuba")

kihCol <- houseDemographics(Kihihi)
walCol <- houseDemographics(Walukuba)
nagCol <- houseDemographics(Nagongera)
homeDemLabels <- c("Households", "Modern", "Wealth Index: Richest","Wealth Index: Middle",
                  "Wealth Index: Poorest","Earth Floors","Average Meals Per Day")
homeDemTable <- data.frame(homeDemLabels,kihCol,walCol,nagCol)
colnames(homeDemTable) <- c("Household Demographics", "Kihihi", "Walukuba", "Nagongera")

remove(homeDemLabels,kihCol,walCol,nagCol,Kihihi,Walukuba,Nagongera)
# -------------------------------------------------------------------------
        # Stage II - Investigating Household Type #
# -------------------------------------------------------------------------
# EXPLORATORY ANALYSIS            ----

# checking that within groupings there is still enough data with either outcome
master %>% select(Subcounty,Malaria) %>% table()
master %>% select(Type,Malaria) %>% table()

caseData <- table(master$Malaria,master$Subcounty,master$Type)
caseData <- as.data.frame(caseData)
colnames(caseData) <- c("malaria","subcounty","Type","freq")
prevData <- as.data.frame(group_by(caseData,Type,subcounty) %>% mutate(percent = freq/sum(freq)*100))
graphData <- subset(prevData, malaria == "Yes")
graphData$sample <- as.data.frame(table(master$Subcounty,master$Type))$Freq
graphData$SE <- sqrt(graphData$percent*(100-graphData$percent)/graphData$sample)

s2explorePlot <- ggplot(graphData, aes(subcounty, percent, color = Type))
s2explorePlot <- s2explorePlot +
  geom_pointrange(aes(ymin=percent-1.96*SE, ymax=percent+1.96*SE)) +
  xlab("Subcounty") +
  ylab("Percent of Observations with Positive Malaria Diagnosis")
s2explorePlot
remove(caseData,prevData,graphData)

# HOUSEHOLD MODELS (Kihihi)       ----
Kihihi <- master %>% filter(Subcounty == "Kihihi")

# adding type - PRIMARY MODEL
m5 <- glmer(Malaria ~ Age + Sex + Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m5) # custom output function created in line 287

# null model
m0 <- glmer(Malaria ~ 1 + (1 | HID/PID), family = binomial(link="logit"), data = Kihihi)
modOutputs(m0)

# Baseline Model - three levels, fixed effects included for age and sex
m1 <- glmer(Malaria ~ Age + Sex + (1 | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m1)

# adding Wealth
m2 <- glmer(Malaria ~ Age + Sex + Wealth + (1 | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m2)

# adding meals
m3 <- glmer(Malaria ~ Age + Sex + Meals.per.day + (1 | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m3)

# adding floors
m4 <- glmer(Malaria ~ Age + Sex + Earth.floor + (1 | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m4)

m6 <- glmer(Malaria ~ Age + Sex + Wealth + Earth.floor+ Meals.per.day + Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"), nAGQ = 0)
modOutputs(m6)

# starting to fit random effects
m7 <- glmer(Malaria ~ Age + Sex + Wealth + Type + (Type | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"), nAGQ = 0)
modOutputs(m7)

m8 <- glmer(Malaria ~ Age + Sex + Wealth + Type + (Age | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"), nAGQ = 0)
modOutputs(m8)

m9 <- glmer(Malaria ~ Age + Sex + Wealth + Type + (Wealth | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"), nAGQ = 0)
modOutputs(m9)

m10 <- glmer(Malaria ~ Age + Sex + Wealth + Type + (Wealth + Type | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"), nAGQ = 0)
modOutputs(m10)

m11 <- glmer(Malaria ~ Age + Sex + Wealth + Type + Type:Age + (1 | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"), nAGQ = 0)
modOutputs(m11)

m12 <- glmer(Malaria ~ Age + Sex + Wealth + Type + Type:Wealth + (1 | HID/PID), 
             family = binomial(link="logit"), data = Kihihi,
             glmerControl(optimizer = "bobyqa"), nAGQ = 0)
modOutputs(m12)

# HOUSEHOLD MODELS (Walukuba)     ----
Walukuba <-  master %>% filter(Subcounty == "Walukuba")

# including type - PRIMARY MODEL
m5 <- glmer(Malaria ~ Age + Sex + Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m5)

# null model
m0 <- glmer(Malaria ~ 1 + (1 | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m0)

# baseline model, includes age and sex
m1 <- glmer(Malaria ~ Age + Sex + (1 | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m1)

# including wealth
m2 <- glmer(Malaria ~ Age + Sex + Wealth + (1 | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m2)

# including meals
m3 <- glmer(Malaria ~ Age + Sex + Meals.per.day + (1 | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m3)

# including earth
m4 <- glmer(Malaria ~ Age + Sex + Earth.floor + (1 | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m4)                                   # Age and earth floors are significant, sex is not
fixef(m4) %>% exp() %>% round(.,digits = 3)   # Age OR: 0.981   Earth OR: 1.684
r2(m4)                                        # Conditional R2: 0.260, Marginal R2: 0.035
calcVPCs(m4)                                  # 0.766 0.070 0.164

# including all household variables
m6 <- glmer(Malaria ~ Age + Sex + Type + Earth.floor + Wealth + Meals.per.day + (1 | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m6)                                   # Age and earth are significant, all others are not
fixef(m6) %>% exp() %>% round(.,digits = 3)   # Age OR: 0.981   Earth OR: 1.901
r2(m6)                                        # Conditional R2: 0.263, Marginal R2: 0.047
calcVPCs(m6)                                  # 0.773 0.068 0.159

# leaving in Type, trying some random effects
m7 <- glmer(Malaria ~ Age + Sex + Type + (Type | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa")) # failed to converge: degenerate  Hessian with 1 negative eigenvalues
modOutputs(m7)                                   # Model is nearly unidentifiable: large eigenvalue ratio
fixef(m7) %>% exp() %>% round(.,digits = 3)   # -
r2(m7)                                        # Conditional R2: 0.260, Marginal R2: 0.022
calcVPCs(m7)                                  # 0.776 0.107 0.117

# leaving in Type, trying some random effects
m8 <- glmer(Malaria ~ Age + Sex + Type + (Age | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m8)                                   # Age is significant
fixef(m8) %>% exp() %>% round(.,digits = 3)   # Age OR: 0.979
r2(m8)                                        # Conditional R2: NA, Marginal R2: NA
calcVPCs(m8)                                  # 0.734 0.038 0.228

m9 <- glmer(Malaria ~ Age + Sex + Type + (Earth.floor | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa")) # large eigenvalue ratio
modOutputs(m9)                                   # Age is significant
fixef(m9) %>% exp() %>% round(.,digits = 3)   # Age OR: 0.980
r2(m9)                                        # Conditional R2: 0.276, Marginal R2: 0.025
calcVPCs(m9)                                  # 0.741 0.045 0.214

m10 <- glmer(Malaria ~ Age + Sex + Type + Type:Age + (1 | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m10)                                   # Age is significant
fixef(m10) %>% exp() %>% round(.,digits = 3)   # Age OR: 0.979
r2(m10)                                        # Conditional R2: 0.256, Marginal R2: 0.023
calcVPCs(m10)                                  # 0.761 0.070 0.169

m11 <- glmer(Malaria ~ Age + Sex + Type + Type:Wealth + (1 | HID/PID), 
             family = binomial(link="logit"), data = Walukuba,
             glmerControl(optimizer = "bobyqa"))
modOutputs(m11)                                   # Age is significant
fixef(m11) %>% exp() %>% round(.,digits = 3)   # Age OR: 0.980
r2(m11)                                        # Conditional R2: 0.261, Marginal R2: 0.036
calcVPCs(m11)                                  # 0.766 0.069 0.165

# HOUSEHOLD MODELS (Nagongera)    ----
Nagongera <- master %>% filter(Subcounty == "Nagongera")

# adding type - no other household vars were significant and none improved the models fit, so im not adjusting for them
m5 <- glmer(Malaria ~ Age + Sex + Type + (1 | HID/PID), 
               family = binomial(link="logit"), data = Nagongera,
               glmerControl(optimizer = "bobyqa"))
modOutputs(m5)                                   # age significant, type not significant
r2(m5)                                        # Conditional R2: 0.369, Marginal R2: 0.284
calcVPCs(m5)                                  # 0.880 0.098 0.022

# null model
m0 <- glmer(Malaria ~ 1 + (1 | HID/PID), family = binomial(link="logit"), data = Nagongera)
modOutputs(m0)                                   # Singularity... uh oh
fixef(m0) %>% exp() %>% round(.,digits = 3)   # Age OR: 0.980
confint(m0,method="Wald") %>% exp() %>% round(.,digits = 3)
r2(m0)                                        # Conditional R2: NA
calcVPCs(m0)                                  # 0.798 0.202 0.000

# baseline model
m1 <- glmer(Malaria ~ Age + Sex + (1 | HID/PID), family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m1)                                   # age significant
r2(m1)                                        # Conditional R2: 0.369, Marginal R2: 0.284
calcVPCs(m1)                                  # 0.881 0.097 0.022

# adding wealth
m2 <- glmer(Malaria ~ Age + Sex + Wealth + (1 | HID/PID), family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m2)                                   # age significant
r2(m2)                                        # Conditional R2: 0.370, Marginal R2: 0.286
calcVPCs(m2)                                  # 0.882 0.098 0.020

# adding meals
m3 <- glmer(Malaria ~ Age + Sex + Meals.per.day + (1 | HID/PID), family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m3)                                   # age significant
r2(m3)                                        # Conditional R2: 0.370, Marginal R2: 0.284
calcVPCs(m3)                                  # 0.880 0.098 0.022

# adding earth
m4 <- glmer(Malaria ~ Age + Sex + Earth.floor + (1 | HID/PID), family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m4)                                   # age significant
r2(m4)                                        # Conditional R2: 0.369, Marginal R2: 0.284
calcVPCs(m4)                                  # 0.881 0.097 0.022

# all household variables 
m6 <- glmer(Malaria ~ Age + Sex + Type + Wealth + Earth.floor + Meals.per.day + (1 | HID/PID), family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m6)                                   # age significant,
r2(m6)                                        # Conditional R2: 0.371, Marginal R2: 0.286
calcVPCs(m6)                                  # 0.884 0.089 0.027

# random slope - type
m7 <- glmer(Malaria ~ Age + Sex + Type + (Type | HID/PID), family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m7)                                   # age significant, singular fit (no variance in slopes wrt age, but singularity coming from HID level i believe)
r2(m7)                                        # Conditional R2: NA, Marginal R2: 0.654 woah
calcVPCs(m7)                                  # 0.646 0.343 0.011

# random slope - age at PID level only
m8 <- glmer(Malaria ~ Age + Sex + Type + (Age | HID/PID), family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m8)                                   # age significant, type not significant
r2(m8)                                        # Conditional R2: 0.792, Marginal R2: 0.391
calcVPCs(m8)                                  # 0.645 0.341 0.014 WOAH

# interaction between age:Type? 
m9 <- glmer(Malaria ~ Age + Sex + Type + Age:Type + (1 | HID/PID), family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m9)                                   # age significant, age:type is significant as well
r2(m9)                                        # Conditional R2: 0.372, Marginal R2: 0.287
calcVPCs(m9)                                  # 0.881 0.097 0.022

fixef(m) %>% exp() %>% round(.,digits = 3)
confint(m,method="Wald") %>% exp() %>% round(.,digits = 3)


# PLOTTING RESIDUALS              ----
# I will do this once for each subcounty using the main m5 model
# redefine the main models for use here:
kihMod <- glmer(Malaria ~ Age + Sex + Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"))
simulationOutput <- simulateResiduals(fittedModel = kihMod, plot = F)
plot(simulationOutput,title = "DHARMa residual: Kihihi")

walMod <- glmer(Malaria ~ Age + Sex + Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
simulationOutput <- simulateResiduals(fittedModel = walMod, plot = F)
plot(simulationOutput,title = "DHARMa residual: Walukuba")

nagMod <- glmer(Malaria ~ Age + Sex + Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
simulationOutput <- simulateResiduals(fittedModel = nagMod, plot = F)
plot(simulationOutput,title = "DHARMa residual: Nagongera")
remove(kihMod,walMod,nagMod,simulationOutput)

# -------------------------------------------------------------------------
        # Stage III - Investigating Season #
# -------------------------------------------------------------------------
# EXPLORATORY ANALYSIS            ----
# exploratory analysis of season - off peak season
masterOff <- master %>% subset(Peak.season == "No")
caseData <- table(masterOff$Malaria,masterOff$Subcounty,masterOff$Type)
caseData <- as.data.frame(caseData)
colnames(caseData) <- c("malaria","subcounty","Type","freq")
prevData <- as.data.frame(group_by(caseData,Type,subcounty) %>% 
                                mutate(percent = freq/sum(freq)*100))
graphData <- subset(prevData, malaria == "Yes")
graphData$sample <- as.data.frame(table(masterOff$Subcounty,masterOff$Type))$Freq
graphData$SE <- sqrt(graphData$percent*(100-graphData$percent)/graphData$sample)

s3explorePlot1 <- ggplot(graphData, aes(subcounty, percent, color = Type))
s3explorePlot1 <- s3explorePlot1 +
  geom_pointrange(aes(ymin=percent-1.96*SE, ymax=percent+1.96*SE)) +
  xlab("Subcounty") +
  ylab("Percent of Observations with Positive Malaria Diagnosis: Off Peak")
s3explorePlot1

# exploratory analysis of season - during peak season
masterOn <- master %>% subset(Peak.season == "Yes")
caseData <- table(masterOn$Malaria,masterOn$Subcounty,masterOn$Type)
caseData <- as.data.frame(caseData)
colnames(caseData) <- c("malaria","subcounty","Type","freq")
prevData <- as.data.frame(group_by(caseData,Type,subcounty) %>% 
                                mutate(percent = freq/sum(freq)*100))
graphData <- subset(prevData, malaria == "Yes")
graphData$sample <- as.data.frame(table(masterOn$Subcounty,masterOn$Type))$Freq
graphData$SE <- sqrt(graphData$percent*(100-graphData$percent)/graphData$sample)

s3explorePlot2 <- ggplot(graphData, aes(subcounty, percent, color = Type))
s3explorePlot2 <- s3explorePlot2 +
  geom_pointrange(aes(ymin=percent-1.96*SE, ymax=percent+1.96*SE)) +
  xlab("Subcounty") +
  ylab("Percent of Observations with Positive Malaria Diagnosis: Peak Season")
s3explorePlot2

# check that there's enough data in each group
masterOff %>% select(Type,Malaria) %>% table()
masterOn  %>% select(Type,Malaria) %>% table()

# clean up
remove(caseData,prevData,graphData,masterOn,masterOff)

# SEASON MODELS (Kihihi)          ----
m1 <- glmer(Malaria ~ Age + Sex + Peak.season + (1 | HID/PID), 
                  family = binomial(link="logit"), data = Kihihi,
                  glmerControl(optimizer = "bobyqa"))
modOutputs(m1)

m2 <- glmer(Malaria ~ Age + Sex + Peak.season + Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m2)

m3 <- glmer(Malaria ~ Age + Sex + Peak.season + Type + Type:Peak.season + (1 | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m3)

m4 <- glmer(Malaria ~ Age + Sex + Peak.season + Type + (Type | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m4)

m5 <- glmer(Malaria ~ Age + Sex + Peak.season + Type + (Peak.season | HID/PID), 
            family = binomial(link="logit"), data = Kihihi,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m5)

# SEASON MODELS (Walukuba)        ----
m1 <- glmer(Malaria ~ Age + Sex + Peak.season + (1 | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m1)

m2 <- glmer(Malaria ~ Age + Sex + Peak.season + Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m2)

m3 <- glmer(Malaria ~ Age + Sex + Peak.season + Type + Type:Peak.season + (1 | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m3)

m4 <- glmer(Malaria ~ Age + Sex + Peak.season + Type + (Type | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m4)

m5 <- glmer(Malaria ~ Age + Sex + Peak.season + Type + (Peak.season | HID/PID), 
            family = binomial(link="logit"), data = Walukuba,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m5)

# SEASON MODELS (Nagongera)       ----
m1 <- glmer(Malaria ~ Age + Sex + Peak.season + (1 | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m1)

m2 <- glmer(Malaria ~ Age + Sex + Peak.season + Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m2)

m3 <- glmer(Malaria ~ Age + Sex + Peak.season + Type + Type:Peak.season + (1 | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m3)

m4 <- glmer(Malaria ~ Age + Sex + Peak.season + Type + (Type | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m4)

m5 <- glmer(Malaria ~ Age + Sex + Peak.season + Type + (Peak.season | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m5)

# -------------------------------------------------------------------------
        # Stage IV - Sensitivity Analysis: IRS #
# -------------------------------------------------------------------------
# EXPLORATORY ANALYSIS ----
# pre-IRS is before Feb 1st 2015 (from Rek et al paper)
# creating Nagongera subsets (incase this isn't run directly after household analysis)
Nagongera <- master %>% filter(Subcounty == "Nagongera")
cutoff <- "2015-02-01" %>% as.Date()
preIRS <- master %>% subset(Subcounty == "Nagongera" & Observation.date < cutoff)
postIRS <- master %>% subset(Subcounty == "Nagongera" & Observation.date >= cutoff)

# checking that within groupings there is still enough data with either outcome
summary(preIRS$Malaria)
summary(postIRS$Malaria)
preIRS %>% select(Type,Malaria) %>% table()
postIRS %>% select(Type,Malaria) %>% table()

# creating a combined nagongera subset for a dotplot
Nagongera <- Nagongera %>% 
  mutate(IRS = case_when(
      Observation.date < cutoff ~ "Pre IRS",
      Observation.date >= cutoff ~ "Post IRS"),
    IRS = factor(IRS),
    IRS = fct_relevel(IRS,"Pre IRS"))

# dotplot like before
caseData <- table(Nagongera$Malaria,Nagongera$IRS,Nagongera$Type)
caseData <- as.data.frame(caseData)
colnames(caseData) <- c("malaria","IRS","Type","freq")
prevData <- as.data.frame(group_by(caseData,Type,IRS) %>% 
                               mutate(percent = freq/sum(freq)*100))
graphData <- subset(prevData, malaria == "Yes")
graphData$sample <- as.data.frame(table(Nagongera$IRS,Nagongera$Type))$Freq
graphData$SE <- sqrt(graphData$percent*(100-graphData$percent)/graphData$sample)

s4explorePlot <- ggplot(graphData, aes(IRS, percent, color = Type))
s4explorePlot <- s4explorePlot +
  geom_pointrange(aes(ymin=percent-1.96*SE, ymax=percent+1.96*SE)) +
  xlab("IRS Status") +
  ylab("Percent of Observations with Positive Malaria Diagnosis")
s4explorePlot
remove(caseData,prevData,graphData,preIRS,postIRS)

# IRS MODELS ----
# baseline model is m1 from household investigation
# adding IRS
m1 <- glmer(Malaria ~ Age + Sex + IRS + (1 | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m1)

# adding type
m2 <- glmer(Malaria ~ Age + Sex + IRS + Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m2)

m3 <- glmer(Malaria ~ Age + Sex + IRS + Peak.season + (1 | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m3)

m4 <- glmer(Malaria ~ Age + Sex + IRS + Peak.season + Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m4)

m5 <- glmer(Malaria ~ Age + Sex + IRS + Peak.season + IRS:Peak.season + (1 | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m5)

m6 <- glmer(Malaria ~ Age + Sex + IRS + Peak.season + IRS:Peak.season + Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m6)

m7 <- glmer(Malaria ~ Age + Sex + IRS + Type + IRS:Type + (1 | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m7)

m8 <- glmer(Malaria ~ Age + Sex + IRS + Type + Peak.season + (IRS | HID/PID), 
            family = binomial(link="logit"), data = Nagongera,
            glmerControl(optimizer = "bobyqa"))
modOutputs(m8) 
# ----



