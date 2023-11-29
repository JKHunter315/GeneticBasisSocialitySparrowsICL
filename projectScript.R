                                        ### PROJECT SCRIPT ###

rm(list=ls())
setwd("~/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Summer research project/MscProject")
# libraries
require(igraph)
require(MCMCglmm)
require(lme4)
require(dplyr)

# Load arrival time  function
source("Data/ArrivalNetwork.R")

                                   ### SOCIAL NETWORK ANALYSIS ###

### Load social data:
raw_socialData <- read.csv("Data/social_dat0821.csv")
str(raw_socialData)
head(raw_socialData)

raw_socialData %>% group_by(birdID) %>% summarise(n=n()) # 500 individuals


    ## 2016-17 winter months

# get 2016 data
sDat16 <- raw_socialData %>% filter(year == 2016 & quarter == "NB")
# remove month 1
sDat16 <- sDat16 %>% filter(month != 1)
# add 2017 data
sDat17 <- raw_socialData %>% filter(year == 2017 & quarter == "NB")
sDat167 <- rbind(sDat16, sDat17)
rm(sDat16, sDat17)

sDat167 %>% group_by(month) %>% summarise(n=n()) # 3 months, all winter
sDat167 %>% group_by(week) %>% summarise(n=n()) # 14 weeks
sDat167 %>% group_by(birdID) %>% summarise(n=n()) # 69 individuals

# make weeks consecutive
sDat167$week <- as.factor(sDat167$week)
levels(sDat167$week) <- 1:length(unique(sDat167$week))
str(sDat167)

                          ## 2015-16 winter months

# get 2015 data
sDat15 <- raw_socialData %>% filter(year == 2015 & quarter == "NB")
# remove month 1
sDat15 <- sDat15 %>% filter(month != 1)

# add 2016 data, month 1 only
sDat16_2 <- raw_socialData %>% filter(year == 2016 & quarter == "NB")
sDat16_2 <- sDat16_2 %>% filter(month == 1)
sDat156 <- rbind(sDat15, sDat16_2)
rm(sDat15, sDat16_2)

sDat156 %>% group_by(month) %>% summarise(n=n()) # 3 months
sDat156 %>% group_by(week) %>% summarise(n=n()) # 15 weeks 
sDat156 %>% group_by(birdID) %>% summarise(n=n()) # 118 individuals

##

# use dataset separated into weeks, combine into both years, z-scale and trivariate model

# add different years - as many as possible

sDat20 <- raw_socialData %>% filter(year == 2020 & quarter == "NB") # no 2019, only 2020
sDat20 %>% group_by(month) %>% summarise(n = n()) # only november
sDat20 %>% group_by(week) %>% summarise(n = n())
sDat13 <- raw_socialData %>% filter(year == 2013 & quarter == "NB") 
sDat13 %>% group_by(week) %>% summarise(n = n())
sDat13 %>% group_by(month) %>% summarise(n = n()) # jan and november - no 2014 so its fine

# convert weeks to consecutive
count <- 1
sDat13$week <- as.factor(sDat13$week)
levels(sDat13$week) <- 1:length(unique(sDat13$week))
levels(sDat13$week) 
count <- 3
sDat156$week <- as.factor(sDat156$week)
levels(sDat156$week) <- (1:length(unique(sDat156$week))) + count
levels(sDat156$week)
count <- 18
sDat167$week <- as.factor(sDat167$week)
levels(sDat167$week) <- (1:length(unique(sDat167$week))) + count
levels(sDat167$week)
count <- 32
sDat20$week <- as.factor(sDat20$week)
levels(sDat20$week) <- (1:length(unique(sDat20$week))) + count
levels(sDat20$week)

# combine into 1 df
allYrDF <- rbind(sDat13, sDat156, sDat167, sDat20)
allYrDF$animal <- allYrDF$birdID
str(allYrDF)
1:length(unique(allYrDF$week))
unique(allYrDF$week)

# create social network for each week

# initialise starting DF
allYr_Socs <- data.frame(birdID = numeric(), degrees = numeric(),
                           strength = numeric(), betweenness = numeric(),
                           density = numeric(), week = numeric())

for (i in 1:length(unique(allYrDF$week))) {
  list <- BuildNetworkArrive(allYrDF[allYrDF$week == i,], t = 150)
  graph <- graph_from_data_frame(list, directed = F)
  degrees <- data.frame(igraph::degree(graph))
  strength <- data.frame(igraph::strength(graph))
  between <- data.frame(betweenness(graph))
  density <- data.frame(igraph::graph.density(graph))

  tempDF <- data.frame(birdID = row.names(degrees), degrees = degrees,
                       strength = strength, betweenness = between, 
                       density = density, week = i)
  colnames(tempDF) <- c("birdID","degrees", "strength", "betweenness",
                        "density", "week")
  allYr_Socs <- rbind(allYr_Socs, tempDF)
}

plot.igraph(graph)

# remove rows with NA
allYr_Socs <- allYr_Socs %>% filter(birdID != "NA")
# rename rows to 1:length of df
row.names(allYr_Socs) <- 1:nrow(allYr_Socs)

allYr_Socs$birdID <- as.factor(allYr_Socs$birdID) # change birdID to factor
allYr_Socs$animal <- allYr_Socs$birdID
allYr_Socs$z.deg <- scale(allYr_Socs$degrees)
allYr_Socs$z.bet <- scale(allYr_Socs$betweenness)
length(unique(allYr_Socs_clean$birdID))

# set priors
priors <- list(R = list(V = diag(3), nu = 0.002), 
               G = list(G1 = list(V = diag(3), nu = 0.002), 
                        G2 = list(V = diag(3), nu = 0.002)))

## genetics part ##

library(pedigree)
library(pedantics)
pd <- read.csv("Data/EP_Final")

# create new df with just columns for id and parents
pd3 <- pd[,1:3]
# convert to MCMC friendly names: animal
colnames(pd3) <- c("animal", "dam", "sire")

# pedigree with dams and sires not in offspring added as founders
pd_add <- pd3
invalid_sires <- !(pd3$sire %in% pd3$animal)
invalid_dams <- !(pd3$dam %in% pd3$animal)
missing_sires <- data.frame(animal = unique(pd3$sire[invalid_sires]), 
                            dam = NA, sire = NA)
missing_dams <- data.frame(animal = unique(pd3$dam[invalid_dams]), 
                           dam = NA, sire = NA)
missing_parents <- rbind(missing_dams, missing_sires)
pd_add <- rbind(pd_add, missing_parents)
pd_add <- fixPedigree(pd_add) # order ped
colnames(pd_add) <- c("animal", "dam", "sire") # rename to animal for mcmcglmm
rm(missing_dams, missing_sires, missing_parents)
rm(pd)
rm(pd3)

row.names(pd_add) <- 1:nrow(pd_add)

# plot pedigree
drawPedigree(pd_add, dots="y")

# indicate which birds are measured in social network
measuredPed <- pd_add
measuredPed$Measured <- ifelse(measuredPed$animal %in% allYrDF$birdID == "TRUE",
                               1,0)
head(measuredPed)

drawPedigree(measuredPed, dots = "y", dat = measuredPed$Measured,
             retain = "informative")
dev.off()


# no covariance in G-structure - try model that has restricted G-structure
# covariance in residuals so keep residuals unstructured
# keep animal unstructured to see covariance

# remove rows with birdIds not in pedigree: throws error
missing_animals <- setdiff(unique(allYr_Socs$animal), unique(pd_add$animal))
allYr_Socs_clean <- subset(allYr_Socs, !animal %in% missing_animals)

allMCMC <- MCMCglmm(c(z.deg, strength, z.bet) ~ trait + density, 
                    random = ~us(trait):animal + idh(trait):birdID, 
                    rcov = ~us(trait):units, 
                    family = c("gaussian", "gaussian", "gaussian"),
                    data = allYr_Socs_clean, ped = pd_add, nitt = 1000000, 
                    burnin = 300000, prior = priors)
summary(allMCMC)
autocorr(allMCMC$VCV)
plot(allMCMC$VCV)



                                ### ANALYSIS: matrices ###

mu <-  c(mean(allMCMC[["Sol"]][ , "(Intercept)"]),
         mean(allMCMC[["Sol"]][ , "traitstrength"]),
         mean(allMCMC[["Sol"]][ , "traitz.bet"]))
G1 <-
  matrix(c(mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.animal"]),
           mean(allMCMC[["VCV"]][ , "traitstrength:traitz.deg.animal"]),
           0,
           mean(allMCMC[["VCV"]][ , "traitz.deg:traitstrength.animal"]),
           mean(allMCMC[["VCV"]][ , "traitstrength:traitstrength.animal"]),
           0,
           0,
           0,
           mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.bet.animal"])),
         ncol = 3)
R <-
  matrix(c(mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.units"]),
           mean(allMCMC[["VCV"]][ , "traitstrength:traitz.deg.units"]),
           mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.deg.units"]),
           mean(allMCMC[["VCV"]][ , "traitz.deg:traitstrength.units"]),
           mean(allMCMC[["VCV"]][ , "traitstrength:traitstrength.units"]),
           mean(allMCMC[["VCV"]][ , "traitz.bet:traitstrength.units"]),
           mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.bet.units"]),
           mean(allMCMC[["VCV"]][ , "traitstrength:traitz.bet.units"]),
           mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.bet.units"])),
         ncol = 3)

G2 <- c(mean(allMCMC[["VCV"]][ , "traitz.deg.birdID"]),
        mean(allMCMC[["VCV"]][ , "traitstrength.birdID"]),
        mean(allMCMC[["VCV"]][ , "traitz.bet.birdID"]))
  
P <- G1 + G2 + R

G1 / P
G2 + G1 /P
G2 / P

GCorr <- allMCMC$VCV[,"traitstrength:traitz.deg.animal"]/
  (sqrt(allMCMC$VCV[,"traitz.deg:traitz.deg.animal"])*
     sqrt(allMCMC$VCV[,"traitstrength:traitstrength.animal"]))
posterior.mode(GCorr)

## percentage variance
totalVar <- mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.animal"]) +
  mean(allMCMC[["VCV"]][ , "traitstrength:traitz.deg.animal"]) +
  mean(allMCMC[["VCV"]][ , "traitstrength:traitstrength.animal"]) +
  mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.bet.animal"]) + 
  mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.units"]) + 
  mean(allMCMC[["VCV"]][ , "traitstrength:traitz.deg.units"]) + 
  mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.deg.units"]) + 
  mean(allMCMC[["VCV"]][ , "traitstrength:traitstrength.units"]) +
  mean(allMCMC[["VCV"]][ , "traitz.bet:traitstrength.units"]) + 
  mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.bet.units"]) + 
  mean(allMCMC[["VCV"]][ , "traitz.deg.birdID"]) +
  mean(allMCMC[["VCV"]][ , "traitstrength.birdID"]) + 
  mean(allMCMC[["VCV"]][ , "traitz.bet.birdID"])

mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.animal"])/totalVar * 100
HPDinterval(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.animal"]/totalVar * 100)
mean(allMCMC[["VCV"]][ , "traitstrength:traitstrength.animal"])/totalVar * 100
HPDinterval(allMCMC[["VCV"]][ , "traitstrength:traitstrength.animal"]/totalVar * 100)
mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.bet.animal"])/totalVar * 100
HPDinterval(allMCMC[["VCV"]][ , "traitz.bet:traitz.bet.animal"]/totalVar * 100)

mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.units"])/totalVar * 100
HPDinterval(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.units"])/totalVar * 100
mean(allMCMC[["VCV"]][ , "traitstrength:traitstrength.units"])/totalVar * 100
HPDinterval(allMCMC[["VCV"]][ , "traitstrength:traitstrength.units"])/totalVar * 100
mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.bet.units"])/totalVar * 100
HPDinterval(allMCMC[["VCV"]][ , "traitz.bet:traitz.bet.units"])/totalVar * 100

mean(allMCMC[["VCV"]][ , "traitz.deg.birdID"])/totalVar * 100
HPDinterval(allMCMC[["VCV"]][ , "traitz.deg.birdID"])/totalVar * 100
mean(allMCMC[["VCV"]][ , "traitstrength.birdID"])/totalVar * 100
HPDinterval(allMCMC[["VCV"]][ , "traitstrength.birdID"])/totalVar * 100
mean(allMCMC[["VCV"]][ , "traitz.bet.birdID"])/totalVar * 100
HPDinterval(allMCMC[["VCV"]][ , "traitz.bet.birdID"])/totalVar * 100

# genetic correlations
  mean(allMCMC[["VCV"]][ , "traitstrength:traitz.deg.animal"]) / 
  (sqrt(mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.animal"])) + 
     sqrt(mean(allMCMC[["VCV"]][ , "traitstrength:traitstrength.animal"])))
  HPDinterval(allMCMC[["VCV"]][ , "traitstrength:traitz.deg.animal"]) / 
    (sqrt(mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.animal"])) + 
       sqrt(mean(allMCMC[["VCV"]][ , "traitstrength:traitstrength.animal"])))

  mean(allMCMC[["VCV"]][ , "traitstrength:traitz.deg.units"]) / 
  (sqrt(mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.units"])) + 
      sqrt(mean(allMCMC[["VCV"]][ , "traitstrength:traitstrength.units"])))
  HPDinterval(allMCMC[["VCV"]][ , "traitstrength:traitz.deg.units"]) / 
    (sqrt(mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.units"])) + 
       sqrt(mean(allMCMC[["VCV"]][ , "traitstrength:traitstrength.units"])))
  
mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.deg.units"]) /
  (sqrt(mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.units"])) + 
     sqrt(mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.bet.units"])))
HPDinterval(allMCMC[["VCV"]][ , "traitz.bet:traitz.deg.units"]) /
  (sqrt(mean(allMCMC[["VCV"]][ , "traitz.deg:traitz.deg.units"])) + 
     sqrt(mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.bet.units"])))

mean(allMCMC[["VCV"]][ , "traitz.bet:traitstrength.units"]) / 
  (sqrt(mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.bet.units"])) + 
     sqrt(mean(allMCMC[["VCV"]][ , "traitstrength:traitstrength.units"])))
HPDinterval(allMCMC[["VCV"]][ , "traitz.bet:traitstrength.units"]) / 
  (sqrt(mean(allMCMC[["VCV"]][ , "traitz.bet:traitz.bet.units"])) + 
     sqrt(mean(allMCMC[["VCV"]][ , "traitstrength:traitstrength.units"])))
