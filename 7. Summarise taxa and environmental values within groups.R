##====== GRADIENT FOREST MODELS IN NEW ZEALAND WATERS ========================##

# SCRIPT 7. SUMMARISE TAXA AND ENVIRONMENTAL VALUES WITHIN GROUPS

##------ Authors: Fabrice Stephenson, Tom Brough, John Leathwick
##------ Start date : 01/07/2019
##------ End date : 01/07/2023

##============================================================================##
# This code was produced as part of a NIWA project funded by the New Zealand 
# Department of Conservation (DOC19208) to develop the New Zealand Seafloor 
# Community Classification (NZSCC) and associated biodiversity layers described
# in:
# 1.  Stephenson et al. (2022). Development of a Seafloor Community 
#     Classification for the New Zealand Region Using a Gradient Forest Approach
#     Frontiers in Marine Science 8.
# 2.  Stephenson et al. (2023). A seafloor bioregionalisation for New Zealand. 
#     Ocean & Coastal Management 242, 106688.
# 3.  Stephenson et al. (2024). Independent statistical validation of the New 
#     Zealand Seafloor Community Classification. Aquatic Conservation: Marine 
#     and Freshwater Ecosystems 34, e4114.

##============================================================================##

# DESCRIPTION: Summarise taxa (occurrence, richness, etc) and environmental values
# (min, mediam, max) within groups

# 1.  Load files and packages
# 2.  Summary of environmental values for a given classification level
# 3.  Summary of biotic values for a given classification level

####==========    1. LOAD  PACKAGES & DATA  ================================####
require(raster); require(cluster); require(devEMF); require(tidyverse);
library(data.table); library(plyr)

rm(list = ls())
imp.vars <- c("Bathy","BedDist", "BotOxy", "BotNi", "BotPhos","BotSal",
              "BotSil", "BotTemp","BPI_broad", "BPI_fine", "ChlAGrad",
              "DET", "PB555nm","SeasTDiff", "Slope", "SSTGrad", "sed.class", 
              "TC", "POCFlux","Ebed")

# Load MPI projection
MPIproj<- CRS("+proj=aea +lat_1=-30 +lat_2=-50 +lat_0=-40 +lon_0=175 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")

# load master predictor stack 
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# load 1km predictors for whole area; extract for imp.var
setwd(paste0(dir, "/Test data"))
load("Pred_1km.CMB.Source") # (1 km resolution)
Pred_1km <- Pred_1km.CMB; rm(Pred_1km.CMB)

Pred_1km$Bathy[Pred_1km$Bathy < 0] <- 0

# load Gradient forest model with appropriate classification level (i.e, as 
# defined in Script in 5)
setwd(paste0(dir, "/Classification"))
GF75 <- raster("CMB.GF75_EEZ_1km.tif")
plot (GF75)

# Load DF; Species names
setwd(paste0(dir, "/Test data"))
load("DF.source") # DF bio 
load("species.DF.source") # species list

# Load RF; Species names
load("RF.source") # RF bio
load("species.RF.source") # species list

# Load BI; Species names
load("BI_LLG.LMG.source") # BI bio
load("species.BI_LLG.LMG.source") # species list

load("BI_MMG.source") # BI bio
load("species.BI_MMG.source") # species list

load("BI_SMG.source") # BI bio
load("species.BI_SMG.source") # species list

load("BI_SSG.source") # BI bio
load("species.BI_SSG.source") # species list

# Load MA; Species names
load("MA.source") # MA bio
load("species.MA.source") # species list

####==========    2. SUMMARY OF ENVIRONMENTAL VALUES    ====================####
### ENVIRONMENTAL SUMMARY
# Tag the environmental data with GF class
GF75_env_1km <- as.data.frame(extract(GF75,Pred_1km[,c(1,2)]))
GF75_env_1km <- as.factor(GF75_env_1km[,1])
GF75_pred_1km <- cbind(GF75_env_1km, Pred_1km)

#Remove groups with less than 5 observations  
#GF75_pred_1km  <- GF75_pred_1km[!(as.numeric(GF75_pred_1km$GF75_env) %in% which(table(GF75_pred_1km$GF75_env)<5)),]  

# loop through and create means by factor classifcation for summary per group
GF75_pred_1km_MEAN <- aggregate(GF75_pred_1km[,imp.vars], list(GF75_pred_1km$GF75_env), median)
colnames(GF75_pred_1km_MEAN)[1] <- "GF_Class"

names(GF75_pred_1km)[1] <- "GF75_Group"
Group_n <- count(GF75_pred_1km$GF75_Group)
Group_n <- subset(Group_n, x!="NA")

GF75_sd <- aggregate(GF75_pred_1km[,imp.vars], list(GF75_pred_1km$GF75_Group), sd)
GF75_95CI <- cbind(GF75_sd[1], 2*(GF75_sd[c(2:21)]/sqrt(Group_n$freq)))
GF75_quant25 <- aggregate(GF75_pred_1km[,imp.vars], list(GF75_pred_1km$GF75_Group), FUN = 'quantile', probs=c(0.25))
GF75_quant75 <- aggregate(GF75_pred_1km[,imp.vars], list(GF75_pred_1km$GF75_Group), FUN = 'quantile', probs=c(0.75))

setwd(paste0(dir, "/Classification"))
write.csv(GF75_MEAN, file = "GF75_MEAN.csv")
write.csv(GF75_sd, file = "GF75_sd.csv")
write.csv(GF75_95CI, file = "GF75_95CI.csv")
write.csv(GF75_quant25, file = "GF75_quant25.csv")
write.csv(GF75_quant75, file = "GF75_quant75.csv")

#### THIS APPROACH CAN BE REPEATED FOR UNCERTAINTY LAYERS TO QUANTIFY THE 
#### MEDIAN ± 25-75% PERCENTILE DISTRIBUTION WITHIN EACH GROUP

####==========    3. SUMMARY OF BIOTIC VALUES    ===========================####
####    3.1 Demersal fish     --------------------------------------------------
# tag the Trawl data - change to factor
GF75_fac1km <- as.data.frame(extract(GF75, DF[,c(1,2)]))
GF75_fac1km <- as.factor(GF75_fac1km[,1])
GF75_samp1km.DF <- cbind(GF75_fac1km, DF)
names(GF75_samp1km.DF)[1] <- "GF75_Group"
GF75_samp.DF <- na.omit(GF75_samp1km.DF)

setwd(paste0(dir, "/Classification"))
write.csv(GF75_samp.DF, file = "GF75_samp_DF.csv")

# loop through and create means by factor classifcation for summary per group
GF75_spe <- GF75_samp.DF[,species.DF]
indx <- sapply(GF75_spe, is.factor)
GF75_spe[indx] <- lapply(GF75_spe[indx], function(x) as.numeric(as.character(x)))
GF75_preds <- GF75_samp.DF[,imp.vars]
GF75_spe <- cbind(GF75_samp.DF[1], GF75_spe)

# calculate the mean species occurrence
mean.Pres <- GF75_spe %>%   #
  group_by(GF75_Group) %>%    #Should only leave GF
  summarise_all(funs(mean), na.rm=T)

setwd(paste0(dir, "/Classification"))
write.csv(mean.Pres, file = "GF75_MEAN.Pres_DF.csv")

# count number of samples per group
FishFullGF.cnt <- as.data.frame(GF75_samp.DF[,1])
names(FishFullGF.cnt) <- "GF_grp"
FishFullGF.cnt$n.samples <- 1

mean.cnt1 <- FishFullGF.cnt %>%   #
  group_by(GF_grp) %>%    #Should only leave GF
  summarise_all(sum, na.rm=T)

# count number of cells per GF group (extent)
FishFullGF.cnt <- as.data.frame(GF75_pred_merge[,1])
names(FishFullGF.cnt) <- "GF_grp"
FishFullGF.cnt$extent <- 1

mean.cnt2 <- FishFullGF.cnt %>%   #
  group_by(GF_grp) %>%    #Should only leave GF
  summarise_all(sum, na.rm=T)

mean.cnt_F<- merge(mean.cnt1,mean.cnt2, all = TRUE)

setwd(paste0(dir, "/Classification"))
write.csv(mean.cnt_F, file = "mean.cnt_DF.csv")

# count number of raster cells in each group  (extent)
GF75_FINAL <- cbind(GF75_preds, GF75_spe)
GF75_MEAN <- aggregate(GF75_FINAL[,1:ncol(GF75_FINAL)], list(GF75_samp.DF$GF75_Group), mean)
colnames(GF75_MEAN)[1] <- "GF_Class"

setwd(paste0(dir, "/Classification"))
write.csv(GF75_MEAN, file = "GF75_MEAN.DF.csv")

####    3.2 Benthic invertebrates     ------------------------------------------
#Benthic Inverts - change for each gear type to avoid very long script
GF75_fac1km <- as.data.frame(extract(GF75, BI_LLG.LMG[,c(1,2)]))
GF75_fac1km <- as.factor(GF75_fac1km[,1])
GF75_samp1km.BI_LLG.LMG <- cbind(GF75_fac1km, BI_LLG.LMG)

names(GF75_samp1km.BI_LLG.LMG)[1] <- "GF75_Group"
GF75_samp.BI_LLG.LMG <- na.omit(GF75_samp1km.BI_LLG.LMG) #Merge locations in TS and offshore, delete NA's
setwd(paste0(dir, "/Classification"))
write.csv(GF75_samp.BI_LLG.LMG, file = "GF75_samp_BI_LLG.LMG.csv")

# loop through and create means by factor classifcation for summary per group
GF75_spe <- GF75_samp.BI_LLG.LMG[,species.BI_LLG.LMG]
indx <- sapply(GF75_spe, is.factor)
GF75_spe[indx] <- lapply(GF75_spe[indx], function(x) as.numeric(as.character(x)))
GF75_preds <- GF75_samp.BI_LLG.LMG[,imp.vars]
GF75_spe <- cbind(GF75_samp.BI_LLG.LMG[1], GF75_spe)


# calculat the mean species occurrence
mean.Pres <- GF75_spe %>%   #
  group_by(GF75_Group) %>%    #Should only leave GF
  summarise_all(funs(mean), na.rm=T)

setwd(paste0(dir, "/Classification"))
write.csv(mean.Pres, file = "GF75_MEAN.Pres_BI_LLG.LMG.csv")

# count number of samples per group
BIFullGF.cnt <- as.data.frame(GF75_samp.BI_LLG.LMG[,1])
names(BIFullGF.cnt) <- "GF_grp"
BIFullGF.cnt$n.samples <- 1

mean.cnt1 <- BIFullGF.cnt %>%   #
  group_by(GF_grp) %>%    #Should only leave GF
  summarise_all(sum, na.rm=T)

# count number of cells per GF group (extent)
BIGF.cnt <- as.data.frame(GF75_pred_merge[,1])
names(BIGF.cnt) <- "GF_grp"
BIGF.cnt$extent <- 1

mean.cnt2 <- BIGF.cnt %>%   #
  group_by(GF_grp) %>%    #Should only leave GF
  summarise_all(sum, na.rm=T)

mean.cnt_LLG.LMG<- merge(mean.cnt1,mean.cnt2, all = TRUE)
setwd(paste0(dir, "/Classification"))
write.csv(mean.cnt_LLG.LMG, file = "mean.cnt_BI_LLG.LMG.csv")

# count number of raster cells in each group  (extent)
GF75_FINAL <- cbind(GF75_preds, GF75_spe)
GF75_MEAN <- aggregate(GF75_FINAL[,1:ncol(GF75_FINAL)], list(GF75_samp.BI_SMG$GF75_Group), mean)
colnames(GF75_MEAN)[1] <- "GF_Class"

setwd(paste0(dir, "/Classification"))
write.csv(GF75_MEAN, file = "GF75_MEAN.BI_SMG.csv")

####    3.3 Reef fish     -----------------------------------------------------
#Reef Fish - ###Make sure using 250m TS classification 
rm(list = ls())

GF75_fac1km <- as.data.frame(extract(GF75, RF[,c(1,2)]))
GF75_fac1km <- as.factor(GF75_fac1km[,1])
GF75_samp1km.RF <- cbind(GF75_fac1km, RF)

names(GF75_samp1km.RF)[1] <- "GF75_Group"
GF75_samp.RF <- na.omit(GF75_samp1km.RF) #Merge locations in TS and offshore, delete NA's
setwd(paste0(dir, "/Classification"))
write.csv(GF75_samp.RF, file = "GF75_samp_RF.csv")

# loop through and create means by factor classifcation for summary per group
GF75_spe <- GF75_samp.RF[,species.RF]
indx <- sapply(GF75_spe, is.factor)
GF75_spe[indx] <- lapply(GF75_spe[indx], function(x) as.numeric(as.character(x)))
GF75_preds <- GF75_samp.RF[,imp.vars]
GF75_spe <- cbind(GF75_samp.RF[1], GF75_spe)

# calculat the mean species occurrence
mean.Pres <- GF75_spe %>%   #
  group_by(GF75_Group) %>%    #Should only leave GF
  summarise_all(funs(mean), na.rm=T)

setwd(paste0(dir, "/Classification"))
write.csv(mean.Pres, file = "GF75_MEAN.Pres_RF.csv")

# count number of samples per group
RFFullGF.cnt <- as.data.frame(GF75_samp.RF[,1])
names(RFFullGF.cnt) <- "GF_grp"
RFFullGF.cnt$n.samples <- 1

mean.cnt1 <- RFFullGF.cnt %>%   #
  group_by(GF_grp) %>%    #Should only leave GF
  summarise_all(sum, na.rm=T)

# count number of cells per GF group (extent)
RFGF.cnt <- as.data.frame(GF75_pred_merge[,1])
names(RFGF.cnt) <- "GF_grp"
RFGF.cnt$extent <- 1

mean.cnt2 <- RFGF.cnt %>%   #
  group_by(GF_grp) %>%    #Should only leave GF
  summarise_all(sum, na.rm=T)

mean.cnt_RF<- merge(mean.cnt1,mean.cnt2, all = TRUE)
setwd(paste0(dir, "/Classification"))
write.csv(mean.cnt_RF, file = "mean.cnt_RF.csv")

# count number of raster cells in each group  (extent)
GF75_FINAL <- cbind(GF75_preds, GF75_spe)
GF75_MEAN <- aggregate(GF75_FINAL[,1:ncol(GF75_FINAL)], list(GF75_samp.RF$GF75_Group), mean)
colnames(GF75_MEAN)[1] <- "GF_Class"

setwd(paste0(dir, "/Classification"))
write.csv(GF75_MEAN, file = "GF75_MEAN.RF.csv")

####    3.4 Macroalgae     -----------------------------------------------------
#Macro-algae make sure using 250m TS classification 
GF75_fac1km <- as.data.frame(extract(GF75, MA[,c(1,2)]))
GF75_fac1km <- as.factor(GF75_fac1km[,1])
GF75_samp1km.MA <- cbind(GF75_fac1km, MA)

names(GF75_samp1km.MA)[1] <- "GF75_Group"
names(GF75_samp250m.MA)[1] <- "GF75_Group"
GF75_samp.MA <- na.omit(GF75_samp1km.MA) #Merge locations in TS and offshore, delete NA's
setwd(paste0(dir, "/Classification"))
write.csv(GF75_samp.MA, file = "GF75_samp_MA.csv")

# loop through and create means by factor classifcation for summary per group
GF75_spe <- GF75_samp.MA[,species.MA]
indx <- sapply(GF75_spe, is.factor)
GF75_spe[indx] <- lapply(GF75_spe[indx], function(x) as.numeric(as.character(x)))
GF75_preds <- GF75_samp.MA[,imp.vars]
GF75_spe <- cbind(GF75_samp.MA[1], GF75_spe)

# calculat the mean species occurrence
mean.Pres <- GF75_spe %>%   #
  group_by(GF75_Group) %>%    #Should only leave GF
  summarise_all(funs(mean), na.rm=T)

setwd(paste0(dir, "/Classification"))
write.csv(mean.Pres, file = "GF75_MEAN.Pres_MA.csv")

# count number of samples per group
MAFullGF.cnt <- as.data.frame(GF75_samp.MA[,1])
names(MAFullGF.cnt) <- "GF_grp"
MAFullGF.cnt$n.samples <- 1

mean.cnt1 <- MAFullGF.cnt %>%   #
  group_by(GF_grp) %>%    #Should only leave GF
  summarise_all(sum, na.rm=T)

# count number of cells per GF group (extent)
MAGF.cnt <- as.data.frame(GF75_pred_merge[,1])
names(MAGF.cnt) <- "GF_grp"
MAGF.cnt$extent <- 1

mean.cnt2 <- MAGF.cnt %>%   #
  group_by(GF_grp) %>%    #Should only leave GF
  summarise_all(sum, na.rm=T)

mean.cnt_MA<- merge(mean.cnt1,mean.cnt2, all = TRUE)
setwd(paste0(dir, "/Classification"))
write.csv(mean.cnt_MA, file = "mean.cnt_MA.csv")

# count number of raster cells in each group  (extent)
GF75_FINAL <- cbind(GF75_preds, GF75_spe)
GF75_MEAN <- aggregate(GF75_FINAL[,1:ncol(GF75_FINAL)], list(GF75_samp.MA$GF75_Group), mean)
colnames(GF75_MEAN)[1] <- "GF_Class"

setwd(paste0(dir, "/Classification"))
write.csv(GF75_MEAN, file = "GF75_MEAN.MA.csv")