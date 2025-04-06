##====== GRADIENT FOREST MODELS IN NEW ZEALAND WATERS ========================##

# SCRIPT 4 - Summarising bootstrapped Gradient Forest models

##------ Authors: Fabrice Stephenson, Tom Brough, John Leathwick
##------ Start date : 01/07/2019
##------ End date : 01/07/2020

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

# DESCRIPTION: summarising the bootstrapped Gradient Forest models for each 
# biotic group. Given the size of the ojbects, bootstrapped models needed to be 
# saved in seperate lists

# 1.  Load final biotic and environmental dataframes 
# 2.  Extract and summarise bootstrapped Gradient Forest model outputs
# 3.  Plot and save results: compositional turnover across environmental 
#     gradients, in PCA space and geographic space
#       3.1 Demersal fish
#       3.2 Benthic invertebrates
#       3.3 Reef fish
#       3.4 Macroalgae
#       3.5 All biotic groups combined

####==========    LOAD  PACKAGES    ========================================####
library("raster");library("rstudioapi");require(dismo);require(gbm); 
require(devEMF); library(gdata); library(rgdal)
gc() # garbage collection to get rid of any residual files.

####==========    LOAD FILES AND PACKAGES   ================================####
# Load MPI projection
MPIproj<- CRS("+proj=aea +lat_1=-30 +lat_2=-50 +lat_0=-40 +lon_0=175 +x_0=0 
              +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# load master predictor stack 
setwd(paste0(dir, "/Test data"))
load("Pred_1km.CMB.Source") # (1 km resolution)
load("Pred_250m.CMB.Source") # 250 m data within the TS
load("Pred_250m.RF.Source") # 250 m data within the TS and on predicted rocky reefs

####==========    1. LOAD FINAL BIOTIC AND ENV DF   ========================####
# Save objects in folder: model object
setwd(paste0(dir, "/Test data"))
load("DF.source") # DF bio 
load("species.DF.source") # species list
load("EnvRanges.DF.source") # DF env ranges

load("RF.source") # RF bio
load("species.RF.source") # species list
load("EnvRanges.RF.source") # RF env ranges

load("BI_SSG.source") # BI_SSG bio
load("species.BI_SSG.source") # species list
load("EnvRanges.BI.source") # RF env ranges

load("BI_SMG.source") # BI_SMG bio
load("species.BI_SMG.source") # species list
load("EnvRanges.BI.source") # RF env ranges

load("BI_MMG.source") # BI_MMG bio
load("species.BI_MMG.source") # species list
load("EnvRanges.BI.source") # RF env ranges

load("BI_LLG.LMG.source") # BI_MMG bio
load("species.BI_LLG.LMG.source") # species list
load("EnvRanges.BI.source") # RF env ranges

load("MA.source") # MA bio
load("species.MA.source") # species list
load("EnvRanges.MA.source") # RF env ranges

load("imp.vars.CMB.source")
load("EnvRanges.CMB.source") # CMB env ranges

imp.vars <- c("Bathy","BedDist", "BotOxy", "BotNi", "BotPhos","BotSal",
              "BotSil", "BotTemp","BPI_broad", "BPI_fine", "ChlAGrad",
              "DET", "PB555nm","SeasTDiff", "Slope", "SSTGrad", "sed.class", 
              "TC", "POCFlux","Ebed")
n.boot <- 4

# function for changing file names from loaded data
loadRData <- function(fileName){
     #loads an RData file, and returns it
     load(fileName)
     get(ls()[ls() != "fileName"])
}

# PREPARE FOR BOOTSTRAPPING THROUGH ITERATIONS
#list of files
ls <- list("1" = "tmp1","2" = "tmp2","3" = "tmp3","4" = "tmp4", "5" = "tmp5")

####===========   2. EXTRACT AND SUMMARISE GF RESULTS   =====================####

start <- Sys.time()
a = 5 # number of bootstrapped lists from script 3.
for (a in 1:length(ls)){
     # load output from first set of bootstraps
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     tmp3 <- loadRData(paste(ls[[a]],".source", sep = ""))

     ###### DEMERSAL FISH #####
     # IMPORTANT VARIABLES
     boot_array_ImpVar.DF <- array(0, dim = c(length(imp.vars), n.boot))
     rownames(boot_array_ImpVar.DF) <- imp.vars
     
     for (i in 1:n.boot){
          boot_array_ImpVar.DF[,i] <- tmp3[[i]]$ImpVar.DF[imp.vars]
     }
    # loop through all versions to save full set of boot array
     if (a == 1){boot_array_ImpVar.DF_F <- boot_array_ImpVar.DF
     } else {boot_array_ImpVar.DF_F <- cbind(boot_array_ImpVar.DF_F, boot_array_ImpVar.DF)}
     
     # save final object 
     if (a == length(ls)){preds_influences.DF <-t(apply(boot_array_ImpVar.DF_F, 1, function(x) c(Mean = mean(x), SD = sd(x))))   #Calculate mean and standard error of the relative influecne of each preds)
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     write.csv(preds_influences.DF, file="preds.influences_DF.csv")}
     
     # SPECIES
     boot_array_SpeR2.DF <- array(0, dim = c(length(DF[,species.DF]), n.boot)) # refer back to the  full dataset
     rownames(boot_array_SpeR2.DF) <- sort(species.DF)
     
     for (i in 1:n.boot){
          boot_array_SpeR2.DF[,i] <- tmp3[[i]]$SpeR2.DF[species.DF]
     }
     
     # Number of species modelled, min, max and mean fit per GF boot
     SpeR2_sum.DF <- array(0, dim = c(n.boot, 4)) # refer back to the  full dataset
     colnames(SpeR2_sum.DF) <- c("length", "min", "mean","max")
     
     for (i in 1:n.boot){
          SpeR2_sum.DF[i,1] <- length(boot_array_SpeR2.DF[,i]) - sum(is.na(boot_array_SpeR2.DF[,i]))
          SpeR2_sum.DF[i,2] <- min(boot_array_SpeR2.DF[,i], na.rm = T)
          SpeR2_sum.DF[i,3] <- mean(boot_array_SpeR2.DF[,i], na.rm = T)
          SpeR2_sum.DF[i,4] <- max(boot_array_SpeR2.DF[,i], na.rm = T)
     }
    
     # loop through all versions to save full set of boot array
     if (a == 1){SpeR2_sum.DF_F <- SpeR2_sum.DF
     } else {SpeR2_sum.DF_F <- rbind(SpeR2_sum.DF_F, SpeR2_sum.DF)}
     # save final object 
     if (a == length(ls)){SpeR2_sum.DF_F <- round(apply(SpeR2_sum.DF_F, 2, function(x) c(Mean = mean(x), SD = sd(x))),2)
         setwd(paste(dir, "/bootstrap_results", sep = ""))
          write.csv(SpeR2_sum.DF_F, file="Sp_meanR2_DF.csv")}
     
     # species model fit and variability
     # loop through all versions to save full set of boot array
     boot_array_SpeR2.DF[is.na(boot_array_SpeR2.DF)] <- 0
     if (a == 1){boot_array_SpeR2.DF_F <- boot_array_SpeR2.DF
     } else {boot_array_SpeR2.DF_F <- cbind(boot_array_SpeR2.DF_F, boot_array_SpeR2.DF)}
     
     # save final object 
     if (a == length(ls)){speciesR2.DF <- t(apply(boot_array_SpeR2.DF_F, 1, function(x) c(Mean.R2 = mean(x), SD = sd(x))))
          setwd(paste(dir, "/bootstrap_results", sep = ""))
          write.csv(speciesR2.DF, file="Sp_R2_DF.csv")}
     
     # PREDICTED ENVIRONMENTAL TRANSFORMS
     boot_array_EnvTran.DF <- array(0,c(length(EnvRanges.DF[[1]]),length(EnvRanges.DF),n.boot))
     dimnames(boot_array_EnvTran.DF)[[2]] <- imp.vars
     dimnames(boot_array_EnvTran.DF)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")
     
     for (i in 1:n.boot){
          for (j in 1:length(EnvRanges.DF)){
               boot_array_EnvTran.DF[,j,i] <- tmp3[[i]]$EnvTran.DF[,j]
          }
     }
     
     if (a == 1){boot_array_EnvTran.DF_F <- boot_array_EnvTran.DF
     } else {boot_array_EnvTran.DF_F <- abind::abind(boot_array_EnvTran.DF_F, boot_array_EnvTran.DF, along = 3)
     }
     # save final object 
     if (a == length(ls)){setwd(paste(dir, "/bootstrap_results", sep = ""))
          save(boot_array_EnvTran.DF_F, file="boot_array_EnvTran_DF.source")}
     
     # FRED
     Fred.DF <- array(0,c(length(Pred_1km.CMB.CMB[[1]]),length(EnvRanges.DF)))
     dimnames(Fred.DF)[[2]] <- imp.vars
     
     # CALCULATE MEAN 
     for (i in 1:n.boot){Fred.DF <- Fred.DF + tmp3[[i]]$TurnOvEEZ.DF}
     Fred.DF1 <- Fred.DF/n.boot
     
     if (a == 1){Fred.DF_F <- Fred.DF1} else {Fred.DF_F <- Fred.DF_F+Fred.DF1}
     # save final object 
     if (a == length(ls)){Fred.DF_mean <- Fred.DF_F / a
          setwd(paste(dir, "/bootstrap_results", sep = ""))
          save(Fred.DF_mean, file="Pred_EEZ_DF.source")}
     
     # FRED 2 - calculate SD 
     Fred2.DF <- array(0,c(length(Pred_1km.CMB[[1]]),length(EnvRanges.DF)))
     dimnames(Fred2.DF)[[2]] <- imp.vars
     
     # Calculate the mean 
     for (i in 1:n.boot){
          if (i == 1){
               Fred2.DF <- (Fred.DF1 - tmp3[[1]]$TurnOvEEZ.DF)^2
          } else {Fred2.DF <- Fred2.DF + (Fred.DF1 - tmp3[[i]]$TurnOvEEZ.DF)^2
          }
     }
     
     Fred2.DF1 <- sqrt(Fred2.DF)
     
     if (a == 1){Fred2.DF_F <- Fred2.DF1} else {Fred2.DF_F <- Fred2.DF_F + Fred2.DF1}
     # save final object 
     if (a == length(ls)){Fred2.DF_mean <- Fred2.DF_F / a
         Fred2.DF_mean <- as.vector(apply(Fred2.DF_mean, 1, mean))
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     save(Fred2.DF_mean, file="SD_pred_EEZ_DF.source")}
     
     ###### REEF FISH #####
     # IMPORTANT VARIABLES
     boot_array_ImpVar.RF <- array(0, dim = c(length(imp.vars), n.boot))
     rownames(boot_array_ImpVar.RF) <- imp.vars
     
     for (i in 1:n.boot){
         boot_array_ImpVar.RF[,i] <- tmp3[[i]]$ImpVar.RF[imp.vars]
     }
     # loop through all versions to save full set of boot array
     if (a == 1){boot_array_ImpVar.RF_F <- boot_array_ImpVar.RF
     } else {boot_array_ImpVar.RF_F <- cbind(boot_array_ImpVar.RF_F, boot_array_ImpVar.RF)}
     
     # save final object 
     if (a == length(ls)){preds_influences.RF <-t(apply(boot_array_ImpVar.RF_F, 1, function(x) c(Mean = mean(x), SD = sd(x))))   #Calculate mean and standard error of the relative influecne of each preds)
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     write.csv(preds_influences.RF, file="preds.influences_RF.csv")}
     
     # SPECIES
     boot_array_SpeR2.RF <- array(0, dim = c(length(RF[,species.RF]), n.boot)) # refer back to the  full dataset
     rownames(boot_array_SpeR2.RF) <- sort(species.RF)
     
     for (i in 1:n.boot){
         boot_array_SpeR2.RF[,i] <- tmp3[[i]]$SpeR2.RF[species.RF]
     }
     
     # Number of species modelled, min, max and mean fit per GF boot
     SpeR2_sum.RF <- array(0, dim = c(n.boot, 4)) # refer back to the  full dataset
     colnames(SpeR2_sum.RF) <- c("length", "min", "mean","max")
     
     for (i in 1:n.boot){
         SpeR2_sum.RF[i,1] <- length(boot_array_SpeR2.RF[,i]) - sum(is.na(boot_array_SpeR2.RF[,i]))
         SpeR2_sum.RF[i,2] <- min(boot_array_SpeR2.RF[,i], na.rm = T)
         SpeR2_sum.RF[i,3] <- mean(boot_array_SpeR2.RF[,i], na.rm = T)
         SpeR2_sum.RF[i,4] <- max(boot_array_SpeR2.RF[,i], na.rm = T)
     }
     
     # loop through all versions to save full set of boot array
     if (a == 1){SpeR2_sum.RF_F <- SpeR2_sum.RF
     } else {SpeR2_sum.RF_F <- rbind(SpeR2_sum.RF_F, SpeR2_sum.RF)}
     # save final object 
     if (a == length(ls)){SpeR2_sum.RF_F <- round(apply(SpeR2_sum.RF_F, 2, function(x) c(Mean = mean(x), SD = sd(x))),2)
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     write.csv(SpeR2_sum.RF_F, file="Sp_meanR2_RF.csv")}
     
     # species model fit and variability
     # loop through all versions to save full set of boot array
     boot_array_SpeR2.RF[is.na(boot_array_SpeR2.RF)] <- 0
     if (a == 1){boot_array_SpeR2.RF_F <- boot_array_SpeR2.RF
     } else {boot_array_SpeR2.RF_F <- cbind(boot_array_SpeR2.RF_F, boot_array_SpeR2.RF)}
     
     # save final object 
     if (a == length(ls)){speciesR2.RF <- t(apply(boot_array_SpeR2.RF_F, 1, function(x) c(Mean.R2 = mean(x), SD = sd(x))))
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     write.csv(speciesR2.RF, file="Sp_R2_RF.csv")}
     
     # PREDICTED ENVIRONMENTAL TRANSFORMS
     boot_array_EnvTran.RF <- array(0,c(length(EnvRanges.RF[[1]]),length(EnvRanges.RF),n.boot))
     dimnames(boot_array_EnvTran.RF)[[2]] <- imp.vars
     dimnames(boot_array_EnvTran.RF)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")
     
     for (i in 1:n.boot){
         for (j in 1:length(EnvRanges.RF)){
             boot_array_EnvTran.RF[,j,i] <- tmp3[[i]]$EnvTran.RF[,j]
         }
     }
     
     if (a == 1){boot_array_EnvTran.RF_F <- boot_array_EnvTran.RF
     } else {boot_array_EnvTran.RF_F <- abind::abind(boot_array_EnvTran.RF_F, boot_array_EnvTran.RF, along = 3)
     }
     # save final object 
     if (a == length(ls)){setwd(paste(dir, "/bootstrap_results", sep = ""))
         save(boot_array_EnvTran.RF_F, file="boot_array_EnvTran_RF.source")}
     
     # FRED
     Fred.RF <- array(0,c(length(Pred_250m.RF[[1]]),length(EnvRanges.RF)))
     dimnames(Fred.RF)[[2]] <- imp.vars
     
     # CALCULATE MEAN 
     for (i in 1:n.boot){Fred.RF <- Fred.RF + tmp3[[i]]$TurnOvTS.RF}
     Fred.RF1 <- Fred.RF/n.boot
     
     if (a == 1){Fred.RF_F <- Fred.RF1} else {Fred.RF_F <- Fred.RF_F+Fred.RF1}
     # save final object 
     if (a == length(ls)){Fred.RF_mean <- Fred.RF_F / a
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     save(Fred.RF_mean, file="Pred_EEZ_RF.source")}
     
     # FRED 2 - calculate SD 
     Fred2.RF <- array(0,c(length(Pred_250m.RF[[1]]),length(EnvRanges.RF)))
     dimnames(Fred2.RF)[[2]] <- imp.vars
     
     # Calculate the mean 
     for (i in 1:n.boot){
         if (i == 1){
             Fred2.RF <- (Fred.RF1 - tmp3[[1]]$TurnOvTS.RF)^2
         } else {Fred2.RF <- Fred2.RF + (Fred.RF1 - tmp3[[i]]$TurnOvTS.RF)^2
         }
     }
     
     Fred2.RF1 <- sqrt(Fred2.RF)
     
     if (a == 1){Fred2.RF_F <- Fred2.RF1} else {Fred2.RF_F <- Fred2.RF_F + Fred2.RF1}
     # save final object 
     if (a == length(ls)){Fred2.RF_mean <- Fred2.RF_F / a
     Fred2.RF_mean <- as.vector(apply(Fred2.RF_mean, 1, mean))
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     save(Fred2.RF_mean, file="SD_pred_EEZ_RF.source")}
     
     
     ###### MACROALGAE #####
     # IMPORTANT VARIABLES
     boot_array_ImpVar.MA <- array(0, dim = c(length(imp.vars), n.boot))
     rownames(boot_array_ImpVar.MA) <- imp.vars
     
     for (i in 1:n.boot){
         boot_array_ImpVar.MA[,i] <- tmp3[[i]]$ImpVar.MA[imp.vars]
     }
     # loop through all versions to save full set of boot array
     if (a == 1){boot_array_ImpVar.MA_F <- boot_array_ImpVar.MA
     } else {boot_array_ImpVar.MA_F <- cbind(boot_array_ImpVar.MA_F, boot_array_ImpVar.MA)}
     
     # save final object 
     if (a == length(ls)){preds_influences.MA <-t(apply(boot_array_ImpVar.MA_F, 1, function(x) c(Mean = mean(x), SD = sd(x))))   #Calculate mean and standard error of the relative influecne of each preds)
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     write.csv(preds_influences.MA, file="preds.influences_MA.csv")}
     
     # SPECIES
     boot_array_SpeR2.MA <- array(0, dim = c(length(MA[,species.MA]), n.boot)) # refer back to the  full dataset
     rownames(boot_array_SpeR2.MA) <- sort(species.MA)
     
     for (i in 1:n.boot){
         boot_array_SpeR2.MA[,i] <- tmp3[[i]]$SpeR2.MA[species.MA]
     }
     
     # Number of species modelled, min, max and mean fit per GF boot
     SpeR2_sum.MA <- array(0, dim = c(n.boot, 4)) # refer back to the  full dataset
     colnames(SpeR2_sum.MA) <- c("length", "min", "mean","max")
     
     for (i in 1:n.boot){
         SpeR2_sum.MA[i,1] <- length(boot_array_SpeR2.MA[,i]) - sum(is.na(boot_array_SpeR2.MA[,i]))
         SpeR2_sum.MA[i,2] <- min(boot_array_SpeR2.MA[,i], na.rm = T)
         SpeR2_sum.MA[i,3] <- mean(boot_array_SpeR2.MA[,i], na.rm = T)
         SpeR2_sum.MA[i,4] <- max(boot_array_SpeR2.MA[,i], na.rm = T)
     }
     
     # loop through all versions to save full set of boot array
     if (a == 1){SpeR2_sum.MA_F <- SpeR2_sum.MA
     } else {SpeR2_sum.MA_F <- rbind(SpeR2_sum.MA_F, SpeR2_sum.MA)}
     # save final object 
     if (a == length(ls)){SpeR2_sum.MA_F <- round(apply(SpeR2_sum.MA_F, 2, function(x) c(Mean = mean(x), SD = sd(x))),2)
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     write.csv(SpeR2_sum.MA_F, file="Sp_meanR2_MA.csv")}
     
     # species model fit and variability
     # loop through all versions to save full set of boot array
     boot_array_SpeR2.MA[is.na(boot_array_SpeR2.MA)] <- 0
     if (a == 1){boot_array_SpeR2.MA_F <- boot_array_SpeR2.MA
     } else {boot_array_SpeR2.MA_F <- cbind(boot_array_SpeR2.MA_F, boot_array_SpeR2.MA)}
     
     # save final object 
     if (a == length(ls)){speciesR2.MA <- t(apply(boot_array_SpeR2.MA_F, 1, function(x) c(Mean.R2 = mean(x), SD = sd(x))))
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     write.csv(speciesR2.MA, file="Sp_R2_MA.csv")}
     
     # PREDICTED ENVIRONMENTAL TRANSFORMS
     boot_array_EnvTran.MA <- array(0,c(length(EnvRanges.MA[[1]]),length(EnvRanges.MA),n.boot))
     dimnames(boot_array_EnvTran.MA)[[2]] <- imp.vars
     dimnames(boot_array_EnvTran.MA)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")
     
     for (i in 1:n.boot){
         for (j in 1:length(EnvRanges.MA)){
             boot_array_EnvTran.MA[,j,i] <- tmp3[[i]]$EnvTran.MA[,j]
         }
     }
     
     if (a == 1){boot_array_EnvTran.MA_F <- boot_array_EnvTran.MA
     } else {boot_array_EnvTran.MA_F <- abind::abind(boot_array_EnvTran.MA_F, boot_array_EnvTran.MA, along = 3)
     }
     # save final object 
     if (a == length(ls)){setwd(paste(dir, "/bootstrap_results", sep = ""))
         save(boot_array_EnvTran.MA_F, file="boot_array_EnvTran_MA.source")}
     
     # FRED
     Fred.MA <- array(0,c(length(Pred_250m.RF[[1]]),length(EnvRanges.MA)))
     dimnames(Fred.MA)[[2]] <- imp.vars
     
     # CALCULATE MEAN 
     for (i in 1:n.boot){Fred.MA <- Fred.MA + tmp3[[i]]$TurnOvTS.MA}
     Fred.MA1 <- Fred.MA/n.boot
     
     if (a == 1){Fred.MA_F <- Fred.MA1} else {Fred.MA_F <- Fred.MA_F+Fred.MA1}
     # save final object 
     if (a == length(ls)){Fred.MA_mean <- Fred.MA_F / a
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     save(Fred.MA_mean, file="Pred_EEZ_MA.source")}
     
     # FRED 2 - calculate SD 
     Fred2.MA <- array(0,c(length(Pred_250m.RF[[1]]),length(EnvRanges.MA)))
     dimnames(Fred2.MA)[[2]] <- imp.vars
     
     # Calculate the mean 
     for (i in 1:n.boot){
         if (i == 1){
             Fred2.MA <- (Fred.MA1 - tmp3[[1]]$TurnOvTS.MA)^2
         } else {Fred2.MA <- Fred2.MA + (Fred.MA1 - tmp3[[i]]$TurnOvTS.MA)^2
         }
     }
     
     Fred2.MA1 <- sqrt(Fred2.MA)
     
     if (a == 1){Fred2.MA_F <- Fred2.MA1} else {Fred2.MA_F <- Fred2.MA_F + Fred2.MA1}
     # save final object 
     if (a == length(ls)){Fred2.MA_mean <- Fred2.MA_F / a
     Fred2.MA_mean <- as.vector(apply(Fred2.MA_mean, 1, mean))
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     save(Fred2.MA_mean, file="SD_pred_EEZ_MA.source")}
     
     ###### BENTHIC INVERTEBRATES #####
     # SPECIES
     species.BI <- c(paste("LLG.LMG.", species.BI_LLG.LMG, sep =""),
                     paste("MMG.",species.BI_MMG, sep = ""),
                     paste("SMG.",species.BI_SMG, sep = ""),
                     paste("SSG.",species.BI_SSG, sep = ""))
     boot_array_SpeR2.BI <- array(0, dim = c(length(species.BI), n.boot)) # refer back to the  full dataset
     rownames(boot_array_SpeR2.BI) <- sort(species.BI)
     
     for (i in 1:n.boot){
         boot_array_SpeR2.BI[,i] <- tmp3[[i]]$SpeR2.BI[species.BI]
     }
     
     # Number of species modelled, min, max and mean fit per GF boot
     SpeR2_sum.BI <- array(0, dim = c(n.boot, 4)) # refer back to the  full dataset
     colnames(SpeR2_sum.BI) <- c("length", "min", "mean","max")
     
     for (i in 1:n.boot){
         SpeR2_sum.BI[i,1] <- length(boot_array_SpeR2.BI[,i]) - sum(is.na(boot_array_SpeR2.BI[,i]))
         SpeR2_sum.BI[i,2] <- min(boot_array_SpeR2.BI[,i], na.rm = T)
         SpeR2_sum.BI[i,3] <- mean(boot_array_SpeR2.BI[,i], na.rm = T)
         SpeR2_sum.BI[i,4] <- max(boot_array_SpeR2.BI[,i], na.rm = T)
     }
     
     # loop through all versions to save full set of boot array
     if (a == 1){SpeR2_sum.BI_F <- SpeR2_sum.BI
     } else {SpeR2_sum.BI_F <- rbind(SpeR2_sum.BI_F, SpeR2_sum.BI)}
     # save final object 
     if (a == length(ls)){SpeR2_sum.BI_F <- round(apply(SpeR2_sum.BI_F, 2, function(x) c(Mean = mean(x), SD = sd(x))),2)
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     write.csv(SpeR2_sum.BI_F, file="Sp_meanR2_BI.csv")}
     
     # species model fit and variability
     # loop through all versions to save full set of boot array
     boot_array_SpeR2.BI[is.na(boot_array_SpeR2.BI)] <- 0
     if (a == 1){boot_array_SpeR2.BI_F <- boot_array_SpeR2.BI
     } else {boot_array_SpeR2.BI_F <- cbind(boot_array_SpeR2.BI_F, boot_array_SpeR2.BI)}
     
     # save final object 
     if (a == length(ls)){speciesR2.BI <- t(apply(boot_array_SpeR2.BI_F, 1, function(x) c(Mean.R2 = mean(x), SD = sd(x))))
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     write.csv(speciesR2.BI, file="Sp_R2_BI.csv")}
     
     # PREDICTED ENVIRONMENTAL TRANSFORMS
     boot_array_EnvTran.BI <- array(0,c(length(EnvRanges.BI[[1]]),length(EnvRanges.BI),n.boot))
     dimnames(boot_array_EnvTran.BI)[[2]] <- imp.vars
     dimnames(boot_array_EnvTran.BI)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")
     
     for (i in 1:n.boot){
         for (j in 1:length(EnvRanges.BI)){
             boot_array_EnvTran.BI[,j,i] <- tmp3[[i]]$EnvTran.BI[,j]
         }
     }
     
     if (a == 1){boot_array_EnvTran.BI_F <- boot_array_EnvTran.BI
     } else {boot_array_EnvTran.BI_F <- abind::abind(boot_array_EnvTran.BI_F, boot_array_EnvTran.BI, along = 3)
     }
     # save final object 
     if (a == length(ls)){setwd(paste(dir, "/bootstrap_results", sep = ""))
         save(boot_array_EnvTran.BI_F, file="boot_array_EnvTran_BI.source")}
     
     # FRED
     Fred.BI <- array(0,c(length(Pred_1km.CMB[[1]]),length(EnvRanges.BI)))
     dimnames(Fred.BI)[[2]] <- imp.vars
     
     # CALCULATE MEAN 
     for (i in 1:n.boot){Fred.BI <- Fred.BI + tmp3[[i]]$TurnOvEEZ.BI}
     Fred.BI1 <- Fred.BI/n.boot
     
     if (a == 1){Fred.BI_F <- Fred.BI1} else {Fred.BI_F <- Fred.BI_F+Fred.BI1}
     # save final object 
     if (a == length(ls)){Fred.BI_mean <- Fred.BI_F / a
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     save(Fred.BI_mean, file="Pred_EEZ_BI.source")}
     
     # FRED 2 - calculate SD 
     Fred2.BI <- array(0,c(length(Pred_1km.CMB[[1]]),length(EnvRanges.BI)))
     dimnames(Fred2.BI)[[2]] <- imp.vars
     
     # Calculate the mean 
     for (i in 1:n.boot){
         if (i == 1){
             Fred2.BI <- (Fred.BI1 - tmp3[[1]]$TurnOvEEZ.BI)^2
         } else {Fred2.BI <- Fred2.BI + (Fred.BI1 - tmp3[[i]]$TurnOvEEZ.BI)^2
         }
     }
     
     Fred2.BI1 <- sqrt(Fred2.BI)
     
     if (a == 1){Fred2.BI_F <- Fred2.BI1} else {Fred2.BI_F <- Fred2.BI_F + Fred2.BI1}
     # save final object 
     if (a == length(ls)){Fred2.BI_mean <- Fred2.BI_F / a
     Fred2.BI_mean <- as.vector(apply(Fred2.BI_mean, 1, mean))
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     save(Fred2.BI_mean, file="SD_pred_EEZ_BI.source")}
     
     ###### COMBINED #####
     # PREDICTED ENVIRONMENTAL TRANSFORMS
     boot_array_EnvTran.CMB <- array(0,c(length(EnvRanges.CMB[[1]]),length(EnvRanges.CMB),n.boot))
     dimnames(boot_array_EnvTran.CMB)[[2]] <- imp.vars
     dimnames(boot_array_EnvTran.CMB)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")
     
     for (i in 1:n.boot){
          for (j in 1:length(EnvRanges.CMB)){
               boot_array_EnvTran.CMB[,j,i] <- tmp3[[i]]$EnvTran.CMB[,j]
          }
     }
     
     if (a == 1){boot_array_EnvTran.CMB_F <- boot_array_EnvTran.CMB
     } else {boot_array_EnvTran.CMB_F <- abind::abind(boot_array_EnvTran.CMB_F, boot_array_EnvTran.CMB, along = 3)}
     # save final object 
     if (a == length(ls)){setwd(paste(dir, "/bootstrap_results", sep = ""))
          save(boot_array_EnvTran.CMB_F, file="boot_array_EnvTran_CMB.source")}
     
     # FRED - 1 km
     Fred.CMB <- array(0,c(length(Pred_1km.CMB[[1]]),length(EnvRanges.CMB)))
     dimnames(Fred.CMB)[[2]] <- imp.vars
     
     # CALCULATE MEAN 
     for (i in 1:n.boot){Fred.CMB <- Fred.CMB + tmp3[[i]]$TurnOvEEZ.CMB}
     Fred.CMB1 <- Fred.CMB/n.boot
     
     if (a == 1){Fred.CMB_F <- Fred.CMB1} else {Fred.CMB_F <- Fred.CMB_F+Fred.CMB1}
     # save final object 
     if (a == length(ls)){Fred.CMB_mean <- Fred.CMB_F / a
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     save(Fred.CMB_mean, file="Pred_EEZ_CMB.source")}
     
     # FRED 2 - 1 km 
     Fred2.CMB <- array(0,c(length(Pred_1km.CMB[[1]]),length(EnvRanges.CMB)))
     dimnames(Fred2.CMB)[[2]] <- imp.vars
     
     # Calculate the mean 
     for (i in 1:n.boot){
          if (i == 1){
               Fred2.CMB <- (Fred.CMB1 - tmp3[[1]]$TurnOvEEZ.CMB)^2
          } else {Fred2.CMB <- Fred2.CMB + (Fred.CMB1 - tmp3[[i]]$TurnOvEEZ.CMB)^2
          }
     }
     
     Fred2.CMB1 <- sqrt(Fred2.CMB)
     
     if (a == 1){Fred2.CMB_F <- Fred2.CMB1} else {Fred2.CMB_F <- Fred2.CMB_F + Fred2.CMB1}
     # save final object 
     if (a == length(ls)){Fred2.CMB_mean <- Fred2.CMB_F / a
         Fred2.CMB_mean <- as.vector(apply(Fred2.CMB_mean, 1, mean))
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     save(Fred2.CMB_mean, file="SD_pred_EEZ_CMB.source")}
     
     # FRED - 250m
     Fred.250.CMB <- array(0,c(length(Pred_250m.CMB[[1]]),length(EnvRanges.CMB)))
     dimnames(Fred.CMB)[[2]] <- imp.vars
     
     # CALCULATE MEAN 
     for (i in 1:n.boot){ Fred.250.CMB <-  Fred.250.CMB + tmp3[[i]]$TurnOvTS.CMB}
     Fred.250.CMB1 <- Fred.250.CMB/n.boot
     
     if (a == 1){Fred.250.CMB_F <- Fred.250.CMB1} else {Fred.250.CMB_F <- Fred.250.CMB_F + Fred.250.CMB1}
     # save final object 
     if (a == length(ls)){Fred.250.CMB_mean <- Fred.250.CMB_F / a
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     save(Fred.250.CMB_mean, file="Pred_TS_CMB.source")}
     
     # FRED 2 - 250 m
     Fred2.250.CMB <- array(0,c(length(Pred_250m.CMB[[1]]),length(EnvRanges.CMB)))
     dimnames(Fred2.CMB)[[2]] <- imp.vars
     
     # Calculate the mean 
     for (i in 1:n.boot){
          if (i == 1){
               Fred2.250.CMB <- (Fred.250.CMB1 - tmp3[[1]]$TurnOvTS.CMB)^2
          } else {Fred2.250.CMB <- Fred2.250.CMB + (Fred.250.CMB1 - tmp3[[i]]$TurnOvTS.CMB)^2
          }
     }
     
     Fred2.250.CMB1 <- sqrt(Fred2.250.CMB)
     
     if (a == 1){Fred2.250.CMB_F <- Fred2.250.CMB1} else {Fred2.250.CMB_F <- Fred2.250.CMB_F + Fred2.250.CMB1}
     # save final object 
     if (a == length(ls)){Fred2.250.CMB_mean <- Fred2.250.CMB_F / a
         Fred2.250.CMB_mean <- as.vector(apply(Fred2.250.CMB_mean, 1, mean))
     setwd(paste(dir, "/bootstrap_results", sep = ""))
     save(Fred2.250.CMB_mean, file="SD_pred_TS_CMB1.source")}
}

end <- Sys.time()
end - start

#####==========   3. PLOT RESULTS   ========================================####
#####==========   3.1 DEMERSAL FISH ========================================####
# PREDICTED ENVIRONMENTAL TRANSFORMS
setwd(paste(dir, "/bootstrap_results", sep = ""))
boot_array_EnvTran.DF <- loadRData("boot_array_EnvTran_DF.source")

# plot means + 95CI and save to file
emf(file = "EnvRanges_fullmodel_DF.emf", emfPlus = FALSE)
par(mar=c(4, 2, 1, 1))
par(mfrow=c(5,4))

# change this value to fit all env range or change tp have order of var importance
for (i in c(1:20)) {
    plot(EnvRanges.DF[,i],apply(boot_array_EnvTran.DF[,i,],1,mean), col = 2,type='l',xlab = imp.vars[i], 
         ylab = '',
         ylim = c(min(boot_array_EnvTran.DF[,,]), max(boot_array_EnvTran.DF[,,]))) #max(boot_array_EnvTran[,,])))
    # SD
    lines(EnvRanges.DF[,i],apply(boot_array_EnvTran.DF[,i,],1,mean) + apply(boot_array_EnvTran.DF[,i,],1, sd), lty = 'dashed')
    lines(EnvRanges.DF[,i],apply(boot_array_EnvTran.DF[,i,],1,mean) - apply(boot_array_EnvTran.DF[,i,],1, sd), lty = 'dashed')

    # 95% PI
    # lines(EnvRanges.DF[,i],apply(boot_array_EnvTran.DF[,i,],1, quantile, probs= c(0.1)), lty = 'dashed')
    # lines(EnvRanges.DF[,i],apply(boot_array_EnvTran.DF[,i,],1, quantile, probs= c(0.9)), lty = 'dashed')
    # 95% CI
    # lines(EnvRanges.DF[,i],apply(boot_array_EnvTran.DF[,i,],1,mean) + 2 * (sqrt(apply(boot_array_EnvTran.DF[,i,],1,var))/sqrt(5)), lty = 'dashed')
    # lines(EnvRanges.DF[,i],apply(boot_array_EnvTran.DF[,i,],1,mean) - 2 * (sqrt(apply(boot_array_EnvTran.DF[,i,],1,var))/sqrt(5)), lty = 'dashed')
    # lines(EnvRanges[,i],apply(boot_array_EnvTran[,i,],1,mean),col = 2)
    rug(quantile(DF[,imp.vars[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
    # if (i == 1) title('Sample size - 1000')
}
dev.off()

# SPATIAL PREDICTIONS OF TRANSFORMED ENVIRONMENTAL SPACE
Fred.DF <- loadRData("Pred_EEZ_DF.source")
# SETUP PCA 
PCs <- prcomp(na.omit(Fred.DF))
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
# Setup plotting parameters
nvs <- dim(PCs$rotation)[1]
vec <- c("BotOxy", "Bathy","BotSal","DET","TC") 
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <- 20
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * + 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * + 1.1

#  PLOT IN PCA SPACE 
# windows()
jpeg(filename = "PCA_DF.jpeg", width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
par(mar=c(2,2,2,2))
plot((PCs$x[, c(1,2)]), xlim = xrng, ylim = yrng, pch = ".", cex = 2, col = rgb(r, g, b, max = 255), asp = 1)
# points(PCs$rotation[!vind, c(1,2)]/scal, pch = "+")
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec,1]/scal, PCs$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec, 1]), 
     PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec, 2]), labels = vec)
title('DF GF model - EEZ')
dev.off()

# PLOT IN GEOGRAPHIC SPACE
jpeg(filename = "Spatial_Turnover_DF.jpeg", width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
par(mar=c(2,2,2,2))
plot(Pred_1km.CMB[, c("x", "y")], pch = ".", cex = 1, asp = 1, col = rgb(r, g, b, max = 255))
title('DF GF model - EEZ')
dev.off()

# save files
GF_ras_R <- rasterFromXYZ(data.frame(x = Pred_1km.CMB[, c("x")], 
                                     y = Pred_1km.CMB[, c("y")], 
                                     z = r),
                          crs = MPIproj)
GF_ras_G <- rasterFromXYZ(data.frame(x = Pred_1km.CMB[, c("x")], 
                                     y = Pred_1km.CMB[, c("y")], 
                                     z = g),
                          crs = MPIproj)
GF_ras_B <- rasterFromXYZ(data.frame(x = Pred_1km.CMB[, c("x")], 
                                     y = Pred_1km.CMB[, c("y")],
                                     z = b),
                          crs = MPIproj)

### SD 
Fred3.DF <- loadRData("SD_pred_EEZ_DF.source")
SD.ras.DF <- rasterFromXYZ(data.frame(x = Pred_1km.CMB[, c("x")],
                                      y = Pred_1km.CMB[, c("y")],
                                      z = Fred3.DF),
                           crs = MPIproj)

# write the geotiff
setwd(paste(dir, "/bootstrap_results/rasters", sep = ""))
writeRaster(GF_ras_R,"DF_R.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_G,"DF_G.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_B,"DF_B.tif","GTiff", overwrite=TRUE)
writeRaster(SD.ras.DF, file = "SD.ras_DF.tif", overwrite = T)

#####==========   3.2 BENTHIC INVERTEBRATES     ============================####
# PREDICTED ENVIRONMENTAL TRANSFORMS
setwd(paste(dir, "/bootstrap_results", sep = ""))
boot_array_EnvTran.BI <- loadRData("boot_array_EnvTran_BI.source")

# plot means + 95CI and save to file
emf(file = "EnvRanges_fullmodel_BI.emf", emfPlus = FALSE)
par(mar=c(4, 2, 1, 1))
par(mfrow=c(5,4))

# change this value to fit all env range or change tp have order of var importance
for (i in c(1:20)) {
    plot(EnvRanges.BI[,i],apply(boot_array_EnvTran.BI[,i,],1,mean), col = 2,type='l',xlab = imp.vars[i], 
         ylab = '',
         ylim = c(min(boot_array_EnvTran.BI[,,]), max(boot_array_EnvTran.BI[,,]))) #max(boot_array_EnvTran[,,])))
    # SD
    lines(EnvRanges.BI[,i],apply(boot_array_EnvTran.BI[,i,],1,mean) + apply(boot_array_EnvTran.BI[,i,],1, sd), lty = 'dashed')
    lines(EnvRanges.BI[,i],apply(boot_array_EnvTran.BI[,i,],1,mean) - apply(boot_array_EnvTran.BI[,i,],1, sd), lty = 'dashed')

    # 95% PI
    # lines(EnvRanges.BI[,i],apply(boot_array_EnvTran.BI[,i,],1, quantile, probs= c(0.1)), lty = 'dashed')
    # lines(EnvRanges.BI[,i],apply(boot_array_EnvTran.BI[,i,],1, quantile, probs= c(0.9)), lty = 'dashed')
    # 95% CI
    # lines(EnvRanges.BI[,i],apply(boot_array_EnvTran.BI[,i,],1,mean) + 2 * (sqrt(apply(boot_array_EnvTran.BI[,i,],1,var))/sqrt(5)), lty = 'dashed')
    # lines(EnvRanges.BI[,i],apply(boot_array_EnvTran.BI[,i,],1,mean) - 2 * (sqrt(apply(boot_array_EnvTran.BI[,i,],1,var))/sqrt(5)), lty = 'dashed')
    # lines(EnvRanges[,i],apply(boot_array_EnvTran[,i,],1,mean),col = 2)
    # rug(quantile([,imp.vars[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
    # if (i == 1) title('Sample size - 1000')
}
dev.off()

# SPATIAL PREDICTIONS OF TRANSFORMED ENVIRONMENTAL SPACE
Fred.BI <- loadRData("Pred_EEZ_BI.source")
# SETUP PCA 
PCs <- prcomp(na.omit(Fred.BI))
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
# Setup plotting parameters
nvs <- dim(PCs$rotation)[1]
vec <- c("BotOxy","BPI_broad", "DET","TC","Ebed") 
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <- 20
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * + 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * + 1.1

#  PLOT IN PCA SPACE 
# windows()
jpeg(filename = "PCA_BI.jpeg", width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
par(mar=c(2,2,2,2))
plot((PCs$x[, c(1,2)]), xlim = xrng, ylim = yrng, pch = ".", cex = 2, col = rgb(r, g, b, max = 255), asp = 1)
# points(PCs$rotation[!vind, c(1,2)]/scal, pch = "+")
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec,1]/scal, PCs$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec, 1]), 
     PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec, 2]), labels = vec)
title('BI GF model - EEZ')
dev.off()

# PLOT IN GEOGRAPHIC SPACE
jpeg(filename = "Spatial_Turnover_BI.jpeg", width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
par(mar=c(2,2,2,2))
plot(Pred_1km.CMB[, c("x", "y")], pch = ".", cex = 1, asp = 1, col = rgb(r, g, b, max = 255))
title('BI GF model - EEZ')
dev.off()

# save files
GF_ras_R <- rasterFromXYZ(data.frame(x = Pred_1km.CMB[, c("x")], 
                                     y = Pred_1km.CMB[, c("y")], 
                                     z = r),
                          crs = MPIproj)
GF_ras_G <- rasterFromXYZ(data.frame(x = Pred_1km.CMB[, c("x")], 
                                     y = Pred_1km.CMB[, c("y")], 
                                     z = g),
                          crs = MPIproj)
GF_ras_B <- rasterFromXYZ(data.frame(x = Pred_1km.CMB[, c("x")], 
                                     y = Pred_1km.CMB[, c("y")],
                                     z = b),
                          crs = MPIproj)

Fred3.BI <- loadRData("SD_pred_EEZ_BI.source")
SD.ras.BI <- rasterFromXYZ(data.frame(x = Pred_1km.CMB[, c("x")],
                                      y = Pred_1km.CMB[, c("y")],
                                      z = Fred3.BI),
                           crs = MPIproj)

# write the geotiff
setwd(paste(dir, "/bootstrap_results/rasters", sep = ""))
writeRaster(GF_ras_R,"BI_R.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_G,"BI_G.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_B,"BI_B.tif","GTiff", overwrite=TRUE)
### SD 
writeRaster(SD.ras.BI, file = "SD.ras_BI.tif", overwrite = T)

#####==========   3.3 REEF FISH   ==========================================####
# PREDICTED ENVIRONMENTAL TRANSFORMS
setwd(paste(dir, "/bootstrap_results", sep = ""))
boot_array_EnvTran.RF <- loadRData("boot_array_EnvTran_RF.source")

# plot means + 95CI and save to file
emf(file = "EnvRanges_fullmodel_RF.emf", emfPlus = FALSE)
par(mar=c(4, 2, 1, 1))
par(mfrow=c(5,4))

# change this value to fit all env range or change tp have order of var importance
for (i in c(1:20)) {
    plot(EnvRanges.RF[,i],apply(boot_array_EnvTran.RF[,i,],1,mean), col = 2,type='l',xlab = imp.vars[i], 
         ylab = '',
         ylim = c(min(boot_array_EnvTran.RF[,,]), max(boot_array_EnvTran.RF[,,]))) #max(boot_array_EnvTran[,,])))
    # SD
    lines(EnvRanges.RF[,i],apply(boot_array_EnvTran.RF[,i,],1,mean) + apply(boot_array_EnvTran.RF[,i,],1, sd), lty = 'dashed')
    lines(EnvRanges.RF[,i],apply(boot_array_EnvTran.RF[,i,],1,mean) - apply(boot_array_EnvTran.RF[,i,],1, sd), lty = 'dashed')

    # 95% PI
    # lines(EnvRanges.RF[,i],apply(boot_array_EnvTran.RF[,i,],1, quantile, probs= c(0.05)), lty = 'dashed')
    # lines(EnvRanges.RF[,i],apply(boot_array_EnvTran.RF[,i,],1, quantile, probs= c(0.95)), lty = 'dashed')
    # 95% CI
    # lines(EnvRanges.RF[,i],apply(boot_array_EnvTran.RF[,i,],1,mean) + 2 * (sqrt(apply(boot_array_EnvTran.RF[,i,],1,var))/sqrt(5)), lty = 'dashed')
    # lines(EnvRanges.RF[,i],apply(boot_array_EnvTran.RF[,i,],1,mean) - 2 * (sqrt(apply(boot_array_EnvTran.RF[,i,],1,var))/sqrt(5)), lty = 'dashed')
    # lines(EnvRanges[,i],apply(boot_array_EnvTran[,i,],1,mean),col = 2)
    rug(quantile(RF[,imp.vars[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
    # if (i == 1) title('Sample size - 1000')
}
dev.off()

# SPATIAL PREDICTIONS OF TRANSFORMED ENVIRONMENTAL SPACE
Fred.RF <- loadRData("Pred_EEZ_RF.source")
# SETUP PCA 
PCs <- prcomp(na.omit(Fred.RF))
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
# Setup plotting parameters
nvs <- dim(PCs$rotation)[1]
vec <- c("BotOxy", "BotSal", "BotTemp", "DET", "PB555nm","POCFlux")
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <- 20
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * + 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * + 1.1

#  PLOT IN PCA SPACE 
# windows()
jpeg(filename = "PCA_RF.jpeg", width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
par(mar=c(2,2,2,2))
plot((PCs$x[, c(1,2)]), xlim = xrng, ylim = yrng, pch = ".", cex = 2, col = rgb(r, g, b, max = 255), asp = 1)
# points(PCs$rotation[!vind, c(1,2)]/scal, pch = "+")
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec,1]/scal, PCs$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec, 1]), 
     PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec, 2]), labels = vec)
title('RF GF model - EEZ')
dev.off()

# PLOT IN GEOGRAPHIC SPACE
jpeg(filename = "Spatial_Turnover_RF.jpeg", width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
par(mar=c(2,2,2,2))
plot(Pred_250m.RF[, c("x", "y")], pch = ".", cex = 1, asp = 1, col = rgb(r, g, b, max = 255))
title('RF GF model - EEZ')
dev.off()

# save files
GF_ras_R <- rasterFromXYZ(data.frame(x = Pred_250m.RF[, c("x")], 
                                     y = Pred_250m.RF[, c("y")], 
                                     z = r),
                          crs = MPIproj)
GF_ras_G <- rasterFromXYZ(data.frame(x = Pred_250m.RF[, c("x")], 
                                     y = Pred_250m.RF[, c("y")], 
                                     z = g),
                          crs = MPIproj)
GF_ras_B <- rasterFromXYZ(data.frame(x = Pred_250m.RF[, c("x")], 
                                     y = Pred_250m.RF[, c("y")],
                                     z = b),
                          crs = MPIproj)
### SD 
Fred3.RF <- loadRData("SD_pred_EEZ_RF.source")
SD.ras.RF <- rasterFromXYZ(data.frame(x = Pred_250m.RF[, c("x")],
                                      y = Pred_250m.RF[, c("y")],
                                      z = Fred3.RF),
                           crs = MPIproj)

# write the geotiff
setwd(paste(dir, "/bootstrap_results/rasters", sep = ""))
# setwd(paste(dir, "/bootstrap_results/final_runs/Rasters", sep = ""))
writeRaster(GF_ras_R,"RF_R.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_G,"RF_G.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_B,"RF_B.tif","GTiff", overwrite=TRUE)
writeRaster(SD.ras.RF, file = "SD.ras_RF.tif", overwrite = T)

#####==========   3.4 MACROALGAE    ========================================####
# PREDICTED ENVIRONMENTAL TRANSFORMS
setwd(paste(dir, "/bootstrap_results", sep = ""))
# setwd(paste(dir, "/bootstrap_results", sep = ""))
boot_array_EnvTran.MA <- loadRData("boot_array_EnvTran_MA.source")

# plot means + 95CI and save to file
emf(file = "EnvRanges_fullmodel_MA.emf", emfPlus = FALSE)
par(mar=c(4, 2, 1, 1))
par(mfrow=c(5,4))

# change this value to fit all env range or change tp have order of var importance
for (i in c(1:20)) {
    plot(EnvRanges.MA[,i],apply(boot_array_EnvTran.MA[,i,],1,mean), col = 2,type='l',xlab = imp.vars[i], 
         ylab = '',
         ylim = c(min(boot_array_EnvTran.MA[,,]), max(boot_array_EnvTran.MA[,,]))) #max(boot_array_EnvTran[,,])))
   # SD
    lines(EnvRanges.MA[,i],apply(boot_array_EnvTran.MA[,i,],1,mean) + apply(boot_array_EnvTran.MA[,i,],1, sd), lty = 'dashed')
    lines(EnvRanges.MA[,i],apply(boot_array_EnvTran.MA[,i,],1,mean) - apply(boot_array_EnvTran.MA[,i,],1, sd), lty = 'dashed')

    # 95% PI
    # lines(EnvRanges.MA[,i],apply(boot_array_EnvTran.MA[,i,],1, quantile, probs= c(0.05)), lty = 'dashed')
    # lines(EnvRanges.MA[,i],apply(boot_array_EnvTran.MA[,i,],1, quantile, probs= c(0.95)), lty = 'dashed')
    # 95% CI
    # lines(EnvRanges.MA[,i],apply(boot_array_EnvTran.MA[,i,],1,mean) + 2 * (sqrt(apply(boot_array_EnvTran.MA[,i,],1,var))/sqrt(5)), lty = 'dashed')
    # lines(EnvRanges.MA[,i],apply(boot_array_EnvTran.MA[,i,],1,mean) - 2 * (sqrt(apply(boot_array_EnvTran.MA[,i,],1,var))/sqrt(5)), lty = 'dashed')
    # lines(EnvRanges[,i],apply(boot_array_EnvTran[,i,],1,mean),col = 2)
    rug(quantile(MA[,imp.vars[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
    # if (i == 1) title('Sample size - 1000')
}
dev.off()

# SPATIAL PREDICTIONS OF TRANSFORMED ENVIRONMENTAL SPACE
Fred.MA <- loadRData("Pred_EEZ_MA.source")
# SETUP PCA 
PCs <- prcomp(na.omit(Fred.MA))
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
# Setup plotting parameters
nvs <- dim(PCs$rotation)[1]
vec <- c("BotSil","DET","BotNi","TC")
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <- 20
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * + 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * + 1.1

#  PLOT IN PCA SPACE 
# windows()
jpeg(filename = "PCA_MA.jpeg", width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
par(mar=c(2,2,2,2))
plot((PCs$x[, c(1,2)]), xlim = xrng, ylim = yrng, pch = ".", cex = 2, col = rgb(r, g, b, max = 255), asp = 1)
# points(PCs$rotation[!vind, c(1,2)]/scal, pch = "+")
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec,1]/scal, PCs$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec, 1]), 
     PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec, 2]), labels = vec)
title('MA GF model - EEZ')
dev.off()

# PLOT IN GEOGRAPHIC SPACE
jpeg(filename = "Spatial_Turnover_MA.jpeg", width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
par(mar=c(2,2,2,2))
plot(Pred_250m.RF[, c("x", "y")], pch = ".", cex = 1, asp = 1, col = rgb(r, g, b, max = 255))
title('MA GF model - EEZ')
dev.off()

# save files
GF_ras_R <- rasterFromXYZ(data.frame(x = Pred_250m.RF[, c("x")], 
                                     y = Pred_250m.RF[, c("y")], 
                                     z = r),
                          crs = MPIproj)
GF_ras_G <- rasterFromXYZ(data.frame(x = Pred_250m.RF[, c("x")], 
                                     y = Pred_250m.RF[, c("y")], 
                                     z = g),
                          crs = MPIproj)
GF_ras_B <- rasterFromXYZ(data.frame(x = Pred_250m.RF[, c("x")], 
                                     y = Pred_250m.RF[, c("y")],
                                     z = b),
                          crs = MPIproj)
### SD 
Fred3.MA <- loadRData("SD_pred_EEZ_MA.source")
SD.ras.MA <- rasterFromXYZ(data.frame(x = Pred_250m.RF[, c("x")],
                                      y = Pred_250m.RF[, c("y")],
                                      z = Fred3.MA),
                           crs = MPIproj)

# write the geotiff
setwd(paste(dir, "/bootstrap_results/rasters", sep = ""))
# setwd(paste(dir, "/bootstrap_results/final_runs/Rasters", sep = ""))
writeRaster(GF_ras_R,"MA_R.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_G,"MA_G.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_B,"MA_B.tif","GTiff", overwrite=TRUE)
writeRaster(SD.ras.MA, file = "SD.ras_MA.tif", overwrite = T)

#####==========   3.5 ALL BIOTIC GROUPS COMBINED    ========================####
# PREDICTED ENVIRONMENTAL TRANSFORMS
setwd(paste(dir, "/bootstrap_results", sep = ""))
boot_array_EnvTran.CMB <- loadRData("boot_array_EnvTran_CMB.source")

# Environemtnal predictor contributions
tmp <- array(0,c(20,2))
dimnames(tmp)[[1]] <- imp.vars
dimnames(tmp)[[2]] <- c("Mean","SD")

for (i in c(1:20)){
     tmp[i,1] <- apply(boot_array_EnvTran.CMB[,i,],1,mean)[[200]]
     tmp[i,2] <- apply(boot_array_EnvTran.CMB[,i,],1, sd)[[200]]
}

# plot means + 95CI and save to file
emf(file = "EnvRanges_fullmodel_CMB.emf", emfPlus = FALSE)
par(mar=c(4, 2, 1, 1))
par(mfrow=c(5,4))

# change this value to fit all env range or change tp have order of var importance
for (i in c(1:20)) {
    plot(EnvRanges.CMB[,i],apply(boot_array_EnvTran.CMB[,i,],1,mean), col = 2,type='l',xlab = imp.vars[i], 
         ylab = '',
         ylim = c(min(boot_array_EnvTran.CMB[,,]), max(boot_array_EnvTran.CMB[,,]))) #max(boot_array_EnvTran[,,])))
    # SD
    lines(EnvRanges.CMB[,i],apply(boot_array_EnvTran.CMB[,i,],1,mean) + apply(boot_array_EnvTran.CMB[,i,],1, sd), lty = 'dashed')
    lines(EnvRanges.CMB[,i],apply(boot_array_EnvTran.CMB[,i,],1,mean) - apply(boot_array_EnvTran.CMB[,i,],1, sd), lty = 'dashed')

    # 95% PI
    # lines(EnvRanges.CMB[,i],apply(boot_array_EnvTran.CMB[,i,],1, quantile, probs= c(0.05)), lty = 'dashed')
    # lines(EnvRanges.CMB[,i],apply(boot_array_EnvTran.CMB[,i,],1, quantile, probs= c(0.95)), lty = 'dashed')
    # 95% CI
    # lines(EnvRanges.CMB[,i],apply(boot_array_EnvTran.CMB[,i,],1,mean) + 2 * (sqrt(apply(boot_array_EnvTran.CMB[,i,],1,var))/sqrt(5)), lty = 'dashed')
    # lines(EnvRanges.CMB[,i],apply(boot_array_EnvTran.CMB[,i,],1,mean) - 2 * (sqrt(apply(boot_array_EnvTran.CMB[,i,],1,var))/sqrt(5)), lty = 'dashed')
    # lines(EnvRanges[,i],apply(boot_array_EnvTran[,i,],1,mean),col = 2)
    # rug(quantile(CMB[,imp.vars[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
    # if (i == 1) title('Sample size - 1000')
}
dev.off()

# plot means for each predictor
emf(file = "EnvRanges_allM_CMB.emf", emfPlus = FALSE)
par(mar=c(4, 2, 1, 1))
par(mfrow=c(5,4))

# change this value to fit all env range or change tp have order of var importance
for (i in c(1:20)) {
    plot(EnvRanges.CMB[,imp.vars[i]],apply(boot_array_EnvTran.CMB[,imp.vars[i],],1,mean), col = "black",type='l',xlab = imp.vars[i], 
         ylab = '',
         ylim = c(min(boot_array_EnvTran.CMB[,,]), max(boot_array_EnvTran.CMB[,,]))) #max(boot_array_EnvTran[,,])))
    lines(EnvRanges.DF[,imp.vars[i]],apply(boot_array_EnvTran.DF[,imp.vars[i],],1,mean), col = "blue")
    lines(EnvRanges.RF[,imp.vars[i]],apply(boot_array_EnvTran.RF[,imp.vars[i],],1,mean), col = "coral")
    lines(EnvRanges.BI[,imp.vars[i]],apply(boot_array_EnvTran.BI[,imp.vars[i],],1,mean), col = "yellow")
    lines(EnvRanges.MA[,imp.vars[i]],apply(boot_array_EnvTran.MA[,imp.vars[i],],1,mean), col = "green")
    # for (j in c(1:length(EnvRanges.RF[]))){
    # if(names(EnvRanges.RF[j]) == imp.vars[i]){
    # lines(EnvRanges.RF[,imp.vars[i]],apply(boot_array_EnvTran.RF[,imp.vars[i],],1,mean), col = "coral")}
    # else {next}
    # }
}
dev.off()

# SPATIAL PREDICTIONS OF TRANSFORMED ENVIRONMENTAL SPACE
Fred.CMB <- loadRData("Pred_EEZ_CMB.source")
# SETUP PCA 
PCs <- prcomp(na.omit(Fred.CMB))
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
# Setup plotting parameters
nvs <- dim(PCs$rotation)[1]
vec <- c("Bathy", "BotOxy","BotSil","BPI_broad","PB555nm","TC")
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <- 20
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * + 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * + 1.1

#  PLOT IN PCA SPACE 
# windows()
jpeg(filename = "PCA_CMB.jpeg", width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
par(mar=c(2,2,2,2))
plot((PCs$x[, c(1,2)]), xlim = xrng, ylim = yrng, pch = ".", cex = 2, col = rgb(r, g, b, max = 255), asp = 1)
# points(PCs$rotation[!vind, c(1,2)]/scal, pch = "+")
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec,1]/scal, PCs$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec, 1]), 
     PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec, 2]), labels = vec)
title('CMB GF model - EEZ')
dev.off()

# PLOT IN GEOGRAPHIC SPACE
jpeg(filename = "Spatial_Turnover_CMB.jpeg", width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
par(mar=c(2,2,2,2))
plot(Pred_1km.CMB[, c("x", "y")], pch = ".", cex = 1, asp = 1, col = rgb(r, g, b, max = 255))
title('CMB GF model - EEZ')
dev.off()

# save files
GF_ras_R <- rasterFromXYZ(data.frame(x = Pred_1km.CMB[, c("x")], 
                                     y = Pred_1km.CMB[, c("y")], 
                                     z = r),
                          crs = MPIproj)
GF_ras_G <- rasterFromXYZ(data.frame(x = Pred_1km.CMB[, c("x")], 
                                     y = Pred_1km.CMB[, c("y")], 
                                     z = g),
                          crs = MPIproj)
GF_ras_B <- rasterFromXYZ(data.frame(x = Pred_1km.CMB[, c("x")], 
                                     y = Pred_1km.CMB[, c("y")],
                                     z = b),
                          crs = MPIproj)
### SD 
Fred3.CMB <- loadRData("SD_pred_EEZ_CMB.source")
SD.ras.CMB <- rasterFromXYZ(data.frame(x = Pred_1km.CMB[, c("x")],
                                       y = Pred_1km.CMB[, c("y")],
                                       z = Fred3.CMB),
                            crs = MPIproj)


# write the geotiff
setwd(paste(dir, "/bootstrap_results/rasters", sep = ""))
# setwd(paste(dir, "/bootstrap_results/final_runs/Rasters", sep = ""))
writeRaster(GF_ras_R,"CMB_R.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_G,"CMB_G.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_B,"CMB_B.tif","GTiff", overwrite=TRUE)
writeRaster(SD.ras.CMB, file = "SD.ras_CMB.tif", overwrite = T)






