##====== GRADIENT FOREST MODELS IN NEW ZEALAND WATERS ========================##

# SCRIPT 3 - Bootstrapping and aggregating Gradient Forest models

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

# DESCRIPTION: bootstrapping and aggregation of Gradient Forest models for each 
# biotic group. Key decisions were number of bootstraps and number of trees
# to allow models to run in parallel without running out memory. 

# 1. Prepare final biotic and enironmental data using demersal fish as example
# 2. Load final dataframes for modelling
# 3. Bootstrap and combine Gradient Forest models for all biotic groups

####==========    LOAD  PACKAGES    ========================================####
library("raster");library("rstudioapi");require(dismo);require(gbm); 
require(devEMF); require(extendedForest); require(gradientForest)
library(parallel); library(foreach); library(bigstatsr)
gc() # garbage collection to get rid of any residual files.

####==========    LOAD FILES AND PACKAGES   ================================####
# For each taxa, provide Presence/Absence dataframe; dataframe of predictors; 
# Load MPI projection
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(dir, "/Test data"))

MPIproj <- CRS("+proj=aea +lat_1=-30 +lat_2=-50 +lat_0=-40 +lon_0=175 +x_0=0 
              +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")

# load("Pred_1km.CMB.Source") # spatial data for prediction at 1 km resolution for exclusive economic zone
# load("Pred_250m.CMB.Source") # spatial data for prediction at 250 m resolution for territorial sea
# load("DF.Source") # species data for tuning model
# load("BI.Source")
# load("RF.Source")
# load("MA.Source")
# 
# ####========    1. EXAMPLE DATA - DEMERSAL FISH (DF)    ==================####
# # imp vars for demersal fish
# imp.vars.DF <- c("Bathy","BedDist", "BotOxy", "BotNi", "BotPhos","BotSal",
#                  "BotSil", "BotTemp","BPI_broad", "BPI_fine","ChlAGrad",
#                  "DET", "PB555nm","SeasTDiff", "Slope", "SSTGrad", "sed.class", 
#                  "TC", "POCFlux", "Ebed")
# 
# # remove excess preds in DF and env DF
# Pred_1km.DF <- cbind(Pred_1km[,c(1:2)],Pred_1km[, imp.vars.DF])
# Pred_1km.DF <- na.omit(Pred_1km.DF)
# Pred_250m.DF <- cbind(Pred_250m[,c(1:2)],Pred_250m[, imp.vars.DF])
# Pred_250m.DF <- na.omit(Pred_250m.DF)
# 
# DF <- cbind(DF[,c(1:2)],DF[, imp.vars.DF],DF[33:ncol(DF)])
# head(DF[c(1:22)])
# # 20 if not including the new change
# species.DF <- names(DF[23:ncol(DF)])
# 
# # create environemtnal gradient files for each taxa
# PredMins <- apply(Pred_1km[,c(3:ncol(Pred_1km))],2,min)
# PredMaxs <- apply(Pred_1km[,c(3:ncol(Pred_1km))],2,max)
# 
# EnvRanges.DF <- as.data.frame(array(0,c(200,length(imp.vars.DF))))
# names(EnvRanges.DF) <- imp.vars.DF
# 
# for (i in c(1:length(imp.vars.DF))) {
#      EnvRanges.DF[,i] <- seq(PredMins[imp.vars.DF[i]],PredMaxs[imp.vars.DF[i]],length = 200)
# }
# 
# EnvRanges.DF$Bathy <- abs(EnvRanges.DF$Bathy)
# 
# #### Save final dataframes for modelling
# # Save objects in folder: model object
# setwd(paste0(dir, "/Test data"))
# save(DF, file = "DF.source") # DF bio 
# save(imp.vars.DF, file = "imp.vars.DF.source")
# save(species.DF, file = "species.DF.source") # species list
# save(Pred_1km.DF, file ="Pred_1km.DF.source") # DF Env
# save(Pred_250m, file ="Pred_250m.DF.source") # DF Env
# save(EnvRanges.DF,file = "EnvRanges.DF.source") # DF env ranges
# 
# Repeat for other taxa groups  ------------------------------------------------

# Environmental data            ------------------------------------------------
# # create imp.vars.CMB (combination of all imp.vars)
# imp.vars.CMB <- imp.vars.DF
# 
# # create pred_1km.CMB 
# Pred_1km.CMB <- cbind(Pred_1km[,c(1:2)],Pred_1km[, imp.vars.CMB])
# Pred_1km.CMB <- na.omit(Pred_1km.CMB)
# 
# # create pred_250m.CMB 
# Pred_250m.CMB <- cbind(Pred_250m[,c(1:2)],Pred_250m[, imp.vars.CMB])
# Pred_250m.CMB <- na.omit(Pred_250m.CMB)
# 
# # create EnvRanges.CMB 
# EnvRanges.CMB <- as.data.frame(array(0,c(200,length(imp.vars.CMB))))
# names(EnvRanges.CMB) <- imp.vars.CMB
# 
# for (i in c(1:length(imp.vars.CMB))) {
#      EnvRanges.CMB[,i] <- seq(PredMins[imp.vars.CMB[i]],PredMaxs[imp.vars.CMB[i]],length = 200)
# }
# 
# #### Save final dataframes for modelling
# # Save objects in folder: model object
# setwd(paste0(dir, "/Test data"))
# save(imp.vars.CMB, file = "imp.vars.CMB.source")
# save(Pred_1km.CMB, file ="Pred_1km.CMB.source") # CMB Env
# save(Pred_250m.CMB, file ="Pred_250m.CMB.source") # CMB Env
# save(EnvRanges.CMB,file = "EnvRanges.CMB.source") # CMB env ranges

####==========    2. LOAD FINAL DATAFRAMES FOR MODELLING    ================####
# Save objects in folder: model object
setwd(paste0(dir, "/Test data"))
load("DF.source") # DF bio 
load("imp.vars.DF.source")
load("species.DF.source") # species list
load("Pred_1km.DF.source") # DF Env
load("Pred_250m.DF.source") # DF Env
load("EnvRanges.DF.source") # DF env ranges

load("RF.source") # RF bio
load("imp.vars.RF.source")
load("species.RF.source") # species list
# load("Pred_1km.RF.source") # RF Env
load("Pred_250m.RF.source") # RF Env
load("EnvRanges.RF.source") # RF env ranges

load("BI_SSG.source") # BI_SSG bio
load("imp.vars.BI.source")
load("species.BI_SSG.source") # species list
load("Pred_1km.BI.source") # RF Env
#load("Pred_250m.RF.source") # RF Env
load("EnvRanges.BI.source") # RF env ranges

load("BI_SMG.source") # BI_SMG bio
load("imp.vars.BI.source")
load("species.BI_SMG.source") # species list
load("Pred_1km.BI.source") # RF Env
#load("Pred_250m.RF.source") # RF Env
load("EnvRanges.BI.source") # RF env ranges

load("BI_MMG.source") # BI_MMG bio
load("imp.vars.BI.source")
load("species.BI_MMG.source") # species list
load("Pred_1km.BI.source") # RF Env
#load("Pred_250m.RF.source") # RF Env
load("EnvRanges.BI.source") # RF env ranges

load("BI_LLG.LMG.source") # BI_MMG bio
load("imp.vars.BI.source")
load("species.BI_LLG.LMG.source") # species list
load("Pred_1km.BI.source") # RF Env
#load("Pred_250m.RF.source") # RF Env
load("EnvRanges.BI.source") # RF env ranges

load("MA.source") # MA bio
load("imp.vars.MA.source")
load("species.MA.source") # species list
# load("Pred_1km.RF.source") # RF Env
load("Pred_250m.MA.source") # RF Env
load("EnvRanges.MA.source") # RF env ranges

load("imp.vars.CMB.source")
load("Pred_1km.CMB.source") # CMB Env
load("Pred_250m.CMB.source") # CMB Env
load("EnvRanges.CMB.source") # CMB env ranges

####===========   3. BOOTSTRAP AND COMBINE GRANDIENT FOREST MODELS    ======####

detectCores()
n.boot <- 4 # number of boots - test this first to the number of cores minus 4

# SETTING UP FILE STRUCTURE FOR SAVING OUTPUTS FROM PARRALEL LOOP
# 1) Variables of importance (list of var, percent contrib)
# 2) Species and R2 (number of species modelled, min, mean and max R2)
# 3) Environmental transforms (for samples)
# 4) Spatial prediction of turnover for EEZ at 1km grid resolution
# 5) Spatial prediction of turnover for TS at 250m grid resolution
multiResultClass <- function(ImpVar.DF=NULL,SpeR2.DF=NULL, EnvTran.DF=NULL, TurnOvEEZ.DF=NULL,
                             ImpVar.RF=NULL,SpeR2.RF=NULL, EnvTran.RF=NULL, TurnOvTS.RF=NULL,
                             ImpVar.BI=NULL,SpeR2.BI=NULL, EnvTran.BI=NULL, TurnOvEEZ.BI=NULL,
                             ImpVar.MA=NULL,SpeR2.MA=NULL, EnvTran.MA=NULL, TurnOvTS.MA=NULL,
                             EnvTran.CMB=NULL, TurnOvEEZ.CMB=NULL, TurnOvTS.CMB=NULL) {me <- list(
                               ImpVar.DF = ImpVar.DF, SpeR2.DF = SpeR2.DF, EnvTran.DF = EnvTran.DF, TurnOvEEZ.DF =  TurnOvEEZ.DF,
                               ImpVar.RF = ImpVar.RF, SpeR2.RF = SpeR2.RF, EnvTran.RF = EnvTran.RF, TurnOvTS.RF =  TurnOvTS.RF,
                               ImpVar.BI = ImpVar.BI, SpeR2.BI = SpeR2.BI, EnvTran.BI = EnvTran.BI, TurnOvEEZ.BI =  TurnOvEEZ.BI,
                               ImpVar.MA = ImpVar.MA, SpeR2.MA = SpeR2.MA, EnvTran.MA = EnvTran.MA, TurnOvTS.MA =  TurnOvTS.MA,
                               EnvTran.CMB = EnvTran.CMB, TurnOvEEZ.CMB =  TurnOvEEZ.CMB, TurnOvTS.CMB =  TurnOvTS.CMB
                             )
     
## Set the name for the class
class(me) <- append(class(me),"multiResultClass")
return(me)
}

# MAIN PARRALLEL LOOP
cl <- parallel::makeCluster(n.boot) # number of cores used - ALWAYS 1 LESS THAN AVAILABLE
# showConnections()
doParallel::registerDoParallel(cl)
packages <- c("extendedForest", "gradientForest") # need to call packages within each parralel loop

start <- Sys.time() # recording the time

tmp5 <- foreach(i = 1:n.boot) %dopar% {
  
  # change this 250 for final runs (but takes a lot longer)
  n.tree = 100
  
     # DEMERSAL FISH
     nSites <- 5000
     lev <- floor(log2(nSites * 0.368/2))
     # create training dataset
     train_ind <- sample(seq_len(nrow(DF)), size = nSites) # index of rows for 5000 samples
     DF_train <- na.omit(DF[train_ind, ])
     
     # Model
     lapply(packages, require, character.only = TRUE)
     DF.gf <- gradientForest(DF_train,
                             predictor.vars = imp.vars.DF,
                             response.vars = species.DF,
                             ntree = n.tree, 
                             transform = NULL,
                             compact = T,
                             nbin = 401,
                             maxLevel = lev,
                             corr.threshold = 0.5,
                             trace = T)
  result <- multiResultClass()
  # Imp of envir variables (sorted into alphabetical order)
  result$ImpVar.DF <- sort(round(DF.gf$overall.imp,5))
  # Species model fits 
  result$SpeR2.DF <- DF.gf$result
  # Env transforms
  result$EnvTran.DF <- predict(DF.gf, EnvRanges.DF, extrap=F) # CHECK  IF THIS NEEDS TO BE FOR MEANS OF ALL OTHER ENV PREDS
  result$TurnOvEEZ.DF <- predict(DF.gf, Pred_1km.DF[,imp.vars.DF], extrap=F)
  # result$TurnOvTS.DF <- predict(DF.gf, Pred_250m.DF[,imp.vars.DF], extrap=F)
  
  # BENTHIC INVERTEBRATES
  # GF for each geartype
  # LLG.LMG
  nSites <- 5000
  lev <- floor(log2(nSites * 0.368/2))
  # create training dataset
  train_ind <- sample(seq_len(nrow(BI_LLG.LMG)), size = nSites) # index of rows for 5000 samples
  BI_LLG.LMG_train <- na.omit(BI_LLG.LMG[train_ind, ])
  # Model
  lapply(packages, require, character.only = TRUE)
  BI_LLG.LMG.gf <- gradientForest(BI_LLG.LMG_train,
                                  predictor.vars = imp.vars.BI,
                                  response.vars = species.BI_LLG.LMG,
                                  ntree = n.tree, 
                                  transform = NULL,
                                  compact = T,
                                  nbin = 401,
                                  maxLevel = lev,
                                  corr.threshold = 0.5,
                                  trace = T)
  # SMG
  nSites <- nrow(BI_SMG)
  lev <- floor(log2(nSites * 0.368/2))
  # create training dataset
  train_ind <- sample(seq_len(nrow(BI_SMG)), size = nSites*0.75, replace = F) # index of rows for 3000 samples
  BI_SMG_train <- na.omit(BI_SMG[train_ind, ])
  # Model
  lapply(packages, require, character.only = TRUE)
  BI_SMG.gf <- gradientForest(BI_SMG_train,
                              predictor.vars = imp.vars.BI,
                              response.vars = species.BI_SMG,
                              ntree = n.tree, 
                              transform = NULL,
                              compact = T, 
                              nbin = 401, 
                              maxLevel = lev,
                              corr.threshold = 0.5,
                              trace = T)
  
  # MMG
  nSites <- nrow(BI_MMG)
  lev <- floor(log2(nSites * 0.368/2))
  # create training dataset
  train_ind <- sample(seq_len(nrow(BI_MMG)), size = nSites*0.75, replace = F) # index of rows for 3000 samples
  BI_MMG_train <- na.omit(BI_MMG[train_ind, ])
  # Model
  lapply(packages, require, character.only = TRUE)
  BI_MMG.gf <- gradientForest(BI_MMG_train,
                              predictor.vars = imp.vars.BI,
                              response.vars = species.BI_MMG,
                              ntree = n.tree, 
                              transform = NULL,
                              compact = T, 
                              nbin = 401, 
                              maxLevel = lev,
                              corr.threshold = 0.5,
                              trace = T)
  
  # SSG
  nSites <- nrow(BI_SSG)
  lev <- floor(log2(nSites * 0.368/2))
  # create training dataset
  train_ind <- sample(seq_len(nrow(BI_SSG)), size = nSites*0.75, replace = F) # index of rows for 3000 samples
  BI_SSG_train <- na.omit(BI_SSG[train_ind, ])
  # Model
  lapply(packages, require, character.only = TRUE)
  BI_SSG.gf <- gradientForest(BI_SSG_train,
                              predictor.vars = imp.vars.BI,
                              response.vars = species.BI_SSG,
                              ntree = n.tree, 
                              transform = NULL,
                              compact = T, 
                              nbin = 401, 
                              maxLevel = lev,
                              corr.threshold = 0.5,
                              trace = T)
  ##### combine BI GFs
  CBI.gf <- combinedGradientForest(LLG.LMG = BI_LLG.LMG.gf, SMG = BI_SMG.gf, MMG = BI_MMG.gf,
                                   SSG = BI_SSG.gf, nbin = 401, 
                                   method=2,
                                   standardize = "before")

  # Imp of envir variables (sorted into alphabetical order)
  # result$ImpVar.BI <- sort(round(CBI.gf$overall.imp,5))
  # Species model fits 
  result$SpeR2.BI <- CBI.gf$rsq
  # Env transforms
  result$EnvTran.BI <- predict(CBI.gf, EnvRanges.BI, extrap=F) # CHECK  IF THIS NEEDS TO BE FOR MEANS OF ALL OTHER ENV PREDS
  result$TurnOvEEZ.BI <- predict(CBI.gf, Pred_1km.BI[,imp.vars.BI], extrap=F)
  
  # MACRO ALGAE
  nSites <- nrow(MA)
  lev <- floor(log2(nSites * 0.368/2))
  # create training dataset
  train_ind <- sample(seq_len(nrow(MA)), size = nSites*0.75, replace = F) # index of rows for 3000 samples
  MA_train <- na.omit(MA[train_ind, ])
  # Model
  lapply(packages, require, character.only = TRUE)
  MA.gf <- gradientForest(MA_train,
                          predictor.vars = imp.vars.MA,
                          response.vars = species.MA,
                          ntree = n.tree, 
                          transform = NULL,
                          compact = T, 
                          nbin = 401, 
                          maxLevel = lev,
                          corr.threshold = 0.5,
                          trace = T)
  
  # Imp of envir variables (sorted into alphabetical order)
  result$ImpVar.MA <- sort(round(MA.gf$overall.imp,5))
  # Species model fits 
  result$SpeR2.MA <- MA.gf$result
  # Env transforms
  result$EnvTran.MA <- predict(MA.gf, EnvRanges.MA[,imp.vars.MA], extrap=F) # CHECK  IF THIS NEEDS TO BE FOR MEANS OF ALL OTHER ENV PREDS
  result$TurnOvTS.MA <- predict(MA.gf, Pred_250m.RF[,imp.vars.MA], extrap=F)
  
  # Reef FISH
  nSites <- nrow(RF)
  lev <- floor(log2(nSites * 0.368/2))
  # create training dataset
  train_ind <- sample(seq_len(nrow(RF)), size = nSites*0.75, replace = F) # index of rows for 3000 samples
  RF_train <- na.omit(RF[train_ind, ])
  # Model
  lapply(packages, require, character.only = TRUE)
  RF.gf <- gradientForest(RF_train,
                          predictor.vars = imp.vars.RF,
                          response.vars = species.RF,
                          ntree = n.tree, 
                          transform = NULL,
                          compact = T, 
                          nbin = 401, 
                          maxLevel = lev,
                          corr.threshold = 0.5,
                          trace = T)
  
  # Imp of envir variables (sorted into alphabetical order)
  result$ImpVar.RF <- sort(round(RF.gf$overall.imp,5))
  # Species model fits 
  result$SpeR2.RF <- RF.gf$result
  # Env transforms
  result$EnvTran.RF <- predict(RF.gf, EnvRanges.RF[,imp.vars.RF], extrap=F) # CHECK  IF THIS NEEDS TO BE FOR MEANS OF ALL OTHER ENV PREDS
  result$TurnOvTS.RF <- predict(RF.gf, Pred_250m.RF[,imp.vars.RF], extrap=F)
  
  # COMBINE MODELS
  CMB.gf <- combinedGradientForest(DF = DF.gf, RF = RF.gf, 
                                   BI_LLG.LMG = BI_LLG.LMG.gf, BI_MMG = BI_MMG.gf,
                                   BI_SMG = BI_SMG.gf, BI_SSG = BI_SSG.gf,
                                   MA = MA.gf, nbin = 401, 
                                   method=2,
                                   standardize = "before")
  
  # save turnover and spatial prediction
  result$EnvTran.CMB <- predict(CMB.gf, EnvRanges.CMB, extrap=F) # CHECK  IF THIS NEEDS TO BE FOR MEANS OF ALL OTHER ENV PREDS
  result$TurnOvEEZ.CMB <- predict(CMB.gf, Pred_1km.CMB[,imp.vars.CMB], extrap=F)
  result$TurnOvTS.CMB <- predict(CMB.gf, Pred_250m.CMB[,imp.vars.CMB], extrap=F)
  
  # Return results files
  return(result)
}

parallel::stopCluster(cl)

end <- Sys.time()
end - start

# 5000 DF; 339 RF; 5000 BI LLG; 1105 BI MMG; 1660 BI SMG; 312 BI SMG; 3041 MA
# 10 trees bootstrapped 4 times using 4 cores; predicted to EEZ data & TS: 
# 4 mins & file size = 1 Gb 
# 100 trees bootstrapped 4 times using 4 cores; predicted to EEZ data & TS: 
# 30 mins & file size = 10 Gb 
object.size(tmp5)

setwd(paste0(dir, "/bootstrap_results"))
save(tmp5, file="tmp5.source") # 
# load("tmp3.source")

# Based on size of objects, number of cores, and memory limits of computer, this 
# approach might need to be repeated multiple times, and objects combined (see
# script 4.
