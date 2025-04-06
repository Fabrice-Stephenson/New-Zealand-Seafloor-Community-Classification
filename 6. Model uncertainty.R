##====== GRADIENT FOREST MODELS IN NEW ZEALAND WATERS ========================##

# SCRIPT 6 - SPATIAL UNCERTAINTY OF GRADIENT FOREST PREDICTIONS

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

# DESCRIPTION: Assessment of spatial uncertainty of Gradient Forest models: 
# uncertainty in compositional turnover (Standard deviation of the mean
# compositional turnove) & coverage of samples of the environmental space. 

# 1.  Load files and packages
# 2.  Uncertainty in compositional turnover
# 3.  Coverage of the environmental space

####==========    1. LOAD FILES AND PACKAGES   =============================####
require(raster); require(dismo);require(gbm); require(devEMF)

# Load MPI projection
MPIproj<- CRS("+proj=aea +lat_1=-30 +lat_2=-50 +lat_0=-40 +lon_0=175 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")

# load master predictor stack 
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(dir, "/Test data"))
load("Pred_1km.CMB.CMB.Source") # (1 km resolution)

R <- raster("Template_1km.tif")

# loa biological data
load("DF.source") # DF bio 
load("RF.source") # RF bio
load("BI_SSG.source") # BI_SSG bio
load("BI_SMG.source") # BI_SMG bio
load("BI_MMG.source") # BI_MMG bio
load("BI_LLG.LMG.source") # BI_MMG bio
load("MA.source") # MA bio

# Combine dataframes
cmb_smp <- rbind(DF[,c(1:22)], RF[,c(1:22)], MA[,c(1:22)], 
                 BI_LLG.LMG[,c(1:22)],BI_MMG[,c(1:22)], BI_SMG[,c(1:22)], BI_SSG[,c(1:22)])
DF <- DF[,c(1:22)]
DF$sample.sites <- 1
DF$FID <- cellFromXY(R, DF[,c("X","Y")])
test <- cmb_smp$FID == unique(cmb_smp$FID)
DF <- DF[DF$FID == unique(DF$FID),]

####==========    2. Uncertainty in compositional turnover    ==============####
setwd(paste0(dir, "/bootstrap_results"))
load("SD_pred_EEZ_CMB.source") # (1 km resolution)

CMB.UC <- rasterFromXYZ(data.frame(x = Pred_1km.CMB.CMB[,1],
                                   y = Pred_1km.CMB.CMB[,2],
                                   z = Fred2.CMB_mean),
                              crs = MPIproj)
plot(CMB.UC)
setwd(paste0(dir, "/bootstrap_results/rasters"))
writeRaster(CMB.UC, filename= "CMB_SD.tif", 
            format = "GTiff", 
            overwrite = TRUE)

####==========    3. COVERAGE OF THE ENV SPACE    ==========================####

# coverage by absences 
xy <- as.matrix(cbind(Pred_1km.CMB$x, Pred_1km.CMB$y))
Pred_1km.CMB_U <-  Pred_1km.CMB
Pred_1km.CMB_U$FID <-  cellFromXY(R, xy)
Pred_1km.CMB_U <- Pred_1km.CMB[!Pred_1km.CMB_U$FID %in% DF$FID,] # no overlap with presences
set.seed(5)
train_ind <- sample(seq_len(nrow(Pred_1km.CMB_U)), size = nrow(DF)*2)
preddat_EC <- Pred_1km.CMB_U[train_ind, ]
preddat_EC <- preddat_EC[,-c(1:2)] 
# preddat_EC <- preddat_EC[,c(1:10,12,11)]
preddat_EC$sample.sites <- 0 # absent

# bind by row
DF_ES <- rbind(preddat_EC,DF[,c(3:23)])

DF_BRT_M2 <- gbm.step(data=DF_ES, gbm.x = 1:20, # environmental varibale columns
                      gbm.y = 21, # presence absence of samples
                      family = "bernoulli", tree.complexity = 2,
                      learning.rate = 0.05, bag.fraction = 0.6, n.folds=10, 
                      max.trees = 5000, plot.main = T, tolerance.method = "fixed",
                      tolerance = 0.01, verbose = T)

pred.map <- predict.gbm(DF_BRT_M2, Pred_1km.CMB, 
                        n.trees = DF_BRT_M2$gbm.call$best.trees, type = "response")
# export the map
ES.map.mean <- cbind(Pred_1km.CMB[, c("x", "y")],pred.map)

# convert to raster
BRT_ES.mean <- rasterFromXYZ(data.frame(x = ES.map.mean[,1], 
                                        y = ES.map.mean[,2], 
                                        z = ES.map.mean[,3]),
                             crs = MPIproj) # same proj as orginal env variables tiff files
plot(BRT_ES.mean)
setwd(paste0(dir, "/bootstrap_results/rasters"))
writeRaster(BRT_ES.mean,filename = "Env_cov_DF_EEZ.tif", overwrite=T) # write raster file
