##====== GRADIENT FOREST MODELS IN NEW ZEALAND WATERS ========================##

# SCRIPT 1. EXAMPLE FILE PREPARATION

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

# DESCRIPTION: Preparation of environmental and biological data for use in 
# Gradient Forest modelling. 

# 1. Species data 
# 2. Environmental Variables
# 3. Linking species and environmental data (final datasets)


####==========    LOAD  PACKAGES    ========================================####
library("raster");library("tidyverse");library("ggplot2");library("rstudioapi")

####    1. SPECIES DATA     ================================================####
# import pre prepared taxa data with projected (AEA -41 <MPI proj) X Y coords 
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(dir, "/Test data"))
DF_allyears <- read.csv("TRAWL_demersal_fish_1979_2018_proj_SP.csv", header = T)
sp.mat <-reshape2::dcast(DF_allyears, X + Y ~ SpeciesName, length)

# AGGREGATE RECORDS TO 1KM GRID CELLS (p/a)
r <- raster::raster("EEZ_NZ_coast_ras1km.tif") 
# raster::plot(r)
xy <- as.matrix(cbind(sp.mat$X, sp.mat$Y))
sp.mat$FID <- raster::cellFromXY(r, xy) ### FID from raster to aggregate records
rm(xy) # tidy workspace
head(sp.mat[,450:462])

sp.mat.sum<-sp.mat[,3:ncol(sp.mat)] %>%   # sum species cols by unique FID
     group_by(FID) %>%    
     summarise_all(funs(sum), na.rm=T)

# FID and X Y values
xy <- cbind(X = sp.mat$X, Y = sp.mat$Y)
FID <- as.data.frame(cbind(xy, FID = cellFromXY(r, xy)))
rm(xy)

# link X Y variables by FID
sp.mat.distinct <-  plyr::join(sp.mat.sum, FID, by = "FID", match = "first") 
sp.mat.distinct <- sp.mat.distinct[,-1] # remove FID and 
rm(FID); rm(sp.mat.sum)
# reshuffle X Y coords to be at start
sp.mat.distinct <- sp.mat.distinct[,c(460,461,1:459)] 

cutoff <- 9 # remove samples with < 10 unique occurrences 
xy <- sp.mat.distinct[,(1:2)]
sp.mat.distinct <- sp.mat.distinct[,-c(1:2)]
sp.mat.distinct[sp.mat.distinct > 1] <- 1 # turn into presence / absence 
# max(sp.mat.distinct$Zeus) # check it has worked
sp.mat.cutoff <- sp.mat.distinct[,colSums(sp.mat.distinct[,1:length(sp.mat.distinct)]) > cutoff]
sp.mat.FINAL <- cbind(xy,sp.mat.cutoff)

# clip to EEZ using a mask
sp.mat.FINAL <- cbind(raster::extract(r, sp.mat.FINAL[,c("X","Y")]),sp.mat.FINAL)
sp.mat.FINAL <- sp.mat.FINAL[!is.na(sp.mat.FINAL[,1]),]
sp.mat.FINAL <- sp.mat.FINAL[,-1]

# SAVE SP MATRIX
write.csv(sp.mat.FINAL, "DF.matrix_all_years_unique_SP.csv")
rm(list = ls()) # clear workspace and load only final DF
sp.mat.FINAL.EEZ <- read.csv("DF.matrix_all_years_unique_SP.csv", header =  T)
sp.mat.FINAL.EEZ <- sp.mat.FINAL.EEZ[,-1]

####    2. ENVIRONMENTAL VARIABLES           ================================####
#### FIRST AT 1KM RESOLUTION
# LOAD RASTERS AS TIFF FILES
f1km <- list.files(getwd())
ras1km <- lapply(f1km,raster) # load as raster
plot(ras1km[[25]]) # plot raster to check
predStack1km <- stack(ras1km) # convert list of rasters to stack
rm(ras1km) # tidy up workspace

# RENAME FILES (ABR.) - consitent names across all rasters and dataframes
names(predStack1km) # complete full list later
names(predStack1km) <- c("Bathy","BedDist","BotNi","BotOxy","BotOxySat",
                         "BotPhos","BotSal","BotSil","BotTemp","BPI_broad",
                         "BPI_fine","Carbonate","ChlA","ChlAGrad","DET",
                         "DynOc","Gravel","KPAR","MLD","Mud","PAR","PB555nm",
                         "Rough","Sand","Slope","SST","SSTGrad","TC","TempRes")

# SAVE 
save(predStack1km, file = "predStack1km.source") 
load("predStack1km.source")
names(predStack1km)

# CONVERT TO DATAFRAME FOR PREDICTION
Pred_1km <- as.data.frame(predStack1km, xy = T)
Pred_1km <- na.omit(Pred_1km)
Pred_1km$Bathy <- Pred_1km$Bathy*-1
save(Pred_1km, file = "Pred_1km.CMB.source")
load("Pred_1km.CMB.source")

####    3. LINKING SPECIES & ENVIRONMENT      ==============================####
# extract variables and link to species data
DF <- cbind(sp.mat.FINAL.EEZ[,c(1:2)], 
            raster::extract(predStack1km, sp.mat.FINAL.EEZ[,c("X","Y")], df =T),sp.mat.FINAL.EEZ[,c(3:ncol(sp.mat.FINAL.EEZ))])
# plot(predStack1km[[1]])
# points(DF[,c(1,2)],type="p", pch=3, col="black", cex = 0.2)

# convert species data to factor to force GF into doing a classification tree
for (i in seq(33,ncol(DF))) DF[,i] <- as.factor(DF[,i])

save(DF, file = "DF.Source")
load("DF.Source")
