##====== GRADIENT FOREST MODELS IN NEW ZEALAND WATERS ========================##

# SCRIPT 5 - Classifying mean bootstrapped compositional turnover into groups

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

# DESCRIPTION: Classification of Gradient Forest models to different levels  
# and exploration of classification strength across biotic groups. 

# 1.  Load compositionl turnover and other files
# 2.  Heirarchical classsification of compositional turnover into 5 - 150 groups
# 3.  Export summary results for each classification (5 - 150 groups) including 
#     plotting PCA, plot and export raster, and tag species' data with group num
# 4.  Analyse classification strength for each biotic group (anomsim, pairwise 
#     differences, R2, p-values).
#       4.1 Reef fish
#       4.2 Macroalgae
#       4.3 Demersal fish  
#       4.4 Benthic invertebrates

####==========    1. LOAD FILES AND PACKAGES   =============================####
# For each taxa, provide Presence Absence dataframe; dataframe of predictors; 
# list of species and env preds
rm(list = ls())
library(cluster); require(raster); require(tidyverse); require(cluster); require(devEMF)
require(vegan); require(ecodist); require(parallel)

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

imp.vars <- c("Bathy","BedDist", "BotOxy", "BotNi", "BotPhos","BotSal",
              "BotSil", "BotTemp","BPI_broad", "BPI_fine", "ChlAGrad",
              "DET", "PB555nm","SeasTDiff", "Slope", "SSTGrad", "sed.class", 
              "TC", "POCFlux","Ebed")

# template raster in case 
setwd(paste0(dir, "/Test data"))
r <- raster::raster("Template_1km.tif") ### Getting FID from raster to sum sightings
plot(r)

# load 1km predictors for whole area; extract for imp.var
MPIproj<- CRS("+proj=aea +lat_1=-30 +lat_2=-50 +lat_0=-40 +lon_0=175 +x_0=0 
              +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# load master predictor stack 
setwd(paste0(dir, "/Test data"))
load("Pred_1km.CMB.Source") # (1 km resolution)

# Load biological data for taxa of interest
load("DF.source") # DF bio 
load("RF.source") # RF bio
load("MA.source") # MA bio
load("BI.source") # BI bio

# generate a depth mask if needed
Pred_1km.CMB$Bathy[Pred_1km.CMB$Bathy < 0] <- 0

# bring in the turnover object for multi-level classification 
setwd(paste0(dir, "/bootstrap_results"))
load(file = "Pred_EEZ_CMB.source")

####==========    2. HEIRARCHICAL CLASSIFICATION    ========================####
# CLASSIFY 
StartTime <- Sys.time()
EEZClaraClassification <- clara(x = Pred_EEZ_CMB[,imp.vars],    
                                   k = 500,
                                   samples = 50,
                                   metric = 'manhattan',   #manhattan distance appropriate because transformation takes care of scaling differences
                                   trace = 2,
                                   pamLike = T)
# FOR FINAL RUN CHANGE SAMPLES = 50
EndTime <- Sys.time()
print(EndTime - StartTime)

# Classification with 4 mil pred points for 20 env vars with 10 samples: 16 mins

setwd(paste0(dir, "/Classification"))
save('EEZClaraClassification',file = 'EEZClaraClassification.source')
load("EEZClaraClassification.source")

# loop through Transformed environmental space and 
# extract each class at 5 - 150 classes

# now create a 500 row dataset summarising the transformed envrionmental attributes for the groups from clara
EEZMedoidMeans <- matrix(0, nrow = 500, ncol = length(imp.vars))

dimnames(EEZMedoidMeans)[[1]] <- paste('Grp_',c(1:500),sep='')
dimnames(EEZMedoidMeans)[[2]] <- imp.vars

for (i in c(1:length(imp.vars))) EEZMedoidMeans[,i] <- tapply(Pred_EEZ_CMB[,imp.vars[[i]]],
                                                              EEZClaraClassification[[4]],mean)

summary(EEZMedoidMeans)

# and apply a hierarchical classification to it using agnes
EEZMedoidAgnesClassification <- agnes(EEZMedoidMeans, 
                                      metric = 'manhattan',
                                      method = 'gaverage',
                                      par.method = -0.1)

# plot(EEZMedoidAgnesClassification, which.plots = 2)
# rect.hclust(EEZMedoidAgnesClassification, k=25, border="red")
# rect.hclust(EEZMedoidAgnesClassification, k=75, border="blue")

# now reduce to a smaller number of groups based on the clara results
ClaraGroupExpansion <- cutree(EEZMedoidAgnesClassification,1:500)

# head(ClaraGroupExpansion)
i <- match(EEZClaraClassification$clustering,ClaraGroupExpansion[,500])
# summary(i)

EEZClaraClassification <- EEZClaraClassification[-c(11:length(EEZClaraClassification))]

Group.Num <- seq(5,150,5)
Group.Name <- paste("Grp_",Group.Num, sep = "")

for (j in 1:length(Group.Num)){
  index <- j + 10
  EEZClaraClassification[[index]] <- ClaraGroupExpansion[i,Group.Num[j]]
  names(EEZClaraClassification)[index] <- Group.Name[j]
}

setwd(paste0(dir, "/Classification"))
save(EEZClaraClassification,file = 'EEZClaraClassification_5_150.source')

####==========    3. EXPORT CLASSIFICATIONS RESULTS    =====================####
setwd(paste0(dir, "/Classification"))
load('EEZClaraClassification_5_150.source')

Group.Num <- seq(5,150,5)
Group.Name <- paste("Grp_",Group.Num, sep = "")

for (k in 1:length(Group.Num)){
  Means <- matrix(0, nrow = Group.Num[k], ncol = length(imp.vars))
  dimnames(Means)[[1]] <- paste('Grp_', c(1:Group.Num[k]),sep='')
  dimnames(Means)[[2]] <- imp.vars
  
  for (i in c(1:length(imp.vars))) Means[,i] <- tapply(Pred_EEZ_CMB[,imp.vars[[i]]],
                                                               EEZClaraClassification[[Group.Name[k]]],mean)
  EEZClusterPCA <- prcomp(Means)
  
  # set up colours using the same PCA space
  a1 <- EEZClusterPCA$x[, 1]
  a2 <- EEZClusterPCA$x[, 2]
  a3 <- EEZClusterPCA$x[, 3]
  r <- a1 + a2
  g <- -a2
  b <- a3 + a2 - a1
  r <- (r - min(r))/(max(r) - min(r)) * 255
  g <- (g - min(g))/(max(g) - min(g)) * 255
  b <- (b - min(b))/(max(b) - min(b)) * 255
  
  # PLOT PCA ------------------------------------------------------------------
  nvs <- dim(EEZClusterPCA$rotation)[1]
  vec <- imp.vars[c(1,2,3,4,6,8,9,16,18,19)]
  lv <- length(vec)
  vind <- rownames(EEZClusterPCA$rotation) %in% vec
  
  scal <- 15
  xrng <- range(EEZClusterPCA$x[, 1], EEZClusterPCA$rotation[, 1]/scal) * + 1.1
  yrng <- range(EEZClusterPCA$x[, 2], EEZClusterPCA$rotation[, 2]/scal) * + 1.1
  
  variableCEX <- (as.numeric(table(EEZClaraClassification[[Group.Name[k]]]))^0.110) * 2
  
  # save as EMF
  dir.create(paste(dir, "/classification/", Group.Name[k],sep =""))
  setwd(paste(dir, "/classification/", Group.Name[k],sep =""))
  emf(file = paste(Group.Name[k], "_EEZ_PCA.emf"), emfPlus = FALSE)
  plot((EEZClusterPCA$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = variableCEX, 
       col = rgb(r, g, b, max = 255), asp = 1)
  arrows(rep(0, lv), rep(0, lv), EEZClusterPCA$rotation[vec,1]/scal, EEZClusterPCA$rotation[vec, 2]/scal, length = 0.06210)
  jit <- 0.00110
  text(EEZClusterPCA$rotation[vec, 1]/scal + jit * sign(EEZClusterPCA$rotation[vec, 1]), 
       EEZClusterPCA$rotation[vec, 2]/scal + jit * sign(EEZClusterPCA$rotation[vec, 2]), labels = vec)
  text(EEZClusterPCA$x[,1], EEZClusterPCA$x[,2], seq(1,Group.Num[k],1),
       cex = variableCEX/7, #0.8,
       pos = 4,
       adj = c(0.2,-0.1))
  dev.off() # finishes the plot and saves
  
  # SAVE MAP. ------------------------------------------------------------------
  rFull <- r[EEZClaraClassification[[Group.Name[k]]]] #clustering]
  gFull <- g[EEZClaraClassification[[Group.Name[k]]]]
  bFull <- b[EEZClaraClassification[[Group.Name[k]]]]
  
  jpeg(filename = paste(Group.Name[k], "_EEZ_Spatial.jpeg", sep = ""), width = 11000, height = 11000, quality = 100, bg = "white", res = 1000)
  par(mar=c(2,2,2,2))
  plot(Pred_EEZ_CMB[,1:2],
       cex = 0.7, 
       col = rgb(rFull, gFull, bFull, max = 255),
       pch = '.',
       asp = 1)
  title(paste(Group.Name[k]))
  dev.off()
  
  # RASTER ---------------------------------------------------------------------
  # create a dataframe with spatial information for th classification 
  CMB_DF <- cbind(Pred_EEZ_CMB[,c(1:2)],EEZClaraClassification[[Group.Name[k]]])
  
  # export as raster
  CMB_DF.1km.R <- rasterFromXYZ(data.frame(x = CMB_DF[,1],
                                           y = CMB_DF[,2],
                                           z = CMB_DF[,3]),
                                       crs = MPIproj)
  # plot(CMB_DF.1km.R)
  GF <- CMB_DF.1km.R
  writeRaster(GF, filename= paste(Group.Name[k], ".tif", sep = ""), 
              format = "GTiff", 
              overwrite = TRUE)
  
  # TAG BIOLOGICAL DATA --------------------------------------------------------
  # DF
  DF.class <- as.data.frame(raster::extract(GF,DF[,c("X","Y")]))
  colnames(DF.class) <- Group.Name[k]
  DF <- cbind(DF, DF.class)
  # RF
  RF.class <- as.data.frame(raster::extract(GF,RF[,c("X","Y")]))
  colnames(RF.class) <- Group.Name[k]
  RF <- cbind(RF, RF.class)
  # MA
  MA.class <- as.data.frame(raster::extract(GF,MA[,c("X","Y")]))
  colnames(MA.class) <- Group.Name[k]
  MA <- cbind(MA, MA.class)
  # DF
  BI.class <- as.data.frame(raster::extract(GF,BI[,c("X","Y")]))
  colnames(BI.class) <- Group.Name[k]
  BI <- cbind(BI, BI.class)
}

setwd(paste0(dir, "/Test data"))
save(DF, file = "DF.class.source")
save(RF, file = "RF.class.source")
save(MA, file = "MA.class.source")
save(BI, file = "BI.class.source")

####==========    4. ANALYSIS OF SPECIES DATA     ==========================####
setwd(dir)
source("./Pairwise_adonis.R")

# biological groups with class number tags
setwd(paste0(dir, "/Test data"))
load("DF.class.source") # DF bio 
load("RF.class.source") # RF bio
load("MA.class.source") # MA bio
load("BI.class.source") # BI bio

####    4.1 REEF FISH   ----------------------------------------------------------
RF <- na.omit(RF)
Summ.table <- data.frame(Num.Classes =integer(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                         R2 = double(), p.val = double(), stringsAsFactors=FALSE)
# j = 1
for (j in 1:length(Group.Num)){
  grps.sum <- as.data.frame(table(RF[,Group.Name[j]]))
  grps.sum <- grps.sum[grps.sum[,2] > 4,1]
  
  RF.cut <- RF[RF[,Group.Name[j]] %in% grps.sum,]
  
  RF.spe <- RF.cut[,c(24:(ncol(RF.cut)-32))]
  for (i in 1:ncol(RF.spe)){RF.spe[,i] <- as.numeric(RF.spe[,i])}
  RF.spe[RF.spe==1] <- 0
  RF.spe[RF.spe==2] <- 1
  RF.spe <- RF.spe[,colSums(RF.spe[,1:length(RF.spe)]) > 0]
  RF.spe <- RF.spe[rowSums(RF.spe[1:nrow(RF.spe),]) > 0,]
  
  RF.class <- RF.cut[,c((ncol(RF.cut)-31):ncol(RF.cut))]
  RF.class <- RF.class[rownames(RF.class) %in% rownames(RF.spe), ]
  
  RF.dist <- vegdist(RF.spe, distance="jaccard", binary = T)
  
  # ANOSIM
  perm <- anosim(RF.dist, RF.class[,Group.Name[j]], permutations = 100)
  # PERMANOVA
  # perm <- adonis(RF.dist ~ RF.class[,Group.Name[j]], data = RF.class, permutations = 100)
  
  Summ.table[j,1] <- Group.Num[j]
  Summ.table[j,2] <- length(grps.sum)/Group.Num[j]
  # if (any(j == c(1,5,10,15,25))){
  #   Perm.pair <- pairwise.adonis(RF.dist,
  #                                RF.class[,Group.Name[j]], perm = 25, parallel = 5)
  #   gc()
  #   Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
  # } else {
  #   Summ.table[j,3] <- 0
  # }
  # Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
  
  Summ.table[j,4] <- perm$statistic
  # Summ.table[j,4] <- perm$aov.tab$R2[1]
  Summ.table[j,5] <- perm$signif
  # Summ.table[j,5] <- perm$aov.tab$`Pr(>F)`[1]
  
  print(paste("Iteration finished for ", Group.Num[j], sep =""))
}

setwd(paste0(dir, "/Classification_assessment"))
write.csv(Summ.table, file = "RF_scores.csv")

# dispersion if you need to look at this
disp.grp <- betadisper(RF.dist, RF.class$Grp_5, type = "median")
pairwise <- permutest(disp.grp, pairwise=TRUE, permutations=100)
pairwise$pairwise[1]
plot(disp.grp, ellipse = F, hull = T, conf = 0.40) # 90% data ellips
boxplot(disp.grp)

# plot figure
plot(Summ.table[,1], Summ.table[,4])

####    4.2 MACROALGAE    ------------------------------------------------------
MA <- na.omit(MA)
Summ.table <- data.frame(Num.Classes =integer(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                         R2 = double(), p.val = double(), stringsAsFactors=FALSE)

for (j in 1:length(Group.Num)){
  grps.sum <- as.data.frame(table(MA[,Group.Name[j]]))
  grps.sum <- grps.sum[grps.sum[,2] > 4,1]
  
  MA.cut <- MA[MA[,Group.Name[j]] %in% grps.sum,]
  
  MA.spe <- MA.cut[,c(24:(ncol(MA.cut)-32))]
  for (i in 1:ncol(MA.spe)){MA.spe[,i] <- as.numeric(MA.spe[,i])}
  MA.spe[MA.spe==1] <- 0
  MA.spe[MA.spe==2] <- 1
  MA.spe <- MA.spe[,colSums(MA.spe[,1:length(MA.spe)]) > 0]
  MA.spe <- MA.spe[rowSums(MA.spe[1:nrow(MA.spe),]) > 0,]
  
  MA.class <- MA.cut[,c((ncol(MA.cut)-31):ncol(MA.cut))]
  MA.class <- MA.class[rownames(MA.class) %in% rownames(MA.spe), ]
  
  MA.dist <- bcdist(MA.spe)
  perm <- anosim(MA.dist, MA.class[,Group.Name[j]], permutations = 100, parallel = 5)
  # perm <- adonis(MA.dist ~ MA.class[,Group.Name[j]], data = MA.class, permutations = 100, parallel = 10)
 
  
  Summ.table[j,1] <- Group.Num[j]
  Summ.table[j,2] <- length(grps.sum)/Group.Num[j]
  # if (any(j == c(1,5,10,15,25))){
  #   Perm.pair <- pairwise.adonis(MA.dist,
  #                                MA.class[,Group.Name[j]], perm = 25, parallel = 5)
  #   gc()
  #   Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
  # } else {
  #   Summ.table[j,3] <- 0
  # }
  # Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
  # Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)

  Summ.table[j,4] <- perm$statistic
  # Summ.table[j,4] <- perm$aov.tab$R2[1]
  Summ.table[j,5] <- perm$signif
  # Summ.table[j,5] <- perm$aov.tab$`Pr(>F)`[1]
  
  print(paste("Iteration finished for ", Group.Num[j], sep =""))
}
ssetwd(paste0(dir, "/Classification_assessment"))
write.csv(Summ.table, file = "MA_scores.csv")
# plot(Summ.table[,1],Summ.table[,4])

####    4.3 DEMERSAL FISH   -----------------------------------------------------
DF <- na.omit(DF)
Summ.table <- data.frame(Num.Classes =integer(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                         R2 = double(), p.val = double(), stringsAsFactors=FALSE)
start <- Sys.time() # recording the time

for (j in 1:length(Group.Num)){
  grps.sum <- as.data.frame(table(DF[,Group.Name[j]]))
  grps.sum <- grps.sum[grps.sum[,2] > 4,1]
  
  DF.cut <- DF[DF[,Group.Name[j]] %in% grps.sum,]
  
  ##IF SUBSAMPLING IS NEEDED
  # for (k in 1:length(grps.sum)){
  #   grp <- grps.sum[k]
  #   count <- as.data.frame(table(DF[,Group.Name[j]]))
  #   samp <- DF[DF[,Group.Name[j]] == grp,]
  #   if(nrow(samp) > 100){
  #     sub.samp <- samp[sample(nrow(samp), 100, replace = F), ]
  #   } else {
  #     sub.samp <-samp
  #   }
  #   if (k == 1){DF.samp <-sub.samp} else {DF.samp <- rbind(DF.samp, sub.samp)}
  # }
  # 
  # DF.cut <- DF.samp
  # 
  # DF.spe <- DF.cut[,c(24:(ncol(DF.cut)-32))]
  # for (i in 1:ncol(DF.spe)){DF.spe[,i] <- as.numeric(DF.spe[,i])}
  # DF.spe[DF.spe==1] <- 0
  # DF.spe[DF.spe==2] <- 1
  # DF.spe <- DF.spe[,colSums(DF.spe[,1:length(DF.spe)]) > 0]
  # DF.spe <- DF.spe[rowSums(DF.spe[1:nrow(DF.spe),]) > 0,]
  # 
  # DF.class <- DF.cut[,c((ncol(DF.cut)-31):ncol(DF.cut))]
  # DF.class <- DF.class[rownames(DF.class) %in% rownames(DF.spe),]
  
  DF.spe <- DF.cut[,c(24:(ncol(DF.cut)-32))]
  for (i in 1:ncol(DF.spe)){DF.spe[,i] <- as.numeric(DF.spe[,i])}
  DF.spe[DF.spe==1] <- 0
  DF.spe[DF.spe==2] <- 1
  DF.spe <- DF.spe[,colSums(DF.spe[,1:length(DF.spe)]) > 0]
  DF.spe <- DF.spe[rowSums(DF.spe[1:nrow(DF.spe),]) > 0,]
  
  DF.class <- DF.cut[,c((ncol(DF.cut)-31):ncol(DF.cut))]
  DF.class <- DF.class[rownames(DF.class) %in% rownames(DF.spe), ]
  
  DF.dist <- bcdist(DF.spe, rmzero = T)
  gc()
  perm <- anosim(DF.dist, DF.class[,Group.Name[j]], permutations = 100, parallel = 5)
  # perm <- adonis(DF.dist ~ DF.class[,Group.Name[j]], permutations = 100, parallel = 5)
  gc()
  
  Summ.table[j,1] <- Group.Num[j]
  Summ.table[j,2] <- length(grps.sum)/Group.Num[j]
  
  # if (any(j == c(1,5,10,15,25))){
  #   Perm.pair <- pairwise.adonis(DF.dist,
  #                                DF.class[,Group.Name[j]], perm = 25, parallel = 5)
  #   gc()
  #   Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
  # } else {
  #   Summ.table[j,3] <- 0
  # }

  # Summ.table[j,4] <- perm$aov.tab$R2[1]
  Summ.table[j,4] <- perm$statistic
  Summ.table[j,5] <- perm$signif
  # Summ.table[j,5] <- perm$aov.tab$`Pr(>F)`[1]
  
  print(paste("Iteration finished for ", Group.Num[j], sep =""))
}
end <- Sys.time()
end - start

# # IF SUBSAMPLING IS NEEDED
# # 
# start <- Sys.time()
# for (j in c(10,15,20,30)){
#   grps.sum <- as.data.frame(table(DF[,Group.Name[j]]))
#   grps.sum <- grps.sum[grps.sum[,2] > 4,1]
# 
#   DF <- DF[DF[,Group.Name[j]] %in% grps.sum,]
# 
#   for (k in 1:length(grps.sum)){
#     grp <- grps.sum[k]
#     count <- as.data.frame(table(DF[,Group.Name[j]]))
#     samp <- DF[DF[,Group.Name[j]] == grp,]
#     if(nrow(samp) > 100){
#       sub.samp <- samp[sample(nrow(samp), 100, replace = F), ]
#     } else {
#       sub.samp <-samp
#     }
#     if (k == 1){DF.samp <-sub.samp} else {DF.samp <- rbind(DF.samp, sub.samp)}
#     }
# 
#   DF.cut <- DF.samp
# 
#   DF.spe <- DF.cut[,c(24:(ncol(DF.cut)-32))]
#   for (i in 1:ncol(DF.spe)){DF.spe[,i] <- as.numeric(DF.spe[,i])}
#   DF.spe[DF.spe==1] <- 0
#   DF.spe[DF.spe==2] <- 1
#   DF.spe <- DF.spe[,colSums(DF.spe[,1:length(DF.spe)]) > 0]
#   DF.spe <- DF.spe[rowSums(DF.spe[1:nrow(DF.spe),]) > 0,]
# 
#   DF.class <- DF.cut[,c((ncol(DF.cut)-31):ncol(DF.cut))]
#   DF.class <- DF.class[rownames(DF.class) %in% rownames(DF.spe), ]
# 
#   DF.dist <- bcdist(DF.spe, rmzero = T)
#   gc()
#   Perm.pair <- pairwise.adonis(DF.dist,DF.class[,Group.Name[j]], perm = 100)
#   gc()
#   Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
# 
#   print(paste("Iteration finished for ", Group.Num[j], sep =""))
# }
# 
# end <- Sys.time()
# end - start

setwd(paste0(dir, "/Classification_assessment"))
write.csv(Summ.table, file = "DF_scores.csv")

####    4.4 BENTHIC INVERTEBRATES   --------------------------------------------
BI <- na.omit(BI)
Summ.table <- data.frame(Num.Classes =integer(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                         R2 = double(), p.val = double(), stringsAsFactors=FALSE)
start <- Sys.time() # recording the time

for (j in 1:length(Group.Num)){
  grps.sum <- as.data.frame(table(BI[,Group.Name[j]]))
  grps.sum <- grps.sum[grps.sum[,2] > 4,1]
  
  BI.cut <- BI[BI[,Group.Name[j]] %in% grps.sum,]
  
  BI.spe <- BI.cut[,c(3:(ncol(BI.cut)-32))]
  for (i in 1:ncol(BI.spe)){BI.spe[,i] <- as.numeric(BI.spe[,i])}
  # BI.spe[BI.spe==1] <- 0
  # BI.spe[BI.spe==2] <- 1
  BI.spe <- BI.spe[,colSums(BI.spe[,1:length(BI.spe)]) > 0]
  BI.spe <- BI.spe[rowSums(BI.spe[1:nrow(BI.spe),]) > 0,]
  
  BI.class <- BI.cut[,c((ncol(BI.cut)-29):ncol(BI.cut))]
  BI.class <- BI.class[rownames(BI.class) %in% rownames(BI.spe), ]
  
  BI.dist <- bcdist(BI.spe, rmzero = T)
  perm <- anosim(BI.dist, BI.class[,Group.Name[j]], permutations = 100, parallel = 5)
  # perm <- adonis(BI.dist ~ BI.class[,Group.Name[j]], permutations = 100, parallel = 5)
  
  Summ.table[j,1] <- Group.Num[j]
  Summ.table[j,2] <- length(grps.sum)/Group.Num[j]
  Summ.table[j,4] <- perm$statistic
  Summ.table[j,5] <- perm$signif
  
  # Summ.table[j,4] <- perm$aov.tab$R2[1]
  # Summ.table[j,5] <- perm$aov.tab$`Pr(>F)`[1]
  # 
  # if (any(j == c(1,2,3,4,5,10,15,20,30))){
  #   for (k in 1:length(grps.sum)){
  #     grp <- grps.sum[k]
  #     count <- as.data.frame(table(BI.cut[,Group.Name[j]]))
  #     samp <- BI.cut[BI.cut[,Group.Name[j]] == grp,]
  #     if(nrow(samp) > 100){
  #       sub.samp <- samp[sample(nrow(samp), 100, replace = F), ]
  #     } else {
  #       sub.samp <-samp
  #     }
  #     if (k == 1){BI.samp <-sub.samp} else {BI.samp <- rbind(BI.samp, sub.samp)}
  #   }
    
    # BI.cut <- BI.samp
    
  #   BI.spe <- BI.cut[,c(3:(ncol(BI.cut)-32))]
  #   for (i in 1:ncol(BI.spe)){BI.spe[,i] <- as.numeric(BI.spe[,i])}
  #   # BI.spe[BI.spe==1] <- 0
  #   # BI.spe[BI.spe==2] <- 1
  #   BI.spe <- BI.spe[,colSums(BI.spe[,1:length(BI.spe)]) > 0]
  #   BI.spe <- BI.spe[rowSums(BI.spe[1:nrow(BI.spe),]) > 0,]
  #   
  #   BI.class <- BI.cut[,c((ncol(BI.cut)-29):ncol(BI.cut))]
  #   BI.class <- BI.class[rownames(BI.class) %in% rownames(BI.spe), ]
  #   
  #   BI.dist <- bcdist(BI.spe, rmzero = T)
  #   Perm.pair <- pairwise.adonis(BI.dist,BI.class[,Group.Name[j]], perm = 100)
  #   Summ.table[j,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
  #   } else {
  #   Summ.table[j,3] <- 0
  # }
  
  print(paste("Iteration finished for ", Group.Num[j], sep =""))
}
end <- Sys.time()
end - start

setwd(paste0(dir, "/Classification_assessment"))
write.csv(Summ.table, file = "BI_scores.csv")
