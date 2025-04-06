##====== GRADIENT FOREST MODELS IN NEW ZEALAND WATERS ========================##

# SCRIPT 2. GRADIENT FOREST MODEL TUNING

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

# DESCRIPTION: example of Gradient Forest model parameterisation and outputs
# exemplified using the Demersal fish data. Key decisions included which 
# environmental variables to use and appropriate sample number.  

# 1. Correlation of environmental predictors at sampling locations
# 2. Exploration of sample size for GF modelling
# 3. Gradient Forest modelling and paramaterisation
# 4. Exploration of the transformed environmental space (compositional turnover)
# 5. Exploration of classification of the transformed environmental space
# 6. Example of processing and outputs of classification using 30-groups

####==========    LOAD  PACKAGES & DATA  ===================================####
library("raster");library("rstudioapi");library(cluster); library(corrplot)
require(extendedForest); require(gradientForest)

set.seed(1987)

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(dir, "/Test data"))
load("Pred_1km.CMB.Source") # spatial data for prediction
load("DF.Source") # species data for tuning model

MPIproj <- CRS("+proj=aea +lat_1=-30 +lat_2=-50 +lat_0=-40 +lon_0=175 +x_0=0 
              +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")

# Function for correlation testing
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

####    1. CORRELATION IN ENV PREDICTORS   ==================================####

######### correlations in the data
head(DF[,c(3:22)]) # Environmental variables only

# Environmental variables
var.names <- names(DF)[c(3:22)]

M <- DF[,var.names]
M <- na.omit(M)
M<-cor(M)
round(M,2)

# matrix of the p-value of the correlation
p.mat <- cor.mtest(M)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corr.plt <- corrplot::corrplot(M, method="color", col=col(200),  
                     type="upper", order="hclust", 
                     addCoef.col = "black", # Add coefficient of correlation
                     tl.col="black", tl.srt=45, #Text label color and rotation
                     # Combine with significance
                     number.cex =0.5,
                     p.mat = p.mat, sig.level = 0.05, insig = "blank", 
                     # hide correlation coefficient on the principal diagonal
                     diag=FALSE)

####    2. REDUCING DATA SIZE FOR MODEL TUNING       ========================####
# reduce number or rows = 5000 samples
set.seed(1987)
train_ind <- sample(seq_len(nrow(DF)), size = 5000) # index of rows for 3000 samples
DF_train <- DF[train_ind, ]
DF_train <- na.omit(DF_train)

####    3. GRADIENT FOREST MODEL   ==========================================####
# setup the lev parameter
nSites <- nrow(DF_train)
lev <- floor(log2(nSites * 0.368/2))
# Species
spe.name <- names(DF_train)[c(23:ncol(DF_train))]

# Specify the model
# NOTE - using compact gives a model that is about 2% of the size of the 
# non-compact, with no diff in performance / almost identical with respect 
# to both species R2 and predictor contribution given the size of the dataset, 
# have specified a large number of bins all final models will be fitted 
# with 500 trees

StartTime <- Sys.time()
FishFull.gf <- gradientForest(DF_train,
                              predictor.vars = var.names,
                              response.vars = spe.name,  # can apply filter here if species occ low with small sample
                              ntree = 20, # can vary number of trees if need to speed it up / for bootstrapping
                              transform = NULL, 
                              compact = T,
                              nbin = 401, 
                              maxLevel = lev, 
                              corr.threshold = 0.5, 
                              trace = T)

EndTime <- Sys.time() # 5000 samples 100 trees for 320 species on laptop = 12.5 mins
print(EndTime - StartTime)

setwd(paste0(dir, "/Outputs"))
save(FishFull.gf, file = "Fish.5k.gf.source")
load("Fish.5k.gf.source")

# EXPLORE MODEL OUTPUTS
# use this code if errors crop up for plotting of GF plots with incdividual species curves
# plot.window.orig <- plot.window
# plot.window <- function(xlim, ylim, log="", asp=NA, ...) {
#      if (!all(is.finite(xlim))) xlim <- c(0,1)
#      if (!all(is.finite(ylim))) ylim <- c(0,1)
#      plot.window.orig(xlim, ylim, log="", asp=NA, ...)
# }
# assignInNamespace("plot.window", plot.window, "graphics")

plot(FishFull.gf, plot.type = "O")

# change final bit to reflect the number of env preds
most_important <- names(importance(FishFull.gf))[1:8] 

plot(FishFull.gf, plot.type = "S", imp.vars = most_important,
     leg.posn = "topright", cex.legend = 1, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9, 
     par.args = list(mgp = c(1.5,0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))

plot(FishFull.gf, plot.type = "Cumulative.Importance", imp.vars = most_important,
     show.species = F, show.overall=T, common.scale = F, cex.axis = 1.2 ,
     cex.lab = 1.2, line.ylab = 0.9, 
     par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,0.3, 0, 0)))

plot(FishFull.gf, plot.type = "P", show.names = T, horizontal = F,
     cex.axis = 1, cex.labels = 0.7, line = 2.5)

speciesR2 <- FishFull.gf$result
length(speciesR2);min(speciesR2); mean(speciesR2); max(speciesR2)

####    4. TRANSFORMED ENVIRONMENTAL SPACE    ==============================####
# check that variable names match and edit if required
match(names(DF_train)[var.names],names(Pred_1km.CMB)[var.names])
# if mismatches correct here, otherwise ignore

# get the most important predictors
imp.vars <- importance(FishFull.gf)
imp.vars <- sort(names(imp.vars[imp.vars>0]))
# only use variables included in the model with positive R2
Pred_1km <- cbind(Pred_1km.CMB[,c(1:2)], Pred_1km.CMB[, imp.vars]) 

# and create the predictions
# extrapolate = F is important assumption - check package help for further information
FullEEZTransPreds <- cbind(Pred_1km[ , c("x", "y")],
                               + predict(FishFull.gf, Pred_1km[, imp.vars], extrap=F))

# Fit a principal components model 
PCs <- prcomp(na.omit(FullEEZTransPreds[, imp.vars]))
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
vec <- names(importance(FishFull.gf))[1:8] 
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <- 20
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * + 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * + 1.1

#  PLOT IN PCA SPACE
quartz()
plot((PCs$x[, c(1,2)]), xlim = xrng, ylim = yrng, pch = ".", cex = 2, col = rgb(r, g, b, max = 255), asp = 1)
# points(PCs$rotation[!vind, c(1,2)]/scal, pch = "+")
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec,1]/scal, PCs$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec, 1]), 
     PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec, 2]), labels = vec)
title('Fish GF model')

# PLOT IN GEOGRAPHIC SPACE
quartz()
plot(FullEEZTransPreds[, c("x", "y")], pch = ".", cex = 1, asp = 1, col = rgb(r, g, b, max = 255))
title('Fish GF model')
# Just around the cook strait
plot(FullEEZTransPreds[, c("x", "y")], pch = ".", cex = 3, asp = 1, col = rgb(r, g, b, max = 255),
     xlim=c(-500000, 500000), ylim=c(-500000, 500000)) # Adjust window size to focus on areas of importance
# points(FullEEZTransPreds[, c("x", "y")], pch = ".", cex = 3, asp = 1,
#        xlim=c(5000000, 7500000), ylim=c(-5500000, -2300000))

# SAVE TURNOVER AS RASTER
#expoting as RGB bands for plotting in arc gis
GF_ras_R <- rasterFromXYZ(data.frame(x = FullEEZTransPreds[,1], 
                                     y = FullEEZTransPreds[,2], 
                                     z = r),
                          crs = MPIproj)
plot(GF_ras_R)
GF_ras_G <- rasterFromXYZ(data.frame(x = FullEEZTransPreds[,1], 
                                     y = FullEEZTransPreds[,2], 
                                     z = g),
                          crs = MPIproj)
plot(GF_ras_G)
GF_ras_B <- rasterFromXYZ(data.frame(x = FullEEZTransPreds[,1], 
                                     y = FullEEZTransPreds[,2], 
                                     z = b),
                          crs = MPIproj)
plot(GF_ras_B)

# write the geotiff - combining multiband rasters in r
writeRaster(GF_ras_R,"GF_R.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_G,"GF_G.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_B,"GF_B.tif","GTiff", overwrite=TRUE)

####    5. CLASSIFICATION OF TRANSFORMED ENV SPACE    ======================####
# CLASSIFY 
StartTime <- Sys.time()
DF.EEZClaraClassification <- clara(x = FullEEZTransPreds[,imp.vars],    
                                   k = 500,
                                   samples = 20,
                                   metric = 'manhattan',   
                                   #manhattan distance appropriate because transformation takes care of scaling differences
                                   trace = 2,
                                   pamLike = T)
# FOR FINAL RUN CHANGE SMAPLES = 50
EndTime <- Sys.time()
print(EndTime - StartTime)
# Classification with 4.1 mil pred points for 20 env vars with 20 samples: 30 mins 
# Classification with 4.1 mil pred points for 20 env vars with 50 samples:

save('DF.EEZClaraClassification',file = 'DF.EEZClaraClassification.source')
load('DF.EEZClaraClassification.source')

imp.vars <- var.names

# now create a 500 row dataset summarising the transformed envrionmental attributes for the groups from clara
DF.EEZMedoidMeans <- matrix(0, nrow = 500, ncol = length(imp.vars))

dimnames(DF.EEZMedoidMeans)[[1]] <- paste('Grp_',c(1:500),sep='')
dimnames(DF.EEZMedoidMeans)[[2]] <- imp.vars

for (i in c(1:length(imp.vars))) DF.EEZMedoidMeans[,i] <- tapply(FullEEZTransPreds[,i+2],
                                                   DF.EEZClaraClassification[[4]],mean)
# summary(DF.EEZMedoidMeans)
# and apply a hierarchical classification to it using agnes
DF.EEZMedoidAgnesClassification <- agnes(DF.EEZMedoidMeans, 
                                         metric = 'manhattan',
                                         method = 'gaverage',
                                         par.method = -0.1)
plot(DF.EEZMedoidAgnesClassification)

# now reduce to a smaller number of groups based on the clara results
ClaraGroupExpansion <- cutree(DF.EEZMedoidAgnesClassification,1:500)
i <- match(DF.EEZClaraClassification$clustering,ClaraGroupExpansion[,500])
# summary(i)

DF.EEZClaraClassification$FiveGroups       <- ClaraGroupExpansion[i,5]
DF.EEZClaraClassification$ThirtyGroups     <- ClaraGroupExpansion[i,30]
DF.EEZClaraClassification$FiftyGroups      <- ClaraGroupExpansion[i,50]
DF.EEZClaraClassification$HundredGroups    <- ClaraGroupExpansion[i,100]

setwd(paste0(dir, "/Outputs"))
save('DF.EEZClaraClassification',file = 'DF.EEZClaraClassification.source')
load('DF.EEZClaraClassification.source')

####    6. EXTRACTING AND PLOTTING 30-GROUP CLASS   ========================####
DF.EEZThirtyMeans <- matrix(0, nrow = 30, ncol = length(imp.vars))

dimnames(DF.EEZThirtyMeans)[[1]] <- paste('Grp_',c(1:30),sep='')
dimnames(DF.EEZThirtyMeans)[[2]] <- imp.vars

for (i in c(1:length(imp.vars))) DF.EEZThirtyMeans[,i] <- tapply(FullEEZTransPreds[,i+2],
                                                  DF.EEZClaraClassification$ThirtyGroups,mean) # this is the bit you could replace the preds to see the env values

head(round(DF.EEZThirtyMeans,4))
summary(DF.EEZThirtyMeans)

# dendrogram
dist <- dist(DF.EEZThirtyMeans, diag=TRUE)
hc <- hclust(dist)
# Plot the result
plot(hc)

# now calculate a PCA object for the 30-group means
DF.EEZClusterPCA <- prcomp(DF.EEZThirtyMeans)
summary(DF.EEZClusterPCA)
round(DF.EEZClusterPCA[[2]],4)

# set up colours using the same PCA space
a1 <- DF.EEZClusterPCA$x[, 1]
a2 <- DF.EEZClusterPCA$x[, 2]
a3 <- DF.EEZClusterPCA$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255

cbind(r,g,b)

## and plot groups in pca space
nvs <- dim(DF.EEZClusterPCA$rotation)[1]
vec <- names(importance(FishFull.gf))[1:6] 
lv <- length(vec)
vind <- rownames(DF.EEZClusterPCA$rotation) %in% vec

scal <- 15
xrng <- range(DF.EEZClusterPCA$x[, 1], DF.EEZClusterPCA$rotation[, 1]/scal) * + 1.1
yrng <- range(DF.EEZClusterPCA$x[, 2], DF.EEZClusterPCA$rotation[, 2]/scal) * + 1.1

variableCEX <- (as.numeric(table(DF.EEZClaraClassification$ThirtyGroups))^0.15) * 2

quartz()
plot((DF.EEZClusterPCA$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = variableCEX, 
     col = rgb(r, g, b, max = 255), asp = 1)
arrows(rep(0, lv), rep(0, lv), DF.EEZClusterPCA$rotation[vec,1]/scal, DF.EEZClusterPCA$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015
text(DF.EEZClusterPCA$rotation[vec, 1]/scal + jit * sign(DF.EEZClusterPCA$rotation[vec, 1]), 
     DF.EEZClusterPCA$rotation[vec, 2]/scal + jit * sign(DF.EEZClusterPCA$rotation[vec, 2]), labels = vec)
text(DF.EEZClusterPCA$x[,1], DF.EEZClusterPCA$x[,2], seq(1,30,1),
     cex = variableCEX/10, #0.8,
     pos = 4,
     adj = c(0.2,-0.1))
title('Thirty DF. Gradient Forest groups')

# save as EMF
emf(file = "Thirty_group.emf", emfPlus = FALSE)
plot((DF.EEZClusterPCA$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = variableCEX, 
     col = rgb(r, g, b, max = 255), asp = 1)
arrows(rep(0, lv), rep(0, lv), DF.EEZClusterPCA$rotation[vec,1]/scal, DF.EEZClusterPCA$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015
text(DF.EEZClusterPCA$rotation[vec, 1]/scal + jit * sign(DF.EEZClusterPCA$rotation[vec, 1]), 
     DF.EEZClusterPCA$rotation[vec, 2]/scal + jit * sign(DF.EEZClusterPCA$rotation[vec, 2]), labels = vec)
text(DF.EEZClusterPCA$x[,1], DF.EEZClusterPCA$x[,2], seq(1,50,1),
     cex = variableCEX/10, #0.8,
     pos = 4,
     adj = c(0.2,-0.1))
dev.off() # finishes the plot and saves

# now in geographic space
# first expand out the colours to all locations
rFull <- r[DF.EEZClaraClassification$ThirtyGroups] #clustering]
gFull <- g[DF.EEZClaraClassification$ThirtyGroups]
bFull <- b[DF.EEZClaraClassification$ThirtyGroups]

plot(FullEEZTransPreds[,1:2],
     cex = 0.5, 
     col = rgb( rFull, gFull, bFull, max = 255),
     pch = '.',
     asp = 1)

title("DF. GF model with PAM clustering - 30 groups")

# this next used to add some number labels to the map showing larger groups
identify(FullEEZTransPreds[,1],
         FullEEZTransPreds[,2],
         DF.EEZClaraClassification$ThirtyGroups,
         cex = 0.7,
         offset = 0,
         n = 50)

# save as a TIFF
DF_ThirtyGroups <- rasterFromXYZ(data.frame(x = FullEEZTransPreds[,1], 
                                                y = FullEEZTransPreds[,2], 
                                                z = DF.EEZClaraClassification$ThirtyGroups))

plot(DF_ThirtyGroups)

writeRaster(DF_ThirtyGroups, filename="DF.GF30.tif", 
            format = "GTiff", 
            overwrite = TRUE)

# this creates hexadecimal colour descriptions to allow transfer of pca based RGB colours to QGis
# values need to be manually edited into QGis legend file colours
colour_code <- cbind(seq(1,30),rgb(r, g, b, max = 255))

# for arc we can use the rgb values and save as a layer file
colour_code_arc_30 <- round(cbind(r,g,b),0)