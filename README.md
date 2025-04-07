README
The following readme describes the structure of the R code used to develop the New Zealand Seafloor Community Classification (NZSCC) and associated biodiversity layers described in: Stephenson, F., Rowden, A.A., Brough, T., Petersen, G., Bulmer, R.H., Leathwick, J.R., Lohrer, A.M., Ellis, J.I., Bowden, D.A., Geange, S.W., Funnell, G.A., Freeman, D.J., Tunley, K., Tellier, P., Clark, D.E., Lundquist, C.J., Greenfield, B.L., Tuck, I.D., Mouton, T.L., Neill, K.F., Mackay, K.A., Pinkerton, M.H., Anderson, O.F., Gorman, R.M., Mills, S., Watson, S., Nelson, W.A., and Hewitt, J.E. (2022). Development of a Seafloor Community Classification for the New Zealand Region Using a Gradient Forest Approach. Frontiers in Marine Science 8.

R CODE:
1. Data Preparation.R
Description: Preparation of environmental and biological data for use in  Gradient Forest modelling.
2. Gradient Forest Model Tuning.R
Description: example of Gradient Forest model parameterisation and outputs exemplified using the Demersal fish data. Key decisions included which  environmental variables to use and appropriate sample number.
3. Bootstrapping GF Models.R
Description: bootstrapping and aggregation of Gradient Forest models for each biotic group. Key decisions were number of bootstraps and number of trees to allow models to run in parallel without running out memory.
4. Combining results of GF models.R
Description: summarising the bootstrapped Gradient Forest models for each biotic group. Given the size of the objects, bootstrapped models needed to be  saved in separate lists
5. Classification.R
Description: Classification of Gradient Forest models to different levels (5 Ð 150 in increments of 5) and exploration of classification strength across biotic groups.
6. Model uncertainty.R
Description: Assessment of spatial uncertainty of Gradient Forest models: uncertainty in compositional turnover (Standard deviation of the mean compositional turnover) and coverage of samples of the environmental space.
7. Summarise taxa and environmental values within groups.R
DESCRIPTION: Summarise taxa (occurrence, richness, etc) and environmental values (min, mediam, max) within groups.

DATA:
DF.source - data for demersal fish (for years 2070 - 2022)
Pred_1km.CMB.source Ð environmental data used at a 1km grid resolution saved as a dataframe
