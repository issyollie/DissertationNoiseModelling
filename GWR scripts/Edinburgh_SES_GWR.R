
#Edinburgh SES GWR 
#install.packages(sf)
library(sf)
library(tmap)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(GWmodel)
library(AICcmodavg)

#----------Importing the data previously cleaned in QGIS----------
getwd()
setwd("C:/Users/ibk1/NoiseModelling/Edinburgh/Analysis/summer")

dataGWR <- sf::read_sf("edi_GWR_Input_summer_datazones_SIMD_extra_final.shp") # Reading as an sf 

colnames(dataGWR)

# Plotting the average noise levels for each datazone with equal intervals 
head(dataGWR)

#--------Correlation Matrix for the Identified Variables-------------------

# First I need to convert dataGWR into a dataframe, not a spatial dataframe 

dataGWR_df <- as.data.frame(dataGWR) # Converting this to a dataframe to avoid issues with the geometry column 
names(dataGWR_df)

#install.packages("openxlsx")
library(openxlsx)

correlation_matrix <- cor(dataGWR_df[, c("mean_noise", "incnorm", "empnorm" , "hlthnorm", "urban_pct","depratnorm" ,"pctmineth", 
                                         "HlthCIF"  ,  "HlthAlcSR",  "HlthDrugSR" ,"HlthSMR"   , "HlthDprsPc", "HlthLBWTPc", "HlthEmergS" ,
                                         "HlthRank" ,"HouseOCrat", "HouseNCrat", "HouseRank",  "Vigintilv2", "CrimeRate", "GAccRank", "EduRank"
                                         
)])

correlation_matrix <- cor(dataGWR_df[, c("mean_noise",  "hlthnorm",
                                         "depratnorm" ,"pctmineth",
                                       "HouseOCrat", 
                                         "HouseNCrat", "HouseRank",  
                                         "Vigintilv2", "CrimeRate", "GAccRank"
                                         
)])

library(corrplot)
corrplot(correlation_matrix, method = "circle", type = "upper")


#--------Variables selected from the correlation analysis-------------
#"mean_noise", "incnorm","HouseOCrat","depratnorm",  "pctmineth", "HouseNCrat", "HlthSMR", urban_pct

#----------Mapping the dependent variable----------
#------------------------2. AICc Optimisation with Selected Variables----------------
library(GWmodel)
library(tmap)
library(sf)
library(dplyr)

# Convert to spatial object
dataGWRspatial <- as_Spatial(dataGWR)

# Define dependent and independent variables (excluding GAccRank and CrimeRate)
DeVar <- "mean_noise"
InDeVar <- c("hlthnorm", "depratnorm", "pctmineth", "HouseOCrat", "HouseNCrat")

# Find optimal bandwidth using AICc
optimalBW <- bw.gwr(mean_noise ~ hlthnorm + depratnorm + pctmineth + HouseOCrat + HouseNCrat,
                    data = dataGWRspatial, 
                    approach = "AICc", 
                    kernel = "bisquare", 
                    adaptive = TRUE)

# Model selection
modelSel <- model.selection.gwr(DeVar, InDeVar, data=dataGWRspatial, kernel="bisquare", adaptive=TRUE, bw=optimalBW)

# Sort models by AICc
sortedModels <- model.sort.gwr(modelSel, numVars=length(InDeVar), ruler.vector=modelSel[[2]][,2])
modelList <- sortedModels[[1]]

# View radial model
model.view.gwr(DeVar, InDeVar, model.list=modelList)

#------------------------3. AICc Plot Prep and Visualisation-------------------------

# Extract AICc values
n <- length(InDeVar)
AICcList <- sortedModels[[2]][,2]
indices <- rep(n, n)

for (i in 2:n) {
  indices[i] <- indices[i-1] + ((n - i) + 1)
}

AICcBestModelValues <- AICcList[indices]
BestModels <- sortedModels[[1]][indices]

# Adjust this to reflect actual variable order
variablesAsAdded <- c("depratnorm","hlthnorm", "pctmineth", "HouseNCrat", "HouseOCrat")

# Plot AICc values
par(mar = c(8, 4, 4, 2))
plot(cbind(1:5, AICcBestModelValues), col = "black", pch = 20, lty = 5,
     main = "AICc Optimisation", ylab = "AICc", type = "b", axes = FALSE)
par(las = 2)
axis(1, at = 1:length(variablesAsAdded), labels = variablesAsAdded)
axis(2, labels = TRUE)
dev.off()

#------------------------4. AICc Differences & Variable Impact------------------------

AICcDifference <- AICcBestModelValues[1:(n-1)] - AICcBestModelValues[2:n]
print(AICcDifference)

#------------------------5. Final GWR Model Using Optimal Variables-------------------

# Define new model formula without GAccRank and CrimeRate
finalFormula <- mean_noise ~ depratnorm+ hlthnorm+ pctmineth + HouseNCrat + HouseOCrat 

# Recalculate optimal bandwidth
optimalBW_final <- bw.gwr(formula = finalFormula,
                          data = dataGWRspatial,
                          approach = "AICc",
                          kernel = "bisquare",
                          adaptive = TRUE)

# Run GWR model
gwrmodel_final <- gwr.basic(formula = finalFormula,
                            data = dataGWRspatial,
                            bw = optimalBW_final,
                            kernel = "bisquare",
                            adaptive = TRUE)

# Print summary
print(gwrmodel_final)

#----------------------------global residuals---------------------
# note the chANGE 
dataGWR$predictedmean_noise <- 5.325486 +0.320662  * dataGWR$pctmineth  +  0.014732 * dataGWR$HouseOCrat  +  -5.295056   * dataGWR$depratnorm      +   0.771391   * dataGWR$hlthnorm  + 0.014732   * dataGWR$HouseNCrat        

# look at the global regression, and the intercept is the first value, followed by the estimate for the second variable and second value  
# Check what this did:
head(dataGWR)

# Step 2: Calculate global residuals by subtracting the predicted value from the actual value 
dataGWR$globalRes <- dataGWR$mean_noise - dataGWR$predictedmean_noise
# Check what this did:
head(dataGWR)

#Step 3: Rescale the global residuals to the 0-1 range using mean and sd 
m <- mean(dataGWR$globalRes)
sd <- sd(dataGWR$globalRes)

#Calculating standardised global residuals 
dataGWR$stGlobalRes <- (dataGWR$globalRes-m)/sd

#------------------ Map Prep --------------------------------------------------------
results <- as.data.frame(gwrmodel_final$SDF)
mapGWR <- cbind(dataGWR, as.matrix(results))

#--------------local coefs--------------------
# Extract the local coefficients
local_coefs <- gwrmodel_final$SDF@data

# Check column names (optional)
names(local_coefs)

# Calculate mean and sd for each coefficient
summary_stats <- sapply(local_coefs, function(x) {
  if (is.numeric(x)) c(mean = mean(x), sd = sd(x)) else NULL
})

# Transpose and convert to a data frame
summary_df <- as.data.frame(t(summary_stats))

# Optional: Keep only the coefficients (filter out other columns like diagnostics)
coef_names <- c("(Intercept)", "depratnorm", "hlthnorm", "pctmineth", "HouseNCrat", "HouseOCrat")
summary_df <- summary_df[rownames(summary_df) %in% coef_names, ]

# View result
print(summary_df)


# Bounding box for plots
bbox1 <- st_bbox(mapGWR)
xrange <- bbox1$xmax - bbox1$xmin
yrange <- bbox1$ymax - bbox1$ymin
bbox1[3] <- bbox1[3] + (0.25 * xrange)
bbox1[2] <- bbox1[2] - (0.25 * yrange)
bbox1 <- bbox1 %>% st_as_sfc()

# Rename coefficient and t-value columns
all_vars <- list(
  Intercept = "Greys",
  pctmineth = "PuOr",
  HouseNCrat = "PuOr",
  HouseOCrat = "PuOr",
  depratnorm = rev(brewer.pal(9, "Reds")),
  hlthnorm = "PuOr"
)

for (var in names(all_vars)) {
  beta_col <- paste0(var, ".1")
  tval_col <- paste0(var, ".1_TV")
  
  if (beta_col %in% names(mapGWR)) {
    names(mapGWR)[which(names(mapGWR) == beta_col)] <- paste0(var, "_beta")
  }
  if (tval_col %in% names(mapGWR)) {
    names(mapGWR)[which(names(mapGWR) == tval_col)] <- paste0(var, "_TV")
  }
}

# Identify significant coefficients
for (var in names(all_vars)) {
  tval_col <- paste0(var, "_TV")
  beta_col <- paste0(var, "_beta")
  beta_sig_col <- paste0(var, "_beta_sig")
  
  if (tval_col %in% names(mapGWR) && beta_col %in% names(mapGWR)) {
    whereNonSig <- which(mapGWR[[tval_col]] > -1.96 & mapGWR[[tval_col]] < 1.96)
    mapGWR[[beta_sig_col]] <- mapGWR[[beta_col]]
    mapGWR[[beta_sig_col]][whereNonSig] <- NA
  }
}

# Create maps
maps_all <- list()
for (var in names(all_vars)) {
  beta_sig_col <- paste0(var, "_beta_sig")
  palette <- all_vars[[var]]
  
  if (beta_sig_col %in% names(mapGWR)) {
    maps_all[[var]] <- tm_shape(mapGWR, bbox = bbox1) +
      tm_fill(beta_sig_col, style = "equal", palette = palette,
              colorNA = "lightgray", textNA = "Non-significant") +
      tm_layout(
        legend.outside = FALSE,
        legend.position = c("right", "bottom"),
        legend.bg.color = "white",
        legend.bg.alpha = 0.7,
        frame = FALSE,
        outer.margins = c(0.05, 0.05, 0.05, 0.05)
      ) +
      tm_borders()
  }
}


# Plot all maps (without GAccRank and CrimeRate)
tmap_mode("plot")
do.call(tmap_arrange, c(maps_all, ncol = 3))
tmap_mode("plot")

# Plot individual maps
maps_all[["Intercept"]]
maps_all[["pctmineth"]]
maps_all[["HouseNCrat"]]
maps_all[["HouseOCrat"]]
maps_all[["depratnorm"]]
maps_all[["hlthnorm"]]

#-------------2.5 Mapping Local R 2------------
tm_shape(mapGWR) + tm_fill("Local_R2", style="equal", n=7, palette="Greens")+tm_borders()
# These values were already calculated when we ran the gwrmodel


# ________________2.6 Mapping Local Residual values with Global Residuals_________________________ 
ming <- min(mapGWR$stGlobalRes)
maxg <- max(mapGWR$stGlobalRes) # Finding the min and max of the global standardised residuals
ming
maxg

breaksGlRes <- c(ming, -2.58, -1.96, 0, 1.96, 2.58, maxg)

head(mapGWR$stGlobalRes) 
head(mapGWR) 

# Sort the breaks in ascending order
breaksGlRes <- sort(breaksGlRes, na.last = TRUE)  # This sorts the breaks and removes NAs at the end

globRes<-tm_shape(mapGWR, bbox=bbox1) + tm_fill("stGlobalRes", style="fixed", breaks=breaksGlRes, palette="RdBu") +tm_borders+tm_borders()

# Check if breaks are sorted and contain no NAs
print(breaksGlRes)  # Make sure it is sorted and no NAs

print(globRes)

#Local Residual, which is already calculated for us in the GWRMODEL
minl <- min(mapGWR$Stud_residual)
maxl <- max(mapGWR$Stud_residual)

breaksLRes <- c(minl, -2.58, -1.96, 0, 1.96, 2.58, maxl)

locRes <- tm_shape(mapGWR, bbox=bbox1) + tm_fill("Stud_residual", style="fixed", breaks=breaksLRes, palette="RdBu") +tm_borders()

grid.newpage() # Plotting on a grid to have them side by side 
pushViewport(viewport(layout=grid.layout(1,2)))
print(globRes, vp=viewport(layout.pos.col = 1, layout.pos.row =1))
print(locRes, vp=viewport(layout.pos.col = 2, layout.pos.row =1))
dev.off()
#------Saving to then calculate local residuals in GeoDa-----------
st_write(mapGWR,"edi_GWR_Results_rq3_final_FIVE.shp", append=FALSE)



