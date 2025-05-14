

# Reading the necessary libraries
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
setwd("C:/Users/ibk1/NoiseModelling/Glasgow/Analysis/summer")

dataGWR <- sf::read_sf("gg_gwr_input_50m.shp") # Reading as an sf 

colnames(dataGWR)

head(dataGWR)

# Exploring the data 
head(dataGWR$green_coun) #

#--------Correlation Matrix for the Identified Variables-------------------

# First I need to convert dataGWR into a dataframe, not a spatial dataframe 

dataGWR_df <- as.data.frame(dataGWR) # Converting this to a dataframe to avoid issues with the geometry column 
names(dataGWR_df)

# checking the data types 
str(dataGWR_df[, c("noise_mean", "ndvi_mean", "tree_count", "urban_coun",  "green_coun", "forest_cou")])
# showing if there are nas
summary(dataGWR_df$noise_mean)

# Now I can run the correlation test with all the potential variables 
correlation_matrix <- cor(dataGWR_df[, c("noise_mean", "ndvi_mean", "tree_count", "urban_coun",  "green_coun", "forest_cou")], use = "complete.obs")

print(correlation_matrix)

library(corrplot)
corrplot(correlation_matrix, method = "circle", type = "upper")

#--------Variables selected from the correlation analysis-------------
#"noise_mean", "ndvi_mean", "tree_count",  "green_coun", "forest_cou"

#----------Mapping the dependent variable and prepping for Geoda----------
colnames(dataGWR)

tm_shape(dataGWR) + tm_fill("noise_mean", palette = "Reds", n = 7, style = "equal")

dataGWR$POLY_ID <- seq_len(nrow(dataGWR))# This creates a sequence that starts at 1 and ends at a specified value -I set this specified value to the number of rows in dataGWR

# Did this work? 
colnames(dataGWR)
# Great, now we can save to bring into GeoDa
#st_write(dataGWR, "Lden_ScotGov_SIMD_Combined_PolyID.shp", append=FALSE)# Now I am ready to calculate weights and conduct the global spatial autocorrelation analysis


#------------------Inspecting the data--------------------------
# Check how many cells are not urban, 
unclassified_cells <- dataGWRspatial[
  dataGWRspatial$green_coun == 0 &
    dataGWRspatial$forest_cou == 0 &
    dataGWRspatial$urban_coun == 0, ]

# Count and calculate percentage
n_unclassified <- nrow(unclassified_cells)
n_total <- nrow(dataGWRspatial)

percentage_unclassified <- (n_unclassified / n_total) * 100

cat("Unclassified cells:", n_unclassified, "\n")
cat("Total cells:", n_total, "\n")
cat("Percentage unclassified:", round(percentage_unclassified, 2), "%\n")

#the total of unclassified cells is 4.84. Thus they will be classed as if they are urban, a limitation of the study. 

#------------AICc Optimisation-----------------
# First, as with any GWmodel call, we need to convert data into a spatial data frame
dataGWRspatial <- as_Spatial(dataGWR)
names(dataGWRspatial)

dataGWRspatial$veg_coun <- pmax(dataGWRspatial$green_coun, dataGWRspatial$forest_cou)

DeVar <- "noise_mean"

# Defining independent variables 
InDeVar <- c("ndvi_mean" , "tree_count" ,"veg_coun")

# Run GWR with the specified independent variables
optimalBW <- bw.gwr(noise_mean ~  ndvi_mean + tree_count +veg_coun, 
                    data = dataGWRspatial, 
                    approach = "AICc", 
                    kernel = "bisquare", 
                    adaptive = TRUE)

# Running the Model Selection 
modelSel <- model.selection.gwr(DeVar, InDeVar, data=dataGWRspatial, kernel="bisquare", adaptive=TRUE, bw=565)

# Extract list of models from the results, this creates a list of the order in which they were generated
sortedModels <- model.sort.gwr(modelSel, numVars <- length(InDeVar), ruler.vector = modelSel[[2]][,2])
modelList <- sortedModels[[1]]

model.view.gwr(DeVar, InDeVar, model.list=modelList) # Viewing the radial model 

# Export in a figure
png(filename="ModelSelection_RadialView_Final.png", width=800, height = 800)
model.view.gwr(DeVar, InDeVar, model.list=modelList)
dev.off()

#----------1.3 Determining the impact of variables on AICc--------------------------------------- 
#Prepping the graph that shows the variables with the lowest impact on AICc

# number of independent variables
n <- length(InDeVar)

# Export list of AICc values from the sorted models
AICcList <- sortedModels[[2]][,2]

indices <- rep(n, n) # initialise a list of indicies as n values of n, the first position is already correct

# Based on the lab: for each position we will take the number from previous step (i-1) and add (n-i), then correct 
# by adding another 1, because we started at i=2
for (i in 2:n) {
  indices[i]=indices[i-1]+((n-i)+1)
}
# Checking what this looks like, seems to match with radial plot: 
indices

# Now let's find AICc values for models at these indices
AICcBestModelValues <- AICcList[indices]

AICcBestModelValues
# Result: [1] 6405.122 6379.580 6354.013 6339.046 6322.967 6313.754 6306.834 6298.617 6294.297 6285.674

# To be able to plot the AICc optimisation plot, we need to find out which variable was the one that was added to each best model - we do this by selecting model descriptions from the model selection result.
BestModels <- sortedModels[[1]][indices]
BestModels

# The last model has variables listed in the order of addition, and we will need this for our AICc plot.
BestModels[n]

# With the correct order 
variablesAsAdded <- c("ndvi_mean", "veg_coun", "tree_count" )

#-----------1.31 Plotting the graph with the AICc impact---------------------------------
par(mar = c(8, 4, 4, 2)) # Increase bottom margin (first value)
# Documentation for this margin increase was found here: 
# https://www.r-bloggers.com/2010/06/setting-graph-margins-in-r-using-the-par-function-and-lots-of-cow-milk/
plot(cbind(1:3,AICcBestModelValues), col = "black", pch = 20, lty = 10, 
     main = "AICc optimisation", ylab = "AICc", type = "b", axes=FALSE)
par(las=2) # This will rotate labels on x axis for 90 degrees
axis(1, at=1:3, labels=variablesAsAdded) # this plots variable names as labels on x axis, changed to 9
axis(2, at=NULL, labels=TRUE) # this plots numbers on y axis
dev.off()

#------------1.32 Calculating the AICc differences between the variables----------------------
#Difference between two consecutive AICc: Takes the first n-1 elements (1:(n-1)) and the last n-1 elements (2:n) and subtracts the second list from the first list
AICcDifference <- AICcBestModelValues[1:(n-1)]-AICcBestModelValues[2:n]
# Check how this looks
AICcDifference

#----------------2.1 Creating the new optimised GWR model using 8 variables------------------
# Optimised bandwidth: we take the bisquare adaptive kernel and use AICc for identification of the optimal bandwith

optimalBW7 <- bw.gwr(noise_mean ~ ndvi_mean+veg_coun+tree_count, data=dataGWRspatial, approach="AICc", kernel="bisquare", adaptive=TRUE)
optimalBW7 

# Run the GWR model (basic)
gwrmodel <- gwr.basic(noise_mean ~ ndvi_mean+veg_coun+tree_count, data=dataGWRspatial, bw=optimalBW7, kernel="bisquare", adaptive=TRUE) 
# As a note, it is not necessary to run the global linear model: the GWR result below gives both the local and global gwr results. 

print(gwrmodel)

#capture.output(gwrmodel, file="GWRmodel_descriptiveResult_50_BW7.txt", append = TRUE)

#----------------2.2 Calculating Global Standardised Residuals ----------------------------------------

#----------------2.21 Calculate global residuals in three steps--------------------------------

# Step 1: calculate predicted mean noise
dataGWR$veg_coun <- dataGWRspatial$veg_coun
all.equal(dataGWR$tree_count, dataGWRspatial$tree_count)

# Here is the actual equation with the values taken from the results of the GWR 
dataGWR$predictednoise_mean <- 7.3084651 + -0.0071409 * dataGWR$tree_count  + -1.7401261 * dataGWR$veg_coun + -2.8929985 * dataGWR$ndvi_mean
# look at the global regression, and the intercept is the first value, followed by the estimate for the second variable and second value  
# Check what this did:
head(dataGWR)

# Step 2: Calculate global residuals by subtracting the predicted value from the actual value 
dataGWR$globalRes <- dataGWR$noise_mean - dataGWR$predictednoise_mean
# Check what this did:
head(dataGWR)

#Step 3: Rescale the global residuals to the 0-1 range using mean and sd 
m <- mean(dataGWR$globalRes)
sd <- sd(dataGWR$globalRes)

#Calculating standardised global residuals 
dataGWR$stGlobalRes <- (dataGWR$globalRes-m)/sd

#------------------2.3 Creating the map of GWR by joining with the data--------------------------------------
# Turning the gwrmodel into a dataframe 
results <- as.data.frame(gwrmodel$SDF)

names(results)
head(results)

#To get the stats on the local regression
summary(results)

# Appending the GWR data and the map data and then renaming variables 
mapGWR <- cbind(dataGWR, as.matrix(results))
head(mapGWR)
names(mapGWR)

# Renaming them from .1 to _beta
names(mapGWR)[which(names(mapGWR)=="Intercept")] <- "Intercept_beta"
names(mapGWR)[which(names(mapGWR)=="tree_count.1")] <- "tree_count_beta"
names(mapGWR)[which(names(mapGWR)=="veg_coun.1")] <- "veg_count_beta"
names(mapGWR)[which(names(mapGWR)=="ndvi_mean.1")] <- "ndvi_mean_beta"

names(mapGWR)
names(mapGWR) <- gsub("veg_coun", "veg_count", names(mapGWR))
names(mapGWR)[names(mapGWR) == "veg_countt_beta"] <- "veg_count_beta"

#------------------------2.31 Setting Up the Bounding Boxes----------------------------------------
#--------------getting local GWR coefs--------------------
# Extract the local coefficients
local_coefs <- gwrmodel$SDF@data

# Check column names (optional)
names(local_coefs)

# Calculate mean and sd for each coefficient
summary_stats <- sapply(local_coefs, function(x) {
  if (is.numeric(x)) c(mean = mean(x), sd = sd(x)) else NULL
})

# Transpose and convert to a data frame
summary_df <- as.data.frame(t(summary_stats))

# View result
print(summary_df)

#----------creating bounding box--------------------------
# Get the bounding boxes of the parameter est maps 
bbox1 <- st_bbox(mapGWR)
# Range of x and y values
xrange <- bbox1$xmax - bbox1$xmin # range of x values
yrange <- bbox1$ymax - bbox1$ymin # range of y values
# Extend the right dimension by 25% more space
bbox1[3] <- bbox1[3] + (0.25 * xrange)
# Extend the bottom dimension by 25% more space
bbox1[2] <- bbox1[2] - (0.25 * yrange)
# Convert bounding box to a simple feature geometry
bbox1 <- bbox1 %>% st_as_sfc()

#---------------2.4 Mapping parameter estimates and T values for each of the 7 variables----------------

variables <- list(
  ndvi_mean = "PRGn",
  tree_count = "PRGn",
  veg_count = "PRGn"
)


for (var in names(variables)) { # names(variables) means that this loop will iterate over the names of the variables 
  # that I just defined (so leaves out the colormap)
  # Dynamically create column names
  t_value_col <- paste0(var, "_TV") # This creates an empty t value column that has the name of the variable with _TV appended to it (which we already defined in the GWR model) 
  # I was able to figure out how to use paste0 from this documentation
  #https://www.digitalocean.com/community/tutorials/paste-in-r
  beta_col <- paste0(var, "_beta") # This does the same for beta 
  beta_sig_col <- paste0(var, "_beta_sig") # Now for beta sig 
  
  # Identify non-significant areas --> because we defined t value col as an empty list that will be filled in iteratively for each variable name, we only need to write this once 
  whereNonSig <- which(mapGWR[[t_value_col]] > -1.96 & mapGWR[[t_value_col]] < 1.96) # 1.96 is the cut off four significant t values 
  
  # Copy original parameter estimates
  mapGWR[[beta_sig_col]] <- mapGWR[[beta_col]]
  # --> it's the same step as here: 
  # mapGWR$Pct_LTH_prob_beta_sig <- mapGWR$Pct_LTH_prob_beta, but instead it has double brackets in order to reference something in a list 
  # the documentation for this double brackets was found here: 
  #https://www.dataquest.io/blog/for-loop-in-r/
  
  # Set non-significant areas to NA
  mapGWR[[beta_sig_col]][whereNonSig] <- NA
  # Same as this line but for the empty list of beta sig col and with double brackets as above:
  # mapGWR$Pct_LTH_prob_beta_sig[whereNonSig_Pct_LTH_prob] <- NA  # Set non-significant areas to NA 
  
} # End for 
maps <- list() 

for (var in names(variables)) { 
  beta_sig_col <- paste0(var, "_beta_sig")
  palette <- variables[[var]] 
  
  # Find symmetric max for breaks
  # Find symmetric max for breaks
  max_abs <- max(abs(mapGWR[[beta_sig_col]]), na.rm = TRUE)
  
  # Define custom breaks, e.g., 6 breaks centered at 0
  breaks <- seq(-max_abs, max_abs, length.out = 7)  # Adjust the number as needed
  
  maps[[var]] <- tm_shape(mapGWR) +
    tm_fill(
      beta_sig_col,
      style = "fixed",
      breaks = breaks,
      palette = palette,
      colorNA = "gray80",
      textNA = "Non-significant"
    ) +
    tm_layout(
      legend.outside = FALSE,
      legend.position = c(0.8, 0.2),
      legend.bg.color = "white",
      legend.bg.alpha = 0.7,
      frame = FALSE,
      outer.margins = c(0.05, 0.05, 0.05, 0.05)
    )
}

Map_tree_count_beta_final <- maps[["tree_count"]]# I've renamed every map in the maps list by assigning the individual map name to its own variable
# This helps with plotting them on a grid 
Map_ndvi_mean_beta_final <- maps[["ndvi_mean"]]
Map_veg_count_beta_final <- maps[["veg_count"]]

# Set tmap to plot mode (important!)
tmap_mode("plot")

# ✅ Plot each map one by one
print(Map_tree_count_beta_final)
print(Map_ndvi_mean_beta_final)
print(Map_veg_count_beta_final)

# Create a 2x2 grid for the first three maps
grid.newpage()
pushViewport(viewport(layout=grid.layout(2, 2)))  # 2 rows, 2 columns

# Place the maps in the grid layout
print(Map_tree_count_beta_final, vp=viewport(layout.pos.col=1, layout.pos.row=1))
print(Map_ndvi_mean_beta_final, vp=viewport(layout.pos.col=2, layout.pos.row=1))
print(Map_veg_count_beta_final, vp=viewport(layout.pos.col=1, layout.pos.row=2))

dev.off()
#-------------2.5 Mapping Local R 2------------
tm_shape(mapGWR, bbox=bbox2) + tm_fill("Local_R2", style="equal", n=7, palette="Greens")#+tm_borders()
# These values were already calculated when we ran the gwrmodel

# ________________2.6 Mapping Local Residual values with Global Residuals_________________________ 
ming <- min(mapGWR$stGlobalRes)
maxg <- max(mapGWR$stGlobalRes) # Finding the min and max of the global standardised residuals
ming
maxg

breaksGlRes <- c(ming, -2.58, -1.96, 0, 1.96, 2.58, maxg)

head(mapGWR$stGlobalRes) 

# Sort the breaks in ascending order
breaksGlRes <- sort(breaksGlRes, na.last = TRUE)  # This sorts the breaks and removes NAs at the end

globRes<-tm_shape(mapGWR, bbox=bbox1) + tm_fill("stGlobalRes", style="fixed", breaks=breaksGlRes, palette="RdBu")# +tm_borders+tm_borders()

# Check if breaks are sorted and contain no NAs
print(breaksGlRes)  # Make sure it is sorted and no NAs

print(globRes)

#Local Residual, which is already calculated for us in the GWRMODEL
minl <- min(mapGWR$Stud_residual)
maxl <- max(mapGWR$Stud_residual)

breaksLRes <- c(minl, -2.58, -1.96, 0, 1.96, 2.58, maxl)

locRes <- tm_shape(mapGWR, bbox=bbox1) + tm_fill("Stud_residual", style="fixed", breaks=breaksLRes, palette="RdBu")# +tm_borders()

grid.newpage() # Plotting on a grid to have them side by side 
pushViewport(viewport(layout=grid.layout(1,2)))
print(globRes, vp=viewport(layout.pos.col = 1, layout.pos.row =1))
print(locRes, vp=viewport(layout.pos.col = 2, layout.pos.row =1))
dev.off()

head(mapGWR)
#------Saving to then calculate local residuals in GeoDa-----------
st_write(mapGWR,"GG_GWR_50m_results_rq1.shp", append=FALSE)

