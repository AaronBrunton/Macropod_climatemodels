library(raster)

# Set the directory containing ESRI Grid files
grid_dir <- "E:/EGK/Sdm_project/rdata/1990"

# Get the list of files with .flt extension
flt_files <- list.files(grid_dir, pattern = "\\.flt$", full.names = TRUE)

# Create an empty list to store raster layers
grid_rasters <- list()

# Loop through the .flt files, read each raster, and add to the list
for (flt_file in flt_files) {
  # Construct the file names based on convention
  hdr_file <- sub("\\.flt$", ".hdr", flt_file)
  prj_file <- sub("\\.flt$", ".prj", flt_file)
  
  # Read the raster using raster() function
  raster_layer <- raster(flt_file, header = hdr_file, proj4string = prj_file)
  
  # Add to the list with a suitable name
  grid_rasters[[basename(flt_file)]] <- raster_layer
}

# Stack the rasters together to create a raster stack
grid_stack <- stack(grid_rasters)

# Check the structure and summary of the raster stack
print(grid_stack)

# 5_ TNM - mean annual minimum temperature
# 8 - TXM - mean annual maximum temperature
# 9 - TXX - mean maximum monthly maximum temperature
# 7 - TXI - mean minimum monthly maximum temperature
# 4- TNI - mean minimum monthly minimum temperature
# 6 - TNX - mean maximum monthly minimum temperature
# 1- PTA - Average total annual rainfall
# 3 - PTX - mean maximum monthly rainfall
# 2 - PTI - mean minimum monthly rainfall


curr_bionames <- c("Annual_Precip",
                   "Min_Mth_Precip",
                   "Max_Mth_Precip",
                   "Min_Mth_Tmp",
                   "Annual_Min_Temp",
                   "Min_Mth_Min_Temp",
                   "Min_Mth_Max_Temp",
                   "Annual_Max_Temp",
                   "Max_Mth_Max_Temp"
)
bioclim_stack <- grid_stack
names(bioclim_stack) <- curr_bionames
names(bioclim_stack)
plot(bioclim_stack[[1]])

# Clip the raster stack to the specified extent
# Define the extent
clip_extent <- extent(151.0961, 153.7035, -28.38, -25.47475)
clipped_bioclim_stack <- crop(bioclim_stack, clip_extent)
clipped_bioclim_stack
plot(clipped_bioclim_stack[[1]], main = names(clipped_bioclim_stack[[1]]))

# Check the structure and summary of the clipped raster stack
print(clipped_bioclim_stack)
summary(clipped_bioclim_stack)
# Export the clipped raster stack
writeRaster(clipped_bioclim_stack, "clippedcsirovars_file.tif", format = "GTiff", overwrite = TRUE)

# Start importing other rasters for predictions
moist <- rast("bioclim_28.tif") # moisture
moist
# Convert SpatRaster to RasterLayer
moist_raster <- raster(moist)
rm(moist_raster_resampled)
# Resample moist_raster_matched to match the resolution of bio raster stack
moist_raster_resampled <- resample(moist_raster, bioclim_stack)
moist_raster_resampled
moist <- crop(moist_raster_resampled, clip_extent)
moist
plot(moist)

# Elevation raster
file_path <- "E:/EGK/Sdm_project/rdata/GEODATA-9-second-DEM-and-D8_-Digital-Elevation-Model-Version-3-and-Flow-Direction-Grid-2008/GEODATA-9-second-DEM-and-D8_-Digital-Elevation-Model-Version-3-and-Flow-Direction-Grid-2008"
# Load the raster
dem_raster <- raster(file.path(file_path, "dem.tif"))
dem_raster
dem <- crop(dem_raster, clip_extent)
dem <- resample(dem,clipped_bioclim_stack)
dem
rm(dem_raster)
# NDVI - Define the file path
file_path <- "E:/EGK/Sdm_project/rdata/NDVI/Australian-NDVI-(Normalised-Difference-Vegetation-Index)-October-2018---March-2019/average"
# Load the raster
ndvi_raster <- raster(file.path(file_path, "ndvi.tif"))
ndvi_raster
plot(ndvi_raster)
clipped_ndvi <- crop(ndvi_raster, clip_extent)
ndvi <- resample(clipped_ndvi, clipped_bioclim_stack)
ndvi
clipped_ndvi
rm(ndvi_raster)
plot(ndvi)
# tree cover - file path
file_path <- "E:/EGK/Sdm_project/rdata/"
# Load the raster
tcover_raster <- raster(file.path(file_path, "treecover.tif"))
clipped_tree <- crop(tcover_raster, clip_extent)
clipped_tree

# Export the clipped raster stack
writeRaster(clipped_tree, "clipped_tree.tif", format = "GTiff", overwrite = TRUE)

# Soil - Define the file path
file_path <- "E:/EGK/Sdm_project/rdata/soil/Soil-and-Landscape-Grid-National-Soil-Attribute-Maps---Soil-Colour-(3-arcsec-(~90m))---Release-2-([195---202)/datasets/environmental/soil_classification_v2/layers"
# Load the raster
soil_raster <- raster(file.path(file_path, "soil.tif"))

soil <- crop(soil_raster, clip_extent)
soil_resampled <- resample(soil,  clipped_bioclim_stack)
soil_resampled
## stack our raster layers of climate and habitat
#------ first match extent for both raster stacks

# Define the common extent (based on one of the rasters or a custom extent)
common_extent <- extent(clipped_bioclim_stack)  # Using 'moist' as the reference
common_extent
# Ensure all extents are identical
#common_extent <- extent(151.095, 153.6, -28.035, -25.475)

# Adjust the extent of clipped_tree_cropped to match the desired extent
extent(clipped_tree) <- common_extent 

# Crop each raster to the common extent
moist_cropped <- crop(moist, common_extent)
#clipped_ndvi_cropped <- crop(clipped_ndvi, common_extent)
soil_resampled_cropped <- crop(soil_resampled, common_extent)
dem_cropped <- crop(dem, common_extent)
ndvi_cropped <- crop(ndvi, common_extent)
clipped_tree_cropped <- crop(clipped_tree, common_extent)
clipped_bioclim_stack_cropped <- crop(clipped_bioclim_stack, common_extent)

# make sure resolutions are all matching
# Define the reference raster (e.g., moist_cropped)
reference_raster <- clipped_bioclim_stack_cropped

# Resample all rasters to match the reference raster
moist_resampled <- resample(moist_cropped, reference_raster, method = "bilinear")
ndvi_resampled <- resample(ndvi_cropped, reference_raster, method = "bilinear")
soil_resampled <- resample(soil_resampled_cropped, reference_raster, method = "bilinear")
dem_resampled <- resample(dem_cropped, reference_raster, method = "bilinear")
tree_resampled <- resample(clipped_tree_cropped, reference_raster, method = "bilinear")
bioclim_resampled <- resample(clipped_bioclim_stack_cropped, reference_raster, method = "bilinear")

# Combine the resampled rasters into a single stack
bio_stack <- stack(moist_resampled, 
                   ndvi_resampled, 
                   soil_resampled, 
                   dem_resampled, 
                   tree_resampled, 
                   bioclim_resampled)

# View the combined raster stack
plot(bio_stack)
bio_stack
# change names
# Check current names
current_names <- names(bio_stack)
print(current_names)

## Rename layers 1, 3, and 4
new_names <- current_names
new_names[1] <- "moist"
new_names[3] <- "soil"
new_names[4] <- "elevation"

# Assign the new names to the raster stack
names(bio_stack) <- new_names
# Verify changes
print(names(bio_stack))
# Save the combined raster stack with new layer names
writeRaster(bio_stack, "csiro_clim_vars.tif", format = "GTiff", overwrite = TRUE)

# Import the stacked raster
mac_preds <- stack("E:/EGK/Sdm_project/rdata/csiro_clim_vars.tif")
mac_preds <- bio_stack
mac_preds
rm(mac_preds2)
# Define new names
#current_names <- c("ndvi", "soil", "elevation", "treecover", "Max_Mth_Max_Temp")

# Assign new names to the raster stack
names(mac_preds) <- current_names

# Check the names to confirm they have been updated
print(names(mac_preds))
#========================read in again
# # Import the stacked raster
# bio_stack <- stack("D:/EGK/project/csiro_clim_vars.tif")
# bio_stack
# # Current layer names
 current_names <- names(mac_preds) #bio_stack
# 
# # Change names for specific layers
 current_names[1] <- "moisture"   
 current_names[2] <- "ndvi"
 current_names[3] <- "soil"       
 current_names[4] <- "elevation"  
 current_names[5] <- "treecover"
 current_names[6] <- "Annual_Precip" 
 current_names[7] <- "Min_Mth_Precip"
 current_names[8] <- "Max_Mth_Precip"  
current_names[9] <- "Min_Mth_Tmp"
 current_names[10] <- "Annual_Min_Temp"
 current_names[11] <- "Min_Mth_Min_Temp"
current_names[12] <- "Min_Mth_Max_Temp"
 current_names[13] <- "Annual_Max_Temp"
 current_names[14] <-"Max_Mth_Max_Temp"    
# 
# 
# # Assign the new names to the raster stack
 names(mac_preds) <- current_names
# 
# # Verify the changes
names(mac_preds)


##-------------- check colinearity of variables

#Extract the values from the raster stack as a data frame
raster_values <- as.data.frame(values(mac_preds))
#install.packages("usdm")
library(usdm)
# Compute VIFs
vif_result <- vifcor(raster_values, th = 0.7)

# View the result
print(vif_result)

# Compute the pairwise plot
pairs(raster_values, type = "l")
library(ggplot2)
library(GGally) 
ggpairs(raster_values)

# remove collinear variables from rasterstack of the variables
# List of layers to drop
layers_to_drop <- c("Annual_Precip", "Annual_Min_Temp", "moisture", 
                    "Min_Mth_Precip", "Min_Mth_Max_Temp", "Min_Mth_Min_Temp",
                    "Annual_Max_Temp", "Min_Mth_Tmp")

# Drop the layers
macpreds <- dropLayer(mac_preds, layers_to_drop)
macpreds
# Check the names of the remaining layers
names(my_raster_stack)
# or 
mac_preds <- exclude (mac_preds,vif_result)
mac_preds
names(mac_preds)
writeRaster(mac_preds, filename = "macropod_preds.tif", overwrite = TRUE)


