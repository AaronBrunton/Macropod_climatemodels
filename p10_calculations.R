# Calculate p10 threshold from current models to transfer to future binary maps

### binary maps
library(terra)
library(raster)
library(sp)

# Convert SpatRaster to RasterLayer
egk_rf_accr45_rast <- raster(egk_rf_accr45)

# Extract occurrence values
occPredVals <- raster::extract(egk_rf_r45, egk_data)
egk_rf_r45
egk_rf
# convert spat to rast 
egk_rf_rast <- raster(egk_rf)
egk_rf_rast
occPredVals <- raster::extract(egk_rf_rast, egk_data)
# Swamp
sw_rf_rast <- raster(sw_rf_pred)
occPredVals <- raster::extract(sw_rf_rast, sw_data)
# RNW
rnw_rf_rast <- raster(rnw_rf_pred)
occPredVals <- raster::extract(rnw_rf_rast, rnw_data)
# POt
pot_rf_rast <- raster(pot_rf_pred)
occPredVals <- raster::extract(pot_rf_rast , pot_data)
# Red leg pad
rlpad_rf_rast <- raster(rlpad_rf_pred)
occPredVals <- raster::extract(rlpad_rf_rast , rlpad_data)
# red neck pad
rnpad_rf_rast <- raster(rnpad_rf_pred)
occPredVals <- raster::extract(rnpad_rf_rast , redneck_data)
# black stripe 
bsw_rf_rast <- raster(bsw_rf_pred)
occPredVals <- raster::extract(bsw_rf_rast , bsw_data)
 
# Calculate p10 threshold
if(length(occPredVals) < 10) {
  p10 <- floor(length(occPredVals) * 0.9)
} else {
  p10 <- ceiling(length(occPredVals) * 0.9)
}
thresh <- rev(sort(occPredVals))[p10]

print(paste("p10 Threshold:", thresh))
# egk = 0.474617
# swamp = 0.454053
# rnw = 0.4612
# pot =  0.419937
# rlpad = p10 Threshold: 0.53589
# rnpad = p10 Threshold: 0.476586666666667
# bsw = "p10 Threshold: 0.452833333333333"


