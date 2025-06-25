# Process a catalog to create all the things
#ctg
ctg<-readLAScatalog("day3/start")
crs(ctg) <- "EPSG:25832"
plot(ctg)

# Optional: create smaller tiles
#cut in smaller areas
# opt_output_files(ctg) <- paste0("day3/small/tile_{XLEFT}_{YBOTTOM}") # label outputs based on coordinates
# opt_chunk_buffer(ctg) <- 0
# opt_chunk_size(ctg) <- 250 # retile to 250 m
# small <- catalog_retile(ctg) # apply retile
# plot(small) # some plotting

#create the following folders: dtm, nlas, chm, chm2,nlas_segmented, segmented

#2.dtm
#tin
opt_output_files(ctg) <- paste0("day3/dtm/{*}_dtm") # set folder output
ctg_dtm <- rasterize_terrain(ctg, res = 5, algorithm = tin()) # normalize

#3.create nlas
opt_output_files(ctg) <- paste0("day3/nlas/", "/{*}_nlas") # set folder output
ctg_nlas <- normalize_height(ctg, tin())
# How to clip this to >= 0?

# Clipping function for map
clip_catalog <- function(las, ...) {
   output <- filter_poi(las, Z>=0)
   return(output)
}

opt_output_files(ctg_nlas) <- ''
ctg_nlas_clipped = catalog_map(ctg_nlas, clip_catalog)

#4. create CHM
ctg_nlas <- readLAScatalog("day3/nlas")
opt_filter(ctg_nlas) <- "-drop_z_below 0" # Doesn't seem to do anything??
ctg_nlas@data$Min.Z

opt_output_files(ctg_nlas) <- paste0("day3/chm/", "/{*}_chm")
crs(ctg_nlas) <- "EPSG:25832"
ctg_chm <- rasterize_canopy(ctg_nlas, res = 5, algorithm = dsmtin())

#5.locate trees
f <- function(x) {pmax(x * 0.1 + 3, 3)} # Window size from Z hacked for negative Z
f <- function(x) {(x * 0.1 + 3, 3)} # Window size from Z
heights <- seq(0,40,5)
ws <- f(heights)
plot(heights, ws, type = "l", ylim = c(0,6))
crs(ctg_nlas) <- "EPSG:25832"
opt_output_files(ctg_nlas) <- paste0("day3/chm/", "/{*}_ttops") # set folder output
ctg_ttops_function <- locate_trees(ctg_nlas, lmf(f))

#6.tree segmentation and tree metric in one file
# Not sure what the intent is here...
opt_output_files(ctg_nlas) <- paste0("day3/chm2", "/chm_{*}") #set folder output
chm <- rasterize_canopy(ctg_nlas, 1, dsmtin())

opt_output_files(ctg_nlas) <- paste0("day3/chm2", "/ttops_{*}")
# opt_output_files(ctg_nlas) <- ''
# I get errors here using bitmerge or gpstime
ttops <- locate_trees(ctg_nlas, lmf(f), uniqueness = "bitmerge")

opt_output_files(ctg_nlas) <- paste0("day3/nlas_segmented", "/{*}_segmented")
algo <- dalponte2016(chm, ttops)
ctg_segmented <- segment_trees(ctg_nlas, algo)

#
ctg_nlas_seg <- readLAScatalog("day3/nlas_segmented")
opt_output_files(ctg_nlas) <- paste0("day3/segmented", "/{*}_segmented") #set folder output
crowns_dalponte <- crown_metrics(ctg_nlas_seg, func = .stdmetrics_z, geom = "convex")
st_crs(crowns_dalponte) <- 25832
setwd("day3/segmented") #settare cartella output
st_write(
  crowns_dalponte,
  "trees_crown_and_metrics.shp",
  delete_layer = TRUE
)

#7. metrics in one file
#at raster level
ctg_nlas_seg <- readLAScatalog("day3/nlas_segmented")
stmetric_pixel <- pixel_metrics(ctg_nlas_seg, func = .stdmetrics_z, 10) # calculate all the z metrics at 10 m
crs(stmetric_pixel) <- "EPSG:25832"
setwd("day3/metrics/") #settare cartella output
writeRaster(stmetric_pixel, "metriche_z_raster.tiff")