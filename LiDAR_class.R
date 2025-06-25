# Day 1

library(tidyverse)
library(lidR)

# Read the point cloud
las = readLAS(here::here('data/las_area1.laz'))

# Check for problems
las_check(las)

# Interactive 3D plot via rgl
plot(las, bg='white', legend=TRUE)

# Color by attribute
plot(las, bg='white', color='Classification')
plot(las, bg='white', color='Intensity', legend=TRUE)

# Classification == 2 is ground
gnd <- filter_ground(las)
plot(gnd, size = 3, bg = "white") 

# mapview plot of the extent
plot(las, mapview=TRUE)

# Read just specific parts
readLAS(filter = "-help") # Lists the filter options
#-keep_first	Keep only the first return of each pulse
#-keep_last	Keep only the last return
#-keep_single	Keep returns from pulses with only one return
#-keep_intermediate	Keep returns that are neither first nor last
#-drop_first	Exclude first returns
#-drop_last	Exclude last returns
#-drop_z_below x
#-drop_z_above 900 NB is m above sea level

las = readLAS(here::here('data/las_area1.laz'), filter='-last_only')
las = readLAS(here::here('data/las_area1.laz'), filter='-keep_class 2')
plot(las)

#exercise 1
#open the las2
las = readLAS(here::here('data/las_area2.laz'))

#plot the point cloud
plot(las, bg='white')
plot(las, color="Intensity", bg="white", legend=TRUE ) #plot intensity
plot(las, color="ScanAngleRank", bg="white", legend=TRUE )

#how is the pulse density?
las
# density      : 10.69 points/m²
# density      : 8.94 pulses/m²

#check if the point cloud is fine
las_check(las)

#re-load the point cloud and plot it using only points below  1750
las = readLAS(here::here('data/las_area2.laz'), filter='-drop_z_above 1759')
plot(las, bg='white')

#plot specific information of the point cloud of the area above 1900m: try to plot the number of returns - never done but you can try! do it with a blue background and adding a legend
las = readLAS(here::here('data/las_area2.laz'), filter='-drop_z_below 1900')
plot(las, bg='white', color='NumberOfReturns', legend=TRUE)

las = readLAS(here::here('data/las_area2.laz'), filter='-keep_class 1 2')

# Use area1 for this
p1 <- c(679368, 5155974)
p2 <- c(679468, 5156074)
transect<-clip_transect(las, p1, p2, width = 5) #plot point cloud transect
plot(transect)
plot(transect, color='Classification')

# Plotting a transect with ggplot gives a cross-sectional view of the data
ggplot(payload(transect), aes(X,Y, color=Z  ) )+ 
  #you can plot also other infos color=Intensity, color=Classification
  geom_point()+
  theme_minimal()+
  scale_color_viridis_c(option="viridis")

ggplot(payload(transect), aes(X,Z, color=factor(Classification)) )+ 
  geom_point()+
  theme_minimal()

###
#exercise 2 to do at home
#create a transect with the new las2
#show the Classification of a two given points of transect of 10m
las2<-readLAS("C:/Users/MicTorresani/OneDrive - Scientific Network South Tyrol/varie/corsi/physalia/LiDAR 2025/data/starting_data/las/area2/DownloadService/las_area2.laz" ) 
p1 <- c(620207, 5156764)
p2 <- c(620257, 5156764)
transect<-clip_transect(las2, p1, p2, width = 5) #plot point cloud transect

ggplot(payload(transect), aes(X,Z, color=Z  ) )+
  geom_point()+
  theme_minimal()+
  scale_color_viridis_c(option="viridis")

#clip the las with the shapefile of a given polygon
library(sf)
area_test<-st_read(here::here("data/area_shp/area_shp_test.shp"))

plot(area_test)
clipped_las<-clip_roi(las, area_test)
ctr = st_centroid(area_test) |> st_coordinates()
plot(clipped_las)

#clip my las with given gps point and given radius
# Use QGIS to get point coordinates (right-click on point)
area_circle<-clip_circle(las, x=679368, y=5155974, radius = 50   )
p1 <- c(679368, 5155974)
p2 <- c(679468, 5156074)
area_circle<-clip_circle(las, x=p1[1], y=p1[2], radius = 50   )
area_circle<-clip_circle(las, x=ctr[1, 'X'], y=ctr[1, 'Y'], radius = 50   )
plot(area_circle)

# Classification - if not classified or badly classified
#classify point cloud, we use the function classify_ground ()
#PMF - not so good
classification1<-classify_ground(las, algorithm = pmf(ws=5, th=3) )

transect<-clip_transect(classification1, p1, p2, width = 5)
ggplot(payload(transect), aes(X,Z, color=Classification  ) )+
  geom_point()+
  ggtitle("pmf")+
  theme_minimal()+
  scale_color_viridis_c(option="viridis")

# PMF with tuning works better
ws<-seq(3, 12, 3)
ws
th<-seq(0.1, 2, length.out= length(ws))

classification2<-classify_ground(las, algorithm = pmf(ws=ws, th=th) )
transect<-clip_transect(classification2, p1, p2, width = 5)
ggplot(payload(transect), aes(X,Z, color=Classification  ) )+
  geom_point()+
  ggtitle("pmf_tuned")+
  theme_minimal()+
  scale_color_viridis_c(option="viridis")

#CSF - better classification
library(RCSF)

classification_csf<-classify_ground(las, algorithm = csf()) #based on Zang et al 2016. 
transect<-clip_transect(classification_csf, p1, p2, width = 5) #we can change the transect width if needed
ggplot(payload(transect), aes(X,Z, color=Classification  ) )+
  geom_point()+
  ggtitle("csf")+
  theme_minimal()+
  scale_color_viridis_c(option="viridis")

ggplot(payload(transect), aes(X,Z, color=factor(Classification)  ) )+
  geom_point()+
  ggtitle("csf")+
  theme_minimal()

ground = filter_ground(las)
plot(ground, bg='white')

# Digital Terrain Model (DTM) with TIN - preferred method
# DTM is a raster of the terrain height (e.g. ground)
dtm_tin = rasterize_terrain(las, res=1, algorithm=tin())
plot_dtm3d(dtm_tin)
plot(dtm_tin)

# Plot of points and DTM
x = plot(filter_first(las))
add_dtm3d(x, dtm_tin)

# Save for import to QGIS
library(raster)
writeRaster(dtm_tin, 'data/dtm_tin.tiff')

# Terrain with IDW - alternate method
dtm_idw = rasterize_terrain(las, res=1, algorithm=knnidw())
plot_dtm3d(dtm_idw)
plot(dtm_idw)

writeRaster(dtm_idw, 'data/dtm_idw.tiff')

#dtm_kriging no tuning
dtm_kriging<-rasterize_terrain(las, res = 1, algorithm = kriging())
#res is in meters, use classification_csf if you don't trust your original classification
plot_dtm3d(dtm_kriging)
plot(dtm_kriging)
#dtm_kriging k=10
dtm_kriging_k10<-rasterize_terrain(las, res = 1, algorithm = kriging(k=10))
plot_dtm3d(dtm_kriging_k10)

#dtm_kriging k=150
dtm_kriging_k150<-rasterize_terrain(las, res = 1, algorithm = kriging(k=150))
writeRaster(dtm_kriging_k150, "data/dtm_kriging_150.tiff")
plot_dtm3d(dtm_kriging_k150)
plot(dtm_kriging_k150)
plot(dtm_kriging_k150-dtm_tin)

#this function allows to create raster of slope and aspect 
dtm_products<-terrain(dtm_tin, v=c("slope", "aspect"), unit="degree" )
plot(dtm_products) #plot both
plot(dtm_products$slope, main="slope")
writeRaster(dtm_products$slope, 'data/slope.tiff')

# Normalization - point cloud of distance above ground
# Wrong way to do it - creates lots of negative values
nlas_wrong = las - dtm_tin
plot(nlas_wrong, legend=TRUE)

# Better way to do it
nlas = normalize_height(las, algorithm=tin())
plot(nlas, legend=TRUE)

nlas_filtered = filter_poi(nlas, Z>=0, Z<=40 ) # Remove negative and too large
plot(nlas_filtered, legend=TRUE)

#### Day 2
# Digital Surface Model (DSM) - raster of canopy height above sea level
# i.e. not normalized
dsm = rasterize_canopy(las, res=1, algorithm=dsmtin())
plot(dsm)

# Canopy Height Model (CHM) - raster of canopy height above ground
# Wrong way to make canopy model - subtract terrain from canopy
# OK to do this if you don't have a point cloud
chm_wrong = dsm - dtm_tin
plot(chm_wrong)

# Right way - use normalized point cloud
chm_right = rasterize_canopy(nlas_filtered, res=1, algorithm=dsmtin())
plot(chm_right)
writeRaster(chm_right, 'data/chm_right.tiff')

# Segment trees - find tree top points
# Use normalized point cloud
library(sf)
tree_top_5 = locate_trees(nlas_filtered, algorithm=lmf(ws=5))
plot(tree_top_5)
st_write(tree_top_5, 'data/tree_top_5.gpkg')

# Changing window size changes the trees found
tree_top_10 = locate_trees(nlas_filtered, algorithm=lmf(ws=10))
plot(tree_top_10)
st_write(tree_top_10, 'data/tree_top_10.gpkg')

# Use a function to determine window size
#create a function that makes a window size that increase from 3 (whan height is 0)
# and increase of 0.1 for each meter. you can change the 0.1 and increase etc 0.2...
f <- function(x) {x * 0.1 + 3} 

# Show the function
heights <- seq(0,40,5)
ws <- f(heights)
plot(heights, ws, type = "l", ylim = c(0,6))

# Find treetops with variable window size
tree_top_function_right<-locate_trees(nlas_filtered, lmf(f))
st_write(tree_top_function_right, 'data/tree_top_function.gpkg')

x <- plot(nlas_filtered, bg = "white")
add_treetops3d(x, tree_top_function_right)


#exercise 6
#normalize the las2 using the tin() algorithm. 
#filter the points as done before. 
las2 = readLAS(here::here('data/las_area2.laz'))
nlas2 = normalize_height(las2, algorithm=tin())
plot(nlas2, legend=TRUE)

nlas_filtered2 = filter_poi(nlas2, Z>=0, Z<=45 ) # Remove negative and too large
plot(nlas_filtered2, legend=TRUE)

dtm_tin2 = rasterize_terrain(las2, res=2, algorithm=tin())
plot(dtm_tin2)
writeRaster(dtm_tin2, 'data/dtm_tin2.tiff')

#exercise 7
#create a 2 new CHM with the area2 with a spatial resolution of 2m. one using the dsmtin() and the other using the p2r(subcircle = 0.2,na.fill = tin() ) 
chm_right2 = rasterize_canopy(nlas_filtered2, res=2, algorithm=dsmtin())
plot(chm_right2)
writeRaster(chm_right2, 'data/chm_right2.tiff', overwrite=TRUE)

chm_right_p2r_2 = rasterize_canopy(nlas_filtered2, res=2, 
                                   algorithm=p2r(subcircle = 0.2,na.fill = tin()))
plot(chm_right_p2r_2)
writeRaster(chm_right2, 'data/chm_right2.tiff', overwrite=TRUE)

#exercise 8
#load the area2
#create a nlas_filtered (with filter 0, 45) in the area 2
#locate trees using  the f function setting height seq(0,45,5)
tree_top_function_right2<-locate_trees(nlas_filtered2, lmf(f))
st_write(tree_top_function_right2, 'data/tree_top_function2.gpkg')

# Tree segmentation - find tree outlines
# DalPonte 2016
# Dalponte, M. and Coomes, D. A. (2016), Tree-centric mapping of forest carbon density from airborne laser scanning and hyperspectral data. Methods Ecol Evol, 7: 1236–1245. doi:10.1111/2041-210X.12575.
algorithm1 = dalponte2016(chm_right, tree_top_function_right)
trees_segmented_dalponte = segment_trees(nlas_filtered, algorithm1)
plot(trees_segmented_dalponte, color='treeID')

# Look at one tree (ID from QGIS)
tree1 = filter_poi(trees_segmented_dalponte, treeID==1985)
plot(tree1)

# Li 2012
# Li, W., Guo, Q., Jakubowski, M. K., & Kelly, M. (2012). A new method for segmenting individual trees from the lidar point cloud. Photogrammetric Engineering & Remote Sensing, 78(1), 75-84.
algorithm2 = li2012()
trees_segmented_li= segment_trees(nlas_filtered, algorithm2)
plot(trees_segmented_li, color='treeID')

# Crown metrics
# Shape of crowns
crown_dalponte = crown_metrics(trees_segmented_dalponte, func=.stdtreemetrics, 
                               geom='convex')
st_write(crown_dalponte, 'data/crown_dalponte.gpkg')


##### Day 3
crown_dalponte_all_metrics = crown_metrics(
  trees_segmented_dalponte, func=.stdmetrics_z, geom='convex')
st_write(crown_dalponte_all_metrics, 'data/crown_dalponte_all_metrics.gpkg')


#cloud information at plot level with a given point, or adding a polygon
# Finds z statistics for the cloud points at given locations
# In this case, cloud points within a radius of a location
test_points<-st_read("data/punctual_plot/punctual_plot.shp")
metric_plot<-plot_metrics(nlas_filtered, func=.stdmetrics_z, radius=20, test_points) #set radius of your plot

### Volume
#volume, use with caution!! This is simplistic
#volume at crown level
crown_dalponte_new<-crown_metrics(trees_segmented_dalponte, func = .stdmetrics_z, geom="convex") #you can use also concave
crown_dalponte_new$volume<-0.05*(crown_dalponte_new$zmax^2.5)/100  
#use with caution!!!! #Some sources suggest values between 0.03 and 0.07 depending on species.
plot(crown_dalponte_new["volume"])

#volume at pixel level
crown_dalponte_new<-pixel_metrics(trees_segmented_dalponte, func = .stdmetrics_z, 5) 
#you can use also concave
crown_dalponte_new$volume<-0.05*(crown_dalponte_new$zmax^2.5)/100  
#use with caution!!!! #Some sources suggest values between 0.03 and 0.07 depending on species.
plot(crown_dalponte_new["volume"])

### segment trees by only CHM
# tree location using a CHM and not a point cloud
# Use if point cloud is not available
chm_right<-rasterize_canopy(nlas_filtered, res=1, algorithm = dsmtin()) 
# re-create the CHM with 1m spatial resolution, the higher the better

tree_top_chm <- locate_trees(chm_right, lmf(f)) 
# locate tree using the chm and the f function

agorithm_tree_by_chm <- dalponte2016(chm_right, tree_top_chm) 
# use of the dalponte2016 algorithm to detect the trees

crowns <- agorithm_tree_by_chm()
plot(crowns, col = pastel.colors(2500))

