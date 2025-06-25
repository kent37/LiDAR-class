# Try Northampton data

noho = readLAS('data/2015_ME_MA_lidar_Job1198982/Job1198982_42072_30_64.laz')
noho = noho |> filter_duplicates() |> remove_noise()
las_check(noho)

plot(noho, bg='white', color='ReturnNumber')
ctr = c(348531.5,2940843.7)

circle = clip_circle(noho, x=ctr[1], y=ctr[2], radius = 500   )
plot(circle, bg='white', color='Classification', size=5, pal=rainbow)
plot(filter_ground(noho), bg='white')

noho_dtm_tin = rasterize_terrain(noho, res=2, algorithm=tin())
plot_dtm3d(noho_dtm_tin) # Crashes

# CHM - not giving good result
noho_nlas = normalize_height(noho, algorithm=tin())
noho_nlas_filtered = filter_poi(noho_nlas, Z>=0 ) # Remove negative
noho_chm_right = rasterize_canopy(noho_nlas_filtered, res=1, algorithm=dsmtin())
plot(noho_chm_right)

