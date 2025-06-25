#biodiversity
chm_right<-rasterize_canopy(nlas_right_filtered, res=20, algorithm = dsmtin())
dtm_tin<-rasterize_terrain(las, res = 1, algorithm = tin()) #res is in meters, use classification_csf if yoo dont trust your original classification

#rao's Q index
mat_values <- values(chm_right)
mat_values_clean <- na.omit(mat_values[, "Z"])

dist_matrix <- as.matrix(dist(mat_values_clean))
#rao's Q result
rao_q <- sum(dist_matrix) / (length(mat_values_clean)^2)

mean_val <- mean(mat_values_clean, na.rm = TRUE)
sd_val <- sd(mat_values_clean, na.rm = TRUE)
#coefficient of variation 
cv_val <- (sd_val / mean_val) * 100

#do the same with dtm for geodiversity