// google earth engine code for the download of gedi global chm 
// https://lpdaac.usgs.gov/documents/998/GEDI02_UserGuide_V21.pdf
// (potapov 30m, lang10m, meta1m): 
// Go to code.earthengine.google.com and draw a rectangle ROI then Run
// =======================
// 1. Define AOI: draw or use a fixed polygon
// =======================
// You can also use Map.drawingTools() to draw it interactively
var aoi = geometry;
Map.centerObject(aoi, 13);
// =======================
// 2. Load Lang et al. CHM (10m) - 2020
// https://langnico.github.io/globalcanopyheight/
// =======================
var langCHM = ee.Image('users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1')
                .clip(aoi);
Map.addLayer(langCHM, {min: 0, max: 50, palette: ['white', 'green', 'darkgreen']}, 'Lang CHM (10m)');
// =======================
// 3. Load Potapov et al. CHM (GEDI V2.7) (30m) - 2019
// https://glad.umd.edu/dataset/gedi
// =======================
var potapovCHM = ee.ImageCollection('users/potapovpeter/GEDI_V27')
                    .mosaic()
                    .clip(aoi);
Map.addLayer(potapovCHM, {min: 0, max: 50, palette: ['white', 'orange', 'brown']}, 'Potapov CHM (30m)');
// =======================
// 4. Load Meta CHM (1m) if available - 2023
// https://gee-community-catalog.org/projects/meta_trees/
// https://sustainability.atmeta.com/blog/2024/04/22/using-artificial-intelligence-to-map-the-earths-forests/
// =======================
var metaCHM = ee.ImageCollection("projects/meta-forest-monitoring-okw37/assets/CanopyHeight")
                .filterBounds(aoi)
                .mosaic()
                .clip(aoi);
Map.addLayer(metaCHM, {min: 0, max: 50, palette: ['white', 'blue', 'purple']}, 'Meta CHM (1m)');
// =======================
// 5. Export to Google Drive
// =======================
// Lang export (10 m)
Export.image.toDrive({
  image: langCHM,
  description: 'Lang_CHM_10m',
  folder: 'GEE_exports',
  region: aoi,
  scale: 10,
  maxPixels: 1e13
});
// Potapov export (30 m)
Export.image.toDrive({
  image: potapovCHM,
  description: 'Potapov_CHM_30m',
  folder: 'GEE_exports',
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});
// Meta export (1 m)
Export.image.toDrive({
  image: metaCHM,
  description: 'Meta_CHM_1m',
  folder: 'GEE_exports',
  region: aoi,
  scale: 1,
  maxPixels: 1e13
});






