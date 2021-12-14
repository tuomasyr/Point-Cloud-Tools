# Load required packages
require(raster); require(ForestTools); require(rgdal);require(sp); 

# Read CHM
chm <- raster('./tmp/raster/chm_max.tif')
projection(chm) <- sp::CRS("+init=epsg:3067")

# Detect trees.
# Apply the Variable Window Filter algorithm to detect tree tops. https://doi.org/10.14358/PERS.70.5.589
lin <- function(x){x * 0.05 + 0.6}
ttops <- vwf(CHM = chm, winFun = lin, minHeight = 3); crs(ttops) <- crs(chm)

# Apply Marker-controlled Watershed Segmentation to segment tree crowns. https://doi.org/10.1016/1047-3203(90)90014-M
crowns_above<- mcws(treetops = ttops, CHM = chm,  minHeight = 2, verbose = FALSE)

# Polygonize raster output
crowns_above_Poly <- rasterToPolygons(crowns_above,dissolve = T); crs(crowns_above_Poly) <- crs(ttops)

# Export data
png(filename = './tmp/figures/chm_treeseg_above.png',width = 20, height = 20, units = "cm", res = 400); plot(chm); plot(ttops, col = 'red',add = T); plot(crowns_above_Poly, add = T); dev.off()
png('./tmp/figures/chm.png'); plot(chm); dev.off()
writeRaster(chm,'./tmp/raster/chm.tif',format = "GTiff", overwrite = T)
writeOGR(crowns_above_Poly, "./tmp/crownseg", "crownseg_above",driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(ttops, "./tmp/treeloc", "treelocs_above",driver="ESRI Shapefile", overwrite_layer = TRUE)