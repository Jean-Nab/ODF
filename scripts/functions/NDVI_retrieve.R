

# data source earthdata.nasa.gov
# citation : Didan, K. (2015). MOD13Q1 MODIS/Terra Vegetation Indices 16-Day L3 Global 250m SIN Grid V006 [Data set].
#            NASA EOSDIS Land Processes DAAC. Accessed 2021-09-14 from https://doi.org/10.5067/MODIS/MOD13Q1.006

# merging des data pour formation du raster aggrege & projection sur jeu de data ? (ou need de projection des points sur crs)
## 
## library(raster)
## library(gdalUtils)
## library(rgdal)
## # library(terra)
## 
## # recup' chemin des .hdf
##   hdf_src <- list.files(path = "C:/Users/Travail/Desktop/Ressource QGis/NDVI/raw_data/",pattern = ".hdf",full.names = T)
## 
## # conversion .hdf -> Gtiff uniquement du layer NDVI
##   for (i in 1:length(hdf_src)){
##     
##     sds <- get_subdatasets(hdf_src[i]) # recup' des layer du hdf
##     gdal_translate(sds[grep("NDVI$",sds)], dst_dataset = paste0("C:/Users/Travail/Desktop/Ressource QGis/NDVI/processed/", # selection du layer NDVI
##                                                                 gsub(".{1,}/raw_data/","",gsub(".hdf",".tif",hdf_src[i]))) # save en .tif
##                    )
##     
##   }
## 
## # merging des tiles de raster
##   tif_src <- list.files(path = "C:/Users/Travail/Desktop/Ressource QGis/NDVI/processed/",pattern = ".tif",full.names = T)
##   mosaic_rasters(gdalfile = tif_src, dst_dataset = "C:/Users/Travail/Desktop/Ressource QGis/france/NDVI_France.tif")
##   
##   j <- raster("C:/Users/Travail/Desktop/Ressource QGis/france/NDVI_France.tif")
##   
##   
## # reprojection / unscaling ?
##   j_l93 <- projectRaster(from = j,
##   crs = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
## ")
##   j_l93 <- j_l93*0.0001
## 
##   writeRaster(j_l93, filename = "C:/Users/Travail/Desktop/Ressource QGis/france/NDVI_France_l93.tif")
## 
## 
## 



# dsnTable = "C:/git/ODF/data/EPOC-ODF2021_barycentre.csv"
# names_coord = c("X_barycentre_L93","Y_barycentre_L93")
# buffer_medium = 500
# prefixe_fichier = "EPOC-ODF2021_"
# dsnRaster = "C:/Users/Travail/Desktop/Ressource QGis/france/NDVI_France_l93.tif"

extract_NDVI = function(dsnTable, dsnRaster,  names_coord, buffer_medium, prefixe_fichier){
  
# packages
  library(sp)
  library(sf)
  library(raster)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  
# recuperation des donnees
  cat("\t-------\tRecuperation des donnees geolocalisees\t-------\n")
  
  Point <- read.csv(dsnTable)
  Coords <- names_coord
  BuffMed <- buffer_medium
  
  # Retrait des NAs
  testCoords=match(Coords,names(Point))
  Point=subset(Point,!is.na(as.data.frame(Point)[,testCoords[1]]))
  Point=subset(Point,!is.na(as.data.frame(Point)[,testCoords[2]]))
  
  Point$id <- c(1:nrow(Point))
  Point$ID_extract <- c(rep(1:nrow(Point)))
  
  # formation de l'objet spatial
  Point_sf <- st_as_sf(x = Point, coords = Coords, crs = 2154)
  
  PointBuffMed <- st_buffer(Point_sf, dist = BuffMed)
  
  
# recuperation des donnes d'altitude BD - Alti (multiple .asc --> 1 raster) ----
  NDVI <- raster(dsnRaster)
  
  
# extraction de l'altitude/pente/orientation -----
  cat("\t-------\tDebut extraction du NDVI\t-------\n")
  
# buffer medium
  cat(paste0(c("\t---\tNDVI moyen dans un buffer de :",BuffMed," m","\t---\n")))
  
  SpNDVIm <- exactextractr::exact_extract(x=NDVI, y=PointBuffMed, fun="mean")
  
  SpNDVIm <- as.data.frame(SpNDVIm)
  colnames(SpNDVIm) <- c("SpNDVIm")
  SpNDVIm$id <- c(1:nrow(SpNDVIm))
  
  
  SpNDVIm <- dplyr::left_join(SpNDVIm, Point[,c("id",Coords,"ID_liste")]) # add d'informations sur les listes
  
  save(SpNDVIm, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"NDVI_envEPOC.RData")) # securite
  write.csv(SpNDVIm, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"NDVI_envEPOC.csv"), row.names = F)
  
  cat("\n\n\t-------\tExtract NDVI done - check /function_output \t-------\n\n")
  
  
}

