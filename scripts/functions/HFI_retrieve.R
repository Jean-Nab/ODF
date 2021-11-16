## # source : https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-2009-human-footprint/maps
## 
## # cropping du raster la france
## library(raster)
## hii_rast <- raster("C:/Users/Travail/Desktop/Ressource QGis/hii-global-geo-grid/wildareas-v3-2009-human-footprint.tif")
## 
## # load d'un extent sous wgs84
##   desired_extent <- extent(projectExtent(raster("C:/Users/Travail/Desktop/Ressource QGis/france/NDVI_France_l93.tif"), # extent france/espagne/UK/suisse
##                             "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") # crs wgs84
##               )
##   hii_rast_france <- crop(hii_rast, desired_extent) # cropping sous wgs84
## 
## # reprojection sous l93
##   hii_rast_france_l93 <- projectRaster(from = hii_rast_france,
##                                        crs = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
## 
## # save
##   writeRaster(hii_rast_france_l93, filename = "C:/Users/Travail/Desktop/Ressource QGis/france/HFI_France_l93.tif")
## 




# dsnTable = "C:/git/ODF/data/EPOC-ODF2021_barycentre.csv"
# names_coord = c("X_barycentre_L93","Y_barycentre_L93")
# buffer_medium = 500
# prefixe_fichier = "EPOC-ODF2021_"
# dsnRaster = "C:/Users/Travail/Desktop/Ressource QGis/france/HFI_France_l93.tif"





extract_HFI = function(dsnTable, dsnRaster,  names_coord, buffer_medium, prefixe_fichier){

# packages
  library(sp)
  library(sf)
  library(raster)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(data.table)
  
  # recuperation des donnees
  cat("\t-------\tRecuperation des donnees geolocalisees\t-------\n")
  
  Point <- as.data.frame(fread(dsnTable,header=T,stringsAsFactors = F,encoding="UTF-8"))
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
  HFI <- raster(dsnRaster)
  
  
  # extraction de l'altitude/pente/orientation -----
  cat("\t-------\tDebut extraction du HFI (Human footprint index\t-------\n")
  
  # buffer medium
  cat(paste0(c("\t---\tHFI moyen dans un buffer de :",BuffMed," m","\t---\n")))
  
  SpHFIm <- exactextractr::exact_extract(x=HFI, y=PointBuffMed, fun="mean")
  
  SpHFIm <- as.data.frame(SpHFIm)
  colnames(SpHFIm) <- c("SpHFIm")
  SpHFIm$id <- c(1:nrow(SpHFIm))
  
  
  SpHFIm <- dplyr::left_join(SpHFIm, Point[,c("id",Coords,"ID_liste")]) # add d'informations sur les listes
  
  save(SpHFIm, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"HFI_envEPOC.RData")) # securite
  write.csv(SpHFIm, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"HFI_envEPOC.csv"), row.names = F)
  
  cat("\n\n\t-------\tExtract HFI done - check /function_output \t-------\n\n")
  
  
  
  
  
}

























