# Fonctions basees selon le modele des fonctions d'extraction de Yves Bas (https://github.com/cesco-lab/Vigie-Chiro_scripts/tree/master/functions/extractGI) 

##  # pre-requis : crop bioclim selon france + stacking/save sur disque
#######################################################################
##  library(sf)
##  library(raster)
##  
##  
##  fra.adm <- st_read(dsn = "C:/Users/Travail/Desktop/Ressource QGis/france/adm/FRA_adm0.shp")
##  bioclim <- stack(list.files(path = "C:/Users/Travail/Desktop/Ressource QGis/Bioclim_world/", pattern = "wc2.1_30s",full.names = T))
##   
##  fra.buff <- st_transform(fra.adm,crs=2154) ; fra.buff <- st_buffer(fra.buff, 25000)
##  fra.buff <- st_transform(fra.buff, crs = 4326)
##  
##  fra.bioclim <- crop(x = bioclim, y = as_Spatial(fra.buff))
##  
##  fra.bioclim2 <- projectRaster(from = fra.bioclim,
##  crs = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
##  method = "ngb"
##  )
##  
##  
##  writeRaster(fra.bioclim2, filename = "C:/Users/Travail/Desktop/Ressource QGis/france/bioclim/bioclim_30s_degrees.grd", format = "raster")



 dsnTable = "C:/git/epoc/DS.v2/epoc_barycentre_liste_density_add.csv"
 names_coord = c("X_barycentre_L93","Y_barycentre_L93")
 buffer_medium = 500
 buffer_large = 5000
 dsnRaster = "C:/Users/Travail/Desktop/Ressource QGis/france/bioclim/bioclim_30s_degrees.grd"



extract_Bioclim = function(dsnTable, names_coord, buffer_medium, buffer_large, dsnRaster){
  # packages
  library(sp)
  library(sf)
  library(raster)
  library(dplyr)

  
  
  
  # recuperation des donnees
  cat("\t-------\tRecuperation des donnees geolocalisees\t-------\n")
  
  Point <- read.csv(dsnTable)
  Coords <- names_coord
  BuffMed <- buffer_medium
  BuffLarg <- buffer_large
  
  # Retrait des NAs
  testCoords=match(Coords,names(Point))
  Point=subset(Point,!is.na(as.data.frame(Point)[,testCoords[1]]))
  Point=subset(Point,!is.na(as.data.frame(Point)[,testCoords[2]]))
  
  Point$id <- c(1:nrow(Point))
  Point$ID_extract <- c(rep(1:nrow(Point)))
  
  # formation de l'objet spatial
  Point_sf <- st_as_sf(x = Point, coords = Coords, crs = 2154)
  
  PointBuffMed <- st_buffer(Point_sf, dist = BuffMed)
  PointBuffLarg <- st_buffer(Point_sf, dist = BuffLarg)
  
  
  # recuperation des donnes bioclim de france ----
  BIOC <- stack(dsnRaster)
  
  
  # extraction des variables bioclim -----
  cat("\t-------\tDebut extraction des variables bioclim\t-------\n")
  
  # buffer medium
    cat(paste0(c("\t---\tVariables bioclim moyenne dans un buffer de :",BuffMed," m","\t---\n")))
    
    SpBiom <- exactextractr::exact_extract(x=BIOC, y=PointBuffMed, fun="mean")
    
    # changement nom de colonnes
      colnames(SpBiom) <- gsub("(mean.wc2.1_30s_bio_)([0-9]{1,})","SpBiom_\\2",colnames(SpBiom))
      SpBiom$id <- c(1:nrow(SpBiom))
      SpBiom <- SpBiom[,order(colnames(SpBiom))]
    
  
  # buffer large
    cat(paste0(c("\t---\tVariables bioclim moyenne dans un buffer de :",BuffLarg," m","\t---\n")))
    
    SpBiol <- exactextractr::exact_extract(x=BIOC, y=PointBuffLarg, fun="mean")
    
    # changement nom de colonnes
      colnames(SpBiol) <- gsub("(mean.wc2.1_30s_bio_)([0-9]{1,})","SpBiol_\\2",colnames(SpBiol))
      SpBiol$id <- c(1:nrow(SpBiol))
      SpBiol <- SpBiol[,order(colnames(SpBiol))]

  # Jointure des tables ----
    SPBIO <- dplyr::left_join(SpBiom,SpBiol, by = "id")
  
    write.csv(SPBIO, file = "C:/git/ODF/output/function_output/BIOCLIM_envEPOC.csv", row.names = F)
    
    
    cat("\n\n\t-------\tExtract BIOCLIM done - check /function_output \t-------\n\n")
  
}
























