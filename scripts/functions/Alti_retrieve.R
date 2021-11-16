# Fonctions basees selon le modele des fonctions d'extraction de Yves Bas (https://github.com/cesco-lab/Vigie-Chiro_scripts/tree/master/functions/extractGI) 


##  # pre-requis : formation du raster altitude/pente/aspect de la pente
########################################################################

##  library( raster)
##  library(sf)
##  library(rasterVis)
##  library(ggplot2)
##  
##  
##  fra.adm <- st_transform(st_read(dsn = "C:/Users/Travail/Desktop/Ressource QGis/france/adm/FRA_adm2.shp"),crs=2154)
##  dsnTable = "C:/git/epoc/DS.v2/epoc_barycentre_liste_density_add.csv"
##  Point_sf <- st_as_sf(x = Point, coords = Coords, crs = 2154)
##  asc_list <- list.files(path = "C://Users/Travail/Desktop/Ressource QGis/france/altitude/BDALTIV2_2-0_75M_ASC_LAMB93-IGN69_FRANCE_2020-09-30",
##                         pattern='asc', full.names=TRUE, recursive = T)
##  
##  
##  # formation du raster altitude a partir des fichiers .ascii
##  rast_list <- list()
##  for(i in 1:length(asc_list)) { rast_list[i] <- raster(asc_list[i]) }
##  rast_list$fun <- mean
##  
##  # collage des dalles raster
##  ALTI <- do.call(mosaic,rast_list)
##  crs(ALTI) <- crs(Point_sf)
##  
##
##  ALTI_pente <- raster::terrain(ALTI, opt=c("slope","aspect"),unit="degrees",neighbors=8)
##
##  h <- stack(c(ALTI,ALTI_pente))
##
##
##  # sauvegarde sur disque
##  # writeRaster(ALTI, filename = "C:/Users/Travail/Desktop/Ressource QGis/france/Altitude.tif", format = "GTiff")
##  # writeRaster(h, filename = "C:/Users/Travail/Desktop/Ressource QGis/france/Altitude_pente_orientation.tif", format = "raster")



## dsnTable = "C:/git/epoc/DS.v2/epoc_barycentre_liste_density_add.csv"
## names_coord = c("X_barycentre_L93","Y_barycentre_L93")
## buffer_medium = 500
## buffer_large = 5000




extract_Alti = function(dsnTable, names_coord, buffer_medium, buffer_large, prefixe_fichier)
{
  # packages
  library(sp)
  library(sf)
  library(raster)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(progress)
  library(data.table)
  
  
    
  # recuperation des donnees
    cat("\t-------\tRecuperation des donnees geolocalisees\t-------\n")
    
    Point <- as.data.frame(fread(dsnTable,header=T,stringsAsFactors = F,encoding="UTF-8"))
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
    
    
    # recuperation des donnes d'altitude BD - Alti (multiple .asc --> 1 raster) ----
    ALTI <- stack("C:/Users/Travail/Desktop/Ressource QGis/france/Altitude_pente_orientation.tif")

    
  # extraction de l'altitude/pente/orientation -----
  cat("\t-------\tDebut extraction de l'altitude\t-------\n")
    
    # locale
    cat("\t---\tAltitude/pente/orientation locale\t---\n")
    
    SpAltis <- raster::extract(ALTI, Point_sf,df=T,along=T,nl=3)
    colnames(SpAltis) <- c("id","SpAltis","SpPens","SpPenOrs")
    
    
    # buffer medium
    cat(paste0(c("\t---\tAltitude/pente/orientation moyenne dans un buffer de :",BuffMed," m","\t---\n")))

    SpAltim <- exactextractr::exact_extract(x=ALTI, y=PointBuffMed, fun="mean")
    colnames(SpAltim) <- c("SpAltim","SpPenm","SpPenOrm")
    SpAltim$id <- c(1:nrow(SpAltim))
    
    
    # buffer large
    cat(paste0(c("\t---\tAltitude/pente/orientation moyenne dans un buffer de :",BuffLarg," m","\t---\n")))
    
    SpAltil <- exactextractr::exact_extract(x=ALTI, y=PointBuffLarg, fun="mean")
    colnames(SpAltil) <- c("SpAltil","SpPenl","SpPenOrl")
    SpAltil$id <- c(1:nrow(SpAltil))
  
  # calcul des pentes (conservion degrees en metre)
    cat(paste0(c("\t---\tMise en commun des tables et calcul des pentes en mètres\t---\n")))
    # mise en commun des tables
    SpAlti_all <- merge(SpAltis,SpAltim,all=T)
    SpAlti_all <- merge(SpAlti_all,SpAltil,all=T)
  
    SpAlti_all[,grep("Pen[a-z]{1}$",colnames(SpAlti_all))] <- apply(X = SpAlti_all[,grep("Pen[a-z]{1}$",colnames(SpAlti_all))],
                                                           MARGIN = 2,
                                                           FUN = function(X){ res(ALTI)[1] * tan(X * pi/180)}) # (75m : resolution du raster) * hauteur determine par l'angle en degrees avec une conversion degrees -> radians
    SpAlti_all <- dplyr::left_join(SpAlti_all, Point[,c("id",Coords,"ID_liste")]) # add d'informations sur les listes
    
    
    save(SpAlti_all, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"ALTI_envEPOC.RData")) # securite
    write.csv(SpAlti_all, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"ALTI_envEPOC.csv"), row.names = F)
    
    cat("\n\n\t-------\tExtract ALTI done - check /function_output \t-------\n\n")
    
    
}









