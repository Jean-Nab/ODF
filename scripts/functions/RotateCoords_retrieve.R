

extract_multiCoords = function(dsnTable, names_coord, prefixe_fichier, NB_rotation)
{
  # packages
  library(sp)
  library(sf)
  library(spdep)
  library(dplyr)
  library(tidyr)
  library(progress)
  library(data.table)
  
  
  
  # recuperation des donnees
  cat("\t-------\tRecuperation des donnees geolocalisees\t-------\n")
  
  Point <- as.data.frame(fread(dsnTable,header=T,stringsAsFactors = F,encoding="UTF-8"))
  Coords <- names_coord
  
  # Retrait des NAs
  testCoords=match(Coords,names(Point))
  Point=subset(Point,!is.na(as.data.frame(Point)[,testCoords[1]]))
  Point=subset(Point,!is.na(as.data.frame(Point)[,testCoords[2]]))
  
  Point$id <- c(1:nrow(Point))
  Point$ID_extract <- c(rep(1:nrow(Point)))
  
  # formation de l'objet spatial
  Point_sf <- as_Spatial(st_transform(st_as_sf(x = Point, coords = Coords, crs = 2154), crs = 4326))
  
  # recup des coordonne en wgs84
  CoordX <- Point_sf@coords[,1]
  CoordY <- Point_sf@coords[,2]
  CoordXY=as.matrix(cbind(CoordX,CoordY))
  
  Coord_recueil <- Point[,c("ID_liste","id","ID_extract")]
  
  # Calcul des coordonnees de l'axe X apres rotation
  cat(paste0(c("\t---\tCalcul des coordonnées de longitude après rotation\t---\n")))
  
  for (a in 0:(as.numeric(NB_rotation)-1)){
    Coordi=Rotation(CoordXY,angle=pi*a/as.numeric(NB_rotation))
    #print(plot(Coordi[,1],CoordDS[,1],main=as.character(a)))
    #print(plot(Coordi[,1],CoordDS[,2],main=as.character(a)))
    Coord_recueil=cbind(Coord_recueil,Coordi[,1])
    names(Coord_recueil)[ncol(Coord_recueil)]=paste0("SpCoordXrotation_",a)
  }
  
  
  save(Coord_recueil, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"Coords_envEPOC.RData")) # securite
  write.csv(Coord_recueil, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"Coords_envEPOC.csv"), row.names = F)
  
  cat("\n\n\t-------\tExtract Coordonnees done - check /function_output \t-------\n\n")
  
  
}













