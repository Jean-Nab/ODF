library(sp)
library(raster)
library(dplyr)
library(ggplot2)
library(maptools)


# data ----
  land <- raster("C:/git/epoc/data/clc2018_clc2018_v2018_20_raster100m/CLC2018_CLC2018_V2018_20.tif")
  fra.adm <- st_transform(st_read(dsn = "C:/Users/Travail/Desktop/Ressource QGis/france/adm/FRA_adm2.shp"),crs=2154)
  

  fra.adm <- st_transform(fra.adm,crs=crs(land))
# crop/mask CLC to france ----
  CLC_fra <- crop(land,fra.adm)
  CLC_fra1 <- mask(CLC_fra, st_buffer(fra.adm,dist=20000))


# projection en L93
  CLC_fra2 <- projectRaster(CLC_fra1, crs = CRS("+init=epsg:2154"))



writeRaster(CLC_fra2, filename = "C:/Users/Travail/Desktop/Ressource QGis/france/CLC2018_raster100m",format = "GTiff")


# V2 ----------
  source("http://raw.githubusercontent.com/brry/misc/master/shp2raster.R")
  shp_ras <- shp2raster("C:/Users/Travail/Desktop/Ressource QGis/france/CLC18_SHP__FRA_2019-08-21/1_DONNEES_LIVRAISON_2019-08-21/CLC2018_FR/CLC18_FR.shp", column="CODE_18",cellsize = 100)

# V3 ------
  # passage par Qgis --> need changement class de $CODE_18
    land.shp <- st_read("C:/Users/Travail/Desktop/Ressource QGis/france/CLC18_SHP__FRA_2019-08-21/1_DONNEES_LIVRAISON_2019-08-21/CLC2018_FR/CLC18_FR.shp",crs=2154)

    land.shp$CODE_18 <- as.integer(as.character(land.shp$CODE_18))

    st_write(land.shp, "C:/Users/Travail/Desktop/Ressource QGis/france/CLC18_SHP__FRA_2019-08-21/1_DONNEES_LIVRAISON_2019-08-21/CLC2018_FR/CLC18_FR_V2.shp", driver = "ESRI Shapefile")





