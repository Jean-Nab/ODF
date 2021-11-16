# Fonctions basees selon le modele des fonctions d'extraction de Yves Bas (https://github.com/cesco-lab/Vigie-Chiro_scripts/tree/master/functions/extractGI) 



## extract_OSO(dsnTable = "C:/git/epoc/DS.v2/epoc_barycentre_liste_density_add.csv",
##             names_coord = c("X_barycentre_L93","Y_barycentre_L93"),
##             buffer_medium = 50,
##             buffer_large = 500,
##             dsnRasterOSO = "C:/Users/Travail/Desktop/Ressource QGis/france/OCS_2018_CESBIO.tif")
## 
## 
## 
## dsnTable = "C:/git/epoc/DS.v2/epoc_barycentre_liste_density_add.csv"
## names_coord = c("X_barycentre_L93","Y_barycentre_L93")
## buffer_medium = 50
## buffer_large = 500
## dsnRasterOSO = "C:/Users/Travail/Desktop/Ressource QGis/france/OCS_2018_CESBIO.tif"



extract_OSO = function(dsnTable, names_coord, buffer_medium, buffer_large, dsnRasterOSO, prefixe_fichier)
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
    
    # recup' data/raster ----
    cat("\t----\tDebut du chargement des donnees\t----\n")
    
    Point <- as.data.frame(fread(dsnTable,header=T,stringsAsFactors = F,encoding="UTF-8"))
    Coords <- names_coord
    BuffMed <- buffer_medium
    BuffLarg <- buffer_large
    
    OSO <- raster(dsnRasterOSO)
    ##### test_stack <- stack(c(dsnRasterCLC,dsnRasterOSO))
    
    # retrait des NAs
    testCoords=match(Coords,names(Point))
    Point=subset(Point,!is.na(as.data.frame(Point)[,testCoords[1]]))
    Point=subset(Point,!is.na(as.data.frame(Point)[,testCoords[2]]))
    
    Point$id <- c(1:nrow(Point))
    Point$ID_extract <- c(rep(1:nrow(Point)))
    
    # formation de l'objet spatial
    Point_sf <- st_as_sf(x = Point, coords = Coords, crs = 2154)
    
    PointBuffMed <- st_buffer(Point_sf, dist = BuffMed)
    PointBuffLarg <- st_buffer(Point_sf, dist = BuffLarg)
  
    
    cat("\t----\tFin du chargement des donnees\t----\n")
  
  # EXTRACTION ----
    # Buffer medium ----
      cat(paste("\t----\tDebut extraction des donnees OSO selon le buffer :",BuffMed,"m\t----\n"))
      
      vec.iteration_med <- seq(from = 1, to = length(PointBuffMed$ID_extract), by = 305) # argument suplementaire ? / add condition : NA/all -> tous
      # vec.iteration_med <- vec.iteration_med[1:10]
      
      OSO_hab_wide_BM <- as.data.frame(matrix(nrow=0,ncol = 1))
      colnames(OSO_hab_wide_BM) <- "id"
      
      # initialisation de la barre de progression
      pb_m <- progress::progress_bar$new(total = length(vec.iteration_med),
                                         format = "Extraction des donnees OSO buffer medium [:bar] :percent    Tps ecoule = :elapsedfull",
                                         clear = F)
    
      
      pb_m$tick(0)
      Sys.sleep(1/20)
      
      # boucle de 305 en 305: ----
      for(i in vec.iteration_med){
        
        
        PointBuffMed.tmp <- PointBuffMed[i:min((i+304),length(PointBuffMed$ID_extract)),] # si add condition --> condition - 1
        
        
        OSO_hab.tmp <- raster::extract(x = OSO, 
                                       y = PointBuffMed.tmp,
                                       df=T,
                                       along=T)
        
        
        # formatage des donnees
        colnames(OSO_hab.tmp)[1] <- "ID_extract"
        OSO_hab.tmp$ID_extract <- OSO_hab.tmp$ID_extract + (i - 1)
        colnames(OSO_hab.tmp)[2] <- "OSO_hab"
        
        
        # gestion dtf des resultats preliminaires -----
        
        OSO_hab.tmp <- left_join(OSO_hab.tmp,st_drop_geometry(PointBuffMed.tmp[,c("id","ID_extract")]), by = "ID_extract")
        
        OSO_hab_wide.tmp <- OSO_hab.tmp %>% reshape2::dcast(ID_extract + id ~ OSO_hab, fun.aggregate = length, value.var = "id")
       
        OSO_hab_wide_BM <- merge(OSO_hab_wide_BM,OSO_hab_wide.tmp,all=T)
        
        
        # actualisation de la progression
        pb_m$tick()
        Sys.sleep(1 / 100)
        
      }
    
      ### OSO_hab_wide_BM <- OSO_hab_wide_BM[,-grep("FALSE|TRUE|,",colnames(OSO_hab_wide_BM))] # retrait d'un artefact de merge
      
      # correction des NAs (= habitats absent lors de l'extract du point)
      OSO_hab_wide_BM[is.na(OSO_hab_wide_BM)] <- 0
      
      # retrait des habitats 0 ==> zone non couverte par le raster [condition : detection de la colonne "0" du dtf]
      OSO_hab_wide_BM <- OSO_hab_wide_BM[,if(length(grep("^0$",names(OSO_hab_wide_BM)) > 0)){ -grep("^0$",names(OSO_hab_wide_BM))}else{names(OSO_hab_wide_BM)}]
      
      
      # calcul des proportions d'habitats 
      OSO_hab_wide_BM[,grep("[0-9]",colnames(OSO_hab_wide_BM))] <- ( OSO_hab_wide_BM[,grep("[0-9]",colnames(OSO_hab_wide_BM))] / 
                                                                       rowSums(OSO_hab_wide_BM[,grep("[0-9]",colnames(OSO_hab_wide_BM))])
      ) * 100
      
      # rearragement de l'ordre des colonnes
      OSO_hab_wide_BM <- OSO_hab_wide_BM %>%
        tidyr::pivot_longer(cols = colnames(OSO_hab_wide_BM)[grep("[0-9]{1,}",colnames(OSO_hab_wide_BM))]) %>%
        mutate(name = as.numeric(name)) %>%
        arrange(name) %>%
        tidyr::pivot_wider(names_from = name, values_from = value)
      
      colnames(OSO_hab_wide_BM) <- gsub("([0-9]{1,})","SpOSOm_\\1",colnames(OSO_hab_wide_BM))    
    
      
  
    # Buffer large -------
      cat(paste("\n\t----\tDebut extraction des donnees OSO selon le buffer :",BuffLarg,"m\t----\n"))
      
      vec.iteration_larg <- seq(from = 1, to = length(PointBuffLarg$ID_extract), by = 200) # argument suplementaire ? / add condition : NA/all -> tous
      # vec.iteration_larg <- vec.iteration_larg[1:10]
      
      OSO_hab_wide_BL <- as.data.frame(matrix(nrow=0,ncol = 1))
      colnames(OSO_hab_wide_BL) <- "id"
      
      # initialisation de la barre de progression
      pb_l <- progress::progress_bar$new(total = length(vec.iteration_larg),
                                         format = "Extraction des donnees OSO buffer large [:bar] :percent    Tps ecoule = :elapsedfull",
                                         clear = F)
      
      pb_l$tick(0)
      Sys.sleep(1/20)
      
      
      # boucle ----
      for(i in vec.iteration_larg){
        
        PointBuffLarg.tmp <- PointBuffLarg[i:min((i+199),length(PointBuffLarg$ID_extract)),] # si add condition --> condition - 1
        
        
        OSO_hab.tmp <- raster::extract(x = OSO, 
                                       y = PointBuffLarg.tmp,
                                       df=T,
                                       along=T)
        
        
        # formatage des donnees
        colnames(OSO_hab.tmp)[1] <- "ID_extract"
        OSO_hab.tmp$ID_extract <- OSO_hab.tmp$ID_extract + (i - 1)
        colnames(OSO_hab.tmp)[2] <- "OSO_hab"
        
        
        # gestion dtf des resultats preliminaires -----
        
        OSO_hab.tmp <- left_join(OSO_hab.tmp,st_drop_geometry(PointBuffLarg.tmp[,c("id","ID_extract")]), by = "ID_extract")
        
        OSO_hab_wide.tmp <- OSO_hab.tmp %>% reshape2::dcast(ID_extract + id ~ OSO_hab, fun.aggregate = length, value.var = "id")
        
        OSO_hab_wide_BL <- merge(OSO_hab_wide_BL,OSO_hab_wide.tmp,all=T)
        
        # actualisation de la progression
        pb_l$tick()
        Sys.sleep(1 / 100)
        
      }
      
      
      # correction des NAs (= habitats absent lors de l'extract du point)
      OSO_hab_wide_BL[is.na(OSO_hab_wide_BL)] <- 0
      
      # retrait des habitats 0 ==> zone non couverte par le raster [condition : detection de la colonne "0" du dtf]
      OSO_hab_wide_BL <- OSO_hab_wide_BL[,if(length(grep("^0$",names(OSO_hab_wide_BL)) > 0)){ -grep("^0$",names(OSO_hab_wide_BL))}else{names(OSO_hab_wide_BL)}]
      
      
      # calcul des proportions d'habitats 
      OSO_hab_wide_BL[,grep("[0-9]",colnames(OSO_hab_wide_BL))] <- ( OSO_hab_wide_BL[,grep("[0-9]",colnames(OSO_hab_wide_BL))] / 
                                                                       rowSums(OSO_hab_wide_BL[,grep("[0-9]",colnames(OSO_hab_wide_BL))])
      ) * 100
      
      # rearragement de l'ordre des colonnes
      OSO_hab_wide_BL <- OSO_hab_wide_BL %>%
        tidyr::pivot_longer(cols = colnames(OSO_hab_wide_BL)[grep("[0-9]{1,}",colnames(OSO_hab_wide_BL))]) %>%
        mutate(name = as.numeric(name)) %>%
        arrange(name) %>%
        tidyr::pivot_wider(names_from = name, values_from = value)
      
      colnames(OSO_hab_wide_BL) <- gsub("([0-9]{1,})","SpOSOl_\\1",colnames(OSO_hab_wide_BL))    
 
  
  
  # jointure des donnes CLC selon les 2 buffers
  
  OSO_hab_BMBL <- dplyr::left_join(OSO_hab_wide_BM,OSO_hab_wide_BL, by = c("id", "ID_extract"))
  
  OSO_hab_BMBL <- dplyr::left_join(Point[,c("id","ID_extract","ID_liste",names_coord)],
                                   OSO_hab_BMBL , by = c("id", "ID_extract"))
  
  OSO_hab_BMBL <- dplyr::left_join(OSO_hab_BMBL, Point[,c("id",Coords,"ID_liste")]) # add d'informations sur les listes
  
  # sauvegarde sur disque
  
  save(OSO_hab_wide_BM, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"OSO_BM_envEPOC.RData")) # securite
  save(OSO_hab_wide_BL, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"OSO_BL_envEPOC.RData")) # securite
  save(OSO_hab_BMBL, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"OSO_envEPOC.RData")) # securite
  write.csv(OSO_hab_BMBL, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"OSO_envEPOC.csv"), row.names = F)
  
  
  cat("\n\n\t-------\tExtract OSO done - check /function_output \t-------\n\n")
  
  
  }


