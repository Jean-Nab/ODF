extract_eauSurface = function(dsnTable, names_coord,buffer_small, buffer_medium, buffer_large, dsnSurf,
                               prefixe_fichier, size_batch_small, size_batch_med, size_batch_large){
  
  
  # packages
  library(sp)
  library(sf)
  library(raster)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(progress)
  # requiert library(data.table)
  
  
  # recuperation des donnees
  cat("\t-------\tRecuperation des donnees geolocalisees\t-------\n")
  
  Point <- read.csv(dsnTable)
  Coords <- names_coord
  BuffSmal <- buffer_small
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
  
  PointBuffSmal <- st_buffer(Point_sf, dist = BuffSmal)
  PointBuffMed <- st_buffer(Point_sf, dist = BuffMed)
  PointBuffLarg <- st_buffer(Point_sf, dist = BuffLarg)
  
  # load des shp
  cat("\t-------\tChargement des shapes SURFACES d'eau\t-------\n")

  SURF <- st_read(dsn = dsnSurf)
  SURF$Nature <- as.character(SURF$Nature)
  
  
  
  # exttraction proportion de surface des zones d'eaux -----
  
  #############################
  # Buffer locaux : small -----
  
  cat(paste0(c("\t---\tSurfaces d'eau dans un buffer de :",BuffSmal," m","\t---\n")))
  
  vec.iteration_smal <- seq(from = 1, to = length(PointBuffSmal$ID_extract), by = size_batch_small) # argument suplementaire ? / add condition : NA/all -> tous
  
  SurfSmal <- as.data.frame(matrix(nrow=0,ncol = 1))
  colnames(SurfSmal) <- c("id")
  class(SurfSmal$Nature) <- class(SURF$Nature)
  
  Rebus <- data.table::data.table()
  
  # initialisation de la barre de progression
  pb_s <- progress::progress_bar$new(total = length(vec.iteration_smal),
                                     format = "Extraction Aire des surfaces d'eau buffer small [:bar] :percent    Tps ecoule = :elapsedfull",
                                     clear = F)
  
  pb_s$tick(0)
  Sys.sleep(1/20)
  
  
  for(i in vec.iteration_smal){
    Point.tmp <- PointBuffSmal[i:min((i+size_batch_small-1),length(PointBuffSmal$ID_extract)),]
    
    SurfSmal.tmp <- st_intersection(SURF, Point.tmp) # intersection en block
    
    for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
      
      if(nrow(SurfSmal.tmp[SurfSmal.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
        SurfSmal.tmp[SurfSmal.tmp$id == h,"Aire_intersect"] <- as.numeric(st_area(SurfSmal.tmp[SurfSmal.tmp$id == h,]))
        
      }else{ # en cas d'absence de route --> formation table de donnees nulles
        
        Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Aire_intersect = 0, Nature = NA)))
      }
      
    }
    
    # jointure des tables de la condition presence de routes
    SurfSmal <- rbind(SurfSmal,st_drop_geometry(SurfSmal.tmp[,c("id","Nature","Aire_intersect")]))
    
    # actualisation de la progression
    pb_s$tick()
    Sys.sleep(1 / 100)
  }
  
  # gestion de l'extractions des longueurs de Surfs intersectees
  SurfSmal <- rbind(SurfSmal,Rebus)
  
  # calcul de la proportion des surfaces
  Aire_tot <- BuffSmal^2*pi
  SurfSmal$Aire_intersect <- SurfSmal$Aire_intersect / Aire_tot
  
  # conversion long -> wide
  SurfSmal.dt <- data.table::data.table(SurfSmal)
  SurfSmal_wide <- data.table::dcast(SurfSmal.dt,id + Aire_intersect ~ Nature, value.var = "Aire_intersect", fill = 0,fun.aggregate = sum) # conversion long -> wide
  
  
  # somme des longueurs de Surfs intersectes par id selon le type
  Extr_Smal <-  SurfSmal_wide %>%
    group_by(id) %>%
    mutate(SpCESs_1 = if(length(grep("Eau salée permanente",colnames(.))) > 0){sum(`Eau salée permanente`)}else{sum(`NA`)}) %>%
    mutate(SpCESs_2 = if(length(grep("Eau douce permanente",colnames(.))) > 0){sum(`Eau douce permanente`)}else{sum(`NA`)}) %>%
    mutate(SpCESs_3 = if(length(grep("Eau salée non permanente",colnames(.))) > 0){sum(`Eau salée non permanente`)}else{sum(`NA`)}) %>%
    mutate(SpCESs_4 = if(length(grep("Eau douce non permanente",colnames(.))) > 0){sum(`Eau douce non permanente`)}else{sum(`NA`)}) %>%
    mutate(SpCESs_5 = if(length(grep("Névé, glacier",colnames(.))) > 0){sum(`Névé, glacier`)}else{sum(`NA`)})
  
  
  Extr_Smal <- unique(Extr_Smal[,grep("id|SpCES",colnames(Extr_Smal))])
  
  
  
  #############################
  # Buffer locaux : medium -----
  
  cat(paste0(c("\t---\tSurfaces d'eau dans un buffer de :",BuffMed," m","\t---\n")))
  
  vec.iteration_med <- seq(from = 1, to = length(PointBuffMed$ID_extract), by = size_batch_med) # argument suplementaire ? / add condition : NA/all -> tous
  
  SurfMed <- as.data.frame(matrix(nrow=0,ncol = 1))
  colnames(SurfMed) <- c("id")
  class(SurfMed$Nature) <- class(SURF$Nature)
  
  Rebus <- data.table::data.table()
  
  # initialisation de la barre de progression
  pb_m <- progress::progress_bar$new(total = length(vec.iteration_med),
                                     format = "Extraction Aire des surfaces d'eau buffer medium [:bar] :percent    Tps ecoule = :elapsedfull",
                                     clear = F)
  
  pb_m$tick(0)
  Sys.sleep(1/20)
  
  
  for(i in vec.iteration_med){
    Point.tmp <- PointBuffMed[i:min((i+size_batch_med-1),length(PointBuffMed$ID_extract)),]
    
    SurfMed.tmp <- st_intersection(SURF, Point.tmp) # intersection en block
    
    for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
      
      if(nrow(SurfMed.tmp[SurfMed.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
        SurfMed.tmp[SurfMed.tmp$id == h,"Aire_intersect"] <- as.numeric(st_area(SurfMed.tmp[SurfMed.tmp$id == h,]))
        
      }else{ # en cas d'absence de route --> formation table de donnees nulles
        
        Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Aire_intersect = 0, Nature = NA)))
      }
      
    }
    # jointure des tables de la condition presence de routes
    SurfMed <- rbind(SurfMed,st_drop_geometry(SurfMed.tmp[,c("id","Nature","Aire_intersect")]))
    
    # actualisation de la progression
    pb_m$tick()
    Sys.sleep(1 / 100)
  }
  
  # gestion de l'extractions des longueurs de Surfs intersectees
  SurfMed <- rbind(SurfMed,Rebus)
  
  # calcul de la proportion des surfaces
  Aire_tot <- BuffMed^2*pi
  SurfMed$Aire_intersect <- SurfMed$Aire_intersect / Aire_tot
  
  
  # conversion long -> wide
  SurfMed.dt <- data.table::data.table(SurfMed)
  SurfMed_wide <- data.table::dcast(SurfMed.dt,id + Aire_intersect ~ Nature, value.var = "Aire_intersect", fill = 0,fun.aggregate = sum) # conversion long -> wide
  
  
  # somme des longueurs de Surfs intersectes par id selon le type
  Extr_Med <-  SurfMed_wide %>%
    group_by(id) %>%
    mutate(SpCESm_1 = if(length(grep("Eau salée permanente",colnames(.))) > 0){sum(`Eau salée permanente`)}else{sum(`NA`)}) %>%
    mutate(SpCESm_2 = if(length(grep("Eau douce permanente",colnames(.))) > 0){sum(`Eau douce permanente`)}else{sum(`NA`)}) %>%
    mutate(SpCESm_3 = if(length(grep("Eau salée non permanente",colnames(.))) > 0){sum(`Eau salée non permanente`)}else{sum(`NA`)}) %>%
    mutate(SpCESm_4 = if(length(grep("Eau douce non permanente",colnames(.))) > 0){sum(`Eau douce non permanente`)}else{sum(`NA`)}) %>%
    mutate(SpCESm_5 = if(length(grep("Névé, glacier",colnames(.))) > 0){sum(`Névé, glacier`)}else{sum(`NA`)})
  
  
  Extr_Med <- unique(Extr_Med[,grep("id|SpCES",colnames(Extr_Med))])
  
  
  
  #############################
  # Buffer locaux : large -----
  
  cat(paste0(c("\t---\tSurfaces d'eau dans un buffer de :",BuffLarg," m","\t---\n")))
  
  vec.iteration_larg <- seq(from = 1, to = length(PointBuffLarg$ID_extract), by = size_batch_large) # argument suplementaire ? / add condition : NA/all -> tous
  
  SurfLarg <- as.data.frame(matrix(nrow=0,ncol = 1))
  colnames(SurfLarg) <- c("id")
  class(SurfLarg$Nature) <- class(SURF$Nature)
  
  Rebus <- data.table::data.table()
  
  # initialisation de la barre de progression
  pb_l <- progress::progress_bar$new(total = length(vec.iteration_larg),
                                     format = "Extraction Aire des surfaces d'eau buffer large [:bar] :percent    Tps ecoule = :elapsedfull",
                                     clear = F)
  
  pb_l$tick(0)
  Sys.sleep(1/20)
  
  
  for(i in vec.iteration_larg){
    Point.tmp <- PointBuffLarg[i:min((i+size_batch_large-1),length(PointBuffLarg$ID_extract)),]
    
    SurfLarg.tmp <- st_intersection(SURF, Point.tmp) # intersection en block
    
    for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
      
      if(nrow(SurfLarg.tmp[SurfLarg.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
        SurfLarg.tmp[SurfLarg.tmp$id == h,"Aire_intersect"] <- as.numeric(st_area(SurfLarg.tmp[SurfLarg.tmp$id == h,]))
        
        
      }else{ # en cas d'absence de route --> formation table de donnees nulles
        
        Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Aire_intersect = 0, Nature = NA)))
      }
      
    }
    
    # jointure des tables de la condition presence de routes
    SurfLarg <- rbind(SurfLarg,st_drop_geometry(SurfLarg.tmp[,c("id","Nature","Aire_intersect")]))
    
    # actualisation de la progression
    pb_l$tick()
    Sys.sleep(1 / 100)
  }
  
  # gestion de l'extractions des longueurs de Surfs intersectees
  SurfLarg <- rbind(SurfLarg,Rebus)
  
  # calcul de la proportion des surfaces
  Aire_tot <- BuffLarg^2*pi
  SurfLarg$Aire_intersect <- SurfLarg$Aire_intersect / Aire_tot
  
  
  # conversion long -> wide
  SurfLarg.dt <- data.table::data.table(SurfLarg)
  SurfLarg_wide <- data.table::dcast(SurfLarg.dt,id + Aire_intersect ~ Nature, value.var = "Aire_intersect", fill = 0,fun.aggregate = sum) # conversion long -> wide
  
  
  # somme des longueurs de Surfs intersectes par id selon le type
  Extr_Larg <-  SurfLarg_wide %>%
    group_by(id) %>%
    mutate(SpCESl_1 = if(length(grep("Eau salée permanente",colnames(.))) > 0){sum(`Eau salée permanente`)}else{sum(`NA`)}) %>%
    mutate(SpCESl_2 = if(length(grep("Eau douce permanente",colnames(.))) > 0){sum(`Eau douce permanente`)}else{sum(`NA`)}) %>%
    mutate(SpCESl_3 = if(length(grep("Eau salée non permanente",colnames(.))) > 0){sum(`Eau salée non permanente`)}else{sum(`NA`)}) %>%
    mutate(SpCESl_4 = if(length(grep("Eau douce non permanente",colnames(.))) > 0){sum(`Eau douce non permanente`)}else{sum(`NA`)}) %>%
    mutate(SpCESl_5 = if(length(grep("Névé, glacier",colnames(.))) > 0){sum(`Névé, glacier`)}else{sum(`NA`)})
  
  
  Extr_Larg <- unique(Extr_Larg[,grep("id|SpCES",colnames(Extr_Larg))])
  
  
  
  # Jointure des tables d'extractions des differents buffers
  SpSURF <- left_join(Extr_Smal,Extr_Med, by="id")
  SpSURF <- left_join(SpSURF,Extr_Larg, by="id")
  
  SpSURF <- dplyr::left_join(SpSURF, Point[,c("id",Coords,"ID_liste")],by="id") # add d'informations sur les listes
  
  
  save(SpSURF, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"EauSURFACES_envEPOC.RData")) # securite
  write.csv(SpSURF, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"EauSURFACES_envEPOC.csv"), row.names = F)
  
  
  cat("\n\n\t-------\tExtract SURFACES d'eau done - check /function_output \t-------\n\n")
  
  
  
  
 
  
  
  
  
  
  
  
  
  
  
   
}