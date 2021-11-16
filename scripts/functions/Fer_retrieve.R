extract_fer = function(dsnTable, names_coord,buffer_small, buffer_medium, buffer_large, dsnFer,
                         prefixe_fichier, size_batch_med, size_batch_large){
  
  
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
  cat("\t-------\tChargement des shapes ROUTES/VOIES FERRES\t-------\n")
  
  FER <- st_read(dsn = dsnFer)
  FER$NATURE <- as.character(FER$NATURE)
  
  
  # Extraction FER ------
  cat("\t-------\tDebut de l'extract des voies ferres \t-------\n\n")
  
  #############################
  # Buffer locaux : small -----
  
  cat(paste0(c("\n\n\t---\tFER dans un buffer de :",BuffSmal," m","\t---\n")))
  
  
  FerSmal <- as.data.frame(matrix(nrow=0,ncol = 1))
  colnames(FerSmal) <- c("id")
  class(FerSmal$NATURE) <- class(FER$NATURE)
  
  Rebus <- data.table::data.table()
  
  
  # pas d'extract en batch : pb de gestion dtf <=> faible densite du reseau ferroviaire
  Point.tmp <- PointBuffSmal[,]
  
  FerSmal.tmp <- st_intersection(FER, Point.tmp) # intersection en block
  
  for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
    
    if(nrow(FerSmal.tmp[FerSmal.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
      FerSmal.tmp[FerSmal.tmp$id == h,"Longueur_intersect"] <- as.numeric(st_length(FerSmal.tmp[FerSmal.tmp$id == h,]))
      
    }else{ # en cas d'absence de Fer --> formation table de donnees nulles
      
      Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Longueur_intersect = 0, NATURE = NA)))
    }
    
  }
  
  # jointure des tables de la condition presence de routes
  FerSmal <- rbind(FerSmal,st_drop_geometry(FerSmal.tmp[,c("id","NATURE","Longueur_intersect")]))
  
  
  # ajout des donnes rebus (0 lignes de fers intersecte)
  FerSmal <- rbind(FerSmal,Rebus)
  
  
  #somme des extracts par id
  ExtrF_Smal <- FerSmal %>%
    group_by(id) %>%
    mutate(SpFes = sum(Longueur_intersect))
  
  ExtrF_Smal <- unique(ExtrF_Smal[,grep("id|SpFe",colnames(ExtrF_Smal))])
  
  
  
  #############################
  # Buffer locaux : medium -----
  
  cat(paste0(c("\n\n\t---\tFER dans un buffer de :",BuffMed," m","\t---\n")))
  
  vec.iteration_med <- seq(from = 1, to = length(PointBuffMed$ID_extract), by = size_batch_med) # argument suplementaire ? / add condition : NA/all -> tous
  
  FerMed <- as.data.frame(matrix(nrow=0,ncol = 1))
  colnames(FerMed) <- c("id")
  class(FerMed$NATURE) <- class(FER$NATURE)
  
  Rebus <- data.table::data.table()
  
  # initialisation de la barre de progression
  pb_m <- progress::progress_bar$new(total = length(vec.iteration_med),
                                     format = "Extraction longueur des voies ferrees buffer medium [:bar] :percent    Tps ecoule = :elapsedfull",
                                     clear = F)
  
  
  pb_m$tick(0)
  Sys.sleep(1/20)
  
  
  for(i in vec.iteration_med){
    Point.tmp <- PointBuffMed[i:min((i+size_batch_med-1),length(PointBuffMed$ID_extract)),]
    
    FerMed.tmp <- st_intersection(FER, Point.tmp) # intersection en block
    
    for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
      
      if(nrow(FerMed.tmp[FerMed.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
        FerMed.tmp[FerMed.tmp$id == h,"Longueur_intersect"] <- as.numeric(st_length(FerMed.tmp[FerMed.tmp$id == h,]))
        
      }else{ # en cas d'absence de Fer --> formation table de donnees nulles
        
        Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Longueur_intersect = 0, NATURE = NA)))
      }
      
    }
    
    # jointure des tables de la condition presence de routes
    FerMed <- rbind(FerMed,st_drop_geometry(FerMed.tmp[,c("id","NATURE","Longueur_intersect")]))
    
    # actualisation de la progression
    pb_m$tick()
    Sys.sleep(1 / 100)
  }
  
  # ajout des donnes rebus (0 lignes de fers intersecte)
  FerMed <- rbind(FerMed,Rebus)
  
  
  #somme des extracts par id
  ExtrF_Med <- FerMed %>%
    group_by(id) %>%
    mutate(SpFem = sum(Longueur_intersect))
  
  ExtrF_Med <- unique(ExtrF_Med[,grep("id|SpFe",colnames(ExtrF_Med))])
  
  
  
  #############################
  # Buffer locaux : large -----
  
  cat(paste0(c("\n\n\t---\tFER dans un buffer de :",BuffLarg," m","\t---\n")))
  
  vec.iteration_larg <- seq(from = 1, to = length(PointBuffLarg$ID_extract), by = size_batch_large) # argument suplementaire ? / add condition : NA/all -> tous
  
  FerLarg <- as.data.frame(matrix(nrow=0,ncol = 1))
  colnames(FerLarg) <- c("id")
  class(FerLarg$NATURE) <- class(FER$NATURE)
  
  Rebus <- data.table::data.table()
  
  # initialisation de la barre de progression
  pb_l <- progress::progress_bar$new(total = length(vec.iteration_larg),
                                     format = "Extraction longueur des voies ferrees buffer large [:bar] :percent    Tps ecoule = :elapsedfull",
                                     clear = F)
  
  
  pb_l$tick(0)
  Sys.sleep(1/20)
  
  
  for(i in vec.iteration_larg){
    Point.tmp <- PointBuffLarg[i:min((i+size_batch_large-1),length(PointBuffLarg$ID_extract)),]
    
    FerLarg.tmp <- st_intersection(FER, Point.tmp) # intersection en block
    
    for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
      
      if(nrow(FerLarg.tmp[FerLarg.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
        FerLarg.tmp[FerLarg.tmp$id == h,"Longueur_intersect"] <- as.numeric(st_length(FerLarg.tmp[FerLarg.tmp$id == h,]))
        
      }else{ # en cas d'absence de Fer --> formation table de donnees nulles
        
        Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Longueur_intersect = 0, NATURE = NA)))
      }
      
    }
    
    # jointure des tables de la condition presence de routes
    FerLarg <- rbind(FerLarg,st_drop_geometry(FerLarg.tmp[,c("id","NATURE","Longueur_intersect")]))
    
    # actualisation de la progression
    pb_l$tick()
    Sys.sleep(1 / 100)
  }
  
  # ajout des donnes rebus (0 lignes de fers intersecte)
  FerLarg <- rbind(FerLarg,Rebus)
  
  
  #somme des extracts par id
  ExtrF_Larg <- FerLarg %>%
    group_by(id) %>%
    mutate(SpFel = sum(Longueur_intersect))
  
  ExtrF_Larg <- unique(ExtrF_Larg[,grep("id|SpFe",colnames(ExtrF_Larg))])
  
  
  # Jointure des tables d'extractions des differents buffers
  SpFER <- left_join(ExtrF_Smal,ExtrF_Med, by="id")
  SpFER <- left_join(SpFER,ExtrF_Larg, by="id")
  
  SpFER <- dplyr::left_join(SpFER, Point[,c("id",Coords,"ID_liste")]) # add d'informations sur les listes
  
  
  save(SpFER, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"FER_envEPOC.RData")) # securite
  write.csv(SpFER, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"FER_envEPOC.csv"), row.names = F)
  
  
  cat("\n\n\t-------\tExtract FER done - check /function_output \t-------\n\n")
  
  
  
  
  
  
  
  
  
  
  
}