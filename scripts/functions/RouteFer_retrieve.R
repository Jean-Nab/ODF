
# Fonctions basees selon le modele des fonctions d'extraction de Yves Bas (https://github.com/cesco-lab/Vigie-Chiro_scripts/tree/master/functions/extractGI) 
## Sources couches shp : BD ROUTE 500 Edition 201 / https://geoservices.ign.fr/documentation/diffusion/telechargement-donnees-libres.html#route-500


### Warning : extraction route large  differents : gestion du dtf d'extract dans la loop --> vs big dtf


##  # telechargement du shape route
##    library(curl)
##    url = "ftp://ROUTE_500_ext:UqueemievaiDah3k@ftp3.ign.fr/ROUTE500_3-0__SHP_LAMB93_FXX_2020-08-04.7z.001"
##    h = new_handle(dirlistonly=TRUE)
##    con = curl(url, "r", h)
##    close(con)
##  
##    curl_download(url,destfile = "C://Users/Travail/Desktop/TEST.7z") # fichier en .7z

## dsnTable = "C:/git/test_epoc1500.csv"
## names_coord = c("X_barycentre_L93","Y_barycentre_L93")
## buffer_small = 50
## buffer_medium = 500
## buffer_large = 5000
## dsnRoute = "C:/Users/Travail/Desktop/Ressource QGis/france/ROUTE500_3-0__SHP_LAMB93_FXX_2020-08-04/ROUTE500/1_DONNEES_LIVRAISON_2020-08-00223/R500_3-0_SHP_LAMB93_FXX-ED201/RESEAU_ROUTIER/TRONCON_ROUTE.shp"
## dsnFer = "C:/Users/Travail/Desktop/Ressource QGis/france/ROUTE500_3-0__SHP_LAMB93_FXX_2020-08-04/ROUTE500/1_DONNEES_LIVRAISON_2020-08-00223/R500_3-0_SHP_LAMB93_FXX-ED201/RESEAU_FERRE/TRONCON_VOIE_FERREE.shp"

# extract_route(dsnTable,names_coord,buffer_small,buffer_medium,buffer_large,dsnRoute,dsnFer)


extract_route = function(dsnTable, names_coord,buffer_small, buffer_medium, buffer_large, dsnRoute, dsnFer){
  
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
  cat("\t-------\tChargement des shapes ROUTES/VOIES FERRES\t-------\n")
  
  ROUTE <- st_read(dsn = dsnRoute)
  ROUTE$VOCATION <- as.character(ROUTE$VOCATION)
  
  FER <- st_read(dsn = dsnFer)
  FER$NATURE <- as.character(FER$NATURE)
  
  
  # Extraction ROUTE ------
  
  #############################
  # Buffer locaux : small -----
  
    cat(paste0(c("\n\n\t---\tROUTE dans un buffer de :",BuffSmal," m","\t---\n")))
  
    vec.iteration_smal <- seq(from = 1, to = length(PointBuffSmal$ID_extract), by = 305) # argument suplementaire ? / add condition : NA/all -> tous
    
    RouteSmal <- as.data.frame(matrix(nrow=0,ncol = 1))
    colnames(RouteSmal) <- c("id")
    class(RouteSmal$VOCATION) <- class(ROUTE$VOCATION)
    
    Rebus <- data.table::data.table()
    
    # initialisation de la barre de progression
    pb_s <- progress::progress_bar$new(total = length(vec.iteration_smal),
                                       format = "Extraction longueur des routes buffer small [:bar] :percent    Tps ecoule = :elapsedfull",
                                       clear = F)
    
    
    pb_s$tick(0)
    Sys.sleep(1/20)
    
    
    for(i in vec.iteration_smal){
      Point.tmp <- PointBuffSmal[i:min((i+304),length(PointBuffSmal$ID_extract)),]
      
      RouteSmal.tmp <- st_intersection(ROUTE, Point.tmp) # intersection en block
      
      for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
        
        if(nrow(RouteSmal.tmp[RouteSmal.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
          RouteSmal.tmp[RouteSmal.tmp$id == h,"Longueur_intersect"] <- as.numeric(st_length(RouteSmal.tmp[RouteSmal.tmp$id == h,]))
          
        }else{ # en cas d'absence de route --> formation table de donnees nulles
          
          Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Longueur_intersect = 0, VOCATION = NA)))
        }
        
      }
      
      # jointure des tables de la condition presence de routes
      RouteSmal <- rbind(RouteSmal,st_drop_geometry(RouteSmal.tmp[,c("id","VOCATION","Longueur_intersect")]))
      
      # actualisation de la progression
      pb_s$tick()
      Sys.sleep(1 / 100)
    }

    
    # gestion de l'extract
      # jointure des tables forme pendant l'extract selon la condition if-else
      RouteSmal <- rbind(RouteSmal,Rebus)
      
    RouteSmal_wide <- data.table::data.table(RouteSmal)
    RouteSmal_wide <- data.table::dcast(RouteSmal_wide,id + Longueur_intersect ~ VOCATION, value.var = "Longueur_intersect", fill = 0) # conversion long -> wide
    RouteSmal_wide$`NA` <- 0
    
    
    #somme des extracts par id selon le type de routes
    Extr_Smal <-  RouteSmal_wide %>%
      group_by(id) %>%
      mutate(SpRoLLs = if(length(grep("locale",colnames(.))) > 0){sum(`Liaison locale`)}else{sum(`NA`)}) %>%
      mutate(SpRoLRs = if(length(grep("régionale",colnames(.))) > 0){sum(`Liaison régionale`)}else{sum(`NA`)}) %>%
      mutate(SpRoLPs = if(length(grep("principale",colnames(.))) > 0){sum(`Liaison principale`)}else{sum(`NA`)}) %>%
      mutate(SpRoBs = if(length(grep("Bretelle",colnames(.))) > 0){sum(`Bretelle`)}else{sum(`NA`)}) %>%
      mutate(SpRoTAPs = if(length(grep("autoroutier",colnames(.))) > 0){sum(`Type autoroutier`)}else{sum(`NA`)})
    
    
    Extr_Smal <- unique(Extr_Smal[,grep("id|SpRo",colnames(Extr_Smal))])
    
  
    
  ##############################
  # Buffer locaux : medium -----
    
    cat(paste0(c("\t---\tROUTE dans un buffer de :",BuffMed," m","\t---\n")))
    
    vec.iteration_med <- seq(from = 1, to = length(PointBuffMed$ID_extract), by = 305) # argument suplementaire ? / add condition : NA/all -> tous
    
    RouteMed <- as.data.frame(matrix(nrow=0,ncol = 1))
    colnames(RouteMed) <- c("id")
    class(RouteMed$VOCATION) <- class(ROUTE$VOCATION)
    
    Rebus <- data.table::data.table()
    
    # initialisation de la barre de progression
    pb_m <- progress::progress_bar$new(total = length(vec.iteration_med),
                                       format = "Extraction longueur des routes buffer medium [:bar] :percent    Tps ecoule = :elapsedfull",
                                       clear = F)
    
    
    pb_m$tick(0)
    Sys.sleep(1/20)
    
    for(i in vec.iteration_med){
      Point.tmp <- PointBuffMed[i:min((i+304),length(PointBuffMed$ID_extract)),]
      
      RouteMed.tmp <- st_intersection(ROUTE, Point.tmp) # intersection en block
      
      for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
        
        if(nrow(RouteMed.tmp[RouteMed.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
          RouteMed.tmp[RouteMed.tmp$id == h,"Longueur_intersect"] <- as.numeric(st_length(RouteMed.tmp[RouteMed.tmp$id == h,]))
          
        }else{ # en cas d'absence de route --> formation table de donnees nulles
          
          Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Longueur_intersect = 0, VOCATION = NA)))
        }
        
      }
      
      # jointure des tables de la condition presence de routes
      RouteMed <- rbind(RouteMed,st_drop_geometry(RouteMed.tmp[,c("id","VOCATION","Longueur_intersect")]))
      
      # actualisation de la progression
      pb_m$tick()
      Sys.sleep(1 / 100)
    }
    
    # gestion de l'extract en parallele
      # jointure des tables forme pendant l'extract selon la condition if-else
      RouteMed <- rbind(RouteMed,Rebus)
    
    RouteMed_wide <- data.table::data.table(RouteMed)
    RouteMed_wide <- data.table::dcast(RouteMed_wide,id + Longueur_intersect ~ VOCATION, value.var = "Longueur_intersect", fill = 0) # conversion long -> wide
    RouteMed_wide$`NA` <- 0
    
    
    #somme des extracts par id selon le type de routes
    Extr_Med <-  RouteMed_wide %>%
      group_by(id) %>%
      mutate(SpRoLLm = if(length(grep("locale",colnames(.))) > 0){sum(`Liaison locale`)}else{sum(`NA`)}) %>%
      mutate(SpRoLRm = if(length(grep("régionale",colnames(.))) > 0){sum(`Liaison régionale`)}else{sum(`NA`)}) %>%
      mutate(SpRoLPm = if(length(grep("principale",colnames(.))) > 0){sum(`Liaison principale`)}else{sum(`NA`)}) %>%
      mutate(SpRoBm = if(length(grep("Bretelle",colnames(.))) > 0){sum(`Bretelle`)}else{sum(`NA`)}) %>%
      mutate(SpRoTAPm = if(length(grep("autoroutier",colnames(.))) > 0){sum(`Type autoroutier`)}else{sum(`NA`)})
    
    
    Extr_Med <- unique(Extr_Med[,grep("id|SpRo",colnames(Extr_Med))])
  
    
  #############################
  # Buffer locaux : large -----
    
    cat(paste0(c("\t---\tROUTE dans un buffer de :",BuffLarg," m","\t---\n")))
    
    vec.iteration_larg <- seq(from = 1, to = length(PointBuffLarg$ID_extract), by = 150) # argument suplementaire ? / add condition : NA/all -> tous
    
    RouteLarg <- as.data.frame(matrix(nrow=0,ncol = 1))
    colnames(RouteLarg) <- c("id")
    class(RouteLarg$VOCATION) <- class(ROUTE$VOCATION)
    
    Rebus <- data.table::data.table()
    
    # initialisation de la barre de progression
    pb_l <- progress::progress_bar$new(total = length(vec.iteration_larg),
                                       format = "Extraction longueur des routes buffer large [:bar] :percent    Tps ecoule = :elapsedfull",
                                       clear = F)
    
    
    pb_l$tick(0)
    Sys.sleep(1/20)
    
    
    for(i in vec.iteration_larg){
      Point.tmp <- PointBuffLarg[i:min((i+149),length(PointBuffLarg$ID_extract)),]
      
      RouteLarg.tmp <- st_intersection(ROUTE, Point.tmp) # intersection en block
      
      for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
        
        if(nrow(RouteLarg.tmp[RouteLarg.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
          RouteLarg.tmp[RouteLarg.tmp$id == h,"Longueur_intersect"] <- as.numeric(st_length(RouteLarg.tmp[RouteLarg.tmp$id == h,]))
          
        }else{ # en cas d'absence de route --> formation table de donnees nulles
          
          Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Longueur_intersect = 0, VOCATION = NA)))
        }
        
      }
      
      # jointure des tables de la condition presence de routes
     RouteLarg <- rbind(RouteLarg,st_drop_geometry(RouteLarg.tmp[,c("id","VOCATION","Longueur_intersect")]))
      
      # actualisation de la progression
      pb_l$tick()
      Sys.sleep(1 / 100)
    }
    
    
      
    # gestion de l'extract en parallele
      # jointure des tables forme pendant l'extract selon la condition if-else
    RouteLarg <- rbind(RouteLarg,Rebus)
    
    RouteLarg_wide <- data.table::data.table(RouteLarg)
    RouteLarg_wide <- data.table::dcast(RouteLarg_wide,id + Longueur_intersect ~ VOCATION, value.var = "Longueur_intersect", fill = 0) # conversion long -> wide
    RouteLarg_wide$`NA` <- 0
    
    
    #somme des extracts par id selon le type de routes
    Extr_Larg <-  RouteLarg_wide %>%
      group_by(id) %>%
      mutate(SpRoLLl = if(length(grep("locale",colnames(.))) > 0){sum(`Liaison locale`)}else{sum(`NA`)}) %>%
      mutate(SpRoLRl = if(length(grep("régionale",colnames(.))) > 0){sum(`Liaison régionale`)}else{sum(`NA`)}) %>%
      mutate(SpRoLPl = if(length(grep("principale",colnames(.))) > 0){sum(`Liaison principale`)}else{sum(`NA`)}) %>%
      mutate(SpRoBl = if(length(grep("Bretelle",colnames(.))) > 0){sum(`Bretelle`)}else{sum(`NA`)}) %>%
      mutate(SpRoTAPl = if(length(grep("autoroutier",colnames(.))) > 0){sum(`Type autoroutier`)}else{sum(`NA`)})
    
    
    Extr_Larg <- unique(Extr_Larg[,grep("id|SpRo",colnames(Extr_Larg))])
      
    
    
  # Jointure des tables d'extractions des differents buffers
    SpROUTE <- left_join(Extr_Smal,Extr_Med, by="id")
    SpROUTE <- left_join(SpROUTE,Extr_Larg, by="id")
    
    SpROUTE <- dplyr::left_join(SpROUTE, Point[,c("id",Coords,"ID_liste")]) # add d'informations sur les listes
    
    
    save(SpROUTE, file = "C:/git/ODF/output/function_output/Rimage_ROUTE_envEPOC.RData") # securite
    write.csv(SpROUTE, file = "C:/git/ODF/output/function_output/ROUTE_envEPOC.csv", row.names = F)
    
    
    cat("\n\n\t-------\tExtract ROUTE done - check /function_output \t-------\n\n")
  
    
    
    
  # Extraction FER ------
    cat("\t-------\tDebut de l'extract des voies ferres \t-------\n\n")
    
    #############################
    # Buffer locaux : small -----
    
    cat(paste0(c("\n\n\t---\tFER dans un buffer de :",BuffSmal," m","\t---\n")))
    
    vec.iteration_smal <- seq(from = 1, to = length(PointBuffSmal$ID_extract), by = 305) # argument suplementaire ? / add condition : NA/all -> tous
    
    FerSmal <- as.data.frame(matrix(nrow=0,ncol = 1))
    colnames(FerSmal) <- c("id")
    class(FerSmal$NATURE) <- class(FER$NATURE)
    
    Rebus <- data.table::data.table()
    
    # initialisation de la barre de progression
    pb_s <- progress::progress_bar$new(total = length(vec.iteration_smal),
                                       format = "Extraction longueur des voies ferrees buffer small [:bar] :percent    Tps ecoule = :elapsedfull",
                                       clear = F)
    
    
    pb_s$tick(0)
    Sys.sleep(1/20)
    
    
    for(i in vec.iteration_smal){
      Point.tmp <- PointBuffSmal[i:min((i+304),length(PointBuffSmal$ID_extract)),]
      
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
      
      # actualisation de la progression
      pb_s$tick()
      Sys.sleep(1 / 100)
    }
    
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
    
    vec.iteration_med <- seq(from = 1, to = length(PointBuffMed$ID_extract), by = 305) # argument suplementaire ? / add condition : NA/all -> tous
    
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
      Point.tmp <- PointBuffMed[i:min((i+304),length(PointBuffMed$ID_extract)),]
      
      FerMed.tmp <- st_intersection(FER, Point.tmp) # intersection en block
      
      for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
        
        if(nrow(FerMed.tmp[FerMed.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
          FerMed.tmp[FerMed.tmp$id == h,"Longueur_intersect"] <- as.numeric(st_length(FerMed.tmp[FerMed.tmp$id == h,]))
          
        }else{ # en cas d'absence de Fer --> formation table de donnees nulles
          
          Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Longueur_intersect = 0, NATURE = NA)))
        }
        
      }
      
      # jointure des tables de la condition presence de routes
      FerMed.tmp <- st_drop_geometry(FerMed.tmp)
      
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
    
    vec.iteration_larg <- seq(from = 1, to = length(PointBuffLarg$ID_extract), by = 150) # argument suplementaire ? / add condition : NA/all -> tous
    
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
      Point.tmp <- PointBuffLarg[i:min((i+149),length(PointBuffLarg$ID_extract)),]
      
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
    
    
    save(SpFER, file = "C:/git/ODF/output/function_output/Rimage_FER_envEPOC.RData") # securite
    write.csv(SpFER, file = "C:/git/ODF/output/function_output/FER_envEPOC.csv", row.names = F)
    
    
    cat("\n\n\t-------\tExtract FER done - check /function_output \t-------\n\n")
    
    
    
}
























