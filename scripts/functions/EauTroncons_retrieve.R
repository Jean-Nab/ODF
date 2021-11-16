
# Fonctions basees selon le modele des fonctions d'extraction de Yves Bas (https://github.com/cesco-lab/Vigie-Chiro_scripts/tree/master/functions/extractGI) 
## Sources couches shp : BD CARTHAGE 2017 tronçons & surfaces hydographique 




## dsnTable = "C:/git/epoc/DS.v2/epoc_barycentre_liste_density_add.csv"
## names_coord = c("X_barycentre_L93","Y_barycentre_L93")
## buffer_small = 50
## buffer_medium = 500
## buffer_large = 5000
## dsnTronc = "C:/Users/Travail/Desktop/Ressource QGis/france/hydrographie/TronconHydrograElt_FXX.shp"
## dsnSurf = "C:/Users/Travail/Desktop/Ressource QGis/france/hydrographie/EltHydroSurface_FXX.shp"



extract_eauTroncons = function(dsnTable, names_coord,buffer_small, buffer_medium, buffer_large, dsnTronc,
                               prefixe_fichier){
  
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
  cat("\t-------\tChargement des shapes TRONCONS d'eau\t-------\n")
  
  TRONC <- st_read(dsn = dsnTronc)
  TRONC$Etat <- as.character(TRONC$Etat)

  
  # exttraction longueur des troncons -----
  
  #############################
  # Buffer locaux : small -----
  
  cat(paste0(c("\t---\tTroncons d'eau dans un buffer de :",BuffSmal," m","\t---\n")))
  
  vec.iteration_smal <- seq(from = 1, to = length(PointBuffSmal$ID_extract), by = 305) # argument suplementaire ? / add condition : NA/all -> tous
  
  TroncSmal <- as.data.frame(matrix(nrow=0,ncol = 1))
  colnames(TroncSmal) <- c("id")
  class(TroncSmal$Etat) <- class(TRONC$Etat)
  
  Rebus <- data.table::data.table()
  
  # initialisation de la barre de progression
  pb_s <- progress::progress_bar$new(total = length(vec.iteration_smal),
                                     format = "Extraction longueur des Troncons d'eau buffer small [:bar] :percent    Tps ecoule = :elapsedfull",
                                     clear = F)
  
  pb_s$tick(0)
  Sys.sleep(1/20)
  
  
  for(i in vec.iteration_smal){
    Point.tmp <- PointBuffSmal[i:min((i+304),length(PointBuffSmal$ID_extract)),]
    
    TroncSmal.tmp <- st_intersection(TRONC, Point.tmp) # intersection en block
    
    for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
      
      if(nrow(TroncSmal.tmp[TroncSmal.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
        TroncSmal.tmp[TroncSmal.tmp$id == h,"Longueur_intersect"] <- as.numeric(st_length(TroncSmal.tmp[TroncSmal.tmp$id == h,]))
        
      }else{ # en cas d'absence de route --> formation table de donnees nulles
        
        Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Longueur_intersect = 0, Etat = NA)))
      }
      
    }
    
    # jointure des tables de la condition presence de routes
    TroncSmal <- rbind(TroncSmal,st_drop_geometry(TroncSmal.tmp[,c("id","Etat","Longueur_intersect")]))
    
    # actualisation de la progression
    pb_s$tick()
    Sys.sleep(1 / 100)
  }
  
  # gestion de l'extractions des longueurs de Troncs intersectees
  TroncSmal <- rbind(TroncSmal,Rebus)
  
  # conversion long -> wide
  TroncSmal.dt <- data.table::data.table(TroncSmal)
  TroncSmal_wide <- data.table::dcast(TroncSmal.dt,id + Longueur_intersect ~ Etat, value.var = "Longueur_intersect", fill = 0,fun.aggregate = sum) # conversion long -> wide
  
  
  # somme des longueurs de Troncs intersectes par id selon le type
  Extr_Smal <-  TroncSmal_wide %>%
    group_by(id) %>%
    mutate(SpCELs_1 = if(length(grep("Permanent",colnames(.))) > 0){sum(`Permanent`)}else{sum(`NA`)}) %>%
    mutate(SpCELs_2 = if(length(grep("Fictif",colnames(.))) > 0){sum(`Fictif`)}else{sum(`NA`)}) %>%
    mutate(SpCELs_3 = if(length(grep("Intermittent",colnames(.))) > 0){sum(`Intermittent`)}else{sum(`NA`)}) %>%
    mutate(SpCELs_4 = if(length(grep("Inconnu",colnames(.))) > 0){sum(`Inconnu`)}else{sum(`NA`)}) %>%
    mutate(SpCELs_5 = if(length(grep("En attente de mise à jour",colnames(.))) > 0){sum(`En attente de mise à jour`)}else{sum(`NA`)}) %>%
    mutate(SpCELs_6 = if(length(grep("A sec",colnames(.))) > 0){sum(`A sec`)}else{sum(`NA`)})
  
  
  Extr_Smal <- unique(Extr_Smal[,grep("id|SpCEL",colnames(Extr_Smal))])
  
  
  
  #############################
  # Buffer locaux : medium -----
  
  cat(paste0(c("\t---\tTroncons d'eau dans un buffer de :",BuffMed," m","\t---\n")))
  
  vec.iteration_med <- seq(from = 1, to = length(PointBuffMed$ID_extract), by = 305) # argument suplementaire ? / add condition : NA/all -> tous
  
  TroncMed <- as.data.frame(matrix(nrow=0,ncol = 1))
  colnames(TroncMed) <- c("id")
  class(TroncMed$Etat) <- class(TRONC$Etat)
  
  Rebus <- data.table::data.table()
  
  # initialisation de la barre de progression
  pb_m <- progress::progress_bar$new(total = length(vec.iteration_med),
                                     format = "Extraction longueur des Troncons d'eau buffer medium [:bar] :percent    Tps ecoule = :elapsedfull",
                                     clear = F)
  
  pb_m$tick(0)
  Sys.sleep(1/20)
  
  
  for(i in vec.iteration_med){
    Point.tmp <- PointBuffMed[i:min((i+304),length(PointBuffMed$ID_extract)),]
    
    TroncMed.tmp <- st_intersection(TRONC, Point.tmp) # intersection en block
    
    for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
      
      if(nrow(TroncMed.tmp[TroncMed.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
        TroncMed.tmp[TroncMed.tmp$id == h,"Longueur_intersect"] <- as.numeric(st_length(TroncMed.tmp[TroncMed.tmp$id == h,]))
        
      }else{ # en cas d'absence de route --> formation table de donnees nulles
        
        Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Longueur_intersect = 0, Etat = NA)))
      }
      
    }
    
    # jointure des tables de la condition presence de routes
    TroncMed <- rbind(TroncMed,st_drop_geometry(TroncMed.tmp[,c("id","Etat","Longueur_intersect")]))
    
    # actualisation de la progression
    pb_m$tick()
    Sys.sleep(1 / 100)
  }
  
  # gestion de l'extractions des longueurs de Troncs intersectees
  TroncMed <- rbind(TroncMed,Rebus)
  
  # conversion long -> wide
  TroncMed.dt <- data.table::data.table(TroncMed)
  TroncMed_wide <- data.table::dcast(TroncMed.dt,id + Longueur_intersect ~ Etat, value.var = "Longueur_intersect", fill = 0,fun.aggregate = sum) # conversion long -> wide
  
  
  # somme des longueurs de Troncs intersectes par id selon le type
  Extr_Med <-  TroncMed_wide %>%
    group_by(id) %>%
    mutate(SpCELm_1 = if(length(grep("Permanent",colnames(.))) > 0){sum(`Permanent`)}else{sum(`NA`)}) %>%
    mutate(SpCELm_2 = if(length(grep("Fictif",colnames(.))) > 0){sum(`Fictif`)}else{sum(`NA`)}) %>%
    mutate(SpCELm_3 = if(length(grep("Intermittent",colnames(.))) > 0){sum(`Intermittent`)}else{sum(`NA`)}) %>%
    mutate(SpCELm_4 = if(length(grep("Inconnu",colnames(.))) > 0){sum(`Inconnu`)}else{sum(`NA`)}) %>%
    mutate(SpCELm_5 = if(length(grep("En attente de mise à jour",colnames(.))) > 0){sum(`En attente de mise à jour`)}else{sum(`NA`)}) %>%
    mutate(SpCELm_6 = if(length(grep("A sec",colnames(.))) > 0){sum(`A sec`)}else{sum(`NA`)})
  
  
  Extr_Med <- unique(Extr_Med[,grep("id|SpCEL",colnames(Extr_Med))])
  
  
  
  #############################
  # Buffer locaux : large -----
  
  cat(paste0(c("\t---\tTroncons d'eau dans un buffer de :",BuffLarg," m","\t---\n")))
  
  vec.iteration_larg <- seq(from = 1, to = length(PointBuffLarg$ID_extract), by = 305) # argument suplementaire ? / add condition : NA/all -> tous
  
  TroncLarg <- as.data.frame(matrix(nrow=0,ncol = 1))
  colnames(TroncLarg) <- c("id")
  class(TroncLarg$Etat) <- class(TRONC$Etat)
  
  Rebus <- data.table::data.table()
  
  # initialisation de la barre de progression
  pb_l <- progress::progress_bar$new(total = length(vec.iteration_larg),
                                     format = "Extraction longueur des Troncons d'eau buffer large [:bar] :percent    Tps ecoule = :elapsedfull",
                                     clear = F)
  
  pb_l$tick(0)
  Sys.sleep(1/20)
  
  
  for(i in vec.iteration_larg){
    Point.tmp <- PointBuffLarg[i:min((i+304),length(PointBuffLarg$ID_extract)),]
    
    TroncLarg.tmp <- st_intersection(TRONC, Point.tmp) # intersection en block
    
    for(h in Point.tmp$id){ # boucle afin d'appliquer la condition de verification (presence ou absence de Routes ?)
      
      if(nrow(TroncLarg.tmp[TroncLarg.tmp$id == h,]) > 0){ # Si presence --> calcul longueur
        TroncLarg.tmp[TroncLarg.tmp$id == h,"Longueur_intersect"] <- as.numeric(st_length(TroncLarg.tmp[TroncLarg.tmp$id == h,]))
        
      }else{ # en cas d'absence de route --> formation table de donnees nulles
        
        Rebus <- data.table::rbindlist(list(Rebus,data.table::data.table(id = h, Longueur_intersect = 0, Etat = NA)))
      }
      
    }
    
    # jointure des tables de la condition presence de routes
    TroncLarg <- rbind(TroncLarg,st_drop_geometry(TroncLarg.tmp[,c("id","Etat","Longueur_intersect")]))

    # actualisation de la progression
    pb_l$tick()
    Sys.sleep(1 / 100)
  }
  
  # gestion de l'extractions des longueurs de Troncs intersectees
  TroncLarg <- rbind(TroncLarg,Rebus)
  
  # conversion long -> wide
  TroncLarg.dt <- data.table::data.table(TroncLarg)
  TroncLarg_wide <- data.table::dcast(TroncLarg.dt,id + Longueur_intersect ~ Etat, value.var = "Longueur_intersect", fill = 0,fun.aggregate = sum) # conversion long -> wide
  
  
  # somme des longueurs de Troncs intersectes par id selon le type
  Extr_Larg <-  TroncLarg_wide %>%
    group_by(id) %>%
    mutate(SpCELl_1 = if(length(grep("Permanent",colnames(.))) > 0){sum(`Permanent`)}else{sum(`NA`)}) %>%
    mutate(SpCELl_2 = if(length(grep("Fictif",colnames(.))) > 0){sum(`Fictif`)}else{sum(`NA`)}) %>%
    mutate(SpCELl_3 = if(length(grep("Intermittent",colnames(.))) > 0){sum(`Intermittent`)}else{sum(`NA`)}) %>%
    mutate(SpCELl_4 = if(length(grep("Inconnu",colnames(.))) > 0){sum(`Inconnu`)}else{sum(`NA`)}) %>%
    mutate(SpCELl_5 = if(length(grep("En attente de mise à jour",colnames(.))) > 0){sum(`En attente de mise à jour`)}else{sum(`NA`)}) %>%
    mutate(SpCELl_6 = if(length(grep("A sec",colnames(.))) > 0){sum(`A sec`)}else{sum(`NA`)})
  
  
  Extr_Larg <- unique(Extr_Larg[,grep("id|SpCEL",colnames(Extr_Larg))])
  
  
  # Jointure des tables d'extractions des differents buffers
  SpTRONC <- left_join(Extr_Smal,Extr_Med, by="id")
  SpTRONC <- left_join(SpTRONC,Extr_Larg, by="id")
  
  SpTRONC <- dplyr::left_join(SpTRONC, Point[,c("id",Coords,"ID_liste")],by="id") # add d'informations sur les listes
  
  
  save(SpTRONC, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"TRONCONS_envEPOC.RData")) # securite
  write.csv(SpTRONC, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"TRONCONS_envEPOC.csv"), row.names = F)
  
  cat("\n\n\t-------\tExtract TRONCONS done - check /function_output \t-------\n\n")
  

}























