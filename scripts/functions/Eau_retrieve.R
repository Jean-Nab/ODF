
# Fonctions basees selon le modele des fonctions d'extraction de Yves Bas (https://github.com/cesco-lab/Vigie-Chiro_scripts/tree/master/functions/extractGI) 
## Sources couches shp : BD CARTHAGE 2017 tronçons & surfaces hydographique 




dsnTable = "C:/git/epoc/DS.v2/epoc_barycentre_liste_density_add.csv"
names_coord = c("X_barycentre_L93","Y_barycentre_L93")
buffer_small = 50
buffer_medium = 500
buffer_large = 5000
dsnTronc = "C:/Users/Travail/Desktop/Ressource QGis/france/hydrographie/TronconHydrograElt_FXX.shp"
dsnSurf = "C:/Users/Travail/Desktop/Ressource QGis/france/hydrographie/EltHydroSurface_FXX.shp"



extract_eau = function(dsnTable, names_coord,buffer_small, buffer_medium, buffer_large, dsnTronc, dsnSurf){
  
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
  cat("\t-------\tChargement des shapes TRONCONS/SURFACES d'eau\t-------\n")
  
  TRONC <- st_read(dsn = dsnTronc)
  TRONC$Etat <- as.character(TRONC$Etat)
  
  SURF <- st_read(dsn = dsnSurf)
  SURF$Nature <- as.character(SURF$Nature) 
  
  
  # exttraction longueur des troncons -----
  
    #############################
    # Buffer locaux : small -----
    
    cat(paste0(c("\t---\tTroncons d'eau dans un buffer de :",BuffSmal," m","\t---\n")))
    
    Extr_Smal <- as.data.frame(matrix(nrow=0,ncol = 1))
    colnames(Extr_Smal) <- c("id")
    
    # initialisation de la barre de progression
    pb_s <- progress::progress_bar$new(total = nrow(PointbuffSmal),
                                       format = "Extraction longueur des Troncons d'eau buffer small [:bar] :percent    Tps ecoule = :elapsedfull",
                                       clear = F)
    
    
    pb_s$tick(0)
    Sys.sleep(1/20)
    
    for(i in 1:nrow(PointbuffSmal)){ # nrow(PointbuffSmal)
      
      Point.tmp <- PointBuffSmal[i,c("ID_liste","id","Maille","Jour_de_l_annee")]
      
      TroncSmal.tmp <- st_intersection(TRONC, Point.tmp)
      
      if(nrow(TroncSmal.tmp) > 0){
        TroncSmal.tmp$Longueur_intersect <- as.numeric(st_length(TroncSmal.tmp))
        
      }else{ # en cas d'absence de Tronc --> recuperation des infos  
        
        TroncSmal.tmp[1,"id"] <- Point.tmp$id
        TroncSmal.tmp$Longueur_intersect <- 0
        TroncSmal.tmp$Etat <- NA
        
      }
      
      
      # gestion de l'extractions des longueurs de Troncs intersectees
      TroncSmal.tmp <- st_drop_geometry(TroncSmal.tmp[,c("id","Etat","Longueur_intersect")])
      
      # conversion long -> wide
      TroncSmal.tmp.dt <- data.table::data.table(TroncSmal.tmp)
      TroncSmal_wide <- data.table::dcast(TroncSmal.tmp.dt,id + Longueur_intersect ~ Etat, value.var = "Longueur_intersect", fill = 0) # conversion long -> wide
      TroncSmal_wide$`NA` <- 0 # rajout d'une colonne pour condition nulle
      
      
      # somme des longueurs de Troncs intersectes par id selon le type
      TroncSmal_wide_uniq <-  TroncSmal_wide %>%
        group_by(id) %>%
        mutate(SpCELs_1 = if(length(grep("Permanent",colnames(.))) > 0){sum(`Permanent`)}else{sum(`NA`)}) %>%
        mutate(SpCELs_2 = if(length(grep("Fictif",colnames(.))) > 0){sum(`Fictif`)}else{sum(`NA`)}) %>%
        mutate(SpCELs_3 = if(length(grep("Intermittent",colnames(.))) > 0){sum(`Intermittent`)}else{sum(`NA`)}) %>%
        mutate(SpCELs_4 = if(length(grep("Inconnu",colnames(.))) > 0){sum(`Inconnu`)}else{sum(`NA`)}) %>%
        mutate(SpCELs_5 = if(length(grep("En attente de mise à jour",colnames(.))) > 0){sum(`En attente de mise à jour`)}else{sum(`NA`)}) %>%
        mutate(SpCELs_6 = if(length(grep("A sec",colnames(.))) > 0){sum(`A sec`)}else{sum(`NA`)})
        
      
      TroncSmal_wide_uniq <- unique(TroncSmal_wide_uniq[,grep("id|SpCEL",colnames(TroncSmal_wide_uniq))])
      
      
      # jointure des tables
      Extr_Smal <- merge(Extr_Smal, TroncSmal_wide_uniq,all=T)
      
      
      # actualisation de la progression
      pb_s$tick()
      Sys.sleep(1 / 100)
      
    }
    
    
    #############################
    # Buffer locaux : medium -----
    
    cat(paste0(c("\t---\tTroncons d'eau dans un buffer de :",BuffMed," m","\t---\n")))
    
    Extr_Med <- as.data.frame(matrix(nrow=0,ncol = 1))
    colnames(Extr_Med) <- c("id")
    
    # initialisation de la barre de progression
    pb_m <- progress::progress_bar$new(total = nrow(PointbuffMed),
                                       format = "Extraction longueur des Troncons d'eau buffer medium [:bar] :percent    Tps ecoule = :elapsedfull",
                                       clear = F)
    
    
    pb_m$tick(0)
    Sys.sleep(1/20)
    
    for(i in 1:nrow(PointbuffMed)){ # nrow(PointbuffMed)
      
      Point.tmp <- PointBuffMed[i,c("ID_liste","id","Maille","Jour_de_l_annee")]
      
      TroncMed.tmp <- st_intersection(TRONC, Point.tmp)
      
      if(nrow(TroncMed.tmp) > 0){
        TroncMed.tmp$Longueur_intersect <- as.numeric(st_length(TroncMed.tmp))
        
      }else{ # en cas d'absence de Tronc --> recuperation des infos  
        
        TroncMed.tmp[1,"id"] <- Point.tmp$id
        TroncMed.tmp$Longueur_intersect <- 0
        TroncMed.tmp$Etat <- NA
        
      }
      
      
      # gestion de l'extractions des longueurs de Troncs intersectees
      TroncMed.tmp <- st_drop_geometry(TroncMed.tmp[,c("id","Etat","Longueur_intersect")])
      
      # conversion long -> wide
      TroncMed.tmp.dt <- data.table::data.table(TroncMed.tmp)
      TroncMed_wide <- data.table::dcast(TroncMed.tmp.dt,id + Longueur_intersect ~ Etat, value.var = "Longueur_intersect", fill = 0) # conversion long -> wide
      TroncMed_wide$`NA` <- 0 # rajout d'une colonne pour condition nulle
      
      
      # somme des longueurs de Troncs intersectes par id selon le type
      TroncMed_wide_uniq <-  TroncMed_wide %>%
        group_by(id) %>%
        mutate(SpCELm_1 = if(length(grep("Permanent",colnames(.))) > 0){sum(`Permanent`)}else{sum(`NA`)}) %>%
        mutate(SpCELm_2 = if(length(grep("Fictif",colnames(.))) > 0){sum(`Fictif`)}else{sum(`NA`)}) %>%
        mutate(SpCELm_3 = if(length(grep("Intermittent",colnames(.))) > 0){sum(`Intermittent`)}else{sum(`NA`)}) %>%
        mutate(SpCELm_4 = if(length(grep("Inconnu",colnames(.))) > 0){sum(`Inconnu`)}else{sum(`NA`)}) %>%
        mutate(SpCELm_5 = if(length(grep("En attente de mise à jour",colnames(.))) > 0){sum(`En attente de mise à jour`)}else{sum(`NA`)}) %>%
        mutate(SpCELm_6 = if(length(grep("A sec",colnames(.))) > 0){sum(`A sec`)}else{sum(`NA`)})
      
      
      TroncMed_wide_uniq <- unique(TroncMed_wide_uniq[,grep("id|SpCEL",colnames(TroncMed_wide_uniq))])
      
      
      # jointure des tables
      Extr_Med <- merge(Extr_Med, TroncMed_wide_uniq,all=T)
      
      
      # actualisation de la progression
      pb_m$tick()
      Sys.sleep(1 / 100)
      
    }
    
    
    
    #############################
    # Buffer locaux : large -----
    
    cat(paste0(c("\t---\tTroncons d'eau dans un buffer de :",BuffLarg," m","\t---\n")))
    
    Extr_Larg <- as.data.frame(matrix(nrow=0,ncol = 1))
    colnames(Extr_Larg) <- c("id")
    
    # initialisation de la barre de progression
    pb_l <- progress::progress_bar$new(total = nrow(PointbuffLarg),
                                       format = "Extraction longueur des Troncons d'eau buffer large [:bar] :percent    Tps ecoule = :elapsedfull",
                                       clear = F)
    
    
    pb_l$tick(0)
    Sys.sleep(1/20)
    
    for(i in 1:nrow(PointbuffLarg)){ # nrow(PointbuffLarg)
      
      Point.tmp <- PointBuffLarg[i,c("ID_liste","id","Maille","Jour_de_l_annee")]
      
      TroncLarg.tmp <- st_intersection(TRONC, Point.tmp)
      
      if(nrow(TroncLarg.tmp) > 0){
        TroncLarg.tmp$Longueur_intersect <- as.numeric(st_length(TroncLarg.tmp))
        
      }else{ # en cas d'absence de Tronc --> recuperation des infos  
        
        TroncLarg.tmp[1,"id"] <- Point.tmp$id
        TroncLarg.tmp$Longueur_intersect <- 0
        TroncLarg.tmp$Etat <- NA
        
      }
      
      
      # gestion de l'extractions des longueurs de Troncs intersectees
      TroncLarg.tmp <- st_drop_geometry(TroncLarg.tmp[,c("id","Etat","Longueur_intersect")])
      
      # conversion long -> wide
      TroncLarg.tmp.dt <- data.table::data.table(TroncLarg.tmp)
      TroncLarg_wide <- data.table::dcast(TroncLarg.tmp.dt,id + Longueur_intersect ~ Etat, value.var = "Longueur_intersect", fill = 0) # conversion long -> wide
      TroncLarg_wide$`NA` <- 0 # rajout d'une colonne pour condition nulle
      
      
      # somme des longueurs de Troncs intersectes par id selon le type
      TroncLarg_wide_uniq <-  TroncLarg_wide %>%
        group_by(id) %>%
        mutate(SpCELl_1 = if(length(grep("Permanent",colnames(.))) > 0){sum(`Permanent`)}else{sum(`NA`)}) %>%
        mutate(SpCELl_2 = if(length(grep("Fictif",colnames(.))) > 0){sum(`Fictif`)}else{sum(`NA`)}) %>%
        mutate(SpCELl_3 = if(length(grep("Intermittent",colnames(.))) > 0){sum(`Intermittent`)}else{sum(`NA`)}) %>%
        mutate(SpCELl_4 = if(length(grep("Inconnu",colnames(.))) > 0){sum(`Inconnu`)}else{sum(`NA`)}) %>%
        mutate(SpCELl_5 = if(length(grep("En attente de mise à jour",colnames(.))) > 0){sum(`En attente de mise à jour`)}else{sum(`NA`)}) %>%
        mutate(SpCELl_6 = if(length(grep("A sec",colnames(.))) > 0){sum(`A sec`)}else{sum(`NA`)})
      
      
      TroncLarg_wide_uniq <- unique(TroncLarg_wide_uniq[,grep("id|SpCEL",colnames(TroncLarg_wide_uniq))])
      
      
      # jointure des tables
      Extr_Larg <- merge(Extr_Larg, TroncLarg_wide_uniq,all=T)
      
      
      # actualisation de la progression
      pb_l$tick()
      Sys.sleep(1 / 100)
      
    }
    
  # Jointure des tables d'extractions des differents buffers
    SpTRONC <- left_join(Extr_Smal,Extr_Med, by="id")
    SpTRONC <- left_join(SpTRONC,Extr_Larg, by="id")
    
    SpTRONC <- dplyr::left_join(SpTRONC, Point[,c("id","X_barycentre_L93","Y_barycentre_L93","ID_liste")],by="id") # add d'informations sur les listes

    
    write.csv(SpTRONC, file = "C:/git/ODF/output/function_output/TRONCONS_envEPOC.csv", row.names = F)
    
    
    cat("\n\n\t-------\tExtract TRONCONS done - check /function_output \t-------\n\n")
    
    
    
  
  # exttraction proportion de surface des zones d'eaux -----
    
    #############################
    # Buffer locaux : small -----
    
    cat(paste0(c("\t---\tSurfaces d'eau dans un buffer de :",BuffSmal," m","\t---\n")))
    
    Extr_Smal <- as.data.frame(matrix(nrow=0,ncol = 1))
    colnames(Extr_Smal) <- c("id")
    
    # initialisation de la barre de progression
    pb_s <- progress::progress_bar$new(total = nrow(PointbuffSmal),
                                       format = "Extraction Aire des surfaces d'eau buffer small [:bar] :percent    Tps ecoule = :elapsedfull",
                                       clear = F)
    
    
    pb_s$tick(0)
    Sys.sleep(1/20)
    
    for(i in 1:nrow(PointbuffSmal)){ # nrow(PointbuffSmal)
      
      Point.tmp <- PointBuffSmal[i,c("ID_liste","id","Maille","Jour_de_l_annee")]
      
      SurfSmal.tmp <- st_intersection(SURF, Point.tmp)
      
      if(nrow(SurfSmal.tmp) > 0){
        SurfSmal.tmp$Aire_intersect <- as.numeric(st_area(SurfSmal.tmp))
        
      }else{ # en cas d'absence de Surf --> recuperation des infos  
        
        SurfSmal.tmp[1,"id"] <- Point.tmp$id
        SurfSmal.tmp$Aire_intersect <- 0
        SurfSmal.tmp$Nature <- NA
        
      }
      
      # calcul de la proportion des surfaces
      Aire_tot <- BuffSmal^2*pi
      SurfSmal.tmp$Aire_intersect <- SurfSmal.tmp$Aire_intersect / Aire_tot
      
      # gestion de l'extractions des surfaces intersectees
      SurfSmal.tmp <- st_drop_geometry(SurfSmal.tmp[,c("id","Nature","Aire_intersect")])
      
      # conversion long -> wide
      SurfSmal.tmp.dt <- data.table::data.table(SurfSmal.tmp)
      SurfSmal_wide <- data.table::dcast(SurfSmal.tmp.dt,id + Aire_intersect ~ Nature, value.var = "Aire_intersect", fill = 0) # conversion long -> wide
      SurfSmal_wide$`NA` <- 0 # rajout d'une colonne pour condition nulle
      
      
      # somme des Aires de Surfs intersectes par id selon le type
      
      SurfSmal_wide_uniq <-  SurfSmal_wide %>%
        group_by(id) %>%
        mutate(SpCESs_1 = if(length(grep("Eau salée permanente",colnames(.))) > 0){sum(`Eau salée permanente`)}else{sum(`NA`)}) %>%
        mutate(SpCESs_2 = if(length(grep("Eau douce permanente",colnames(.))) > 0){sum(`Eau douce permanente`)}else{sum(`NA`)}) %>%
        mutate(SpCESs_3 = if(length(grep("Eau salée non permanente",colnames(.))) > 0){sum(`Eau salée non permanente`)}else{sum(`NA`)}) %>%
        mutate(SpCESs_4 = if(length(grep("Eau douce non permanente",colnames(.))) > 0){sum(`Eau douce non permanente`)}else{sum(`NA`)}) %>%
        mutate(SpCESs_5 = if(length(grep("Névé, glacier",colnames(.))) > 0){sum(`Névé, glacier`)}else{sum(`NA`)}) 
      
      
      SurfSmal_wide_uniq <- unique(SurfSmal_wide_uniq[,grep("id|SpCES",colnames(SurfSmal_wide_uniq))])
      
      
      # jointure des tables
      Extr_Smal <- merge(Extr_Smal, SurfSmal_wide_uniq,all=T)
      
      
      # actualisation de la progression
      pb_s$tick()
      Sys.sleep(1 / 100)
      
    }
    
    
  
    #############################
    # Buffer locaux : medium -----
    
    cat(paste0(c("\t---\tSurfaces d'eau dans un buffer de :",BuffMed," m","\t---\n")))
    
    Extr_Med <- as.data.frame(matrix(nrow=0,ncol = 1))
    colnames(Extr_Med) <- c("id")
    
    # initialisation de la barre de progression
    pb_m <- progress::progress_bar$new(total = nrow(PointbuffMed),
                                       format = "Extraction Aire des surfaces d'eau buffer medium [:bar] :percent    Tps ecoule = :elapsedfull",
                                       clear = F)
    
    
    pb_m$tick(0)
    Sys.sleep(1/20)
    
    for(i in 1:nrow(PointbuffMed)){ # nrow(PointbuffMed)
      
      Point.tmp <- PointBuffMed[i,c("ID_liste","id","Maille","Jour_de_l_annee")]
      
      SurfMed.tmp <- st_intersection(SURF, Point.tmp)
      
      if(nrow(SurfMed.tmp) > 0){
        SurfMed.tmp$Aire_intersect <- as.numeric(st_area(SurfMed.tmp))
        
      }else{ # en cas d'absence de Surf --> recuperation des infos  
        
        SurfMed.tmp[1,"id"] <- Point.tmp$id
        SurfMed.tmp$Aire_intersect <- 0
        SurfMed.tmp$Nature <- NA
        
      }
      
      # calcul de la proportion des surfaces
      Aire_tot <- BuffMed^2*pi
      SurfMed.tmp$Aire_intersect <- SurfMed.tmp$Aire_intersect / Aire_tot
      
      # gestion de l'extractions des surfaces intersectees
      SurfMed.tmp <- st_drop_geometry(SurfMed.tmp[,c("id","Nature","Aire_intersect")])
      
      # conversion long -> wide
      SurfMed.tmp.dt <- data.table::data.table(SurfMed.tmp)
      SurfMed_wide <- data.table::dcast(SurfMed.tmp.dt,id + Aire_intersect ~ Nature, value.var = "Aire_intersect", fill = 0) # conversion long -> wide
      SurfMed_wide$`NA` <- 0 # rajout d'une colonne pour condition nulle
      
      
      # somme des Aires de Surfs intersectes par id selon le type
      
      SurfMed_wide_uniq <-  SurfMed_wide %>%
        group_by(id) %>%
        mutate(SpCESm_1 = if(length(grep("Eau salée permanente",colnames(.))) > 0){sum(`Eau salée permanente`)}else{sum(`NA`)}) %>%
        mutate(SpCESm_2 = if(length(grep("Eau douce permanente",colnames(.))) > 0){sum(`Eau douce permanente`)}else{sum(`NA`)}) %>%
        mutate(SpCESm_3 = if(length(grep("Eau salée non permanente",colnames(.))) > 0){sum(`Eau salée non permanente`)}else{sum(`NA`)}) %>%
        mutate(SpCESm_4 = if(length(grep("Eau douce non permanente",colnames(.))) > 0){sum(`Eau douce non permanente`)}else{sum(`NA`)}) %>%
        mutate(SpCESm_5 = if(length(grep("Névé, glacier",colnames(.))) > 0){sum(`Névé, glacier`)}else{sum(`NA`)}) 
      
      
      SurfMed_wide_uniq <- unique(SurfMed_wide_uniq[,grep("id|SpCES",colnames(SurfMed_wide_uniq))])
      
      
      # jointure des tables
      Extr_Med <- merge(Extr_Med, SurfMed_wide_uniq,all=T)
      
      
      # actualisation de la progression
      pb_m$tick()
      Sys.sleep(1 / 100)
      
    }
  
  
    
    #############################
    # Buffer locaux : large -----
    
    cat(paste0(c("\t---\tSurfaces d'eau dans un buffer de :",BuffLarg," m","\t---\n")))
    
    Extr_Larg <- as.data.frame(matrix(nrow=0,ncol = 1))
    colnames(Extr_Larg) <- c("id")
    
    # initialisation de la barre de progression
    pb_m <- progress::progress_bar$new(total = nrow(PointbuffLarg),
                                       format = "Extraction Aire des surfaces d'eau buffer large [:bar] :percent    Tps ecoule = :elapsedfull",
                                       clear = F)
    
    
    pb_m$tick(0)
    Sys.sleep(1/20)
    
    for(i in 1:nrow(PointbuffLarg)){ # nrow(PointbuffLarg)
      
      Point.tmp <- PointBuffLarg[i,c("ID_liste","id","Maille","Jour_de_l_annee")]
      
      SurfLarg.tmp <- st_intersection(SURF, Point.tmp)
      
      if(nrow(SurfLarg.tmp) > 0){
        SurfLarg.tmp$Aire_intersect <- as.numeric(st_area(SurfLarg.tmp))
        
      }else{ # en cas d'absence de Surf --> recuperation des infos  
        
        SurfLarg.tmp[1,"id"] <- Point.tmp$id
        SurfLarg.tmp$Aire_intersect <- 0
        SurfLarg.tmp$Nature <- NA
        
      }
      
      # calcul de la proportion des surfaces
      Aire_tot <- BuffLarg^2*pi
      SurfLarg.tmp$Aire_intersect <- SurfLarg.tmp$Aire_intersect / Aire_tot
      
      # gestion de l'extractions des surfaces intersectees
      SurfLarg.tmp <- st_drop_geometry(SurfLarg.tmp[,c("id","Nature","Aire_intersect")])
      
      # conversion long -> wide
      SurfLarg.tmp.dt <- data.table::data.table(SurfLarg.tmp)
      SurfLarg_wide <- data.table::dcast(SurfLarg.tmp.dt,id + Aire_intersect ~ Nature, value.var = "Aire_intersect", fill = 0) # conversion long -> wide
      SurfLarg_wide$`NA` <- 0 # rajout d'une colonne pour condition nulle
      
      
      # somme des Aires de Surfs intersectes par id selon le type
      
      SurfLarg_wide_uniq <-  SurfLarg_wide %>%
        group_by(id) %>%
        mutate(SpCESl_1 = if(length(grep("Eau salée permanente",colnames(.))) > 0){sum(`Eau salée permanente`)}else{sum(`NA`)}) %>%
        mutate(SpCESl_2 = if(length(grep("Eau douce permanente",colnames(.))) > 0){sum(`Eau douce permanente`)}else{sum(`NA`)}) %>%
        mutate(SpCESl_3 = if(length(grep("Eau salée non permanente",colnames(.))) > 0){sum(`Eau salée non permanente`)}else{sum(`NA`)}) %>%
        mutate(SpCESl_4 = if(length(grep("Eau douce non permanente",colnames(.))) > 0){sum(`Eau douce non permanente`)}else{sum(`NA`)}) %>%
        mutate(SpCESl_5 = if(length(grep("Névé, glacier",colnames(.))) > 0){sum(`Névé, glacier`)}else{sum(`NA`)}) 
      
      
      SurfLarg_wide_uniq <- unique(SurfLarg_wide_uniq[,grep("id|SpCES",colnames(SurfLarg_wide_uniq))])
      
      
      # jointure des tables
      Extr_Larg <- merge(Extr_Larg, SurfLarg_wide_uniq,all=T)
      
      
      # actualisation de la progression
      pb_m$tick()
      Sys.sleep(1 / 100)
      
    }
  
  
  
    # Jointure des tables d'extractions des differents buffers
    SpSURF <- left_join(Extr_Smal,Extr_Med, by="id")
    SpSURF <- left_join(SpSURF,Extr_Larg, by="id")
    
    SpSURF <- dplyr::left_join(SpSURF, Point[,c("id","X_barycentre_L93","Y_barycentre_L93","ID_liste")],by="id") # add d'informations sur les listes
    
    
    write.csv(SpSURF, file = "C:/git/ODF/output/function_output/EauSURFACES_envEPOC.csv", row.names = F)
    
    
    cat("\n\n\t-------\tExtract TRONCONS done - check /function_output \t-------\n\n")
  
  
  
  
  
}























