# chemmin
  setwd("C:/git/ODF/data")

# packages
  library(sf)
  library(tmap)
  library(dplyr)
  library(data.table)
  library(ggplot2)

# upload data

  fra.adm <- st_transform(st_read(dsn = "C:/Users/Travail/Desktop/Ressource QGis/france/adm/FRA_adm2.shp"),crs=2154)
  
  mailles2x2_france <- st_read("C:/git/ODF/data/couche_SSU_initiale_id_maille10km_centroides_intersection_France.shp")
  colnames(mailles2x2_france)[grep("grille_mai",colnames(mailles2x2_france))] <- "NUMÉRO"
  

# fonctions
  

# verification de la distribution des aires des mailles de 2km
  # summary(mailles$area)
  # quantile(mailles$area, c(1:99)/100)  
  
  
# Retrait des SSU incomplets/coupes --> selectionnant les SSU w/ centroides conservees ----
  mailles <- st_buffer(mailles2x2_france, dist = 1000, endCapStyle = "SQUARE")

  
# formation d'un identifiant SSU ----

  mailles_dt <- data.table(mailles)
  mailles_dt <- mailles_dt[, ID_SSU := 1:.N, by = "CODE10KM"] # incrementation des mailles 2x2 selon le code maille
  mailles_dt <- as.data.frame(mailles_dt)
  
  mailles_dt <- mailles_dt %>%
    mutate(ID_PSU_SSU = paste(CODE10KM,ID_SSU, sep="_")) # identifiant PSU_SSU
  
  mailles2 <- left_join(mailles, mailles_dt[,c("ID","ID_SSU","ID_PSU_SSU")]) # jointure des tables (apport de l'information de l'id_SSU)
 

# recuperation du nombre de SSU par PSU ----
  PSU.dtf <- mailles2[,c("ID","CODE10KM","ID_SSU","ID_PSU_SSU")] %>%
    group_by(CODE10KM) %>%
    mutate(nb_SSU = n())

  PSU.dtf <- st_drop_geometry(PSU.dtf)


# selection des PSU ayant au moins 10 SSU complets ----
  PSU.dtf$CODE10KM <- as.character(PSU.dtf$CODE10KM)
  id.PSU10SSU <- which(PSU.dtf$nb_SSU >= 10)
  id.PSU10SSU <- unique(PSU.dtf[id.PSU10SSU,"CODE10KM"])
  id.PSU10SSU <- id.PSU10SSU$CODE10KM

  mailles3 <- mailles2[which(mailles2$CODE10KM %in% id.PSU10SSU),]


# ajout des departements -----
  
  # recuperration de l'information : appartenance des mailles dans les departements
    mailles10x10_1 <- st_crop(mailles10x10,fra.adm)
    mailles10x10_1 <- mailles10x10_1[which(mailles10x10_1$ESPACE != "Marin"),]
  
    dep <- st_intersection(mailles10x10_1,fra.adm)
    dep <- dep[,c("CODE10KM","NAME_2")] ; colnames(dep)[grep("NAME",colnames(dep))] <- "Departements"
    dep <- unique(st_drop_geometry(dep))
    dep$Departements  <- as.character(dep$Departements)
    
    dep <- dep[which(dep$CODE10KM %in% id.PSU10SSU),]
    
  # ajout de cette information sur les PSU/SSU
    mailles3 <- left_join(mailles3, dep)
    mailles3 <- mailles3[-which(duplicated(mailles3$ID_PSU_SSU) == T),] # rm des mailles SSU dupliquees (frontieres des departements)
    
  

  
# tirage aleatoire des mailles 10x10 au sein des departements -----
  dep.dtf <- st_drop_geometry(mailles3[,c("CODE10KM","ID_PSU_SSU","Departements")]) #"ID_SSU","ID_PSU_SSU"

  # selection des departements possedant au moins 10 PSU ----
    nb_mailles10x10_bydep <- unique(dep.dtf[,c("Departements","CODE10KM")]) %>%
      group_by(Departements) %>%
      mutate(Nb_mailles = n())
    nb_mailles10x10_bydep$Departements <- as.character(nb_mailles10x10_bydep$Departements)
    
    nb_mailles10x10_bydep <- nb_mailles10x10_bydep[which(duplicated(nb_mailles10x10_bydep$Departements) == F),]
    ## Paris / Hauts de seine / territoire de belfort / Seine-st-Denis / Val-de-Marne--> departements avec - de 10 mailles 10x10 completes
    
    dep_moins10mailles <- nb_mailles10x10_bydep[which(nb_mailles10x10_bydep$Nb_mailles <= 10),] ; dep_moins10mailles <- dep_moins10mailles$Departements
  
    dep.dtf2 <- dep.dtf[-which(dep.dtf$Departements %in% dep_moins10mailles),] # subset avec les departements possedant assez de mailles pour faire une selection aleaoire
    dep.dtf2a <- dep.dtf[which(dep.dtf$Departements %in% dep_moins10mailles),]
    
    
  # tirage aleatoire 10 PSU par departements ----
    dep.dtf2bis <- unique(dep.dtf2[,c("Departements","CODE10KM")])
    
    mailles_selected <- dep.dtf2bis %>%
      group_by(Departements) %>%
      sample_n(10,replace = F) %>%
      mutate(PSU_tires = 1)
    
  
  # cas par cas pour les departements avec ~ - de 10 mailles ----
    
    dep.dtf2abis <- unique(dep.dtf2a[,c("Departements","CODE10KM")])
    
    # Tirage aleatoire sans remise de 50% des mailles (arrondi au supérieure)
    dep10_small <- dep.dtf2abis %>%
      group_by(Departements) %>%
      sample_n(round(n()*0.51),replace = F) %>% # 0.51 : cas de belfort pour que l'arrondi soit au superieur
      mutate(PSU_tires = 1)
    
    mailles_selected <- rbind(mailles_selected,dep10_small)
    
  
  # jointure des mailles choisis aleatoire -----
      
    dep.dtf3 <- left_join(dep.dtf,mailles_selected)
    dep.dtf3[which(is.na(dep.dtf3$PSU_tires)),"PSU_tires"] <- 0 
    
  
    mailles4 <- left_join(mailles3, dep.dtf3[,c("ID_PSU_SSU","PSU_tires")])
    mailles4[which(is.na(mailles4$PSU_tires)),"PSU_tires"] <- 0 
    
  # visualisation des mailles 10x10 choisis aleatoirement -----
   ' 
    ggplot() +
      geom_sf(data=mailles4, aes(fill = as.factor(PSU_tires), color = as.factor(PSU_tires))) +
      geom_sf(data=fra.adm, alpha=0.5) +
      scale_color_discrete(guide=F) +
      labs(fill = "Mailles selectionnes")
   '
    
# Tirage aelatoire des SSU a echantillonnees -----
  # creation table (COde10/ID_SSU/ID_PSU_SSU)
    m2x2_selected <- st_drop_geometry(mailles4[which(mailles4$PSU_tires == 1),c("Departements","CODE10KM","ID_PSU_SSU","ID_SSU")]) # selection des mailles 10x10 choisis avant de proceder aux tirages des mailles 2x2
    
    
  # Tirage aleatoire des SSu officiels -----
    
    m2x2_selected_off <- m2x2_selected %>%
      group_by(CODE10KM) %>%
      sample_n(5) %>%
      mutate(Statut_SSU = "Officiel")
    
    m2x2_selected <- left_join(m2x2_selected,m2x2_selected_off)
    
    
  # Tirage aleatoire des SSU de reserve -----
    
    m2x2_selected_resv <- m2x2_selected[which(is.na(m2x2_selected$Statut_SSU)== T),] %>%
      group_by(CODE10KM) %>%
      sample_n(5) %>%
      mutate(Statut_SSU = "Reserve")
    
    m2x2_selected_resv_dt <- data.table(m2x2_selected_resv)
    m2x2_selected_resv_dt <- m2x2_selected_resv_dt[, Rang_reserve := 1:.N, by = "CODE10KM"] # incrementation des mailles de reserve 2x2 selon le code maille
    m2x2_selected_resv <- as.data.frame(m2x2_selected_resv_dt)
      
    # m2x2_selected_resv <- m2x2_selected_resv %>%
    #   mutate(Statut_SSU = paste(Statut_SSU,Rang_reserve, sep="_"))
    
      
  # jointure des deux tables ----
    
    m2x2_selected <- merge(m2x2_selected, m2x2_selected_resv, all=T)
    m2x2_selected <- m2x2_selected[-which(duplicated(m2x2_selected$ID_PSU_SSU) == T),] # retrait des mailles dupliquees
    

# formation de l'objet spatial regroupant l'information sur les mailles 2x2 a echantillonner -----
  mailles5 <- left_join(mailles4, m2x2_selected)
  mailles5 <- mailles5[,-grep("ID_SSU",colnames(mailles5))]  
  
  # sauvegarde des shape sur disque
  st_write(mailles5, dsn = "C:/git/ODF/output/0X_tirages_aleatoire_2x2/Mailles2x2_tirees_aleatoirement_v5.geojson", driver = "GEOJSON")
  st_write(mailles5, dsn = "C:/git/ODF/output/0X_tirages_aleatoire_2x2/Mailles2x2_tirees_aleatoirement_v5_ESRI", driver = "ESRI Shapefile")
    
  # visualisation des mailles 2x2 tiree aleatoirement (ex : )
    ggplot() +
      geom_sf(data = mailles5[which(mailles5$Departements == "Charente-Maritime"),], aes(fill = as.factor(Statut_SSU),col = as.factor(PSU_tires))) +
      labs(fill = "Mailles 2x2 a echantillonner") +
      scale_color_manual(values = c(NA,"snow3"),guide=F) +
      ggtitle("Zoom sur le département de Charente-Maritime")
    
    ggplot() +
      geom_sf(data = mailles5[which(mailles5$Departements == "Hautes-Alpes"),], aes(fill = as.factor(Statut_SSU),col = as.factor(PSU_tires))) +
      labs(fill = "Mailles 2x2 a echantillonner") +
      scale_color_manual(values = c(NA,"snow3"),guide=F) +
      ggtitle("Zoom sur le département de Hautes-Alpes")

    ggplot() +
      geom_sf(data = mailles5[which(mailles5$Departements == "Territoire de Belfort"),], aes(fill = as.factor(Statut_SSU),col = as.factor(PSU_tires))) +
      labs(fill = "Mailles 2x2 a echantillonner") +
      scale_color_manual(values = c(NA,"snow3"),guide=F) +
      ggtitle("Zoom sur le département des Hauts-de-Seine")    
    
    
    
    
  # sauvegarde du shape par departement ----
    vec.dep <- unique(mailles5$Departements)
    
    for(i in vec.dep){
      
      mailles.tmp <- mailles5[which(mailles5$Departements == i),]
      
      st_write(mailles.tmp, dsn = paste0("C:/git/ODF/output/0X_tirages_aleatoire_2x2/Departements/",i,".geojson"), driver = "GEOJSON")
      st_write(mailles.tmp, dsn = paste0("C:/git/ODF/output/0X_tirages_aleatoire_2x2/Departements/",i,"_ESRI"), driver = "ESRI Shapefile")
      
    }
    
    
    
    
    
    
    
    
    
  
  
  
  
  
  
  











