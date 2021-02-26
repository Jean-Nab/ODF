# chemmin
  setwd("C:/git/ODF/data")

# packages
  library(sf)
  library(sp)
  library(raster)
  library(dplyr)
  library(ggplot2)
  library(maptools)
  library(reshape2)
  library(data.table)
  

# upload data -----
  asc_files<- list.files("C:/Users/Travail/Desktop/Ressource QGis/france/altitude/BDALTIV2_2-0_75M_ASC_LAMB93-IGN69_FRANCE_2020-09-30/BDALTIV2/1_DONNEES_LIVRAISON_2020-10-00295/BDALTIV2_MNT_75M_ASC_LAMB93_IGN69_FRANCE"
                         ,pattern =".asc$",full.names=T)
  mailles2x2_france <- st_read("C:/git/ODF/data/couche_SSU_initiale_id_maille10km_centroides_intersection_France.shp")
  mailles10x10 <-  st_transform(st_read(dsn = "C:/git/epoc/data/maille_atlas/L93_10x10_TerreMer.shp"),crs=2154)
  
  fra.adm <- st_transform(st_read(dsn = "C:/Users/Travail/Desktop/Ressource QGis/france/adm/FRA_adm2.shp"),crs=2154)
  land <- raster("C:/Users/Travail/Desktop/Ressource QGis/france/CLC2018_raster100m_L93.tif")
  
  
# formatage data ----
  # buffer 500m pour extraction donnes habitats
    mailles2x2_500mbuf <- st_buffer(mailles2x2_france,dist=500)
  
  # conversion sf -> sp
    mailles2x2_france.sp <- as(mailles2x2_france,Class = "Spatial")
  
  # formation dtf resume
    #mailles2x2_env <- as.data.frame(coordinates(mailles2x2_france.sp))
    #mailles2x2_env <- SpatialPointsDataFrame(data=mailles2x2_env,coords = mailles2x2_env[,c(1,2)], 
    #                                         proj4string = CRS("+init=epsg:2154")) #attribution du crs L93
  
  # rasterisation donnes alti
    rast.list <- list()
    for(i in 1:length(asc_files)) { rast.list[i] <- raster(asc_files[i]) }
    rast.list$fun <- mean

    AltiTot <- do.call(mosaic,rast.list) # 3:30 min
    
    
# extraction data pour l'ensemble des centroides SSU ----
  # extraction des donnees alti
    # altitude locale
      SpAltiS=extract(AltiTot,mailles2x2_france.sp)
      mailles2x2_france.sp <- spCbind(mailles2x2_france.sp,SpAltiS)
      
    # altitude moyenne : rayon de 500m
      #SpAltiM=extract(AltiTot,mailles2x2_env,buffer=500,fun=mean) # 0.01 sec / points
      #mailles2x2_env=spCbind(mailles2x2_env,SpAltiM)


    
  # extract : manipulation des cellules de lextract ------
    
    mailles2x2_500mbuf$ID_extract <- c(rep(1:nrow(mailles2x2_500mbuf)))
    
    vec.iteration <- seq(from = 1, to = length(mailles2x2_500mbuf$ID_extract), by = 305) # 305
    
    CLC_hab_wide <- as.data.frame(matrix(nrow=0,ncol = 1))
    colnames(CLC_hab_wide) <- "ID"

    
    # boucle de 305 en 305: ----
    for(i in vec.iteration){
      
      mailles2x2_500mbuf.tmp <- mailles2x2_500mbuf[i:min((i+304),length(mailles2x2_500mbuf$ID_extract)),] #304
      
      beginCluster()
      CLC_hab.tmp <- extract(x = land, 
                             y = mailles2x2_500mbuf.tmp,
                             df=T,
                             along=T)
      endCluster()
      
      # formatage des donnees
        colnames(CLC_hab.tmp)[1] <- "ID_extract"
        CLC_hab.tmp$ID_extract <- CLC_hab.tmp$ID_extract + (i - 1)
      # extract selon CLC niv3 --> recuperation du niveau 2
        CLC_hab.tmp$CLC2018_raster100m_L93 <- trunc(CLC_hab.tmp$CLC2018_raster100m_L93/10)

      
      
      # gestion dtf des resultats preliminaires -----
      
      CLC_hab.tmp <- left_join(CLC_hab.tmp,st_drop_geometry(mailles2x2_500mbuf.tmp[,c("ID","ID_extract")]))
      
      #CLC_hab_wide.tmp <- CLC_hab.tmp %>% dcast(ID_extract + ID ~ CLC2018_raster100m_L93)
        CLC_hab_wide.tmp <- CLC_hab.tmp %>% dcast(ID_extract + ID ~ CLCniv3)
        
        CLCniv3_hab_wide.tmp <- CLC_hab.tmp %>% reshape2::dcast(ID_extract + ID ~ CLCniv3)
        CLCniv2_hab_wide.tmp <- CLC_hab.tmp %>% reshape2::dcast(ID_extract + ID ~ CLCniv2)
        CLCniv1_hab_wide.tmp <- CLC_hab.tmp %>% reshape2::dcast(ID_extract + ID ~ CLCniv1)
        
        CLC_hab_wide.tmp <- left_join(CLCniv1_hab_wide.tmp, left_join(CLCniv2_hab_wide.tmp,CLCniv3_hab_wide.tmp))
        
        
      CLC_hab_wide <- merge(CLC_hab_wide,CLC_hab_wide.tmp,all=T)
    
          
     # avancement boucle
      cat(i,"/ ",vec.iteration[length(vec.iteration)],"\n")
      
    }
    
    #save.image("C:/git/ODF/CLC_hab_centroide_SSU.RData")
    
    # correction des NAs (= habitats absent lors de l'extract du point)
      CLC_hab_wide[is.na(CLC_hab_wide)] <- 0
    
    
    # calcul des proportions d'habitats -----
      CLC_hab_wide[,grep("[0-9]",colnames(CLC_hab_wide))] <- ( CLC_hab_wide[,grep("[0-9]",colnames(CLC_hab_wide))] / 
                                                                 rowSums(CLC_hab_wide[,grep("[0-9]",colnames(CLC_hab_wide))])
                                                              ) * 100

    # mise en forme des donnes d'habitats ----
      
      # rearragement de l'ordre des colonnes
        CLC_hab_wide <- CLC_hab_wide[,order(colnames(CLC_hab_wide))] 
        CLC_hab_wide <- CLC_hab_wide[,c(grep("[a-zA-Z]",colnames(CLC_hab_wide)) , 
                                        grep("[0-9]{2,}",colnames(CLC_hab_wide)) , 
                                        grep("0",colnames(CLC_hab_wide)))] # ordre : IDentifiants / code CLC niv 2 / 0 = extraction de cellule NA (a retirer si jamais)
        
        colnames(CLC_hab_wide) <- gsub("([0-9]{1,})","SpCLC_\\1",colnames(CLC_hab_wide))    
      
    # jointure des donnees altitudes / CLC -----
      mailles2x2_france.sp@data <- left_join(mailles2x2_france.sp@data,CLC_hab_wide)
    
  # sauvegarde de la table / sf des altitudes-CLC extraites au centroides des SSU ----
      save_2x2.sf <- st_as_sf(mailles2x2_france.sp[,grep("ID|Sp",colnames(mailles2x2_france.sp@data))])
      save_2x2.spData <- mailles2x2_france.sp@data
      
    #save("save_2x2.spData","save_2x2.sf", file ="C:/git/ODF/SSU_centroide_alt_CLChab.Rdata")
    
 
# tirage aléatoire des PSU / SSU sur les 3 annees (= 3 tirages aleatoire sans remise ) -----
      
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
          
      
      # ajout des departements ----
      
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
        ## Paris / Hauts de seine / territoire de belfort  --> departements avec - de 30 mailles 10x10 completes
        ## Rhône / Val-d'Oise / Yvelines / Tarn-et-Garonne
        
        
        dep_moins10mailles <- nb_mailles10x10_bydep[which(nb_mailles10x10_bydep$Nb_mailles <= 30),] ; dep_moins10mailles <- dep_moins10mailles$Departements
        
        dep.dtf2 <- dep.dtf[-which(dep.dtf$Departements %in% dep_moins10mailles),] # subset avec les departements possedant suffisame de mailles pour faire une selection aleaoire
        dep.dtf2a <- dep.dtf[which(dep.dtf$Departements %in% dep_moins10mailles),] # subset avec les departements ne possedant pas assez de mailles
        
    
        # tirage aleatoire 10 PSU par departements ----
          dep.dtf2bis <- unique(dep.dtf2[,c("Departements","CODE10KM")])
          
          mailles_selected <- dep.dtf2bis %>%
            group_by(Departements) %>%
            sample_n(30,replace = F) %>%
            mutate(Maille_choisis = 1)
        
        #mailles_selected <- left_join(dep.dtf2,mailles_selected)
        
        
        # cas par cas pour les departements avec ~ - de 10 mailles ----
        
          dep.dtf2abis <- unique(dep.dtf2a[,c("Departements","CODE10KM")])
    
        # Tirage aleatoire sans remise de 75% des mailles (arrondi au supérieure)
          dep10_small <- dep.dtf2abis %>%
            group_by(Departements) %>%
            sample_n(round(n()*0.75),replace = F) %>%
            mutate(Maille_choisis = 1)
          
          mailles_selected <- rbind(mailles_selected,dep10_small)
          
          
        # jointure des mailles choisis aleatoire -----
          
          dep.dtf3 <- left_join(dep.dtf,mailles_selected)
          dep.dtf3[which(is.na(dep.dtf3$Maille_choisis)),"Maille_choisis"] <- 0 
          
          
          mailles4 <- left_join(mailles3, dep.dtf3[,c("ID_PSU_SSU","Maille_choisis")])
          mailles4[which(is.na(mailles4$Maille_choisis)),"Maille_choisis"] <- 0 
    
          # visualisation des mailles 10x10 choisis aleatoirement -----
            k <- st_drop_geometry(mailles4)
            tabl_PSU_bydep <- table(k[which(duplicated(k$CODE10KM) == F),c("Departements","Maille_choisis")]) # tableau croise dynamique nb mailles choisis ou non ; save("tabl_PSU_bydep", file ="C:/git/ODF/output/0x1_table_nb_mailles_tire.Rdata")
            
          # recherche de mailles 10x10 possiblement doublons (~ en frontiere)
            k1 <- which(duplicated(mailles_selected$CODE10KM)==T)
            
           ' 
            ggplot() +
              geom_sf(data=mailles4, aes(fill = as.factor(Maille_choisis), color = as.factor(Maille_choisis))) +
              geom_sf(data=fra.adm, alpha=0.5) +
              scale_color_discrete(guide=F) +
              labs(fill = "Mailles selectionnes")
           '
            
            
          # Tirage aelatoire des SSU a echantillonnees -----
            # creation table (COde10/ID_SSU/ID_PSU_SSU)
            m2x2_selected <- st_drop_geometry(mailles4[which(mailles4$Maille_choisis == 1),c("Departements","CODE10KM","ID_PSU_SSU","ID_SSU")]) # selection des mailles 10x10 choisis avant de proceder aux tirages des mailles 2x2
            
          
          # Tirage aleatoire des SSu officiels -----
          
            m2x2_selected_off <- m2x2_selected %>%
              group_by(CODE10KM) %>%
              sample_n(5) %>%
              mutate(SSU_echantillonnees = "Officiel")
            
            m2x2_selected <- left_join(m2x2_selected,m2x2_selected_off)
            
          
          # Tirage aleatoire des SSU de reserve -----
          
            m2x2_selected_resv <- m2x2_selected[which(is.na(m2x2_selected$SSU_echantillonnees)== T),] %>%
              group_by(CODE10KM) %>%
              sample_n(5) %>%
              mutate(SSU_echantillonnees = "Reserve")
            
            m2x2_selected_resv_dt <- data.table(m2x2_selected_resv)
            m2x2_selected_resv_dt <- m2x2_selected_resv_dt[, ID_SSU_de_reserve := 1:.N, by = "CODE10KM"] # incrementation des mailles de reserve 2x2 selon le code maille
            m2x2_selected_resv <- as.data.frame(m2x2_selected_resv_dt)
            
            m2x2_selected_resv <- m2x2_selected_resv %>%
              mutate(SSU_echantillonnees = paste(SSU_echantillonnees,ID_SSU_de_reserve, sep="_"))
            
          
          # jointure des deux tables ----
          
            m2x2_selected <- merge(m2x2_selected, m2x2_selected_resv, all=T)
            m2x2_selected <- m2x2_selected[-which(duplicated(m2x2_selected$ID_PSU_SSU) == T),] # retrait des mailles dupliquees
            
            
     # formation de l'objet spatial regroupant l'information sur les mailles 2x2 a echantillonner -----
       mailles5 <- left_join(mailles4, m2x2_selected)
       mailles5 <- mailles5[,-grep("ID_SSU",colnames(mailles5))]  
       
       st_write(mailles5, dsn = "C:/git/ODF/output/0X1_tirages_aleatoire_2x2_30PSU/Mailles2x2_tirees_aleatoirement_v4", driver = "GEOJSON")
       st_write(mailles5, dsn = "C:/git/ODF/output/0X1_tirages_aleatoire_2x2_30PSU/Mailles2x2_tirees_aleatoirement_v4_ESRI", driver = "ESRI Shapefile")
     
     # visualisation des mailles 2x2 tiree aleatoirement (ex : )
       ggplot() +
         geom_sf(data = mailles5[which(mailles5$Departements == "Charente-Maritime"),], aes(fill = as.factor(SSU_echantillonnees),col = as.factor(Maille_choisis))) +
         labs(fill = "Mailles 2x2 a echantillonner") +
         scale_color_manual(values = c(NA,"snow3"),guide=F) +
         ggtitle("Zoom sur le département de Charente-Maritime")
       
       ggplot() +
         geom_sf(data = mailles5[which(mailles5$Departements == "Hautes-Alpes"),], aes(fill = as.factor(SSU_echantillonnees),col = as.factor(Maille_choisis))) +
         labs(fill = "Mailles 2x2 a echantillonner") +
         scale_color_manual(values = c(NA,"snow3"),guide=F) +
         ggtitle("Zoom sur le département de Hautes-Alpes")
       
            
# Comparaison entre habitats/altitudes (france metropolitaine - SSU a echantilloner) -----
  # save.image("C:/git/ODF/0X1_avant_comparaison.RData") ; load("C:/git/ODF/0X1_avant_comparaison.RData")

       
  # recuperation du nombre cellule etraites
  # formation d'une table ID / SpCLC | ID / SP_alti
       
       
  # identification des SSU a echantillonne (Officiel / reserv / off + resev)
    mailles5_nonsf <- st_drop_geometry(mailles5)
       
    ID_SSU_off <- mailles5_nonsf[which(mailles5_nonsf$SSU_echantillonnees %in% "Officiel"),"ID"]
    ID_SSU_res <- mailles5_nonsf[which(mailles5_nonsf$SSU_echantillonnees %in% c("Reserve_1","Reserve_2","Reserve_3","Reserve_4","Reserve_5")),"ID"]
    ID_SSU_offres <-  mailles5_nonsf[-which(is.na(mailles5_nonsf$SSU_echantillonnees)),"ID"]
    
  # formation de deux table (alti & CLC)
    # altitude
      centroid_alt <- save_2x2.spData[,grep("ID|SpAlti",colnames(save_2x2.spData))] ; centroid_alt <- centroid_alt[,-grep("ID_extract",colnames(centroid_alt))]
      
    # CLC
      centroid_env <- save_2x2.spData[,grep("ID|SpCLC",colnames(save_2x2.spData))] ; centroid_env <- centroid_env[,-grep("ID_extract",colnames(centroid_env))]
      # recuperation du nombre de cellule extraites
        centroid_env[,grep("Sp",colnames(centroid_env))] <- centroid_env[,grep("Sp",colnames(centroid_env))]*76/100

    
  # Visualisation effort d'echantillonnage ALTITUDE ----
    # Identification des centroides (prospecte ou non) -----
      centroid_alt[which(centroid_alt$ID %in% ID_SSU_off),"Color_differentiel"] <- "Officiel"
      centroid_alt[which(centroid_alt$ID %in% ID_SSU_res),"Color_differentiel"] <- "Reserve"          
      centroid_alt[which(is.na(centroid_alt$Color_differentiel)),"Color_differentiel"] <- "Non echantillonne" 
      
      centroid_alt[which(centroid_alt$ID %in% ID_SSU_offres),"Color_differentiel1"] <- "Officiel_Reserve"
      centroid_alt[which(is.na(centroid_alt$Color_differentiel1)),"Color_differentiel1"] <- "Non echantillonne"    
        
    # plot -----
      ggplot(data = centroid_alt) +
        geom_boxplot(aes(x = Color_differentiel, y = SpAltiS),
                     outlier.shape=NA) +
        scale_y_continuous(limits = quantile(centroid_alt$SpAltiS, c(0, 0.95),na.rm=T))+
        coord_cartesian(ylim = quantile(centroid_alt$SpAltiS, c(0, 0.95),na.rm=T)) +
        xlab("Statut des SSU") + ylab("Altitude au centroïde") +
        ggtitle("Altitude au centroide des SSU",subtitle = paste(" Nb SSU non echantillonne :",length(which(centroid_alt$Color_differentiel == "Non echantillonne")),"\n",
                                                                 "Nb SSU officiel :",length(which(centroid_alt$Color_differentiel == "Officiel")),"\n",
                                                                 "Nb SSU de reserve :",length(which(centroid_alt$Color_differentiel == "Reserve")),"\n",
                                                                 "Representation excluant 5% des donnees de tres hautes altitudes")
                )
      
      ggplot(data = centroid_alt) +
        geom_boxplot(aes(x = Color_differentiel, y = SpAltiS),
                     outlier.shape=NA) +
        scale_y_continuous(limits = quantile(centroid_alt$SpAltiS, c(0, 1),na.rm=T))+
        coord_cartesian(ylim = quantile(centroid_alt$SpAltiS, c(0, 1),na.rm=T)) +
        xlab("Statut des SSU") + ylab("Altitude au centroïde") +
        ggtitle("Altitude au centroide des SSU",subtitle = paste(" Nb SSU non echantillonne :",length(which(centroid_alt$Color_differentiel == "Non echantillonne")),"\n",
                                                                 "Nb SSU officiel :",length(which(centroid_alt$Color_differentiel == "Officiel")),"\n",
                                                                 "Nb SSU de reserve :",length(which(centroid_alt$Color_differentiel == "Reserve")))
        )

      
      ggplot(data = centroid_alt) +
        geom_boxplot(aes(x = Color_differentiel1, y = SpAltiS),
                     outlier.shape=NA) +
        scale_y_continuous(limits = quantile(centroid_alt$SpAltiS, c(0, 0.95),na.rm=T))+
        coord_cartesian(ylim = quantile(centroid_alt$SpAltiS, c(0, 0.95),na.rm=T)) +
        xlab("Statut des SSU") + ylab("Altitude au centroïde") +
        ggtitle("Altitude au centroide des SSU",subtitle = paste(" Nb SSU non echantillonne :",length(which(centroid_alt$Color_differentiel1 == "Non echantillonne")),"\n",
                                                                 "Nb SSU officiel + reserve :",length(which(centroid_alt$Color_differentiel1 == "Officiel_Reserve")),"\n",
                                                                 "Representation excluant 5% des donnees de tres hautes altitudes")
        )
      
      
      # Histogrammes -----
        # 0 sample des SSU non echantilonne
          ggplot(data = centroid_alt,aes(x = SpAltiS,fill=Color_differentiel1)) +
            scale_x_continuous(limits = quantile(centroid_alt$SpAltiS, c(0,1),na.rm=T))+
            geom_histogram(position="identity",binwidth = 25) +
            xlab("Altitude") + ylab("Nombre de SSU") +
            labs(fill = "Statut des SSU") + 
            ggtitle("Nombre de SSU en fonction de l'altitude observée au centroïde par tranche de 25",
                    subtitle = paste(" Nb SSU non echantillonne :",length(which(centroid_alt$Color_differentiel1 == "Non echantillonne")),"\n",
                                     "Nb SSU officiel + reserve :",length(which(centroid_alt$Color_differentiel1 == "Officiel_Reserve")))
                    )
            
      
      
        # w/ sample des SSU non echantilonne
          centroid_alt_noech <- centroid_alt[-which(centroid_alt$ID %in% ID_SSU_offres),] %>%
            sample_n(size = length(ID_SSU_offres),replace=F)
          
          centroid_alt_equivlt <- rbind(centroid_alt_noech,centroid_alt[which(centroid_alt$ID %in% ID_SSU_offres),])
      
          ggplot(data = centroid_alt_equivlt,aes(x = SpAltiS,fill=Color_differentiel1)) +
            scale_x_continuous(limits = quantile(centroid_alt$SpAltiS, c(0,1),na.rm=T))+
            scale_fill_manual(values = c("red","grey"))+
            geom_histogram(position="identity",binwidth = 25,alpha=0.75) +
            xlab("Altitude") + ylab("Nombre de SSU") +
            labs(fill = "Statut des SSU") + 
            ggtitle("Nombre de SSU en fonction de l'altitude observée au centroïde par tranche de 25",
                    subtitle = paste(" Nb SSU non echantillonne après tirage aléatoire :",length(which(centroid_alt_equivlt$Color_differentiel1 == "Non echantillonne")),"\n",
                                     "Nb SSU officiel + reserve :",length(which(centroid_alt_equivlt$Color_differentiel1 == "Officiel_Reserve")))
            )
      
      
        
  # Visualisation effort d'echantillonnage CLC ----     
      # Identification des centroides (prospecte ou non) -----
        centroid_env[which(centroid_env$ID %in% ID_SSU_off),"Color_differentiel"] <- "Officiel"
        centroid_env[which(centroid_env$ID %in% ID_SSU_res),"Color_differentiel"] <- "Reserve"          
        centroid_env[which(is.na(centroid_env$Color_differentiel)),"Color_differentiel"] <- "Non echantillonne" 
        
        centroid_env[which(centroid_env$ID %in% ID_SSU_offres),"Color_differentiel1"] <- "Officiel_Reserve"
        centroid_env[which(is.na(centroid_env$Color_differentiel1)),"Color_differentiel1"] <- "Non echantillonne" 
        
    # Somme des habitats dans un rayon de 500m autour du centroide ----
        
      centroid_tot_hab <- reshape2::melt(centroid_env[,-grep("Color_differentiel1",colnames(centroid_env))],id.vars = c("ID","Color_differentiel"))    
        
      centroid_tot_hab1 <- centroid_tot_hab %>%
        group_by(variable,Color_differentiel) %>%
        mutate(total = sum(as.numeric(value)))
      
      centroid_tot_hab2 <- unique(centroid_tot_hab1[,c("Color_differentiel","variable","total")])
      
      # plot 
        ggplot(centroid_tot_hab2,aes(fill = Color_differentiel,y=total,x=variable)) +
          geom_bar(position = "dodge", stat = "identity") +
          xlab("Habitats CLC niveau 2") + ylab("Surface en ha") +
          labs(fill = "Statut des SSU") +
          ggtitle("Surface totale des habitats CLC niv2 selon un buffer de 500m au centroide des SSU",
                  subtitle = paste(" Nb SSU non echantillonne :",length(which(centroid_env$Color_differentiel == "Non echantillonne")),"\n",
                                   "Nb SSU officiel :",length(which(centroid_env$Color_differentiel == "Officiel")),"\n",
                                   "Nb SSU de reserve :",length(which(centroid_env$Color_differentiel == "Reserve")),"\n")
                  )
    
      # version commant officiel et reserve
        centroid_tot_haba <- reshape2::melt(centroid_env[,-19],id.vars = c("ID","Color_differentiel1"))  # retrait de la colonne w/ distinction entre officiel et reserve  
        
        centroid_tot_haba1 <- centroid_tot_haba %>%
          group_by(variable,Color_differentiel1) %>%
          mutate(total = sum(as.numeric(value)))
        
        centroid_tot_haba2 <- unique(centroid_tot_haba1[,c("Color_differentiel1","variable","total")])
      
        ggplot(centroid_tot_haba2,aes(fill = Color_differentiel1,y=total,x=variable)) +
          geom_bar(position = "dodge", stat = "identity") +
          xlab("Habitats CLC niveau 2") + ylab("Surface en ha") +
          labs(fill = "Statut des SSU") +
          ggtitle("Surface totale des habitats CLC niv2 selon un buffer de 500m au centroide des SSU",
                  subtitle = paste(" Nb SSU non echantillonne :",length(which(centroid_env$Color_differentiel1 == "Non echantillonne")),"\n",
                                   "Nb SSU officiel + reserve :",length(which(centroid_env$Color_differentiel1 == "Officiel_Reserve")))
          )
        
        
    # tirage aleatoire des SSU non echantillonne --> equivalence pour la representation ------
      # Officiel / Reserve / Non echantillonnee ----
        tbl_resum_orne <- as.data.frame(matrix(nrow=0,ncol=2)) ; colnames(tbl_resum_orne) <- c("Color_differentiel","variable") # initialisation de la table resume finale
        
        for(i in 1:30){
          # detection des SSU non echantillonne
            SSU_non_echan <- centroid_env[-which(centroid_env$ID %in% ID_SSU_offres),] %>%
              sample_n(size = length(ID_SSU_off),replace=F) # tirage aleatoire
          
          # creation de la table regroupant tout les SSU officiels/reserve + un tirage des SSU non echantillonnes de taille equivalente
            SSU_env_orne_ta <- rbind(SSU_non_echan,centroid_env[which(centroid_env$ID %in% ID_SSU_offres),])
          
          # formation du tableau croise dynamique
          SSU_hab_orne_ta <- reshape2::melt(SSU_env_orne_ta[,-grep("Color_differentiel1",colnames(SSU_env_orne_ta))],id.vars = c("ID","Color_differentiel"))    
          
          SSU_hab_orne_ta1 <- SSU_hab_orne_ta %>%
            group_by(variable,Color_differentiel) %>%
            mutate(total = sum(as.numeric(value)))
          
          SSU_hab_orne_ta2 <- unique(SSU_hab_orne_ta1[,c("Color_differentiel","variable","total")])
          colnames(SSU_hab_orne_ta2)[3] <- paste0("Tirage_",i)
          
          tbl_resum_orne <- merge(tbl_resum_orne,SSU_hab_orne_ta2,all=T) # jointure w/ table resume finale
        
          
          cat(i,"\n")
        }
        
        # calcul de la surface moyenne apres tirage aleatoire des SSU non echantillonnee
          tbl_resum_orne$Tirage_mean <- apply(tbl_resum_orne[,grep("Tirage",colnames(tbl_resum_orne))],1,mean)
          
      # plot ----
        ggplot(tbl_resum_orne,aes(fill = Color_differentiel,y=Tirage_mean,x=variable)) +
          geom_bar(position = "dodge", stat = "identity") +
          xlab("Habitats CLC niveau 2") + ylab("Surface en ha") +
          labs(fill = "Statut des SSU") +
          ggtitle(paste("Surface totale des habitats CLC niv2 selon un buffer de 500m au centroide des SSU\nSSU non echantillonnés tirés aleatoirement",i,"fois"),
                  subtitle = paste(" Nb SSU non echantillonne :",nrow(SSU_non_echan),"\n",
                                   "Nb SSU officiel :",length(which(centroid_env$Color_differentiel == "Officiel")),"\n",
                                   "Nb SSU de reserve :",length(which(centroid_env$Color_differentiel == "Reserve")),"\n")
          )
          
      # Officiel & Reserve / non echantillonne -----
          
        tbl_resum_o2rne <- as.data.frame(matrix(nrow=0,ncol=2)) ; colnames(tbl_resum_o2rne) <- c("Color_differentiel1","variable") # initialisation de la table resume finale
        
        for(i in 1:30){
          # detection des SSU non echantillonne
          SSU_non_echan <- centroid_env[-which(centroid_env$ID %in% ID_SSU_offres),] %>%
            sample_n(size = length(ID_SSU_offres),replace=F) # tirage aleatoire
          
          # creation de la table regroupant tout les SSU officiels/reserve + un tirage des SSU non echantillonnes de taille equivalente
          SSU_env_o2rne_ta <- rbind(SSU_non_echan,centroid_env[which(centroid_env$ID %in% ID_SSU_offres),])
          
          # formation du tableau croise dynamique
          SSU_hab_o2rne_ta <- reshape2::melt(SSU_env_o2rne_ta[,-19],id.vars = c("ID","Color_differentiel1")) # rm de la colonne a 3 niveaux
          
          SSU_hab_o2rne_ta1 <- SSU_hab_o2rne_ta %>%
            group_by(variable,Color_differentiel1) %>%
            mutate(total = sum(as.numeric(value)))
          
          SSU_hab_o2rne_ta2 <- unique(SSU_hab_o2rne_ta1[,c("Color_differentiel1","variable","total")])
          colnames(SSU_hab_o2rne_ta2)[3] <- paste0("Tirage_",i)
          
          tbl_resum_o2rne <- merge(tbl_resum_o2rne,SSU_hab_o2rne_ta2,all=T) # jointure w/ table resume finale
          
          
          cat(i,"\n")
        }
        
        # calcul de la surface moyenne apres tirage aleatoire des SSU non echantillonnee
        tbl_resum_o2rne$Tirage_mean <- apply(tbl_resum_o2rne[,grep("Tirage",colnames(tbl_resum_o2rne))],1,mean)
        
        # plot ----
        ggplot(tbl_resum_o2rne,aes(fill = Color_differentiel1,y=Tirage_mean,x=variable)) +
          geom_bar(position = "dodge", stat = "identity") +
          xlab("Habitats CLC niveau 2") + ylab("Surface en ha") +
          labs(fill = "Statut des SSU") +
          ggtitle(paste("Surface totale des habitats CLC niv2 selon un buffer de 500m au centroide des SSU\nSSU non echantillonnés tirés aleatoirement",i,"fois"),
                  subtitle = paste(" Nb SSU non echantillonne :",nrow(SSU_non_echan),"\n",
                                   "Nb SSU officiel + reserve :",length(which(centroid_env$Color_differentiel1 == "Officiel_Reserve")))
          )
      

  
       
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
              
    
    
    
    