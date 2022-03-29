
# Fonctions basees selon le modele des fonctions d'extraction de Yves Bas (https://github.com/cesco-lab/Vigie-Chiro_scripts/tree/master/functions/extractGI) 


  #extract_CLC(dsnTable = "C:/git/epoc/DS.v2/epoc_barycentre_liste_density_add.csv",
  #            names_coord = c("X_barycentre_L93","Y_barycentre_L93"),
  #            buffer_medium = 500,
  #            buffer_large = 5000,
  #            dsnRasterCLC = "C:/Users/Travail/Desktop/Ressource QGis/france/CLC2018_raster100m_L93.tif")
 
  
  
  # dsnTable = "C:/git/epoc/DS.v2/epoc_barycentre_liste_density_add.csv"
  # names_coord = c("X_barycentre_L93","Y_barycentre_L93")
  # buffer_medium = 500
  # buffer_large = 5000
  # dsnRasterCLC = "C:/Users/Travail/Desktop/Ressource QGis/france/CLC2018_raster100m_L93.tif"



extract_CLC = function(dsnTable, names_coord, buffer_medium, buffer_large, dsnRasterCLC,
                       prefixe_fichier)
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
    
    CLC <- raster(dsnRasterCLC)
    
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
    #coordinates(Point) <- Coords
    #proj4string(Point) <- CRS("+init=epsg:2154") # L93 /// CRS("+init=epsg:4326") # WGS 84
    
    
    # Si coords en wgs84
    ## Point_L93=spTransform(Point,proj4string(CLC))
    ##
    ##  # Retrait des points en dehors du raster
    ##    PointL_L93=spTransform(PointL,proj4string(Hab))
    ##    PointL_L93=crop(PointL_L93,Hab)
    ##    PointL=subset(PointL,PointL$id %in% PointL_L93$id)
    
    cat("\t----\tFin du chargement des donnees\t----\n")
  
  # EXTRACTION ----
    # Buffer medium ----
      cat(paste("\t----\tDebut extraction des donnees CLC selon le buffer :",BuffMed,"m\t----\n"))
    
      vec.iteration_med <- seq(from = 1, to = length(PointBuffMed$ID_extract), by = 1000) # argument suplementaire ? / add condition : NA/all -> tous
      
      CLC_hab_wide_BM <- as.data.frame(matrix(nrow=0,ncol = 1))
      colnames(CLC_hab_wide_BM) <- "id"
      
      # initialisation de la barre de progression
      pb_m <- progress::progress_bar$new(total = length(vec.iteration_med),
                                       format = "Extraction des donnees CLC buffer medium [:bar] :percent    Tps ecoule = :elapsedfull",
                                       clear = F)
     
      pb_m$tick(0)
      Sys.sleep(1/20)
      
      
      # boucle de 305 en 305: ----
      for(i in vec.iteration_med){
        
        PointBuffMed.tmp <- PointBuffMed[i:min((i+999),length(PointBuffMed$ID_extract)),] # si add condition --> condition - 1
        
        CLC_hab.tmp <- raster::extract(x = CLC, 
                               y = PointBuffMed.tmp,
                               df=T,
                               along=T)
        
        
        # formatage des donnees
        colnames(CLC_hab.tmp)[1] <- "ID_extract"
        CLC_hab.tmp$ID_extract <- CLC_hab.tmp$ID_extract + (i - 1)
        colnames(CLC_hab.tmp)[2] <- "CLCniv3"
        # extract selon CLC niv3 --> recuperation du niveau 2
        CLC_hab.tmp$CLCniv2 <- trunc(CLC_hab.tmp$CLCniv3/10)
        CLC_hab.tmp$CLCniv1 <- trunc(CLC_hab.tmp$CLCniv3/100)
        
        
        # gestion dtf des resultats preliminaires -----
        
        CLC_hab.tmp <- left_join(CLC_hab.tmp,st_drop_geometry(PointBuffMed.tmp[,c("id","ID_extract")]), by = "ID_extract")
        
        CLCniv3_hab_wide.tmp <- CLC_hab.tmp %>% reshape2::dcast(ID_extract + id ~ CLCniv3, fun.aggregate = length, value.var = "id")
        CLCniv3_hab_wide.tmp <- CLCniv3_hab_wide.tmp %>%
          dplyr::select(-contains("0"))  # retrait des habitats 0 ==> zone non couverte par le raster
            w <- as.data.frame(t(apply(X = CLCniv3_hab_wide.tmp[,grep("[0-9]",colnames(CLCniv3_hab_wide.tmp))],
                       MARGIN = 1,
                       FUN = function(X){t(X <- X/sum(X)*100)})))
            names(w) <- names(CLCniv3_hab_wide.tmp)[grep("[0-9]",colnames(CLCniv3_hab_wide.tmp))]
            CLCniv3_hab_wide.tmp[,grep("[0-9]",colnames(CLCniv3_hab_wide.tmp))] <- w
            
        
        CLCniv2_hab_wide.tmp <- CLC_hab.tmp %>% reshape2::dcast(ID_extract + id ~ CLCniv2, fun.aggregate = length, value.var = "id")
        CLCniv2_hab_wide.tmp <- CLCniv2_hab_wide.tmp %>%
          dplyr::select(-contains("0")) 
            w <- as.data.frame(t(apply(X = CLCniv2_hab_wide.tmp[,grep("[0-9]",colnames(CLCniv2_hab_wide.tmp))],
                                       MARGIN = 1,
                                       FUN = function(X){t(X <- X/sum(X)*100)})))
            names(w) <- names(CLCniv2_hab_wide.tmp)[grep("[0-9]",colnames(CLCniv2_hab_wide.tmp))]
            CLCniv2_hab_wide.tmp[,grep("[0-9]",colnames(CLCniv2_hab_wide.tmp))] <- w
        
        
        CLCniv1_hab_wide.tmp <- CLC_hab.tmp %>% reshape2::dcast(ID_extract + id ~ CLCniv1, fun.aggregate = length, value.var = "id")
        CLCniv1_hab_wide.tmp <- CLCniv1_hab_wide.tmp %>%
          dplyr::select(-contains("0")) 
            w <- as.data.frame(t(apply(X = CLCniv1_hab_wide.tmp[,grep("[0-9]",colnames(CLCniv1_hab_wide.tmp))],
                                       MARGIN = 1,
                                       FUN = function(X){t(X <- X/sum(X)*100)})))
            names(w) <- names(CLCniv1_hab_wide.tmp)[grep("[0-9]",colnames(CLCniv1_hab_wide.tmp))]
            CLCniv1_hab_wide.tmp[,grep("[0-9]",colnames(CLCniv1_hab_wide.tmp))] <- w
        
        
        
        CLC_hab_wide_BM.tmp <- left_join(CLCniv1_hab_wide.tmp, left_join(CLCniv2_hab_wide.tmp,CLCniv3_hab_wide.tmp,by = c("ID_extract", "id")),by = c("ID_extract", "id"))
        
        CLC_hab_wide_BM <- merge(CLC_hab_wide_BM,CLC_hab_wide_BM.tmp,all=T)
        CLC_hab_wide_BM[is.na(CLC_hab_wide_BM)] <- 0 # corrections des NA a la suite du merge (habitats absent de certaine listes)
        
        # actualisation de la progression
        pb_m$tick()
        Sys.sleep(1 / 100)
      }
      
      ### CLC_hab_wide_BM <- CLC_hab_wide_BM[,-grep("FALSE|TRUE|,",colnames(CLC_hab_wide_BM))] # retrait d'un artefact de merge
      
      ### correction des NAs (= habitats absent lors de l'extract du point)
      ##CLC_hab_wide_BM[is.na(CLC_hab_wide_BM)] <- 0
      ##
      ##CLC_hab_wide_BM <- CLC_hab_wide_BM %>%
      ##  dplyr::select(-contains("0"))  # retrait des habitats 0 ==> zone non couverte par le raster
      ##
      ##
      ### calcul des proportions d'habitats 
      ##CLC_hab_wide_BM[,grep("[0-9]",colnames(CLC_hab_wide_BM))] <- ( CLC_hab_wide_BM[,grep("[0-9]",colnames(CLC_hab_wide_BM))] / 
      ##                                                           rowSums(CLC_hab_wide_BM[,grep("[0-9]",colnames(CLC_hab_wide_BM))])
      ##) * 100
    
      colnames(CLC_hab_wide_BM) <- gsub("([0-9]{1,})","SpCLCm_\\1",colnames(CLC_hab_wide_BM))    
      
      
    # Buffer large -------
      
      cat(paste("\n\t----\tDebut extraction des donnees CLC selon le buffer :",BuffLarg,"m\t----\n"))
      
      vec.iteration_larg <- seq(from = 1, to = length(PointBuffLarg$ID_extract), by = 400) # argument suplementaire ? / add condition : NA/all -> tous
      
      CLC_hab_wide_BL <- as.data.frame(matrix(nrow=0,ncol = 1))
      colnames(CLC_hab_wide_BL) <- "id"
      
      # initialisation de la barre de progression
      pb_l <- progress::progress_bar$new(total = length(vec.iteration_larg),
                                       format = "Extraction des donnees CLC buffer large [:bar] :percent    Tps ecoule = :elapsedfull",
                                       clear = F)
      
      pb_l$tick(0)
      Sys.sleep(1/20)
      
      
      # boucle de 305 en 305: ----
      for(i in vec.iteration_larg){
        
        PointBuffLarg.tmp <- PointBuffLarg[i:min((i+399),length(PointBuffLarg$ID_extract)),] # si add condition --> condition - 1
        
        
        CLC_hab.tmp <- raster::extract(x = CLC, 
                                       y = PointBuffLarg.tmp,
                                       df=T,
                                       along=T)
        
        
        # formatage des donnees
        colnames(CLC_hab.tmp)[1] <- "ID_extract"
        CLC_hab.tmp$ID_extract <- CLC_hab.tmp$ID_extract + (i - 1)
        colnames(CLC_hab.tmp)[2] <- "CLCniv3"
        # extract selon CLC niv3 --> recuperation du niveau 2
        CLC_hab.tmp$CLCniv2 <- trunc(CLC_hab.tmp$CLCniv3/10)
        CLC_hab.tmp$CLCniv1 <- trunc(CLC_hab.tmp$CLCniv3/100)
        
        
        # gestion dtf des resultats preliminaires -----
        
        CLC_hab.tmp <- left_join(CLC_hab.tmp,st_drop_geometry(PointBuffLarg.tmp[,c("id","ID_extract")]), by = "ID_extract")
        
        CLCniv3_hab_wide.tmp <- CLC_hab.tmp %>% reshape2::dcast(ID_extract + id ~ CLCniv3, fun.aggregate = length, value.var = "id")
        CLCniv3_hab_wide.tmp <- CLCniv3_hab_wide.tmp %>%
          dplyr::select(-contains("0"))  # retrait des habitats 0 ==> zone non couverte par le raster
        w <- as.data.frame(t(apply(X = CLCniv3_hab_wide.tmp[,grep("[0-9]",colnames(CLCniv3_hab_wide.tmp))],
                                   MARGIN = 1,
                                   FUN = function(X){t(X <- X/sum(X)*100)})))
        names(w) <- names(CLCniv3_hab_wide.tmp)[grep("[0-9]",colnames(CLCniv3_hab_wide.tmp))]
        CLCniv3_hab_wide.tmp[,grep("[0-9]",colnames(CLCniv3_hab_wide.tmp))] <- w
        
        
        CLCniv2_hab_wide.tmp <- CLC_hab.tmp %>% reshape2::dcast(ID_extract + id ~ CLCniv2, fun.aggregate = length, value.var = "id")
        CLCniv2_hab_wide.tmp <- CLCniv2_hab_wide.tmp %>%
          dplyr::select(-contains("0")) 
        w <- as.data.frame(t(apply(X = CLCniv2_hab_wide.tmp[,grep("[0-9]",colnames(CLCniv2_hab_wide.tmp))],
                                   MARGIN = 1,
                                   FUN = function(X){t(X <- X/sum(X)*100)})))
        names(w) <- names(CLCniv2_hab_wide.tmp)[grep("[0-9]",colnames(CLCniv2_hab_wide.tmp))]
        CLCniv2_hab_wide.tmp[,grep("[0-9]",colnames(CLCniv2_hab_wide.tmp))] <- w
        
        
        CLCniv1_hab_wide.tmp <- CLC_hab.tmp %>% reshape2::dcast(ID_extract + id ~ CLCniv1, fun.aggregate = length, value.var = "id")
        CLCniv1_hab_wide.tmp <- CLCniv1_hab_wide.tmp %>%
          dplyr::select(-contains("0")) 
        w <- as.data.frame(t(apply(X = CLCniv1_hab_wide.tmp[,grep("[0-9]",colnames(CLCniv1_hab_wide.tmp))],
                                   MARGIN = 1,
                                   FUN = function(X){t(X <- X/sum(X)*100)})))
        names(w) <- names(CLCniv1_hab_wide.tmp)[grep("[0-9]",colnames(CLCniv1_hab_wide.tmp))]
        CLCniv1_hab_wide.tmp[,grep("[0-9]",colnames(CLCniv1_hab_wide.tmp))] <- w
        
        
        CLC_hab_wide_BL.tmp <- left_join(CLCniv1_hab_wide.tmp, left_join(CLCniv2_hab_wide.tmp,CLCniv3_hab_wide.tmp,by = c("ID_extract", "id")),by = c("ID_extract", "id"))
        
        CLC_hab_wide_BL <- merge(CLC_hab_wide_BL,CLC_hab_wide_BL.tmp,all=T)
        
        
        # actualisation de la progression
        pb_l$tick()
        Sys.sleep(1 / 100)
        
      }
      
      
      ### correction des NAs (= habitats absent lors de l'extract du point)
      ##CLC_hab_wide_BL[is.na(CLC_hab_wide_BL)] <- 0
      ##
      ##CLC_hab_wide_BL <- CLC_hab_wide_BL %>%
      ##  dplyr::select(-contains("0"))  # retrait des habitats 0 ==> zone non couverte par le raster
      ##
      ##
      ### calcul des proportions d'habitats 
      ##CLC_hab_wide_BL[,grep("[0-9]",colnames(CLC_hab_wide_BL))] <- ( CLC_hab_wide_BL[,grep("[0-9]",colnames(CLC_hab_wide_BL))] / 
      ##                                                                 rowSums(CLC_hab_wide_BL[,grep("[0-9]",colnames(CLC_hab_wide_BL))])
      ##) * 100

      
      colnames(CLC_hab_wide_BL) <- gsub("([0-9]{1,})","SpCLCl_\\1",colnames(CLC_hab_wide_BL))  
    

  # jointure des donnes CLC selon les 2 buffers
    CLC_hab_BMBL <- dplyr::left_join(CLC_hab_wide_BM,CLC_hab_wide_BL, by = c("id", "ID_extract"))
    
    CLC_hab_BMBL <- dplyr::left_join(Point[,c("id","ID_extract","ID_liste",names_coord)],
                                     CLC_hab_BMBL , by = c("id", "ID_extract"))
  
    CLC_hab_BMBL <- dplyr::left_join(CLC_hab_BMBL, Point[,c("id",Coords,"ID_liste")]) # add d'informations sur les listes
    
  # sauvegarde sur disque
    save(CLC_hab_BMBL, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"CLC_envEPOC.RData")) # securite
    write.csv(CLC_hab_BMBL, file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"CLC_envEPOC.csv"), row.names = F)
  
  cat("\n\n\t-------\tExtract CLC done - check /function_output \t-------\n\n")
  
  }














