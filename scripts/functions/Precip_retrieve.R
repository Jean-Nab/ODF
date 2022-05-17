##############################################################################################################
# This script allows to extract the precipitations sum (in mm) at sites for given dates.
# Those sites have to be in metropolitan France (and Corsica)
#
# It returns the data frame of each day and site (given in arguments) with a new column for the precipitations
# of the day (rr_day), and as many columns as the number of days before the survey specified in arguments 
# (e.g. precipitations of the day before: rr_day_before1, precipitations two days before: rr_day_before2, etc).
#
##############################################################################################################
#
# Meteorological data (from E-OBS) as to be download at this link : 
# https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php#datafiles
# Choose either  : Ensemble Mean RR from the 0.25 deg. regular grid 
# or the 0.1 deg. regular grid (more precise but heavier)
#
##############################################################################################################
#
# Arguments : 
#
# path_to_meteo_nc : path to access to the meteo .nc file
# 
# tableDaysSites : a table with (minimum):
# -> a column named "Latitude" with the site latitude in WGS84
# -> a column named "Longitude" with the site longitude in WGS84
# -> a column with the sites unique IDs
# -> a column named "Date" with the date of the survey for each sites
#       Dates should be formatted this way : "2020-05-26" if characters
#       They might also be c("POSIXct","POSIXt") objects
#
# nbrPreviousDays : number of previous days for which T° should be returned
# (for example "3", to obtain the temperature for the survey date but for also the three previous days)

##############################################################################################################

# path_to_meteo_nc <- "C:/Users/jean.nabias/Desktop/Ressource QGis/weather/rr_ens_spread_0.1deg_reg_v25.0e.nc"

extract_precipitations <- function(path_to_meteo_nc,tableDaysSites,nbrPreviousDays, prefixe_fichier){
  
  ####Preliminary####
  
  #Install and open required packages
  load <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  } 
  packages <- c("ncdf4","sf","stringr","lubridate","ggplot2","tidyr","data.table")
  load(packages)
  
  # Conversion tabl
    #lon/lat
    tableDaysSites <- st_transform(st_as_sf(tableDaysSites, coords = c("X_barycentre_L93","Y_barycentre_L93"), crs = 2154,remove = F), crs = 4326)
    wxc <- st_coordinates(tableDaysSites)
    tableDaysSites$Longitude <- wxc[,1]
    tableDaysSites$Latitude <- wxc[,2]
    tableDaysSites <- st_drop_geometry(tableDaysSites)
    
    #date
    if(length(grep("^date$",names(tableDaysSites))) == 0){
      tableDaysSites$Date <- paste0(tableDaysSites$Annee,"-",tableDaysSites$Mois,"-",tableDaysSites$Jour)
      tableDaysSites$Date <- as.POSIXct(tableDaysSites$Date)
    }else{ # cas STOC
      names(tableDaysSites)[grep("^date$",names(tableDaysSites))] <- "Date"
    }
    
    
  
  
  #creation of a site layer
  sites_SHP <- st_as_sf(tableDaysSites, coords = c("Longitude","Latitude"), crs = 4326)
  sites_SHP <- st_transform(sites_SHP,crs=2154)
  
  if (all(class(tableDaysSites[,which(colnames(tableDaysSites)=="Date")])!=c("POSIXct","POSIXt"))){
    tableDaysSites[,which(colnames(tableDaysSites)=="Date")] <- parse_date_time(
      as.character(tableDaysSites[,which(colnames(tableDaysSites)=="Date")]),"%Y-%m-%d", tz="Europe/Paris")
  }
  
  ####I. Meteo data formatting####
  
  #Open meteo file
  rr.ncdf <- nc_open(path_to_meteo_nc)
  
  #This file has three dimensions : longitude, latitude, date
  #Extraction of the values taken by each of these dimensions
  lon <- rr.ncdf$dim$longitude$vals
  lat <- rr.ncdf$dim$latitude$vals
  date <- rr.ncdf$dim$time$vals
  date <- as_date(date,origin="1950-01-01") #conversion in dates object
  
  #lim date : from 1980-01-01 to today
  start_date = which(date=="1980-01-01")
  end_date = length(date)
  len_date <- end_date-start_date+1
  
  #lim longitude (for France)
  start_lon <- which(lon+6==min(abs(lon+6)))
  stop_lon <- which(lon-10==min(abs(lon-10)))
  len_lon <- stop_lon-start_lon+1
  
  #lim longitude (for France)
  start_lat <- which(lat-41==min(abs(lat-41)))
  stop_lat <- which(lat-52==min(abs(lat-52)))
  len_lat <- stop_lat-start_lat+1
  
  #data extraction from the .nc file into an array
  rr_array <- ncvar_get(rr.ncdf,"rr",start=c(start_lon,start_lat,start_date),count=c(len_lon,len_lat,len_date))
  dimnames(rr_array) <- list(c(lon[start_lon:stop_lon]),c(lat[start_lat:stop_lat]),c(as.character(date[start_date:end_date])))
  
  print("1. Data formatting OK")
  
  ####II. Link between sites and meteo points####
  
  ####____A. Detect longitude * latitude for which data are available####
  
  #Take one date not to recent (data has to exist) but not old either
  #Here arbitrarily : 2019-05-15
  #And extract temperature for each longitude * latitude at this date (array with 2 dimensions)
  
  date_map = which(date=="2019-05-15")
  rr_link_array <- ncvar_get(rr.ncdf,"rr",start=c(start_lon,start_lat,date_map),count=c(len_lon,len_lat,1)) #extraction
  dimnames(rr_link_array) <- list(c(lon[start_lon:stop_lon]),c(lat[start_lat:stop_lat])) #longitude and latitude coord
  rr_link_arrayDT <- as.data.frame(rr_link_array) #array to data.frame
  rr_link_arrayDT$lon <- rownames(rr_link_arrayDT) #longitude in a column
  rr_link_arrayDT <- gather(rr_link_arrayDT,key="lat",value="rr",c(1:((dim(rr_link_arrayDT)[2])-1))) #reshape the data frame in a data frame with a column longitude, a column latitude and a column with the mean temperature
  rr_link_arrayDT <- rr_link_arrayDT[which(!(is.na(rr_link_arrayDT$rr))),] #remove longitude * latitude with NA
  rr_link_arrayDT$coord <- paste0(rr_link_arrayDT$lon,"_",rr_link_arrayDT$lat) #create a column with the concatenation of longitude and latitude
  
  ####____B. Create a shapefile with the points of the meteo data####
  
  coord_meteo_SHP <- st_as_sf(rr_link_arrayDT, coords = c("lon","lat"), crs = 4326)
  coord_meteo_SHP <- st_transform(coord_meteo_SHP,crs=2154)
  
  ####____C. create the link between meteo data and sites####
  
  #Creation of a new column in the shapefile with the sites that indicate, for each site, the closest point
  #of the meteo data
  
  sites_SHP$corresp_EObs <- coord_meteo_SHP$coord[st_nearest_feature(sites_SHP,coord_meteo_SHP)]
  
  print("2. Link between sites and meteo points OK")
  
  ####III. Temperature and anomaly extraction for each night####
  
  tableDaysSites$corresp_EObs <- as.data.frame(sites_SHP$corresp_EObs)
  
  tableDaysSites$rr_day <- rep(NA,dim(tableDaysSites)[1])
  
  if (nbrPreviousDays !=0){
    for (i in 1:nbrPreviousDays){
      tableDaysSites[,dim(tableDaysSites)[2]+1] <- rep(NA,dim(tableDaysSites)[1])
      colnames(tableDaysSites)[dim(tableDaysSites)[2]] <- paste0("rr_day_before",i)
    }
  }
  
  for (index in c(1:dim(tableDaysSites)[1])){
    
    date_day <- tableDaysSites[index,which(colnames(tableDaysSites)=="Date")]
    
    if (as.character(date_day) %in% (dimnames(rr_array)[[3]])==T){
      
      coord_meteo <- sites_SHP$corresp_EObs[index]
      coord_lon <- str_split(coord_meteo,"_",simplify = T)[1]
      coord_lat <- str_split(coord_meteo,"_",simplify = T)[2]
      
      #raw temperature
      tableDaysSites$rr_day[index] <- rr_array[coord_lon,coord_lat,as.character(date_day)]
      
      if (nbrPreviousDays!=0){
        for (i in (1:nbrPreviousDays)){
          tableDaysSites[index,which(colnames(tableDaysSites)==paste0("rr_day_before",i))] <- rr_array[coord_lon,coord_lat,as.character(date_day-days(i))]
        }
      }
      
    }
    
    cat(paste0(round(index/(dim(tableDaysSites)[1])*100,digits=1),"%\r"))
  }
  
  print("3. Precipitations extraction for each night OK")
  
  names(tableDaysSites)[ncol(tableDaysSites)] <- "Sites_precipiationDay"
  write.csv(tableDaysSites[,c("ID_liste","Sites_precipiationDay")],
            file = paste0("C:/git/ODF/output/function_output/",prefixe_fichier,"Precip_envEPOC.csv"), row.names = F)
}