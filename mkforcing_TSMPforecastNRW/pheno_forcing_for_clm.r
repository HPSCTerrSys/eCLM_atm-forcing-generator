## Program to prepare atmospheric forcing from CLM ENS2 data
## 
## to load R at JUWELS just type:
## module load Intel ParaStationMPI R
## module load Silo/4.10.2
##
## Stefan Poll aug 2021
## source("pheno_forcing_for_clm.r")

rm(list=ls())

#### settings ####
require(ncdf4)
require(lubridate) # for date manipulation


# switch
l_calc = TRUE
l_100m = TRUE #FALSE
l_10m=TRUE

#directories
path_forc  = "/p/largedata/slts/shared_data/imod_TSMP-NRW_IBG3/o.data.tsmpm_nrw/"
lonl = 300
latl = 300

# setting specific ens
if (l_100m){
#directories
path_out = "/p/project/cslts/poll1/data/pheno/kleinaltendorf_100m/clm/atm_forcing/"
path_clm  = "/p/project/cslts/poll1/data/pheno/kleinaltendorf_100m/clm/" 
}

if (l_10m){
#directories
path_out2 = "/p/project/cslts/poll1/data/pheno/kleinaltendorf_10m/clm_new/atm_forcing/"
path_clm2  = "/p/project/cslts/poll1/data/pheno/kleinaltendorf_10m/clm_new/" 
}

year <- 2019
monthv <- c(2,3,4,5,6)
#monthv <- c(1,2,3,4,5,6)
#monthv <- c(7,8,9,10,11,12)

#### function #####



#### start program ####

# 
 rad2deg=180/pi
 deg2rad=pi/180

#

for (month in monthv) {  
  
monthstr <- sprintf("%02i", month)
mdate <- as.POSIXct(paste0(year,monthstr,"01"),format="%Y%m%d")

day_no = as.numeric(days_in_month(mdate))
hrs_sta = 0
hrs_end = 23 

  print("#####################")
  print(paste0("Start programm: ",paste0(year,monthstr)))
  print("#####################")
  
  datefail <- ""

if (l_calc) {
  
  icount <- 1 # reset counter
  
  hrl  <- (hrs_end-hrs_sta+1)*day_no
  tbotfld  <- windfld  <- tprecfld <- fsdsfld  <- fldsfld  <- thbotfld <- zbotfld <- qbotfld <- array(NA,c(lonl,latl,hrl))

for (iday in 1:day_no) {

 print(paste0("Process day: ",iday))

 daystr = sprintf("%02i", iday)
  
 for (ihr in 0:23) {

  hrstr = sprintf("%05i", ihr*3600) 
 
# read data
# change directory for the 0th hour
  if (ihr==0){
#     fname <-  paste0(path_forc,"/TSMPForecastNRW",as.character(mdate+(iday-1-1)*3600*24),"/clm/clmoas.clm2.h0.",year,"-",monthstr,"-",daystr,"-",hrstr,".nc")
    fname <-  paste0(path_forc,"/TSMPForecastNRW",ymd(mdate)+(iday-1-1),"/clm/clmoas.clm2.h0.",year,"-",monthstr,"-",daystr,"-",hrstr,".nc")
  }else{
    fname <-  paste0(path_forc,"/TSMPForecastNRW",ymd(mdate)+(iday-1),"/clm/clmoas.clm2.h0.",year,"-",monthstr,"-",daystr,"-",hrstr,".nc")
  } # if ihr
# 
    # workaround if file does not exist
    if (!(file.exists(fname))){
      print(paste0("File does not exit for day ",iday," hour ",ihr,""))
      datefail <- c(datefail,paste0(" ",year,monthstr,daystr, sprintf("%02i", ihr)))
      dshift <- 0 
      # take same forecast hour from different run
      while ((!(file.exists(fname))) & (!(dshift>14))){
        fname <- paste0(path_forc,"/TSMPForecastNRW",ymd(mdate)+(iday-1-1-dshift),"/clm/clmoas.clm2.h0.",year,"-",monthstr,"-",daystr,"-",hrstr,".nc")
        dshift <- dshift + 1 
      }
      # to avoid infinity loop take a "similar" hour to complete dataset
      if (dshift>14){
        print(paste0("!!!!!! File does not exit at all for day ",iday," hour ",ihr,""))
        print(paste0("Take file instead: ",fname_h))
        fname <- fname_h
      }
      # tmp variable for missing 
      fname_h <- fname
    }

   ff     <- nc_open(fname)

#   } #l_ens1
# 
#  if (ihr==hrs_sta & iday==1 & iens==1) { 
  if (iday==1) { 
    lon <- ncvar_get(ff,varid="longxy")
    lat <- ncvar_get(ff,varid="latixy")
  } # if
#
  tbot  <- ncvar_get(ff,varid="TBOT")
  wind  <- ncvar_get(ff,varid="WIND")
  qbot  <- ncvar_get(ff,varid="QBOT")
#   tprec <- ncvar_get(ff,varid="RAIN")
  tprec <- ncvar_get(ff,varid="QFLX_RAIN_GRND")
  fsds  <- ncvar_get(ff,varid="FSDS")
  flds  <- ncvar_get(ff,varid="FLDS")
  thbot <- ncvar_get(ff,varid="THBOT")
#  zbot  <- ncvar_get(ff,varid="ZBOT")
#
# ESTIMATION of lowermost layer height! 
   zbot <- array(9.5,dim(tbot))
#
# filevarattdef(fout, "LONGXY"  , longxy)
# filevarattdef(fout, "LATIXY"  , latixy)
# filevarattdef(fout, "TBOT"  , tbot)
# filevarattdef(fout, "WIND"  , wind)
# filevarattdef(fout, "QBOT"  , qbot)
# filevarattdef(fout, "PRECTmms"  , prectmms)
# filevarattdef(fout, "FSDS"  , fsds)
# filevarattdef(fout, "FLDS"  , flds)
# filevarattdef(fout, "PSRF"  , psrf)
#
  nc_close(ff)
# 
#####
# check if data have a time dimension
 if (length(dim(tbot)<3)){
   timeslt <- 1
 }else{
   timeslt <- dim(tbot)[3]
 }
# 
 tbotfld[,,icount:(icount+timeslt-1)] <- tbot
 windfld[,,icount:(icount+timeslt-1)] <- wind
 qbotfld[,,icount:(icount+timeslt-1)] <- qbot
 tprecfld[,,icount:(icount+timeslt-1)]<- tprec
 fsdsfld[,,icount:(icount+timeslt-1)] <- fsds
 fldsfld[,,icount:(icount+timeslt-1)] <- flds
 thbotfld[,,icount:(icount+timeslt-1)]<- thbot
 zbotfld[,,icount:(icount+timeslt-1)] <- zbot
# 
 icount <- icount + timeslt
#####
#
 } # for ihr
} # for iday

### save file to 

#  
  subdir <- paste0("o.data")
  ifelse(!dir.exists(file.path(path_out, subdir)), dir.create(file.path(path_out, subdir)), FALSE)

  fnameout <- paste0(path_out,"/",subdir,"/",year,"-",monthstr,".nc")

  print(paste0("Save nc file: ",fnameout))

  # definition of dimensions
  dimlon  <- ncdim_def("lon","",1:dim(tbotfld)[1])
  dimlat  <- ncdim_def("lat","",1:dim(tbotfld)[2])
  dimtime <- ncdim_def("time","",1:dim(tbotfld)[3], unlim=TRUE)

  # Variables
  mv      <- 9.96920996838687e+36
  var01   <- ncvar_def("LONGXY"     , "degrees"  , list(dimlon, dimlat), prec="float",longname="lon")
  var02   <- ncvar_def("LATIXY"     , "degrees"  , list(dimlon, dimlat), prec="float",longname="lat") 
  var1    <- ncvar_def("TBOT"       , "K"        , list(dimlon, dimlat, dimtime), mv, prec="float",longname="temperature at lowest atm level")
  var12   <- ncvar_def("THBOT"      , "K"        , list(dimlon, dimlat, dimtime), mv, prec="float",longname="pot. temperature at lowest atm level")
  var2    <- ncvar_def("WIND"       , "m/s"      , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Wind speed at lowest atm level")
  var23   <- ncvar_def("QBOT"       , "kg/kg"    , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Sp. humidity at lowest atm level")
  var3    <- ncvar_def("PRECTmms"   , "mm/s"     , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Surface precipitation")
  var4    <- ncvar_def("FSDS"       , "W/m2"     , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Incident Solar radiation")
  var5    <- ncvar_def("FLDS"       , "W/m2"     , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Incident Long Wave radiation")
  var6    <- ncvar_def("ZBOT"       , "m"        , list(dimlon, dimlat, dimtime), mv, prec="float",longname="reference height")

  nc      <- nc_create(fnameout,list(var01,var02,var1,var12,var2,var23,var3,var4,var5,var6))

  # add source info metadata to file
  ncatt_put( nc, 0, "content", "Atmospheric Forcing for CLM 3.5 extracted from TSMP NRW")
  ncatt_put( nc, 0, "creator", "s.poll")
  ncatt_put( nc, 0, "date_fail", datefail)
  ncatt_put( nc, 0, "date_of_creation", paste0(Sys.time()))

  # put value to var 
  ncvar_put(nc,var01,lon)
  ncvar_put(nc,var02,lat)
  ncvar_put(nc,var1,tbotfld)
  ncvar_put(nc,var12,thbotfld)
  ncvar_put(nc,var2,windfld)
  ncvar_put(nc,var23,qbotfld)
  ncvar_put(nc,var3,tprecfld)
  ncvar_put(nc,var4,fsdsfld)
  ncvar_put(nc,var5,fldsfld)
  ncvar_put(nc,var6,zbotfld)

  nc_close(nc)

  print("nc file written")

}else { # if l_calc

   fnameout <- paste0(path_out,"/",subdir,"/",year,"-",monthstr,".nc")

} # if l_calc
  
## Interpolate NRW CLM data to CKA CLM
  
if (l_100m){

   print(paste0("Do nearest neighbor interpolation for 100m"))
#   library(raster)
  
  filename_clm <- paste0(path_clm,"/griddata_0250x0250.nc")
  nc  <- nc_open(filename_clm)
  lon_clm <- ncvar_get(nc,"LONGXY")
  lat_clm <- ncvar_get(nc,"LATIXY")
  nc_close(nc)

  # find nearest neighbor
  indnn <- array(NA,c(dim(lon_clm),2))
  for (ilon in 1:dim(lon_clm)[1]){
   for (ilat in 1:dim(lon_clm)[2]){
     indnn[ilon,ilat,] <- which((abs(lon_clm[ilon,ilat]-lon)+abs(lat_clm[ilon,ilat]-lat))==min(abs(lon_clm[ilon,ilat]-lon)+abs(lat_clm[ilon,ilat]-lat)),arr.ind=TRUE)
   }   
  }
  # image(thbotfld[191:244,140:191]) # test image
  
  # do nearest neighbor interpolation
  tbotint <- thbotint <- windint <- qbotint <- tprecint <- fsdsint <- fldsint <- zbotint <- array(NA,c(dim(lon_clm),dim(tbotfld)[3]))
  for (ilon in 1:dim(lon_clm)[1]){
   for (ilat in 1:dim(lon_clm)[2]){
    tbotint[ilon,ilat,] <- tbotfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    thbotint[ilon,ilat,] <- thbotfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    windint[ilon,ilat,] <- windfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    qbotint[ilon,ilat,] <- qbotfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    tprecint[ilon,ilat,] <- tprecfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    fsdsint[ilon,ilat,] <- fsdsfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    fldsint[ilon,ilat,] <- fldsfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    zbotint[ilon,ilat,] <- zbotfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
   }   
  }

  subdir <- paste0("m.data")
  ifelse(!dir.exists(file.path(path_out, subdir)), dir.create(file.path(path_out, subdir)), FALSE)

  fnameout <- paste0(path_out,"/",subdir,"/",year,"-",monthstr,".nc")

  print(paste0("Save nc file: ",fnameout))

  # definition of dimensions
  dimlon  <- ncdim_def("lon","",1:dim(lon_clm)[1])
  dimlat  <- ncdim_def("lat","",1:dim(lon_clm)[2])
  dimtime <- ncdim_def("time","",1:dim(tbotfld)[3], unlim=TRUE)

  # Variables
  mv      <- 9.96920996838687e+36
  var01   <- ncvar_def("LONGXY"     , "degrees"  , list(dimlon, dimlat), prec="float",longname="lon")
  var02   <- ncvar_def("LATIXY"     , "degrees"  , list(dimlon, dimlat), prec="float",longname="lat") 
  var1    <- ncvar_def("TBOT"       , "K"        , list(dimlon, dimlat, dimtime), mv, prec="float",longname="temperature at lowest atm level")
  var12   <- ncvar_def("THBOT"      , "K"        , list(dimlon, dimlat, dimtime), mv, prec="float",longname="pot. temperature at lowest atm level")
  var2    <- ncvar_def("WIND"       , "m/s"      , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Wind speed at lowest atm level")
  var23   <- ncvar_def("QBOT"       , "kg/kg"    , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Sp. humidity at lowest atm level")
  var3    <- ncvar_def("PRECTmms"   , "mm/s"     , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Surface precipitation")
  var4    <- ncvar_def("FSDS"       , "W/m2"     , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Incident Solar radiation")
  var5    <- ncvar_def("FLDS"       , "W/m2"     , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Incident Long Wave radiation")
  var6    <- ncvar_def("ZBOT"       , "m"        , list(dimlon, dimlat, dimtime), mv, prec="float",longname="reference height")

  nc      <- nc_create(fnameout,list(var01,var02,var1,var12,var2,var23,var3,var4,var5,var6))

  # add source info metadata to file
  ncatt_put( nc, 0, "content", "Atmospheric Forcing for CLM 3.5 extracted from TSMP NRW and interpolated to CKA domain")
  ncatt_put( nc, 0, "creator", "s.poll")
  ncatt_put( nc, 0, "date_of_creation", paste0(Sys.time()))

  # put value to var 
  ncvar_put(nc,var01,lon_clm)
  ncvar_put(nc,var02,lat_clm)
  ncvar_put(nc,var1,tbotint)
  ncvar_put(nc,var12,thbotint)
  ncvar_put(nc,var2,windint)
  ncvar_put(nc,var23,qbotint)
  ncvar_put(nc,var3,tprecint)
  ncvar_put(nc,var4,fsdsint)
  ncvar_put(nc,var5,fldsint)
  ncvar_put(nc,var6,zbotint)

  nc_close(nc)

  print("nc file written")
  
} # if l100m 

if (l_10m){

#   library(raster)
  print(paste0("Do nearest neighbor interpolation for 10m"))
  
  filename_clm <- paste0(path_clm2,"/griddata_0300x0300.nc")
  nc  <- nc_open(filename_clm)
  lon_clm <- ncvar_get(nc,"LONGXY")
  lat_clm <- ncvar_get(nc,"LATIXY")
  nc_close(nc)

  # find nearest neighbor
  indnn <- array(NA,c(dim(lon_clm),2))
  for (ilon in 1:dim(lon_clm)[1]){
   for (ilat in 1:dim(lon_clm)[2]){
     ind_h <- which((abs(lon_clm[ilon,ilat]-lon)+abs(lat_clm[ilon,ilat]-lat))==min(abs(lon_clm[ilon,ilat]-lon)+abs(lat_clm[ilon,ilat]-lat)),arr.ind=TRUE)
     if (dim(ind_h)[1]==1){
       indnn[ilon,ilat,] <- ind_h
     }else{
       indnn[ilon,ilat,] <- ind_h[1,]
     } 
   }   
  }
  
  # do nearest neighbor interpolation
  tbotint <- thbotint <- windint <- qbotint <- tprecint <- fsdsint <- fldsint <- zbotint <- array(NA,c(dim(lon_clm),dim(tbotfld)[3]))
  for (ilon in 1:dim(lon_clm)[1]){
   for (ilat in 1:dim(lon_clm)[2]){
    tbotint[ilon,ilat,] <- tbotfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    thbotint[ilon,ilat,] <- thbotfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    windint[ilon,ilat,] <- windfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    qbotint[ilon,ilat,] <- qbotfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    tprecint[ilon,ilat,] <- tprecfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    fsdsint[ilon,ilat,] <- fsdsfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    fldsint[ilon,ilat,] <- fldsfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
    zbotint[ilon,ilat,] <- zbotfld[indnn[ilon,ilat,1],indnn[ilon,ilat,2],]
   }   
  }

  subdir <- paste0("m.data")
  ifelse(!dir.exists(file.path(path_out2, subdir)), dir.create(file.path(path_out2, subdir)), FALSE)

  fnameout <- paste0(path_out2,"/",subdir,"/",year,"-",monthstr,".nc")

  print(paste0("Save nc file: ",fnameout))

  # definition of dimensions
  dimlon  <- ncdim_def("lon","",1:dim(lon_clm)[1])
  dimlat  <- ncdim_def("lat","",1:dim(lon_clm)[2])
  dimtime <- ncdim_def("time","",1:dim(tbotfld)[3], unlim=TRUE)

  # Variables
  mv      <- 9.96920996838687e+36
  var01   <- ncvar_def("LONGXY"     , "degrees"  , list(dimlon, dimlat), prec="float",longname="lon")
  var02   <- ncvar_def("LATIXY"     , "degrees"  , list(dimlon, dimlat), prec="float",longname="lat") 
  var1    <- ncvar_def("TBOT"       , "K"        , list(dimlon, dimlat, dimtime), mv, prec="float",longname="temperature at lowest atm level")
  var12   <- ncvar_def("THBOT"      , "K"        , list(dimlon, dimlat, dimtime), mv, prec="float",longname="pot. temperature at lowest atm level")
  var2    <- ncvar_def("WIND"       , "m/s"      , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Wind speed at lowest atm level")
  var23   <- ncvar_def("QBOT"       , "kg/kg"    , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Sp. humidity at lowest atm level")
  var3    <- ncvar_def("PRECTmms"   , "mm/s"     , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Surface precipitation")
  var4    <- ncvar_def("FSDS"       , "W/m2"     , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Incident Solar radiation")
  var5    <- ncvar_def("FLDS"       , "W/m2"     , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Incident Long Wave radiation")
  var6    <- ncvar_def("ZBOT"       , "m"        , list(dimlon, dimlat, dimtime), mv, prec="float",longname="reference height")

  nc      <- nc_create(fnameout,list(var01,var02,var1,var12,var2,var23,var3,var4,var5,var6))

  # add source info metadata to file
  ncatt_put( nc, 0, "content", "Atmospheric Forcing for CLM 3.5 extracted from TSMP NRW and interpolated to CKA domain")
  ncatt_put( nc, 0, "creator", "s.poll")
  ncatt_put( nc, 0, "date_of_creation", paste0(Sys.time()))

  # put value to var 
  ncvar_put(nc,var01,lon_clm)
  ncvar_put(nc,var02,lat_clm)
  ncvar_put(nc,var1,tbotint)
  ncvar_put(nc,var12,thbotint)
  ncvar_put(nc,var2,windint)
  ncvar_put(nc,var23,qbotint)
  ncvar_put(nc,var3,tprecint)
  ncvar_put(nc,var4,fsdsint)
  ncvar_put(nc,var5,fldsint)
  ncvar_put(nc,var6,zbotint)

  nc_close(nc)

  print("nc file written")
  
} # if l_100m


} # for month

# end of program

#   # creat spatial points
#   sp<- SpatialPoints(coords=cbind(as.vector(lon),as.vector(lat)))
#   proj4string(sp) <- CRS('+proj=longlat +datum=WGS84')
#   
#   sp_clm <- SpatialPoints(coords=cbind(as.vector(lon_clm),as.vector(lat_clm)))
#   proj4string(sp_clm) <- CRS('+proj=longlat +datum=WGS84')
#   
#   
#   # neigherst neighbor regridding
#   varsdf = SpatialPointsDataFrame(sp, as.data.frame(t(tbotfld))
#   sp_extract <- extract(varsdf,sp)
#   r_mat_sp <- matrix(sp_extract,nrow=dim(lon_clm)[1],ncol=dim(lon_clm)[2])
# 
#   
# #   d <- pointDistance(sp_clm,sp,lonlat = F)
# #   closest <- apply(d, 1, function(x) which(x == min(x)))
# #   tbotint <- tbotfld[closest,]
