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
l_intp = FALSE

#directories
path_forc  = "/p/scratch/cslts/poll1/test_simulations/idealcordex11_30d_110722-110606"
lonl = 24
latl = 14

# setting specific ens
#directories
path_out = "/p/project/cslts/poll1/data//idealcordex/clm/atm_forcing/"
path_clm  = "/p/project/cslts/poll1/data//idealcordex/clm/"

year <- 2008
monthv <- c(5)
#monthv <- c(1,2,3,4,5,6,7,8,9,10,11,12)
#monthv <- c(1,2,3,4,5,6)
#monthv <- c(7,8,9,10,11,12)


l_grid = TRUE
grid_forc = "/p/project/cslts/poll1/data//idealcordex/clm/griddata_0014x0024.nc"
grid_intp = ""

#### functions #####

##--------------------------------------------------------------------------------------------
## Function to calculate distance on the globe (Wikipedia,Orthodrome)
## INPUT
## b1 lat station
## b2 lat grid point
## l1 lon station
## l2 lon grid point
## OUTPUT
## s distance in [km]
## 
##--------------------------------------------------------------------------------------------

get_nearest <- function(b1,b2,l1,l2){

f<-1/298.257223563				# Abplattung Erde
a<-6378137/1000					# Aequatorradius Erde

F<-(b1+b2)/2
G<-(b1-b2)/2
l<-(l1-l2)/2
F<-pi/180*F
G<-pi/180*G
l<-pi/180*l

S<-(sin(G))^2 *cos(l)^2 + cos(F)^2 *sin(l)^2
C<-(cos(G))^2 * cos(l)^2 + sin(F)^2 *sin(l)^2

w<-atan(sqrt(S/C))
D<-2*w*a

R<-sqrt(S*C)/w

H1<-(3*R -1)/(2*C)
H2<-(3*R+1)/(2*C)

s<-D*(1 + f*H1*(sin(F))^2*(cos(G))^2 - f*H2*(cos(F)^2*(sin(G))^2))

return(s)

}

#### start program ####

# 
 rad2deg=180/pi
 deg2rad=pi/180
 l_first=TRUE
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
  fname_h <- "test.nc"

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
  fname <-  paste0(path_forc,"/clmoas.clm2.h0.",year,"-",monthstr,"-",daystr,"-",hrstr,".nc")
#
  print(fname)
  ff     <- nc_open(fname)
# 
  if (iday==1) { 
    lon <- ncvar_get(ff,varid="longxy")
    lat <- ncvar_get(ff,varid="latixy")
    if (l_grid){
       fg    <- nc_open(grid_forc)
       lone  <- ncvar_get(fg,varid="LONE")
       lonw  <- ncvar_get(fg,varid="LONW")
       lats  <- ncvar_get(fg,varid="LATS")
       latn  <- ncvar_get(fg,varid="LATN")
       edgen <- ncvar_get(fg,varid="EDGEN")
       edgee <- ncvar_get(fg,varid="EDGEE")
       edges <- ncvar_get(fg,varid="EDGES")
       edgew <- ncvar_get(fg,varid="EDGEW")
       nc_close(fg)      
    } # if l_grid

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
  zbot  <- ncvar_get(ff,varid="ZBOT")
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
print(paste0("pot temp range ",range(thbotfld)))
#  
  subdir <- paste0("o.data")
  ifelse(!dir.exists(file.path(path_out, subdir)), dir.create(file.path(path_out, subdir)), FALSE)

  fnameout <- paste0(path_out,"/",subdir,"/",year,"-",monthstr,".nc")

  print(paste0("Save nc file: ",fnameout))

  # definition of dimensions
  dimlon  <- ncdim_def("lon","",1:dim(tbotfld)[1])
  dimlat  <- ncdim_def("lat","",1:dim(tbotfld)[2])
  dimtime <- ncdim_def("time","",1:dim(tbotfld)[3], unlim=TRUE)
#  if (l_grid){ # l_grid
#  dimedge <- ncdim_def("edge","",1)
#  }

  # Variables
  mv      <- 9.96920996838687e+36
  var01   <- ncvar_def("LONGXY"     , "degrees_east"   , list(dimlon, dimlat), prec="float",longname="longitude")
  var02   <- ncvar_def("LATIXY"     , "degrees_north"  , list(dimlon, dimlat), prec="float",longname="latitude") 
  var1    <- ncvar_def("TBOT"       , "K"              , list(dimlon, dimlat, dimtime), mv, prec="float",longname="temperature at lowest atm level")
  var12   <- ncvar_def("THBOT"      , "K"              , list(dimlon, dimlat, dimtime), mv, prec="float",longname="pot. temperature at lowest atm level")
  var2    <- ncvar_def("WIND"       , "m/s"            , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Wind speed at lowest atm level")
  var23   <- ncvar_def("QBOT"       , "kg/kg"          , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Sp. humidity at lowest atm level")
  var3    <- ncvar_def("PRECTmms"   , "mm/s"           , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Surface precipitation")
  var4    <- ncvar_def("FSDS"       , "W/m2"           , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Incident Solar radiation")
  var5    <- ncvar_def("FLDS"       , "W/m2"           , list(dimlon, dimlat, dimtime), mv, prec="float",longname="Incident Long Wave radiation")
  var6    <- ncvar_def("ZBOT"       , "m"              , list(dimlon, dimlat, dimtime), mv, prec="float",longname="reference height")
  if (l_grid){ # l_grid
  var101    <- ncvar_def("LONE"     , "degrees_east"   , list(dimlon, dimlat), prec="float",longname="longitude of east edge")
  var102    <- ncvar_def("LONW"     , "degrees_east"   , list(dimlon, dimlat), prec="float",longname="longitude of west edge")
  var103    <- ncvar_def("LATS"     , "degrees_north"  , list(dimlon, dimlat), prec="float",longname="latitude of south edge")
  var104    <- ncvar_def("LATN"     , "degrees_north"  , list(dimlon, dimlat), prec="float",longname="latitude of north edge")
  }


  if (l_grid){
  nc      <- nc_create(fnameout,list(var01,var02,var1,var12,var2,var23,var3,var4,var5,var6,var101,var102,var103,var104))
  } else {
  nc      <- nc_create(fnameout,list(var01,var02,var1,var12,var2,var23,var3,var4,var5,var6))
  } # if l_grid

  # add source info metadata to file
  ncatt_put( nc, 0, "content", "Atmospheric Forcing for CLM 3.5 extracted from CLM forcing output")
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
  if (l_grid){ # l_grid
  ncvar_put(nc,var101,lone)
  ncvar_put(nc,var102,lonw)
  ncvar_put(nc,var103,lats)
  ncvar_put(nc,var104,latn)
  }

  nc_close(nc)

  print("nc file written")

}else { # if l_calc

  fnameout <- paste0(path_out,"/o.data/",year,"-",monthstr,".nc")

  nc  <- nc_open(fnameout)
  lon <- ncvar_get(nc,"LONGXY")
  lat <- ncvar_get(nc,"LATIXY")
  tbotfld  <- ncvar_get(nc,"TBOT")
  thbotfld <- ncvar_get(nc,"THBOT")
  windfld  <- ncvar_get(nc,"WIND")
  qbotfld  <- ncvar_get(nc,"QBOT")
  tprecfld <- ncvar_get(nc,"PRECTmms")
  fsdsfld  <- ncvar_get(nc,"FSDS") 
  fldsfld  <- ncvar_get(nc,"FLDS")
  zbotfld  <- ncvar_get(nc,"ZBOT")
  nc_close(nc)

} # if l_calc
  
## Interpolate org CLM -> int CLM
  
if (l_intp){

   print(paste0("Do nearest neighbor interpolation"))
#   library(raster)
  
  filename_clm <- paste0(path_clm,"/griddata_0250x0250_rot.nc")
  nc  <- nc_open(filename_clm)
  lon_clm <- ncvar_get(nc,"LONGXY")
  lat_clm <- ncvar_get(nc,"LATIXY")
  nc_close(nc)

  # find nearest neighbor
  if (l_first ) {
  indnn <- array(NA,c(dim(lon_clm),2))
  for (ilon in 1:dim(lon_clm)[1]){
   for (ilat in 1:dim(lon_clm)[2]){
     aux_nearest <- get_nearest(lat,lat_clm[ilon,ilat],lon,lon_clm[ilon,ilat])
     indnn[ilon,ilat,] <- which(aux_nearest == min(aux_nearest),arr.ind=T)
#     indnn[ilon,ilat,] <- which((abs(lon_clm[ilon,ilat]-lon)+abs(lat_clm[ilon,ilat]-lat))==min(abs(lon_clm[ilon,ilat]-lon)+abs(lat_clm[ilon,ilat]-lat)),arr.ind=TRUE)
   }   
  }
  l_first=FALSE
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
  ncatt_put( nc, 0, "content", "Atmospheric Forcing for CLM 3.5 extracted from CLM")
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
  
} # if l_intp


} # for month

# end of program

