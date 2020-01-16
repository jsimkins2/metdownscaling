# attempting to replicate christy's envelope plot 
#

library(ncdf4)
library(gridExtra)
library(grid)
library(scales)
library(grid)
library(gtable)
library(ggplot2)
library(lubridate)
library(corrplot)
library(viridis)
library(dplyr)
library(tibbletime)
read.wcr = function(fname){
  fullname = strsplit(fname, "/")
  dataset_str = fullname[[1]][length(fullname[[1]])]
  datname = strsplit(dataset_str, "_")[[1]][1]
  data.year = substr(dataset_str, nchar(dataset_str)-6, nchar(dataset_str)-3)
  data.date = seq(from=as.POSIXct(paste0(data.year,"-1-1 0:00", tz="UTC")),to=as.POSIXct(paste0(data.year,"-12-31 23:00", tz="UTC")),by="hour")
  vars.info <- data.frame(CF.name = c("date", "air_temperature", "precipitation_flux", "surface_downwelling_shortwave_flux_in_air",
                                      "specific_humidity", "surface_downwelling_longwave_flux_in_air", "air_pressure",
                                      "eastward_wind", "northward_wind", "wind_speed"))
  df <- list()
  tem <- ncdf4::nc_open(fname)
  dim <- tem$dim
  for (j in seq_along(vars.info$CF.name)) {
    if (exists(as.character(vars.info$CF.name[j]), tem$var)) {
      df[[j]] <- ncdf4::ncvar_get(tem, as.character(vars.info$CF.name[j]))
    } else {
      df[[j]] = NA
    }
  }
  names(df) <- vars.info$CF.name
  df <- data.frame(df)
  nc_close(tem)
  if(all(is.na(df$date))){
    df$date = data.date
  }
  if(all(is.na(df$wind_speed))){
    df$wind_speed = sqrt(df$eastward_wind^2 + df$northward_wind^2)
  }
  return(df)
}
#############################################################################
# Willow Creek Training Data
#############################################################################
year_seq = 1999:2008
year_seq = append(year_seq, 2010:2014)
for (y in year_seq){
  fname = paste0("/Users/james/Documents/metdownscaling/data_wcr/raw/hourly_wcr/WCr_1hr.", y, ".nc")
  tem = read.wcr(fname = fname)
  tem$northward_wind = NULL
  tem$eastward_wind = NULL
  tem$date = NULL
  tem$precipitation_flux = tem$precipitation_flux*3600
  prec = tem
  assign(paste0("prec",y), prec)
}
wcr.df = as.vector(rbind(prec1999, prec2000, prec2001, prec2002, prec2003, prec2004, prec2005, prec2006, prec2007, 
                         prec2008, prec2010, prec2011, prec2012, prec2013, prec2014))
date.seq = seq(as.POSIXct("1999-01-01 00:00:00"), as.POSIXct("2008-12-31 23:00:00"), by="hour")
date.seq = append(date.seq, seq(as.POSIXct("2010-01-01 00:00:00"), as.POSIXct("2014-12-31 23:00:00"), by="hour"))
date.seq = strftime(date.seq)
years = year(date.seq)
wcr.df$year = years

# group the wcr.df by year, summarise it by the mean, tjhen append the sum of precip by year
grouped = group_by(wcr.df,year)
wcr.annual = summarise_each(grouped,list(mean),air_temperature,surface_downwelling_shortwave_flux_in_air,
                 specific_humidity, surface_downwelling_longwave_flux_in_air, air_pressure, wind_speed)

wcr.annual = append(wcr.annual, summarise_each(grouped,list(sum), precipitation_flux))

wcr.annual = data.frame(wcr.annual)
wcr.annual$year.1 = NULL

#############################################################################
# now read in NLDAS raw daily data 
#############################################################################

read.nldas.raw = function(fname){
  fullname = strsplit(fname, "/")
  dataset_str = fullname[[1]][length(fullname[[1]])]
  datname = strsplit(dataset_str, "_")[[1]][1]
  data.year = substr(dataset_str, nchar(dataset_str)-6, nchar(dataset_str)-3)
  data.date = seq(from=as.POSIXct(paste0(data.year,"-1-1 0:00", tz="UTC")),to=as.POSIXct(paste0(data.year,"-12-31 23:00", tz="UTC")),by="day")
  vars.info <- data.frame(CF.name = c("air_temperature_minimum","air_temperature_maximum","precipitation_flux", "surface_downwelling_shortwave_flux_in_air",
                                      "specific_humidity", "surface_downwelling_longwave_flux_in_air", "air_pressure",
                                      "eastward_wind", "northward_wind", "wind_speed"))
  df <- list()
  tem <- ncdf4::nc_open(fname)
  dim <- tem$dim
  for (j in seq_along(vars.info$CF.name)) {
    if (exists(as.character(vars.info$CF.name[j]), tem$var)) {
      df[[j]] <- ncdf4::ncvar_get(tem, as.character(vars.info$CF.name[j]))
    } else {
      df[[j]] = NA
    }
  }
  names(df) <- vars.info$CF.name
  df <- data.frame(df)
  nc_close(tem)
  df$air_temperature = (df$air_temperature_minimum + df$air_temperature_maximum)/2
  df$air_temperature_minimum=NULL
  df$air_temperature_maximum=NULL
  df$date=NA
  if(all(is.na(df$date))){
    df$date = data.date
  }
  if(all(is.na(df$wind_speed))){
    df$wind_speed = sqrt(df$eastward_wind^2 + df$northward_wind^2)
  }
  return(df)
}
year_seq = 1999:2008
year_seq = append(year_seq, 2010:2014)
for (y in year_seq){
  fname = paste0("/Users/james/Documents/metdownscaling/data_wcr/raw/NLDAS_day/NLDAS_day.", y, ".nc")
  tem = read.nldas.raw(fname = fname)
  tem$northward_wind = NULL
  tem$eastward_wind = NULL
  tem$date = NULL
  tem$precipitation_flux = tem$precipitation_flux/24*3600
  prec = tem
  assign(paste0("prec",y), prec)
}
nldas.df = as.vector(rbind(prec1999, prec2000, prec2001, prec2002, prec2003, prec2004, prec2005, prec2006, prec2007, 
                           prec2008, prec2010, prec2011, prec2012, prec2013, prec2014))
date.seq = seq(as.POSIXct("1999-01-01 00:00:00"), as.POSIXct("2008-12-31 23:00:00"), by="day")
date.seq = append(date.seq, seq(as.POSIXct("2010-01-01 00:00:00"), as.POSIXct("2014-12-31 23:00:00"), by="day"))
years = year(date.seq)
nldas.df$year = years
nldas.cor = cor(nldas.df)

grouped = group_by(nldas.df,year)
nldas.annual = summarise_each(grouped,list(mean),air_temperature,surface_downwelling_shortwave_flux_in_air,
                            specific_humidity, surface_downwelling_longwave_flux_in_air, air_pressure, wind_speed)

nldas.annual = append(nldas.annual, summarise_each(grouped,list(sum), precipitation_flux))

nldas.annual = data.frame(nldas.annual)
nldas.annual$year.1 = NULL

#############################################################################
# now read in NLDAS DOWNSCALED
#############################################################################
ens_seq = sprintf('%0.2d', seq(1,10))
year_seq = 1999:2014
debias_seq = sprintf('%0.3d', seq(1,10))
path = "~/Documents/metdownscaling/data_wcr/downscaled_v2/hourly/NLDAS/"
for (d in debias_seq){
  for (e in ens_seq){
    for (y in year_seq){
      fname = paste0(path, "NLDAS_downscaled_", d, ".", e, "/NLDAS_downscaled_", d, ".", e, ".", y, ".nc")
      tem = read.wcr(fname = fname)
      tem$northward_wind = NULL
      tem$eastward_wind = NULL
      tem$date = NULL
      tem$precipitation_flux = tem$precipitation_flux*3600
      prec = tem
      assign(paste0("prec",y), prec)
    }
    ens.df = as.vector(rbind(prec1999, prec2000, prec2001, prec2002, prec2003, prec2004, prec2005, prec2006, prec2007, 
                             prec2008, prec2009, prec2010, prec2011, prec2012, prec2013, prec2014))
    date.seq = seq(as.POSIXct("1999-01-01 00:00:00"), as.POSIXct("2014-12-31 23:00:00"), by="hour")
    years = year(date.seq)
    ens.df$year = years
    
    grouped = group_by(ens.df,year)
    ens.annual = summarise_each(grouped,list(mean),air_temperature,surface_downwelling_shortwave_flux_in_air,
                                  specific_humidity, surface_downwelling_longwave_flux_in_air, air_pressure, wind_speed)
    
    ens.annual = append(ens.annual, summarise_each(grouped,list(sum), precipitation_flux))
    
    ens.annual = data.frame(ens.annual)
    ens.annual$year.1 = NULL
    assign(paste0('debmem', d, 'ensmem', e, 'annual'), ens.annual)
  }
}

##############################
# place all 1999 tempertures in the same list
temp1999 = list()
for (d in debias_seq){
  for (e in ens_seq){
    df = eval(parse(text=paste0(('debmem', d, 'ensmem', e, 'annual')))
    bar <- subset(df, year == y)
    temp1999 = append(temp1999,df$air_temperature)
    #rm(eval(parse(text=paste0(('debmem', d, 'ensmem', e, 'annual'))))
# need to plot everything together but most importantly the ensemble with error bars
p<-ggplot(data=debmem010ensmem01annual, aes(x=year, y=air_temperature, colour='blue')) + geom_point() + geom_line()
p<-p+geom_ribbon(aes(ymin=data$lower, ymax=data$upper), linetype=2, alpha=0.1)




