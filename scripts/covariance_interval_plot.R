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

# function for calculating the average correlation of a correlation matrix
avgcor <- function(x){mean(abs(x[lower.tri(x)]))}

# load in the datafiles
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
hours = hour(date.seq)
days = yday(date.seq)
months = month(date.seq)
wcr.df$hour_of_day = hours
wcr.df$day = days
wcr.df$month = months
wcr.df$year = years
# group the wcr.df by year, summarise it by the mean, tjhen append the sum of precip by year
grouped = group_by(wcr.df,year)
wcr.annual = summarise_each(grouped,list(mean),air_temperature,surface_downwelling_shortwave_flux_in_air,
                            specific_humidity, surface_downwelling_longwave_flux_in_air, air_pressure, wind_speed,
                            precipitation_flux)
# if we want to take the sum of precipitation flux
#wcr.annual = append(wcr.annual, summarise_each(grouped,list(sum), precipitation_flux))

wcr.annual = data.frame(wcr.annual)
wcr.annual$year.1 = NULL
wcr.annual$year = NULL
wcr.cor = cor(wcr.annual)

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
date.seq = strftime(date.seq)
years = year(date.seq)
hours = hour(date.seq)
days = yday(date.seq)
months = month(date.seq)
nldas.df$hour_of_day = hours
nldas.df$day = days
nldas.df$month = months
nldas.df$year = years

grouped = group_by(nldas.df,year)
nldas.annual = summarise_each(grouped,list(mean),air_temperature,surface_downwelling_shortwave_flux_in_air,
                              specific_humidity, surface_downwelling_longwave_flux_in_air, air_pressure, wind_speed)

nldas.annual = append(nldas.annual, summarise_each(grouped,list(sum), precipitation_flux))

nldas.annual = data.frame(nldas.annual)
nldas.annual$year.1 = NULL
nldas.annual$year = NULL
nldas.cor = cor(nldas.annual)
#############################################################################
# now read in NLDAS DOWNSCALED
#############################################################################
ens_seq = sprintf('%0.2d', seq(1,10))
year_seq = 1999:2014
debias_seq = sprintf('%0.3d', seq(1,10))
path = "~/Documents/metdownscaling/data_wcr/downscaled_v2/hourly/NLDAS/"
ens.cor.vec = vector()
ens.cor.index = vector()
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
    date.seq = strftime(date.seq)
    years = year(date.seq)
    hours = hour(date.seq)
    days = yday(date.seq)
    months = month(date.seq)
    ens.df$hour_of_day = hours
    ens.df$day = days
    ens.df$month = months
    ens.df$year = years
    grouped = group_by(ens.df,year)
    ens.annual = summarise_each(grouped,list(mean),air_temperature,surface_downwelling_shortwave_flux_in_air,
                                  specific_humidity, surface_downwelling_longwave_flux_in_air, air_pressure, wind_speed, precipitation_flux)
    
    ens.annual = data.frame(ens.annual)
    ens.annual$year.1 = NULL
    ens.annual$year = NULL
    ens.cor = cor(ens.annual)
    ens.cor.vec = append(ens.cor.vec, avgcor(ens.cor))
    ens.cor.index = append(ens.cor.index, paste0(d,'.',e))
  }
}

ens.cor.mean = which(abs(ens.cor.vec - mean(ens.cor.vec)) == min(abs(ens.cor.vec - mean(ens.cor.vec))))
ens.cor.pos.sd = which(abs(ens.cor.vec - (mean(ens.cor.vec) + sd(ens.cor.vec))) == min(abs(ens.cor.vec - (mean(ens.cor.vec) + sd(ens.cor.vec)))))
ens.cor.neg.sd = which(abs(ens.cor.vec - (mean(ens.cor.vec) - sd(ens.cor.vec))) == min(abs(ens.cor.vec - (mean(ens.cor.vec) - sd(ens.cor.vec)))))

ind.list = list(ens.cor.index[ens.cor.mean], ens.cor.index[ens.cor.pos.sd], ens.cor.index[ens.cor.neg.sd])
name_enscor = c('ens.cor.mean', 'ens.cor.pos.sd', 'ens.cor.neg.sd')
name_ensdf = c('ens.df.mean', 'ens.df.pos.sd', 'ens.df.neg.sd')
for (i in seq_len(length(ind.list))){
  d = substr(ind.list[[i]],1,3)
  e = substr(ind.list[[i]],5,6)
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
  date.seq = strftime(date.seq)
  years = year(date.seq)
  hours = hour(date.seq)
  days = yday(date.seq)
  months = month(date.seq)
  ens.df$hour_of_day = hours
  ens.df$day = days
  ens.df$month = months
  ens.df$year = years
  grouped = group_by(ens.df,year)
  ens.annual = summarise_each(grouped,list(mean),air_temperature,surface_downwelling_shortwave_flux_in_air,
                              specific_humidity, surface_downwelling_longwave_flux_in_air, air_pressure, wind_speed, precipitation_flux)
  
  ens.annual = data.frame(ens.annual)
  ens.annual$year.1 = NULL
  ens.annual$year = NULL
  assign(name_enscor[i], cor(ens.annual))
  assign(name_ensdf[i], ens.annual)
}


png("~/Documents/metdownscaling/plots/validate/annual_met_corrplots.png", width = 4, height = 16, units = "in",res = 200)
par(mfrow=c(5,1), mar=c(1,1,1,1))
newnames = c("AirT", "Precip", "SW", "SHum", "LW", "Pres", "WS")

# Plot Willow Creek Correlations
res1 <- cor.mtest(wcr.df, conf.level = .95)
res2 <- cor.mtest(wcr.df, conf.level = .99)
rownames(wcr.cor) = newnames
colnames(wcr.cor) = newnames
corrplot(wcr.cor, p.mat = res1$p, method = 'color', number.cex = .7, type = 'lower',
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         sig.level = 0.05, insig = "blank")
title("Willow Creek Aggregated Daily Correlations")

# Plot the NLDAS Daily Correlations
res1 <- cor.mtest(nldas.df, conf.level = .95)
res2 <- cor.mtest(nldas.df, conf.level = .99)
rownames(nldas.cor) = newnames
colnames(nldas.cor) = newnames
corrplot(nldas.cor, p.mat = res1$p, method = 'color', number.cex = .7, type = 'lower',
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         sig.level = 0.05, insig = "blank")
title("Raw NLDAS Daily Correlations")

# Plot the mean ensemble correlations
res1 <- cor.mtest(ens.df, conf.level = .95)
res2 <- cor.mtest(ens.df, conf.level = .99)
rownames(ens.cor.mean) = newnames
colnames(ens.cor.mean) = newnames
corrplot(ens.cor.mean, p.mat = res1$p, method = 'color', number.cex = .7, type = 'lower',
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         sig.level = 0.05, insig = "blank")
title("Mean Ensemble Aggregated Daily Correlations")

# Plot the +sd ensemble correlations
res1 <- cor.mtest(ens.df, conf.level = .95)
res2 <- cor.mtest(ens.df, conf.level = .99)
rownames(ens.cor.pos.sd) = newnames
colnames(ens.cor.pos.sd) = newnames
corrplot(ens.cor.pos.sd, p.mat = res1$p, method = 'color', number.cex = .7, type = 'lower',
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         sig.level = 0.05, insig = "blank")
title("+1 SD Ensemble Aggregated Daily Correlations")

# Plot the -sd ensemble correlations
res1 <- cor.mtest(ens.df, conf.level = .95)
res2 <- cor.mtest(ens.df, conf.level = .99)
rownames(ens.cor.neg.sd) = newnames
colnames(ens.cor.neg.sd) = newnames
corrplot(ens.cor.neg.sd, p.mat = res1$p, method = 'color', number.cex = .7, type = 'lower',
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         sig.level = 0.05, insig = "blank")
title("-1 SD Ensemble Aggregated Daily Correlations")
dev.off()


