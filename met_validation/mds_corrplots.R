
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

# first things first, let's compare the hourly resolution nldas and wcr observed
# nldas downscaled
ens_seq = sprintf('%0.2d', 5)
year_seq = 1999:2014
debias_seq = sprintf('%0.3d', 5)
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
    hours = hour(date.seq)
    months = month(date.seq)
    ens.df$hour_of_day = hours
    ens.df$month = months
    Mens = cor(ens.df)
  }
}
year_seq = 1999:2008
year_seq = append(year_seq, 2010:2014)
for (y in year_seq){
  fname = paste0("Documents/metdownscaling/data_wcr/raw/hourly_wcr/WCr_1hr.", y, ".nc")
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
hours = hour(date.seq)
months = month(date.seq)
wcr.df$hour_of_day = hours
wcr.df$month = months
Mwcr = cor(wcr.df)






png("~/Documents/metdownscaling/plots/validate/met_corrplots.png", width = 10, height = 6, units = "in",res = 200)
par(mfrow=c(1,2), mar=c(1,1,1,1))
newnames = c("T", "Precip", "SW", "SH", "LW", "Pres", "WS", "Hour", "Month")
res1 <- cor.mtest(wcr.df, conf.level = .95)
res2 <- cor.mtest(wcr.df, conf.level = .99)
rownames(Mwcr) = newnames
colnames(Mwcr) = newnames
corrplot(Mwcr, p.mat = res1$p, method = 'color', number.cex = .7, type = 'lower',
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         sig.level = 0.05, insig = "blank")
title("Willow Creek Observed Correlations", line=-3)
## now the ensemble
res1 <- cor.mtest(ens.df, conf.level = .95)
res2 <- cor.mtest(ens.df, conf.level = .99)
rownames(Mens) = newnames
colnames(Mens) = newnames
corrplot(Mens, p.mat = res1$p, method = 'color', number.cex = .7, type = 'lower',
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         sig.level = 0.05, insig = "blank")
title("Random Ensemble Member Correlations", line=-3)

dev.off()



Msub = Mwcr - Mens
newnames = c("T", "Precip", "SW", "SH", "LW", "Pres", "WS")
rownames(Msub) = newnames
colnames(Msub) = newnames
corrplot(Msub, method = 'color', number.cex = .7, type = 'lower',
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90)