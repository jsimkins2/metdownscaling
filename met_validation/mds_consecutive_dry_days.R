
library(ncdf4)
library(gridExtra)
library(grid)
library(scales)
library(grid)
library(gtable)
library(ggplot2)
library(lubridate)
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

cumul_zeros <- function(x)  {
  x <- !x
  rl <- rle(x)
  len <- rl$lengths
  v <- rl$values
  cumLen <- cumsum(len)
  z <- x
  # replace the 0 at the end of each zero-block in z by the 
  # negative of the length of the preceding 1-block....
  iDrops <- c(0, diff(v)) < 0
  z[ cumLen[ iDrops ] ] <- -len[ c(iDrops[-1],FALSE) ]
  # ... to ensure that the cumsum below does the right thing.
  # We zap the cumsum with x so only the cumsums for the 1-blocks survive:
  x*cumsum(z)
}


# EXTREME TEMPERATURES 
year_seq = 1999:2014

for (y in year_seq){
  nldas = nc_open(paste0("~/Documents/metdownscaling/data_wcr/raw/NLDAS_day/NLDAS_day.",y,".nc"))
  prec = ncvar_get(nldas, "precipitation_flux")
  tair = (ncvar_get(nldas, "air_temperature_maximum") + ncvar_get(nldas, "air_temperature_minimum"))/2
  tair = tair[1:365]
  prec = prec[1:365]
  assign(paste0("tair",y), tair)
  assign(paste0("prec",y), prec)
}

x = seq(1,15*365)
nldas_temp = data.frame(x)
nldas_prec = data.frame(x)
nldas_temp= as.vector(rbind(tair1999, tair2000, tair2001, tair2002, tair2003, tair2004, tair2005, tair2006, tair2007, 
                            tair2008, tair2009, tair2010, tair2011, tair2012, tair2013, tair2014))
nldas_prec= as.vector(rbind(prec1999, prec2000, prec2001, prec2002, prec2003, prec2004, prec2005, prec2006, prec2007, 
                            prec2008, prec2009, prec2010, prec2011, prec2012, prec2013, prec2014))

nltemp = numeric()
nlprec = numeric()
for (i in nldas_prec){
  tem = rep(i, 24)
  nlprec=append(nlprec, tem)
  
  tem = rep(i, 24)
  nltemp=append(tem,nltemp)
}
# now for willow creek
for (y in year_seq){
  wcr = nc_open(paste0("~/Documents/metdownscaling/data_wcr/raw/hourly_wcr/WCr_1hr.",y,".nc"))
  prec = ncvar_get(wcr, "precipitation_flux")
  tair = ncvar_get(wcr, "air_temperature")
  tair = tair[1:8760]
  prec = prec[1:8760]
  assign(paste0("tair",y), tair)
  assign(paste0("prec",y), prec)
}

x = seq(1,140160)
wcr_temp = data.frame(x)
wcr_prec = data.frame(x)
wcr_temp= as.vector(rbind(tair1999, tair2000, tair2001, tair2002, tair2003, tair2004, tair2005, tair2006, tair2007, 
                          tair2008, tair2009, tair2010, tair2011, tair2012, tair2013, tair2014))
wcr_prec= as.vector(rbind(prec1999, prec2000, prec2001, prec2002, prec2003, prec2004, prec2005, prec2006, prec2007, 
                          prec2008, prec2009, prec2010, prec2011, prec2012, prec2013, prec2014))



#ens_seq = sprintf('%0.2d', 1:10)
ens_seq = sprintf('%0.2d', 2)
year_seq = 1999:2014
debias_seq = sprintf('%0.3d', 2)
path = "~/Documents/metdownscaling/data_wcr/downscaled_v2/hourly/NLDAS/"
x = seq(1,140160)
temp.df = data.frame(x)
prec.df = data.frame(x)
for (d in debias_seq){
  for (e in ens_seq){
    for (y in year_seq){
      fname = paste0(path, "NLDAS_downscaled_", d, ".", e, "/NLDAS_downscaled_", d, ".", e, ".", y, ".nc")
      tem = read.wcr(fname = fname)
      tair = tem$air_temperature
      prec = tem$precipitation_flux
      tair = tair[1:8760]
      prec = prec[1:8760]
      prec = prec*3600
      #tair = mean(tair)
      assign(paste0("tair",y), tair)
      assign(paste0("prec",y), prec)
    }
    temp.df[paste0("tair",d, ".", e)]= as.vector(rbind(tair1999, tair2000, tair2001, tair2002, tair2003, tair2004, tair2005, tair2006, tair2007, 
                                                       tair2008, tair2009, tair2010, tair2011, tair2012, tair2013, tair2014))
    prec.df[paste0("prec",d, ".", e)]= as.vector(rbind(prec1999, prec2000, prec2001, prec2002, prec2003, prec2004, prec2005, prec2006, prec2007, 
                                                       prec2008, prec2009, prec2010, prec2011, prec2012, prec2013, prec2014))
    
  }
}
# convert all to mm/hr
prec.df$x = NULL
temp.df$x = NULL
nlprec = nlprec/24*3600
wcr_prec = wcr_prec*3600
# the year 2006 density plots 


# ----------------------------------------------
# GROWING SEASOND ROUGHT DROUGHT DROUGHT
date.seq = seq(as.POSIXct("1999-01-01 00:00:00"), as.POSIXct("2014-12-31 23:00:00"), by="hour")
doy = yday(date.seq)
doy = doy[doy < 366]
wcr.grow = data.frame(doy, wcr_prec, wcr_temp)
dwnsc.grow = data.frame(doy, prec.df, temp.df)
nldas.grow = data.frame(doy, nlprec, nltemp)

#growing season May 1 to September 1...roughly
wcr.grow = subset(wcr.grow, wcr.grow$doy >121 & wcr.grow$doy <244)
dwnsc.grow = subset(dwnsc.grow, dwnsc.grow$doy >121 & dwnsc.grow$doy <244)
nldas.grow = subset(nldas.grow, nldas.grow$doy >121 & nldas.grow$doy <244)


d.days = 7
# wcr.grow < .1 mm/hr is considered the cutoff
wcr.grow$wcr_prec[wcr.grow$wcr_prec < .1] = 0
wcr.zeros = cumul_zeros(wcr.grow$wcr_prec)
wcr.drought = wcr.zeros[wcr.zeros >= d.days*24]
wcr.drought = wcr.drought[wcr.drought == d.days*24]

nldas.grow$nlprec[nldas.grow$nlprec < .1] = 0
nldas.zeros = cumul_zeros(nlraw.grow$nlprec)
nldas.drought = nldas.zeros[nldas.zeros >= d.days*24]
nldas.drought = nldas.drought[nldas.drought == d.days*24]

dwnsc.grow$prec002.02[dwnsc.grow$prec002.02 < .1] = 0
dwnsc.zeros = cumul_zeros(dwnsc.grow$prec002.02)
dwnsc.drought = dwnsc.zeros[dwnsc.zeros >= d.days*24]
dwnsc.drought = dwnsc.drought[dwnsc.drought == d.days*24]


which(nldas.zeros == d.days*24)
which(wcr.zeros == d.days*24)

png("~/Documents/metdownscaling/plots/validate/consecutive_drought.png", width=12, height=6,res = 200, units = "in")
par(mfrow=c(1,1), cex.axis=.5)
plot(dwnsc.zeros, col=c("dodgerblue4"),main = "Consecutive Dry Hours During Growing Season (May 1st - Sep 1st)", 
     xlab = "Growing Season Year", type = "l", ylim = c(0, max(dwnsc.zeros)), xaxt = "n", lwd = 2, ylab = "Hours")
axis(1, at = seq(0,46848 - 46848/16,46848/16), labels = seq(1999,2014,1))
lines(nldas.zeros, col=c("orange"), lwd = 2)
lines(wcr.zeros, col=c("black"), lwd = 2)
legend("topright", legend = c("Random Ensemble Member", "Observed", "NLDAS Raw"), lty = c(1,1), col = c("dodgerblue4", "black", "orange"), lwd = 4,cex=0.8 )
abline(a = 0, b = 0, h = d.days*24, lty = c(2), col= 'purple')
text(-900, d.days*24 + 4, paste0(d.days," days"), col = 'purple',cex = 0.6)
dev.off()


#during x summers, we have 2 periods of 7 day drought

#--------------------------------------
# COUNTING AVERAGE NUMBER OF DROUGHT PERIODS FOR ALL ENSEMBLE MEMBERS, MAKE A DISTRIBUTION
#--------------------------------------
ens_drought = numeric()
ens_seq = sprintf('%0.2d', 1:10)
year_seq = 1999:2014
debias_seq = sprintf('%0.3d', 1:10)
path = "~/Documents/metdownscaling/data_wcr/downscaled_v2/hourly/NLDAS/"
x = seq(1,140160)
prec.df = data.frame(x)
for (d in debias_seq){
  for (e in ens_seq){
    for (y in year_seq){
      fname = paste0(path, "NLDAS_downscaled_", d, ".", e, "/NLDAS_downscaled_", d, ".", e, ".", y, ".nc")
      tem = read.wcr(fname = fname)
      #tair = tem$air_temperature
      prec = tem$precipitation_flux
      #nc_close(tem)
      #tair = tair[1:8760]
      prec = prec[1:8760]
      prec = prec*3600
      #tair = mean(tair)
      assign(paste0("prec",y), prec)
    }
    #temp.df[paste0("tair",d, ".", e)]= as.vector(rbind(tair1999, tair2000, tair2001, tair2002, tair2003, tair2004, tair2005, tair2006, tair2007, 
    #tair2008, tair2009, tair2010, tair2011, tair2012, tair2013, tair2014))
    prec.df= as.vector(rbind(prec1999, prec2000, prec2001, prec2002, prec2003, prec2004, prec2005, prec2006, prec2007, 
                             prec2008, prec2009, prec2010, prec2011, prec2012, prec2013, prec2014))
    
    rm(prec1999, prec2000, prec2001, prec2002, prec2003, prec2004, prec2005, prec2006, prec2007, 
       prec2008, prec2009, prec2010, prec2011, prec2012, prec2013, prec2014)
    
    date.seq = seq(as.POSIXct("1999-01-01 00:00:00"), as.POSIXct("2014-12-31 23:00:00"), by="hour")
    doy = yday(date.seq)
    doy = doy[doy < 366]
    dwnsc.grow = data.frame(doy, prec.df)
    
    #growing season May 1 to September 1...roughly
    dwnsc.grow = subset(dwnsc.grow, dwnsc.grow$doy >121 & dwnsc.grow$doy <244)
    
    
    d.days = 7
    # wcr.grow < .1 mm/hr is considered the cutoff
    dwnsc.grow$prec.df[dwnsc.grow$prec.df < .1] = 0
    dwnsc.zeros = cumul_zeros(dwnsc.grow$prec.df)
    dwnsc.drought = dwnsc.zeros[dwnsc.zeros >= d.days*24]
    dwnsc.drought = dwnsc.drought[dwnsc.drought == d.days*24]
    
    ens_drought= append(ens_drought, length(dwnsc.drought))
  }
}


png("~/Documents/metdownscaling/plots/validate/ensemble_drought_growing_season.png", width = 10, height = 6, units = "in",res = 200)
p = ggplot(ens_drought) + 
  geom_histogram(aes(x=ens_drought$`Random Ensemble Member`,color='darkblue'), color = 'darkblue', fill = 'lightblue',   size=1.5, bins=20)+ 
  ggtitle("Ensemble Member Growing Season Drought Periods (1999-2014)") + xlab(" Number of Drought Periods (7 days)") + ylab("Count")
p = p + geom_vline(xintercept = mean(ens_drought$`Random Ensemble Member`), color = 'seagreen', linetype="dashed", size=1.5)
p = p + annotate("text", label = paste0("Mean = ", mean(ens_drought$`Random Ensemble Member`)),color = 'seagreen',size=4.5,  x = 20, y = 14)
print(p)
dev.off()

#--------------------------------------
#--------------------------------------


#--------------------------------------
#--------------------------------------





















# -------------------------------------
# EXTREME PRECIPITATION 


`GFDL-ESM2Gprec` = `GFDL-ESM2Gprec`[`GFDL-ESM2Gprec` > 0.001111111] # 4 mm/hr
`HadGEM2-CC365prec` = `HadGEM2-CC365prec`[`HadGEM2-CC365prec` > 0.001111111]
bins = seq(min(`GFDL-ESM2Gprec`), max(`GFDL-ESM2Gprec`), l=20)
h3 = hist(`GFDL-ESM2Gprec`, breaks = bins, main = "GFDL-ESM2G : Carbon Sink", col = "seagreen")
h4 = hist(`HadGEM2-CC365prec`, breaks = bins, main = "HadGEM2-CC365 : Carbon Source",col = "moccasin")


png("MACA/plots/x.prec.hist.png", width=12, height=6,res = 200, units = "in")
par(oma=c(0,0,2,0))
par(mfrow=c(1,2))
plot(h3, col=c("seagreen"), xlab = "Precipitation Flux (kg m-2 s-1)", main = "", ylim = c(0,5), cex = .5, xlim = c(.0012, .002))

mtext("GFDL-ESM2G : Carbon Sink")
text(.00124,5, labels = "Count = 24", cex = .7)
plot(h4, col=c("moccasin"), xlab = "Precipitation Flux (kg m-2 s-1)", main = "" ,ylim = c(0,5), cex = .5, xlim = c(.0012, .002))
mtext("HadGEM2-CC365 : Carbon Source")
text(.00124,5, labels = "Count = 2", cex = .7)
title(main="Precipitation Flux > 4 mm/hr",outer=T)
dev.off()
