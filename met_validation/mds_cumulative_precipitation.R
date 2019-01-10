
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



# EXTREME TEMPERATURES 
year_seq = 1999:2014

for (y in year_seq){
  nldas = nc_open(paste0("~/Documents/metdownscaling/data_wcr/raw/NLDAS_day/NLDAS_day.",y,".nc"))
  prec = ncvar_get(nldas, "precipitation_flux")
  prec = prec[1:365]
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
  hrly_prec = i/24
  tem = rep(hrly_prec, 24)
  nlprec=append(nlprec, tem)
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

# convert all to mm/hr, note that we converted nlprec to hourly so now we aren't dividing by 24 any more
nlprec = nlprec*3600
wcr_prec = wcr_prec*3600
# the year 2006 density plots 
# ----------------------------------------------
# DATE STRINGS
date.seq = seq(as.POSIXct("1999-01-01 00:00:00"), as.POSIXct("2014-12-31 23:00:00"), by="hour")
date.seq <- date.seq[!grepl(x = date.seq, pattern = "-02-29")]
years = year(date.seq)
df = data.frame(as.Date(date.seq), years, nlprec, wcr_prec)


##### plotting 
ens_seq = sprintf('%0.2d', 1:10)
year_seq = 1999:2014
debias_seq = sprintf('%0.3d', 1:10)
path = "~/Documents/metdownscaling/data_wcr/downscaled_v2/hourly/NLDAS/"
x = seq(1,140160)
prec.df = data.frame(x)



png("~/Documents/metdownscaling/plots/validate/cumulative_precip.png", width = 10, height = 8, units = "in",res = 200)

p = ggplot(data = df, aes(x = df$as.Date.date.seq., y = cumsum(wcr_prec), color='Observed')) + geom_line() + 
    theme(axis.text.x = element_text(angle=90, hjust = 1)) + ggtitle("Cumulative Precipitation 1999-2015") + 
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") + xlab("Year") + ylab("Cumulative Precipitation (mm)") + 
  scale_color_manual(values = c('Observed' = 'black', 'Debiased Ensemble Members'='dodgerblue4', 'NLDAS Raw'='orange')) +
  scale_shape_manual(labels = c("Observed", "Debiased Ensemble Members", "NLDAS Raw"),
                     values = c(1, 1, 1))
for (d in debias_seq){
  for (e in ens_seq){
    for (y in year_seq){
      fname = paste0(path, "NLDAS_downscaled_", d, ".", e, "/NLDAS_downscaled_", d, ".", e, ".", y, ".nc")
      tem = read.wcr(fname = fname)
      prec = tem$precipitation_flux
      prec = prec[1:8760]
      prec = prec*3600
      #tair = mean(tair)
      assign(paste0("prec",y), prec)
    }
    prec.df = as.vector(rbind(prec1999, prec2000, prec2001, prec2002, prec2003, prec2004, prec2005, prec2006, prec2007, 
                              prec2008, prec2009, prec2010, prec2011, prec2012, prec2013, prec2014))
    
    
    tem = data.frame(prec.df)
    colnames(tem) = "v1"
    p = p + geom_line(data = tem, aes(x = df$as.Date.date.seq., y = cumsum(v1), color = 'Debiased Ensemble Members'), size=0.5, alpha=0.5)
  }
}
p = p + geom_line(data = df, aes(x = df$as.Date.date.seq., y = cumsum(df$nlprec), color = 'NLDAS Raw'), size=1)
p = p +   theme(legend.title = element_blank(), 
                legend.position = c(.2, .89), 
                legend.text = element_text(size = 11)) + 
  guides(colour = guide_legend(override.aes = list(size=2, stroke=1))) 
print(p)
dev.off()