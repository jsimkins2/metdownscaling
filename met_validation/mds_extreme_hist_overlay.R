# 3,3 histograms of high low average temps/precipitation amounts
# highlight the frequency of high, low average temperatures
library(forecast)
library(ncdf4)
library(gridExtra)
library(grid)
library(scales)
library(grid)
library(gtable)
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
  
  if(all(is.na(df$date))){
    df$date = data.date
  }
  if(all(is.na(df$wind_speed))){
    df$wind_speed = sqrt(df$eastward_wind^2 + df$northward_wind^2)
  }
  return(df)
}
# EXTREME TEMPERATURES 

nldas = read.wcr(fname = "~/Documents/metdownscaling/data_wcr/downscaled_v2/hourly/NLDAS/NLDAS_downscaled_001.01/NLDAS_downscaled_001.01.1999.nc")
wcr = read.wcr(fname = "~/Documents/metdownscaling/data_wcr/raw/hourly_wcr/WCr_1hr.1999.nc")


n_hot = nldas$air_temperature[nldas$air_temperature > 300]
w_hot = wcr$air_temperature[wcr$air_temperature > 300]

minval = ifelse(min(w_hot) < min(n_hot), min(w_hot), min(n_hot))
maxval = ifelse(max(w_hot) > max(n_hot), max(w_hot), max(n_hot))
bins = seq(minval, maxval, l=20)




# Now plot on a 2 x 2


n_hot = nldas$air_temperature[nldas$air_temperature > 300]
w_hot = wcr$air_temperature[wcr$air_temperature > 300]

minval = ifelse(min(w_hot) < min(n_hot), min(w_hot), min(n_hot))
maxval = ifelse(max(w_hot) > max(n_hot), max(w_hot), max(n_hot))
hotbins = seq(minval, maxval, l=20)

n_cool = nldas$air_temperature[nldas$air_temperature < 250]
w_cool = wcr$air_temperature[wcr$air_temperature < 250]

minval = ifelse(min(w_cool) < min(n_cool), min(w_cool), min(n_cool))
maxval = ifelse(max(w_cool) > max(n_cool), max(w_cool), max(n_cool))
coolbins = seq(minval, maxval, l=20)

# now for precipitation 
nldas$precipitation_flux = nldas$precipitation_flux*3600
wcr$precipitation_flux = wcr$precipitation_flux*3600

n_dry = nldas$precipitation_flux
w_dry = wcr$precipitation_flux

minval = ifelse(min(w_dry) < min(n_dry), min(w_dry), min(n_dry))
maxval = ifelse(max(w_dry) > max(n_dry), max(w_dry), max(n_dry))
drybins = seq(minval, maxval, l=20)

# set precipitation bins to greater than 5 mm/hr
n_wet = nldas$precipitation_flux[nldas$precipitation_flux > 15]
w_wet = wcr$precipitation_flux[wcr$precipitation_flux > 15]

minval = ifelse(min(w_wet) < min(n_wet), min(w_wet), min(n_wet))
maxval = ifelse(max(w_wet) > max(n_wet), max(w_wet), max(n_wet))
wetbins = seq(minval, maxval, l=20)



#png("~/Documents/metdownscaling/plots/validate/histoverlay_temp_precip.png", width = 10, height = 6, units = "in",res = 200)
png("~/Documents/metdownscaling/plots/validate/histoverlay_extreme_temp_prec.png", width = 10, height = 10, units = "in",res = 200)

par(mfcol=c(2,2), mar=c(4,4,1,0.5), oma=c(3,2,3,1))

# now for the cool weather
hist.obs = hist(w_cool, breaks = coolbins, 
                col = alpha("deepskyblue", 1), main = "", xlab = "Temperature < 250 Kelvin")
hist.ens = hist(n_cool, breaks = coolbins, 
                col = alpha("indianred1", 1), add = T, main = "")
diff = list()
diffhist = list()
for (i in seq_along(hist.obs$counts)){
  if (hist.obs$counts[[i]] >= hist.ens$counts[[i]]){
    diffhist[[i]] = hist.ens$counts[[i]]
  } 
  if (hist.ens$counts[[i]] > hist.obs$counts[[i]]){
    diffhist[[i]] = hist.obs$counts[[i]]
  }
}

diff = as.integer(diffhist)

diffhist = hist.ens
diffhist$counts = diff
plot(diffhist, add = T, col = "mediumpurple4")


# dry weather
hist.obs = hist(w_dry, breaks = drybins, 
                col = alpha("deepskyblue", 1), main = "", xlab = "Precipitation (mm/hr)")
hist.ens = hist(n_dry, breaks = drybins, 
                col = alpha("indianred1", 1), add = T, main = "")
diff = list()
diffhist = list()
for (i in seq_along(hist.obs$counts)){
  if (hist.obs$counts[[i]] >= hist.ens$counts[[i]]){
    diffhist[[i]] = hist.ens$counts[[i]]
  } 
  if (hist.ens$counts[[i]] > hist.obs$counts[[i]]){
    diffhist[[i]] = hist.obs$counts[[i]]
  }
}

diff = as.integer(diffhist)

diffhist = hist.ens
diffhist$counts = diff
plot(diffhist, add = T, col = "mediumpurple4")


# hot weather
hist.obs = hist(w_hot, breaks = hotbins, 
                col = alpha("deepskyblue", 1), main = "", xlab = "Temperature > 300 Kelvin")
hist.ens = hist(n_hot, breaks = hotbins, 
                col = alpha("indianred1", 1), add = T, main = "")
diff = list()
diffhist = list()
for (i in seq_along(hist.obs$counts)){
  if (hist.obs$counts[[i]] >= hist.ens$counts[[i]]){
    diffhist[[i]] = hist.ens$counts[[i]]
  } 
  if (hist.ens$counts[[i]] > hist.obs$counts[[i]]){
    diffhist[[i]] = hist.obs$counts[[i]]
  }
}

diff = as.integer(diffhist)

diffhist = hist.ens
diffhist$counts = diff
plot(diffhist, add = T, col = "mediumpurple4")
#legend("topright", legend = c("Observations", "Ensemble Member", "Overlap"), col = c(alpha("deepskyblue", 1), alpha("indianred1", 1), "mediumpurple4"),
#lty = c(1,1,1), lwd = 5)
#title("1999 Observations and Random Ensemble")


# wet weather
hist.obs = hist(w_wet, breaks = wetbins, 
                col = alpha("deepskyblue", 1), main = "", xlab = "Precipitation > 15 mm/hr", xlim = range(0,250))
hist.ens = hist(n_wet, breaks = wetbins, 
                col = alpha("indianred1", 1), add = T, main = "", xlim = range(0,250))
diff = list()
diffhist = list()
for (i in seq_along(hist.obs$counts)){
  if (hist.obs$counts[[i]] >= hist.ens$counts[[i]]){
    diffhist[[i]] = hist.ens$counts[[i]]
  } 
  if (hist.ens$counts[[i]] > hist.obs$counts[[i]]){
    diffhist[[i]] = hist.obs$counts[[i]]
  }
}

diff = as.integer(diffhist)

diffhist = hist.ens
diffhist$counts = diff
plot(diffhist, add = T, col = "mediumpurple4")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(x = "bottom",inset = 0.01,
       legend = c("Observations", "Ensemble Member", "Overlap"), col = c(alpha("deepskyblue", 1), alpha("indianred1", 1), "mediumpurple4"), 
       lwd=3, cex=1, horiz = TRUE)
title("1999 Observations vs. Random Ensemble Member", line = -2)
dev.off()













png("~/Documents/metdownscaling/plots/validate/histoverlay_maxtemp.png", width = 10, height = 6, units = "in",res = 200)
hist.obs = hist(w_hot, breaks = bins, 
col = alpha("deepskyblue", 1), main = "", xlab = "Temperature (Kelvin)")
hist.ens = hist(n_hot, breaks = bins, 
col = alpha("indianred1", 1), add = T, main = "")
diff = list()
diffhist = list()
for (i in seq_along(hist.obs$counts)){
if (hist.obs$counts[[i]] >= hist.ens$counts[[i]]){
diffhist[[i]] = hist.ens$counts[[i]]
} 
if (hist.ens$counts[[i]] > hist.obs$counts[[i]]){
diffhist[[i]] = hist.obs$counts[[i]]
}
}

diff = as.integer(diffhist)

diffhist = hist.ens
diffhist$counts = diff
plot(diffhist, add = T, col = "mediumpurple4")
legend("topright", legend = c("Observations", "Ensemble Member", "Overlap"), col = c(alpha("deepskyblue", 1), alpha("indianred1", 1), "mediumpurple4"),
lty = c(1,1,1), lwd = 5)
title("1999 Observations and Random Ensemble Member")
dev.off()
