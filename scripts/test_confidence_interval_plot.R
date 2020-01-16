library(ncdf4)
library(ggplot2)
library(cowplot)

# Ensemble directories
wd.base <- "/Users/james/Documents/metdownscaling/data_wcr/"
path.pecan <- "~/Desktop/Research/pecan"

# Site name for indexing
site.name = "WILLOWCREEK"
vers=".v2"
site.lat  =  45.805822 # 45°48′21″N
site.lon  = -90.079722 # 90°04′47″W

# Setting up some file paths, etc
path.raw.base <- file.path(wd.base, "raw/hourly_wcr/")
path.nldas.base <- file.path(wd.base, "raw/NLDAS_day/")
path.tdm.base <- file.path(wd.base, "downscaled_v2/hourly/NLDAS/")
# -----------------------------------
# 1. Read in met data
#    1.1. Raw 
#    1.2. Bias-corrected (summarize)
# -----------------------------------
cols.meta <- c("type", "source", "Date", "Year", "DOY", "Hour")
vars.CF <- c("air_temperature", "precipitation_flux", "surface_downwelling_shortwave_flux_in_air", "surface_downwelling_longwave_flux_in_air", "air_pressure", "specific_humidity", "wind_speed")
vars.short <- c("tair", "precip", "swdown", "lwdown", "press", "qair", "wind")

extract.met <- function(train.path){
  files.train <- dir(train.path, ".nc")
  
  yrs.file <- strsplit(files.train, "[.]")
  yrs.file <- matrix(unlist(yrs.file), ncol=length(yrs.file[[1]]), byrow=T)
  yrs.file <- as.numeric(yrs.file[,ncol(yrs.file)-1]) # Assumes year is always last thing before the file extension
  
  met.out <- data.frame()
  for(i in 1:length(files.train)){
    yr.now <- yrs.file[i]
    
    ncT <- ncdf4::nc_open(file.path(train.path, files.train[i]))
    
    # Set up the time data frame to help index
    nday <- ifelse(lubridate::leap_year(yr.now), 366, 365)
    ntime <- length(ncT$dim$time$vals)
    step.day <- nday/ntime
    step.hr  <- step.day*24
    stamps.hr <- seq(step.hr/2, by=step.hr, length.out=1/step.day) # Time stamps centered on period
    
    # Create a data frame with all the important time info
    # center the hour step
    df.tmp <- data.frame(Year=yr.now, DOY=rep(1:nday, each=1/step.day), Hour=rep(stamps.hr, length.out=ntime))
    df.tmp$Date <- strptime(paste(df.tmp$Year, df.tmp$DOY, df.tmp$Hour, sep="-"), format=("%Y-%j-%H"), tz="UTC")
    
    # Extract the met info, making matrices with the appropriate number of ensemble members
    for(v in names(ncT$var)){
      df.tmp[,v] <- ncdf4::ncvar_get(ncT, v)
    }
    ncdf4::nc_close(ncT)
    
    met.out <- rbind(met.out, df.tmp)
  } # End looping through training data files
  
  return(met.out)
} # End function

dat.wcr <- extract.met(train.path = file.path(path.raw.base))
dat.wcr$wind_speed <- sqrt(dat.wcr$northward_wind^2 + dat.wcr$eastward_wind^2)
dat.wcr$source <- "Ameriflux"
dat.wcr$type <- "raw-tower"
summary(dat.wcr)

dat.nldas <- extract.met(train.path = file.path(path.nldas.base))
dat.nldas$air_temperature <- (dat.nldas$air_temperature_minimum + dat.nldas$air_temperature_maximum)/2
dat.nldas$source <- "NLDAS"
dat.nldas$type <- "raw-gridded"
summary(dat.nldas)

ens.mems <- dir(file.path(path.tdm.base), "NLDAS_downscaled_")
dat.ens <- data.frame()
for(ENS in ens.mems){
  ens.id <- stringr::str_split(ENS, "_")[[1]][3]
  dat.mem <- extract.met(train.path = file.path(path.tdm.base, ENS))
  dat.mem$source <- paste("ensemble", ens.id, sep="-")
  # summary(dat.mem)
  
  dat.ens <- rbind(dat.ens, dat.mem)
}
dat.ens$type <- "downscaled"

# Binding everything together
# cols.meta[!cols.meta %in% names(dat.wcr)]
# vars.CF[!vars.CF %in% names(dat.wcr)]
dat.all <- rbind(dat.wcr[,c(cols.meta, vars.CF)], dat.nldas[,c(cols.meta, vars.CF)], dat.ens[,c(cols.meta, vars.CF)])
summary(dat.all)
dim(dat.all)

dat.long <- stack(dat.all[,vars.CF])
dat.long[,cols.meta] <- dat.all[,cols.meta]
summary(dat.long)
# -----------------------------------

# -----------------------------------
# Aggregating to different scales and graphing
# -----------------------------------
# -------------------
# Annual means
# -------------------
# Aggregate to annual means, preserving each ensemble member (source) along the way 
dat.ann <- aggregate(dat.long$values, by=dat.long[,c("type", "source", "Year", "ind")], FUN=mean, na.rm=F)
names(dat.ann)[which(names(dat.ann)=="x")] <- "values"

# Aggregate again to get the ensemble mean & CI
dat.ann2 <- aggregate(dat.ann$values, by=dat.ann[,c("type", "Year", "ind")], FUN=mean, na.rm=F)
names(dat.ann2)[which(names(dat.ann2)=="x")] <- "val.mean"
dat.ann2$val.025 <- aggregate(dat.ann$values, by=dat.ann[,c("type", "Year", "ind")], FUN=quantile, 0.025, na.rm=F)$x
dat.ann2$val.975 <- aggregate(dat.ann$values, by=dat.ann[,c("type", "Year", "ind")], FUN=quantile, 0.975, na.rm=F)$x
dat.ann2$scale <- "annual"
dat.ann2$type <- factor(dat.ann2$type, levels=c("raw-tower", "raw-gridded", "downscaled"))
summary(dat.ann2)

ggplot(data=dat.ann2) +
  facet_grid(ind ~ scale, scales="free_y", switch="y") +
  geom_ribbon(aes(x=Year, ymin=val.025, ymax=val.975, fill=type), alpha=0.5) +
  geom_line(aes(x=Year, y=val.mean, color=type)) +
  scale_color_manual(values=c("black", "red2", "blue2")) + 
  scale_fill_manual(values=c("black", "red2", "blue2")) + 
  theme_bw() +
  theme(strip.placement = "outside",
        axis.title.y=element_blank(),
        legend.position="top")