
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

nldas = nc_open("~/Documents/metdownscaling/data_wcr/raw/NLDAS_day/NLDAS_day.2006.nc")
nldas_prec = ncvar_get(nldas, "precipitation_flux")
nldas = data.frame(nldas_prec)
wcr = read.wcr(fname = "~/Documents/metdownscaling/data_wcr/raw/hourly_wcr/WCr_1hr.2006.nc")




vars.info <- data.frame(CF.name = c("air_temperature", "precipitation_flux", "surface_downwelling_shortwave_flux_in_air", 
                                    "specific_humidity", "wind_speed"))

#ens_seq = sprintf('%0.2d', 1:10)
ens_seq = sprintf('%0.2d', 1:10)
year_seq = 2006
debias_seq = sprintf('%0.3d', 1:10)
path = "~/Documents/metdownscaling/data_wcr/downscaled_v2/hourly/NLDAS/"
x = seq(1,8760)
temp.df = data.frame(x)
prec.df = data.frame(x)
for (d in debias_seq){
  for (y in year_seq){
    for (e in ens_seq){
      fname = paste0(path, "NLDAS_downscaled_", d, ".", e, "/NLDAS_downscaled_", d, ".", e, ".", y, ".nc")
      tem = read.wcr(fname = fname)
      tair = tem$air_temperature
      prec = tem$precipitation_flux
      tair = tair[1:8760]
      prec = prec[1:8760]
      prec = prec*3600
      #tair = mean(tair)
      #assign(paste0("tair",d, ".", e), tair)
      #assign(paste0("prec",d, ".", e), prec)
      if (exists("temp.df") == FALSE){
        temp.df = data.frame(tair)
        prec.df = data.frame(prec)
        colnames(temp.df) = paste0("tair",d, ".", e)
        colnames(prec.df) = paste0("prec",d, ".", e)
      } else {
        temp.df[paste0("tair",d, ".", e)]= tair
        prec.df[paste0("prec",d, ".", e)]= prec
      }
    }
  }
}
# convert all to mm/hr

nldas$nldas_prec = nldas$nldas_prec/24*3600
wcr$precipitation_flux = wcr$precipitation_flux*3600
# the year 2006 density plots 

library(ggplot2)

png("~/Documents/metdownscaling/plots/validate/debens_prec_densities.png", width = 10, height = 6, units = "in",res = 200)

p = ggplot(nldas) + 
  geom_density(data = nldas, aes(x=nldas_prec, y=..density.., color = 'NLDAS Raw'), size=0.5)  + xlab("Precipitation (mm/hr)") + ggtitle("2006 Gaussian Precipitation Densities") + 
  scale_color_manual(values = c('Observed' = 'black', 'Debiased Ensemble Members'='dodgerblue4', 'NLDAS Raw'='orange')) +
  scale_shape_manual(labels = c("Observed", "Debiased Ensemble Members", "NLDAS Raw"),
                     values = c(1, 1, 1))
                    
for (d in debias_seq){
  for (e in ens_seq){
    tem = prec.df[paste0("prec",d, ".", e)]
    colnames(tem) = "v1"
    p = p + geom_density(data = tem, aes(x=v1, y=..density.., color='Debiased Ensemble Members'), size=0.7, fill="lightblue", alpha=0.5) # + xlim(240,315) + ylim(0, 0.05)
  }
}
p = p + geom_density(data = wcr,aes(x=precipitation_flux, y=..density.., color = 'Observed'), size=1) 
p = p +   theme(legend.title = element_blank(), 
                legend.position = c(.85, .85), 
                legend.text = element_text(size = 11),
                legend.key = element_rect(fill = 'springgreen')) + 
  guides(colour = guide_legend(override.aes = list(size=2, stroke=1, fill=c('lightblue', 'orange', 'black')))) 
print(p)
dev.off()



png("~/Documents/metdownscaling/plots/validate/combined_prec_densities.png", width = 10, height = 6, units = "in",res = 200)

p1 = ggplot(nldas) + 
  geom_density(data = nldas, aes(x=nldas_prec, y=..density.., color = 'NLDAS Raw'), adjust = 10, size=0.5)  + xlab("Precipitation (mm/hr)") + ggtitle("2006 Gaussian Precipitation Densities") + 
  scale_color_manual(values = c('Observed' = 'black', 'Debiased Ensemble Members'='dodgerblue4', 'NLDAS Raw'='orange')) +
  scale_shape_manual(labels = c("Observed", "Debiased Ensemble Members", "NLDAS Raw"),
                     values = c(1, 1, 1))

for (d in debias_seq){
  for (e in ens_seq){
    tem = prec.df[paste0("prec",d, ".", e)]
    colnames(tem) = "v1"
    p1 = p1 + geom_density(data = tem, aes(x=v1, y=..density.., color='Debiased Ensemble Members'), size=0.5, adjust = 10, fill="lightblue", alpha=0.5) # + xlim(240,315) + ylim(0, 0.05)
  }
}
p1 = p1 + geom_density(data = wcr, aes(x=precipitation_flux, y=..density.., color = 'Observed'),adjust=10, size=0.5) 
p1 = p1 +   theme(legend.title = element_blank(), 
                legend.position = c(.78, .9), 
                legend.text = element_text(size = 10)) + 
  guides(colour = guide_legend(override.aes = list(size=1, stroke=1, fill=c('lightblue', 'orange', 'black')),
                               keywidth = .7, keyheight = .7)) 

#### HOT HOT HOT HOT HOT HOT HOT HOT HOT HOT HOT 
p2 = ggplot(nldas) + 
  geom_density(data = nldas, aes(x=nldas_prec, y=..density.., color = 'orange'), color='orange', adjust=1/10, size=0.5)  + xlab("Precipitation (mm/hr)")+ coord_cartesian(xlim=c(.1,1), ylim=c(0,2))# + 
  #scale_color_manual(values = c('Observed' = 'black', 'Debiased Ensemble Members'='dodgerblue4', 'NLDAS Raw'='orange'))
for (d in debias_seq){
  for (e in ens_seq){
    tem = prec.df[paste0("prec",d, ".", e)]
    colnames(tem) = "v1"
    p2 = p2 + geom_density(data = tem, aes(x=v1, y=..density.., color='dodgerblue4'),  color='dodgerblue4', adjust=1/10, fill="lightblue", alpha=0.5) # + xlim(240,315) + ylim(0, 0.05)
  }
}
p2 = p2 + geom_density(data = wcr, aes(x=precipitation_flux, y=..density.., colour="black"), color='black',  adjust=1/10, size=1)
p2 = p2 + geom_density(data = nldas, aes(x=nldas_prec, y=..density.., color = 'orange'), color='orange', adjust=1/10,size=1) + theme(legend.position="none")
# COLD COLD COLD COLD COLD COLD
p3 = ggplot(nldas) + 
  geom_density(data = nldas, aes(x=nldas_prec, y=..density.., color = 'orange'), color='orange', adjust=1/10, size=0.5)  + xlab("Precipitation (mm/hr)")+ coord_cartesian(xlim=c(1,14), ylim=c(0,0.4))# + 
#scale_color_manual(values = c('Observed' = 'black', 'Debiased Ensemble Members'='dodgerblue4', 'NLDAS Raw'='orange'))
for (d in debias_seq){
  for (e in ens_seq){
    tem = prec.df[paste0("prec",d, ".", e)]
    colnames(tem) = "v1"
    p3 = p3 + geom_density(data = tem, aes(x=v1, y=..density.., color='dodgerblue4'),  color='dodgerblue4', adjust=1/10, fill="lightblue", alpha=0.5, size=0.5) # + xlim(240,315) + ylim(0, 0.05)
  }
}
p3 = p3 + geom_density(data = wcr, aes(x=precipitation_flux, y=..density.., colour="black"), color='black',  adjust=1/10, size=0.5)
p3 = p3 + geom_density(data = nldas, aes(x=nldas_prec, y=..density.., color = 'orange'), color='orange', adjust=1/10,size=0.5) + theme(legend.position="none")

lay <- rbind(c(1,1,1,2,2),
             c(1,1,1,3,3))
grid.arrange(p1,p2,p3, layout_matrix = lay)
dev.off()

