# Stephanie.Clay@dfo-mpo.gc.ca
# Aug 2025

# Read filtered fits, map each metric, and write to netCDF.

rm(list=ls())
library(dplyr)
library(oceancolouR)
library(terra)
library(fst)
library(sf)
library(pbapply)
library(lubridate)
library(ncdf4)
source("write_netcdf_marea.R")

years <- 1998:2024

input_file <- "gaussian_fits/fit_data_combined_and_filtered_with_clim.fst"

output_file <- "gaussian_fits/output_nc/gaussian_fit_metrics"


#*******************************************************************************

dir.create(dirname(output_file), showWarnings=FALSE)

df <- read_fst(input_file)

# without this line you'll have a curve at the bottom of the plate-caree projection
sf_use_s2(FALSE)

# cut out Hudson Bay/Strait, Northwestern Passages, Baffin Bay
coords <- data.frame(latitude=c(39,46,60,63,65,66,66,39,39),
                     longitude=c(-76,-76,-65,-65,-64,-62.5,-42,-42,-76)) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    summarise(geometry = st_combine(geometry)) %>%
    st_cast("POLYGON")

bbox <- st_bbox(coords)
xlim <- bbox[c(1,3)]
ylim <- bbox[c(2,4)]

bins <- binlatlon(lonlim=xlim,latlim=ylim)
binsf <- bins %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(coords)) %>% 
    st_as_sf()

coords$bins <- list(st_intersection(binsf,coords[1,])$bin)

df <- df %>% dplyr::filter(bin %in% coords$bins[[1]])

ncinfo <- list(t_start=list(units="",
                            long_name="Day of year of the start of the spring phytoplankton bloom",
                            range=c(1,176)), # ti_limits boundaries used in gaussian fitting
               t_duration=list(units="days",
                               long_name="Duration of the spring phytoplankton bloom",
                               range=c(1,220)), # t_range boundaries used in gaussian fitting
               amplitude_real=list(units="mg/m3",
                                   long_name="Maximum concentration during the spring phytoplankton bloom period",
                                   range=c(0,100)),
               magnitude_real=list(units="days*mg/m3",
                                   long_name="Total chlorophyll-a produced during the spring phytoplankton bloom period",
                                   range=c(0,22000)), # if t.start=1, t.end=220, and chla=100 every day
               annual_mean=list(units="mg/m3",
                                long_name="Average chlorophyll-a over the year",
                                standard_name="mass_concentration_of_chlorophyll_a_in_sea_water",
                                range=c(0,100)),
               NRMSE_bloom=list(units="mg/m3",
                                long_name="Root mean square error between the fitted Gaussian and real chl-a values during the bloom period, normalized to amplitude_real",
                                range=c(0,1)),
               percent_dineof=list(units="",
                                   long_name="Percentage of days with pixel values estimated using DINEOF",
                                   range=c(0,100)))


#**************************************************
# Layers: t.start, t.duration, amplitude.real, magnitude.real, annual.mean, error metrics

metrics <- c("t_start","t_duration","amplitude_real","magnitude_real","annual_mean")
error_metrics <- c("NRMSE_bloom","percent_dineof")

for (year in years) {
    
    cat(year,"...\n")
    
    tmp <- df %>% dplyr::filter(year==!!year)
    
    # rasterize the metrics
    r <- pblapply(1:length(metrics), function(i) {
        tmpr <- terra::rast(bin_to_raster(tmp %>% dplyr::filter(metric==metrics[i]) %>% dplyr::distinct(bin,value), ext=c(xlim,ylim)))
        tmpr[!is.finite(tmpr)] <- NA
        # add the time to the spatraster
        time(tmpr) <- rep(as_datetime(paste0(year,"0101")),nlyr(tmpr))
        # time(tmpr) <- rep(year,nlyr(tmpr))
        return(tmpr)
    }) %>% do.call(what=c)
    
    re <- pblapply(1:length(error_metrics), function(i) {
        tmpr <- terra::rast(bin_to_raster(tmp %>% dplyr::select(bin,all_of(error_metrics[i])) %>% dplyr::distinct(), ext=c(xlim,ylim)))
        tmpr[!is.finite(tmpr)] <- NA
        # add the time to the spatraster
        time(tmpr) <- rep(as_datetime(paste0(year,"0101")),nlyr(tmpr))
        # time(tmpr) <- rep(year,nlyr(tmpr))
        return(tmpr)
    }) %>% do.call(what=c)
    
    r <- c(r,re)
    names(r) <- c(metrics,error_metrics)
    
    dim_lon <- ncdim_def(name="longitude",
                         longname="Longitude in decimal degrees",
                         units="degrees_east",
                         vals=xFromCol(r),
                         create_dimvar=TRUE,
                         unlim=FALSE)
    dim_lat <- ncdim_def(name="latitude",
                         longname="Latitude in decimal degrees",
                         units="degrees_north",
                         vals=yFromRow(r),
                         create_dimvar=TRUE,
                         unlim=FALSE)
    dim_time <- ncdim_def(name="time",
                          longname="time",
                          units="seconds since 1970-01-01T00:00:00Z",
                          # units="years since 0",
                          calendar="standard",
                          vals=as.numeric(time(r))[1],
                          create_dimvar=TRUE,
                          unlim=FALSE)
    outdims <- list(dim_lon,dim_lat,dim_time) # keep the dimensions in this order (x,y)
    
    ncvals <- lapply(1:nlyr(r), function(i) values(r[[i]])[,1])
    names(ncvals) <- names(r)
    
    new_file <- paste0(output_file,"_",year,".nc")
    
    # make list of variables
    nclist <- lapply(1:nlyr(r), function(i) {
        ncv <- names(r[[i]])
        ncvar_def(name=ncv, units=ncinfo[[ncv]]$units, dim=outdims, missval=NaN, shuffle=FALSE, compression=5, prec="float")
    })
    # create new output netcdf
    ncout <- nc_create(new_file, nclist, force_v4=TRUE)
    # add global attributes
    gatts <- make_att_list_mapped(new_file=basename(new_file),
                                  title=paste0("NWA OCCCI spring bloom metrics ",year),
                                  timeatt=year)
    dummy <- lapply(1:length(gatts), function(g) ncatt_put(nc=ncout,varid=0,attname=names(gatts)[g],attval=gatts[[g]]))
    # add extra attributes for some variables
    dummy <- lapply(1:nlyr(r), function(i) {
        ncv <- names(r[[i]])
        sn <- ncinfo[[ncv]]$standard_name
        if (!is.null(sn)) {ncatt_put(nc=ncout, varid=ncv, attname="standard_name", attval=sn)}
        ncatt_put(nc=ncout, varid=ncv, attname="long_name", attval=ncinfo[[ncv]]$long_name)
        ncatt_put(nc=ncout, varid=ncv, attname="valid_min", attval=ncinfo[[ncv]]$range[1], prec="float")
        ncatt_put(nc=ncout, varid=ncv, attname="valid_max", attval=ncinfo[[ncv]]$range[2], prec="float")
        ncatt_put(nc=ncout, varid=ncv, attname="ioos_category", attval="Ocean Color")
    })
    # add standard names to lat, lon, time
    ncatt_put(nc=ncout, varid="longitude", attname="standard_name", attval="longitude")
    ncatt_put(nc=ncout, varid="latitude", attname="standard_name", attval="latitude")
    ncatt_put(nc=ncout, varid="time", attname="standard_name", attval="time")
    # add ioos_category to lat, lon, time
    ncatt_put(nc=ncout, varid="longitude", attname="ioos_category", attval="Location")
    ncatt_put(nc=ncout, varid="latitude", attname="ioos_category", attval="Location")
    ncatt_put(nc=ncout, varid="time", attname="ioos_category", attval="Time")
    # put variable values in file
    dummy <- lapply(1:length(nclist), function(i) ncvar_put(ncout,nclist[[i]],ncvals[[i]]))
    # close file
    nc_close(ncout)
    
    # nc <- nc_open(paste0(output_file,"_",year,".nc"))
    # nc
    # test1 <- ncvar_get(nc,"t_start")
    # nc_close(nc)
    # test2 <- terra::rast(paste0(output_file,"_",year,".nc"), lyrs="t_start")
    
}
