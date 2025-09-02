# Stephanie.Clay@dfo-mpo.gc.ca
# Aug 2025

# Make pngs of values, climatologies, and standardized anomalies for manual inspection.

rm(list=ls())
library(dplyr)
library(oceancolouR)
library(terra)
library(fst)
library(sf)
library(pbapply)
library(lubridate)
library(patchwork)
library(ggplot2)

years <- 1998:2024

input_file <- "gaussian_fits/fit_data_combined_and_filtered_with_clim.fst"

output_file <- "gaussian_fits/output_png/gaussian_fit_metrics"


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

metrics <- c("t_start","t_duration","amplitude_real","magnitude_real","annual_mean")
error_metrics <- c("NRMSE_bloom","percent_dineof")


#*******************************************************************************
# CLIMATOLOGY
# Layers: t_start, t_duration, amplitude_real, magnitude_real, annual_mean, valid_fits

# make the layers across the whole time series
dfclim <- df %>% dplyr::distinct(bin,metric,climatology_mean,climatology_sd)
rclimmean <- pblapply(1:length(metrics), function(i) {
    tmpr <- terra::rast(bin_to_raster(dfclim %>% dplyr::filter(metric==metrics[i]) %>% dplyr::distinct(bin,climatology_mean), ext=c(xlim,ylim)))
    tmpr[!is.finite(tmpr)] <- NA
    return(tmpr)
}) %>% do.call(what=c) %>% setNames(paste0(metrics,"_climatology_mean"))
rclimsd <- pblapply(1:length(metrics), function(i) {
    tmpr <- terra::rast(bin_to_raster(dfclim %>% dplyr::filter(metric==metrics[i]) %>% dplyr::distinct(bin,climatology_sd), ext=c(xlim,ylim)))
    tmpr[!is.finite(tmpr)] <- NA
    return(tmpr)
}) %>% do.call(what=c) %>% setNames(paste0(metrics,"_climatology_sd"))
dfvf <- df %>% dplyr::distinct(bin,valid_fits)
rvf <- terra::rast(bin_to_raster(dfvf, ext=c(xlim,ylim)))
rvf[!is.finite(rvf)] <- NA
r <- c(rclimmean,rclimsd,rvf)

m <- pblapply(r, function(x) {
    if (startsWith(names(x),"amplitude_real")|startsWith(names(x),"magnitude_real")|startsWith(names(x),"annual_mean")) {
        tmpm <- make_raster_map(x, title=names(x), xlim=xlim, ylim=ylim, trans="log10", labels = scales::format_format(scientific=FALSE), na.value="lightgrey", map_fill="darkgrey", map_colour="black", map_linewidth=0.1, hires_land=TRUE)
    } else {
        tmpm <- make_raster_map(x, title=names(x), xlim=xlim, ylim=ylim, labels = scales::format_format(scientific=FALSE), na.value="lightgrey", map_fill="darkgrey", map_colour="black", map_linewidth=0.1, hires_land=TRUE)
    }
    tmpm <- tmpm +
        labs(fill="") +
        theme(legend.position="bottom",legend.direction="horizontal",plot.title=element_text(size=18)) +
        guides(fill = guide_colourbar(title.hjust = 0,
                                      ticks.colour = "black",
                                      barheight = unit(0.6, "cm"),
                                      barwidth = unit(8, "cm"),
                                      frame.colour = "black"))
    return(tmpm)
}) %>% wrap_plots(ncol=5)

ggsave(filename=paste0(output_file,"_climatology.png"),
       plot=m,
       dpi=150,
       units="px",
       width=3200,
       height=2900)



#**************************************************
# INDIVIDUAL YEARS
# Layers: value and standardized_anomaly for each metric, and the error_metrics

for (year in years) {

    cat(year,"...\n")
    
    tmp <- df %>% dplyr::filter(year==!!year)
    
    # rasterize the metrics
    r <- pblapply(1:length(metrics), function(i) {
        tmpr <- terra::rast(bin_to_raster(tmp %>% dplyr::filter(metric==metrics[i]) %>% dplyr::distinct(bin,value), ext=c(xlim,ylim)))
        tmpr[!is.finite(tmpr)] <- NA
        # add the time to the spatraster
        time(tmpr) <- rep(as_datetime(paste0(year,"0101")),nlyr(tmpr))
        return(tmpr)
    }) %>% do.call(what=c)
    
    rsa <- pblapply(1:length(metrics), function(i) {
        tmpr <- terra::rast(bin_to_raster(tmp %>% dplyr::filter(metric==metrics[i]) %>% dplyr::distinct(bin,standardized_anomaly), ext=c(xlim,ylim)))
        tmpr[!is.finite(tmpr)] <- NA
        # add the time to the spatraster
        time(tmpr) <- rep(as_datetime(paste0(year,"0101")),nlyr(tmpr))
        return(tmpr)
    }) %>% do.call(what=c)
    
    re <- pblapply(1:length(error_metrics), function(i) {
        tmpr <- terra::rast(bin_to_raster(tmp %>% dplyr::select(bin,all_of(error_metrics[i])) %>% dplyr::distinct(), ext=c(xlim,ylim)))
        tmpr[!is.finite(tmpr)] <- NA
        # add the time to the spatraster
        time(tmpr) <- rep(as_datetime(paste0(year,"0101")),nlyr(tmpr))
        return(tmpr)
    }) %>% do.call(what=c)
    
    r <- c(r,rsa,re)
    names(r) <- c(metrics,paste0(metrics,"_standardized_anomaly"),error_metrics)
    
    m <- pblapply(r, function(x) {
        if (endsWith(names(x),"standardized_anomaly")) {
            tmpm <- make_raster_map(x, title=names(x), xlim=xlim, ylim=ylim, labels = scales::format_format(scientific=FALSE), col_limits=c(-3,3), cm=c("#0000FF","#5555FF","#AAAAFF","#FFFFFF","#FFAAAA","#FF5555","#FF0000"), set_extremes = TRUE, na.value="lightgrey", map_fill="darkgrey", map_colour="black", map_linewidth=0.1, hires_land=TRUE)
        } else {
            if (names(x) %in% c("amplitude_real","magnitude_real","annual_mean")) {
                tmpm <- make_raster_map(x, title=names(x), xlim=xlim, ylim=ylim, trans="log10", labels = scales::format_format(scientific=FALSE), na.value="lightgrey", map_fill="darkgrey", map_colour="black", map_linewidth=0.1, hires_land=TRUE)
            } else if (names(x)=="NRMSE_bloom") {
                tmpm <- make_raster_map(x, title=names(x), xlim=xlim, ylim=ylim, col_limits=c(0.1,0.4), set_extremes=TRUE, labels = scales::format_format(scientific=FALSE), na.value="lightgrey", map_fill="darkgrey", map_colour="black", map_linewidth=0.1, hires_land=TRUE)
            } else {
                tmpm <- make_raster_map(x, title=names(x), xlim=xlim, ylim=ylim, labels = scales::format_format(scientific=FALSE), na.value="lightgrey", map_fill="darkgrey", map_colour="black", map_linewidth=0.1, hires_land=TRUE)
            }
        }
        tmpm <- tmpm +
            labs(fill="") +
            theme(legend.position="bottom",legend.direction="horizontal",plot.title=element_text(size=18)) +
            guides(fill = guide_colourbar(title.hjust = 0,
                                          ticks.colour = "black",
                                          barheight = unit(0.6, "cm"),
                                          barwidth = unit(8, "cm"),
                                          frame.colour = "black"))
        return(tmpm)
    }) %>% wrap_plots(nrow=3,ncol=5)
    
    ggsave(filename=paste0(output_file,"_",year,".png"),
           plot=m,
           dpi=150,
           units="px",
           width=3400,
           height=3000)
    
}
