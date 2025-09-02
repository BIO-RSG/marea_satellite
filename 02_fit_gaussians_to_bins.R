# Stephanie.Clay@dfo-mpo.gc.ca
# 2025-04-11

# For a given bin and year, fit a gaussian to the daily time series.

rm(list=ls())
library(terra)
library(dplyr)
library(lubridate)
library(minpack.lm)
library(pbapply)
library(oceancolouR)
library(fst)
library(sf)
library(patchwork)
source("functions_modified.R")
source("gaussFit_modified.R")

years <- 1998:2024

input_file <- "smoothed_bins/bin_from_filtered_20km_median"
output_file <- "gaussian_fits/fit_data" # excluding year and extension


#*******************************************************************************

# FIND WHICH BINS ARE IN YOUR REGION OF INTEREST, FOR FILTERING
# without this line you'll have a curve at the bottom of the plate caree projection
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
binsf <- binlatlon(lonlim=xlim,latlim=ylim) %>% st_as_sf(coords=c("longitude","latitude"), crs=st_crs(coords)) %>% st_as_sf()
coords$bins <- list(st_intersection(binsf,coords[1,])$bin)

bindf <- get_bins(region="nwa",resolution="4km",variable=c("bin","longitude","latitude")) %>% 
    dplyr::filter(bin %in% coords$bins[[1]]) %>%
    dplyr::arrange(bin)

# GAUSSIAN FIT SETTINGS
composite <- 1 # daily
fitmethod <- "gauss"
bloomShape <- "symmetric"
smoothMethod <- "loess"
log_chla <- TRUE
loessSpan <- 0.2
t_range <- c(1,220)
ti_limits <- c(1,176)
tm_limits <- c(60,180)
ti_threshold_type <- "percent_thresh"
ti_threshold <- tt_threshold <- 0.2
ti_threshold_constant <- 0.1
threshcoef <- 1.05
tm <- TRUE
beta <- TRUE
use_weights <- TRUE
rm_bkrnd <- TRUE
flag1_lim1 <- 0.75
flag1_lim2 <- 1.25
flag2_lim1 <- 0.85
flag2_lim2 <- 1.15
sv <- list(use_weights,smoothMethod,loessSpan,fitmethod,bloomShape,tm, beta,tm_limits,ti_limits,threshcoef,flag1_lim1,flag1_lim2, flag2_lim1,flag2_lim2,ti_threshold,tt_threshold,rm_bkrnd, ti_threshold_type,ti_threshold_constant)
names(sv) <- c("use_weights","smoothMethod","loessSpan","fitmethod","bloomShape","tm","beta","tm_limits", "ti_limits","threshcoef","flag1_lim1","flag1_lim2","flag2_lim1","flag2_lim2","ti_threshold", "tt_threshold","rm_bkrnd","ti_threshold_type","ti_threshold_constant")
doy_vec <- 1:365
final_colnames <- c("Annual_Mean","Annual_Median","Annual_StDev","t[start]","t[max_real]","t[max_fit]","t[end]","t[duration]","Magnitude[real]","Magnitude[fit]","Amplitude[real]","Amplitude[fit]","Flags","B0", "h", "sigma", "beta", "failure_code","RMSE","RMSLE","RMSE_bloom","RMSLE_bloom")


tmp_gauss_fn <- function(i) {
    chl_to_use <- r[i,]
    
    
    # # remove outliers
    # limit <- quantile(chl_to_use,probs=0.99,na.rm=TRUE)
    # chl_to_use[chl_to_use>limit] <- NA
    
    
    # fit gaussian to time series
    if (sum(is.finite(chl_to_use))==0 | length(unique(chl_to_use[is.finite(chl_to_use)]))==1) {
        final <- matrix(nrow=1,ncol=length(final_colnames))
        colnames(final) <- final_colnames
        return(final)
    }
    if (log_chla) {chl_to_use <- log10(chl_to_use)}
    chl_to_use[!is.finite(chl_to_use)] <- NA
    ind_dayrange <- rep(FALSE,365)
    ind_dayrange[t_range[1]:t_range[2]] <- TRUE
    daily_percov <- rep(0,365)
    daily_percov[is.finite(chl_to_use)] <- 100
    ind_percov <- daily_percov==100
    ind_dayrange_percov <- ind_dayrange & ind_percov
    fitlist <- bf_data_calc(composite, chl_to_use, ind_dayrange_percov, ind_percov, ind_dayrange,
                            daily_percov, t_range, log_chla, doy_vec, sv)
    if (is.null(fitlist$fitparams)) {
        final <- matrix(nrow=1,ncol=length(final_colnames))
        colnames(final) <- final_colnames
    } else {
        # ADD SOME METRICS MEASURED FROM WITHIN THE BLOOM PERIOD
        final <- fitlist$fitparams
        start <- final$`t[start]`
        end <- final$`t[end]`
        bloom_ind <- between(fitlist$df_final$doy,start,end)
        rmsedf <- fitlist$df_final[bloom_ind,]
        final$RMSE_bloom <- sqrt(mean((rmsedf$bfy - rmsedf$model)^2, na.rm = TRUE))
        rmsedf <- rmsedf %>% dplyr::filter(bfy>0&model>0) %>% dplyr::mutate(bfy=log10(bfy),model=log10(model))
        final$RMSLE_bloom <- sqrt(mean((rmsedf$bfy - rmsedf$model)^2, na.rm = TRUE))
        final <- as.matrix(final)
    }
    return(final)
}



for (year in years) {
    
    cat(year,"...\n")
    
    r <- read_fst(paste0(input_file, "_", year, "_binned.fst")) %>%
        dplyr::group_by(bin) %>%
        dplyr::mutate(n=length(unique(value))) %>%
        dplyr::ungroup() %>%
        dplyr::filter(n>1) %>%
        dplyr::right_join(y=expand.grid(bin=bindf$bin,doy=1:ifelse(leap_year(year),366,365)), by=c("bin","doy")) %>%
        dplyr::arrange(doy,bin)
    r <- matrix(r$value,nrow=nrow(bindf))
    r <- r[,1:365]
    
    newdf <- pblapply(1:nrow(r), FUN=tmp_gauss_fn, cl=10) %>%
        do.call(what=rbind) %>%
        as.data.frame()
    
    annstats <- data.frame(annual_mean=rowMeans(r,na.rm=TRUE),
                           annual_sd=apply(r,MARGIN=1,FUN=sd,na.rm=TRUE),
                           annual_geometric_mean=apply(r,MARGIN=1,FUN=geoMean),
                           annual_median=apply(r,MARGIN=1,FUN=median,na.rm=TRUE))
    
    newdf <- dplyr::bind_cols(bindf,newdf,annstats) %>%
        dplyr::select(bin,longitude,latitude,`t[start]`:annual_median)
    
    write_fst(newdf, path=paste0(output_file, "_", year, "_binned.fst"), compress=100)
    
    gc()
    
}


