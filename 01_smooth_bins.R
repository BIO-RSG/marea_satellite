# Stephanie.Clay@dfo-mpo.gc.ca
# 2025-08-10

# For a given bin, calculate a filtered value of all the bins in the surrounding area.

rm(list=ls())
library(terra)
library(dplyr)
library(lubridate)
library(minpack.lm)
library(pbapply)
library(oceancolouR)
library(fst)

years <- 1998:2024

output_file <- "smoothed_bins/bin_from_filtered_20km_median" # excluding year and extension

max_dist <- 20000 # metres from central bin

# these are the default settings used in phytofit azmp/azomp fitting - use this to replace a single bin value with a filtered summary value of all the bins within 20km (max_dist)
percent <- 20 # required percent coverage within 20km (ignored if there aren't enough pixels)
outlier <- "sd3" # sd1, sd2, or sd3
dailystat <- "average"


#*******************************************************************************

dir.create(dirname(output_file), showWarning=FALSE, recursive=TRUE)

bindf <- get_bins(region="nwa",resolution="4km",variable=c("bin","longitude","latitude")) %>% dplyr::arrange(bin)

# function to replace a bin with a filtered summary value of surrounding bins within a radius defined by max_dist
bin_summary_fn <- function(i) {
    alllat <- as.numeric(r[,3])
    alllon <- as.numeric(r[,2])
    midlat <- alllat[i]
    midlon <- alllon[i]
    # reduce to bins within 1.5 degrees for the sake of speed, then calculate distance and reduce to bins within max_dist km
    inds <- abs(alllat-midlat)<1.5 & abs(alllon-midlon)<1.5
    dists <- geodist::geodist(data.frame(longitude=midlon,latitude=midlat),
                              data.frame(longitude=alllon[inds], latitude=alllat[inds]),
                              measure="geodesic")
    tmpr <- r[inds,][dists<max_dist,][,4:368]
    # log the chla data
    tmpr <- log10(tmpr)
    tmpr[!is.finite(tmpr)] <- NA
    # get the average and sd to define outlier boundaries
    tmpr_mean <- colMeans(tmpr,na.rm=TRUE)
    tmpr_sd <- apply(tmpr,MARGIN=2,FUN=sd,na.rm=TRUE)
    sdfactor <- as.numeric(gsub("sd","",outlier))
    lim1 <- tmpr_mean - sdfactor*tmpr_sd
    lim2 <- tmpr_mean + sdfactor*tmpr_sd
    tmpr <- lapply(1:ncol(tmpr), function(j) {
        tmprj <- tmpr[,j]
        tmprj[tmprj < lim1[j] | tmprj > lim2[j]] <- NA
        return(tmprj)
    }) %>% do.call(what=cbind)
    # check percent coverage within the 20km radius on each day
    tmpr_pcov <- 100*colSums(is.finite(tmpr))/nrow(tmpr) >= percent
    tmpr[,!tmpr_pcov] <- NA
    # calculate average or median
    if (dailystat=="average") {
        tmpr <- matrix(colMeans(tmpr,na.rm=TRUE), nrow=1)
    } else if (dailystat=="median") {
        tmpr <- matrix(apply(tmpr,MARGIN=2,FUN=median,na.rm=TRUE), nrow=1)
    }
    # convert back to linear space and return result
    tmpr <- 10^tmpr
    return(tmpr)
}


for (year in years) {
    
    cat(year,"...\n")
    
    if (year==(year(Sys.Date())-1)) {
        r <- read_fst(paste0("/mnt/data3/claysa/ESA_L3b/v6.0/occci/CHL_POLY4/NWA/annual_fst_DINEOFfilled/NWA_occci_filled2yrs_daily_CHL_POLY4_",year,".fst"))
    } else {
        r <- read_fst(paste0("/mnt/data3/claysa/ESA_L3b/v6.0/occci/CHL_POLY4/NWA/annual_fst_DINEOFfilled/NWA_occci_filled3yrs_daily_CHL_POLY4_",year,".fst"))
    }
    
    r <- r %>%
        # first remove any bins where there is only one unique value all year
        dplyr::group_by(bin) %>%
        dplyr::mutate(n=length(unique(value))) %>%
        dplyr::ungroup() %>%
        dplyr::filter(n>1) %>%
        # join and reshape so it's a matrix with row = bin, column = day of year, for easier processing
        # note that this forces the number of days to 365, ignoring leap years
        dplyr::right_join(y=expand.grid(bin=bindf$bin,doy=1:365), by=c("bin","doy")) %>%
        dplyr::arrange(doy,bin)
    r <- matrix(r$value,nrow=nrow(bindf))
    r <- cbind(as.matrix(bindf),r)
    
    newdf <- pblapply(1:nrow(r), FUN=bin_summary_fn, cl=10) %>% do.call(what=rbind)
    
    newdf <- data.frame(bin=rep(bindf$bin,365),
                        doy=rep(1:365,each=nrow(bindf)),
                        value=as.numeric(c(newdf))) %>%
        tidyr::drop_na(value)
    
    write_fst(newdf, path=paste0(output_file, "_", year, "_binned.fst"), compress=100)
    
    gc()
    
}


