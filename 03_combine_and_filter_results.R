# Stephanie.Clay@dfo-mpo.gc.ca
# 2025-08-18

# Combine gaussian fit data from every year and filter results.

rm(list=ls())
library(dplyr)
library(pbapply)
library(fst)
library(data.table)
library(oceancolouR)

years <- 1998:2024

# file with per-pixel % dineof
dineof_file <- "dineof_stats/dineof_percent_reconstructed.fst"

output_file <- "gaussian_fits/fit_data_combined_filtered.fst"

# output columns for metadata and columns specific to the gaussian fit
output_cols1 <- c("bin","longitude","latitude","year","NRMSE_bloom","t[start]","t[duration]","Amplitude[real]","Magnitude[real]","B0","h","sigma","beta","Flags")
output_cols1_newnames <- c("bin","longitude","latitude","year","NRMSE_bloom","t_start","t_duration","amplitude_real","magnitude_real","B0","h","sigma","beta","Flags")

# output columns not affected by the gaussian fit
output_cols2 <- c("bin","year","percent_dineof","annual_mean")


#*******************************************************************************
# MERGE METRICS AND ERROR METRICS FOR ALL YEARS

df <- pblapply(years, function(y) {
    read_fst(paste0("gaussian_fits/fit_data_",y,"_binned.fst")) %>%
        dplyr::mutate(year=y) %>%
        dplyr::select(bin,longitude,latitude,year,Flags,`t[start]`:annual_median) %>%
        dplyr::left_join(y=read_fst(gsub(".fst",paste0("_",y,".fst"),dineof_file)),by=c("bin","year"))
}) %>%
    do.call(what=dplyr::bind_rows) %>%
    dplyr::arrange(year,bin) %>%
    # FILTER BINS THAT EITHER HAVE NO VALID DINEOF DATA, OR ONE VALUE REPEATED OVER THE WHOLE YEAR
    dplyr::filter(is.finite(dineof_range) & dineof_range>0)

# get gaussian fit columns and filter them
df1 <- df %>%
    # valid chl-a in the input poly4 files ranges from 0 to 100 mg/m3, so this is a reasonable cutoff for the following variables
    dplyr::filter(between(`Amplitude[fit]`,0,100) &
                  between(`Amplitude[real]`,0,100) &
                  between(RMSE,0,100) &
                  between(RMSE_bloom,0,100)) %>%
    # add an extra metric for analysis
    dplyr::mutate(NRMSE_bloom = RMSE_bloom/`Amplitude[real]`) %>%
    dplyr::filter(NRMSE_bloom<=1) %>%
    dplyr::select(all_of(output_cols1)) %>%
    setNames(output_cols1_newnames)

# get other columns
df2 <- df %>% dplyr::select(all_of(output_cols2))

df <- dplyr::full_join(df1,df2,by=c("bin","year"))

# WRITE FILTERED RESULTS TO FST
write_fst(df, path=output_file, compress=100)

