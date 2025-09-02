# Stephanie.Clay@dfo-mpo.gc.ca
# 2025-08-18

# Calculate climatologies and standardized anomalies for each metric of interest.

rm(list=ls())
library(dplyr)
library(lubridate)
library(oceancolouR)
library(fst)

# years to use in the climatology
clim_years <- 1998:2024

input_file <- "gaussian_fits/fit_data_combined_filtered.fst"

output_file <- "gaussian_fits/fit_data_combined_and_filtered_with_clim.fst"

metrics_for_clim <- c("t_start","t_duration","amplitude_real","magnitude_real","annual_mean") 


#*******************************************************************************

df <- read_fst(input_file) %>%
    tidyr::pivot_longer(cols=all_of(metrics_for_clim), names_to="metric") %>%
    dplyr::group_by(bin,metric) %>%
    dplyr::mutate(climatology_mean=mean(value[year %in% clim_years],na.rm=TRUE),
                  climatology_sd=sd(value[year %in% clim_years],na.rm=TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(anomaly=value-climatology_mean) %>%
    dplyr::mutate(standardized_anomaly=anomaly/climatology_sd) %>%
    dplyr::group_by(bin) %>%
    dplyr::mutate(valid_fits=sum(is.finite(NRMSE_bloom[metric=="t_start"]))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(year,bin,metric)

write_fst(df, path=output_file, compress=100)

