# Stephanie.Clay@dfo-mpo.gc.ca
# 2025-05-26

# For a given pixel and year, determine the % of days that have DINEOF-estimated data as opposed to real satellite data.

rm(list=ls())
library(dplyr)
library(pbapply)
library(fst)
library(lubridate)

years <- 1998:2024

output_file <- "dineof_stats/dineof_percent_reconstructed.fst"


#*******************************************************************************

# for (y in years) {
#     
#     print(y)
# 
#     if (y<(year(Sys.Date())-1)) {
#         yearstr <- paste0(range(plus_minus(y,1)),collapse="-")
#     } else {
#         yearstr <- paste0((y-1),"-",y)
#     }
#     f <- paste0("/mnt/data1/claysa/DINEOF/OUTPUT_FILES/NWA_daily/forDineof_NWA_occci_CHL_POLY4_",yearstr,"_daily_logged/forDineof_NWA_occci_CHL_POLY4_",yearstr,"_daily_logged_formatted.fst")
#     
#     read_fst(f) %>%
#         dplyr::filter(year==y) %>%
#         dplyr::group_by(bin) %>%
#         dplyr::summarize(percent_dineof=100*sum(is.finite(var_imputed)&!is.finite(var))/n(),
#                          dineof_range=diff(range(var_imputed,na.rm=TRUE))) %>%
#         dplyr::ungroup() %>%
#         dplyr::mutate(year=y) %>%
#         write_fst(path=gsub(".fst",paste0("_",y,".fst"),output_file), compress=100)
#     
#     gc()
#     
# }


files <- paste0(gsub(".fst","_",output_file),years,".fst")

df <- pblapply(files, read_fst) %>% do.call(what=dplyr::bind_rows)

write_fst(df, path=output_file, compress=100)

