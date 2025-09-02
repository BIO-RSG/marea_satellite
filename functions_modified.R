
# EXCERPT FROM PHYTOFIT FUNCTIONS.R

# custom functions from oceancolouR package
plus_minus <- function (x, n) return((x - n):(x + n))
find_line <- function (x1, y1, x2, y2) {
  m <- (y2 - y1)/(x2 - x1)
  b <- y1 - m * x1
  return(list(slope=m, intercept=b))
}

# chl_to_use is the mean or the median (logged or not), depending on user selection - it's not subset
bf_data_calc <- function(composite, chl_to_use, ind_dayrange_percov, ind_percov, ind_dayrange,
                         daily_percov, t_range, log_chla, doy_vec, sv) {
  
  nofit_msg <- NULL
  yfit <- ybkrnd <- rep(NA, sum(ind_dayrange))
  
  # make a dataframe of input vectors
  # roc and thresh background lines, loess, and the scatterplot points (sized by percent coverage) use ind_percov
  # roc and thresh models, and gauss use ind_dayrange_percov (all chl within the day range, with sufficient percent coverage)
  # yfit and ybkrnd use ind_dayrange
  df_input <- data.frame(yday=doy_vec, var=chl_to_use, var_to_fit=chl_to_use) %>% dplyr::mutate(weight=1)
  if (sv$use_weights) {df_input$weight <- daily_percov}
  
  # make a dataframe for output
  df_output <- data.frame(doy=doy_vec, percent_coverage=daily_percov, bfy=df_input$var) %>% dplyr::mutate(model=NA,background=NA,loess=NA)
  
  # use loess smooth on the real data points
  if (sv$smoothMethod == 'loess'){
    mod <- try(loess(var_to_fit ~ yday, data=df_input[ind_percov,], weights=df_input$weight[ind_percov], span=sv$loessSpan, degree=2), silent=TRUE)
    if (!(class(mod)=="try-error" | is.null(mod))) {
      df_output$loess[ind_percov] <- df_input$var_to_fit[ind_percov] <- fitted(mod) # use loess for fitting instead of real values
      df_input$weight <- 1 # reset weights so they aren't used again in a gaussian fit
    }
  }
  
  if (sv$fitmethod == 'gauss') {
    
    if (sv$tm) {tmp_ti_lim <- c(1,365)
    } else {tmp_ti_lim <- sv$ti_limits}
    gauss_res <- gaussFit(dfin=df_input[ind_dayrange_percov,],
                          bloomShape=sv$bloomShape,
                          tm=sv$tm, beta=sv$beta,
                          tm_limits=sv$tm_limits,
                          ti_limits=tmp_ti_lim,
                          t_range=t_range,
                          log_chla=log_chla,
                          composite=composite,
                          flag1_lim1=sv$flag1_lim1,
                          flag1_lim2=sv$flag1_lim2,
                          flag2_lim1=sv$flag2_lim1,
                          flag2_lim2=sv$flag2_lim2,
                          ti_threshold=sv$ti_threshold,
                          tt_threshold=sv$tt_threshold,
                          ydays_dayrange=df_input$yday[ind_dayrange],
                          rm_bkrnd=sv$rm_bkrnd,
                          ti_threshold_type=sv$ti_threshold_type,
                          ti_threshold_constant=sv$ti_threshold_constant)
    fitparams <- gauss_res$values
    # calculate rmse if non-null fit exists
    if (is.null(gauss_res$fit)) {
      fitparams$RMSE <- NA
      fitparams$RMSLE <- NA
      nofit_msg <- gauss_res$nofit_msg
    } else {
      rmsedf <- data.frame(x=df_input$var[ind_dayrange_percov],
                           y=predict(gauss_res$fit, newdata=list(t=df_input$yday[ind_dayrange_percov])))
      if (log_chla) {rmsedf <- rmsedf %>% dplyr::mutate(x=10^x,  y=10^y)}
      fitparams$RMSE <- sqrt(mean((rmsedf$x - rmsedf$y)^2, na.rm = TRUE))
      rmsedf <- rmsedf %>% dplyr::mutate(x=log10(x),  y=log10(y)) %>% dplyr::filter(is.finite(x) & is.finite(y))
      fitparams$RMSLE <- sqrt(mean((rmsedf$x - rmsedf$y)^2, na.rm = TRUE))
      yfit <- gauss_res$yfit
      ybkrnd <- gauss_res$ybkrnd
    }
    # round tmax_fit day
    fitparams[,"t[max_fit]"] <- round(as.numeric(fitparams[,"t[max_fit]"]))
  
  } else {
    
    if (sv$fitmethod == 'roc') {
      fitparams <- rateOfChange(dfin=df_input[ind_dayrange_percov,],
                                dfin_all=df_input[ind_percov,],
                                 bloomShape=sv$bloomShape,
                                 tm_limits=sv$tm_limits,
                                 ti_limits=sv$ti_limits,
                                 log_chla=log_chla,
                                 rm_bkrnd=sv$rm_bkrnd,
                                use_weights=sv$use_weights)
    } else if (sv$fitmethod == "thresh") {
      fitparams <- threshold(dfin=df_input[ind_dayrange_percov,],
                             dfin_all=df_input[ind_percov,],
                              threshcoef=sv$threshcoef, 
                              bloomShape=sv$bloomShape,
                              tm_limits=sv$tm_limits,
                              ti_limits=sv$ti_limits,
                              log_chla=log_chla,
                              rm_bkrnd=sv$rm_bkrnd,
                             use_weights=sv$use_weights)
    }
    
    nofit_msg <- fitparams$nofit_msg
    modelfit <- fitparams$yfit
    if (!is.null(modelfit)) {
      ybkrnd <- predict(modelfit,newdata=data.frame(tall=df_input$yday[ind_dayrange]))
    }
    fitparams <- fitparams$values
    
  }
  
  # round more days
  tvars <- c("t[start]", "t[max_real]", "t[end]", "t[duration]")
  fitparams[,tvars] <- round(as.numeric(fitparams[,tvars]))
  # add extra summary stats 
  bf_extra <- data.frame(Annual_Mean = mean(chl_to_use, na.rm=TRUE),
                         Annual_Median = median(chl_to_use, na.rm=TRUE),
                         Annual_StDev = sd(chl_to_use, na.rm=TRUE))
  fitparams <- dplyr::bind_cols(bf_extra,fitparams)
  
  if (log_chla) {
    df_output$bfy <- 10^df_output$bfy
    df_output$loess <- 10^df_output$loess
    yfit <- 10^yfit
    ybkrnd <- 10^ybkrnd
    fitparams$Annual_Mean <- 10^fitparams$Annual_Mean
    fitparams$Annual_Median <- 10^fitparams$Annual_Median
    fitparams$Annual_StDev <- 10^fitparams$Annual_StDev
  }
  
  # after all the modelling and transformations, add to output dataframe
  df_output$model[ind_dayrange] <- yfit
  df_output$background[ind_dayrange] <- ybkrnd
  
  return(list(fitparams=fitparams, nofit_msg=nofit_msg, df_final=df_output))
  
}
