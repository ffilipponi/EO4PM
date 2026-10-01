# Title         : phenoEO4PM.R
# Description   : Define functions to estimate Key Phenological Stages using Gu method
# Date          : May 2026
# Version       : 2.1
# Copyright name: CC BY
# Licence       : GPL v3
# Authors       : Federico Filipponi
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
# ##############################################################################
#' @title phenoEO4PM
#'
#' @description Calculate phenological metrics from daily time series
#'
#' @param y Numeric vector. Vegetation index time-series.
#' @param t Integer vector or object of class \code{\linkS4class{Date}}. Time of y.
#' @param minValue_ylu Double. Input minimum y value allowed as baseline in cut cycle (default: 0).
#' @param s_lag Integer. Number of lag days for processing single season (default: 5).
#' @param c_win Integer. Window in days corresponding to the minimum duration of increasing and decreasing around peak (default: 20).
#' @param maxExtendMonth Integer. Maximum number of months for extending time series when performing cut cycle (default: 4). Passed to function \code{\link{phenofit::season_mov}}
#' @param ypeak_min Double. Minimum time series value to be considered as peak in cut cycle (local maximum) (default: 0.1). Passed to function \code{\link{phenofit::season_mov}}
#' @param minpeakdistance Integer. Minimum number of days between two consecutive peaks in cut cycle (local maximum) (default: 91). Passed to function \code{\link{phenofit::season_mov}}
#' @param length_min Integer. Minimum number of days for a single season identified in cut cycle (default: 45). Passed to function \code{\link{phenofit::season_mov}}
#' @param r_min Double. Used to eliminate fake peaks and troughs (default: 0.05; in (0, 1)). Passed to function \code{\link{phenofit::season_mov}}
#' @param r_max Double. Used to eliminate fake peaks and troughs (default: 0.2; in (0, 1)). Passed to function \code{\link{phenofit::season_mov}}
#' @param rtrough_max Double. Used to eliminate fake peaks and troughs (default: 0.5; in (0, 1)). Passed to function \code{\link{phenofit::season_mov}}
#' @param nptperyear Integer. Number of observations per year in time series (default: 365).
#' @param max_season Integer. Maximum number of seasons to be processed after cut cycle (default: 3).
#' @param der_method Character. Method to compute derivative. Can be one of: 'none' (default), 'spline' ,'whittaker', 'movavg'.
#' @param force Logical. Force computation of phenological metrics using Asymmetric Gaussian fit in case of failure of Gu extraction method. Default: TRUE.
#' @param reorder Logical. Reorder fields of output data.frame to improve readability.
# #' @param ... Additional arguments to be passed through to function \code{\link{phenoEO4PM}}
#'
#' @return
#' The function returns phenological metrics estimated from daily time series. 
#'
#' @author Federico Filipponi
#'
#' @keywords EO4PM phenology
#' 
# #' @seealso
#'
#' @examples
#' \dontrun{
#' # estimate phenological metrics
#' y <- dnorm(seq(-4, 4, length.out = 500), mean = 0, sd = 1)*4
#' phenometrics <- phenoEO4PM(y = y, t = 1:length(y))
#' 
#' # plot result
#' plotEO4PM(ts=data.frame("y"=y, "t"=1:length(y)), phenometrics)
#' }
#'
#' @import phenofit
#' @importFrom stats smooth.spline
#' @importFrom stats predict
#' @importFrom imputeTS na_interpolation
#' @importFrom phenofit findpeaks
#' @importFrom phenofit season_mov
#' @importFrom phenofit smooth_wWHIT
#' @importFrom phenofit wTSM
#' @importFrom phenofit FitDL.AG
# #' @export

# ###########################################################
# Define function to extract time_set
# ###########################################################

time_set_extract <- function(obs_dates, full=TRUE){
  
  # initialize list
  time_set <- list()
  
  # extract year
  first_jan <- which(as.integer(strftime(as.Date(obs_dates[1:(length(obs_dates)*0.75)], origin="1970-01-01"), format = "%j")) == 1 )[1]
  # force attribution of first january when is outside the observed dates
  if(is.na(first_jan)){
    syear <- median(as.integer(strftime(as.Date(obs_dates, origin="1970-01-01"), format = "%Y")))
    first_jan_date <- as.integer(as.Date(paste(syear, "-01-01", sep=""), origin="1970-01-01"))
  } else {
    first_jan_date <- obs_dates[first_jan]
  }
  year <- as.integer(strftime(as.Date(first_jan_date, origin="1970-01-01"), format = "%Y"))
  
  time_set$first_jan <- first_jan
  time_set$first_jan_date <- first_jan_date
  time_set$year <- year
  time_set$maxline_end <- min(c(which(strftime(as.Date(obs_dates, origin="1970-01-01"), format = "%Y") == (year + 1) ), length(obs_dates)), na.rm=TRUE)
  
  if(full){
    # set temporal position to compute multitemporal metrics
    time_set$solyear <- which(as.integer(strftime(as.Date(obs_dates, origin="1970-01-01"), format = "%Y")) == year)
    # set temporal position of winter (DJF)
    time_set$soldjf1 <- which(as.integer(strftime(as.Date(obs_dates, origin="1970-01-01"), format = "%Y")) == (year-1) & as.integer(strftime(as.Date(obs_dates, origin="1970-01-01"), format = "%m")) == 12)
    time_set$soldjf2 <- which(as.integer(strftime(as.Date(obs_dates, origin="1970-01-01"), format = "%Y")) == year & as.integer(strftime(as.Date(obs_dates, origin="1970-01-01"), format = "%m")) %in% c(1,2))
    # soldjf <- c(soldjf1,soldjf2)
    
    # set temporal position of summer (JJA)
    time_set$solmam <- which(as.integer(strftime(as.Date(obs_dates, origin="1970-01-01"), format = "%m")) %in% c(3,4,5))
    time_set$soljja <- which(as.integer(strftime(as.Date(obs_dates, origin="1970-01-01"), format = "%m")) %in% c(6,7,8))
    time_set$solson <- which(as.integer(strftime(as.Date(obs_dates, origin="1970-01-01"), format = "%m")) %in% c(9,10,11) & as.integer(strftime(as.Date(obs_dates, origin="1970-01-01"), format = "%Y")) %in% time_set$year)
  }
  
  return(time_set)
}

# ###########################################################
# Define smoothing function compatible with phenofit to skip smoothing in cut_cycle
# ###########################################################

smooth_none <- function(y, w, ylu, nptperyear, ...) {
  y <- as.numeric(y)
  w <- as.numeric(w)
  list(
    zs = list(ziter1 = y),
    ws = list(witer1 = w)
  )
}

# ###########################################################
# Define function to compute Phenological Metrics
# ###########################################################

phenoGetEO4PM <- function(pixel, obs_dates, year=NULL, c_win=20, s_lag=0, der_method="none", ...){
  
  if(!(der_method %in% c("none","spline","whittaker","movavg"))){
    stop("Argument 'der_method' must be set to one of the followings: 'none','spline','whittaker','movavg'")
  }
  
  phenokey_result <- as.vector(rep(as.double(NA),37))
  names(phenokey_result) <- c("PM_mask", "SoS_date", "SoS_doy", "SoS_value", "SGS_date", "SGS_doy", "SGS_value", "greenup_date", "greenup_doy", "greenup_value", "greenup_rate",
                              "SMP_date", "SMP_doy", "SMP_value", "PoS_date", "PoS_doy", "PoS_value", "EGS_date", "EGS_doy", "EGS_value",
                              "senescence_date", "senescence_doy", "senescence_value", "senescence_rate", "EoS_date", "EoS_doy", "EoS_value",
                              "seasonal_amplitude", "DOS", "LMP", "maturity_plateau_slope","STI","NSTI", "SMV","PoS_AG_diff","STI_AG_diff","QF")
  
  if(anyNA(pixel) | length(pixel)<(c_win*4+1) | c_win < 1 | anyNA(obs_dates)){
    # if(verbose){
    #   message("Too few temporal observation for phenology key stages extraction. Try using 365 or more observations")
    # }
    phenokey_result[1] <- 0
    phenokey_result[37] <- 3
  } else {
    
    # Initialize Quality Flags
    # qf <- 1
    qf0 <- 1
    qf1 <- 0
    qf2 <- 0
    qf3 <- 0
    qf4 <- 0
    qf5 <- 0
    qf6 <- 0
    qf7 <- 0
    qf8 <- 0
    qf9 <- 0
    qf10 <- 0
    
    #############################################################
    # Extract Phenological Metrics
    #############################################################
    
    ### add checks
    check1 <- 0
    check2 <- 0
    check3 <- 0
    # check if 'pixel' has the same length of 'obs_dates'
    if(length(pixel) != length(obs_dates)){
      check1 <- 1
    }
    # check if obs_dates are ordered and daily
    if(is.unsorted(obs_dates)){
      check2 <- 1
    }
    # check if 'obs_dates' are equally spaced by 1 day
    if(length(obs_dates) != length(obs_dates[1]:obs_dates[length(obs_dates)])){
      check3 <- 1
    }
    
    if(sum(check1,check2,check3)>0){
      metrics <- rep(NA, 10)
      qf1 <- 1
    } else {
      if(!is.null(pixel)){
        
        if(der_method == "none"){
          # der1 <- c(diff(pixel), 0)
          # usa la differenza centrale
          der1 <- c(0, (pixel[3:length(pixel)] - pixel[1:(length(pixel)-2)]) / 2, 0)
          # use Whittaker as fallback
          if(length(which.min(der1)) > 1 | length(which.max(der1)) > 1){
            fit_der1 <- phenofit::smooth_wWHIT(y=pixel, w=rep(1, length(pixel)), ylu=c(min(pixel, na.rm=TRUE), max(pixel, na.rm=TRUE)), nptperyear = 365, iters=1, lambda=1500, second=FALSE)
            der1 <- c(diff(as.double(fit_der1[[1]][[length(fit_der1[[1]])]])), 0)
          }
        }
        
        # der1 using spline
        if(der_method == "spline"){
          spline.eq <- tryCatch({
            stats::smooth.spline(pixel, df = floor(length(pixel) * 0.1))
          }, error = function(e){
            stats::smooth.spline(jitter(pixel, factor=0.1), df = floor(length(pixel) * 0.1))
          })
          der1 <- stats::predict(spline.eq, d = 1)$y
        }
        
        # der1 using Whittaker
        if(der_method == "whittaker"){
          fit_der1 <- phenofit::smooth_wWHIT(y=pixel, w=rep(1, length(pixel)), ylu=c(min(pixel, na.rm=TRUE), max(pixel, na.rm=TRUE)), nptperyear = 365, iters=1, lambda=1500, second=FALSE)
          # fit_der1 <- phenofit::whit2(y=pixel, lambda=1500)
          der1 <- c(diff(as.double(fit_der1[[1]][[length(fit_der1[[1]])]])), 0)
        }
        
        # # der1 using 'pracma::movavg()'
        if(der_method == "movavg"){
          if (!("pracma" %in% (.packages()))) {
            suppressPackageStartupMessages(library("pracma", character.only = TRUE, warn.conflicts = FALSE))
          }
          fit_der1 <- pracma::movavg(pixel, n = 7, type = "t")
          # der1 <- c(0, diff(fit_der1))
          # usa la differenza centrale
          der1 <- c(0, (fit_der1[3:length(fit_der1)] - fit_der1[1:(length(fit_der1)-2)]) / 2, 0)
        }
        
        days <- as.vector(as.integer(obs_dates))
        expectedNA <- 0
      }
      
      metrics <- rep(NA, 10)
      if(class(der1) == "try-error" | length(which(is.na(der1) == TRUE)) != expectedNA | length(which(is.infinite(der1) == TRUE)) != 0){
        qf2 <- 1
      } else {
        
        # get time_set information
        if(is.null(year)){
          time_set <- time_set_extract(obs_dates=obs_dates, full=FALSE)
          first_jan_date <- time_set$first_jan_date
          year <- time_set$year
          maxline_end <- time_set$maxline_end
        } else {
          first_jan_date <- as.integer(as.Date(paste(year, "-01-01", sep=""), origin="1970-01-01"))
          maxline_end <- min(c(which(strftime(as.Date(obs_dates, origin="1970-01-01"), format = "%Y") == (year + 1) ), length(obs_dates)), na.rm=TRUE)
        }
        
        # ##############
        # identify temporal windows where to look for minimum and maximum 1st derivative values
        
        # define frame
        fr <- floor(c_win/(2*0.80))
        positive_frame <- c_win
        negative_frame <- c_win
        
        # define time period to look for max PoS day
        max_pos_min_date <- which(obs_dates >= as.integer(as.Date(paste(year, "-01-01", sep=""), origin="1970-01-01")))[1]
        max_pos_max_date <- rev(which(obs_dates <= as.integer(as.Date(paste(year, "-12-31", sep=""), origin="1970-01-01"))))[1]
        # # in caso non si utilizzi il cut_cycle
        # max_pos_min_date <- which(obs_dates >= as.integer(as.Date(paste(year, "-03-01", sep=""), origin="1970-01-01")))[1]
        # max_pos_max_date <- rev(which(obs_dates <= as.integer(as.Date(paste(year, "-09-30", sep=""), origin="1970-01-01"))))[1]
        
        maxline_end <- min(c(maxline_end, max_pos_max_date), na.rm=TRUE)
        
        # identify max position
        max_pos <- median(which.max(pixel[max_pos_min_date:max_pos_max_date])+(max_pos_min_date - 1))
        max_pos_day <- obs_dates[max_pos]
        
        # find relative minima before and after max
        # minima_pos <- fast_findpeaks_cpp(-pixel, minpeakdistance = 45)$X$pos
        minima_pos <- phenofit::findpeaks(-pixel, minpeakdistance = 45)$X$pos
        if(length(minima_pos) > 0){
          minima_pre_pos <- rev(minima_pos[which(minima_pos < max_pos)])[1]
          minima_post_pos <- minima_pos[which(minima_pos > max_pos)][1]
        }
        
        if(max_pos_day < first_jan_date){
          qf3 <- 1
        } else {
          # compute potential positions of positive der1
          der1pos <- rep(0, length(der1))
          der1_pos <- rep(0, length(der1))
          for(i in 1:length(der1pos)){
            der1pos[i] <- as.integer(sum(der1[max(1, i-fr):min(length(der1), i+fr)] > 0) > positive_frame)
          }
          for(i in 1:length(der1_pos)){
            der1_pos[i] <- as.integer(max(der1pos[max(1, i-fr):min(length(der1pos), i+fr)]))
          }
          # remove positions after max
          der1_pos[max_pos:length(der1_pos)] <- 0
          # remove positions before relative minima before max
          if(exists("minima_pre_pos")){
            if(length(minima_pre_pos) == 1){
              if(!is.na(minima_pre_pos)){
                der1_pos[1:minima_pre_pos] <- 0
              }
            }
          }
          
          # compute potential positions of negative der1
          der1neg <- rep(0, length(der1))
          for(i in 1:length(der1neg)){
            der1neg[i] <- as.integer(sum(der1[max(1, i-fr):min(length(der1), i+fr)] < 0) > negative_frame)
          }
          der1_neg <- rep(0, length(der1))
          for(i in 1:length(der1_neg)){
            der1_neg[i] <- as.integer(max(der1neg[max(1, i-fr):min(length(der1neg), i+fr)]))
          }
          # remove positions before max
          der1_neg[1:max_pos] <- 0
          # remove positions after relative minima after max
          if(exists("minima_post_pos")){
            if(length(minima_post_pos) == 1){
              if(!is.na(minima_post_pos)){
                der1_neg[minima_post_pos:length(der1_neg)] <- 0
              }
            }
          }
          
          # proceed only if there are identified time windows
          if(max(der1_pos)<1){qf4 <- 1}
          if(max(der1_neg)<1){qf5 <- 1}
          if(max(der1_pos)<1 || max(der1_neg)<1){
            metrics <- rep(NA, 10)
          } else {
            
            # ### optimized version
            # Estrazione PRR e PSR
            idx_pos <- which(der1_pos == 1)
            prr <- max(der1[idx_pos], na.rm = TRUE)
            prd_idx <- idx_pos[which.max(der1[idx_pos])]
            prd <- obs_dates[prd_idx]
            
            # Aggiorna maschera negativa per forzare ordine cronologico
            der1_neg[1:prd_idx] <- FALSE
            idx_neg <- which(der1_neg == 1)
            
            if(length(idx_neg) == 0) {
              # qf_vec[7] <- TRUE
              metrics <- rep(NA, 10)
            } else {
              
              psr <- min(der1[idx_neg], na.rm = TRUE)
              psd_idx <- idx_neg[which.min(der1[idx_neg])]
              psd <- obs_dates[psd_idx]
              
              # 5. Calcoli Geometrici (Lineare Regression)
              # Invece di lm() completo (lento), calcoliamo solo i coefficienti necessari
              # Formula retta: y = mx + q -> q = y - mx
              q_rl <- pixel[prd_idx] - prr * prd
              q_sl <- pixel[psd_idx] - psr * psd
              
              Sbaseline <- min(pixel[1:max_pos], na.rm = TRUE)
              Ebaseline <- min(pixel[max_pos:length(pixel)], na.rm = TRUE)
              maxline <- max(pixel[prd_idx:psd_idx], na.rm = TRUE)
              
              UD <- (Sbaseline - q_rl) / prr
              RD <- (Ebaseline - q_sl) / psr
              SD <- (maxline - q_rl) / prr
              DD <- (maxline - q_sl) / psr
              old.DD <- DD
              
              # --- Calcolo del Plateau con Fallback Robusto ---
              sub_idx <- which(obs_dates >= SD & obs_dates <= DD)
              plateau.slope <- NA
              
              if (length(sub_idx) > 3) {
                # Uso di .lm.fit per la massima velocità
                x_matrix <- matrix(c(rep(1, length(sub_idx)), as.numeric(obs_dates[sub_idx])), ncol = 2)
                # plateau.lm <- try(.lm.fit(x_matrix, as.numeric(pixel[sub_idx])), silent = TRUE)
                plateau.lm <- .lm.fit(x_matrix, as.numeric(pixel[sub_idx]))
                
                # estrazione coefficienti
                p_coef        <- plateau.lm$coefficients
                q_plat        <- p_coef[1]
                plateau.slope <- p_coef[2]
                
                # Gestione dei casi in cui il calcolo produce NA numerici (es. collinearità perfetta)
                if (!is.na(plateau.slope)) {
                  denom <- plateau.slope - psr
                  if (abs(denom) > 1e-6 && plateau.slope <= 0) {
                    DD_tmp <- (q_sl - q_plat) / denom
                    if (!is.na(DD_tmp) && DD_tmp > SD && DD_tmp <= max(obs_dates)) {
                      DD <- round(DD_tmp)
                    } else {
                      DD <- round(psd_idx + (obs_dates[1] - 1))
                    }
                  } else {
                    DD <- round(old.DD)
                  }
                } else {
                  DD <- round(max_pos_day)
                }
              } else {
                DD <- round(max_pos_day)
              }
              
              # prr <- suppressWarnings(max(der1[which(der1_pos == 1)], na.rm = T))
              # prd <- days[which(der1_pos == 1)][which.max(der1[which(der1_pos == 1)])]
              #
              # # update 'der1_neg'
              # der1_neg[1:(which(days == prd))] <- 0
              #
              # # minimum value in 1st order derivative should occur after detected maximum value
              # psr <- suppressWarnings(min(der1[which(der1_neg == 1)], na.rm = T))
              # psd <- days[which(der1_neg == 1)][which.min(der1[which(der1_neg == 1)])]
              #
              # # stop process if both 'prd' and 'psd' are not identified
              # if(!(!is.null(prr) & !is.null(prd) & !is.null(psr) & !is.null(psd) & length(psd)>0 & length(prd)>0)){
              #   qf6 <- 1
              # } else {
              #   y.prd <- pixel[which(days == prd)]
              #   rl.y <- prr * days + y.prd - prr * prd
              #   rl.eq <- lm(rl.y ~ days)
              #   y.psd <- pixel[which(days == psd)]
              #   sl.y <- psr * days + y.psd - psr * psd
              #   sl.eq <- lm(sl.y ~ days)
              #   Sbaseline <- min(pixel[1:max_pos], na.rm = T)
              #   Ebaseline <- min(pixel[max_pos:length(pixel)], na.rm = T)
              #   # maxline <- max(pixel[max_pos_min_date:maxline_end], na.rm = T)
              #   maxline <- max(pixel[which(days == psd):which(days == prd)], na.rm = T)
              #   UD <- (Sbaseline - rl.eq$coefficients[1])/rl.eq$coefficients[2]
              #   RD <- (Ebaseline - sl.eq$coefficients[1])/sl.eq$coefficients[2]
              #   SD <- (maxline - rl.eq$coefficients[1])/rl.eq$coefficients[2]
              #   DD <- (maxline - sl.eq$coefficients[1])/sl.eq$coefficients[2]
              #   old.DD <- DD
              #
              #   sub.time <- days[which(days >= SD & days <= DD)]
              #   sub.gcc <- pixel[which(days >= SD & days <= DD)]
              #   if(length(sub.time) > 3){
              #     plateau.lm <- try(lm(sub.gcc ~ sub.time))
              #     if(class(plateau.lm) != "try-error"){
              #       M <- matrix(c(coef(plateau.lm)[2], coef(sl.eq)[2],
              #                     -1, -1), nrow = 2, ncol = 2)
              #       intercepts <- as.matrix(c(coef(plateau.lm)[1],
              #                                 coef(sl.eq)[1]))
              #       interception <- -solve(M) %*% intercepts
              #
              #       DD <- interception[1, 1]
              #     }
              #   }
              #   plateau.slope <- try(plateau.lm$coefficients[2], silent = TRUE)
              #   if(plateau.slope > 0){
              #     DD <- old.DD
              #   }
              #   plateau.intercept <- try(plateau.lm$coefficients[1],
              #                            silent = TRUE)
              #   if(class(plateau.slope) == "try-error"){
              #     plateau.slope <- NA
              #     plateau.intercept <- NA
              #   }
              #   # round values to closest integer
              UD <- round(as.double(UD))
              RD <- round(as.double(RD))
              SD <- round(as.double(SD))
              DD <- round(as.double(DD))
              
              metrics <- c(UD, SD, DD, RD, maxline, Sbaseline, Ebaseline, prr, psr, plateau.slope)
            }
          }
        }
      }
    }
    names(metrics) <- c("UD", "SD", "DD", "RD", "maxline", "Sbaseline", "Ebaseline", "prr", "psr", "plateau.slope")
    
    # perform check before continuing
    check_null <- as.logical(any(is.na(metrics[c(1,2,3,4,6,8,9,10)])))
    if(!any(is.na(metrics[c(1,2,3,4)]))){
      check_low <- as.logical(any(metrics[c(1,2,3,4)] < min(obs_dates)))
    } else {
      check_low <- FALSE
    }
    # check_null <- as.logical(any(is.na(metrics))) ### check all metrics instead of all relevant
    check_ordered <- as.logical(sum(sort(na.omit(metrics[1:4])) == na.omit(metrics[1:4]))<length(na.omit(metrics[1:4])))
    if(check_ordered){
      qf7 <- 1
    }
    
    check_within_year <- as.logical(sum(as.integer(na.omit(metrics[c(2,3)]) %in% (obs_dates[1]:obs_dates[length(obs_dates)])))<2)
    # check_within_year <- as.logical(sum(as.integer(na.omit(metrics[1:4]) %in% seq(obs_dates[1], obs_dates[length(obs_dates)], by=1)))<4)
    
    # update Quality Flags
    if(check_within_year | check_null | check_low){
      qf8 <- 1
    }
    
    # set Quality Flag to output dataset
    
    # Quality Flags meaning
    # qf0: Processed seasonal cycle
    # qf1: Pixel or date input problem
    # qf2: First derivative error
    # qf3: Identified peak date is before reference year (before 1st January)
    # qf4: No increasing period found
    # qf5: No decreasing period found
    # qf6: Gu stages identification failed
    # qf7: Estimated Gu stages are not chronologically ordered
    # qf8: Estimated Gu stages are not within input temporal range or have missing metrics values
    # qf9: Cannot detect relative minimum (SoS) using inflection point
    # qf10 Identified Gu stages are not within the analyzed time period (Identified SoS or SGS or EGS or EoS date falls outside season temporal range)
    # qf11: Failed to cut cycles
    # qf12: Identified more than specified max_season
    # qf13: Incomplete set of Phenological Metrics (e.g. due to out-of-range issues)
    # qf14: Used Asymmetric Gaussian fit after Gu phenometric estimate failure (only applied after successfull cut-cycle)
    
    # Extract Phenological Metrics
    # if(check_null | check_ordered | check_within_year | check_low){ # check_within_year AND check_low utilizzati solo come QF
    if(check_null | check_ordered){
      qf <- as.integer(qf0 * 2^0 + qf1 * 2^1 + qf2 * 2^2 + qf3 * 2^3 + qf4 * 2^4 + qf5 * 2^5 + qf6 * 2^6 + qf7 * 2^7 + qf8 * 2^8)
      phenokey_result[37] <- as.integer(qf)
      phenokey_result[1] <- 0
    } else {
      
      # set mask value to '1'
      phenokey_result[1] <- as.integer(1)
      
      # extract SGS
      sgs <- max(1,which(obs_dates == as.integer(metrics[1])), na.rm=TRUE)
      phenokey_result[5] <- as.integer(metrics[1])
      phenokey_result[6] <- as.integer(1 - first_jan_date + as.integer(metrics[1]))
      phenokey_result[7] <- pixel[sgs]
      # phenokey_result[7] <- ifelse(as.integer(metrics[1]) %in% obs_dates, pixel[sgs], NA)
      
      # extract SoS
      # sos <- max(1, sort(which(pixel[1:which(obs_dates == as.integer(metrics[1]))] == min(pixel[1:which(obs_dates == as.integer(metrics[1]))])), decreasing = TRUE)[1])
      # sos_pos <- as.integer(which(diff(sign(diff(pixel))) == 2) + 1)
      # minima_pos <- phenofit::findpeaks(-pixel, minpeakdistance = 45)$X$pos
      sos_pos <- minima_pos[which(minima_pos <= sgs)]
      sos <- sos_pos[which.min(abs(sos_pos - sgs))]
      # sos <- sos_pos[rev(which(sos_pos < sgs))[1]]
      if(length(sos) == 0){
        sos <- 1 + s_lag
        qf9 <- 1
      }
      phenokey_result[2] <- as.integer(obs_dates[sos])
      phenokey_result[3] <- as.integer(1 - first_jan_date + obs_dates[sos])
      phenokey_result[4] <- pixel[sos]
      
      # ### TODO use phenofit:::check_GS_HeadTail() to Check growing season head and tail minimum values https://rdrr.io/cran/phenofit/man/check_GS_HeadTail.html
      
      # extract SMP
      # smp <- which(obs_dates == as.integer(metrics[2]))
      phenokey_result[12] <- as.integer(metrics[2])
      phenokey_result[13] <- as.integer(1 - first_jan_date + as.integer(metrics[2]))
      phenokey_result[14] <- pixel[which(obs_dates == as.integer(metrics[2]))]
      # phenokey_result[14] <- ifelse(as.integer(metrics[2]) %in% obs_dates, pixel[which(obs_dates == as.integer(metrics[2]))], NA)
      
      # extract EGS
      if(as.integer(metrics[3]) <= obs_dates[length(obs_dates)]){
        egs <- which(obs_dates == as.integer(metrics[3]))
      } else {
        egs <- length(obs_dates) - s_lag
        qf10 <- 1
      }
      phenokey_result[18] <- as.integer(obs_dates[egs])
      phenokey_result[19] <- as.integer(1 - first_jan_date + obs_dates[egs])
      phenokey_result[20] <- pixel[egs]
      
      # extract EoS
      if(min(c(1000000,as.integer(metrics[4])), na.rm=TRUE) <= obs_dates[length(obs_dates)]){
        eos <- which(obs_dates == as.integer(metrics[4]))
      } else {
        eos <- length(obs_dates) - s_lag
        qf10 <- 1
      }
      phenokey_result[25] <- as.integer(obs_dates[eos])
      phenokey_result[26] <- as.integer(1 - first_jan_date + obs_dates[eos])
      phenokey_result[27] <- pixel[eos]
      
      # extract peak date
      # peak <- median(which.max(pixel[sgs:egs]))+sgs
      # peak_value <- pixel[max_pos]
      phenokey_result[15] <- as.integer(max_pos_day)
      phenokey_result[16] <- as.integer(1 - first_jan_date + obs_dates[max_pos])
      if(!(is.na(phenokey_result[16]))){
        if(phenokey_result[16] < 1){
          qf3 <- 1
        }
      }
      phenokey_result[17] <- pixel[max_pos]
      
      # extract amplitude
      amplitude <- as.double(pixel[max_pos]-pixel[sos])
      phenokey_result[28] <- amplitude
      
      # extract maximum growth rate (observation date)
      max_growth_rate <- which(obs_dates == prd)
      phenokey_result[8] <- as.integer(obs_dates[max_growth_rate])
      phenokey_result[9] <- as.integer(1 - first_jan_date + obs_dates[max_growth_rate])
      phenokey_result[10] <- as.double(pixel[max_growth_rate])
      phenokey_result[11] <- as.double(metrics[8])
      
      # extract maximum degrowth rate (observation date)
      max_degrowth_rate <- which(obs_dates == psd)
      phenokey_result[21] <- as.integer(obs_dates[max_degrowth_rate])
      phenokey_result[22] <- as.integer(1 - first_jan_date + obs_dates[max_degrowth_rate])
      phenokey_result[23] <- as.double(pixel[max_degrowth_rate])
      phenokey_result[24] <- as.double(metrics[9])
      
      # extract maturity pleateau slope
      phenokey_result[31] <- as.double(metrics[10])
      
      # extract DoS (Duration of Season) (in days)
      phenokey_result[29] <- as.integer(obs_dates[eos] - obs_dates[sos]) # use SoS
      # phenokey_result[29] <- as.integer(obs_dates[eos] - obs_dates[sgs]) # use SGS
      
      # extract LMP (Length of Maturity Plateau) (in days)
      phenokey_result[30] <- as.integer(metrics[3] - metrics[2])
      
      # extract Seasonal Time Integrated value
      sti <- as.double(sum(pixel[sos:eos]))
      phenokey_result[32] <- sti
      # extract Net Seasonal Time Integrated value
      nsti <- as.integer(sti - (pixel[sos] * phenokey_result[29]))
      phenokey_result[33] <- nsti
      
      # calculate Seasonal Mean Value (SMV)
      phenokey_result[34] <- sti / phenokey_result[29]
      
      # write default values for comparison with Asymmetric Gaussian fit
      phenokey_result[35] <- 0.0
      phenokey_result[36] <- 0.0
      
      # aggregate QF
      qf <- as.integer(qf0 * 2^0 + qf1 * 2^1 + qf2 * 2^2 + qf3 * 2^3 + qf4 * 2^4 + qf5 * 2^5 + qf6 * 2^6 + qf7 * 2^7 + qf8 * 2^8 + qf9 * 2^9 + qf10 * 2^10)
      phenokey_result[37] <- as.integer(qf)
      
    }
  }
  return(phenokey_result)
}

# define function to orchestrate Phenological Metrics estimation from single temporal profile
phenoEO4PM <- function(y, t, minValue_ylu=0, s_lag=5, c_win=20, maxExtendMonth=4, ypeak_min=1.0, minpeakdistance=91, length_min=45, r_min=0.05, r_max=0.2, rtrough_max=0.5, nptperyear=365, max_season=3, der_method="none", force=TRUE, reorder=TRUE){
  
  if (!("phenofit" %in% (.packages()))) {
    suppressPackageStartupMessages(library("phenofit", character.only = TRUE, warn.conflicts = FALSE))
  }
  
  if(anyNA(y)){
    stop("Argument 'y' contains missing values.")
  }
  if(anyNA(t)){
    stop("Argument 't' contains missing values.")
  }
  
  # initialize QFs
  qf0 <- 0
  qf11 <- 0
  qf12 <- 0
  # qf13 <- 0
  qf14 <- 0
  
  # ###########################################################
  # Cut cycles
  # ###########################################################
  
  # define weights to clean baseline-corrected curve artifacts
  w <- rep(as.double(1.0), length(t))
  # w[which(pixel[2:length(pixel)] == pixel[1:(length(pixel)-1)]) + 1] <- 0
  
  # convert to phenofit INPUT object
  ylu <- quantile(as.double(y), c(0.01,0.99), na.rm=TRUE)
  ylu[1] <- pmax(ylu[1], min(0, minValue_ylu, na.rm=TRUE), na.rm=TRUE)
  ylu[2] <- max(y, na.rm=TRUE)
  pfo <- list(t = as.Date(t),
              y = as.double(y),
              w = w,
              ylu = ylu,
              nptperyear = 365,
              south = as.logical(FALSE)
  )
  # detect cycles
  s <- suppressWarnings(phenofit::season_mov(pfo,
                                             options = list(
                                               rFUN = "smooth_wWHIT", 
                                               # rFUN = "smooth_none", 
                                               wFUN = "wTSM",
                                               iters = 1, lambda = 10,
                                               maxExtendMonth = maxExtendMonth,
                                               minpeakdistance = minpeakdistance,
                                               length_min = length_min,
                                               ypeak_min = ypeak_min,
                                               rtrough_max = rtrough_max,
                                               r_min = r_min,
                                               r_max = r_max,
                                               rm.closed=TRUE,
                                               adj.param = FALSE,
                                               verbose = FALSE
                                             )))
  if(is.null(s)){
    nseas <- 0
    qf11 <- 1
  } else {
    sdt <- s[["dt"]]
    nseas <- nrow(sdt)
    # if(nseas > max_season){
    #   sdt <- sdt[sort(order(sdt$y_peak, decreasing = TRUE)[1:max_season]),]
    #   # nseas <- max_season
    #   nseas <- nrow(sdt)
    #   qf12 <- 1
    # }
  }
  
  suppressWarnings(rm(list=c("ylu","pfo","s","w")))
  
  # ###########################################################
  # Estimate phenological metrics
  # ###########################################################
  
  if(nseas > 0){
    
    for(s in 1:nseas){
      
      # set QFs
      qf0 <- 1
      
      s_beg <- as.integer(sdt$beg[s])
      s_pos <- which(t %in% s_beg)
      s_end <- as.integer(sdt$end[s])
      e_pos <- which(t %in% s_end)
      s_year <- as.integer(strftime(sdt$peak[s], format = "%Y"))
      
      # compute phenometrics
      pks_result <- phenoGetEO4PM(pixel=as.double(y[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)]), obs_dates = t[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)], year=s_year, s_lag=s_lag, der_method=der_method, c_win=c_win)
      
      # fit Asymmetric Gaussian curve if Gu Method fails
      if(force){
        if(pks_result[1] == 0){
          
          # fit Asymmetric Gaussian curve
          # pars <- phenofit::FitDL.Beck(y = y[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)], t = t[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)], tout = t[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)], w=ifelse(y[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)] > quantile(y[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)], probs=0.75, na.rm=TRUE), 1.0, 0.25), type = 1, iters = 2)
          # dl <- phenofit::doubleLog.Beck(par = as.double(phenofit::get_param(pars)), t=t[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)])
          dl <- phenofit::FitDL.AG(y = y[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)], t = t[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)], tout = t[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)], w=ifelse(y[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)] > quantile(y[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)], probs=0.75, na.rm=TRUE), 1.0, 0.25), type = 1, iters = 2)
          dl <- dl$zs$iter2
          
          # compute phenometrics
          pks_result <- phenoGetEO4PM(pixel=dl, obs_dates = t[max(1, s_pos-s_lag):min(length(t), e_pos+s_lag)], year=s_year, s_lag=s_lag, der_method=der_method, c_win=c_win)
          
          pks_result[35] <- max(dl[((s_lag+1):(length(dl-s_lag)))], na.rm=TRUE) - max(y[(s_pos:e_pos)], na.rm=TRUE)
          pks_result[36] <- sum(dl[((s_lag+1):(length(dl-s_lag)))], na.rm=TRUE) - sum(y[(s_pos:e_pos)], na.rm=TRUE)
          qf14 <- 1
          suppressWarnings(rm(list=c("dl")))
        }
      }
      
      # update  Quality Flags
      pks_result[37] <- as.integer(pks_result[37] + qf0 * 2^0 + qf11 * 2^11 + qf12 * 2^12 + qf14 * 2^14)
      
      # write values for Start of Cycle (SoC) and End of Cycle (EoC)
      cyc_result <- rep(NA, 6)
      names(cyc_result) <- c("SoC_date", "SoC_doy", "SoC_value", "EoC_date", "EoC_doy", "EoC_value")
      cyc_result[1] <- s_beg
      cyc_result[2] <- as.integer(strftime(as.Date(s_beg, origin="1970-01-01"), format = "%j"))
      cyc_result[3] <- as.double(y[s_pos])
      cyc_result[4] <- s_end
      cyc_result[5] <- as.integer(strftime(as.Date(s_end, origin="1970-01-01"), format = "%j"))
      cyc_result[6] <- as.double(y[e_pos])
      
      # bind result to data.frame
      if(exists("pm_df")){
        pm_df <- rbind(pm_df, data.frame(t(c(cyc_result, pks_result))))
      } else {
        pm_df <- data.frame(t(c(cyc_result, pks_result)))
      }
    }
    
  } else {
    
    # compute phenometrics
    pks_result <- phenoGetEO4PM(pixel=y, obs_dates = t, der_method=der_method, c_win=c_win)
    
    # update  Quality Flags
    pks_result[37] <- as.integer(pks_result[37] + qf0 * 2^0 + qf11 * 2^11 + qf12 * 2^12 + qf14 * 2^14)
    
    # write values for Start of Cycle (SoC) and End of Cycle (EoC)
    cyc_result <- rep(NA, 6)
    names(cyc_result) <- c("SoC_date", "SoC_doy", "SoC_value", "EoC_date", "EoC_doy", "EoC_value")
    rMin_pos <- phenofit::findpeaks(-y, minpeakdistance = 45)$X$pos
    
    if(!is.na(pks_result[15])){
      first_jan_date <- as.integer(as.Date(paste(as.integer(strftime(as.Date(pks_result[15], origin="1970-01-01"), format = "%Y")), "-01-01", sep=""), origin="1970-01-01"))
    } else {
      first_jan_date <- as.integer(as.Date(paste(as.integer(strftime(as.Date(median(t), origin="1970-01-01"), format = "%Y")), "-01-01", sep=""), origin="1970-01-01"))
    }
    
    SoC_pos <- rMin_pos[which(rMin_pos <= which(t == pks_result[5]))]
    SoC_pos <- SoC_pos[which.min(abs(SoC_pos - which(t == pks_result[5])))][1]
    if(length(SoC_pos) == 1 & !is.na(SoC_pos)){
      cyc_result[1] <- t[SoC_pos]
      cyc_result[2] <- as.integer(1 - first_jan_date + t[SoC_pos])
      cyc_result[3] <- as.double(y[SoC_pos])
    }
    EoC_pos <- rMin_pos[which(rMin_pos >= which(t == pks_result[25]))]
    EoC_pos <- EoC_pos[which.min(abs(EoC_pos - which(t == pks_result[25])))][1]
    if(length(EoC_pos) == 1 & !is.na(EoC_pos)){
      cyc_result[4] <- t[EoC_pos]
      cyc_result[5] <- as.integer(1 - first_jan_date + t[EoC_pos])
      cyc_result[6] <- as.double(y[EoC_pos])
    }
    # convert to data frame
    pm_df <- data.frame(t(c(cyc_result, pks_result)))
  }
  
  # add year column
  year_vec <- as.integer(strftime(as.Date(pm_df$PoS_date, origin="1970-01-01"), format = "%Y"))
  if(any(is.na(year_vec))){
    for(i in which(is.na(year_vec))){
      if(!is.na(pm_df$SoC_date[i]) & !is.na(pm_df$EoC_date[i])){
        year_vec[i] <- as.integer(strftime(as.Date(t[which(t %in% pm_df$SoC_date[i]:pm_df$EoC_date[i])][which.max(y[which(t %in% pm_df$SoC_date[i]:pm_df$EoC_date[i])])], origin="1970-01-01"), format = "%Y"))
      }
    }
  }
  pm_df <- cbind(data.frame("year"=year_vec, "cycle"=rep(as.integer(NA), nrow(pm_df)), "season"=rep(as.integer(NA), nrow(pm_df)), pm_df))
  for(i in sort(unique(pm_df$year))){
    if(length(which(pm_df$year == i)) == 1){
      pm_df$cycle[which(pm_df$year == i)] <- 1
      pm_df$season[which(pm_df$year == i)] <- 1
    } else {
      pm_df$season[which(pm_df$year == i)] <- order(pm_df$SoC_date[which(pm_df$year == i)])
      pm_df$cycle[which(pm_df$year == i)] <- order(pm_df$STI[which(pm_df$year == i)])
      # pm_df$cycle[which(pm_df$year == i)] <- order(pm_df$SMV[which(pm_df$year == i)])
    }
  }
  
  # remove cycles for years exceeding max_season
  max_season_rm <- which(pm_df$cycle > max_season)
  if(length(max_season_rm) > 0){
    qf12 <- 1
    pm_df$QF[which(pm_df$year %in% unique(pm_df$year[max_season_rm]))] <- as.integer(pm_df$QF[which(pm_df$year %in% unique(pm_df$year[max_season_rm]))] + qf12 * 2^12)
    pm_df <- pm_df[-max_season_rm,]
  }
  
  # reorder fields
  if(reorder){
    pm_df <- pm_df[,c("year", "season", "cycle", "PM_mask", 
                      "SoC_doy", "SoS_doy", "SGS_doy", "greenup_doy", "SMP_doy", "PoS_doy", "EGS_doy", "senescence_doy", "EoS_doy", "EoC_doy", 
                      "SoC_date", "SoS_date", "SGS_date", "greenup_date", "SMP_date", "PoS_date", "EGS_date", "senescence_date", "EoS_date", "EoC_date", 
                      "SoC_value", "SoS_value", "SGS_value", "greenup_value", "SMP_value", "PoS_value", "EGS_value", "senescence_value", "EoS_value", "EoC_value", 
                      "greenup_rate", "senescence_rate", "seasonal_amplitude", "DOS", "LMP", "maturity_plateau_slope","STI","NSTI", "SMV","PoS_AG_diff","STI_AG_diff","QF")]
  }
  
}
