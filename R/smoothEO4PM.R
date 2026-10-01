# Title         : smoothEO4PM.R
# Description   : Perform weighted smoothing and daily interpolation of time series
# Date          : May 2026
# Version       : 2.1
# Copyright name: CC BY
# Licence       : GPL v3
# Authors       : Federico Filipponi
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
# ##############################################################################
#' @title smoothEO4PM
#'
#' @description Perform smoothing and daily interpolation of time series
#' using weighted iterative Whittaker smoother and optional 
#' weighted iterative Savitski-Golay polynomial fitting (evenly spaced in time)
#'
#' @param y Numeric vector. Vegetation index time-series.
#' @param t Integer vector or object of class \code{\linkS4class{Date}}. Time of y.
#' @param w (optional) Numeric vector. Weights of y. If not specified, weights of all NA values will be 0.0, the others will be 1.0.
#' @param tout (optional) Integer vector or object of class \code{\linkS4class{Date}}. Target time dimension of smoothed interpolated y (daily timesteps required).
#' @param win_noise Integer odd. Temporal window (number of not NA y observations) to be used in noise removal step (small drops and optionally spikes). Set to 0 to skip noise removal.
#' @param win_savgol Integer odd. Lag window to be used for weighted iterative Savitski-Golay polynomial fitting (set to 0 to skip noise removal).
#' @param lambda Integer. Lambda to be used in weighted iterative Whittaker smoother (default: 300).
#' @param minValue Double. Input minimum threshold values for output y values.
#' @param iters Integer. Number of iterations of both Whittaker smoother and optional Savitski-Golay polynomial fitting (default: 3).
#' @param spike_removal Character. Set removal of small spikes found in the noise removal step. Can be one of: 'none' (default), 'sync', 'later'.
#' @param wFUN_type Character. Function to update weights. Can be one of: 'wBisquare' (default), 'wTSM', 'wChen', 'wSELF'.
#' @param fitMax Logical. Set additional Whittaker smoother pass to improve fit of local maxima.
#' @param raw Logical. Store input raw values in output object. Default: TRUE.
#' @param plot Logical. Display plot.
#' @param plot_filename (optional) Character. Output folder where to store plots.
# #' @param ... Additional arguments to be passed through to function \code{\link{smoothEO4PM}}
#'
#' @return
#' The function returns smoothed daily interpolated time series. 
#' Optionally plot can be generate and/or saved to output folder set using plot_filename argument.
#'
#' @author Federico Filipponi
#'
#' @keywords EO4PM smoothing
#'
# #' @seealso
#' 
#' @examples
#' \dontrun{
#' # temporal smoothing
#' y <- jitter(dnorm(seq(-4, 4, length.out = 500), mean = 0, sd = 1)*4, amount = 0.1)
#' tsg <- smoothEO4PM(y = y, t = 1:length(y))
#' 
#' # perform temporal smoothing and plot result
#' tsg <- smoothEO4PM(y = y, t = 1:length(y), plot = TRUE)
#' 
#' # temporal smoothing with additional arguments
#' tsg <- smoothEO4PM(y = y, t = 1:length(y), tout = 1:(length(y)-10), 
#'                   spikes_removal = "none", win_savgol = 35, lambda = 100, 
#'                   minValue = 0.05, iters = 4, wFUN_type = "wTSM", fitMax = TRUE)
#' }
#' 
#' @import imputeTS
#' @import phenofit
#' @import ggplot2
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @importFrom imputeTS na_interpolation
#' @importFrom phenofit sources
#' @importFrom phenofit smooth_wSG
#' @importFrom phenofit smooth_wWHIT
#' @importFrom phenofit wTSM
#' @importFrom phenofit wBisquare
#' @importFrom phenofit wChen
#' @importFrom phenofit wSELF
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ggtitle
# #' @export

# ###########################################################
# Define noise removal function to clean small drops and spikes
# ###########################################################

noiseRemoval <- function(y, w=3, spikes="none", borders=FALSE, max_pos=NULL){
  
  # load required libraries
  if (!("phenofit" %in% (.packages()))) {
    suppressPackageStartupMessages(library("phenofit", character.only = TRUE, warn.conflicts = FALSE))
  }
  
  if(!exists("cpp_med")){
    # define C++ function to compute median (3x faster)
    # Source: https://stackoverflow.com/questions/34771088/why-is-standard-r-median-function-so-much-slower-than-a-simple-c-alternative/34771504#34771504
    Rcpp::sourceCpp(code='
      #include <Rcpp.h>

      // [[Rcpp::export]]
      double cpp_med(Rcpp::NumericVector x){
      std::size_t size = x.size();
      std::sort(x.begin(), x.end());
      if (size  % 2 == 0) return (x[size / 2 - 1] + x[size / 2]) / 2.0;
      return x[size / 2];
      }', cleanupCacheDir=TRUE
    )
    
    # define C++ function (not clear cache in Rtmp folder)
    # Rcpp::cppFunction(
    #   'double cpp_med(Rcpp::NumericVector x){
    # std::size_t size = x.size();
    # std::sort(x.begin(), x.end());
    # if (size  % 2 == 0) return (x[size / 2 - 1] + x[size / 2]) / 2.0;
    # return x[size / 2];
    # }'
    # )
  }
  
  # fallback: define default median function
  if(!exists("cpp_med")){
    cpp_med <- median
  }
  
  if(!(spikes %in% c("none", "sync", "later"))){
    stop("Argument 'spikes' must be set to one of the followings: 'none' (default), 'sync', 'later'.")
  }
  
  # set window size for noise detection
  w <- floor((w - 1) / 2)
  
  pixel_pos <- which(!is.na(y))
  pixel_out <- y
  pixel <- as.vector(y[pixel_pos])
  
  if(!is.null(max_pos) & spikes != "none"){
    max_ppos <- which(pixel_pos %in% max_pos)
  }
  
  # ###################################
  # search the small drops and spikes
  spike_pos <- integer()
  for(i in (1+w):(length(pixel)-w)){
    app <- pixel[((i-w):(i+w))[-(w+1)]]
    med <- cpp_med(app)
    sta <- sd(app)
    
    # identify spikes
    if(abs(pixel[i]-med)>(2*sta)){
      spike_pos <- c(spike_pos, i)
    }
  }
  
  if(!is.null(max_pos) & spikes != "none"){
    spike_pos <- spike_pos[which(!(spike_pos %in% max_ppos))]
  }
  
  # remove small drops (spikes 'sync')
  if(length(spike_pos) > 0){
    for(i in spike_pos){
      # remove small drops
      if(pixel[i] < pixel[i-1] & pixel[i] < pixel[i+1]){
        pixel_out[pixel_pos[i]] <- NA
      }
      # remove small spikes
      if(spikes == "sync"){
        if(pixel[i] > pixel[i-1] & pixel[i] > pixel[i+1]){
          pixel_out[pixel_pos[i]] <- NA
        }
      }
    }
  }
  
  if(borders){
    # check first observations
    for(j in 1:w){
      pos <- (1:max(3, j + (j-1)))[-j]
      if(abs(pixel[j]-cpp_med(pixel[pos]))>(2*sd(pixel[pos]))){
        # identify small drops
        if(pixel[j] < pixel[j+1]){
          pixel_out[pixel_pos[j]] <- NA
        }
        # identify small spikes
        if(spikes %in% c("sync", "later")){
          if(pixel[j] > pixel[j+1]){
            pixel_out[pixel_pos[j]] <- NA
          }
        }
      }
    }
    # check last observations
    for(j in (length(pixel)-w):length(pixel)){
      pos <- (min(length(pixel)-2, j - (length(pixel)-j)):length(pixel))[-(abs(-length(pixel)+j)+1)]
      if(abs(pixel[j]-cpp_med(pixel[pos]))>(2*sd(pixel[pos]))){
        # identify small drops
        if(pixel[j] < pixel[j-1]){
          pixel_out[pixel_pos[j]] <- NA
        }
        # identify small spikes
        if(spikes %in% c("sync", "later")){
          if(pixel[j] > pixel[j-1]){
            pixel_out[pixel_pos[j]] <- NA
          }
        }
      }
    }
  }
  
  # remove small spikes (spikes 'later')
  if(spikes == "later"){
    # update valid positions
    pixel_pos <- which(!is.na(pixel_out))
    pixel <- as.vector(pixel_out[pixel_pos])
    # obs_dates <- as.vector(x[pixel_pos])
    if(!is.null(max_pos) & spikes != "none"){
      max_ppos <- which(pixel_pos %in% max_pos)
    }
    # search the small drops and spikes
    spike_pos <- integer()
    for(i in (1+w+1):(length(pixel)-(w+1))){
      app <- pixel[((i-w):(i+w))[-(w+1)]]
      med <- cpp_med(app)
      sta <- sd(app)
      
      # remove small spikes
      if(!is.null(max_pos) & length(max_pos) > 0){
        if(i %in% max_ppos){
          if(abs(pixel[i]-med)>(3*sta)){
            if(pixel[i] > pixel[i-1] & pixel[i] > pixel[i+1]){
              pixel_out[pixel_pos[i]] <- NA
            }
          }
        } else {
          if(abs(pixel[i]-med)>(2*sta)){
            if(pixel[i] > pixel[i-1] & pixel[i] > pixel[i+1]){
              pixel_out[pixel_pos[i]] <- NA
            }
          }
        }
      } else {
        if(abs(pixel[i]-med)>(2*sta)){
          if(pixel[i] > pixel[i-1] & pixel[i] > pixel[i+1]){
            pixel_out[pixel_pos[i]] <- NA
          }
        }
      }
    }
  }
  
  # return function result
  return(pixel_out)
}

# ###########################################################
# Define findpeaks function to identify local maxima
# ###########################################################

# findpeaks <- function (x, nups = 1, ndowns = nups, zero = "0", peakpat = NULL,
#                        minpeakheight = -Inf, minpeakdistance = 1,
#                        h_min = 0, h_max = 0,
#                        npeaks = 0, sortstr = FALSE,
#                        include_gregexpr = FALSE,
#                        IsPlot = F)
# {
#   stopifnot(is.vector(x, mode = "numeric") ||
#               is.vector(x, mode = "logical") || length(is.na(x)) == 0)
#   
#   # if (minpeakdistance < 1)
#   #     warning("Handling 'minpeakdistance < 1' is logically not possible.")
#   if (!zero %in% c("0", "+", "-"))
#     stop("Argument 'zero' can only be '0', '+', or '-'.")
#   
#   str_replace_midzero <- function(x) {
#     replace = function(x, replacement = "+") {
#       # paste(rep(replacement, nchar(x)), collapse = "")
#       strrep(replacement, nchar(x))
#     }
#     stringr::str_replace_all(x, "\\++0\\++", . %>% replace("+")) %>%
#       stringr::str_replace_all("-+0-+", . %>% replace("-"))
#   }
#   
#   Rcpp::cppFunction(
#     'double cpp_med(Rcpp::NumericVector x){
#     std::size_t size = x.size();
#     std::sort(x.begin(), x.end());
#     if (size  % 2 == 0) return (x[size / 2 - 1] + x[size / 2]) / 2.0;
#     return x[size / 2];
#     }'
#   )
#   
#   # # extend the use of findpeaks:
#   # if (IsDiff) {
#   #     xc <- sign(diff(x))
#   # } else {
#   #     xc <- x
#   # }
#   xc <- sign(diff(x))
#   xc <- paste(as.character(sign(xc)), collapse = "") %>%
#     {gsub("1", "+", gsub("-1", "-", .))}
#   xc <- str_replace_midzero(xc) # in v0.3.4
#   
#   if (zero != "0") xc <- gsub("0", zero, xc)
#   # if (is.null(peakpat)) peakpat <- sprintf("[+]{%d,}[-]{%d,}", nups, ndowns)
#   if (is.null(peakpat)) peakpat <- sprintf("[+]{%d,}[0]{0,}[-]{%d,}", nups, ndowns)
#   
#   # Fatal bug found at 20211114
#   # Because `diff` operation lead to the length of `x` reduced one
#   # Hence x2 should be `rc + attr(rc, "match.length")`, without `-1`.
#   rc <- gregexpr(peakpat, xc)[[1]]
#   
#   if (rc[1] < 0) return(NULL)
#   x1 <- rc
#   x2 <- rc + attr(rc, "match.length") # - 1
#   
#   attributes(x1) <- NULL
#   attributes(x2) <- NULL
#   n <- length(x1)
#   xv <- xp <- numeric(n)
#   for (i in 1:n) {
#     # for duplicated extreme values, get the median
#     vals <- x[x1[i]:x2[i]]
#     maxI <- which(vals == max(vals, na.rm = TRUE))
#     xp[i] <- floor(cpp_med(maxI)) + x1[i] - 1
#     xv[i] <- x[xp[i]]
#   }
#   inds <- which(xv >= minpeakheight &
#                   xv - pmin(x[x1], x[x2]) >= h_max &
#                   xv - pmax(x[x1], x[x2]) >= h_min)
#   X <- cbind(xv[inds], xp[inds], x1[inds], x2[inds])
#   
#   if (length(X) == 0) return(NULL)
#   # remove near point where dist < minpeakdistance
#   rm_near <- function(x){
#     I <- which(diff(x[, 2]) < minpeakdistance)[1]
#     if (is.na(I) | length(x) == 0){
#       return(x)
#     }else{
#       I_del <- I + as.integer(x[I, 1] > x[I+1, 1]) #remove the min
#       x <- x[-I_del, , drop = F]
#       if (nrow(x) <= 1) return(x);
#       rm_near(x)
#     }
#   }
#   
#   X <- X[order(X[, 2]), ,drop = F] # update 20180122; order according to index
#   X <- rm_near(X) # sort index is necessary before `rm_near`
#   
#   if (sortstr) { # || minpeakdistance > 1
#     sl <- sort.list(X[, 1], na.last = NA, decreasing = TRUE)
#     X <- X[sl, , drop = FALSE]
#   }
#   
#   if (npeaks > 0 && npeaks < nrow(X)) {
#     X <- X[1:npeaks, , drop = FALSE]
#   }
#   
#   X <- setNames(as.data.frame(X), c("val", "pos", "left", "right"))
#   
#   return(list(X))
# }

# ###########################################################
# Define 'smoothEO4PM' function
# ###########################################################

smoothEO4PM <- function(y, t, w=NULL, tout=NULL, win_noise=3, win_savgol=0, lambda=300, minValue=0, iters=3, spikes_removal="later", wFUN_type="wBisquare", fitMax=FALSE, raw=TRUE, plot=FALSE, plot_filename=NULL){
  
  # check input arguments
  if(length(y) != length(t)){
    stop("Object 't' should have the same length of object 'y'.")
  }
  if(!is.null(w)){
    if(length(y) != length(w)){
      stop("Object 'w' should have the same length of object 'y'.")
    }
  }  
  if(win_noise > 0){
    if(win_noise %% 2 == 0){
      stop("Temporal window to be used in spike removal set using argument 'win_noise' must be an odd positive integer.")
    }
  } else {
    if(win_noise < 0){
      stop("Temporal window to be used in spike removal set using argument 'win_noise' must be an odd positive integer or zero.")
    }
  }
  if(win_savgol > 0){
    if(win_savgol %% 2 == 0){
      stop("Temporal window to be used in Savistky-Golay smoothing set using argument 'w_savgol' must be an odd positive integer.")
    }
  } else {
    if(win_savgol < 0){
      stop("Temporal window to be used in Savistky-Golay smoothing set using argument 'w_savgol' must be an odd positive integer or zero.")
    }
  }
  if(!(wFUN_type %in% c("wTSM", "wBisquare", "wChen", "wSELF"))){
    stop("Argument 'wFUN_type' should be set to one of the followings: 'wBisquare' (dafault), 'wTSM', 'wChen', 'wSELF'.")
  }
  
  if(!(spikes_removal %in% c("none", "sync", "later"))){
    stop("Argument 'spikes_removal' must be set to one of the followings: 'none' (default), 'sync', 'later'.")
  }
  
  generate_plot <- FALSE
  if(!is.null(plot_filename)){
    generate_plot <- TRUE
    if(!dir.exists(dirname(plot_filename))){
      dir.create(dirname(plot_filename), showWarnings = FALSE, recursive = TRUE)
    }
  }
  if(plot){
    generate_plot <- TRUE
  }
  
  # convert input dates to integer
  t <- as.integer(as.Date(t, tz="UTC"))
  
  if(is.null(tout)){
    tout <- t[1]:t[length(t)]
  } else {
    tout <- sort(as.integer(as.Date(tout, tz="UTC")))
    if(tout[1] < t[1]){
      tout <- tout[which(tout >= t[1])]
      warning("Output dates start before the first available observation. Adjusting 'tout' argument accordingly.")
    }
    if(tout[length(tout)] > t[length(t)]){
      tout <- tout[which(tout <= t[length(t)])]
      warning("Output dates start before the first available observation. Adjusting 'tout' argument accordingly.")
    }
  }
  
  if(length(which(duplicated(t))) > 0){
    dupl_dates <- TRUE
  } else {
    dupl_dates <- FALSE
  }
  
  # rename input objects
  obs_dates_out <- tout
  pix <- y
  
  # define additional arguments
  savgol_method <- "phenofit"
  # savgol_method <- "sen2rts"
  
  # # load required libraries
  if (!("phenofit" %in% (.packages()))) {
    suppressPackageStartupMessages(library("phenofit", character.only = TRUE, warn.conflicts = FALSE))
  }
  if (!("imputeTS" %in% (.packages()))) {
    suppressPackageStartupMessages(library("imputeTS", character.only = TRUE, warn.conflicts = FALSE))
  }
  if(generate_plot){
    if (!("ggplot2" %in% (.packages()))) {
      suppressPackageStartupMessages(library("ggplot2", character.only = TRUE, warn.conflicts = FALSE))
    }
  }
  # if(savgol_method == "sen2rts"){
  #   suppressPackageStartupMessages(library(sen2rts))
  # }

  # wFUN_type <- "wBisquare"
  # wFUN_type <- "wTSM"
  wFUN_fun <- get(wFUN_type)
  
  # #####################################
  
  # remove NAs
  ppix <- c(1,(which(!is.na(pix[2:(length(pix)-1)]))+1), length(pix))
  pix <- pix[ppix]
  obsd <- t[ppix]
  
  # get weights
  if(!is.null(w)){
    w <- w[ppix]
  } else {
    w <- rep(1, length(pix))
  }
  
  if(raw | generate_plot){
    # initialize data frame to store series to plot
    dfp <- data.frame("time"=as.Date(obsd, tz="UTC"), "raw"=pix)
  }
  
  if(generate_plot){
    # initialize data frame to store series to plot
    sdfp <- data.frame("time"=as.Date(obs_dates_out, tz="UTC"))
  }
  
  # set temporal frame
  obs_dates_frame <- sort(obsd)[1]:sort(obsd, decreasing = TRUE)[1]
  
  # ######
  # remove small drops and spikes
  if(win_noise >= 3){
    
    if(spikes_removal == "later"){
      
      # find maxima in time series (to avoid spike removal when spike is a relative maximum)
      if(dupl_dates){
        result <- suppressWarnings(aggregate(pix, by=list(time=obsd), FUN=max, na.rm=TRUE))
        max_pos <- result$time
        result <- as.double(result$x)
        result[which(is.infinite(result))] <- NA
        
        ppix <- rep(as.double(NA), length(obs_dates_frame))
        ppix[which(obs_dates_frame %in% max_pos)] <- result
        ppix <- suppressMessages(imputeTS::na_interpolation(ppix, option = "linear"))
        # max_pos <- phenofit::findpeaks(ppix, minpeakdistance = 90, minpeakheight = quantile(pix, probs=0.75, na.rm=TRUE))
        max_pos <- findpeaks(ppix, minpeakdistance = 90, minpeakheight = quantile(pix, probs=0.75, na.rm=TRUE))
        result <- obs_dates_frame[max_pos[[1]]$pos]
        max_pos <- integer()
        for(i in result){
          ppix <- which(obsd %in% i)
          max_pos <- c(max_pos, which(pix == max(pix[ppix], na.rm=TRUE)))
        }
      } else {
        ppix <- rep(as.double(NA), length(obs_dates_frame))
        ppix[which(obs_dates_frame %in% obsd)] <- pix
        ppix <- suppressMessages(imputeTS::na_interpolation(ppix, option = "linear"))
        # max_pos <- phenofit::findpeaks(ppix, minpeakdistance = 90, minpeakheight = quantile(pix, probs=0.75, na.rm=TRUE))
        max_pos <- findpeaks(ppix, minpeakdistance = 90, minpeakheight = quantile(pix, probs=0.75, na.rm=TRUE))
        max_pos <- which(obsd %in% obs_dates_frame[max_pos[[1]]$pos])
      }
      
      # remove small drops and spikes
      ppix <- noiseRemoval(y=pix, w=win_noise, spikes=spikes_removal, borders = TRUE, max_pos=max_pos)
      
    } else {
      
      # remove small drops and spikes
      ppix <- noiseRemoval(y=pix, w=win_noise, spikes=spikes_removal, borders = TRUE)
      
    }
    
    # update weights
    w[(which(pix != ifelse(is.na(ppix), -999, ppix)))] <- 0.0
    
    # # update max_pos (may be useful for Whittaker 2nd pass)
    # max_pos <- max_pos[which(pix[max_pos] == ifelse(is.na(ppix[max_pos]), -999, ppix[max_pos]))]
    
    pix <- ppix
    
    if(plot){
      dfp <- cbind(dfp, data.frame("despike"=pix))
    }
  }
  
  # ######
  # create daily aggregates
  if(dupl_dates){
    ppix <- suppressWarnings(aggregate(pix, by=list(time=obsd), FUN=mean, na.rm=TRUE))
    pix <- as.double(ppix$x)
    ppix <- suppressWarnings(aggregate(w, by=list(time=obsd), FUN=max, na.rm=TRUE))
    w <- ppix$x
    obsd <- ppix$time
  }
  
  if(all(is.na(pix))){
    result <- rep(NA, length(obs_dates_out))
  } else {
    
    # set weights update function
    wFUN_fun <- get(wFUN_type)
    
    # ######
    # Time series interpolation (using stine)
    smoothed_pixel <- rep(NA, length(obs_dates_out))
    smoothed_pixel[which(obs_dates_out %in% obsd)] <- pix[which(obsd %in% obs_dates_out)]
    
    # interpolate using stine
    ppix <- suppressMessages(imputeTS::na_interpolation(smoothed_pixel, option = "stine"))
    
    if(generate_plot){
      sdfp <- cbind(sdfp, data.frame("stine"=ppix))
    }
    
    # ######
    # Time series local smoothing (using Savitski-Golay)
    if(win_savgol >= 3){
      
      # update weights
      wt <- rep(0.1, length(obs_dates_out))
      wt[which(obs_dates_out %in% obsd)] <- w[which(obsd %in% obs_dates_out)]
      
      # smoothing using Savitsky-Golay (phenofit implementation)
      if(savgol_method == "phenofit"){
        # if(is.null(minValue)){
        #   minValue_ylu <- NA
        # } else {
        #   minValue_ylu <- minValue
        # }
        if(is.null(minValue)){
          ylu <- quantile(smoothed_pixel[which(smoothed_pixel != -999)], c(0.01,0.99))
          # ylu[1] <- pmax(ylu[1], 0, minValue_ylu, na.rm=TRUE)
          ylu[2] <- max(smoothed_pixel, na.rm=TRUE)
        } else {
          ylu <- c(minValue, max(smoothed_pixel, na.rm=TRUE))
        }
        
        ppix <- phenofit::smooth_wSG(y=ppix, w=wt, nptperyear = 365, d = 3, frame = win_savgol, iters=2, ylu=ylu, wFUN=wFUN_fun)
        smoothed_pixel <- as.double(ppix[[1]][[length(ppix[[1]])]])
        wt <- as.double(ppix[[2]][[length(ppix[[2]])]])
        # wt <- phenofit::wTSM(y = ppix, yfit = smoothed_pixel, w = wt, iter = 1, nptperyear = 365, wfact = 0.5)
      }
      
      # smoothing using Savitsky-Golay (sen2rts implementation)
      if(savgol_method == "sen2rts"){
        # # interpolate head and tail
        # smoothed_pixel <- w_savgol(y=ppix, x=obs_dates_out, q=wt, polynom = 3, window = win_savgol, borders = TRUE)
        smoothed_pixel <- sen2rts:::w_savgol(y=ppix, x=obs_dates_out, q=wt, polynom = 3, window = win_savgol)
        smoothed_pixel[which(is.na(smoothed_pixel))] <- pix[which(is.na(smoothed_pixel))]
        wt <- phenofit::wTSM(y = ppix, yfit = smoothed_pixel, w = wt, iter = 1, nptperyear = 365, wfact = 0.5)
      }
      
      if(generate_plot){
        sdfp <- cbind(sdfp, data.frame("savgol"=smoothed_pixel))
      }
      
    } else {
      smoothed_pixel <- ppix
      wt <- rep(0, length(obs_dates_out))
      wt[which(obs_dates_out %in% obsd)] <- 1
    }
    
    # clean values lower than minValue
    if(!is.null(minValue)){
      wt[which(smoothed_pixel < minValue)] <- wt[which(smoothed_pixel < minValue)] * 0.5
      smoothed_pixel[which(smoothed_pixel < minValue)] <- minValue
    }
    
    # #####
    # Global temporal smoothing (using Whittaker)
    if(lambda == 0){
      w_lambda <- NULL
    } else {
      w_lambda <- lambda
    }
    ylu <- quantile(smoothed_pixel, c(0.01,0.99), na.rm=TRUE)
    # if(is.null(minValue)){
    #   minValue_ylu <- NA
    # } else {
    #   minValue_ylu <- minValue
    # }
    
    # update weights to account for low values
    wt[which(smoothed_pixel <= ylu[1])] <- 1
    
    if(is.null(minValue)){
      # ylu[1] <- pmax(ylu[1], 0, minValue_ylu, na.rm=TRUE)
      ylu[2] <- max(smoothed_pixel, na.rm=TRUE)
    } else {
      ylu <- c(minValue, max(smoothed_pixel, na.rm=TRUE))
    }
    
    # # replace missing values in smoothed_pixel
    # wt[which(is.na(smoothed_pixel))] <- 0
    # smoothed_pixel[which(is.na(smoothed_pixel))] <- -999
    
    # wt_zero <- which(wt == 0)
    # result <- smoothed_pixel
    # for(i in 1:iters){
    #   result <- phenofit::smooth_wWHIT(y=result, w=wt, ylu=ylu, nptperyear = 365, d = 3, iters=2, lambda=w_lambda, wFUN=wFUN_fun, second=TRUE)
    #   wt <- as.double(result[[2]][[length(result[[2]])]])
    #   wt[wt_zero] <- 0
    #   result <- as.double(result[[1]][[length(result[[1]])]])
    # }
    
    result <- phenofit::smooth_wWHIT(y=smoothed_pixel, w=wt, ylu=ylu, nptperyear = 365, d = 3, iters=3, lambda=w_lambda, wFUN=wFUN_fun, second=TRUE)
    wt <- as.double(result[[2]][[length(result[[2]])]])
    result <- as.double(result[[1]][[length(result[[1]])]])
    
    if(generate_plot){
      sdfp <- cbind(sdfp, data.frame("whit"=result))
    }
    
    # ######
    # Temporal smoothing pass 2 - Improve fit of maximum values
    
    if(fitMax){
      # identify anomaly in interpolated annual maximum (10% difference compared to original value)
      max_delta_thres <- 0.05
      max_delta_perc_selection <- 1.0
      fix_obsd <- integer()
      pix_year <- as.integer(format(as.Date(obsd, origin="1970-01-01"), "%Y"))
      interp_pix_year <- as.integer(format(as.Date(obs_dates_out, origin="1970-01-01"), "%Y"))
      year_list <- sort(unique(pix_year))
      for(y in year_list){
        yd <- as.integer(max(obsd[which(pix_year %in% y)], na.rm=TRUE) - min(obsd[which(pix_year %in% y)], na.rm=TRUE))
        if(yd >= 180){
          pix_max_value <- suppressWarnings(max(pix[which(pix_year %in% y)], na.rm=TRUE))
          pix_max_interp <- suppressWarnings(max(result[which(interp_pix_year %in% y)], na.rm=TRUE))
          if(abs(pix_max_value - pix_max_interp) > (pix_max_value * max_delta_thres)){
            fix_obsd <- c(fix_obsd, obsd[which(pix_year %in% y)][which(pix[which(pix_year %in% y)] >= (pix_max_value * max_delta_perc_selection))])
          }
        }
      }
      
      # perform smoothing
      if(length(fix_obsd) > 0){
        
        # overwrite annual maximum values resulted in wrong smoothing
        smoothed_pixel[which(obsd %in% fix_obsd)] <- pix[which(obsd %in% fix_obsd)]
        wt <- wt * 0.5
        wt[which(obsd %in% fix_obsd)] <- 1
        
        # result <- smoothed_pixel
        # for(i in 1:iters){
        #   result <- phenofit::smooth_wWHIT(y=result, w=wt, ylu=ylu, nptperyear = 365, d = 3, iters=2, lambda=w_lambda, wFUN=wFUN_fun, second=TRUE)
        #   wt <- as.double(result[[2]][[length(result[[2]])]])
        #   wt[wt_zero] <- 0
        #   result <- as.double(result[[1]][[length(result[[1]])]])
        # }
        
        result <- phenofit::smooth_wWHIT(y=smoothed_pixel, w=wt, ylu=ylu, nptperyear = 365, d = 3, iters=3, lambda=w_lambda, wFUN=wFUN_fun, second=TRUE)
        wt <- as.double(result[[2]][[length(result[[2]])]])
        result <- as.double(result[[1]][[length(result[[1]])]])
        
        if(generate_plot){
          sdfp <- cbind(sdfp, data.frame("whit_fitMax"=result))
        }
      }
    }
    
  }
  
  # clean values lower than minValue
  if(!is.null(minValue)){
    result[which(result < minValue)] <- minValue
  }
  
  # generate plot
  if(generate_plot){
    # initialize plot
    g <- ggplot2::ggplot()
    g <- g + ggplot2::geom_line(data=sdfp, ggplot2::aes(time, stine), linetype="solid", color="cyan3", linewidth=0.4)
    
    if("savgol" %in% names(sdfp)){
      g <- g + ggplot2::geom_line(data=sdfp, ggplot2::aes(time, savgol), linetype="solid", color="blue4", linewidth=0.4)
    }
    
    if("whit_fitMax" %in% names(sdfp)){
      g <- g + ggplot2::geom_line(data=sdfp, ggplot2::aes(time, whit), linetype="solid", color="orange4", linewidth=0.7)
      g <- g + ggplot2::geom_line(data=sdfp, ggplot2::aes(time, whit_fitMax), linetype="solid", color="forestgreen", linewidth=0.6)
    } else {
      g <- g + ggplot2::geom_line(data=sdfp, ggplot2::aes(time, whit), linetype="solid", color="forestgreen", linewidth=0.6)
    }
    
    # plot raw points
    g <- g + ggplot2::geom_point(data=dfp, ggplot2::aes(time, raw), colour="red2", fill="red2", size=0.3, na.rm=TRUE)
    # plot points after spike removal
    if("despike" %in% names(dfp)){
      g <- g + ggplot2::geom_point(data=dfp, ggplot2::aes(time, despike), colour="black", fill="black", size=0.3, na.rm=TRUE)
    }
    
    # add extra graphics
    g <- g + ggplot2::theme_bw() + ggplot2::ylab("") + ggplot2::xlab("Time")
    # ymax <- max(l3a_profile$value[which(as.integer(strftime(l3a_profile$time, format="%Y")) == year)], na.rm=TRUE)*1.05
    # g <- g + ggplot2::ylim(0, ymax)
    g <- g + ggplot2::scale_x_date(limits=c(as.Date(tout[1], tz="UTC"), as.Date(tout[length(tout)], tz="UTC")), date_labels="%b-%y")
    
    # display plot
    if(!is.null(plot_filename)){
      # printo to file
      png(filename=plot_filename, width = 3000, height = 1500, units = "px", pointsize = 12, bg = "white", res = 300)
      suppressWarnings(print(g))
      invisible(dev.off())
    }
    if(plot){
      suppressWarnings(print(g))
    }
  }
  
  result <- data.frame("t"=as.integer(obs_dates_out), "y"=as.double(result), "w"=as.double(wt))
  if(raw){
    result <- cbind(result, data.frame("raw"=rep(as.double(NA), nrow(result))))
    if(dupl_dates){
      ppix <- suppressWarnings(aggregate(as.double(dfp$raw), by=list(time=dfp$time), FUN=mean, na.rm=TRUE))
      dfp <- data.frame("time"=as.integer(ppix$time), "raw"=as.double(ppix$x))
    }
    result$raw[which(result$t %in% as.integer(dfp$time))] <- as.double(dfp$raw[which(as.integer(dfp$time) %in% result$t)])
  }
  
  return(result)
}
