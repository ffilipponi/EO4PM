#!/usr/bin/env Rscript
# ########################################################
# Title         : EO4PM.R
# Description   : Extract Phenological Metrics from daily time series
# Date          : Jul 2026
# Version       : 2.1
# Copyright name: CC BY-NC-SA
# Licence       : GPL v3
# Authors       : Federico Filipponi
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
# ########################################################

# ### TODO
# - integrate custom functions in the script
# - create object to store temporal statistics
# OK valuta utilizzo variabile alternativa per CRS (es. WKT o altro)
# - aggiorna global attributes
# - aggiorna versione GCMD Science Keywords, Version 8.7
### verifica se c'è già questo: set history to get all the attributes from 'multitemporal_stacking_processor' and 'multitemporal_smoothing_processor'
# - aggiungi sistema con units per gestire da solo il formato data in input (es. epoch secondi o epoch giorni)
# - load mask from external file (both need terra::compareGeom() check, e opzionalmente crop mask to input file o generala da vettoriale)
# - aggiungi possibilità clip a bounding box (anche in LatLong)
# - verifica QF e scrivi funzione per leggerne facilmente i valori
# - aggiungi più statistiche stagionali

# load optparse library for argument parsing
options("repos"="http://cran.at.r-project.org/")
invisible(tryCatch(find.package("optparse"), error=function(e) install.packages("optparse", dependencies=TRUE)))
require(optparse, quietly=TRUE)

# read arguments
option_list <- list(
  optparse::make_option(c("-i", "--input"), type="character", default=NULL, 
                        help="Input dataset file name", metavar="file"),
  optparse::make_option(c("-o", "--output"), type="character", default=getwd(), 
                        help="Output path [default = %default]", metavar="folder"),
  optparse::make_option(c("-c", "--calendar"), type="character", default="natural", 
                        help="Calendar type (must be 'crop', 'extended', 'user' or 'natural' [default = %default])", metavar="character"),
  optparse::make_option(c("-d", "--start_date"), type="character", default=NULL, 
                        help="Start date (es. '2017-01-01')", metavar="character"),
  optparse::make_option(c("-f", "--end_date"), type="character", default=NULL, 
                        help="End date (es. '2017-12-31')", metavar="character"),
  optparse::make_option(c("-y", "--year"), type="integer", default=NULL, 
                        help="Year to be used for the analysis", metavar="integer"), 
  optparse::make_option(c("-m", "--mask"), type="character", default=NULL, 
                        help="Input mask file path", metavar="file"),
  optparse::make_option(c("-v", "--variable"), type="character", default=NULL, 
                        help="Variable name to be processed stored in input file", metavar="character"),
  optparse::make_option(c("-s", "--seasons"), type="integer", default=1, 
                        help="Maximum number of seasons to find in a solar year (positive integer in range 1-9) [default = %default]", metavar="integer"),
  optparse::make_option(c("-b", "--min_value"), type="double", default=NULL, 
                        help="Minimum allowed baseline value in cut cycle", metavar="double"),
  optparse::make_option(c("--prefix"), type="character", default=NULL, 
                        help="Prefix to be assigned to output file name", metavar="character"),
  optparse::make_option(c("--chunksize"), type="integer", default=128, 
                        help="Chunksize to optimize RAM usage [default = %default]", metavar="integer"),
  optparse::make_option(c("--deflate"), type="integer", default=6, 
                        help="NetCDF file compression (should be an integer value between '0' and '9' [default = %default])", metavar="integer"),
  optparse::make_option(c("-q","--cores"), type="integer", default=1, 
                        help="Number of cores to use [default = %default]", metavar="integer"),
  # optparse::make_option(c("--compute_tstats"), type="logical", default=FALSE, action="store_true",
  #                       help="Compute temporal statistics (optional)"),
  optparse::make_option(c("--force"), type="logical", default=FALSE, action="store_true",
                        help="Force computation (if there are less than 365 observations in the time window, or force the use asymmetric gaussian fit is Gu estimate fails)"),
  optparse::make_option(c("--recycling"), type="character", default=NULL, 
                        help="Force recycling of an existing output NetCDF file (use in combination with '--mask' argument is adviced)", metavar="file"),
  optparse::make_option(c("--progress"), type="logical", default=FALSE, action="store_true",
                       help="Display progress bar"),
  optparse::make_option(c("--verbose"), type="logical", default=FALSE, action="store_true",
                        help="Verbose mode"),
  optparse::make_option(c("--c_win"), type="integer", default=20, 
                        help="Window in days corresponding to the minumum duration of increasing and decreasing around peak [default = %default]", metavar="integer"),
  optparse::make_option(c("--s_lag"), type="integer", default=5, 
                        help="Head and tail (in days) to extend vegetation cycle in phenology calculation [default = %default]", metavar="integer"),
  optparse::make_option(c("--maxExtendMonth"), type="integer", default=4, 
                        help="Maximum number of months to extend time series in vegetation cycle identification [default = %default]", metavar="integer"),
  optparse::make_option(c("--minPeakValue"), type="double", default=1.0, 
                        help="Minimum peak value (in variable units) [default = %default]", metavar="double"),
  optparse::make_option(c("--minPeakDistance"), type="integer", default=91, 
                        help="Minimum peak distance in days [default = %default]", metavar="integer"),
  optparse::make_option(c("--minSeasonLength"), type="integer", default=45, 
                        help="Minimum season length in days [default = %default]", metavar="integer"),
  optparse::make_option(c("--r_min"), type="double", default=0.05, 
                        help="Minimum relative height difference between trough and peak values [default = %default]", metavar="double"),
  optparse::make_option(c("--r_max"), type="double", default=0.2, 
                        help="Maximum relative height difference between trough and peak values [default = %default]", metavar="double"),
  optparse::make_option(c("--rtrough_max"), type="double", default=0.5, 
                        help="Maximum ratio between trough and peak values [default = %default]", metavar="double")
)
opt_parser <- optparse::OptionParser(option_list=option_list)
opt <- optparse::parse_args(opt_parser)

# check if required arguments are supplied and point to existing files
if(is.null(opt$input)){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (input).", call.=FALSE)
}
if(is.null(opt$recycling)){
  if(is.null(opt$output)){
    optparse::print_help(opt_parser)
    stop("At least one argument must be supplied (output) if argument '--recycling' is not provided.", call.=FALSE)
  }
}
if(!(opt$calendar == "user")){
  if(is.null(opt$year)){
    optparse::print_help(opt_parser)
    stop("At least one argument must be supplied (year).", call.=FALSE)
  }
}

############################################################
# Setup environment
############################################################

# define R objects from options
input_file <- normalizePath(path=opt$input, winslash="/", mustWork=TRUE)
output_path <- normalizePath(path=opt$output, winslash="/", mustWork=FALSE)
variable <- opt$variable
max_season <- opt$seasons
minValue <- opt$min_value
calendar <- opt$calendar
year <- as.integer(opt$year)
nc_mask_arg <- opt$mask
out_fname_prefix <- opt$prefix
chunksize <- opt$chunksize
nc_compression <- opt$deflate
cores <- opt$cores
# compute_tstats <- opt$compute_tstats
recycling <- opt$recycling
force <- opt$force
progress <- opt$progress
verbose <- opt$verbose
# parameters
nptperyear <- 365
c_win <- opt$c_win
s_lag <- opt$s_lag
maxExtendMonth <- opt$maxExtendMonth
ypeak_min <- opt$minPeakValue
minpeakdistance <- opt$minPeakDistance
length_min <- opt$minSeasonLength
r_min <- opt$r_min
r_max <- opt$r_max
rtrough_max <- opt$rtrough_max

if(!is.null(recycling)){
  recycling <- normalizePath(path=opt$recycling, winslash="/", mustWork=TRUE)
} else {
  # check if supplied output path is a file
  if(file.exists(output_path)){
    if(file_test("-f", output_path)){
      stop("Supplied output path is an existing file")
    }
  } else {
    # create output folder if does not exists
    dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  }
}

# check if mask file exists
if(!is.null(nc_mask_arg) && nc_mask_arg != "internal"){
  nc_mask_file <- normalizePath(path=nc_mask_arg, winslash="/", mustWork=FALSE)
  if(!file.exists(nc_mask_file)){
    stop(paste("Mask file '", nc_mask_file,"' does not exists. Set another mask file name", sep=""))
  }
}

# ##########################
# check input arguments

# check calendar type
if(!(calendar %in% c("natural","crop","extended","user"))){
  stop("Calendar type must be either 'crop', 'extended', 'user' or 'natural'")
}

if(calendar == "user"){
  # define function to check date format
  IsDate <- function(mydate, date.format = "%Y-%m-%d"){
    tryCatch(!is.na(as.Date(mydate, date.format)),  
             error = function(err) {FALSE})  
  }
  if(is.null(opt$start_date)){
    optparse::print_help(opt_parser)
    stop("At least one argument must be supplied (start_date).", call.=FALSE)
  } else {
    # check if date has the correct format
    if(!IsDate(opt$start_date)){
      stop("Argument 'start_date' has not the format 'YYYY-MM-DD'.", call.=FALSE)
    }
  }
  if(!is.null(opt$end_date)){
    # check if date has the correct format
    if(!IsDate(opt$end_date)){
      stop("Argument 'end_date' has not the format 'YYYY-MM-DD'.", call.=FALSE)
    }
    if(as.Date(opt$start_date) > as.Date(opt$end_date)){
      stop("Argument 'start_date' must be set to a date preceding argument 'end_date'.", call.=FALSE)
    }
  }
  if((as.integer(as.Date(opt$end_date)) - as.integer(as.Date(opt$start_date))) > 500){
    warning("Temporal range between 'start_date' and 'end_date' is longer than 500 days.")
  }
  start_date <- as.Date(opt$start_date)
  end_date <- as.Date(opt$end_date)
  if(!is.null(year)){
    if(!(year %in% as.integer(strftime(x = seq(start_date, end_date, 'days'), format="%Y", tz = "UTC", usetz = FALSE)))){
      year_rle <- rle(as.integer(strftime(x = seq(start_date, end_date, 'days'), format="%Y", tz = "UTC", usetz = FALSE)))
      year <- as.integer(year_rle$values[order(year_rle$lengths, decreasing = TRUE)[1]])
    }
  } else {
    year_rle <- rle(as.integer(strftime(x = seq(start_date, end_date, 'days'), format="%Y", tz = "UTC", usetz = FALSE)))
    year <- as.integer(year_rle$values[order(year_rle$lengths, decreasing = TRUE)[1]])
  }
}

if(max_season < 1 | max_season > 9){
  stop("Max season argument muste be a positive integer number between 1 and 9.")
}

if(!(nc_compression %in% c(0:9))){
  stop("NetCDF file compression set using '--deflate' argument should be set to an integer value between '0' and '9'.")
} else {
  if(nc_compression == 0){
    nc_compression <- NA
  }
}

if(chunksize < 0){
  stop("Chunksize must be a positive integer number.")
}

# ###################################
# ### for debug
# # input_file <- "/space/T32TQQ_20150704_20240123_S2_L3B_20m_LAI_IdS.nc"
# # input_file <- "/mnt/repository/nr_products/sentinel2/l3a_tsg/T32TQQ/T32TQQ_20150704_20260127_S2_L3A_TSG_10m_LAI.nc"
# input_file <- "/media/data/workspace/T32TQQ_20150704_20251031_S2_L3A_TSG_10m_LAI.nc"
# output_path <- "/space/workspace/EO4PMv2/out01_test72"
# nc_mask_arg <- NULL
# # nc_mask_arg <- "/space/data/vector/ESU_LAI_2021.2022.2023_mask.tif"
# year <- 2021
# calendar <- "extended"
# variable <- NULL
# nc_compression <- 9
# chunksize <- 16
# # chunksize <- 128
# cores <- 24
# compute_tstats <- FALSE
# # force <- FALSE
# force <- TRUE
# # progress <- TRUE
# progress <- FALSE
# verbose <- TRUE
# # parameters
# max_season <- 3
# # max_season <- 1
# minValue <- NA
# s_lag <- 5
# nptperyear <- 365
# maxExtendMonth <- 4
# minpeakdistance <- 91
# length_min <- 45
# # rtrough_max <- 0.6
# rtrough_max <- 0.5
# ypeak_min <- 1.0
# r_min <- 0.05
# r_max <- 0.2
# # # per detection medica
# # r_min <- 0.01
# # r_max <- 0.1
# source("/space/analysis/IREA/auxiliary_files/codes/EO4PM/scripts/phenoGetEO4PM.R")

# ###################################
# load required libraries
invisible(tryCatch(find.package("sf"), error=function(e) install.packages("sf", dependencies=TRUE)))
invisible(tryCatch(find.package("ncdf4"), error=function(e) install.packages("ncdf4", dependencies=TRUE)))
invisible(tryCatch(find.package("terra"), error=function(e) install.packages("terra", dependencies=TRUE)))
invisible(tryCatch(find.package("data.table"), error=function(e) install.packages("data.table", dependencies=TRUE)))
invisible(tryCatch(find.package("phenofit"), error=function(e) install.packages("phenofit", dependencies=TRUE)))
invisible(tryCatch(find.package("progressr"), error=function(e) install.packages("progressr", dependencies=TRUE)))
invisible(tryCatch(find.package("future.callr"), error=function(e) install.packages("future.callr", dependencies=TRUE)))
invisible(tryCatch(find.package("doFuture"), error=function(e) install.packages("doFuture", dependencies=TRUE)))
suppressWarnings(suppressPackageStartupMessages(require(future.callr, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(require(doFuture, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(require(progressr, quietly=TRUE)))
suppressPackageStartupMessages(require(sf, quietly=TRUE))
suppressPackageStartupMessages(require(ncdf4, quietly=TRUE))
suppressPackageStartupMessages(require(terra, quietly=TRUE))
suppressPackageStartupMessages(require(data.table, quietly=TRUE))
suppressPackageStartupMessages(require(phenofit, quietly=TRUE))

terra::terraOptions(todisk=TRUE)
terra::terraOptions(progress=0)
terra::setGDALconfig("IGNORE_XY_AXIS_NAME_CHECKS","YES")

readRast_opts <- NULL

# setup progress bar
if(progress){
  options(progressr.enable = TRUE)
  progressr::handlers(global = TRUE)
  # progressr::handlers("cli")
}

options(warn=-1)

############################
# # get environmental variables
# # PROCESSOR_HOME <- Sys.getenv("PROCESSOR_HOME")
# initial.options <- commandArgs(trailingOnly = FALSE)
# file.arg.name <- "--file="
# script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
# PROCESSOR_HOME <- normalizePath(path=dirname(script.name), winslash = "/", mustWork = TRUE)
# 
# # ### for debug
# # source("/space/analysis/IREA/auxiliary_files/codes/EO4PM/scripts/phenoGetEO4PM.R")
# # PROCESSOR_HOME <- "/space/scripts"
# 
# # Load processing functions from external file
# source(normalizePath(path=paste(PROCESSOR_HOME, "/phenoGetEO4PM.R", sep="")))

tmp_folder <- tempdir()

ptm <- proc.time()
nc_hist_systime <- as.character(format(Sys.time(), "%a %b %d %H:%M:%S %Y"))

# ###########################################################
# Define function to compute Phenological Metrics
# ###########################################################

# v 2.1
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

# define function to extract time_set
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

# define smoothing function compatible with phenofit to skip smoothing in cut_cycle
smooth_none <- function(y, w, ylu, nptperyear, ...) {
  y <- as.numeric(y)
  w <- as.numeric(w)
  list(
    zs = list(ziter1 = y),
    ws = list(witer1 = w)
  )
}

# # check if jemalloc library is used
# system(paste0("grep jemalloc /proc/", Sys.getpid(), "/maps"))

# ###########################################################
# Print parameters
# ###########################################################

if(verbose){
  message("================================================================================")
  message("Running 'EO4PM' processor for Phenological Metrics estimation - Version 2.1")
  message(paste("Processing file: '", input_file, "'", sep=""))
  message(paste("Number of cores used for the processing is set to: '", cores, "'", sep=""))
  message("")
  message(paste("Using mask: ", ifelse(is.null(nc_mask_arg), "NONE", nc_mask_arg), sep=""))
  if(!(calendar == "user")){
    message(paste("Analysis year is: ", year, sep=""))
  }
  message(paste("Calendar type is: '", calendar, "'", ifelse(calendar == "user", paste(" [", start_date, " - ", end_date, "]", sep=""), ""), sep=""))
  # ### TODO add extra parameter print
  message("")
  if(is.null(recycling)){
    message(paste("Results will be stored in path: '", output_path, "'", sep=""))
  } else {
    message(paste("Results will recycle file: '", recycling, "'", sep=""))
  }
  message(paste("Starting Phenological Metrics estimation: ", Sys.time(), sep=""))
}

############################################################
# Setup processor
############################################################

# read NetCDF file
nc_in <- ncdf4::nc_open(filename = input_file, write = FALSE, readunlim = FALSE)
nc_file_name <- basename(input_file)

# check input NetCDF variable names
nc_var_name <- names(nc_in$var)[!(names(nc_in$var) %in% c("crs", "weights","mask","time_obs"))]
if(is.null(variable)){
  if(length(nc_var_name) != 1){
    ncdf4::nc_close(nc_in)
    stop("Input file do not contain a unique variable")
  } else {
    var_name <- nc_var_name
  }
} else {
  if(!(variable %in% nc_var_name)){
    ncdf4::nc_close(nc_in)
    stop(paste("Variable name '", variable, "'defined using argument '--variable' is not available in input NetCDF file.", sep=""))
  } else {
    var_name <- variable
  }
}

# get global attributes
nc_attrs <- ncdf4::ncatt_get(nc_in, 0)
# get variable attributes
nc_var_atts <- ncdf4::ncatt_get(nc_in, var_name)

# read input NetCDF using terra
nc_rast <- suppressWarnings(terra::rast(input_file, opts=readRast_opts))

# check if SpatRaster object has time dimension
if(!terra::has.time(nc_rast)){
  stop("EO4PM requires input file with time dimension.")
}

# get time information from NetCDF file
# tm_int <- ncdf4::ncvar_get(nc_in,"time")
# obs_dates <- as.integer(tm_int)
# get time information from NetCDF file (using terra library)
tm_int <- terra::time(nc_rast)
obs_dates <- as.integer(as.Date(tm_int, tz = "UTC"))

# check if date spacing is daily
if(!(length(seq(from=obs_dates[1], to=obs_dates[length(obs_dates)], by = 1)) == length(obs_dates))){
  stop("EO4PM requires input temporal profiles with daily spacing.")
}

# check if input NetCDF has spatial CRS
if(!(c("crs") %in% names(nc_in$var))){
  if(verbose){
    warning("Input NetCDF has not a variable named 'crs' containing spatial 'Coordinate Reference System'. Going on however.")
  }
} else {
  if(!(ncdf4::ncatt_get(nc_in, "crs", "spatial_ref")$hasatt)){
    warning("Variable 'crs' in input NetCDF file does not contain attribute 'spatial_ref' for the definition of 'Coordinate Reference System'. Going on however.")
  }
}

# set the temporal window to be used for phenological metrics estimation
if(calendar == "natural"){
  start_date <- as.integer(as.Date(strptime(paste(year, "-01-01", sep=""), format = "%Y-%m-%d", tz = "GMT"), origin="1970-01-01"))
  end_date <- as.integer(as.Date(strptime(paste(year+1, "-01-01", sep=""), format = "%Y-%m-%d", tz = "GMT"), origin="1970-01-01"))-1
}
if(calendar == "crop"){
  start_date <- as.integer(as.Date(strptime(paste(year-1, "-10-01", sep=""), format = "%Y-%m-%d", tz = "GMT"), origin="1970-01-01"))
  end_date <- as.integer(as.Date(strptime(paste(year, "-11-01", sep=""), format = "%Y-%m-%d", tz = "GMT"), origin="1970-01-01"))-1
}
if(calendar == "extended"){
  start_date <- as.integer(as.Date(strptime(paste(year-1, "-09-01", sep=""), format = "%Y-%m-%d", tz = "GMT"), origin="1970-01-01"))
  end_date <- as.integer(as.Date(strptime(paste(year+1, "-01-01", sep=""), format = "%Y-%m-%d", tz = "GMT"), origin="1970-01-01"))-1
}
if(calendar == "natural"){
  start_date <- as.integer(start_date)
  end_date <- as.integer(end_date)
}

# define temporal window for the analysis
start_pos <- max(1, which(obs_dates >= start_date)[1])
end_pos <- min(length(obs_dates), sort(which(obs_dates <= end_date), decreasing = TRUE)[1])

if(length(which(obs_dates >= start_date & obs_dates <= end_date)) == 0){
  stop(paste("Temporal observations in NetCDF dataset (", as.Date(obs_dates[1], origin="1970-01-01"), " - ", as.Date(obs_dates[length(obs_dates)], origin="1970-01-01"), ") falls outside analysis year '", year, "'", sep=""))
}

# set import time window
pheno_lag <- 0

if(is.null(start_date) | is.null(end_date)){
  nc_time_start <- 1
  nc_time_count <- -1
} else {
  nc_time_start <- max(1, start_pos - pheno_lag)
  nc_time_count <- min(length(obs_dates), end_pos + pheno_lag) - nc_time_start + 1

  if(nc_time_count < 365){
    if(force){
      if(verbose){
        warning("Input NetCDF has few available obseravations. Continuing in '--force' mode.")
      }
    } else {
      if(start_date < obs_dates[1] | end_date > obs_dates[length(obs_dates)]){
        stop(paste("Temporal observations in NetCDF dataset (", as.Date(obs_dates[1], origin="1970-01-01"), " - ", as.Date(obs_dates[length(obs_dates)], origin="1970-01-01"), ") does not fully include analysis year '", year, "'", sep=""))
      }
      stop(paste("Input NetCDF has few available observations.", sep=""))
    }
  }
}
# define phenodates
phenodates <- obs_dates[start_pos:end_pos]

# # define if phenodates full include season for DJF computation
# if(sum(as.integer(phenodates < as.integer(as.Date(strptime(paste(year, "-03-01", sep=""), format = "%Y-%m-%d", tz = "GMT"), origin="1970-01-01")))) > 1){
#   if(phenodates[1] > as.integer(as.Date(strptime(paste(year-1, "-12-01", sep=""), format = "%Y-%m-%d", tz = "GMT"), origin="1970-01-01"))){
#     full_seasonal_djf_flag <- FALSE
#     seasonal_djf_bound1 <- as.character(as.Date(phenodates[1], origin="1970-01-01"))
#   } else {
#     full_seasonal_djf_flag <- TRUE
#     seasonal_djf_bound1 <- as.character(as.Date(strptime(paste(year-1, "-12-01", sep=""), format = "%Y-%m-%d", tz = "GMT"), origin="1970-01-01"))
#   }
#   seasonal_djf_bound2 <- as.character(as.Date(as.integer(as.Date(strptime(paste(year, "-03-01", sep=""), format = "%Y-%m-%d", tz = "GMT"), origin="1970-01-01"))-1, origin="1970-01-01"))
#   seasonal_djf_bounds_attr <- paste(seasonal_djf_bound1, seasonal_djf_bound2, sep=" - ")
# } else {
#   full_seasonal_djf_flag <- FALSE
#   seasonal_djf_bounds_attr <- "NONE"
# }

# create temporal dimension
nc_out_time_dim <- as.integer(obs_dates)
# set start and end dates
nc_att_date_start <- format(as.Date(nc_out_time_dim[1], origin="1970-01-01"), "%Y-%m-%d")
nc_att_date_end <- format(as.Date(nc_out_time_dim[length(nc_out_time_dim)], origin="1970-01-01"), "%Y-%m-%d")
if(calendar == "natural"){
  nc_file_date_start <- as.integer(paste(year, "0101", sep=""))
  nc_file_date_end <- as.integer(paste(year, "1231", sep=""))
}
if(calendar == "crop"){
  nc_file_date_start <- as.integer(paste(year-1, "1101", sep=""))
  nc_file_date_end <- as.integer(paste(year, "1031", sep=""))
}
if(calendar == "extended"){
  nc_file_date_start <- as.integer(paste(year-1, "0901", sep=""))
  nc_file_date_end <- as.integer(paste(year, "1231", sep=""))
}
if(calendar == "user"){
  nc_file_date_start <- as.integer(paste(strftime(start_date, format="%Y%m%d", tz="UTC", usetz=FALSE), sep=""))
  nc_file_date_end <- as.integer(paste(strftime(end_date, format="%Y%m%d", tz="UTC", usetz=FALSE), sep=""))
}

# get variable dimensions
var_dims_names <- character()
var_dims <- nc_in$var[[var_name]]$dim
for(d in 1:length(var_dims)){
  var_dims_names <- c(var_dims_names, var_dims[[d]]$name)
}
# set dimension names
xname <- "x"
yname <- "y"
if("easting" %in% var_dims_names && "northing" %in% var_dims_names){
  xname <- "easting"
  yname <- "northing"
}
if("longitude" %in% var_dims_names && "latitude" %in% var_dims_names){
  xname <- "longitude"
  yname <- "latitude"
}
if("lon" %in% var_dims_names && "lat" %in% var_dims_names){
  xname <- "lon"
  yname <- "lat"
}
if(!("x" %in% var_dims_names && "y" %in% var_dims_names) && xname == "x" && yname == "y"){
  warning("Neither 'easting and 'northing' nor 'x' and 'y' dimension names found. Going on using 'x' and 'y'.")
}

# get dimensions
# nc_out_x_dim <- tryCatch(expr=ncdf4::ncvar_get(nc_in, "x"), error = function(e){
#   ncdf4::ncvar_get(nc_in, xname)
# })
# nc_out_y_dim <- tryCatch(expr=ncdf4::ncvar_get(nc_in, "y"), error = function(e){
#   ncdf4::ncvar_get(nc_in, yname)
# })
nc_out_x_dim <- tryCatch(expr=ncdf4::ncvar_get(nc_in, xname), error = function(e){
  ncdf4::ncvar_get(nc_in, "x")
})
nc_out_y_dim <- tryCatch(expr=ncdf4::ncvar_get(nc_in, yname), error = function(e){
  ncdf4::ncvar_get(nc_in, "y")
})

# get CRS information from input file
proj4_args <- terra::crs(nc_rast, proj=TRUE)
epsg_code <- terra::crs(nc_rast, describe=TRUE)$code
crs_wkt <- terra::crs(nc_rast)
geoTransform <- as.character(paste(terra:::.geotransform(input_file), collapse=" "))

nc_att_geo_extent <- as.double(as.vector(terra::ext(nc_rast)))
nc_att_geo_res <- terra::res(nc_rast)
resolution <- terra::res(nc_rast)[1]
# convert to EPSG:4326
geo_extent <- as.double(as.vector(terra::project(terra::ext(nc_rast), from=paste("epsg:", epsg_code, sep=""), to="epsg:4326")))

# get history from input NetCDF file
if(ncdf4::ncatt_get(nc_in,0, "history")$hasatt){
  nc_in_hist <- as.character(ncdf4::ncatt_get(nc_in,0, "history")$value)
} else {
  nc_in_hist <- ""
}

# get variable information
if(ncdf4::ncatt_get(nc_in, var_name, "long_name")$hasatt){
  var_longname <- as.character(ncdf4::ncatt_get(nc_in, var_name, "long_name")$value)
} else {
  var_longname <- var_name
}
if(ncdf4::ncatt_get(nc_in, var_name, "units")$hasatt){
  var_units <- as.character(ncdf4::ncatt_get(nc_in, var_name, "units")$value)
} else {
  var_units <- "dl"
  # var_units <- "m2 m-2"
}

# get offset and scale_factor for the processed variable
if(nc_in$var[[which(names(nc_in$var) == var_name)]]$hasAddOffset){
  var_offset <- as.double(nc_in$var[[which(names(nc_in$var) == var_name)]]$addOffset)
} else {
  if(!is.null(ncatt_get(nc_in, 0)$var_add_offset)){
    var_offset <- as.double(ncatt_get(nc_in, 0)$var_add_offset)
  } else {
    var_offset <- as.double(0.0)
  }
}
if(nc_in$var[[which(names(nc_in$var) == var_name)]]$hasScaleFact){
  var_scale_factor <- as.double(nc_in$var[[which(names(nc_in$var) == var_name)]]$scaleFact)
} else {
  if(!is.null(ncatt_get(nc_in, 0)$var_scale_factor)){
    var_scale_factor <- as.double(ncatt_get(nc_in, 0)$var_scale_factor)
  } else {
    var_scale_factor <- as.double(1.0)
  }
}

if(var_offset == 0 && var_scale_factor == 1){
  var_scaled <- FALSE
} else {
  var_scaled <- TRUE
}
# set NetCDF missing value and precision
if(var_scaled){
  # nc_missval <- -32768
  nc_missval <- nc_var_atts$`_FillValue`
  nc_prec <- "short"
} else {
  nc_missval <- NULL
  nc_prec <- "double"
}

# define specific output variables offset and scale factor
sti_var_scale_factor <- as.double(var_scale_factor * 100)
sti_var_offset <- as.double(var_offset + 0)
mps_var_scale_factor <- as.double(var_scale_factor / 100)
mps_var_offset <- as.double(var_offset + 0)

# ##########################
# import mask dataset

# check mask argument
if(is.null(nc_mask_arg)){
  masking <- FALSE
} else {
  masking <- TRUE
  if(nc_mask_arg != "internal"){
    mask_file <- normalizePath(path=nc_mask_arg, winslash = "/", mustWork = TRUE)
  }
}
if(masking){
  if(nc_mask_arg == "internal"){
    if(length(which(names(nc_in$var) == "mask")) != 1){
      if(force){
        warning(paste("Input NetCDF file '", input_file, "' does not contain a variable named 'mask'. Going on however skipping without masking out pixels ...", sep=""))
      } else {
        ncdf4::nc_close(nc_in)
        stop(paste("Input NetCDF file '", input_file, "' does not contain a variable named 'mask'. Remove '-m' argument or set another input file.", sep=""))
      }
    }
  } else {
    mask_rast_check <- tryCatch({terra::rast(mask_file, quiet=TRUE); TRUE}, error = function(e) FALSE)
    if(mask_rast_check){
      rast_mask <- terra::rast(mask_file)
      mask_check <- terra::compareGeom(nc_rast, rast_mask, lyrs=FALSE, crs=TRUE, rowcol=TRUE, res=TRUE, stopOnError=FALSE)
      if(!mask_check){
        if(force){
          warning("Geometry of input mask is not consistent with input NetCDF file. Going on however without masking out pixels ...")
          masking <- FALSE
        } else {
          ncdf4::nc_close(nc_in)
          stop("Geometry of input mask is not consistent with input NetCDF file. Remove '-m' argument or set another input file.")
        }
      }
    } else {
      vect_check <- tryCatch({sf::st_read(mask_file, stringsAsFactors=FALSE, quiet=TRUE); TRUE}, error = function(e) FALSE)
      if(vect_check){
        # import vector mask
        vect_mask <- sf::st_read(mask_file, stringsAsFactors=FALSE, quiet=TRUE)
        vect_mask_epsg <- sf::st_crs(vect_mask)$epsg
        # reproject
        if(vect_mask_epsg != epsg_code){
          vect_mask <- sf::st_transform(vect_mask, crs=sf::st_crs(as.integer(epsg_code)))
        }
        # convert to raster
        rast_mask <- terra::rasterize(vect_mask, nc_rast[[1]], background=0)
        time(rast_mask) <- NULL
      } else {
        if(force){
          warning("Input mask is not a raster file nor a vector file. Going on however without masking out pixels ...")
          masking <- FALSE
        } else {
          ncdf4::nc_close(nc_in)
          stop("Input mask is not a raster file nor a vector file. Remove '-m' argument or set another input file.")
        }
      }
    }
  }
}

# import mask file
if(masking){
  # import mask data
  if(nc_mask_arg == "internal"){
    r_mask <- as.matrix(ncdf4::ncvar_get(nc=nc_in, varid="mask", collapse_degen = FALSE))
  } else {
    r_mask <- as.matrix(rast_mask, wide=TRUE)
  }
  r_mask[which(is.na(r_mask))] <- 0
  # check if mask contains pixels to be processed
  if(max(as.vector(r_mask), na.rm=TRUE) == 0){
    ncdf4::nc_close(nc_in)
    stop("Input mask does not contain pixel to process in the input raster file extent. Remove '-m' argument or set another input file.")
  }
  invisible(gc())
} else {
  # create full raster mask if mask file is not provided
  r_mask <- NULL
}

############################################################
# Setup output NetCDF file
############################################################

# check prefix argument
if(is.null(out_fname_prefix)){
  out_fname_prefix <- ""
} else {
  out_fname_prefix <- paste(gsub(pattern=" ", replacement = "_", x = out_fname_prefix), "_", sep="")
}
# get input NetCDF file attributes
if(ncatt_get(nc_in, 0, attname="S2_tile")$hasatt){
  tile <- as.character(ncatt_get(nc_in, 0, attname="S2_tile")$value)
  out_fname_prefix <- paste(out_fname_prefix, tile, "_", sep="")
}
# ### TODO get attributes to establish if input is Sentinel-2

# set product ID
name_id <- paste(out_fname_prefix, year, "_", "S2", "_", "L4", "_", resolution, "m", "_", "PM", "_", "s", max_season, "_v2.1", sep="")

# set output file name
if(is.null(recycling)){
  output_file <- normalizePath(path=paste(output_path, "/", out_fname_prefix, year, "_", "S2", "_", "L4", "_", resolution, "m", "_", var_name, "_", "PM", "_", "s", max_season, ".nc", sep=""), winslash = "/", mustWork = FALSE)
  if(file.exists(output_file)){
    stop("Output file name already exists in the selected output folder. Provide another output path using argument '-o'.")
  }
} else {
  output_file <- recycling
}

# create output NetCDF lock file
output_file_lck <- normalizePath(paste(tmp_folder, "/", basename(output_file), "_output.lck", sep=""), winslash = "/", mustWork = FALSE)

# create input NetCDF lock file
input_file_lck <- normalizePath(paste(tmp_folder, "/", basename(output_file), "_input.lck", sep=""), winslash = "/", mustWork = FALSE)

# export temporary mask
if(masking & !is.null(r_mask)){
  r_mask <- terra::rast(extent=terra::ext(nc_rast), resolution=terra::res(nc_rast), crs=terra::crs(nc_rast), vals=r_mask)
  writeRaster(r_mask, filename = normalizePath(paste(tmp_folder, "/", basename(output_file), "_rastMask.tif", sep=""), winslash = "/", mustWork = FALSE), filetype="GTiff", datatype="INT1U", gdal=c("COMPRESS=LZW"))
  r_mask <- normalizePath(paste(tmp_folder, "/", basename(output_file), "_rastMask.tif", sep=""), winslash = "/", mustWork = FALSE)
  invisible(gc())
  # create input mask lock file
  input_mask_lck <- normalizePath(paste(tmp_folder, "/", basename(output_file), "_inputMask.lck", sep=""), winslash = "/", mustWork = FALSE)
}

# ###########################
# Create output NetCDF

# setup output NetCDF file
nc_out_xdim <- ncdf4::ncdim_def(name=xname, longname="easting", units="meters", vals=nc_out_x_dim, unlim=FALSE)
nc_out_ydim <- ncdf4::ncdim_def(name=yname, longname="northing", units="meters", vals=nc_out_y_dim, unlim=FALSE)

if(chunksize > length(nc_out_x_dim) | chunksize == 0){
  nc_out_x_chunksize <- length(nc_out_x_dim)
} else {
  if(chunksize < 128){
    nc_out_x_chunksize <- 128
  } else {
    nc_out_x_chunksize <- chunksize
  }
}
if(chunksize > length(nc_out_y_dim) | chunksize == 0){
  nc_out_y_chunksize <- length(nc_out_y_dim)
} else {
  if(chunksize < 128){
    nc_out_y_chunksize <- 128
  } else {
    nc_out_y_chunksize <- chunksize
  }
}
nc_out_chunksize <- c(nc_out_x_chunksize,nc_out_y_chunksize)
nc_out_xy_size <- as.integer(length(nc_out_x_dim) * length(nc_out_y_dim))

utm_number <- 0
if(c("crs") %in% names(nc_in$var)){
  if(ncdf4::ncatt_get(nc_in, "crs", "utm_zone_number")$hasatt){
    utm_number <- as.integer(ncdf4::ncatt_get(nc_in, "crs", "utm_zone_number")$value)
  }
} else {
  if((epsg_num > 32600 & epsg_num <= 32660) | (epsg_num > 32700 & epsg_num <= 32760)){
    utm_number <- substr(epsg_num, 4, 5)
  } else {
    utm_number <- 0
  }
}

nc_miss_value <- as.integer(-32768)

# ###########################
# create variable list
nc_var_crs <- ncdf4::ncvar_def(name="crs", longname="CRS definition", units="dl", prec="short", dim=list())

# nc_var_veg_mask <- ncdf4::ncvar_def(name="veg_mask", longname="Vegetation Mask", units="dl", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
# nc_var_vi_avg <- ncdf4::ncvar_def(name=paste(var_name, "_", "avg", sep=""), longname=paste(var_longname, " ", "average", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
# nc_var_vi_std <- ncdf4::ncvar_def(name=paste(var_name, "_", "std", sep=""), longname=paste(var_longname, " ", "standard deviation", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
# nc_var_vi_min <- ncdf4::ncvar_def(name=paste(var_name, "_", "min", sep=""), longname=paste(var_longname, " ", "minimum", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
# nc_var_vi_max <- ncdf4::ncvar_def(name=paste(var_name, "_", "max", sep=""), longname=paste(var_longname, " ", "maximum", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
# nc_var_vi_delta <- ncdf4::ncvar_def(name=paste(var_name, "_", "delta", sep=""), longname=paste(var_longname, " ", "delta", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
# nc_var_vi_djf_avg <- ncdf4::ncvar_def(name=paste(var_name, "_", "djf_avg", sep=""), longname=paste(var_longname, " ", "DJF average", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
# nc_var_vi_djf_min <- ncdf4::ncvar_def(name=paste(var_name, "_", "djf_min", sep=""), longname=paste(var_longname, " ", "DJF minimum", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
# nc_var_vi_djf_max <- ncdf4::ncvar_def(name=paste(var_name, "_", "djf_max", sep=""), longname=paste(var_longname, " ", "DJF maximum", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
# nc_var_vi_jja_avg <- ncdf4::ncvar_def(name=paste(var_name, "_", "jja_avg", sep=""), longname=paste(var_longname, " ", "JJA average", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
# nc_var_vi_jja_max <- ncdf4::ncvar_def(name=paste(var_name, "_", "jja_max", sep=""), longname=paste(var_longname, " ", "JJA maximum", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)

nc_var_list <- list()
for(s in 1:max_season){
  var_prefix <- paste("VC", s, "_", sep="")
  
  nc_var_season <- ncdf4::ncvar_def(name=paste(var_prefix, "season", sep=""), longname="Seasons number", units="dl", missval=-127, prec="byte", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_pm_mask <- ncdf4::ncvar_def(name=paste(var_prefix, "PM_mask", sep=""), longname="Phenological Metrics mask", units="dl", missval=-127, prec="byte", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  # nc_var_seasons <- ncdf4::ncvar_def(name=paste(var_prefix, "seasons", sep=""), longname="Number of seasons", units="dl", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_qf <- ncdf4::ncvar_def(name=paste(var_prefix, "QF", sep=""), longname="Quality Flags", units="dl", missval=NULL, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  
  nc_var_soc <- ncdf4::ncvar_def(name=paste(var_prefix, "SoC_date", sep=""), longname="Start of Cycle date", units=c("days since 1970-01-01 00:00:00"), missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_soc_doy <- ncdf4::ncvar_def(name=paste(var_prefix, "SoC_doy", sep=""), longname="Start of Cycle (DOY)", units="dl", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_soc_value <- ncdf4::ncvar_def(name=paste(var_prefix, "SoC_value", sep=""), longname=paste("Start of Cycle ", var_name," value", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_eoc <- ncdf4::ncvar_def(name=paste(var_prefix, "EoC_date", sep=""), longname="End of Cycle date", units=c("days since 1970-01-01 00:00:00"), missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_eoc_doy <- ncdf4::ncvar_def(name=paste(var_prefix, "EoC_doy", sep=""), longname="End of Cycle (DOY)", units="dl", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_eoc_value <- ncdf4::ncvar_def(name=paste(var_prefix, "EoC_value", sep=""), longname=paste("End of Cycle ", var_name," value", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  
  nc_var_sos <- ncdf4::ncvar_def(name=paste(var_prefix, "SoS_date", sep=""), longname="Start of Season date", units=c("days since 1970-01-01 00:00:00"), missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_sos_doy <- ncdf4::ncvar_def(name=paste(var_prefix, "SoS_doy", sep=""), longname="Start of Season (DOY)", units="dl", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_sos_value <- ncdf4::ncvar_def(name=paste(var_prefix, "SoS_value", sep=""), longname=paste("Start of Season ", var_name," value", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_sgs <- ncdf4::ncvar_def(name=paste(var_prefix, "SGS_date", sep=""), longname="Start of Growing Season date", units=c("days since 1970-01-01 00:00:00"), missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_sgs_doy <- ncdf4::ncvar_def(name=paste(var_prefix, "SGS_doy", sep=""), longname="Start of Growing Season (DOY)", units="dl", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_sgs_value <- ncdf4::ncvar_def(name=paste(var_prefix, "SGS_value", sep=""), longname=paste("Start of Growing Season ", var_name," value", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_greenup <- ncdf4::ncvar_def(name=paste(var_prefix, "greenup_date", sep=""), longname="Greenup date", units=c("days since 1970-01-01 00:00:00"), missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_greenup_doy <- ncdf4::ncvar_def(name=paste(var_prefix, "greenup_doy", sep=""), longname="Greenup (DOY)", units="dl", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_greenup_value <- ncdf4::ncvar_def(name=paste(var_prefix, "greenup_value", sep=""), longname=paste("Greenup ", var_name, " value", sep=""), units=var_units, missval=nc_missval, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_greenup_rate <- ncdf4::ncvar_def(name=paste(var_prefix, "greenup_rate", sep=""), longname=paste(var_name, " Greenup rate", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_smp <- ncdf4::ncvar_def(name=paste(var_prefix, "SMP_date", sep=""), longname="Start of Maturity Plateau date", units=c("days since 1970-01-01 00:00:00"), missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_smp_doy <- ncdf4::ncvar_def(name=paste(var_prefix, "SMP_doy", sep=""), longname="Start of Maturity Plateau (DOY)", units="dl", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_smp_value <- ncdf4::ncvar_def(name=paste(var_prefix, "SMP_value", sep=""), longname=paste("Start of Maturity Plateau ", var_name," value", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_pos <- ncdf4::ncvar_def(name=paste(var_prefix, "PoS_date", sep=""), longname="Peak of Season date", units=c("days since 1970-01-01 00:00:00"), missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_pos_doy <- ncdf4::ncvar_def(name=paste(var_prefix, "PoS_doy", sep=""), longname="Peak of Season (DOY)", units="dl", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_pos_value <- ncdf4::ncvar_def(name=paste(var_prefix, "PoS_value", sep=""), longname=paste("Peak of Season ", var_name, " value", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_egs <- ncdf4::ncvar_def(name=paste(var_prefix, "EGS_date", sep=""), longname="End of Growing Season date", units=c("days since 1970-01-01 00:00:00"), missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_egs_doy <- ncdf4::ncvar_def(name=paste(var_prefix, "EGS_doy", sep=""), longname="End of Growing Season (DOY)", units="dl", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_egs_value <- ncdf4::ncvar_def(name=paste(var_prefix, "EGS_value", sep=""), longname=paste("End of Growing Season ", var_name," value", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_senescence <- ncdf4::ncvar_def(name=paste(var_prefix, "senescence_date", sep=""), longname="Senescence date", units=c("days since 1970-01-01 00:00:00"), missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_senescence_doy <- ncdf4::ncvar_def(name=paste(var_prefix, "senescence_doy", sep=""), longname="Senescence (DOY)", units="dl", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_senescence_value <- ncdf4::ncvar_def(name=paste(var_prefix, "senescence_value", sep=""), longname=paste("Senescence ", var_name," value", sep=""), units=var_units, missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_senescence_rate <- ncdf4::ncvar_def(name=paste(var_prefix, "senescence_rate", sep=""), longname=paste(var_name, " Senescence rate", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_eos <- ncdf4::ncvar_def(name=paste(var_prefix, "EoS_date", sep=""), longname="End of Season date", units=c("days since 1970-01-01 00:00:00"), missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_eos_doy <- ncdf4::ncvar_def(name=paste(var_prefix, "EoS_doy", sep=""), longname="End of Season (DOY)", units="dl", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_eos_value <- ncdf4::ncvar_def(name=paste(var_prefix, "EoS_value", sep=""), longname=paste("End of Season ", var_name," value", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_amplitude <- ncdf4::ncvar_def(name=paste(var_prefix, "seasonal_amplitude", sep=""), longname=paste(var_name, " ", "seasonal amplitude", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_dos <- ncdf4::ncvar_def(name=paste(var_prefix, "DoS", sep=""), longname="Duration of Season", units="day", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_lmp <- ncdf4::ncvar_def(name=paste(var_prefix, "LMP", sep=""), longname="Length of Maturity Plateau", units="day", missval=nc_miss_value, prec="short", dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_plateau_slope <- ncdf4::ncvar_def(name=paste(var_prefix, "maturity_plateau_slope", sep=""), longname=paste("Maturity Plateau slope ", var_name," value", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_sti <- ncdf4::ncvar_def(name=paste(var_prefix, "STI", sep=""), longname=paste("Seasonal Time Integrated ", var_name," value", sep=""), units=var_units, missval=nc_miss_value, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_nsti <- ncdf4::ncvar_def(name=paste(var_prefix, "NSTI", sep=""), longname=paste("Net Seasonal Time Integrated ", var_name," value", sep=""), units=var_units, missval=nc_miss_value, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_smv <- ncdf4::ncvar_def(name=paste(var_prefix, "SMV", sep=""), longname=paste("Seasonal Mean Value ", var_name," value", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_pos_diff <- ncdf4::ncvar_def(name=paste(var_prefix, "PoS_AG_difference", sep=""), longname=paste("Peak of Season Asymmetric Gaussian difference ", var_name, " value", sep=""), units=var_units, missval=nc_missval, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  nc_var_sti_diff <- ncdf4::ncvar_def(name=paste(var_prefix, "STI_AG_difference", sep=""), longname=paste("Seasonal Time Integrated Asymmetric Gaussian difference ", var_name," value", sep=""), units=var_units, missval=nc_miss_value, prec=nc_prec, dim=list(nc_out_xdim, nc_out_ydim), chunksizes=nc_out_chunksize, compression=nc_compression)
  
  nc_var_list <- c(nc_var_list, list(nc_var_season, nc_var_soc, nc_var_soc_doy, nc_var_soc_value, nc_var_eoc, nc_var_eoc_doy, nc_var_eoc_value, nc_var_pm_mask, 
                                     nc_var_sos, nc_var_sos_doy, nc_var_sos_value, nc_var_sgs, nc_var_sgs_doy, nc_var_sgs_value, nc_var_greenup, nc_var_greenup_doy, nc_var_greenup_value, nc_var_greenup_rate,
                                                                          nc_var_smp, nc_var_smp_doy, nc_var_smp_value, nc_var_pos, nc_var_pos_doy, nc_var_pos_value, nc_var_egs, nc_var_egs_doy, nc_var_egs_value, 
                                                                          nc_var_senescence, nc_var_senescence_doy, nc_var_senescence_value,  nc_var_senescence_rate, nc_var_eos, nc_var_eos_doy, nc_var_eos_value, 
                                                                          nc_var_amplitude, nc_var_dos, nc_var_lmp, nc_var_plateau_slope, nc_var_sti, nc_var_nsti, nc_var_smv, nc_var_pos_diff, nc_var_sti_diff, nc_var_qf))
}

# generate list of variables
# nc_var_list <- list(nc_var_veg_mask, nc_var_vi_avg, nc_var_vi_std, nc_var_vi_min, nc_var_vi_max, nc_var_vi_delta, nc_var_vi_djf_avg, nc_var_vi_djf_min, nc_var_vi_djf_max, nc_var_vi_jja_avg, nc_var_vi_jja_max)
# nc_var_list <- rep(list(nc_var_season, nc_var_pm_mask, nc_var_sos, nc_var_sos_doy, nc_var_sos_value, 
#                     nc_var_sgs, nc_var_sgs_doy, nc_var_sgs_value, nc_var_greenup, nc_var_greenup_doy, nc_var_greenup_value, nc_var_greenup_rate,
#                     nc_var_smp, nc_var_smp_doy, nc_var_smp_value, nc_var_pos, nc_var_pos_doy, nc_var_pos_value, nc_var_egs, nc_var_egs_doy, nc_var_egs_value, 
#                     nc_var_senescence, nc_var_senescence_doy, nc_var_senescence_value,  nc_var_senescence_rate, nc_var_eos, nc_var_eos_doy, nc_var_eos_value, 
#                     nc_var_cos, nc_var_cos_doy, nc_var_cos_value, nc_var_amplitude, nc_var_dos, nc_var_lmp, nc_var_plateau_slope, nc_var_baseline, 
#                     nc_var_sti, nc_var_nsti, nc_var_smv, nc_var_pos_diff, nc_var_sti_diff, nc_var_qf), max_season)
# nc_var_list[[length(nc_var_list)+1]] <- list(nc_var_crs)

# get variable names
nc_var_names <- character()
for(j in 1:44){
  nc_var_names <- c(nc_var_names, gsub(pattern = "VC1_", replacement = "", nc_var_list[[j]]$name))
}

# add CRS variable
nc_var_list_full <- c(nc_var_list, list(nc_var_crs))

# create NetCDF file
if(is.null(recycling)){
  ncout <- ncdf4::nc_create(output_file, vars=nc_var_list_full, force_v4=TRUE)
} else {
  # check if output recycling file was generated using EO4PM processor
  ncout <- ncdf4::nc_open(filename = output_file, write=FALSE, readunlim = FALSE)
  if(ncdf4::ncatt_get(ncout,0, "processor_id")$hasatt){
    processor_id_check <- as.character(ncdf4::ncatt_get(ncout, 0, "processor_id")$value)
  } else {
    stop("Output NetCDF file must be generated by processor 'EO4PM' to be recycled.")
  }
  if(ncdf4::ncatt_get(ncout, 0, "processor_id")$hasatt){
    processor_id_check <- as.character(ncdf4::ncatt_get(ncout, 0, "processor_id")$value)
    if(processor_id_check != "EO4PM"){
      stop("Output NetCDF file must be generated by processor 'EO4PM' to be recycled.")
    }
  } else {
    stop("Output NetCDF file must be generated by processor 'EO4PM' to be recycled.")
  }
  if(ncdf4::ncatt_get(ncout,0, "processor_version")$hasatt){
    processor_version_check <- as.double(ncdf4::ncatt_get(ncout, 0, "processor_version")$value)
    if(processor_version_check < 2.0){
      stop("Output NetCDF file must be generated by processor 'EO4PM' version >= 2.0.")
    }
  }
  if(ncdf4::ncatt_get(ncout, 0, "vegetation_cycles")$hasatt){
    vegetation_cycles_check <- as.integer(ncdf4::ncatt_get(ncout, 0, "vegetation_cycles")$value)
    if(vegetation_cycles_check != max_season){
      if(force){
        warning(paste("Force maximum number of seasons in a solar year to the one of output NetCDF file: '", max_season, "'.", sep=""))
        max_season <- as.integer(vegetation_cycles_check)
      } else {
        stop(paste("Maximum number of seasons in a solar year set using argument '-s' differs from the one in output NetCDF file: '", max_season, "'.", sep=""))
      }
    }
  } else {
    stop("Output NetCDF file must generated by processor 'EO4PM' does not have a global attribute named 'vegetation_cycles'.")
  }
  suppressWarnings(rm(list=c("processor_id_check","processor_version_check","vegetation_cycles_check")))
}

# ###########################
# add attributes

if(is.null(recycling)){
  # CRS to every variable
  nc_list_var <- names(ncout$var)
  for(variable in nc_list_var){
    ncdf4::ncatt_put(ncout,variable,"grid_mapping","crs")
  }
  
  # put attributes to 'crs' variable
  ncdf4::ncatt_put(ncout,"crs","crs_wkt", as.character(crs_wkt))
  ncdf4::ncatt_put(ncout,"crs","spatial_ref", as.character(crs_wkt))
  ncdf4::ncatt_put(ncout,"crs","proj4", as.character(proj4_args))
  ncdf4::ncatt_put(ncout,"crs","epsg_code", as.character(epsg_code))
  ncdf4::ncatt_put(ncout,"crs","geotransform", as.character(geoTransform))
  # ncdf4::ncatt_put(ncout,"crs","grid_mapping_name","transverse_mercator", prec="text")
  
  # ncdf4::ncatt_put(ncout,"crs","long_name","CRS definition", prec="text")
  
  # # GDAL standard
  # ### GDAL fix for UTM longitude coordinates shifted by 360 m: 'gdalinfo --config GDAL_NETCDF_CENTERLONG_180 NO' (ref: https://trac.osgeo.org/gdal/ticket/6759)
  # if(!is_chunk){
  #   ncdf4::ncatt_put(ncout,"crs","grid_mapping_name","transverse_mercator", prec="text")
  #   ncdf4::ncatt_put(ncout,"crs","longitude_of_central_meridiane",crs_LoCM,prec="float")
  #   ncdf4::ncatt_put(ncout,"crs","false_easting",crs_FE,prec="float")
  #   ncdf4::ncatt_put(ncout,"crs","false_northing",crs_FN,prec="float")
  #   ncdf4::ncatt_put(ncout,"crs","latitude_of_projection_origin",crs_LoPO,prec="float")
  #   ncdf4::ncatt_put(ncout,"crs","longitude_of_prime_meridian",0.0,prec="float")
  #   ncdf4::ncatt_put(ncout,"crs","scale_factor_at_central_meridian",crs_SFCM,prec="float")
  #   ncdf4::ncatt_put(ncout,"crs","grid_mapping_name","transverse_mercator", prec="text")
  #   ncdf4::ncatt_put(ncout,"crs","spatial_ref", crs_wkt)
  #   ncdf4::ncatt_put(ncout,"crs","crs_wkt", crs_wkt)
  #   ncdf4::ncatt_put(ncout,"crs","GeoTransform", geoTransform, prec="text")
  #   #ncdf4::ncatt_put(ncout,"crs","_CoordinateTransformType","Projection", prec="text")
  #   #ncdf4::ncatt_put(ncout,"crs","_CoordinateAxisTypes", "GeoX GeoY", prec="text")
  #   ncdf4::ncatt_put(ncout,variable,"proj4", as.character(proj4_args))
  #   # ncdf4::ncatt_put(ncout,variable,"proj4", as.character(crs_proj4))
  # }
  
  # add offset and scale factor for 'Quality Flags - QF' variable
  # for(s in 1:max_season){
  #   ncdf4::ncatt_put(ncout,paste("VC", s, "_", "QF", sep=""),"add_offset", as.integer(32768), prec="integer")
  #   ncdf4::ncatt_put(ncout,paste("VC", s, "_", "QF", sep=""),"scale_factor", as.integer(0), prec="integer")
  # }
  # add offset and scale factor for selected variables
  var_scalist <- c(5,8,12,15,18,19,22,25,28,31,32,35,36,42,43)-1
  nc_list_var_scale <- nc_var_names[var_scalist]
  if(var_scaled){
    for(s in 1:max_season){
      for(variable in nc_list_var_scale){
        ncdf4::ncatt_put(ncout,paste("VC", s, "_", variable, sep=""),"add_offset", as.double(var_offset), prec="float")
        ncdf4::ncatt_put(ncout,paste("VC", s, "_", variable, sep=""),"scale_factor", as.double(var_scale_factor), prec="float")
      }
    }
  }
  # add offset and scale factor for 'STI' variables
  var_stilist <- c(39,40,43)
  nc_list_var_scale <- nc_var_names[var_stilist]
  for(s in 1:max_season){
    for(variable in nc_list_var_scale){
      ncdf4::ncatt_put(ncout,paste("VC", s, "_", variable, sep=""),"add_offset", as.double(sti_var_offset), prec="float")
      ncdf4::ncatt_put(ncout,paste("VC", s, "_", variable, sep=""),"scale_factor", as.double(sti_var_scale_factor), prec="float")
    }
  }
  # add offset and scale factor for 'maturity_plateau_slope' variable
  for(s in 1:max_season){
    ncdf4::ncatt_put(ncout,paste("VC", s, "_", "maturity_plateau_slope", sep=""),"add_offset", as.double(mps_var_offset), prec="float")
    ncdf4::ncatt_put(ncout,paste("VC", s, "_", "maturity_plateau_slope", sep=""),"scale_factor", as.double(mps_var_scale_factor), prec="float")
  }
  
  # set global attributes to output NetCDF-4 file
  ncdf4::ncatt_put(ncout,0,"title",paste("Phenological Metrics from satellite time series", sep=""))
  ncdf4::ncatt_put(ncout,0,"product_version", 2.1, prec="float")
  ncdf4::ncatt_put(ncout,0,"Conventions","CF-1.10")
  ncdf4::ncatt_put(ncout,0,"geospatial_lon_min",as.double(geo_extent[1]), prec="float")
  ncdf4::ncatt_put(ncout,0,"geospatial_lon_max",as.double(geo_extent[2]), prec="float")
  ncdf4::ncatt_put(ncout,0,"geospatial_lat_min",as.double(geo_extent[3]), prec="float")
  ncdf4::ncatt_put(ncout,0,"geospatial_lat_max",as.double(geo_extent[4]), prec="float")
  ncdf4::ncatt_put(ncout,0,"geospatial_x_resolution",as.double(nc_att_geo_res[1]), prec="float")
  ncdf4::ncatt_put(ncout,0,"geospatial_y_resolution",as.double(nc_att_geo_res[2]), prec="float")
  ncdf4::ncatt_put(ncout,0,"geospatial_x_min",as.integer(nc_att_geo_extent[1]), prec="integer")
  ncdf4::ncatt_put(ncout,0,"geospatial_x_max",as.integer(nc_att_geo_extent[2]), prec="integer")
  ncdf4::ncatt_put(ncout,0,"geospatial_y_min",as.integer(nc_att_geo_extent[3]), prec="integer")
  ncdf4::ncatt_put(ncout,0,"geospatial_y_max",as.integer(nc_att_geo_extent[4]), prec="integer")
  ncdf4::ncatt_put(ncout,0,"geospatial_x_units","meters")
  ncdf4::ncatt_put(ncout,0,"geospatial_y_units","meters")
  ncdf4::ncatt_put(ncout,0,"spatial_resolution", paste(nc_att_geo_res[1], " meters", sep=""))
  ncdf4::ncatt_put(ncout,0,"time_coverage_start", nc_file_date_start)
  ncdf4::ncatt_put(ncout,0,"time_coverage_end", nc_file_date_end)
  ncdf4::ncatt_put(ncout,0,"start_date", nc_att_date_start)
  ncdf4::ncatt_put(ncout,0,"end_date", nc_att_date_end)
  ncdf4::ncatt_put(ncout,0,"id",name_id)
  ncdf4::ncatt_put(ncout,0,"processing_level","L4")
  # ncdf4::ncatt_put(ncout,0,"summary","This dataset contains phenological metrics estimated from Sentinel-2 satellite observations distributed in granules")
  ncdf4::ncatt_put(ncout,0,"keywords","EARTH SCIENCE > BIOSPHERE > VEGETATION > PLANT PHENOLOGY")
  ncdf4::ncatt_put(ncout,0,"GCMD_id","3f45aadf-ec7c-43a1-a008-b24ca139837a")
  ncdf4::ncatt_put(ncout,0,"keywords_vocabulary","GCMD Science Keywords, Version 23.4")
  # ncdf4::ncatt_put(ncout,0,"standard_name_vocabulary","CF-1.7, Version 58")
  # ncdf4::ncatt_put(ncout,0,"comment","Generated using EO4PM operator version 2.1 from atmospherically corrected Sentinel-2 L2A data using MAJA algorithm distributed by Theia in MUSCATE format.")
  ncdf4::ncatt_put(ncout,0,"license","free and open access")
  #ncdf4::ncatt_put(ncout,0,"terms_for_use","These data can be used freely for research purposes provided that the following source is acknowledged: ")
  ncdf4::ncatt_put(ncout,0,"disclaimer","This data is made available in the hope that will be useful, but WITHOUT ANY WARRANTY")
  ncdf4::ncatt_put(ncout,0,"processor_id","EO4PM")
  ncdf4::ncatt_put(ncout,0,"processor_name","phenological_metrics_estimation_processor")
  ncdf4::ncatt_put(ncout,0,"processor_version","2.1")
  ncdf4::ncatt_put(ncout,0,"estimation_method","Gu")
  ncdf4::ncatt_put(ncout,0,"calendar_type",calendar)
  # ncdf4::ncatt_put(ncout,0,"full_seasonal_djf",as.character(full_seasonal_djf_flag))
  # ncdf4::ncatt_put(ncout,0,"seasonal_djf_bounds",seasonal_djf_bounds_attr)
  ncdf4::ncatt_put(ncout,0,"phenological_year",year)
  ncdf4::ncatt_put(ncout,0,"vegetation_cycles",max_season, prec="short")
  ncdf4::ncatt_put(ncout,0,"input_filename",nc_file_name)
  ncdf4::ncatt_put(ncout,0,"processor_keywords","phenology")
  ncdf4::ncatt_put(ncout,0,"cdm_data_type","Grid")
  ncdf4::ncatt_put(ncout,0,"grid_type","GeoX GeoY")
  ncdf4::ncatt_put(ncout,0,"proj4_params",proj4_args)
  ncdf4::ncatt_put(ncout,0,"EPSG_code",epsg_code)
  
  # get extra attributes from input file
  if(ncdf4::ncatt_get(nc_in,0, "temporal_coregistration_applied")$hasatt){
    nc_att_coreg <- as.character(ncdf4::ncatt_get(nc_in,0, "temporal_coregistration_applied")$value)
  } else {
    nc_att_coreg <- as.character("Unknown")
  }
  ncdf4::ncatt_put(ncout,0,"temporal_coregistration_applied",nc_att_coreg)
  if(ncdf4::ncatt_get(nc_in,0, "smoothing_type")$hasatt){
    nc_att_smoothing <- as.character(ncdf4::ncatt_get(nc_in,0, "smoothing_type")$value)
  } else {
    nc_att_smoothing <- as.character("Unknown")
  }
  ncdf4::ncatt_put(ncout,0,"smoothing_type",nc_att_smoothing)
  if(ncdf4::ncatt_get(nc_in,0, "smoothed_variable")$hasatt){
    nc_att_smooth_var <- as.character(ncdf4::ncatt_get(nc_in,0, "smoothed_variable")$value)
  } else {
    nc_att_smooth_var <- as.character(var_name)
  }
  ncdf4::ncatt_put(ncout,0,"smoothed_variable",nc_att_smooth_var)
  
  # ### TODO add this lines  only if satellite is Sentinel-2
  # ncdf4::ncatt_put(ncout,0,"source","Copernicus Sentinel-2 data at L2A")
  # ncdf4::ncatt_put(ncout,0,"platform","Sentinel-2A, Sentinel-2B, Sentinel-2C")
  # ncdf4::ncatt_put(ncout,0,"sensor","MSI")
  # ncdf4::ncatt_put(ncout,0,"S2_tile",tile)
  # nc_acknowledgment <- paste("Input Sentinel-2 data available at no cost from ESA Copernicus Open Access Hub. Copernicus Sentinel-2 data processed at level 2A by CNES for THEIA Land data center. Data processing has been performed using ''. This dataset contains modified Copernicus Sentinel data (", format(Sys.Date(), "%Y"), ").", sep="")
  # ncdf4::ncatt_put(ncout,0,"acknowledgment", nc_acknowledgment)
  
  # ncdf4::ncatt_put(ncout,0,"naming_authority","")
  # ncdf4::ncatt_put(ncout,0,"institution","")
  # ncdf4::ncatt_put(ncout,0,"processing_centre","")
  # ncdf4::ncatt_put(ncout,0,"creator_name","")
  # ncdf4::ncatt_put(ncout,0,"creator_url","")
  # ncdf4::ncatt_put(ncout,0,"creator_email","")
  # ncdf4::ncatt_put(ncout,0,"project","")
  
  # # optional attributes
  # ncdf4::ncatt_put(ncout,0,"date_issued","")
  # ncdf4::ncatt_put(ncout,0,"date_modified","")
  # ncdf4::ncatt_put(ncout,0,"publisher_name","")
  # ncdf4::ncatt_put(ncout,0,"publisher_url","")
  # ncdf4::ncatt_put(ncout,0,"publisher_email","")
  # ncdf4::ncatt_put(ncout,0,"contributor_name","")
  # ncdf4::ncatt_put(ncout,0,"contributor_role","")
  
  # ncdf4::ncatt_put(ncout,0,"contact","")
  # ncdf4::ncatt_put(ncout,0,"date_issued","")
  # ncdf4::ncatt_put(ncout,0,"references","")
  
  # setup nc history
  nc_hist_past <- as.character(paste("EO4PM.R --input=", input_file, " --mask=", nc_mask_arg, 
                                     " --recycling=", ifelse(is.null(opt$recycling), paste("NULL", " --output=", output_file, sep=""), as.character(opt$recycling)),
                                     " --year=", year, " --variable=", var_name, " --calendar=", calendar, 
                                     " --start_date=", ifelse(is.null(opt$start_date), "NULL", as.character(opt$start_date)), 
                                     " --end_date=", ifelse(is.null(opt$end_date), "NULL", as.character(opt$end_date)), 
                                     " --season=", max_season, " --min_value=", ifelse(is.null(minValue), "NULL", as.character(minValue)), " --s_lag=", s_lag, " --maxExtendMonth=", maxExtendMonth,
                                     " --minPeakValue=", ypeak_min, " --minPeakDistance=", minpeakdistance, " --minSeasonLength=", length_min, 
                                     " --r_min=", r_min, " --r_max=", r_max, " --rtrough_max=", rtrough_max, 
                                     " --force=", force, " --chunksize=", chunksize, " --deflate=", nc_compression, " --cores=", cores, sep=""))
  
  # fill in QF variables with zeros
  for(s in 1:max_season){
    ncdf4::ncvar_put(ncout, varid=paste("VC", s, "_", "QF", sep=""), vals=rep(as.integer(0), nc_out_xy_size), verbose = FALSE)
  }
  
  ncdf4::nc_sync(ncout)
  ncdf4::nc_close(ncout)
}

# close input NetCDF file
ncdf4::nc_close(nc_in)

# remove NetCDF objects
suppressWarnings(rm(list=c("nc_in","ncout")))

# clean workspace
suppressWarnings(rm(list=c("calendar", "crs_wkt", "end_date", "end_pos", "epsg_code", "file.arg.name", "full_seasonal_djf_flag", 
                           "geo_extent", "geoTransform", "initial.options", "mask_check", "mask_file", "name_id", "nc_acknowledgment", 
                           "nc_att_coreg", "nc_att_date_end", "nc_att_date_start", "nc_att_geo_extent", "nc_att_geo_res", "nc_att_smooth_var", "nc_att_smoothing", 
                           "nc_file_date_end", "nc_file_date_start", "nc_file_name", "nc_mask_arg", "nc_missval", "nc_out_chunksize", "nc_out_time_dim", 
                           "nc_out_x_chunksize", "nc_out_xdim", "nc_out_xy_size", "nc_out_y_chunksize", "nc_out_ydim", "nc_prec", "nc_rast", 
                           "obs_dates", "output_path", "pheno_lag", "PROCESSOR_HOME", "nc_compression", "nc_var_list", "nc_var_list_full", "nc_var_atts", "nc_var_crs", 
                           "proj4_args", "rast_mask", "resolution", "s", "script.name", "seasonal_djf_bound1", "seasonal_djf_bound2", "seasonal_djf_bounds_attr", 
                           "start_date", "start_pos", "tile", "tm_int", "utm_number", "var_longname", "var_units", "variable", "out_fname_prefix", 
                           "nc_var_season", "nc_var_pm_mask", "nc_var_qf", "nc_var_sos", "nc_var_sos_doy", "nc_var_sos_value", "nc_var_sgs", "nc_var_sgs_doy", "nc_var_sgs_value", 
                           "nc_var_greenup", "nc_var_greenup_doy", "nc_var_greenup_value", "nc_var_greenup_rate", "nc_var_smp", "nc_var_smp_doy", "nc_var_smp_value", 
                           "nc_var_pos", "nc_var_pos_doy", "nc_var_pos_value", "nc_var_egs", "nc_var_egs_doy", "nc_var_egs_value", "var_dims", "var_scalist", "var_stilist", 
                           "nc_var_senescence", "nc_var_senescence_doy", "nc_var_senescence_value", "nc_var_senescence_rate", "nc_var_eos", "nc_var_eos_doy", "nc_var_eos_value", 
                           "nc_var_soc", "nc_var_soc_doy", "nc_var_soc_value", "nc_var_eoc", "nc_var_eoc_doy", "nc_var_eoc_value", "nc_var_amplitude", "nc_var_dos", "nc_var_lmp", "nc_var_plateau_slope", 
                           "nc_var_sti", "nc_var_nsti", "nc_var_smv", "nc_var_pos_diff", "nc_var_sti_diff", "xname", "yname", "var_prefix", "var_dims_names", "readRast_opts", 
                           "nc_var_name", "nc_list_var", "nc_list_var_scale", "nc_miss_value", "d", "j")))
invisible(gc())

############################################################
# Perform Phenological Metrics estimation
############################################################

# define minimum value to cut cycles
if(is.null(minValue)){
  minValue_ylu <- NA
} else {
  minValue_ylu <- minValue
}

# use chunks (to optimize RAM usage)
nr <- length(nc_out_x_dim)
tnr <- ceiling(nr/chunksize)
tnrv <- as.vector(1:tnr)
nc <- length(nc_out_y_dim)
tnc <- ceiling(nc/chunksize)
tncv <- as.vector(1:tnc)

rlg <- expand.grid(tnrv,tncv)
names(rlg) <- c("x","y")
rl <- nrow(rlg)
r_list <- 1:rl
# rl <- length(tnrv)*length(tncv)

suppressWarnings(rm(list=c("nc_out_x_dim", "nc_out_y_dim", "minValue")))

# subset list of chunks based on those with pixels to be processed according to input mask
if(masking & !is.null(r_mask)){
  # read mask data
  full_mask <- suppressWarnings(terra::rast(r_mask))
  r_mask_count <- data.frame("CHUNK"=as.integer(r_list), "COUNT"=rep(as.integer(NA), length(r_list)))
  for(l in r_list){
    # define chunk to be imported
    r <- rlg[l,1]
    s <- rlg[l,2]
    rs <- as.integer(((r-1)*chunksize)+1)
    re <- as.integer(r*chunksize)
    if(r == tnr){
      re <- nr
    }
    cs <- as.integer(((s-1)*chunksize)+1)
    ce <- as.integer(s*chunksize)
    if(s == tnc){
      ce <- nc
    }
    rs_mask <- terra::crop(full_mask, terra::ext(c(min(c(terra::xFromCol(full_mask, rs),terra::xFromCol(full_mask, re)))-(terra::res(full_mask)[1]/2),max(c(terra::xFromCol(full_mask, rs),terra::xFromCol(full_mask, re)))+(terra::res(full_mask)[1]/2),min(c(terra::yFromRow(full_mask, cs),terra::yFromRow(full_mask, ce)))-(terra::res(full_mask)[2]/2),max(c(terra::yFromRow(full_mask, cs),terra::yFromRow(full_mask, ce)))+(terra::res(full_mask)[2]/2))))
    # count number of pixel to be processed in each chunk
    n_pix_count <- length(which(as.vector(rs_mask[]) == 1))
    r_mask_count$COUNT[which(r_mask_count$CHUNK == l)] <- n_pix_count
  }
  
  r_mask_count <- r_mask_count[which(r_mask_count$COUNT > 0),]
  
  if(nrow(r_mask_count) == 0){
    stop("No pixels to be processed in input mask.")
  }
  
  # reorder mask_count
  r_mask_order <- as.vector(rep(NA, (ceiling(nrow(r_mask_count)/cores) * cores)))
  r_mask_order[1:nrow(r_mask_count)] <- r_mask_count$CHUNK[order(r_mask_count$COUNT, decreasing = TRUE)]
  r_mask_order <- matrix(data=r_mask_order, ncol = cores)
  r_list <- as.vector(na.omit(as.vector(t(r_mask_order))))
  
  if(nrow(r_mask_count) < length(r_list)){
    if(verbose){
      message(paste("Masked chunks to be processed is: ", nrow(r_mask_count), "/", length(r_list), sep=""))
    }
  }
  suppressWarnings(rm(list=c("full_mask", "rs_mask", "r_mask_count","r_mask_order")))
  invisible(gc(verbose=FALSE))
}

# set doFuture export list
export_list <- c("rl","rlg","chunksize","tnr","tnc","nr","nc","input_file","input_file_lck","output_file","output_file_lck","r_mask","input_mask_lck", 
                 "var_name", "force", "nc_time_start","nc_time_count", "phenodates", "phenoGetEO4PM", "smooth_none", 
                 "year", "time_set_extract", "var_scaled", "var_scale_factor", "var_offset", "std_var_offset", "std_var_scale_factor", "sti_var_offset", "sti_var_scale_factor", 
                 "nc_var_names", "masking", "progress", "prg", "minValue_ylu", "max_season", "c_win", "s_lag", "mps_var_offset", "mps_var_scale_factor", 
                 # "compute_tstats", "tstatsEO4PM", 
                 # "vi_dataset", "r_list", 
                 "nptperyear", "maxExtendMonth", "ypeak_min", "minpeakdistance", "length_min", "r_min", "r_max", "rtrough_max")

analysis_fun <- function(numlist){
  prg <- progressr::progressor(along=numlist)
  foreach::foreach(l=numlist,
                   .combine=c,
                   .inorder=FALSE,
                   .errorhandling = "pass",
                   .options.future = list(globals = export_list, scheduling = FALSE, packages=c("terra","ncdf4","phenofit","data.table","progressr"), seed = TRUE),
                   # .options.future = list(globals = export_list, scheduling = 1.0, packages=c("terra","ncdf4","phenofit","progressr"), seed = TRUE),
                   .verbose=FALSE) %dofuture% {
                     
                     # ### for debug
                     # for(l in r_list){
                     # message(paste("Chunk:", l, sep=""))
                     
                     chunk_counter <- 0
                     
                     # blocco OpenMP e Algebra Lineare tramite variabili d'ambiente
                     Sys.setenv(OMP_THREAD_LIMIT = 1)
                     Sys.setenv(OPENBLAS_NUM_THREADS = 1)
                     Sys.setenv(MKL_NUM_THREADS = 1)
                     Sys.setenv(MATHLIB_NUM_THREADS = 1)
                     
                     # Sys.setenv(R_DATATABLE_NUM_THREADS = 1)
                     suppressPackageStartupMessages(require(data.table, quietly=TRUE))
                     data.table::setDTthreads(1)
                     
                     # define chunk to be imported
                     r <- rlg[l,1]
                     s <- rlg[l,2]
                     rs <- as.integer(((r-1)*chunksize)+1)
                     re <- as.integer(r*chunksize)
                     if(r == tnr){
                       re <- nr
                     }
                     chunks_nrow <- re -rs +1
                     cs <- as.integer(((s-1)*chunksize)+1)
                     ce <- as.integer(s*chunksize)
                     if(s == tnc){
                       ce <- nc
                     }
                     chunks_ncol <- ce - cs +1
                     
                     # ### for debug
                     # message(paste("ROW ", "start: ", rs, " count: ", chunks_nrow, " - COL ", "start: ", cs, " count: ", chunks_ncol, sep=""))
                     
                     # check input mask
                     mskip <- FALSE
                     m_mask <- NULL
                     if(masking & !is.null(r_mask)){
                       lck_mask <- TRUE
                       while(lck_mask == TRUE){
                         if(!file.exists(input_mask_lck)){
                           # lock input file
                           # lck_in <- filelock::lock(input_file, exclusive = FALSE)
                           write.csv(as.character(""), file = input_mask_lck)
                           # read mask data
                           r_mask <- suppressWarnings(terra::rast(r_mask))
                           rs_mask <- terra::crop(r_mask, terra::ext(c(min(c(terra::xFromCol(r_mask, rs),terra::xFromCol(r_mask, re)))-(terra::res(r_mask)[1]/2),max(c(terra::xFromCol(r_mask, rs),terra::xFromCol(r_mask, re)))+(terra::res(r_mask)[1]/2),min(c(terra::yFromRow(r_mask, cs),terra::yFromRow(r_mask, ce)))-(terra::res(r_mask)[2]/2),max(c(terra::yFromRow(r_mask, cs),terra::yFromRow(r_mask, ce)))+(terra::res(r_mask)[2]/2))))
                           m_mask <- as.matrix(rs_mask, wide=TRUE)
                           suppressWarnings(rm(list=c("r_mask","rs_mask")))
                           if(max(as.vector(m_mask), na.rm=TRUE) == 0){
                             mskip <- TRUE
                           }
                           lck_rm <- file.remove(input_mask_lck)
                           lck_mask <- FALSE
                         } else {
                           lck_mask <- TRUE
                           Sys.sleep(5)
                         }
                       }
                     }
                     
                     # skip if mask is void
                     if(!mskip){
                       
                       chunk_counter <- 1
                       
                       # ###########################################################
                       # Import data
                       # ###########################################################
                       
                       lck_in <- TRUE
                       
                       while(lck_in == TRUE){
                         if(!file.exists(input_file_lck)){
                           # lock input file
                           # lck_in <- filelock::lock(input_file, exclusive = FALSE)
                           write.csv(as.character(""), file = input_file_lck)
                           
                           # read data
                           nc_input <- ncdf4::nc_open(filename = input_file, write = FALSE, readunlim = FALSE)
                           nc_data <- ncdf4::ncvar_get(nc=nc_input, varid=var_name, start=c(rs,cs,nc_time_start), count=c(chunks_nrow,chunks_ncol,nc_time_count), collapse_degen = FALSE)
                           ncdf4::nc_close(nc_input)
                           suppressWarnings(rm(list=c("nc_input")))
                           # invisible(gc(verbose=FALSE))
                           
                           # release file lock
                           # flock <- filelock::unlock(lck_in)
                           lck_rm <- file.remove(input_file_lck)
                           lck_in <- FALSE
                         } else {
                           lck_in <- TRUE
                           Sys.sleep(10)
                           
                         }
                       }
                       
                       nc_out_dims <- c(dim(nc_data)[1],dim(nc_data)[2])
                       vi_dataset <- apply(nc_data, 3, c)
                       n_rows <- nrow(vi_dataset)
                       suppressWarnings(rm(list=c("nc_data")))
                       invisible(gc(verbose=FALSE))
                       
                       # subset input dataset according to input mask
                       if(masking & !is.null(m_mask)){
                         # identify pixel to be processed based on input mask
                         g_pixel_list <- as.vector(which(t(m_mask) == 1))
                         vi_dataset <- matrix(vi_dataset[g_pixel_list,], nrow=length(g_pixel_list))
                         suppressWarnings(rm(list=c("m_mask")))
                       } else {
                         g_pixel_list <- seq(1:n_rows)
                       }
                       
                       # ###########################################################
                       # Create pixel list to be processed
                       # ###########################################################
                       
                       # identify null pixel in temporal profiles
                       no_gap_obs <- apply(vi_dataset, MARGIN=1, FUN=function(x) sum(!is.na(x)))
                       
                       if(!is.na(ypeak_min)){
                         y_start <- which(phenodates == max(phenodates[1], as.integer(as.Date(paste(year, "-01-01", sep=""))), na.rm=TRUE))
                         y_end <- which(phenodates == min(phenodates[length(phenodates)], as.integer(as.Date(paste(year, "-12-31", sep=""))), na.rm=TRUE))
                         max_obs <- apply(vi_dataset[,(y_start:y_end)], MARGIN=1, FUN=function(x) max(x, na.rm=TRUE))
                         # skip time series with missing data or low peak values
                         pixel_list <- sort(unique(which(no_gap_obs == ncol(vi_dataset) & max_obs >= ypeak_min)))
                       } else {
                         # skip time series with missing data
                         pixel_list <- sort(unique(which(no_gap_obs == ncol(vi_dataset))))
                       }
                       
                       g_pixel_list <- g_pixel_list[pixel_list]
                       
                       # clean workspace
                       suppressWarnings(rm(list=c("no_gap_obs","max_obs","y_start","y_end")))
                       invisible(gc(verbose=FALSE))
                       
                       # get number of pixels to be processed
                       pl <- length(pixel_list)
                       
                       # continue only if there are pixels to be processed
                       if(pl > 0){
                         
                         # subset input dataset
                         vi_dataset <- matrix(vi_dataset[pixel_list,], nrow=length(pixel_list))
                         p_list <- c(1:length(pixel_list))
                         
                         # create empty object to store results
                         pks_dataset <- list()
                         for(s in 1:max_season){
                           pks_dataset[[s]] <- matrix(data=as.double(NA), nrow=(nc_out_dims[1]*nc_out_dims[2]), ncol=44)
                         }
                         invisible(gc(verbose=FALSE))
                         
                         # # define faster 'di2dt' function
                         # fast_di2dt <- function(di, t, ypred) {
                         #   res_t <- lapply(di, function(idx) t[idx])
                         #   res_y <- lapply(di, function(idx) ypred[idx])
                         #   names(res_y) <- paste0("y_", names(di))
                         #   dt <- data.table::setDT(c(res_t, res_y))
                         #   if (inherits(t, "Date") || inherits(t, "POSIXt")) {
                         #     dt[, `:=`(
                         #       len = as.integer(end - beg) + 1L,
                         #       year = data.table::year(peak) 
                         #     )]
                         #   } else {
                         #     dt[, `:=`(
                         #       len = as.numeric(end - beg) + 1,
                         #       year = NA_integer_
                         #     )]
                         #   }
                         #   return(dt)
                         # }
                         # 
                         # assignInNamespace("di2dt", fast_di2dt, ns = "phenofit")
                         
                         # ###########################################################
                         # Process single pixels
                         # ###########################################################
                         
                         for(p in p_list){
                           
                           # ### for debug
                           # message(p)
                           
                           # get pixel position in target matrix
                           q <- g_pixel_list[p]
                           
                           # ###########################################################
                           # Calculate temporal statistics
                           # ###########################################################
                           
                           tstats_check <- TRUE
                           # if(compute_tstats){
                           #   # ### TODO add percentiles?
                           #   y_tstats <- tstatsEO4PM(pixel = vi_dataset[p,], obs_dates = phenodates)
                           #   if(y_tstats[1] == 0){
                           #     tstats_check <- FALSE
                           #   }
                           # }
                           
                           # only continue processing if all conditions are met (no NA values)
                           if(tstats_check){
                             
                             # initialize QFs
                             qf0 <- 0
                             qf11 <- 0
                             qf12 <- 0
                             # qf13 <- 0
                             qf14 <- 0
                             
                             # ###########################################################
                             # Cut cycles
                             # ###########################################################
                             
                             if(max_season > 1){
                               # define weights to clean baseline-corrected curve artifacts
                               # w <- rep(as.double(1.0), length(phenodates))
                               # w[which(pixel[2:length(pixel)] == pixel[1:(length(pixel)-1)]) + 1] <- 0
                               
                               # convert to phenofit INPUT object
                               ylu <- quantile(as.double(vi_dataset[p,]), c(0.01,0.99))
                               ylu[1] <- pmax(ylu[1], min(0, minValue_ylu, na.rm=TRUE), na.rm=TRUE)
                               ylu[2] <- max(vi_dataset[p,], na.rm=TRUE)
                               pfo <- list(t = as.Date(phenodates),
                                           y = as.double(vi_dataset[p,]),
                                           # w = w,
                                           w = rep(as.double(1.0), length(phenodates)),
                                           ylu = ylu,
                                           nptperyear = 365,
                                           south = as.logical(FALSE)
                               )
                               
                               # detect cycles
                               # s <- suppressWarnings(phenofit::season_mov(pfo,
                               #                                        options = list(
                               #                                          # rFUN = "smooth_wWHIT",
                               #                                          rFUN = "smooth_none",
                               #                                          wFUN = "wTSM",
                               #                                          iters = 1, lambda = 10,
                               #                                          maxExtendMonth = maxExtendMonth,
                               #                                          minpeakdistance = minpeakdistance,
                               #                                          length_min = length_min,
                               #                                          ypeak_min = ypeak_min,
                               #                                          rtrough_max = rtrough_max,
                               #                                          r_min = r_min,
                               #                                          r_max = r_max,
                               #                                          rm.closed=TRUE,
                               #                                          adj.param = FALSE,
                               #                                          verbose = FALSE
                               #                                        )))
                               
                               s <- tryCatch(expr=suppressWarnings(phenofit::season(pfo,
                                                                      # rFUN = "smooth_wWHIT",
                                                                      rFUN = "smooth_none",
                                                                      wFUN = "wTSM",
                                                                      iters = 1, lambda = 10,
                                                                      minpeakdistance = minpeakdistance,
                                                                      length_min = length_min,
                                                                      ypeak_min = ypeak_min,
                                                                      rtrough_max = rtrough_max,
                                                                      MaxPeaksPerYear = max_season,
                                                                      MaxTroughsPerYear = as.integer(max_season + 1),
                                                                      r_min = r_min,
                                                                      r_max = r_max,
                                                                      rm.closed=TRUE,
                                                                      adj.param = FALSE,
                                                                      )), error = function(e){
                                                                        return(NULL)
                                                                      })
                               
                               if(is.null(s)){
                                 nseas <- 0
                                 qf11 <- 1
                               } else {
                                 sdt <- s[["dt"]]
                                 # plot_season(pfo,s)
                                 nseas <- max(0, nrow(sdt), na.rm=TRUE)
                                 if(nseas > 0){
                                   sdt <- sdt[which(as.integer(strftime(sdt$peak, format = "%Y")) == year),]
                                 }
                                 
                                 if(nseas > max_season){
                                   sdt <- sdt[sort(order(sdt$y_peak, decreasing = TRUE)[1:max_season]),]
                                   # nseas <- max_season
                                   nseas <- nrow(sdt)
                                   qf12 <- 1
                                 }
                               }
                               
                               suppressWarnings(rm(list=c("ylu","pfo","s","w")))
                             } else {
                               nseas <- 1
                             }
                             
                             # ###########################################################
                             # Estimate phenological metrics
                             # ###########################################################
                             
                             if(nseas > 0 & max_season != 1){
                               
                               for(s in 1:nseas){
                                 
                                 # set QFs
                                 qf0 <- 1
                                 
                                 s_beg <- as.integer(sdt$beg[s])
                                 s_pos <- which(phenodates %in% s_beg)
                                 s_end <- as.integer(sdt$end[s])
                                 e_pos <- which(phenodates %in% s_end)
                                 s_year <- as.integer(strftime(sdt$peak[s], format = "%Y"))
                                 
                                 # compute phenometrics
                                 pks_result <- phenoGetEO4PM(pixel=as.double(vi_dataset[p,max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)]), obs_dates = phenodates[max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)], year=s_year, s_lag=s_lag, c_win=c_win)
                                 
                                 # fit Asymmetric Gaussian curve if Gu Method fails
                                 if(force){
                                   if(pks_result[1] == 0){
                                     
                                     # fit Asymmetric Gaussian curve
                                     # pars <- phenofit::FitDL.Beck(y = vi_dataset[p,max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)], t = phenodates[max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)], tout = phenodates[max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)], w=ifelse(vi_dataset[p,max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)] > quantile(vi_dataset[p,max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)], probs=0.75, na.rm=TRUE), 1.0, 0.25), type = 1, iters = 2)
                                     # dl <- phenofit::doubleLog.Beck(par = as.double(phenofit::get_param(pars)), t=phenodates[max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)])
                                     dl <- phenofit::FitDL.AG(y = vi_dataset[p,max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)], t = phenodates[max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)], tout = phenodates[max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)], w=ifelse(vi_dataset[p,max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)] > quantile(vi_dataset[p,max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)], probs=0.75, na.rm=TRUE), 1.0, 0.25), type = 1, iters = 2)
                                     dl <- dl$zs$iter2
                                     
                                     # compute phenometrics
                                     pks_result <- phenoGetEO4PM(pixel=dl, obs_dates = phenodates[max(1, s_pos-s_lag):min(length(phenodates), e_pos+s_lag)], year=s_year, s_lag=s_lag, c_win=c_win)
                                     
                                     pks_result[35] <- max(dl[((s_lag+1):(length(dl-s_lag)))], na.rm=TRUE) - max(vi_dataset[p,(s_pos:e_pos)], na.rm=TRUE)
                                     pks_result[36] <- sum(dl[((s_lag+1):(length(dl-s_lag)))], na.rm=TRUE) - sum(vi_dataset[p,(s_pos:e_pos)], na.rm=TRUE)
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
                                 cyc_result[3] <- as.double(vi_dataset[p,s_pos])
                                 cyc_result[4] <- s_end
                                 cyc_result[5] <- as.integer(strftime(as.Date(s_end, origin="1970-01-01"), format = "%j"))
                                 cyc_result[6] <- as.double(vi_dataset[p,e_pos])
                                 
                                 # bind result to data.frame
                                 if(exists("pm_df")){
                                   pm_df <- rbind(pm_df, data.frame(t(c(cyc_result, pks_result))))
                                 } else {
                                   pm_df <- data.frame(t(c(cyc_result, pks_result)))
                                 }
                               }
                               
                             } else {
                               
                               # compute phenometrics
                               pks_result <- phenoGetEO4PM(pixel=vi_dataset[p,], obs_dates = phenodates, year=year, c_win=c_win)
                               
                               # update  Quality Flags
                               pks_result[37] <- as.integer(pks_result[37] + qf0 * 2^0 + qf11 * 2^11 + qf12 * 2^12 + qf14 * 2^14)
                               
                               # write values for Start of Cycle (SoC) and End of Cycle (EoC)
                               cyc_result <- rep(NA, 6)
                               names(cyc_result) <- c("SoC_date", "SoC_doy", "SoC_value", "EoC_date", "EoC_doy", "EoC_value")
                               rMin_pos <- phenofit::findpeaks(-vi_dataset[p,], minpeakdistance = 45)$X$pos
                               
                               if(!is.na(pks_result[15])){
                                 first_jan_date <- as.integer(as.Date(paste(as.integer(strftime(as.Date(pks_result[15], origin="1970-01-01"), format = "%Y")), "-01-01", sep=""), origin="1970-01-01"))
                               } else {
                                 first_jan_date <- as.integer(as.Date(paste(as.integer(strftime(as.Date(median(phenodates), origin="1970-01-01"), format = "%Y")), "-01-01", sep=""), origin="1970-01-01"))
                               }
                               
                               SoC_pos <- rMin_pos[which(rMin_pos <= which(phenodates == pks_result[5]))]
                               SoC_pos <- SoC_pos[which.min(abs(SoC_pos - which(phenodates == pks_result[5])))][1]
                               if(length(SoC_pos) == 1){
                                 if(!is.na(SoC_pos)){
                                   cyc_result[1] <- phenodates[SoC_pos]
                                   cyc_result[2] <- as.integer(1 - first_jan_date + phenodates[SoC_pos])
                                   cyc_result[3] <- as.double(vi_dataset[p,SoC_pos])
                                 }
                               }
                               EoC_pos <- rMin_pos[which(rMin_pos >= which(phenodates == pks_result[25]))]
                               EoC_pos <- EoC_pos[which.min(abs(EoC_pos - which(phenodates == pks_result[25])))][1]
                               if(length(EoC_pos) == 1){
                                 if(!is.na(EoC_pos)){
                                   cyc_result[4] <- phenodates[EoC_pos]
                                   cyc_result[5] <- as.integer(1 - first_jan_date + phenodates[EoC_pos])
                                   cyc_result[6] <- as.double(vi_dataset[p,EoC_pos])
                                 }
                               }
                               # convert to data frame
                               pm_df <- data.frame(t(c(cyc_result, pks_result)))
                             }
                             
                             # sort results according to 'STI' values
                             pm_df <- cbind(data.frame("cycle"=order(pm_df$STI, decreasing = TRUE), "season"=c(1:nrow(pm_df)), pm_df))
                             # # sort results according to 'SMV' values
                             # pm_df <- cbind(data.frame("cycle"=order(pm_df$SMV, decreasing = TRUE), "season"=c(1:nrow(pm_df)), pm_df))
                             if(is.unsorted(pm_df$cycle)){
                               pm_df <- pm_df[pm_df$cycle,]
                             }
                             
                             # #################
                             # ### not applicable to function
                             
                             # apply scale factor and offset
                             if(var_scaled){
                               var_scalist <- c(5,8,12,15,18,19,22,25,28,31,32,35,36,42,43)
                               pm_df[,var_scalist] <- suppressWarnings(round((pm_df[,var_scalist] - var_offset) / var_scale_factor, digits = 0))
                               var_mps <- c(39)
                               pm_df[,var_mps] <- suppressWarnings(round((pm_df[,var_mps] - mps_var_offset) / mps_var_scale_factor, digits = 0))
                               var_stilist <- c(40,41,44)
                               if(sti_var_offset != 0 & sti_var_scale_factor != 1){
                                 pm_df[,var_stilist] <- suppressWarnings(round((pm_df[,var_stilist] - sti_var_offset) / sti_var_scale_factor, digits = 0))
                               } else {
                                 pm_df[,var_stilist] <- suppressWarnings(round(pm_df[,var_stilist], digits = 0))
                               }
                             }
                             
                             # clean values out-of-range
                             for(i in 1:nrow(pm_df)){
                               pm_df[i,which(pm_df[i,] < -32768)] <- NA
                               pm_df[i,which(pm_df[i,] > 32768)] <- NA
                               # pm_df[i,which(is.infinite(pm_df[i,]))] <- NA
                             }
                             # update 'bit13' for incomplete values for all the metrics
                             pm_df[which(!complete.cases(pm_df)),45] <- as.integer(pm_df[which(!complete.cases(pm_df)),45] + 2^13)
                             
                             # apply offset to QFs
                             # pm_df[,37] <- pm_df[,37] -32768
                             
                             # store results to target matrix
                             for(s in 1:max(1, nseas)){
                               pks_dataset[[s]][q,] <- as.integer(pm_df[s,c(2:45)])
                             }
                             # update NA values for QF
                             for(s in 1:max_season){
                               pks_dataset[[s]][which(is.na(pks_dataset[[s]][,44])),44] <- 0
                             }
                           }
                           
                           # # ### DEV function
                           # pm_df <- cbind(data.frame("year"=as.integer(strftime(as.Date(pm_df$PoS_date, origin="1970-01-01"), format = "%Y")), pm_df))
                           
                           # clean workspace
                           suppressWarnings(rm(list=c("pm_df", "pks_result", "cyc_result", "SoC_pos", "EoC_pos", "sdt", "s", "nseas", "first_jan_date")))
                         }
                         
                         # clean workspace
                         suppressWarnings(rm(list=c("vi_dataset", "p_list", "g_pixel_list", "pixel_list")))
                         invisible(gc(verbose = FALSE))
                         
                         # ###########################
                         # write results to output NetCDF file
                         lck_out <- TRUE
                         output_file_chunk <- paste(output_file, sprintf(fmt="%00006i", l), ".nc", sep="")
                         
                         while(lck_out == TRUE){
                           if(!file.exists(output_file_lck)){
                             write.csv(as.character(""), file = output_file_lck)
                             
                             # rename output file
                             mv_outfile <- suppressWarnings(file.rename(from=output_file, to=output_file_chunk))
                             if(mv_outfile){
                               # open output NetCDF
                               ncout <- ncdf4::nc_open(filename = output_file_chunk, write = TRUE)
                               # on.exit(ncdf4::nc_close(ncout), add = TRUE)
                               
                               # write data
                               # three attempts to write file
                               nc_write_status <- NULL
                               nc_write_status <- tryCatch(expr=ncdf4::ncvar_put(ncout, varid=paste("VC1_", nc_var_names[1], sep=""), vals=as.integer(pks_dataset[[1]][,1]), start=c(rs,cs), count=c(chunks_nrow,chunks_ncol), verbose = FALSE), error = function(e){
                                 Sys.sleep(60)
                                 tryCatch(expr=ncdf4::ncvar_put(ncout, varid=paste("VC1_", nc_var_names[1], sep=""), vals=as.integer(pks_dataset[[1]][,1]), start=c(rs,cs), count=c(chunks_nrow,chunks_ncol), verbose = FALSE), error = function(e){
                                   Sys.sleep(60)
                                   tryCatch(expr=ncdf4::ncvar_put(ncout, varid=paste("VC1_", nc_var_names[1], sep=""), vals=as.integer(pks_dataset[[1]][,1]), start=c(rs,cs), count=c(chunks_nrow,chunks_ncol), verbose = FALSE), error = function(e){
                                     # ncdf4::nc_close(ncout)
                                     return(2)
                                   })
                                 })
                               })
                               
                               if(!is.null(nc_write_status)){
                                 message(paste("Failed to write chunk: ", l, sep=""))
                               } else {
                                 
                                 for(s in 1:max_season){
                                   for(n in 1:length(nc_var_names)){
                                     # message(paste("Variable ", n, "/", length(nc_var_names), " - ", paste("VC", s, "_", nc_var_names[n], sep=""), sep="")) # ### for debug
                                     ncdf4::ncvar_put(ncout, varid=paste("VC", s, "_", nc_var_names[n], sep=""), vals=pks_dataset[[s]][,n], start=c(rs,cs), count=c(chunks_nrow,chunks_ncol), verbose = FALSE)
                                   }
                                 }
                               }
                               
                               # close output NetCDF file
                               # ncdf4::nc_sync(ncout)
                               ncdf4::nc_close(ncout)
                               suppressWarnings(rm(list=c("ncout")))
                               invisible(gc(verbose=FALSE))
                               
                               mv_outfile <- suppressWarnings(file.rename(from=output_file_chunk, to=output_file))
                               lck_rm <- file.remove(output_file_lck)
                               lck_out <- FALSE
                               
                             } else {
                               lck_out <- TRUE
                               Sys.sleep(10)
                             }
                           } else {
                             lck_out <- TRUE
                             Sys.sleep(10)
                           }
                         }
                       }
                     }
                     # } # ### for debug
                     # display progress bar
                     if(progress){
                       prg(sprintf("x=%g", l))
                     }
                     
                     # clean workspace
                     suppressWarnings(rm(list=c("pks_dataset")))
                     terra::tmpFiles(remove = TRUE)
                     invisible(gc(full=TRUE, verbose=FALSE))
                     
                     # # cleanup RAM
                     # if (.Platform$OS.type == "unix") {
                     #   
                     #   # # call malloc_trim(0) via system call to force RAM cleaning
                     #   # force_trim <- function() {
                     #   #   # Carica libc e chiama malloc_trim
                     #   #   invisible(.C("malloc_trim", as.integer(0), PACKAGE = "base"))
                     #   # }
                     #   
                     #   # define malloc_trim(0) C function to force RAM cleaning
                     #   Rcpp::cppFunction('
                     #   #include <malloc.h>
                     # 
                     #   int force_trim() {
                     #     malloc_trim(0);
                     #     return 1;
                     #   }'
                     #   )
                     #   
                     #   # force RAM cleaning
                     #   ft <- force_trim()
                     # }
                     
                     return(chunk_counter)
                   }
}

chunk_counter <- 0
invisible(gc(verbose=FALSE))

# ################################
# initialize future for parallel computing (parallelize using doFuture)
# future::plan(multisession, workers=cores, gc=TRUE)
# Enable garbage collection
# options(future.gc = TRUE)
# options(future.scheduling.gc = TRUE)

# future::plan(multicore, workers=cores, gc=TRUE)

future::plan(future.callr::callr, workers=cores)
options(future.gc = TRUE)
# oopts <- options(future.globals.maxSize = 1.5 * 1000000000)  ## 1.5 GB
# on.exit(options(oopts))

# run in parallel
# Rprof(append = TRUE, filename=paste0(output_file, ".out"), memory.profiling = TRUE)
result_par <- analysis_fun(r_list)
# Rprof(NULL)
# stop cluster
future::plan(sequential)

# remove temporary raster mask
if(file.exists(normalizePath(paste(tmp_folder, "/", basename(output_file), "_rastMask.tif", sep=""), winslash = "/", mustWork = FALSE))){
  r_mask_rm <- file.remove(normalizePath(paste(tmp_folder, "/", basename(output_file), "_rastMask.tif", sep=""), winslash = "/", mustWork = FALSE))
}

chunk_counter <- sum(as.integer(result_par))
invisible(gc())

############################################################
# Include additional attributes
############################################################

# set creation date
get_systime <- Sys.time()
nc_attr_isodate <- as.character(ISOdatetime(year = strftime(get_systime, format = "%Y", tz = "GMT"), month = strftime(get_systime, format = "%m", tz = "GMT"), day = strftime(get_systime, format = "%d", tz = "GMT"), hour = strftime(get_systime, format = "%H", tz = "GMT"), min = strftime(get_systime, format = "%M", tz = "GMT"), sec = strftime(get_systime, format = "%S", tz = "GMT"), tz = "GMT"))
# nc_hist_systime <- as.character(format(Sys.time(), "%a %b %d %H:%M:%S %Y"))
get_systime_split <- unlist(strsplit(as.character(get_systime), split=" "))
date_creation <- paste(get_systime_split[1], "T", get_systime_split[2], "+0000", sep="")

Sys.sleep(10)
iter_loop <- 1
while(!file.exists(output_file) && iter_loop < 6) {
  Sys.sleep(5)
  iter_loop <- iter_loop + 1
}
ncout <- ncdf4::nc_open(filename = output_file, write = TRUE, readunlim = FALSE)

ncdf4::ncatt_put(ncout,0,"date_created",date_creation)
ncdf4::ncatt_put(ncout,0,"ISOdate",nc_attr_isodate)

# create history line
nc_hist <- paste(nc_in_hist, "\n", nc_hist_systime, ": ", nc_hist_past, sep="\n")
ncdf4::ncatt_put(ncout, 0,"history", nc_hist)

# close output NetCDF file
ncdf4::nc_sync(ncout)
ncdf4::nc_close(ncout)

if(verbose){
  message(paste("Ending Phenological Metrics estimation: ", Sys.time(), sep=""))
  message(paste("Elapsed time ", as.character(paste(as.integer(as.numeric(proc.time() - ptm)[3]/3600), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/3600 - as.integer(as.numeric(proc.time() - ptm)[3]/3600)) * 60)), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/60 - as.integer(as.numeric(proc.time() - ptm)[3]/60)) * 60)), sep="")), " hours", sep=""))
  message("")
}

# clean workspace and exit
rm(list=ls())
invisible(gc())
q("no")
