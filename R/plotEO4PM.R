# Title         : plotEO4PM.R
# Description   : Plot Phenological Metrics generated using EO4PM
# Date          : Jun 2026
# Version       : 1.2
# Copyright name: CC BY
# Licence       : GPL v3
# Authors       : Federico Filipponi
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
# ##############################################################################
# ChangeLog
# 30/03/2026
# Change scale_x_datetime() to scale_x_date() to fix problems with x-axis
#
#' @title plotEO4PM
#'
#' @description Plot Phenological Metrics generated using EO4PM, along with 
#' interpolated time series and optionally raw time series.
#'
#' @param ts Numeric vector, vegetation index time-series or object of class \code{\linkS4class{data.frame}} generated using function \code{\link{smoothEO4PM}}.
#' @param pm (optional) Object of class \code{\linkS4class{data.frame}} generated using function \code{\link{phenoEO4PM}}.
# #' @param raw (optional) Numeric vector, raw vegetation index time-series.
#' @param start_date (optional) Character. Input start date, either an object of class \code{\linkS4class{Date}} or a character object with format 'YYYY-MM-DD'.
#' @param end_date (optional) Character. Input end date, either an object of class \code{\linkS4class{Date}} or a character object with format 'YYYY-MM-DD'.
#' @param legend Logical. Display a legend on plot bottom. If FALSE ggplot object is returned from function. Default: TRUE.
#' @param title (optional) Character. Name to be used for plot title.
#' @param ylab (optional) Character. Plot y-axis label.
#' @param outfile (optional) Character. Output file name where to store plots.
# #' @param verbose Logical. Set verbose mode to enable warning messages.
# #' @param ... Additional arguments to be passed through to function \code{\link{plotEO4PM}}
#'
#' @return
#' The function returns one plot with phenological metrics over source time series. 
#' Optionally plots can be saved to output folder set using dest argument.
#'
#' @author Federico Filipponi
#'
#' @keywords EO4PM phenology
#'
# #' @seealso
#'
#' @examples
#' \dontrun{
#' # temporal smoothing
#' y <- jitter(dnorm(seq(-4, 4, length.out = 500), mean = 0, sd = 1)*4, amount = 0.1)
#' tsg <- smoothEO4PM(y = y, t = 1:length(y))
#' # estimate phenological metrics
#' phenometrics <- phenoEO4PM(y = tsg$y, t = tsg$t)
#' # generate plot with minimal set of parameters
#' plotEO4PM(ts = tsg, pm = phenometrics)
#' 
#' # generate plot with additional parameters
#' plotEO4PM(ts = tsg, pm = phenometrics, legend = FALSE, ylab = "VI", title = "My plot", subtitle = "of phenology")
#' 
#' # generate plot for specific temporal range
#' plotEO4PM(ts = tsg, pm = phenometrics, start_date="1970-03-01", end_date="1971-04-01")
#' 
#' # save plot to file
#' plotEO4PM(ts = tsg, pm = phenometrics, outfile="plot.png")
#' 
#' # store plot to object and replot
#' EO4PM_ggplot <- plotEO4PM(ts = tsg, pm = phenometrics)
#' plot(EO4PM_ggplot)
#' }
#' 
#' @import ggplot2
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_x_date
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 ggplot_gtable
#' @importFrom gridExtra grid.arrange
# #' @export

plotEO4PM <- function(ts, pm=NULL, start_date=NULL, end_date=NULL, title=NULL, subtitle=NULL, ylab=NULL, legend=TRUE, outfile=NULL){
  
  # check input arguments
  if(!is.data.frame(ts)){
    stop("Input 'ts' object must by of type 'data.frame' with the following fields: 't', 'y'.")
  } else {
    if(sum(as.integer(c('t', 'y') %in% names(ts))) != 2){
      stop("Input 'ts' object must by of type 'data.frame' with the following fields: 't', 'y'.")
    }
  }
  
  # check input arguments
  if(!is.null(pm)){
    if(!is.data.frame(pm)){
      stop("Input 'pm' object must by of type 'data.frame' generated from function 'phenoEO4PM'.")
    } else {
      if(sum(as.integer(c('SoC_date','EoC_date','SoS_date','SGS_date','greenup_date','SMP_date','PoS_date','EGS_date','senescence_date','EoS_date','SoS_value','SGS_value','greenup_value','SMP_value','PoS_value','EGS_value','senescence_value','EoS_value') %in% names(pm))) < 18){
        stop("Input 'pm' object must by of type 'data.frame' generated using function 'phenoEO4PM' with the following fields: 'SoC_date','EoC_date','SoS_date','SGS_date','greenup_date','SMP_date','PoS_date','EGS_date','senescence_date','EoS_date','SoS_value','SGS_value','greenup_value','SMP_value','PoS_value','EGS_value','senescence_value','EoS_value'.")
      }
    }
  }
  
  # define function to check dates
  IsDate <- function(mydate, date.format = "%Y-%m-%d"){
    tryCatch(!is.na(as.Date(mydate, date.format)),  
             error = function(err) {FALSE})  
  }
  
  if(!is.null(start_date)){
    # check if date has the correct format
    if(!IsDate(start_date)){
      stop("Argument 'start_date' has not the format 'YYYY-MM-DD'.", call.=FALSE)
    }
    if(as.integer(as.Date(ts$t[nrow(ts)], tz="UTC")) < as.integer(as.Date(start_date, tz="UTC"))){
      stop("Input time series is not within temporal range set using argument 'start_date'.", call.=FALSE)
    }
  }
  if(!is.null(end_date)){
    # check if date has the correct format
    if(!IsDate(end_date)){
      stop("Argument 'end_date' has not the format 'YYYY-MM-DD'.", call.=FALSE)
    }
    if(!is.null(start_date)){
      if(as.Date(start_date) >= as.Date(end_date)){
        stop("Argument 'start_date' must be set to a date preceding argument 'end_date'.", call.=FALSE)
      }
      # check if temporal window falls outside input time series data
      if(as.integer(as.Date(ts$t[1], tz="UTC")) > as.integer(as.Date(end_date, tz="UTC")) | as.integer(as.Date(ts$t[nrow(ts)], tz="UTC")) < as.integer(as.Date(start_date, tz="UTC"))){
        stop("Input time series is not within temporal range set using arguments 'start_date' and 'end_date'.", call.=FALSE)
      }
    } else {
      if(as.integer(as.Date(ts$t[1], tz="UTC")) > as.integer(as.Date(end_date, tz="UTC"))){
        stop("Input time series is not within temporal range set using argument 'end_date'.", call.=FALSE)
      }
    }
  }
  
  start_date <- as.Date(max(ts$t[1], as.integer(as.Date(start_date)), na.rm=TRUE))
  end_date <- as.Date(min(ts$t[nrow(ts)], as.integer(as.Date(end_date)), na.rm=TRUE))
  
  if(!("w" %in% names(ts))){
    ts$w <- rep(as.double(0), nrow(ts))
  }
  
  ts$t <- as.Date(ts$t, tz="UTC")
  
  # # load ggplot2 library
  # suppressPackageStartupMessages(library(ggplot2))
  
  # build plot
  g <- ggplot2::ggplot()
  
  if(!is.null(pm)){
    # fill areas outside identified seasons in grey
    g <- g + ggplot2::geom_rect(aes(xmin=as.Date(ts$t[1]), xmax=as.Date(pm$SoC_date[1]), ymin=-Inf, ymax=Inf), alpha=0.2, colour="transparent", fill="black")
    if(nrow(pm) > 1){
      for(s in 1:(nrow(pm)-1)){
        if(as.Date(pm$SoC_date[s+1]) > as.Date(pm$EoC_date[s])){
          g <- g + suppressWarnings(ggplot2::geom_rect(aes_(xmin=as.Date(pm$EoC_date[s]), xmax=as.Date(pm$SoC_date[s+1]), ymin=-Inf, ymax=Inf), alpha=0.2, colour="transparent", fill="black"))
        }
      }
    }
    g <- g + ggplot2::geom_rect(aes(xmin=as.Date(pm$EoC_date[nrow(pm)]), xmax=as.Date(ts$t[nrow(ts)]), ymin=-Inf, ymax=Inf), alpha=0.2, colour="transparent", fill="black")
    
    # plot time series later replaced by asymmetric gaussian function
    if("QF" %in% names(pm)){
      bit14 <- floor(pm$QF / 2^14) %% 2
      if(any(bit14 == 1)){
        ts$dotted <- rep(as.double(NA), nrow(ts))
        for(b in which(bit14 == 1)){
          ts_pos <- sort(which(ts$t %in% (pm$SoC_date[b]:pm$EoC_date[b])))
          ts$dotted[ts_pos] <- ts$y[ts_pos]
          dl <- tryCatch(expr = phenofit::FitDL.AG(y = ts$y[ts_pos], t = ts$t[ts_pos], tout = ts$t[ts_pos], w=ifelse(ts$y[ts_pos] > quantile(ts$y[ts_pos], probs=0.75, na.rm=TRUE), 1.0, 0.25), type = 1, iters = 2), error = function(e){return(NULL)})
          if(!is.null(dl)){
            ts$y[ts_pos] <- dl$zs$iter2
          }
        }
      }
    }
  }
  
  # plot time series
  g <- g + ggplot2::geom_line(data=ts, ggplot2::aes(t, y), linetype="solid", color="green4", linewidth=1.2)
  
  # plot original smoothed curve (when asymmetric gaussian was applied)
  if("dotted" %in% names(ts)){
    g <- g + ggplot2::geom_line(data=ts, ggplot2::aes(t, dotted), linetype="dotted", color="darkviolet", linewidth=0.6)
  }
  
  # plot raw time series
  raw <- FALSE
  if("raw" %in% names(ts)){
    raw <- TRUE
    g <- g + ggplot2::geom_point(data=ts, ggplot2::aes(t, raw, alpha=w), size=0.4, fill="black", na.rm=TRUE, show.legend = FALSE)
  }
  
  if(!is.null(pm)){
    # add cycle boundaries
    g <- g + ggplot2::geom_vline(data = pm, ggplot2::aes(xintercept = as.Date(EoC_date)), color="darkorange3", linetype="solid", linewidth=0.4)
    g <- g + ggplot2::geom_vline(data = pm, ggplot2::aes(xintercept = as.Date(SoC_date)), color="blue4", linetype="dotted", linewidth=0.6)
    
    # add PM symbols
    scaled_size <- max(1.5, min(4, (1.5 + (nrow(ts) - 365) * (4 - 1.5) / (5500 - 365))))
    g <- g + ggplot2::geom_point(data=pm, ggplot2::aes(as.Date(SoS_date), SoS_value), shape=23, color="black", size=scaled_size, na.rm=TRUE) + scale_color_manual(labels = c("SoS"))
    g <- g + ggplot2::geom_point(data=pm, ggplot2::aes(as.Date(SGS_date), SGS_value), shape=23, color="cyan3", size=scaled_size, na.rm=TRUE)
    g <- g + ggplot2::geom_point(data=pm, ggplot2::aes(as.Date(greenup_date), greenup_value), shape=21, color="green2", size=scaled_size, na.rm=TRUE)
    g <- g + ggplot2::geom_point(data=pm, ggplot2::aes(as.Date(SMP_date), SMP_value), shape=22, color="chartreuse3", size=scaled_size, na.rm=TRUE)
    g <- g + ggplot2::geom_point(data=pm, ggplot2::aes(as.Date(PoS_date), PoS_value), shape=24, color="blue3", size=scaled_size, na.rm=TRUE)
    g <- g + ggplot2::geom_point(data=pm, ggplot2::aes(as.Date(EGS_date), EGS_value), shape=22, color="gold", size=scaled_size, na.rm=TRUE)
    g <- g + ggplot2::geom_point(data=pm, ggplot2::aes(as.Date(senescence_date), senescence_value), shape=21, color="orangered", size=scaled_size, na.rm=TRUE)
    g <- g + ggplot2::geom_point(data=pm, ggplot2::aes(as.Date(EoS_date), EoS_value), shape=23, color="red3", size=scaled_size, na.rm=TRUE)
  }
  
  # add extra graphics
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(panel.grid.major.x = element_line(color = "black", linewidth = 0.2, linetype="solid"), panel.grid.minor.x = element_line(color = "grey", linewidth = 0.2, linetype="dotted"), panel.grid.minor.y =element_blank())
  ymax <- max(ts$y[which(as.integer(as.Date(ts$t)) %in% (as.integer(start_date:end_date)))], na.rm=TRUE)*1.05
  g <- g + ggplot2::ylim(0, ymax)
  g <- g + ggplot2::scale_x_date(limits=c(as.Date(start_date), as.Date(end_date)),  
                                     date_minor_breaks= "1 months",
                                     date_breaks= "1 year", expand=c(0,0), date_labels="%Y")
  g <- g + ggplot2::xlab("")
  
  if(!is.null(ylab)){
    g <- g + ggplot2::ylab(ylab)
  }
  
  # add title
  if(!is.null(title)){
    if(is.null(subtitle)){
      g <- g + ggplot2::ggtitle(title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    } else {
      g <- g + ggplot2::ggtitle(title, subtitle = subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5, face="italic"))
    }
  }
  
  # generate legend
  if(legend){
    # # load gridExtra library
    # suppressPackageStartupMessages(library(grid))
    # suppressPackageStartupMessages(library(gridExtra))
    
    if(raw & !is.null(pm)){
      g_group <- c("TS","RAW","SoC","EoC")
      scm <- scale_color_manual(
        # name = "Legend",
        breaks= c("RAW","TS","SoC","EoC"),
        values = c("RAW"="black", "TS"="green4", "SoC"="blue4","EoC"="darkorange3"),
        labels = c("RAW"="raw data", "TS"="smoothed time series", "SoC"="SoC","EoC"="EoC")
      )
      g_guides <- guides(color = guide_legend(nrow = 1, 
                                              override.aes = list(
                                                linetype = c("blank","solid","dotted","solid"), # 'blank' nasconde la linea
                                                linewidth = c(NA,1.5,0.8,0.8),
                                                shape    = c(16, NA, NA, NA), # NA nasconde il punto (16=cerchio)
                                                size     = c(2, rep(3, 3))
                                              )
      ))
    }
    if(!raw & !is.null(pm)){
      g_group <- c("TS","SoC","EoC")
      scm <- scale_color_manual(
        # name = "Legend",
        breaks= c("TS","SoC","EoC"),
        values = c("TS"="green4", "SoC"="blue4","EoC"="darkorange3"),
        labels = c("TS"="smoothed time series", "SoC"="SoC","EoC"="EoC")
      )
      g_guides <- guides(color = guide_legend(nrow = 1, 
                                              override.aes = list(
                                                linetype = c("solid","dotted","solid"), # 'blank' nasconde la linea
                                                linewidth = c(1.5,0.8,0.8),
                                                shape    = c(NA, NA, NA), # NA nasconde il punto (16=cerchio)
                                                size     = c(rep(2, 3))
                                              )
      ))
    }
    if(raw & is.null(pm)){
      g_group <- c("TS","RAW")
      scm <- scale_color_manual(
        # name = "Legend",
        breaks= c("RAW","TS"),
        values = c("RAW"="black", "TS"="green4"),
        labels = c("RAW"="raw data", "TS"="smoothed time series")
      )
      g_guides <- guides(color = guide_legend(nrow = 1, 
                                              override.aes = list(
                                                linetype = c("blank","solid"), # 'blank' nasconde la linea
                                                linewidth = c(NA,1.5),
                                                shape    = c(16, NA), # NA nasconde il punto (16=cerchio)
                                                size     = c(2, 3)
                                              )
      ))
    }
    if(!raw & is.null(pm)){
      g_group <- c("TS")
      scm <- scale_color_manual(
        # name = "Legend",
        breaks= c("TS"),
        values = c("TS"="green4"),
        labels = c("TS"="smoothed time series")
      )
      g_guides <- guides(color = guide_legend(nrow = 1, 
                                              override.aes = list(
                                                linetype = c("solid"), # 'blank' nasconde la linea
                                                linewidth = c(1.5),
                                                shape    = c(NA), # NA nasconde il punto (16=cerchio)
                                                size     = c(2)
                                              )
      ))
    }
    
    # create time series legend
    p_legenda_ts <- ggplot(data.frame(x=1, y=1:length(g_group), gruppo=g_group), 
                           aes(x, y, color=gruppo)) +
      geom_point() +
      geom_line(group=1) +
      scm +
      g_guides + 
      theme(legend.background = element_rect(fill="transparent", color=NA), legend.key = element_rect(fill = NA, color = NA), legend.direction = "horizontal", legend.title = element_blank(), legend.text=element_text(size=14), plot.margin = unit(c(0,0,0,0), "pt"))
    
    if(!is.null(pm)){
      # create phenological metrics legend line
      p_legenda_pm <- ggplot2::ggplot(data.frame(x=1, y=1:8, gruppo=c("SoS","SGS","greenup","SMP","PoS","EGS","senescence","EoS")), 
                                      aes(x, y, color=gruppo)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_manual(
          name = "Phenological Metrics",
          breaks= c("SoS","SGS","greenup","SMP","PoS","EGS","senescence","EoS"),
          values = c("SoS"="black","SGS"="cyan3", "greenup"="green2","SMP"="chartreuse3", "PoS"="blue3", "EGS"="gold", "senescence"="orangered", "EoS"="red3"),
          labels = c("SoS"="SoS","SGS"="SGS", "greenup"="greenup","SMP"="SMP", "PoS"="PoS", "EGS"="EGS", "senescence"="senescence", "EoS"="EoS")
        ) +
        
        ggplot2::guides(color = guide_legend(nrow = 1, 
                                             override.aes = list(
                                               linetype = c("blank","blank","blank","blank","blank","blank","blank","blank"), 
                                               shape    = c(23, 23, 21, 22, 24, 22, 21, 23),                
                                               size     = c(rep(4, 8))
                                             )
        )) +
        ggplot2::theme(legend.background = element_rect(fill="transparent", color=NA), legend.key = element_rect(fill = NA, color = NA), legend.title=element_text(size=14), legend.text=element_text(size=14), legend.direction = "horizontal", legend.box.margin = margin(0,0,0,0), plot.margin = unit(c(0,0,0,0), "pt"))
      
      # transform to 'gtable' object
      g_pm <- ggplot_gtable(ggplot_build(p_legenda_pm))
      g_pm_leg_idx <- which(sapply(g_pm$grobs, function(x) x$name) == "guide-box")
      pm_legend <- g_pm$grobs[[g_pm_leg_idx]]
    }
    
    # transform to 'gtable' object
    g_ts <- ggplot2::ggplot_gtable(ggplot_build(p_legenda_ts))
    g_ts_leg_idx <- which(sapply(g_ts$grobs, function(x) x$name) == "guide-box")
    ts_legend <- g_ts$grobs[[g_ts_leg_idx]]
    
    if(!is.null(outfile)){
      # set output file name
      outfile <- normalizePath(path=outfile, winslash = "/", mustWork = FALSE)
      if(!dir.exists(dirname(outfile))){
        dir.create(path=dirname(outfile), showWarnings = FALSE, recursive = TRUE)
      }
      # initialize PNG file
      png(filename=outfile, width = 6000, height = 2000, units = "px", pointsize = 12, bg = "white", res = 300)
      # print plot
      # suppressWarnings(print(g))
      if(!is.null(pm)){
        suppressWarnings(gridExtra::grid.arrange(grobs=list(g, ts_legend, pm_legend), ncol = 1, nrow=3, heights = grid::unit(c(10, 0.8, 0.5), units = "null")))
      } else {
        suppressWarnings(gridExtra::grid.arrange(grobs=list(g, ts_legend), ncol = 1, nrow=3, heights = grid::unit(c(10, 0.8, 0.5), units = "null")))
      }
      # write to PNG file
      invisible(dev.off())
    } else {
      if(!is.null(pm)){
        # gr <- gridExtra::arrangeGrob(g, ts_legend, pm_legend, ncol = 1, nrow=3, heights = grid::unit(x = c(10, 0.8, 0.5), units="null"))
        g <- suppressWarnings(gridExtra::grid.arrange(grobs=list(g, ts_legend, pm_legend), ncol = 1, nrow=3, heights = grid::unit(c(10, 0.8, 0.5), units = "null")))
      } else {
        # gr <- gridExtra::arrangeGrob(g, ts_legend, ncol = 1, nrow=3, heights = grid::unit(x = c(10, 0.8, 0.5), units="null"))
        g <- suppressWarnings(gridExtra::grid.arrange(grobs=list(g, ts_legend), ncol = 1, nrow=3, heights = grid::unit(c(10, 0.8, 0.5), units = "null")))
      }
      
      return(g)
    }
  
  } else {
    if(!is.null(outfile)){
      # set output file name
      outfile <- normalizePath(path=outfile, winslash = "/", mustWork = FALSE)
      if(!dir.exists(dirname(outfile))){
        dir.create(path=dirname(outfile), showWarnings = FALSE, recursive = TRUE)
      }
      # initialize PNG file
      png(filename=outfile, width = 6000, height = 2000, units = "px", pointsize = 12, bg = "white", res = 300)
      # print plot
      suppressWarnings(print(g))
      # write to PNG file
      invisible(dev.off())
    } else {
      # display plot
      suppressWarnings(print(g))
      return(g)
    }
  }
}
