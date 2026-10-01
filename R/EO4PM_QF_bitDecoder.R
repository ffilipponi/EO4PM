# Title         : EO4PM_QF_bitDecoder.R
# Description   : Decode bit Quality Flags generated during estimation of phenological metrics using EO4PM
# Date          : Feb 2026
# Version       : 1.0
# Copyright name: CC BY
# Licence       : GPL v3
# Authors       : Federico Filipponi
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
# ##############################################################################
#' @title QFbitDecoder
#' 
#' @description Decode bit Quality Flags 
#' generated during estimation of phenological metrics using EO4PM.
#'
#' @param x Object of class \code{\linkS4class{data.frame}} generated using function \code{\link{phenoEO4PM}}.
#' @param append Logical. Append fields of decoded bit to input data frame (default: TRUE).
#'
#' @return
#' The function returns data.frame with QF decoded bit. 
#' Optionally new fields are appended to input EO4PM data.frame.
#' 
#' @details
#' Quality Flags meaning
#' bit0: Processed seasonal cycle
#' bit1: Pixel or date input problem
#' bit2: First derivative error
#' bit3: Identified peak date is before reference year (before 1st January)
#' bit4: No increasing period found
#' bit5: No decreasing period found
#' bit6: Gu stages identification failed
#' bit7: Estimated Gu stages are not chronologically ordered
#' bit8: Estimated Gu stages are not within input temporal range or have missing metrics values
#' bit9: Cannot detect relative minimum (SoS) using inflection point
#' bit10 Identified Gu stages are not within the analyzed time period (Identified SoS or SGS or EGS or EoS date falls outside season temporal range)
#' bit11: Failed to cut cycles
#' bit12: Identified more than specified max_season
#' bit13: Incomplete set of Phenological Metrics (e.g. due to out-of-range issues)
#' bit14: Used Asymmetric Gaussian fit after Gu phenometric estimate failure (only applied after successfull cut-cycle)
#'
#' @author Federico Filipponi
#'
#' @keywords EO4PM phenology bit QF decoder
#' 
#' @family EO4PM
#'
# #' @seealso
#'
#' @examples
#' \dontrun{
#' # decode QF bits
#' y <- dnorm(seq(-4, 4, length.out = 500), mean = 0, sd = 1)*4
#' phenometrics <- phenoEO4PM(y = y, t = 1:length(y))
#' # decoded bits
#' EO4PM_QF <- QFbitDecoder(x = phenometrics, append = FALSE)
#' # append decoded bit fields to phenometrics data.frame
#' phenometrics <- QFbitDecoder(x = phenometrics)
#' }
#' 
# #' @export

# ###########################################################
# Define function to decode bits
# ###########################################################

QFbitDecoder <- function(x, append=TRUE){
  
  # check input arguments
  if(!is.data.frame(x)){
    stop("Input 'x' object must by of type 'data.frame' with the following field: 'QF'.")
  } else {
    if(!("QF" %in% names(x))){
      stop("Input 'x' object must by of type 'data.frame' with the following field: 'QF'.")
    }
  }
  
  # check if bit 
  if(append){
    if(sum(as.integer(names(x) %in% paste("QF_bit", c(0:15), sep=""))) > 0){
      stop("Input 'x' object already contains at least one decoded bit field.")
    }
  }
  
  # create data frame to store decoded bits
  bdf <- data.frame(matrix(data = as.integer(NA), nrow = nrow(x), ncol = 16))
  names(bdf) <- paste("QF_bit", c(0:15), sep="")
  
  # decode bits
  for(bit in 0:15){
    bdf[,(bit+1)] <- as.integer(floor(x$QF/(2^bit)) %% 2 == 1)
  }
  
  if(append){
    return(cbind(x,bdf))
  } else {
    return(bdf)
  }
}

# ###########################################################
# Define function to show help in Rstudio viewer
# ###########################################################

source_help <- function(x) {
  f_obj <- get(x, envir = .GlobalEnv)
  src_info <- attr(f_obj, "srcref")
  file_path <- attr(src_info, "srcfile")$filename
  lines <- readLines(file_path)
  func_line <- grep("function\\(", lines)[1]
  comments <- c()
  for(i in 1:(func_line-1)){
    if(grepl("^\\s*#'", lines[i])){
      comments <- c(comments, lines[i])
    }
  }
  if (length(comments) == 0) stop("No comment #' found.")
  clean <- gsub("^\\s*#'\\s*", "", comments)
  clean <- clean[!grepl("^@(import|importFrom)", clean)]
  all_tags_idx <- grep("^@", clean)
  extract_tag_block <- function(tag_name) {
    start <- grep(paste0("^", tag_name), clean)
    if (length(start) == 0) return("")
    next_tags <- all_tags_idx[all_tags_idx > start]
    end <- if (length(next_tags) > 0) next_tags[1] - 1 else length(clean)
    block <- clean[start:end]
    block[1] <- gsub(paste0("^", tag_name, "\\s*"), "", block[1])
    return(paste(block, collapse = "\n\n"))
  }
  title_val <- clean[1]
  first_tag <- if(length(all_tags_idx) > 0) min(all_tags_idx) else length(clean) + 1
  descr_text <- if(first_tag > 2) paste(clean[2:(first_tag-1)], collapse = "\n\n") else ""
  descr_explicit <- extract_tag_block("@description")
  final_description <- if(nchar(descr_explicit) > 0) descr_explicit else descr_text
  final_description <- gsub("\n", " ", final_description)
  details_val  <- extract_tag_block("@details")
  return_val   <- extract_tag_block("@return")
  author_val   <- extract_tag_block("@author")
  keywords_val <- extract_tag_block("@keywords")
  examples_val <- extract_tag_block("@examples")
  family_val   <- extract_tag_block("@family")
  param_lines <- clean[grep("^@param", clean)]
  params_rd <- gsub("@param\\s+(\\w+)\\s+(.*)", "\\\\item{\\1}{\\2}", param_lines)
  usage_text <- capture.output(args(f_obj))
  usage_full <- gsub("^function", x, paste(usage_text[-length(usage_text)], collapse = "\n"))
  rd_txt <- c(
    paste0("\\name{", x, "}"),
    paste0("\\alias{", x, "}"),
    paste0("\\title{", gsub("@title", "", title_val), "}"),
    paste0("\\description{", final_description, "}")
  )
  if (nchar(details_val) > 0) {
    rd_txt <- c(rd_txt, paste0("\\details{", details_val, "}"))
  }
  rd_txt <- c(rd_txt, paste0("\\usage{", usage_full, "}"))
  if (length(params_rd) > 0)   rd_txt <- c(rd_txt, "\\arguments{", params_rd, "}")
  if (nchar(return_val) > 0)   rd_txt <- c(rd_txt, paste0("\\value{", gsub("\n", "", return_val), "}"))
  if (nchar(author_val) > 0)   rd_txt <- c(rd_txt, paste0("\\author{", author_val, "}"))
  if (nchar(keywords_val) > 0) rd_txt <- c(rd_txt, paste0("\\keyword{", keywords_val, "}"))
  if (nchar(examples_val) > 0) {
    ex_clean <- gsub("\n\n", "\n", examples_val)
    rd_txt <- c(rd_txt, paste0("\\examples{", ex_clean, "}"))
  }
  rd_con <- textConnection(rd_txt)
  rd_obj <- tools::parse_Rd(rd_con)
  close(rd_con)
  out_html <- tempfile(fileext = ".html")
  tools::Rd2HTML(rd_obj, out = out_html)
  viewer <- getOption("viewer", utils::browseURL)
  viewer(out_html)
}
