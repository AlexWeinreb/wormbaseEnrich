#' Determine cache directory
#'
#' This internal function determines what directory should be used to store cache
#' files. It can be a user-provided path, set as a global option, or determined
#' on an OS basis.
#' @param cache_dir User-specified cache directory.
#'
#' @return The location of the cache directory, either explicitly user-provided,
#' specified as a global option, or (by default) in the OS cache directory.
#'
get_cache_dir <- function(cache_dir){
  if(is.null(cache_dir)){
    if(!is.null(getOption("WBenrich_cache_dir"))){
      cache_dir <- getOption("WBenrich_cache_dir")
    } else{
      cache_dir <- rappdirs::user_cache_dir("WBenrich", "WBenrich")
    }
  }

  if(!dir.exists(cache_dir)){
    dir.create(cache_dir, recursive = TRUE)
    message("Created cache directory: ", cache_dir)
  }

  cache_dir
}


download_dictionary <- function(filename, cache_dir){

  url_base <- 'http://caltech.wormbase.org/TissueEnrichmentAnalysis/'

  url <- paste0(url_base, filename)

  utils::download.file(url,
                file.path(cache_dir, filename))

}






#' Fetch a dictionary
#'
#' Load a dictionary file from Wormbase
#'
#'
#' @param analysis the type of analysis (tissue, phenotype, go).
#' @param cache If a number, only download the dictionary after this many days.
#' @param cache_dir Specify a cache directory.
#'
#' # Caching
#'
#' Since the Wormbase annotations change relatively rarely, it may be not be desirable
#' to download the dictionary after every session refresh. Instead, on the first call to
#' this function, the dictionary is downloaded and cached (saved to disk). Subsequent
#' calls to this function first attempt to read a cached file.
#'
#'
#' If `cache = FALSE` or `cache = 0`, the dictionary is downloaded, saved in `cache_dir`, and returned.
#' If `cache` is a (positive) number, the content of `cache_dir` is examined, and if
#' a file less than `cache` days is found, it is read and returned. If a file older than `cache` days
#' is found, a new dictionary is downloaded and overwrites the existing one.
#'
#' To keep an analysis reproducible on the long term, set `cache = Inf` and a local (version-controlled) `cache_dir`.
#'
#' ## Cache directory
#'
#' If `cache_dir` is a character, it is taken as the path to a directory to use for caching. Otherwise,
#' if the option `WBenrich_cache_dir` is set, it is used as cache directory. Otherwise,
#' the cache directory is chosen by `rappdirs`.
#'
#' @returns A data.frame containing the dictionary of interest, one row per gene, one column per term
#'
#' @export
fetch_dictionary <- function(analysis = c("tissue", "phenotype", "go"), cache = 2L, cache_dir = NULL) {

  analysis <- tolower(analysis)
  if(length(analysis) != 1L){
    stop("Only provide one type of analysis")
  }
  analysis <- match.arg(analysis)

  cache_dir <- get_cache_dir(cache_dir)

  if(! cache >= 0){
    stop("Value of `cache` incorrect: ", cache)
  }



  filename <- switch (analysis,
                      "tissue" = "anatomy_dict.csv",
                      "phenotype" = "phenotype_dict.csv",
                      "go" = "go_dict.csv"
  )

  filepath <- file.path(cache_dir, filename)

  if(file.exists(filepath)){
    file_age <- as.numeric(difftime(Sys.time(), file.mtime(filepath), units = "days"))
  }

  if( ! file.exists(filepath) || file_age > cache){
    download_dictionary(filename, cache_dir)
  }

  read_dictionary(filepath)

}

read_dictionary <- function(filepath){

  utils::read.csv(filepath, check.names = FALSE)
}





