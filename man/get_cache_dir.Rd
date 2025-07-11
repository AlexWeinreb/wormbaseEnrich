% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_dictionary.R
\name{get_cache_dir}
\alias{get_cache_dir}
\title{Determine cache directory}
\usage{
get_cache_dir(cache_dir)
}
\arguments{
\item{cache_dir}{User-specified cache directory.}
}
\value{
The location of the cache directory, either explicitly user-provided,
specified as a global option, or (by default) in the OS cache directory.
}
\description{
This internal function determines what directory should be used to store cache
files. It can be a user-provided path, set as a global option, or determined
on an OS basis.
}
