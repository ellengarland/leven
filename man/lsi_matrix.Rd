% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leven.R
\name{lsi_matrix}
\alias{lsi_matrix}
\title{Compute the Levenshtein Similarity Index for each possible pair of strings in the input file}
\usage{
lsi_matrix(filename)
}
\arguments{
\item{filename}{The file containing the input data. The input must be comma-separated values, where the columns are:
\itemize{ \item Location 
\item Year 
\item Song type 
\item Theme 
\item Unit 1 
\item Unit 2 
\item Unit 3 
\item ..... 
\item Unit n
}}
}
\value{
A symmetric matrix where the cell in row i, col j is the Levenshtein Similarity Index between line i and line j from the file filename
}
\description{
Compute the Levenshtein Similarity Index for each possible pair of strings in the input file
}
