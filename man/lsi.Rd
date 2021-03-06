% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leven.R
\name{lsi}
\alias{lsi}
\title{Compute the Levenshtein Similary Index between two vectors}
\usage{
lsi(x, y, cost_matrix = NULL, fleven = leven)
}
\arguments{
\item{x}{The first vector}

\item{y}{The second vector}

\item{cost_matrix}{An optional cost matrix to specify the cost of specific substitutions.}

\item{fleven}{The function to use to compute the Levenshtein distance. Defaults to leven. Use fast_leven instead to use the native (fast) code}
}
\value{
The LSI of x and y
}
\description{
Compute the Levenshtein Similary Index between two vectors
}
\examples{
lsi(c("k", "i", "t", "t", "e", "n"), c("s", "i", "t", "t", "i", "n", "g"))
}
\seealso{
\link{cost_matrix_from_file}

\link{fast_leven}

\link{leven}
}
