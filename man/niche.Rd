% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topological.R
\name{niche}
\alias{niche}
\title{Niche Model Food Web}
\usage{
niche(S, C)
}
\arguments{
\item{S}{Number of species in the community.}

\item{C}{The connectance, or fraction of realized links in the food web.}
}
\value{
An adjacency matrix for a niche model food web.
}
\description{
Niche Model Food Web
}
\section{Reference}{

Williams, R. J., and N. D. Martinez. 2000. Simple rules yield complex food webs. Nature 404:180–183.
}

\examples{
niche(20, .1)
}
