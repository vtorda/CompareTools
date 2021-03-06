% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/retrieve_event3.R
\name{retrieve_event3}
\alias{retrieve_event3}
\title{Retrieve event3}
\usage{
retrieve_event3(compare, new_object = FALSE)
}
\arguments{
\item{compare}{a compare object}

\item{partitioning}{logical, does the clusters handled separately (TRUE) or all together (FALSE).
If TRUE and no groups is given, all clusters will be handeled separately.
If a group file is given clusteres will be grouped according to that file.}

\item{groups}{a Tab Separated file containing the cluster names and their position in predefined groups}

\item{part_name}{}
}
\value{
this function returns a compare object with additional event matrices
}
\description{
Retrieve event3
}
