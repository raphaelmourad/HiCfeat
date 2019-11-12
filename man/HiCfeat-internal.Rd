\name{averagePerBin}
\alias{averagePerBin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute average value of a genomic feature for each bin.
}
\description{
This function is for internal use only! 
The function computes the average value of a genomic feature for each bin of size binsize. 
}
\usage{
averagePerBin(x, binsize, mcolname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The GRanges object used to compute the average for each bin.
}
  \item{binsize}{
Size of each bin. 
}
  \item{mcolname}{
The metadata column name used to compute the average for each bin. 
}
}
\value{
A vector object is returned.
}
\author{
Raphael Mourad
}
