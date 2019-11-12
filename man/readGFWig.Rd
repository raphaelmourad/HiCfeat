\name{readGFWig}
\alias{readGFWig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read genomic feature bigwig file
}
\description{
Function to import a genomic feature wig or bigWig file. The file can contain for instance normalized ChIP-seq values (log2(ChIP/Input) values).
}
\usage{
readGFWig(GFWigFile, seqInfoChr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GFWigFile}{
The name of the wig or bigWig file. If it does not contain an absolute path, the file name is relative to the current working directory.
}
  \item{seqInfoChr}{
A Seqinfo object for the corresponding genome. 
}
}
\value{
A GRanges object is returned.
}
\author{
Raphael Mourad
}

\examples{
# Same use as for readGFBed(). An example wig or bigWig file is not provided in the package 
# because these files have large memory sizes.
}
