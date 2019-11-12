\name{readGFBed}
\alias{readGFBed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read genomic feature bed file
}
\description{
Function to import a genomic feature bed file. The bed file must contain coordinates of the genomic features. It can be for instance ChIP-seq peak coordinates or functional element coordinates.
}
\usage{
readGFBed(GFBedFile, seqInfoChr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GFBedFile}{
The name of the bed file. If it does not contain an absolute path, the file name is relative to the current working directory.
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
bed_file <- system.file("inst/extdata", "BEAF32_Kc167_CS_dm3.bed", package="HiCfeat")
data(dataExample)
seqInfoChr=dataExample$SeqInfoChr
GFBed=readGFBed(bed_file,seqInfoChr)
}
