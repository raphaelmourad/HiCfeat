\name{readDomBed}
\alias{readDomBed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read chromatin domain bed file
}
\description{
Function to import a chromatin domain bed file. The bed file must contain coordinates of the chromatin domains. It can be for instance topologically associating domain (TAD) coordinates. In the bed file, the first column is the chromosome, the second column is the starting position of the domain, the third column is the ending position of the domain. The fourth column (optional) is the annotation of the domain (for instance: active, PcG, ...).
}
\usage{
readDomBed(domainBedFile, seqInfoChr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{domainBedFile}{
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
dom_file <- system.file("inst/extdata", "TAD_Corces_1kb_type_p1.bed", package="HiCfeat")
data(dataExample)
seqInfoChr=dataExample$SeqInfoChr
domBed=readDomBed(dom_file,seqInfoChr)
}
