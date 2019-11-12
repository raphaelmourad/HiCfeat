\name{dataExample}
\alias{dataExample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read chromatin domain bed file
}
\description{
This dataExample object comprises three objects: GenomicFeatureList.GR, AnnotNames and Domains.GR. GenomicFeatureList.GR is a list of GRanges objects, one for each genomic feature to analyze. AnnotNames is a vector of the genomic feature names. Domains.GR is a GRanges object for topologically associating domains.

Topologically associating domains (object Domains.GR) were computed from Hi-C experiments on Drosophila melanogaster Kc167 cells from Gene Expression Omnibus (GEO) accession GSE62904. ChIP-seq data (object GenomicFeatureList.GR) on Drosophila melanogaster Kc167 cells were downloaded from Gene Expression Omnibus (GEO) accessions GSE30740, GSE42085 and GSE54529. Only data for chromosome 2L are provided.
}
\usage{
data(dataExample)
}
\references{
Kevin Van Bortle, Michael H. Nichols, Li Li, Chin-Tong Ong, Naomi Takenaka, Zhaohui S. Qin, and Victor G. Corces. Insulator function and topological domain border strength scale with architectural protein occupancy. Genome Biology, 15(5):R82+, June 2014.

Li Li, Xiaowen Lyu, Chunhui Hou, Naomi Takenaka, Huy Q. Nguyen, Chin-Tong Ong, Caelin Cubenas-Potts, Ming Hu, Elissa P. Lei, Giovanni Bosco, Zhaohui S. Qin and Victor G. Corces. Widespread rearrangement of 3D chromatin organization underlies Polycomb-mediated stress-induced silencing. Molecular Cell, 15:S1097-2765, March 2015.
}
\author{
Raphael Mourad
}
