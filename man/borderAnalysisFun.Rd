\name{borderAnalysisFun}
\alias{borderAnalysisFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate enrichment test for chromatin domain borders
}
\description{
This is the main function. This function estimates the influences of genomic features on chromatin domain borders by a multiple logistic regression model. Enrichment test by simple logistic regression is also included in this function for comparisons. Genomic features can be either coordinate data (e.g. ChIP-seq peak or functional element coordinates) or quantitative data (e.g. normalized ChIP-seq values: log2(ChIP/Input) values). The proposed method is flexible and can account for statistical interactions among multiple genomic features.
}
\usage{
borderAnalysisFun(genomicFeatureList.GR, GFDataType, annotNames, domains.GR, 
		  			seqInfoChr, analysisMode, binSize, borderSize, 
					LRT, interactionTerms, verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genomicFeatureList.GR}{
A list of GRanges objects. Each GRanges object has been built from either coordinate data using readGFBed function (for instance ChIP-seq peak coordinates or functional element coordinates), or quantitative data using readGFWig (e.g. normalized ChIP-seq values: log2(ChIP/Input) values). All GRanges object should be either coordinate data or quantitative data. There cannot be a mix of GRanges from coordinate data and quantitative data.
}
  \item{GFDataType}{
GFDataType = "bed" if all GRanges were built from readGFBed function, or GFDataType = "wig" if all GRanges were built from readGFWig function.
}
  \item{annotNames}{
A character vector defining the name of each genomic feature. 
Names should not comprise any special character such as ":/+-*^,;!?" because an internal R formula object is created inside the function borderAnalysisFun.
}
  \item{domains.GR}{
A GRanges for the chromatin domain data that was built from readDomBed function.
}
  \item{seqInfoChr}{
A Seqinfo object of the corresponding genome. 
}
  \item{analysisMode}{
A character vector analysisMode=c("EnrichmentTest", "MLR", "MLRLasso", "MLRInter", "MLRInterLasso"). 
"EnrichmentTest" for enrichment test by simple logistic regression, "MLR" for multiple logistic regression, "MLRLasso" for multiple logistic regression using lasso parameter estimation, and "MLRInter" for multiple logistic regression with interaction terms provided by the interactionTerms parameter, and "MLRInter" for multiple logistic regression with interaction terms using lasso parameter estimation.
}
  \item{binSize}{
The size of the genomic bins (in bp) used for logistic regression computation. By default, binSize = 50.
}
  \item{borderSize}{
The size x 2 of the window (in bp) surrounding the border between two domains. By default, borderSize = 1000.
}
  \item{LRT}{
If TRUE, a likelihood ratio test is used to assess the significance of beta parameters. By default, LRT = FALSE, and Wald's test is computed instead.
}
  \item{interactionTerms}{
A character vector InteractionTerms that contains the interaction terms between the genomic features. For instance, if there are two genomic features "F1" and "F2", the user can add the interaction term "F1:F2" in the same manner as in an lm or glm object.
}
  \item{verbose}{
If verbose = TRUE, then each step of the method is printed. By default, verbose = FALSE.
}
}
\value{
The object contains 8 attributes:
  \item{Enrich}{Enrichment test (simple logistic regression) result for each genomic feature.}
  \item{MLR}{Multiple logistic regression results for all genomic features.}
  \item{MLRLasso}{Multiple logistic regression results for all genomic features. Parameters of the logistic regression model are learned by lasso.}
  \item{MLRInter}{Multiple logistic regression results for all genomic features together with interaction terms that are provided in the character vector interactionTerms.}
  \item{MLRInterLasso}{Multiple logistic regression results for all genomic features together with interaction terms that are provided in the character vector interactionTerms. Parameters of the logistic regression model are learned by lasso.}
  \item{Mat}{Matrix used for all logistic regression computations.}
  \item{MultGLM}{The glm object from the multiple logistic regression model without interaction terms.}
  \item{InterGLM}{The glm object from the multiple logistic regression model with interaction terms.}
  \item{CorGF}{Correlation matrix between genomic features.}
}
\author{
Raphael Mourad
}

\examples{
# Load library HiCfeat
library(HiCfeat)

# Load data example.
data(dataExample) 

# Multiple logistic regression estimated by maximum likelihood
print("Multiple logistic regression estimated by maximum likelihood")
BA_res=borderAnalysisFun(genomicFeatureList.GR=dataExample$GenomicFeatureList.GR,
	GFDataType="bed",annotNames=dataExample$AnnotNames,domains.GR=dataExample$Domains.GR,
	seqInfoChr=dataExample$SeqInfoChr,analysisMode="MLR",binSize=500,borderSize=5e3,
	LRT=FALSE,interactionTerms="",verbose=FALSE)

}
