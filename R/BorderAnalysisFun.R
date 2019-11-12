# Raphael Mourad
# Cuvier team, LBME lab
# 02/04/2015





# MAIN FUNCTION ----------------------------------------------------------
borderAnalysisFun<-function(genomicFeatureList.GR,GFDataType,annotNames,domains.GR,seqInfoChr,analysisMode,binSize=50,borderSize=1000,LRT=FALSE,interactionTerms="",verbose=FALSE){


# CHECK INPUT DATA ------------------------------------------------------------ 

if(verbose){print("DATA CHECKING")}

if(class(genomicFeatureList.GR)!="list"){print("genomicFeatureList.GR is not a list object!"); return(0)}
for(i in 1:length(genomicFeatureList.GR)){
 if(class(genomicFeatureList.GR[[i]])!="GRanges"){print("i-th object of genomicFeatureList.GR is not a GenomicRanges object!"); return(0)}
}
if(class(GFDataType)!="character"){print("GFDataType is not a character object!"); return(0)}
if(class(annotNames)!="character"){print("annotNames is not a character object!"); return(0)}
if(class(domains.GR)!="GRanges"){print("domains.GR is not a GenomicRanges object!"); return(0)}
if(class(seqInfoChr)!="Seqinfo"){print("seqInfoChr is not a seqinfo object!");return(0)}
if(class(binSize)!="integer" & class(binSize)!="numeric"){print("binSize is not an integer or numeric object!");return(0)}
if(class(borderSize)!="integer" & class(borderSize)!="numeric"){print("borderSize is not an integer or numeric object!");return(0)}
if(class(LRT)!="logical"){print("LRT is not a logical object!"); return(0)}
if(class(analysisMode)!="character"){print("analysisMode is not a character object!"); return(0)}
if(class(interactionTerms)!="character"){print("interactionTerms is not a character object!"); return(0)}



# PROCESS DATA ------------------------------------------------------------ 

if(verbose){print("DATA PREPROCESSING")}

chr.V=as.character(seqnames(seqInfoChr))

Borders.GR=NULL
for(i in 1:length(chr.V)){
 if(sum(seqnames(domains.GR)==chr.V[i])){
  domains.GRi=domains.GR[seqnames(domains.GR)==chr.V[i]]

  Borders.GRi=GRanges(seqnames=seqnames(domains.GRi[-1]),IRanges(start=start(domains.GRi[-1])-borderSize,end=start(domains.GRi[-1])+borderSize-1),seqinfo=seqInfoChr)
  if(i==1){
   Borders.GR=Borders.GRi
  }else{
   Borders.GR=c(Borders.GR,Borders.GRi)
  }
 }else{print(paste0("No ",chr.V[i]," in Domains.GR"))}
}

# Binned matrix
binMat=NULL
seqLengthChr=seqlengths(seqInfoChr)
for(i in 1:length(chr.V)){
 BordStarti=seq(1,seqLengthChr[i],by=binSize)
 BordEndi=BordStarti+binSize-1
 BordEndi[length(BordEndi)]=seqLengthChr[i]
 binMat=rbind(binMat,cbind(chr.V[i],BordStarti,BordEndi))
 if(verbose){print(chr.V[i])}
}
binMat.GR=GRanges(seqnames=binMat[,1],IRanges(start=as.numeric(binMat[,2]),end=as.numeric(binMat[,3])))
seqinfo(binMat.GR)=seqInfoChr
olBinBorders=findOverlaps(binMat.GR,Borders.GR)
binMat.GR$Border=rep(0,length(binMat.GR))
#binMat.GR$Border[queryHits(olBinBorders)]=1
binMat.GR$Border[as.data.frame(olBinBorders)[[1]]]=1

# Annotate borders
binMat.Mat=NULL
for(i in 1:length(genomicFeatureList.GR)){
 if(GFDataType=="bed"){
  genomicFeatureList.GR[[i]]$score=rep(1,length(genomicFeatureList.GR[[i]]))
 }

 binMati=averagePerBin(genomicFeatureList.GR[[i]],binSize,"score")
 binMat.Mat=as(cbind(binMat.Mat,binMati),"Matrix")
 if(verbose){print(paste0(annotNames[i]," annotated"))}
}
binMat.Mat=cbind(binMat.GR$Border,binMat.Mat)
colnames(binMat.Mat)<-c("Border",annotNames)

# Compute correlations among genomic features
corGF=cor(as.matrix(binMat.Mat[,-1]))


# ENRICHMENT TEST ------------------------------------------------------------ 

if(verbose){print("DATA ANALYSIS")}

matCoefMargGLM=NULL
matCoefMultGLM=NULL
matCoefInterGLM=NULL
matCoefMultLasso=NULL
matCoefInterLasso=NULL
MultGLM=NULL
InterGLM=NULL

# Binned matrix 
binMat.Mati=binMat.Mat
if(sum(analysisMode%in%c("EnrichmentTest","MLR","MLRInter","MLRInterLasso"))){
 binMat.mati=as.data.frame(as.matrix(binMat.Mati))
}

# Enrichment test
if(sum(analysisMode=="EnrichmentTest")){
 if(verbose){print("Enrichment Test")}
 matCoefMargi=NULL
 pval_LRT=rep(NA,length(annotNames))
 AnalysisMargZero=glm(Border~1,data=binMat.mati,family=binomial())
 for(j in 1:length(annotNames)){
  formGLMMarg=as.formula(paste0("Border~",annotNames[j]))
  AnalysisMarg=glm(formGLMMarg,data=binMat.mati,family=binomial())
  AnalysisMargRes=summary(AnalysisMarg)
  coefj=AnalysisMargRes$coefficients[2,]
  if(LRT){
  #lrtj=lrtest(AnalysisMarg,AnalysisMargZero)
  Dj=(logLik(AnalysisMargZero)[1]-logLik(AnalysisMarg)[1])*-2
  pval_LRT[j]=1-pchisq(Dj,1)
  }
  matCoefMargi=rbind(matCoefMargi,coefj)
  if(verbose){print(annotNames[j])}
 }
 rownames(matCoefMargi)=annotNames
 freqBins=colSums(as.matrix(binMat.mati[,-1]))
 freqPeaks=sapply(genomicFeatureList.GR,length)
 matCoefMargGLM=data.frame(GenomicFeature=annotNames,matCoefMargi,pval_LRT,freqBins,freqPeaks)
}

# Multiple Logistic Regression
if(sum(analysisMode=="MLR")){
 if(verbose){print("Multiple Logistic Regression")}
 Analysisi=glm(Border~.,data=binMat.mati,family=binomial())
 MultGLM=summary(Analysisi)

 pval_LRT=rep(NA,length(annotNames))
 if(LRT){
  for(j in 1:length(annotNames)){
   formGLMMultj=as.formula(paste0("Border~",paste(annotNames[-j],collapse="+")))
   Analysisij=glm(formGLMMultj,data=binMat.mati,family=binomial())
   #lrtj=lrtest(Analysisij,Analysisi)
   Dj=(logLik(Analysisij)[1]-logLik(Analysisi)[1])*-2
   pval_LRT[j]=1-pchisq(Dj,1)
   if(verbose){print(paste0("LRT multi: ",annotNames[j]))}
  }
 }
 freqBins=colSums(binMat.mati[,-1])
 freqPeaks=sapply(genomicFeatureList.GR,length)
 matCoefMultGLM=data.frame(GenomicFeature=rownames(summary(Analysisi)$coefficients[-1,]),summary(Analysisi)$coefficients[-1,],pval_LRT,freqBins,freqPeaks)
}

# Multiple Logistic Regression Estimated by Lasso
if(sum(analysisMode=="MLRLasso")){
 if(verbose){print("Multiple Logistic Regression with Lasso Estimation")}
 CVLasso=cv.glmnet(binMat.Mati[,-(1:2)],binMat.Mati[,1],family="binomial")
 lambda=CVLasso$lambda.min
 CVLassoError=CVLasso$cvm[which(CVLasso$lambda==lambda)]
 coefLasso=CVLasso$glmnet.fit$beta[,which(CVLasso$lambda==lambda)]
 matCoefMultLasso=data.frame(GenomicFeature=names(coefLasso),Estimate=round(coefLasso,5))
}

# Multiple Logistic Regression with Interaction Terms 
if(sum(analysisMode=="MLRInter")){
 if(verbose){print("Multiple Logistic Regression with Interaction Terms")}
 OneWayTerms=paste(annotNames,collapse='+')
 formInterTesti=as.formula(paste0("Border~",paste(OneWayTerms,collapse="+"),"+",paste(interactionTerms,collapse="+")))
 enrichInterTesti=glm(formInterTesti,data=binMat.mati,family=binomial())
 InterGLM=summary(enrichInterTesti)

 pval_LRTInter=rep(NA,length(annotNames)+length(interactionTerms))
 if(LRT){
  for(j in 1:length(interactionTerms)){
   if(length(interactionTerms)>1){
    formGLMInterj=as.formula(paste0("Border~",paste(OneWayTerms,collapse="+"),"+",paste(interactionTerms[-j],collapse="+")))
   }else{
    formGLMInterj=as.formula(paste0("Border~",paste(OneWayTerms,collapse="+")))
   }
   enrichInterTestij=glm(formGLMInterj,data=binMat.mati,family=binomial())
   #lrtInterj=lrtest(enrichInterTestij,enrichInterTesti)
   DInterj=(logLik(enrichInterTestij)[1]-logLik(enrichInterTesti)[1])*-2
   pval_LRTInter[length(annotNames)+j]=1-pchisq(DInterj,1)
   rm(enrichInterTestij)
   if(verbose){print(paste0("LRT multi: ",interactionTerms[j]))}
  }
 }
 matCoefInterGLM=data.frame(GenomicFeature=rownames(summary(enrichInterTesti)$coefficients[-1,]),summary(enrichInterTesti)$coefficients[-1,],pval_LRT=pval_LRTInter)
}

# Multiple Logistic Regression with Interaction Terms Estimated by Lasso
if(sum(analysisMode=="MLRInterLasso")){
 if(verbose){print("Multiple Logistic Regression with Interaction Terms Estimated by Lasso")}
 OneWayTerms=paste(annotNames,collapse='+')
 formInterTesti=as.formula(paste0("Border~",paste(OneWayTerms,collapse="+"),"+",paste(interactionTerms,collapse="+")))
 binMatInter.Mati=sparse.model.matrix(formInterTesti,data=binMat.mati)

 CVLasso=cv.glmnet(binMatInter.Mati[,-1],binMat.Mati[,1],family="binomial")
 lambda=CVLasso$lambda.min
 CVLassoError=CVLasso$cvm[which(CVLasso$lambda==lambda)]
 coefInterLasso=CVLasso$glmnet.fit$beta[,which(CVLasso$lambda==lambda)]
 matCoefInterLasso=data.frame(GenomicFeature=names(coefInterLasso),round(coefInterLasso,5))
}




if(verbose){print("All analyses done")}


list2return=list(Enrich=matCoefMargGLM, MLR=matCoefMultGLM, MLRLasso=matCoefMultLasso,MLRInter=matCoefInterGLM,MLRInterLasso=matCoefInterLasso, Mat=binMat.Mat, MLRGLM=MultGLM, MLRInterGLM=InterGLM,CorGF=corGF)

return(list2return)
}

