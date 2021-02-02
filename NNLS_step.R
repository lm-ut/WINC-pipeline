##This is a preliminary version

##LOAD FILES##

referenceTableName="CASmatrix.txt"
admixedTable <- "parfile_MXL.csv"

#functions from GLOBETROTTER
#GLOBETROTTER website: https://people.maths.bris.ac.uk/~madjl/finestructure/globetrotter.html
#GLOBETROTTER citation: Hellenthal, G. et al (2014), A Genetic Atlas of Human Admixture History, Science 343:747-751.

 getfit=function(predmat,fitdata,restrict=1){
   
   temp=matrix(predmat[-restrict,],ncol=dim(predmat)[2])
   for(i in 1:nrow(temp)) temp[i,]=temp[i,]-predmat[restrict,]
   
   fitdata2=fitdata-predmat[restrict,]
   
   v=nnls(t(temp),fitdata2)
   x=v$x
   newx=1:nrow(predmat)
   newx[!((1:nrow(predmat))==restrict)]=x
   newx[restrict]=1-sum(x)
   v$x=newx
   names(v$x)=rownames(predmat)
   
   return(v)
 }
 
 getoverallfit=function(predmat,fitdata){
   restrict=1
   rep=1
   i=1
   while(rep==1){
     q=getfit(predmat,fitdata,restrict=i)
     
     if(q$x[i]>0) rep=0
     i=i+1
   }
   
   return(q)
 }
 
 ##NNLS.mat
 nnls.mat=function(donors,recipients){
   tmp=t(sapply(1:nrow(recipients),function(X) getoverallfit(predmat=donors,fitdata=recipients[X,])$x))
   rownames(tmp)=rownames(recipients)
   return(tmp)
 }
 

putna <- function(winner,max,thres){
  tmpwin=winner
  tmpwin[max<thres]=NA
  return(tmpwin)
}


referenceTable=read.table(referenceTableName,header = T)

paramTable <- read.table(admixedTable,header = T,stringsAsFactors = F)

threemode="S3" %in% colnames(paramTable)

for (par in 1:nrow(paramTable)){
  S1n=paramTable[par,"S1"]
  S2n=paramTable[par,"S2"]
  ADMn=paramTable[par,"admixedName"]
  donorDir <- paramTable[par,"donorDir"]
  recipientDir <- paramTable[par,"pathPainting"]
  listsamplesFile <- paramTable[par,"listSamplesFile"]
  suffix<- paramTable[par,"suffixWindowsPainted"]
  
  
  listsamples=read.table(listsamplesFile)
  
  subsetS1 <- listsamples[listsamples[,1] %in% c(S1n),2]
  subsetS2 <- listsamples[listsamples[,1] %in% c(S2n),2]
  availsamples=gsub(x=list.files(path = donorDir,pattern = "\\.out"),pattern = suffix,replacement = "")
  subsetS1 <- subsetS1[(subsetS1 %in% availsamples)]
  subsetS2 <- subsetS2[(subsetS2 %in% availsamples)]
  
  subsetADM <- listsamples[listsamples[,1] %in% c(ADMn),2]
  
  
  if (threemode==TRUE) {
    S3n=paramTable[par,"S3"]
    subsetS3 <- listsamples[listsamples[,1] %in% c(S3n),2]
    subsetS3 <- subsetS3[(subsetS3 %in% availsamples)]
    S3 <- lapply(subsetS3,function(X) read.table(paste0(donorDir,"/",X,suffix),header = T,check.names = F))
    winnames=S3[[1]][,2:3]
    colnames(winnames)=c("Start","End")
    S3 <- lapply(S3,function (X) X[,-c(1:3)])
    S3=as.matrix(Reduce("+", S3) )
    S3=S3[1:(nrow(S1)/2),]+S3[((nrow(S3)/2)+1):nrow(S3),]
    S3=t(apply(S3,1,function(X) X/(sum(X)+1e-36)))
    S3=rbind(S3,S3)
    print(sprintf("Analysing %s-%s-%s...",S1n,S2n,S3n))
  } else {
    print(sprintf("Analysing %s-%s...",S1n,S2n))
  }

  ##########
  #process S1
  
  S1 <- lapply(subsetS1,function(X) read.table(paste0(donorDir,"/",X,suffix),header = T,check.names = F))
  winnames=S1[[1]][,2:3]
  colnames(winnames)=c("Start","End")
  S1 <- lapply(S1,function (X) X[,-c(1:3)])
  S1=as.matrix(Reduce("+", S1) )
  S1=S1[1:(nrow(S1)/2),]+S1[((nrow(S1)/2)+1):nrow(S1),]
  S1=t(apply(S1,1,function(X) X/(sum(X)+1e-36)))
  S1=rbind(S1,S1)
  
  S2 <- lapply(subsetS2,function(X) read.table(paste0(donorDir,"/",X,suffix),header = T,check.names = F))
  S2 <- lapply(S2,function (X) X[,-c(1:3)])
  S2=as.matrix(Reduce("+", S2) )
  S2=S2[1:(nrow(S2)/2),]+S2[((nrow(S2)/2)+1):nrow(S2),]
  S2=t(apply(S2,1,function(X) X/(sum(X)+1e-24)))
  S2=rbind(S2,S2)
  ####
  #columns 2 take
  cols=!(apply(S1,2,function(X) (all(X==0))))
  
  #process Recipent
  availRec=gsub(x=list.files(path = recipientDir,pattern = "\\.out"),pattern = suffix,replacement = "")
  subsetADM <- subsetADM[(subsetADM %in% availRec)]
  ADM <- lapply(subsetADM,function(X) read.table(paste0(recipientDir,X,suffix),header = T,check.names = F))
  ADM <- lapply(ADM,function (X) X[,-c(1:3)])
  ADM=lapply(ADM,function (Y) (t(apply(Y,1,function(X) X/(sum(X)+1e-24)))))
  admcols=!apply(ADM[[1]],2,function(X) all(X==0))
  
  allwindows=list()
  allcorr=c()
  for (win in 1:nrow(S1)){
    allwin=t(sapply(ADM,function (X) X[win,]))
    if (all(allwin==0)){usable="N"} else {usable="Y"} 
    if (threemode == TRUE){
    sources=rbind(S1[win,],S2[win,],S3[win,])
    corr <- cor(cbind(sources[1,cols],sources[2,cols],sources[3,cols]))
    corr=max(corr[upper.tri(corr)])}
    else{
      sources=rbind(S1[win,],S2[win,])
      corr <- cor(sources[1,cols],sources[2,cols])}
    #tmp<- nnls.mat(donors = as.matrix(sources[,cols]),recipients = as.matrix(allwin[,cols]))
    tmp= nnls.mat(donors = as.matrix(sources[,cols]),recipients = as.matrix(allwin[,admcols]))
    ref.tmp=referenceTable[referenceTable$rho==ceiling(corr*10)/10,]
    thres2use=sapply(c(.8,.9,.95,.99),function(X) min(ref.tmp[ref.tmp$accuracy>=X,"threshold"]))
    thres2use=c(thres2use,c(seq(.55,.9,by=.05),seq(.91,.99,by=.01)))
    winner <- apply(tmp,1,which.max)-1
    max=apply(tmp,1,max)
    allthres=sapply(thres2use,function(X) putna(winner,max,X))
    tmp <- cbind(tmp,winner,allthres,corr,usable)
    if (threemode == TRUE){
    colnames(tmp)=c("S1","S2","S3","NoThreshold","accuracy.8","accuracy.9","accuracy.95","accuracy.99",
                    paste0("thres",c(seq(.55,.9,by=.05),seq(.91,.99,by=.01))),"rho","Usable")} 
    else{
    colnames(tmp)=c("S1","S2","NoThreshold","accuracy.8","accuracy.9","accuracy.95","accuracy.99",
                    paste0("thres",c(seq(.55,.9,by=.05),seq(.91,.99,by=.01))),"rho","Usable")}
    allwindows[[length(allwindows)+1]] <- as.data.frame(tmp)
  }  

  haps=split(allwindows,rep(c("A","B"),each=length(allwindows)/2))
  hapsA=haps[[1]]
  hapsB=haps[[2]]
  
  allinds=list()
  for (i in (1:length(hapsA))){
    allinds[[length(allinds)+1]] <- hapsA[[i]]
    allinds[[length(allinds)+1]] <- hapsB[[i]]
  }
  
  N=2*nrow(allinds[[1]])
  
  allinds <- (do.call(rbind,allinds))
  
  allinds=split(allinds,rep(1:N,length(hapsA)))
  x=1
  for (ind in subsetADM){
    write.table(file = paste0(ind,".winc.tmp"),x = cbind(winnames,allinds[[x]]),quote=F,row.names = F)
    x=x+1
  }

}

