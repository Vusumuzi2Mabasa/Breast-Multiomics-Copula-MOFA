# Official mixOmics DIABLO, kept separate from unsupervised integration.
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))args[1] else "."
pilot<-"--pilot" %in% args
suppressPackageStartupMessages(library(mixOmics));library(jsonlite)
readmatrix<-function(p)as.matrix(read.csv(p,header=FALSE,check.names=FALSE))
ba<-function(y,p)mean(c(mean(p[y==0]==0),mean(p[y==1]==1)))
prepare<-function(a,b,xa,xb){
  oa<-ob<-list();x<-cbind(1,xa);xt<-cbind(1,xb)
  for(m in names(a)){
    aa<-a[[m]];bb<-b[[m]]
    if(m %in% c("RNA","protein"))for(j in seq_len(ncol(aa))){
      good<-is.finite(aa[,j]);coef<-lm.fit(x[good,,drop=FALSE],aa[good,j])$coefficients;coef[!is.finite(coef)]<-0
      aa[,j]<-aa[,j]-x%*%coef;bb[,j]<-bb[,j]-xt%*%coef
    }
    for(j in seq_len(ncol(aa))){mu<-mean(aa[,j],na.rm=TRUE);aa[!is.finite(aa[,j]),j]<-mu;bb[!is.finite(bb[,j]),j]<-mu}
    good<-apply(aa,2,sd)>1e-7
    aa<-aa[,good,drop=FALSE];bb<-bb[,good,drop=FALSE]
    if(ncol(aa)<2)stop("Insufficient variable features")
    oa[[m]]<-aa;ob[[m]]<-bb
  }
  list(train=oa,test=ob)
}
fitmodel<-function(x,y,k,keep){
  design<-matrix(.1,length(x),length(x));diag(design)<-0
  keepx<-lapply(x,function(a)rep(min(keep,ncol(a)),k))
  withCallingHandlers(block.splsda(x,factor(y,levels=c(0,1)),ncomp=k,keepX=keepx,design=design,scale=TRUE),warning=function(w){currentWarnings<<-c(currentWarnings,conditionMessage(w));invokeRestart("muffleWarning")})
}
predictclass<-function(fit,x,k){
  result<-predict(fit,newdata=x)$WeightedVote$centroids.dist
  as.integer(result[,k])
}
folder<-file.path(root,if(pilot)"v2/work/integration_pilot" else "v2/results/integration")
outdir<-file.path(root,if(pilot)"v2/work/diablo_pilot" else "v2/results/diablo");dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
files<-list.files(folder,pattern="[.]json$",full.names=TRUE)
for(file in files){
  id<-sub("[.]json$","",basename(file));out<-file.path(outdir,paste0(id,".json"));if(file.exists(out))next
  meta<-fromJSON(file,simplifyVector=FALSE);set.seed(meta$seed+181)
  input<-file.path(root,"v2/work/diablo_inputs",id)
  blocks<-c("germline","somatic","RNA","protein")
  readblocks<-function(label){
    d<-lapply(blocks,function(m){a<-readmatrix(file.path(input,paste0(label,"_",m,".csv")));colnames(a)<-paste0(m,"_",seq_len(ncol(a)));a});names(d)<-blocks;d
  }
  a<-readblocks("train");b<-readblocks("test");u<-readblocks("unseen")
  xa<-readmatrix(file.path(input,"train_covariates.csv"));xb<-readmatrix(file.path(input,"test_covariates.csv"));xu<-readmatrix(file.path(input,"unseen_covariates.csv"))
  colnames(xa)<-colnames(xb)<-colnames(xu)<-paste0("clinical_",seq_len(ncol(xa)))
  a$clinical<-xa;b$clinical<-xb;u$clinical<-xu
  ya<-scan(file.path(input,"train_y.csv"),quiet=TRUE);yb<-scan(file.path(input,"test_y.csv"),quiet=TRUE);yu<-scan(file.path(input,"unseen_y.csv"),quiet=TRUE)
  began<-proc.time()[3]
  currentWarnings<-character()
  result<-tryCatch({
    fold<-integer(length(ya));for(g in c(0,1)){idx<-sample(which(ya==g));fold[idx]<-rep(1:3,length.out=length(idx))}
    candidates<-list()
    for(k in c(1,2))for(keep in c(4,8)){
      score<-numeric(3)
      for(f in 1:3){
        tr<-fold!=f;te<-!tr;p<-prepare(lapply(a,function(x)x[tr,,drop=FALSE]),lapply(a,function(x)x[te,,drop=FALSE]),xa[tr,,drop=FALSE],xa[te,,drop=FALSE])
        fit<-fitmodel(p$train,ya[tr],k,keep);score[f]<-ba(ya[te],predictclass(fit,p$test,k))
      }
      candidates[[length(candidates)+1]]<-list(k=k,keep=keep,validation_balanced_accuracy=mean(score))
    }
    best<-candidates[[which.max(sapply(candidates,function(x)x$validation_balanced_accuracy))]]
    p<-prepare(a,b,xa,xb);pu<-prepare(a,u,xa,xu)
    fit<-fitmodel(p$train,ya,best$k,best$keep)
    pred<-predictclass(fit,p$test,best$k);predu<-predictclass(fit,pu$test,best$k)
    list(status="ok",method="DIABLO supervised",scenario=meta$scenario,n=meta$n,replicate=meta$replicate,seed=meta$seed,k=best$k,keep=best$keep,candidates=candidates,
         balanced_accuracy=ba(yb,pred),unseen_balanced_accuracy=ba(yu,predu),test_class=pred,unseen_class=predu,
         warnings=unique(currentWarnings),diagnostic_pass=!any(grepl("converg",currentWarnings,ignore.case=TRUE)),software=list(R=as.character(getRversion()),mixOmics=as.character(packageVersion("mixOmics"))),seconds=unname(proc.time()[3]-began))
  },error=function(e)list(status="failed",scenario=meta$scenario,n=meta$n,replicate=meta$replicate,seed=meta$seed,error=conditionMessage(e)))
  write_json(result,out,auto_unbox=TRUE,pretty=TRUE,digits=16)
  cat(id,result$status,"\n")
}
