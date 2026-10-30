# Official supervised mixOmics DIABLO; two molecular blocks plus clinical adjustment inputs.
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))args[1] else "."
pilot<-'--pilot' %in% args
suppressPackageStartupMessages(library(mixOmics));library(jsonlite)
readmatrix<-function(p)as.matrix(read.csv(p,header=FALSE,check.names=FALSE))
ba<-function(y,p)mean(c(mean(p[y==0]==0),mean(p[y==1]==1)))
prepare<-function(a,b,ca,cb){
  xa<-cbind(1,ca[,3:6,drop=FALSE]);xb<-cbind(1,cb[,3:6,drop=FALSE]);oa<-ob<-list()
  for(name in names(a)){
    coef<-lm.fit(xa,a[[name]])$coefficients;coef[!is.finite(coef)]<-0
    aa<-a[[name]]-xa%*%coef;bb<-b[[name]]-xb%*%coef
    good<-apply(aa,2,sd)>1e-7
    oa[[name]]<-aa[,good,drop=FALSE];ob[[name]]<-bb[,good,drop=FALSE]
  }
  # Same observed clinical variables supplied to MOFA's downstream ridge model.
  oa$clinical<-ca;ob$clinical<-cb
  colnames(oa$clinical)<-colnames(ob$clinical)<-paste0('clinical_',1:ncol(ca))
  list(train=oa,test=ob)
}
fitmodel<-function(x,y,k,keep){
  design<-matrix(.1,length(x),length(x));diag(design)<-0
  keepx<-lapply(x,function(a)rep(min(keep,ncol(a)),k))
  withCallingHandlers(block.splsda(x,factor(y,levels=c(0,1)),ncomp=k,keepX=keepx,design=design,scale=TRUE,max.iter=1000),warning=function(w){currentWarnings<<-c(currentWarnings,conditionMessage(w));invokeRestart('muffleWarning')})
}
predclass<-function(model,x,k)as.integer(predict(model,newdata=x)$WeightedVote$centroids.dist[,k])
files<-list.files(file.path(root,if(pilot)'v3/work/pilot_integration' else 'v3/results/integration'),pattern='[.]json$',full.names=TRUE)
for(file in files){
 id<-sub('[.]json$','',basename(file));out<-file.path(root,if(pilot)'v3/work/pilot_diablo' else 'v3/results/diablo',paste0(id,'.json'));dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE);if(file.exists(out))next
 meta<-fromJSON(file,simplifyVector=FALSE);if(meta$status!='ok')next
 set.seed(meta$seed+181);input<-file.path(root,if(pilot)'v3/work/pilot_diablo_inputs' else 'v3/work/diablo_inputs',id);began<-proc.time()[3];currentWarnings<-character()
 result<-tryCatch({
   readblocks<-function(label){x<-lapply(c('genotype','RNA'),function(m){a<-readmatrix(file.path(input,paste0(label,'_',m,'.csv')));colnames(a)<-paste0(m,'_',1:ncol(a));a});names(x)<-c('genotype','RNA');x}
   a<-readblocks('train');b<-readblocks('test');ca<-readmatrix(file.path(input,'train_covariates.csv'));cb<-readmatrix(file.path(input,'test_covariates.csv'))
   ya<-scan(file.path(input,'train_y.csv'),quiet=TRUE);yb<-scan(file.path(input,'test_y.csv'),quiet=TRUE)
   fold<-integer(length(ya));for(group in c(0,1)){idx<-sample(which(ya==group));fold[idx]<-rep(1:3,length.out=length(idx))}
   candidates<-list()
   for(k in c(1,2))for(keep in c(4,8)){
     score<-numeric(3)
     for(f in 1:3){tr<-fold!=f;te<-!tr;p<-prepare(lapply(a,function(x)x[tr,,drop=FALSE]),lapply(a,function(x)x[te,,drop=FALSE]),ca[tr,,drop=FALSE],ca[te,,drop=FALSE]);model<-fitmodel(p$train,ya[tr],k,keep);score[f]<-ba(ya[te],predclass(model,p$test,k))}
     candidates[[length(candidates)+1]]<-list(k=k,keep=keep,validation_balanced_accuracy=mean(score))
   }
   best<-candidates[[which.max(sapply(candidates,function(x)x$validation_balanced_accuracy))]]
   p<-prepare(a,b,ca,cb);model<-fitmodel(p$train,ya,best$k,best$keep);pred<-predclass(model,p$test,best$k)
   selected<-lapply(model$loadings,function(x)which(rowSums(abs(x))>0));names(selected)<-names(model$loadings)
   list(status='ok',method='DIABLO',scenario=meta$scenario,n=meta$n,replicate=meta$replicate,seed=meta$seed,k=best$k,keep=best$keep,candidates=candidates,balanced_accuracy=ba(yb,pred),test_class=pred,test_y=yb,selected=selected,train_scores=lapply(model$variates,unname),warnings=unique(currentWarnings),diagnostic_pass=!any(grepl('converg',currentWarnings,ignore.case=TRUE)),R=as.character(getRversion()),mixOmics=as.character(packageVersion('mixOmics')),seconds=unname(proc.time()[3]-began))
 },error=function(e)list(status='failed',scenario=meta$scenario,n=meta$n,replicate=meta$replicate,seed=meta$seed,error=conditionMessage(e)))
 write_json(result,out,auto_unbox=TRUE,pretty=TRUE,digits=16);cat(id,result$status,'\n')
}
