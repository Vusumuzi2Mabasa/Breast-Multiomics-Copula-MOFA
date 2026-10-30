# Independent QR nested F arithmetic and held-out prediction checks.
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))args[1] else '.'
library(jsonlite)
fixtures<-fromJSON(file.path(root,'v3/results/R_eqtl_fixtures.json'),simplifyVector=FALSE)
matrixfrom<-function(x)do.call(rbind,lapply(x,unlist))
ftest<-function(x0,x1,y){
 a<-lm.fit(x0,y);b<-lm.fit(x1,y);r0<-sum(a$residuals^2);r1<-sum(b$residuals^2)
 stat<-max(0,(r0-r1)/(b$rank-a$rank))/(r1/(length(y)-b$rank))
 pf(stat,b$rank-a$rank,length(y)-b$rank,lower.tail=FALSE)
}
errors<-numeric();count<-0
for(f in fixtures){
 g<-matrixfrom(f$genotype);rna<-matrixfrom(f$RNA);c<-matrixfrom(f$covariates);y<-unlist(f$y);base<-cbind(1,c,y);pa<-pi<-numeric(12)
 for(j in 1:12){
  idx<-((j-1)%/%2)*4+1:4;gg<-g[,idx,drop=FALSE];common<-cbind(base,gg);full<-cbind(common,gg*y)
  pa[j]<-ftest(base,full,rna[,j]);pi[j]<-ftest(common,full,rna[,j]);e<-f$expected[[j]]
  errors<-c(errors,abs(pa[j]-e$omnibus_p),abs(pi[j]-e$interaction_p));count<-count+1
 }
 qa<-p.adjust(pa,'BH');qi<-p.adjust(pi,'BH')
 for(j in 1:12)errors<-c(errors,abs(qa[j]-f$expected[[j]]$omnibus_q),abs(qi[j]-f$expected[[j]]$interaction_q))
}
metric<-fromJSON(file.path(root,'v3/results/R_metric_fixtures.json'),simplifyVector=FALSE);me<-numeric()
for(f in metric){
 y<-unlist(f$y);p<-unlist(f$p);n1<-sum(y==1);n0<-sum(y==0)
 auc<-(sum(rank(p)[y==1])-n1*(n1+1)/2)/(n1*n0)
 ba<-mean(c(mean((p[y==1]>=.5)),mean((p[y==0]<.5))))
 me<-c(me,abs(auc-f$auc),abs(ba-f$balanced_accuracy),abs(mean((p-y)^2)-f$brier))
}
stopifnot(max(errors)<1e-9,max(me)<1e-10)
write_json(list(status='passed',gene_model_fixtures=count,max_qtl_p_or_q_error=max(errors),metric_fixtures=length(metric),max_metric_error=max(me),R=as.character(getRversion())),file.path(root,'v3/results/R_validation.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat('Independent R checks passed\n')
