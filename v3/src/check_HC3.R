args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))args[1] else '.'
library(jsonlite);fixtures<-fromJSON(file.path(root,'v3/results/R_HC3_fixtures.json'),simplifyVector=FALSE)
matrixfrom<-function(a)do.call(rbind,lapply(a,unlist));errors<-numeric();count<-0
for(d in fixtures){
 g<-matrixfrom(d$genotype);rna<-matrixfrom(d$RNA);c<-matrixfrom(d$covariates);y<-unlist(d$y);base<-cbind(1,c,y)
 for(j in 1:12){
  gg<-g[,((j-1)%/%2)*4+1:4,drop=FALSE];x<-cbind(base,gg,gg*y);fit<-lm.fit(x,rna[,j]);bread<-chol2inv(chol(crossprod(x)));h<-rowSums(qr.Q(qr(x))^2)
  adjusted<-x*as.vector(fit$residuals/(1-h));V<-bread%*%crossprod(adjusted)%*%bread
  for(family in c('omnibus','interaction')){
   idx<-if(family=='omnibus')(ncol(base)+1):ncol(x) else (ncol(base)+5):ncol(x)
   b<-fit$coefficients[idx];w<-as.numeric(crossprod(b,solve(V[idx,idx,drop=FALSE],b)))
   p<-pf(w/length(idx),length(idx),nrow(x)-fit$rank,lower.tail=FALSE)
   errors<-c(errors,abs(p-d$expected[[j]][[paste0(family,'_p')]]));count<-count+1
  }
 }
}
stopifnot(max(errors)<1e-9)
write_json(list(status='passed',hypotheses_checked=count,max_p_error=max(errors)),file.path(root,'v3/results/R_HC3_validation.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat('HC3 R validation passed\n')
