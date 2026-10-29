args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))args[1] else "."
f<-jsonlite::fromJSON(file.path(root,"v2/results/validation/R_integration_fixtures.json"),simplifyVector=FALSE)
auc<-function(y,p){n1<-sum(y==1);n0<-sum(y==0);(sum(rank(p)[y==1])-n1*(n1+1)/2)/(n1*n0)}
checks<-lapply(f,function(a){
 y<-unlist(a$truth);p<-unlist(a$p);yu<-unlist(a$unseen_truth);pu<-unlist(a$unseen_p)
 ba<-mean(c(mean(p[y==1]>=.5),mean(p[y==0]<.5)))
 actual<-c(auc(y,p),mean((p-y)^2),auc(yu,pu),ba);error<-max(abs(actual-unlist(a$expected)))
 if(error>1e-12)stop("Metric mismatch ",a$name)
 list(fixture=a$name,max_error=error)
})
jsonlite::write_json(list(status="passed",fixtures=length(checks),checks=checks),file.path(root,"v2/results/validation/R_integration.json"),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat("Independent R integration metrics passed for",length(checks),"fixtures.\n")
