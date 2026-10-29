# Independent bounded optimisation and efficient-score arithmetic.
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))args[1] else "."
f<-jsonlite::fromJSON(file.path(root,"v2/results/validation/R_conditional_fixtures.json"),simplifyVector=FALSE)
checks<-list()
for(item in f){
 u<-unlist(item$u);v<-unlist(item$v);t<-unlist(item$t);x<-do.call(rbind,lapply(item$x,unlist))
 objective<-function(theta){r<-tanh(drop(x%*%theta));vv<-1-r*r;mean(.5*log(vv)+.5*(u*u+v*v-2*r*u*v)/vv)}
 initial<-rep(0,ncol(x));initial[1]<-atanh(max(-.8,min(.8,cor(u,v))))
 # Numerical finite-difference derivatives intentionally avoid the Python gradient.
 fit<-optim(initial,objective,method="L-BFGS-B",lower=rep(-2,ncol(x)),upper=rep(2,ncol(x)),control=list(maxit=500,factr=1e4,pgtol=1e-7,ndeps=rep(1e-5,ncol(x))))
 rho<-tanh(drop(x%*%fit$par));w<-1+rho*rho
 s<-rho+((1+rho*rho)*u*v-rho*(u*u+v*v))/(1-rho*rho)
 coef<-qr.solve(x*sqrt(w),t*sqrt(w));eff<-t-drop(x%*%coef)
 score<-sum(eff*s)/sqrt(sum((eff*s)^2)+1e-12)
 error<-abs(score-item$expected)
 if(error>0.01)stop("Conditional score mismatch: ",item$scenario," ",error)
 checks[[length(checks)+1]]<-list(scenario=item$scenario,n=item$n,gene=item$gene,score=score,max_absolute_error=error,optim_code=fit$convergence)
}
jsonlite::write_json(list(status="passed",fixtures=length(checks),checks=checks,R_version=as.character(getRversion())),file.path(root,"v2/results/validation/R_conditional.json"),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat("Independent R conditional score checks passed for",length(checks),"fixtures.\n")
