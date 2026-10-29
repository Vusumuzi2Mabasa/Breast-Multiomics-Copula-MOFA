# Independently recompute AUC, average precision, Brier and subspace recovery.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else "."
fixtures <- jsonlite::fromJSON(file.path(root,"results/validation/R_fixtures.json"),simplifyVector=FALSE)
results <- lapply(fixtures,function(f){
 y <- unlist(f$truth); p <- unlist(f$probability)
 n1 <- sum(y==1); n0 <- sum(y==0)
 auc <- (sum(rank(p,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0)
 ord <- order(p,decreasing=TRUE); ys <- y[ord]
 ap <- sum(cumsum(ys)/seq_along(ys)*ys)/n1
 brier <- mean((y-p)^2)
 scores <- do.call(rbind,f$scores); latent <- do.call(rbind,f$latent)
 qa <- qr.Q(qr(scale(scores,scale=FALSE))); qb <- qr.Q(qr(scale(latent,scale=FALSE)))
 recovery <- mean(svd(crossprod(qa,qb),nu=0,nv=0)$d^2)
 err <- max(abs(c(auc,ap,brier,recovery)-unlist(f$expected)))
 if(err>1e-10) stop(f$name,": independent metric mismatch ",err)
 list(fixture=f$name,max_absolute_error=err)
})
jsonlite::write_json(list(status="passed",fixtures=length(results),checks=results),file.path(root,"results/validation/R_metrics.json"),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat("Independent R metric checks passed for",length(results),"fixtures.\n")
