 mlfit <- function(b,p,x,offset,d,h,ds,y,maxiter,tol){

 
 for(nit in 1:maxiter){
	  b0 <- b
	  fd <- matrix(0,p,1)
	  sd <- matrix(0,p,p)
	  if(nit==1){
	  	ebx <- exp(x%*%b)*exp(offset)
	  	l0 <- sum(d*log(h+ebx)-ds-y*ebx)
	  }
	  for(it in 1:p){
		  fd[it,1] <- sum((d/(h+ebx)-y)* x[,it]*ebx)
		  for(jt in 1:p)   sd[it,jt]=sum(  (d/(h+ebx)-d*ebx/(h+ebx)^2-y)*x[,it]*x[,jt]*ebx)
	  }  


	b <- b-solve(sd)%*%fd
	ebx <- exp(x%*%b)*exp(offset)
	l <- sum(d*log(h+ebx)-ds-y*ebx)
	bd <- abs(b-b0)
	if(max(bd)< tol) break()
  }
  
  
out <- list(b=b,sd=sd, nit=nit, loglik=c(l0,l))
out
}
