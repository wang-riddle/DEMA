##Q-aggregation
Q.agg=function(beta.com,X,Y,total.step=10,n.sour)
{
  ##beta.com is a matrix which combines all beta estimator
  ##X,Y is test data(target data)
  L=ncol(beta.com)
  p=nrow(beta.com)#dimension
  lambda=rep(0,L)#weight
  alpha=2/((1:total.step)+1)
  n=length(Y)+n.sour
  for (k in 1:total.step) {
    alpha_k=alpha[k]
    ###Q
    Q.value=rep(0,L)
    for (j in 1:L) {
      e_j=rep(0,L)
      e_j[j]=1
      lambda_j=lambda+alpha_k*(e_j-lambda)
      beta.para=2*sqrt(n)/log(8*sqrt(n))
      Q.value[j]=Q(X,Y,beta.com,lambda_j,v=1/2,beta.para,pi=rep(1/L,L))
      
    }
    idx=which.min(Q.value)
    e_idx=rep(0,L)
    e_idx[idx]=1
    lambda=lambda+alpha_k*(e_idx-lambda)
  }
  beta=beta.com%*%lambda
  return(list(lambda=lambda,beta=beta))
}
####
Q=function(X,Y,beta.com,lambda,v=1/2,beta.para,pi)
{
  MSE1=sum((Y-X%*%(as.numeric(beta.com%*%lambda)))^2)
  MSE2=t(lambda)%*%colSums((Y-X%*%beta.com)^2)
  K=t(lambda)%*%log(1/pi)
  n=length(Y)
  Q.result=(1-v)*MSE1+v*MSE2+(beta.para/n)*K
  return(Q.result)
}


