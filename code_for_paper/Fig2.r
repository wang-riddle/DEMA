DEMA.sp=function(Xt,Yt,Xs.list,Ys.list,M,alpha=.9,c=1){
  ####target:Xt yt ;source:Xs ys;
  nt=length(Yt)
  ns=length(Ys.list[[1]])
  K=length(Xs.list)
  p=ncol(Xt)
  
  ###divide the target sample for aggregation
  slide=3
  I_abs=round(nt/slide)
  idx.set=list()
  ##3 times record
  beta.record=matrix(0, nrow = p, ncol = slide)
  weight.record=matrix(0, nrow = M, ncol = slide)
  for (s in 1:slide) {
    Itil=(I_abs*s-(I_abs-1)):(I_abs*s)
    #target sample for candidate model
    Xt.can=Xt[-Itil,]
    Yt.can=Yt[-Itil]
    #target sample for aggregation
    Xt.agg=Xt[Itil,]
    Yt.agg=Yt[Itil]
    ##randomly divide the variables for candidate model
    #shuffle_p=sample(1:p)
    Xt.model=list()
    Xs.model=list()
    #DC_scan=DC(Xt.can,Yt.can)
    #corr=order(DC_scan,decreasing = T)
    shuffle_p=1:32
    pm=4
    
    for (m in 1:M) {
      idx.set[[m]]=shuffle_p[1:(pm*m)]
      Xt.model[[m]]=Xt.can[,idx.set[[m]]]
    }
    for (k in 1:K) {
      Xs.model[[k]]=list()
      Xs.model.tmp = Xs.list[[k]]
      for (i in 1:M) {
        Xs.model[[k]][[i]]=Xs.model.tmp[,idx.set[[i]]]
      }
    }
    ##calculate estimators
    beta_EST.M=matrix(0,p,M)
    beta_source=list()
    for (m in 1:M) {
      Xt.sub=Xt.model[[m]]
      Xs.sub=list()
      for (k in 1:K) {
        Xs.sub[[k]]=Xs.model[[k]][[m]]
      }
      penalty.res=penalty(Xt.sub,Yt.can,Xs.sub,Ys.list)
      beta_EST=penalty.res$beta.target
      beta_source[[m]]=penalty.res
      ##return to the p-dimension
      beta_EST.M[,m][idx.set[[m]]]=beta_EST
    }
    ##aggregation
    agg.re=Q.agg(beta_EST.M,Xt.agg,Yt.agg)
    weight.record[,s]=agg.re$lambda
    beta.record[,s]=agg.re$beta
    
  }
  beta.DEMA=apply(beta.record,1,mean)
  return(list(beta.DEMA=beta.DEMA,
              beta.record=beta.record,
              weight=weight.record,
              #corr=corr,
              idx.set=idx.set,
              beta_source=beta_source))
}
