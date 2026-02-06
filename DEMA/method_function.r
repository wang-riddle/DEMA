##function of different method
DEMA=function(Xt,Yt,Xs.list,Ys.list,M,alpha=.9,c=1,slide=3,type="LASSO"){
  ####target:Xt yt ;source:Xs ys;
  nt=length(Yt)
  ns=length(Ys.list[[1]])
  K=length(Xs.list)
  p=ncol(Xt)
  ###divide the target sample for aggregation
  I_abs=round(nt/slide)
  if(is.null(c)){
    threshold=floor((nt-I_abs)/log(nt-I_abs))
  }else{
    threshold=c*(floor((nt-I_abs)^alpha))
  }
  
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
    D_scan=Dscan(Xt.can,Yt.can)
    corr=order(D_scan,decreasing = T)
    shuffle_p=corr[1:threshold]
    pm=threshold/M
    
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
      penalty.res=penalty(Xt.sub,Yt.can,Xs.sub,Ys.list,type = type)
      beta_EST=penalty.res$beta.target
      beta_source[[m]]=penalty.res
      ##return to the p-dimension
      beta_EST.M[,m][idx.set[[m]]]=beta_EST
    }
    ##aggregation
    agg.re=Q.agg(beta_EST.M,Xt.agg,Yt.agg,n.sour=K*ns)
    weight.record[,s]=agg.re$lambda
    beta.record[,s]=agg.re$beta
    
  }
  beta.DEMA=apply(beta.record,1,mean)
  return(list(beta.DEMA=beta.DEMA,
              beta.record=beta.record,
              weight=weight.record,
              corr=corr,
              idx.set=idx.set,
              beta_source=beta_source))
}


##########################################################
Dscan=function(Xt,Yt){
  p=ncol(Xt)
  n=nrow(Xt)
  Dscan=rep(0,p)
  Xt=scale(Xt)
  Yt=Yt-mean(Yt)
  for (j in 1:p) {
    Xt_j=Xt[,j]
    Dscan[j]=crossprod(Xt_j,Yt)
  }
  return(Dscan)
}


##########################################################
penalty=function(Xt,Yt,Xs.list,Ys.list,type = "LASSO"){
  ##generate for algorithm.r；penalty for unified
  ##target:Xt,Yt; source:Xs.list,Ys.list
  #generate data
  K=length(Xs.list)
  p=ncol(Xt)
  Xapp=Xs.list
  Xapp[[K+1]]=Xt
  dimensions = lapply(Xapp, dim)
  zeros_list = lapply(dimensions, function(d) matrix(0, nrow = d[1], ncol = d[2]))
  X_final = list()
  for (i in 1:(K + 1)) {
    row_list <- list()
    for (j in 1:(K + 1)) {
      if (j == i) {
        row_list[[j]] <- Xapp[[i]]
      } else if(j == K + 1) {
        row_list[[j]] <- Xapp[[i]]
      } else {
        row_list[[j]] <- zeros_list[[i]]
      }
    }
    X_final[[i]] <- do.call(cbind, row_list)
  }
  X_final <- do.call(rbind, X_final)
  
  Yapp=Ys.list
  Yapp[[K+1]] = Yt
  Y_final <- do.call(c, Yapp)
  #theory value for penalty, low dimension
  N=sum(sapply(Xs.list, nrow)) + nrow(Xt)
  lambda_source=sapply(Xs.list, function(x) nrow(x)/N * sqrt(ncol(x) / nrow(x)))
  lambda=unlist(lapply(1:length(Xs.list), function(i) rep(lambda_source[i], ncol(Xs.list[[i]]))))
  lambda_target=sqrt( ncol(Xt) / N )
  lambda=c(lambda,rep(lambda_target, ncol(Xt)))
  #conduct lasso
  penalty.vec=lambda/max(lambda)
  if(type == "LASSO"){
    penalty_fit=cv.glmnet(X_final, Y_final, alpha=1, intercept=FALSE, standardize=FALSE,
                          penalty.factor=penalty.vec, nfolds = 10,
                          lambda.min.ratio=0.01)
  }else if(type == "SCAD"){
    penalty_fit=cv.ncvreg(X_final, Y_final, 
                          penalty = "SCAD",      # 对应 glmnet 的 alpha=1
                          intercept = FALSE, 
                          standardize = FALSE,
                          penalty.factor = penalty.vec, # 参数名完全一致
                          nfolds = 10,
                          eps = 0.01)
  }
  beta=coef(penalty_fit, s = "lambda.min")[-1]
  beta.mat=matrix(beta, nrow = p, ncol = K+1, byrow = FALSE)
  beta.final=cbind(beta.mat[,1:K] + beta.mat[,K+1],beta.mat[,K+1])
  beta.source=beta.final[,1:K]
  beta.target=beta.final[,K+1]
  
  beta.res.list=list()
  beta.res.list$beta.source=beta.source
  beta.res.list$beta.target=beta.target
  
  return(beta.res.list)
}

##########################################################
###########################################################
mse.fun<- function(beta,est, X.test=NULL){
  pred.err<-NA
  est.err<- sum((beta-est)^2)
  
  if(!is.null(X.test)){
    pred.err<-  mean((X.test%*%(beta-est))^2)
  }
  return(list(est.err=est.err, pred.err= pred.err))
}

###########################################################
Q.tar=function(Xt,Yt,M,cand="chunked"){
  nt=length(Yt)
  p=ncol(Xt)
  ###divide the target sample for aggregation
  I_abs=round(nt/3)
  pm=p/M
  idx.set=list()
  ##3 times record
  beta.record=matrix(0, nrow = p, ncol = 3)
  weight.record=matrix(0, nrow = M, ncol = 3)
  for (s in 1:3) {
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
    if(cand=="all in"){
      Xt.can=scale(Xt.can)
      corr=crossprod(Xt.can,Yt.can)/nrow(Xt.can)
      shuffle_p=order(corr,decreasing = T)
      if(s==1){M=2^M-1}
      weight.record=matrix(0, nrow = M, ncol = 3)
      for (m in 1:M) {
        unnested=1:p
        sub_vectors <- split(unnested, rep(1:5, each = 100))
        all_combinations <- unlist(
          lapply(1:5, function(k) combn(5, k, simplify = FALSE)),
          recursive = FALSE
        )
        idx.set <- lapply(
          all_combinations, 
          function(indices) unlist(sub_vectors[indices])
        )
        Xt.model[[m]]=Xt.can[,idx.set[[m]]]
      }
    }else if(cand=="chunked"){
      weight.record=matrix(0, nrow = M, ncol = 3)
      for (m in 1:M) {
        shuffle_p=1:p
        idx.set[[m]]=shuffle_p[(pm*(m-1)+1):(pm*m)]
        Xt.model[[m]]=Xt.can[,idx.set[[m]]]
      }
    }else if(cand=="nested"){
      weight.record=matrix(0, nrow = M, ncol = 3)
      MSE.chunk=rep(0,M)
      for (m in 1:M) {
        chunk=(1:p)[(pm*(m-1)+1):(pm*m)]
        fit.ols=lm(Yt.can~Xt.can[,chunk]-1)
        MSE.chunk[m]=mean((Yt.can-Xt.can[,chunk]%*%coef(fit.ols))^2)
        
      }
      order.chunk=order(MSE.chunk,decreasing = T)
      shuffle_p=rep(0,p)
      for (m in 1:M) {
        i=order.chunk[m]
        shuffle_p[(pm*(m-1)+1):(pm*m)]=(1:p)[(pm*(i-1)+1):(pm*i)]
      }
      for (m in 1:M) {
        idx.set[[m]]=shuffle_p[1:(pm*m)]
        Xt.model[[m]]=Xt.can[,idx.set[[m]]]
      }
    }
    ##calculate estimators
    beta_EST.M=matrix(0,p,M)
    for (m in 1:M) {
      Xt.sub=Xt.model[[m]]
      cv.las.Q=cv.glmnet(Xt.sub,Yt.can)
      lambda.min=cv.las.Q$lambda.min
      fit.las.Q=glmnet(Xt.sub,Yt.can,lambda = lambda.min)
      beta_EST=coef(fit.las.Q)[-1]
      ##return to the p-dimension
      beta_EST.M[,m][idx.set[[m]]]=beta_EST
    }
    ##aggregation
    agg.re=Q.agg(beta_EST.M,Xt.agg,Yt.agg)
    weight.record[,s]=agg.re$lambda
    beta.record[,s]=agg.re$beta
    
  }
  beta.Qma=apply(beta.record,1,mean)
  return(list(beta.Qma=beta.Qma,
              beta.record=beta.record,
              weight=weight.record))
  
}





