#' Obtain median scale factors based on a adjusted center
#'
#' @export
scale.factor=function(X,average="mean",trim=0.1,adjval=NULL,robust.lm=TRUE) {

  mean.vec = apply(X,1,FUN=function(x){adj.center(x,average=average,trim=trim,adjval=adjval)})
  if (sum(mean.vec^2)>0) {
    if (robust.lm) {
      require(MASS)
      msf = apply(X,2,FUN=function(x){rlm(x ~ mean.vec - 1,psi=psi.bisquare)$coefficients})
    } else {
      msf = as.vector(t(mean.vec)%*%X/(sum(mean.vec^2)));
    }
  } else {
    msf = rep(0,ncol(X))
  }
  NewX = X - mean.vec%*%t(msf)
  return(list(msf=msf,mean.vec=mean.vec,NewX=NewX))
}
