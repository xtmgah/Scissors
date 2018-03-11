#' Shift parameter
getShiftParam = function(X,param.grid=NULL,average="median",adjval=NULL,plot=FALSE) {

  if (is.null(param.grid)){
    param.grid = seq(1,10,by=1);
  }
  res = apply(matrix(param.grid,ncol=1),1,
              FUN=function(x){SlogADstat.hy(X=X,shift.val=x,average=average,adjval=adjval)$ADval})
  optim.idx = which.min(res);
  optim.param = param.grid[optim.idx];
  res.temp = SlogADstat.hy(X=X,shift.val=optim.param,average=average,adjval=adjval);
  optim.mvs = res.temp$medscale; # median vector scores
  optim.stat = res.temp$ADval;
  
  if (is.infinite(optim.stat)) {
    warning("Minimum = Inf")
  }
  
  if (plot){
    plot(param.grid,res,type="o",ylab="A-D statistic",xlab="Shift value");
    abline(v=optim.param,col="red");
    text(x=optim.param+(0.2*(max(param.grid)-min(param.grid))),
         y=(min(res)+0.8*(max(res)-min(res))),labels=round(optim.param,digits=3),col="red");
  }
  return(list(optim.param=optim.param,optim.stat=optim.stat,optim.mvs=optim.mvs,
              optim.idx=optim.idx,stat.evals=res));
}
