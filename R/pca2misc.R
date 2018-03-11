pca2misc = function(scoremat,eigenvalue,dimension,
                   PCnum=NULL,maxPCnum=20,eps=NULL,rm.PCdir=TRUE,ADcutoff=10,
                   filter.dir=TRUE,
                   siglev=1e-10,
                   less.return=FALSE) {
  # Version 3: Add choosing PC directions
  # Computes the Stahel-Donoho outlyingness of every element in x
  # x can be univariate or multivariate
  # Assumes that dim(x) = n x d if multivariate
  #
  n = ncol(scoremat); d = dimension;
  eigenval = eigenvalue;
  projmat = scoremat;

  # Choose the number of spikes
  if (is.null(PCnum)){
    if (is.null(eps)){
      # eps = (1/nrow(Data));
      eps = (1/d);
    }
    result=PCnum.gamma.hy(eigenval=eigenval,d=d,n=n,eps=eps);
    K = min(result$K,maxPCnum);
  } else {
    K = PCnum;
  }
  ## Choose normal directions
  if (is.logical(rm.PCdir)){
    if (rm.PCdir) {
      badPCdir = which(apply(projmat[1:K,],1,ADstatWins.hy)>ADcutoff);
      if (length(badPCdir)>0) {
        PCsubset = c(1:K)[-badPCdir];
      } else {
        PCsubset = 1:K;
      }
    } else {
      PCsubset = 1:K;
      badPCdir = c();
    }
    rm.PCdir=badPCdir;
  } else {
    if (is.null(rm.PCdir) | (length(rm.PCdir)==0)) {
      PCsubset = 1:K;
    } else {
      PCsubset = c(1:K)[-rm.PCdir];
    }
  }
  M = length(PCsubset);

  # Get OS statistic
  if (M==0) {
    # cat(paste("!!!!!!!     ","No directions included","    !!!!!!!"),"\n");
    OS = rep(0,n);
    OSpval = rep(1,n);
    NPS = NULL;  # Projection Score
    if (less.return) {
      NPS = NULL;
    }
    return(list(OS=OS,OSpval=OSpval,NPS=NPS,
                outliers=NULL,outOS=NULL,outliers.sort=NULL,outOS.sort=NULL,
                M=M,PCsubset=PCsubset,rm.PCdir=rm.PCdir,K=K));
  } else {
    x = projmat[PCsubset,];   # projmat with the chosen directions

    if (is.null(dim(x))){
      if (filter.dir) {
        ADstat = ADstatWins.hy(x);
        if (ADstat > ADcutoff) {
          OS = rep(0,n);
          OSpval = rep(1,n);
          NPS = NULL;  # Projection Score
          directions = NULL;
          return(list(OS=OS,OSpval=OSpval,NPS=NPS,directions=directions,
                      outliers=NULL,outOS=NULL,outliers.sort=NULL,outOS.sort=NULL,
                      M=M,PCsubset=PCsubset,rm.PCdir=rm.PCdir,K=K));
        } else {
          OS = abs(x-median(x))/mad(x);
          OSpval = 1-pchisq(q=OS^2,df=M);
          NPS = matrix(rep(OS,n),ncol=n);
          directions = matrix(rep(x,n),ncol=n);
          if (less.return) {
            NPS = directions = NULL;
          }
          out = which(OSpval<siglev); outOS = OS[out];
          out.order = order(outOS,decreasing=TRUE);
          out.sort = out[out.order]; outOS.sort = outOS[out.order];
          return(list(OS=OS,OSpval=OSpval,NPS=NPS,directions=directions,
                      outliers=out,outOS=outOS,outliers.sort=out.sort,outOS.sort=outOS.sort,
                      M=M,PCsubset=PCsubset,rm.PCdir=rm.PCdir,K=K));
        }# if (ADstat > ADcutoff)
      } else {
        OS = abs(x-median(x))/mad(x);
        OSpval = 1-pchisq(q=OS^2,df=M);
        NPS = matrix(rep(OS,n),ncol=n);
        directions = matrix(rep(x,n),ncol=n);
        if (less.return) {
          NPS = directions = NULL;
        }
        out = which(OSpval<siglev); outOS = OS[out];
        out.order = order(outOS,decreasing=TRUE);
        out.sort = out[out.order]; outOS.sort = outOS[out.order];
        return(list(OS=OS,OSpval=OSpval,NPS=NPS,directions=directions,
                    outliers=out,outOS=outOS,outliers.sort=out.sort,outOS.sort=outOS.sort,
                    M=M,PCsubset=PCsubset,rm.PCdir=rm.PCdir,K=K));
      }
    } else {
      x = t(x);
      ndir = 300*M;
      A = generdir(x,ndir=ndir) # generates `ndir' directions (ndir by M)
      # Add more potential directions
      subprojmat = t(x);
      Scov = (subprojmat%*%t(subprojmat))/n;
      B0 = solve(Scov)%*%subprojmat;   # M by n
      B0 = B0[,-which(apply(B0,2,FUN=function(x){sqrt(sum(x^2))})<1e-10)];
      B1 = sweep(B0,2,apply(B0,2,FUN=function(x){sqrt(sum(x^2))}),"/")  # normalized M by n
      A = rbind(A,t(B1))  # ndir by M

      Y = x %*% t(A) # project x onto A (n by ndir)
      if (filter.dir) {
        ADstat = apply(Y,2,ADstatWins.hy);
        indir = which(ADstat<ADcutoff);
        A = A[indir,];
        Y = x %*% t(A) # project x onto A (n by length(indir))
        out_temp = apply(X=Y,MARGIN=2,
                         FUN= function(t) (t-median(t))/mad(t)); # n by length(indir)
        OS = apply(abs(out_temp),1,max);
        OSpval = 1-pchisq(q=OS^2,df=M);
        indexmax=unlist(apply(abs(out_temp),1,
                              function(x) min(which(x==max(x,na.rm=TRUE)))));
        NPS = out_temp[,indexmax];
        directions=t(A[indexmax,]);

        neg.index = which(diag(NPS)<0);
        if (length(neg.index)>0) {
          NPS[,neg.index] = -NPS[,neg.index];
          directions[,neg.index] = -directions[,neg.index];
        }
        if (less.return) {
          NPS = directions = NULL;
        }
        out = which(OSpval<siglev); outOS = OS[out];
        out.order = order(outOS,decreasing=TRUE);
        out.sort = out[out.order]; outOS.sort = outOS[out.order];
        return(list(OS=OS,OSpval=OSpval,NPS=NPS,directions=directions,
                    outliers=out,outOS=outOS,outliers.sort=out.sort,outOS.sort=outOS.sort,
                    M=M,PCsubset=PCsubset,rm.PCdir=rm.PCdir,K=K));
      } else {
        out_temp = apply(X=Y,MARGIN=2,
                         FUN= function(t) (t-median(t))/mad(t)); # n by length(indir)
        OS = apply(abs(out_temp),1,max);
        OSpval = 1-pchisq(q=OS^2,df=M);
        indexmax=unlist(apply(abs(out_temp),1,
                              function(x) min(which(x==max(x,na.rm=TRUE)))))
        NPS = out_temp[,indexmax];
        directions=t(A[indexmax,]);
        neg.index = which(diag(NPS)<0);
        if (length(neg.index)>0) {
          NPS[,neg.index] = -NPS[,neg.index];
          directions[,neg.index] = -directions[,neg.index];
        }
        if (less.return) {
          NPS = directions = NULL;
        }
        out = which(OSpval<siglev); outOS = OS[out];
        out.order = order(outOS,decreasing=TRUE);
        out.sort = out[out.order]; outOS.sort = outOS[out.order];
        return(list(OS=OS,OSpval=OSpval,NPS=NPS,directions=directions,
                    outliers=out,outOS=outOS,outliers.sort=out.sort,outOS.sort=outOS.sort,
                    M=M,PCsubset=PCsubset,rm.PCdir=rm.PCdir,K=K));
      }# if (filter.dir)
    }
  }
}
