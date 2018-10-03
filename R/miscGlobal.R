#' Detect global shape variants in RNA-seq data
#'
#' This function discovers outlying subjects whose RNA-seq have abnormal shapes
#' and provides the most outlying direction for each outlier.
#'
#' @param data precessed RNA-seq data from \link{process_data}
#' @param siglev the significance level of the Chi-squared distribution. Default is 1e-10.
#' @param subt.mean logical, whether to subtract mean before SVD. Default is FALSE.
#' @param PCnum the number of PCs to be used.
#' If NULL (the default) the number of PCs will be estimated.
#' @param maxPCnum the maximum number of PCs to be used. Default is 20.
#' @param eps tuning parameter for estimating the number of PCs.
#' If NULL (default), eps = (1/nrow(data))
#' @param rm.PCdir logical or numeric, the indices of PCs to be excluded from further study.
#' If TRUE (the default), this automatically excludes a set of PCs that deviate from normality.
#' If a sequence is specified, the corresponding PCs will be excluded.
#' @param ADcutoff a cutoff value for checking the normality based on Anderson-Darling test statistic
#' @param filter.dir logical, wether to filter out directions that deviate from normality.
#' Default is TRUE.
#' @param less.return logical, whether to show less results and not to return large matrices.
#' Default is TRUE.
#'
#' @export
miscGlobal = function(data,siglev=1e-10,subt.mean=FALSE,
                    PCnum=NULL,maxPCnum=20,eps=NULL,
                    rm.PCdir=TRUE,ADcutoff=10,
                    filter.dir=TRUE,
                    less.return=TRUE) {

  d = nrow(data); n = ncol(data);
  z.pca = pca.hy(data,subt.mean=subt.mean);
  projmat = z.pca$projmat; eigenval = z.pca$eigenval; dirmat = z.pca$dirmat;

  out = pca2misc(scoremat=projmat,eigenvalue=eigenval,dimension=d,
                  rm.PCdir=rm.PCdir,ADcutoff=ADcutoff,siglev=siglev,
                  PCnum=PCnum,maxPCnum=maxPCnum,eps=eps,
                  filter.dir=filter.dir,
                  less.return=FALSE);
  resmat = data - matrix(z.pca$dirmat[,1:out$K],ncol=out$K)%*%matrix(z.pca$projmat[1:out$K,],nrow=out$K)

  PCsubset = out$PCsubset; M = out$M; NPS = out$NPS;
  if (M==0) {
    MOD = NULL;
  } else if (M==1) {
    MOD = matrix(rep(dirmat[,PCsubset],n),ncol=n);
  } else {
    MOD = dirmat[,PCsubset]%*%out$directions;
  }

  if (less.return) {
    NPS = MOD = z.pca = resmat = NULL;
  }
  cutoff = sqrt(qchisq(p=(1-siglev),df=length(PCsubset)))

  return(list(OS=out$OS,OSpval=out$OSpval,NPS=NPS,MOD=MOD,pca=z.pca,resmat=resmat,
              outliers=out$outliers,outOS=out$outOS,outliers.sort=out$outliers.sort,outOS.sort=out$outOS.sort,
              M=out$M,PCsubset=out$PCsubset,rm.PCdir=out$rm.PCdir,cutoff=cutoff));
}
