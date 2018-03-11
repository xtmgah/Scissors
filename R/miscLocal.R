#' Step 2 outlier detection
#' @export
miscLocal = function(miscGlobalobj,covermat,exonset,
                      winlength=100,readconstr=10,
                      siglev=1e-5,cutoff=NULL) {
  resmat = miscGlobalobj$resmat;
  rawmat = covermat;
  n = ncol(resmat); out1 = miscGlobalobj$outliers;
  if (length(out1)>=1) {
    resmat = resmat[,-out1]; rawmat = covermat[,-out1] ;
  }

  resmat_out = detect_onoffout(resmat=resmat,rawmat=rawmat,exonset=exonset,
                               winlength=winlength,readconstr=readconstr,
                               siglev=siglev,cutoff=cutoff)

  changeID = Relist.hy(1:n,out1);
  out2 = changeID$new[resmat_out$outliers];
  out2on = changeID$new[resmat_out$on_outliers]; out2off = changeID$new[resmat_out$off_outliers];
  ## Step 2 MOD
  MOD = matrix(0,nrow=nrow(resmat),ncol=n)
  MOD[,-out1] = resmat_out$MOD;

  cutoff = resmat_out$cutoff;
  out2stat = resmat_out$onoff_stat;
  out2site = resmat_out$onoff_site[resmat_out$outliers,];
  if (!is.null(out2site)) {
    if (is.null(dim(out2site))) {
      out2site = matrix(out2site,ncol=2)
    }
  }

  return(list(outliers=out2,outliers.on=out2on,outliers.off=out2off,
              out2stat=out2stat,cutoff=cutoff,
              MOD=MOD,out2site=out2site,
              out2res=resmat_out));
}