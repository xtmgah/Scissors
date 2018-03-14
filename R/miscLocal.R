#' Detect local shape abnormality in RNA-seq data
#'
#' This function discovers outlying subjects whose RNA-seq have "local" abnormal shapes
#' and provides the most outlying window-direction for each outlier.
#'
#' @param miscGlobalobj object from \link{miscGlobal}
#' @param overmat raw coverage matrix, dataI from \link{process_data}
#' @param exonset data annotation matrix, exonset from \link{process_data}
#' @param winlength the window length. Default is 100.
#' @param readconstr the minimum mean read-counts. Default is 10.
#' @param siglev the significance level for a robust outlier detection. Default is 1e-5.
#' IF cutoff is specified, siglev is not used.
#' @param cutoff the cutoff value for outlying statistics.
#' If NULL, the cutoff value is computed based on the specified siglev.
#'
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
