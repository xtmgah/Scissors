#' Step 2 outlier detection: get_offstat
#' @export
get_offstat= function(resmat,rawmat,exonset,
                      winlength=100,readconstr=10) {
  ## Step 2 : detect "off" outliers
  require(zoo)
  # Find eligible regions
  n2 = ncol(resmat)
  medvec = apply(rawmat,1,median)

  window_min = rollapply(medvec,width=winlength,FUN=min)
  window_idx = which(window_min>=readconstr)

  # Find exons whose lengths are less than the pre-determined size.
  smallexon_0 = which((exonset[,3]-exonset[,2]+1)<winlength)
  smallexon_min = c()
  for (j in smallexon_0) {
    smallexon_min = min(medvec[c(exonset[j,2]:exonset[j,3])])
  }
  smallexon = smallexon_0[which(smallexon_min>=readconstr)]

  # Determine window size
  winsize_vec = rep(0,length(window_min))
  winsize_vec[window_idx] = winlength;
  winsize_vec[exonset[smallexon,2]] = exonset[smallexon,3]-exonset[smallexon,2]+1

  # Final start_positions for calculating statistics
  window_site = which(winsize_vec>0)
  window_size = winsize_vec[window_site]
  pdrate = matrix(0,nrow=length(window_site),ncol=ncol(resmat))
  if (length(window_site)>0) {
    for (i in 1:length(window_site)) {
      carea = c((window_site[i]):(window_site[i]-1+window_size[i]));
      pdrate[i,] = -pd.rate.hy(apply(resmat[carea,],2,sum));
    }
    off_stat = apply(pdrate,2,max)
    where_on = apply(pdrate,2,which.max)
    start_loc = window_site[where_on]
    end_loc = start_loc + window_size[where_on] - 1
  } else {
    off_stat = rep(0,ncol(resmat))
    start_loc = rep(0,ncol(resmat))
    end_loc = rep(0,ncol(resmat))
  }

  return(list(off_stat=off_stat,off_site=cbind(start_loc,end_loc),
              window_site=window_site,window_size=window_size))
}
