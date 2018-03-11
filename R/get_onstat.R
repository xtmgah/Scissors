#' Step 2 outlier detection: get_onstat
#' @export
get_onstat= function(resmat,rawmat,exonset,readconstr=10) {
  ## Step 2 : detect "on" outliers
  n2 = ncol(resmat)
  medvec = apply(rawmat,1,median)

  window_site = sort(c(exonset[,2],(exonset[,3]+1)))
  window_size = window_site[2:length(window_site)] - window_site[1:(length(window_site)-1)]
  window_site = window_site[1:(length(window_site)-1)]
  cbind(window_site,window_size)

  pdrate = matrix(0,nrow=length(window_site),ncol=n2)
  medianmat = matrix(0,nrow=length(window_site),ncol=n2)
  # read_max = rep(0,length(window_site))
  for (j in 1:length(window_site)) {
    carea = c(window_site[j]:(window_site[j]+window_size[j]-1))
    # read_max[j] = max(medvec[carea])
    medianmat[j,] = apply(rawmat[carea,],2,FUN=function(x){(median(x)>=readconstr)})
    pdrate[j,] = pd.rate.hy(apply(resmat[carea,],2,sum))
  }
  on_stat = apply(pdrate*medianmat,2,max)
  where_on = apply(pdrate*medianmat,2,which.max)
  start_loc = window_site[where_on]
  end_loc = start_loc + window_size[where_on] - 1

  return(list(on_stat=on_stat,on_site=cbind(start_loc,end_loc),
              window_site=window_site,window_size=window_size))
}
