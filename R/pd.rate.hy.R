pd.rate.hy = function(x) {
  # projection depth
  return((x-median(x))/mad(x))
}
