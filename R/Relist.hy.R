Relist.hy = function(a,b){
  # Input  : a = original list; b = number will be removed from a
  # Output : c = new list; d => d[k] = case number in c of a[k]
  if (is.null(b)) {
    c = a; d = a;
  } else {
    c = a[-b]; d = rep(0,length(a));
    for (i in 1:length(a)){
      if (length(which(c==a[i]))>0){
        d[i] = which(c==a[i]);
      }
    }
  }
  return(list(new=c,ord=d));
}
