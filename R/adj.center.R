adj.center = function(x,average="median",adjval=NULL) {
  ##  % Obtain adjusted center (mean or median)
  if (is.null(adjval)) {
    if (average=="median") {
      center.val = median(x);
    } else {
      center.val = mean(x);
    }
  } else {
    nonzero.pos = which(x>adjval)
    if (length(nonzero.pos)>0) {
      if (average=="median") {
        center.val = median(x[nonzero.pos])
      } else {
        center.val = mean(x[nonzero.pos])
      }
    } else {
      center.val = 0
    }
  }
  return(center.val);
}
