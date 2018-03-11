ADstat.hy=function(x){
  # Based on ADStatQF function in MATLAB by Qing Feng
  # This function is used to calculate the Anderson Darling Test Statistics of
  # Standard normal distribution
  
  n = length(x);
  if (n < 7){
    return(print("Sample size must be greater than 7"));
  } else {
    xs = sort(x);
    f = pnorm(xs,mean(xs),sd(xs));
    i = 1:n ;
    S = sum(((2*i-1)/n)*(log(f)+log(1-rev(f))));
    ADval = -n-S;
    return(ADval);
  }
}