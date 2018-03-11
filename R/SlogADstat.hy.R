SlogADstat.hy = function(X,shift.val=1,average="median",adjval=NULL){
  # log transform
  datalog = log10(X + shift.val) - log10(shift.val) ;
  # calculate the median scale factor
  msf = scale.factor(X=datalog,average=average,adjval=adjval)$msf

  ADval = ADstatWins.hy(msf);
  return(list(msf=msf,ADval=ADval));
}
