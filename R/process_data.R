#' Processing data
#' @export
process_data = function(rawdata,exon,input.type="whole_intron",
                        output.type="part_intron",intron.len=NULL,
                        logshift.val=NULL,param.grid=NULL,average="median",adjval=NULL,
                        center.type=1,center.adjval=NULL,
                        smoothness=0.7,draw.plot=FALSE,plot.main="Gene",
                        weight.method="mean",max.weight=2) {

  # Annotate pileup data
  data.annotation = annotate_pileup(pileup=rawdata,exon=exon,
                             input.type=input.type,output.type=output.type,intron.len=intron.len)
  dai = data.annotation$dai # data annotation information
  exonset = dai$epm
  intron.len = dai$intron.len

  # Log transform data
  data.log = logtransform_data(data=data.annotation$out.pileup,logshift.val=logshift.val,
                               param.grid=param.grid,average=average,adjval=adjval)

  # Center and normalize data
  data.normalized = normalize_data(data=data.log$outdata,rawmat=data.annotation$out.pileup,exonset=exonset,
                                   center.type=center.type,smoothness=smoothness,draw.plot=draw.plot,
                                   main=plot.main,average=average,adjval=center.adjval)

  # Weigh intronic region
  data.weighted = weigh_data(data=data.normalized$outdata,dai=dai,method=weight.method,max.weight=max.weight)

  # Save output
  dataI = data.annotation$out.pileup
  datalog = data.log$outdata
  dataS = data.weighted$outdata

  return(list(dataI=dataI,datalog=datalog,dataS=dataS,dai=dai,exonset=exonset,
              logshift.val=data.log$logshift.val,
              msf=data.normalized$msf,g=data.normalized$g,data.center=data.normalized$data.center,
              w=data.weighted$w))
}
