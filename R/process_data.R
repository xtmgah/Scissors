#' Processing data
#' @export
process_data = function(rawdata,exon,input.type="whole_intron",
                        output.type="part_intron",intron.len=NULL,
                        logshift.val=NULL,param.grid=NULL,average="median",adjval=NULL,
                        center.type=1,center.adjval=NULL,
                        weight.method="mean",max.weight=2) {

  # Annotate pileup data
  data.annotation = annotate_pileup(pileup=rawdata,exon=exon,
                             input.type=input.type,output.type=output.type,intron.len=intron.len)
  exonset = data.annotation$exonset; # save exon positions
  dai = data.annotation$dai # data annotation information

  # Log transform data
  data.log = logtransform_data(data=data.annotation$out.pileup,logshift.val=logshift.val,
                               param.grid=param.grid,average=average,adjval=adjval)

  # Center and normalize data
  data.centering = center_data(data=data.log$outdata,type=center.type,average=average,adjval=center.adjval)

  # Weigh intronic region
  data.weighted = weigh_data(data=data.centering$outdata,dai=dai,method=weight.method,max.weight=max.weight)

  # Save output
  dataI = data.annotation$out.pileup
  datalog = data.log$outdata
  dataS = data.weighted$outdata

  return(list(dataI=dataI,datalog=datalog,dataS=dataS,dai=dai,exonset=exonset,
              logshift.val=data.log$logshift.val,
              msf=data.centering$msf,data.center=data.centering$data.center,
              w=data.weighted$w))
}
