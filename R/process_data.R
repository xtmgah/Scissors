#' Processing data
#' @export
process_data = function(rawdata,exon,input.type="whole_intron",
                        output.type="part_intron",intron.len=NULL,
                        logshift.val=NULL,param.grid=NULL,average="mean",trim=0.1,adjval=NULL,
                        center.type=1,center.adjval=NULL,
                        loop=TRUE,smoothness=0.7,draw.plot=FALSE,plot.main="Gene",
                        weigh=FALSE,weight.method="mean",max.weight=2) {

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
                                   loop=loop,center.type=center.type,smoothness=smoothness,draw.plot=draw.plot,
                                   main=plot.main,average=average,trim=trim,adjval=center.adjval)

  # Save output
  dataraw = data.annotation$out.pileup
  datalog = data.log$outdata

  # Weigh intronic region
  if (weigh) {
    data.weighted = weigh_data(data=data.normalized$outdata,dai=dai,method=weight.method,max.weight=max.weight)
    dataS = data.weighted$outdata
  } else {
    data.weighted=as.list(NULL);
    dataS = data.normalized$outdata;
  }

  return(list(dataraw=dataraw,datalog=datalog,dataS=dataS,dai=dai,exonset=exonset,
              logshift.val=data.log$logshift.val,
              msf=data.normalized$msf,
              g1.offset=data.normalized$g1.offset,g2.offset=data.normalized$g2.offset,
              goodcase=data.normalized$goodcase,
              data.center=data.normalized$data.center,
              w=data.weighted$w))
}
