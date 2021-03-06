#' Annotating pileup
#'
#' This function is used to get pileup data based on new transcript
#' for the purpose of the downstream outlier detection analysis.
#' @param pileup raw coverage data
#' @param exon exon
#' @param input.type type of intronic region contained in pileup,
#' with choices "whole_intron", "part_intron", or "only_exon"; the first is the default.
#' @param output.type type of intronic region that will be included in output,
#' with choices "whole_intron", "part_intron", or "only_exon"; the second is the default.
#' @param intron.len a length of intronic region that will be included in output.
#' @keywords
#' @export
#' @examples
annotate_pileup = function(pileup,exon,input.type="whole_intron",
                           output.type="part_intron",intron.len=NULL) {
  # input.type  = "whole_intron", "part_intron", "only_exon"
  # output.ytpe = "whole_intron", "part_intron", "only_exon"

  if (input.type=="whole_intron") {

    if (output.type=="whole_intron") {
      ep.new = find.exon.hy(exon,is.intron=TRUE,num.intron=NULL) ;
      output = pileup;
    } else if (output.type=="part_intron") {
      if (is.null(intron.len)) {
        intron.len = ceiling(len.intron.hy(exon=exon)*0.5);
      }
      ep.new = find.exon.hy(exon,is.intron=TRUE,num.intron=intron.len) ;
      output = pileup[ep.new$coverage.col,] ;    #   Area to be included.
    } else if (output.type=="only_exon") {
      intron.len = 0;
      ep.new = find.exon.hy(exon,is.intron=TRUE,num.intron=intron.len) ;
      output = pileup[ep.new$coverage.col,] ;    #   Area to be included.
    } else {
      stop(output.type," is not an option for output.type.")
    }

  } else if (input.type=="part_intron") {

    if (output.type=="whole_intron") {
      stop(output.type," is not an option when input.type=part_intron.")
    } else if (output.type=="part_intron") {
      if (is.null(intron.len)) {
        intron.len = ceiling(len.intron.hy(exon=exon)*0.5);
      }
      ep.new = find.exon.hy(exon,is.intron=TRUE,num.intron=intron.len) ;
      output = pileup ;    #   Area to be included.
      ##### If the dimension of the given pileup data is not equal to the expected one
      ##### from the intron.len calculated, should display an error. Fix this!
    } else if (output.type=="only_exon") {
      intron.len.temp = ceiling(len.intron.hy(exon=exon)*0.5);
      ep.new.temp = find.exon.hy(exon,is.intron=TRUE,num.intron=intron.len.temp) ;

      exonic.region = c()
      for (i in 1:nrow(ep.new.temp$ep)) {
        exonic.region = c(exonic.region,(ep.new.temp$epl[i]:ep.new.temp$epr[i]));
      }
      output = pileup[exonic.region,];
      rm(intron.len.temp,ep.new.temp,exonic.region);

      intron.len = 0;
      ep.new = find.exon.hy(exon,is.intron=TRUE,num.intron=intron.len) ;
    } else {
      stop(output.type," is not an option for output.type.")
    }

  } else if (input.type=="only_exon") {

    if (output.type=="whole_intron") {
      stop(output.type," is not an option when input.type==only_exon.")
    } else if (output.type=="part_intron") {
      stop(output.type," is not an option when input.type==only_exon.")
    } else if (output.type=="only_exon") {
      intron.len = 0;
      ep.new = find.exon.hy(exon,is.intron=TRUE,num.intron=intron.len) ;
      output = pileup ;
    } else {
      stop(output.type," is not an option for output.type.")
    }

  } else {
    stop(input.type," is not an option for input.type.")
  }

  # Save exon position as a matrix for the new pileup
  epm0 = ep.new$ep; epm = ep.new$ep;

  # Change the directions
  seqdir = strsplit(exon,":")[[1]][3] ;
  if (seqdir=="-") {
    output = apply(output,2,rev); # log counts
    seqlen = epm0[nrow(epm0),4]
    epm[,1] = seqlen - rev(epm0[,4]) + 1;
    epm[,2] = seqlen - rev(epm0[,3]) + 1;
    epm[,3] = seqlen - rev(epm0[,2]) + 1;
    epm[,4] = seqlen - rev(epm0[,1]) + 1;
  }

  # pileup annotation info (dai)
  dai = as.list(NULL)
  dai$epm = epm
  dai$intron.len = intron.len
  return(list(out.pileup=output,dai=dai))
}
##########################################################
# annotate.pileup.v1 = function(pileup,exon,pileup.type="whole_intron",
#                            output.type="part_intron",intron.len=NULL) {
#   # pileup.type = "whole_intron", "part_intron", "only_exon"
#   # output.ytpe = "whole_intron", "part_intron", "only_exon"
#
#   if (pileup.type=="whole_intron") {
#
#     if (output.type=="whole_intron") {
#       ep.new = find.exon.hy(exon,is.intron=TRUE,num.intron=NULL) ;
#       output = pileup;
#     } else if (output.type=="part_intron") {
#       if (is.null(intron.len)) {
#         intron.len = ceiling(len.intron.hy(exon=exon)*0.5);
#       }
#       ep.new = find.exon.hy(exon,is.intron=TRUE,num.intron=intron.len) ;
#       output = pileup[ep.new$coverage.col,] ;    #   Area to be included.
#     } else if (output.type=="only_exon") {
#       intron.len = 0;
#       ep.new = find.exon.hy(exon,is.intron=TRUE,num.intron=intron.len) ;
#       output = pileup[ep.new$coverage.col,] ;    #   Area to be included.
#     } else {
#       stop(output.type,"is not an option for output.type.")
#     }
#
#   } else if (pileup.type=="part_intron") {
#
#     if (output.type=="whole_intron") {
#       stop(output.type,"is not an option when pileup.type=part_intron.")
#     } else if (output.type=="part_intron") {
#       if (is.null(intron.len)) {
#         intron.len = ceiling(len.intron.hy(exon=exon)*0.5);
#       }
#       ep.new = find.exon.hy(exon,is.intron=TRUE,num.intron=intron.len) ;
#       output = pileup ;    #   Area to be included.
#       ##### If the dimension of the given pileup data is not equal to the expected one
#       ##### from the intron.len calculated, should display an error. Fix this!
#     } else if (output.type=="only_exon") {
#       intron.len.temp = ceiling(len.intron.hy(exon=exon)*0.5);
#       ep.new.temp = find.exon.hy(exon,is.intron=TRUE,num.intron=intron.len.temp) ;
#
#       exonic.region = c()
#       for (i in 1:nrow(ep.new.temp$ep)) {
#         exonic.region = c(exonic.region,(ep.new.temp$epl[i]:ep.new.temp$epr[i]));
#       }
#       output = pileup[exonic.region,];
#       rm(intron.len.temp,ep.new.temp,exonic.region);
#
#       len.intron = 0;
#       ep.new = find.exon.hy(exon,is.intron=TRUE,num.intron=len.intron) ;
#     } else {
#       stop(output.type,"is not an option for output.type.")
#     }
#
#   } else if (pileup.type=="only_exon") {
#
#     if (output.type=="whole_intron") {
#       stop(output.type,"is not an option when pileup.type==only_exon.")
#     } else if (output.type=="part_intron") {
#       stop(output.type,"is not an option when pileup.type==only_exon.")
#     } else if (output.type=="only_exon") {
#       intron.len = 0;
#       ep.new = find.exon.hy(exon,is.intron=TRUE,num.intron=intron.len) ;
#       output = pileup ;
#     } else {
#       stop(output.type,"is not an option for output.type.")
#     }
#
#   } else {
#     stop(pileup.type,"is not an option for output.type.")
#   }
#
#   # Save exon position as a matrix for the new pileup
#   epm0 = ep.new$ep; epm = ep.new$ep;
#
#   # Change the directions
#   seqdir = strsplit(exon,":")[[1]][3] ;
#   if (seqdir=="-") {
#     output = apply(output,2,rev); # log counts
#     seqlen = epm0[nrow(epm0),4]
#     epm[,1] = seqlen - rev(epm0[,4]) + 1;
#     epm[,2] = seqlen - rev(epm0[,3]) + 1;
#     epm[,3] = seqlen - rev(epm0[,2]) + 1;
#     epm[,4] = seqlen - rev(epm0[,1]) + 1;
#   }
#
#   # pileup annotation info (dai)
#   dai = as.list(NULL)
#   dai$epm = epm
#   dai$intron.len = intron.len
#   return(list(out.pileup=output,dai=dai))
# }
