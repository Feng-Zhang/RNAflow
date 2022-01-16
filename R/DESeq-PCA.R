##' @title DESeq_PCA
##' @description This function would plot PCA based on DESeqDataSet object and save as pdf and png format
##' @param dds a matrix with row name
##' @param group_name the number of colData to split case and control group
##' @param label A logical value to plot label of sample
##' @param output_dir the path of output
##' @param filename NULL or a character, the filename of output. When value is NULL, the figures are not saved.
##' @return NULL
##' @examples
##' \dontrun{
##' input = createCountPhe()
##' dds = DEseqObj(input[[1]],input[[2]],ref_level="untreated")
##' DESeqObjPCA(dds)
##' }
##' @importFrom ggplot2 geom_text
##' @importFrom grDevices pdf png dev.off
##' @importFrom DESeq2 plotPCA vst varianceStabilizingTransformation
##' @export
##'



DESeq_PCA <- function(dds,group_name="condition",label=FALSE,output_dir=".",filename=NULL){
  #groupNum=1;
  name <- NULL
  stopifnot("DESeqDataSet"%in% class(dds))
  if(nrow(dds)>1000) vsd <- vst(dds, blind=FALSE) else vsd=varianceStabilizingTransformation(dds)
  p <- plotPCA(vsd, intgroup=group_name)
  if(label){p <- p+geom_text(aes(label=name), vjust = 'inward', hjust = 'inward')}

  if(!is.null(filename)){
    ggsave(filename=filename,plot=p,path=output_dir)
  }
  return(p)
}
