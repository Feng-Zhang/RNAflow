##' @title DESeqRes
##' @description This function would get the results based on DESeqDataSet object
##' @param dds a matrix with row name
##' @param fold_change_line a numeric, the threshold to define significant fold change
##' @param adj_pvalue_line a numeric, the threshold of significant log2(fold change)
##' @param ... Aavailable arguments to be passed to results
##' @return a dataframe
##' @examples
##' \dontrun{
##' input = createCountPhe()
##' dds = DEseqObj(input[[1]],input[[2]],ref_level="untreated")
##' res = DESeqRes(dds)
##' }
##' @importFrom  DESeq2 results
##' @importFrom utils modifyList
##' @export
##'


DESeq_res <- function(dds,fold_change_line=1,adj_pvalue_line=0.05,...){
  dotargs <- list(...)
  defargs <- list(object=dds)
  res <- do.call("results",modifyList(defargs,dotargs))
  res <- as.data.frame(res[order(res$padj), ])
  res <- define_de(res,fold_change_line=fold_change_line,adj_pvalue_line=adj_pvalue_line)
  return(res)
}

##' @title define_de
##' @description This function would create a new column to define the differential expression
##' @param DEres a data.frame contains fold_change, adj_p columns
##' @param fold_change a character, the name of column represented log2(fold change)
##' @param adj_pvalue a character, the name of column represented p-value
##' @param de_class a character, the name of created column represented differential expression
##' @param fold_change_line a numeric, the threshold to define significant fold change
##' @param adj_pvalue_line a numeric, the threshold to define significant p-value
##' @return a dataframe
##' @export
define_de <- function(DEres,fold_change="log2FoldChange",adj_pvalue="padj",de_class="regulate",fold_change_line=1,adj_pvalue_line=0.05){
  stopifnot(class(DEres)[1]=="data.frame" & all(c(fold_change,adj_pvalue)%in%colnames(DEres)))
  DEres[,de_class] <- "Normal"
  DEres[DEres[,fold_change] < -fold_change_line & DEres[,adj_pvalue] < adj_pvalue_line & !is.na(DEres[,adj_pvalue]), de_class] <- "Down"
  DEres[DEres[,fold_change] > fold_change_line & DEres[,adj_pvalue] < adj_pvalue_line & !is.na(DEres[,adj_pvalue]), de_class] <- "Up"
  return(DEres)
}

