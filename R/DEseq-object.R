##' @title create_DEseq
##' @description This function would generate a object of DESeq2 class
##' @param count_data a matrix with row name
##' @param col_data a dataframe include group information
##' @param design_names designs with multiple variables, e.g., group + condition, and designs with interactions, e.g., genotype + treatment + genotype:treatment.
##' @param group_name the name of col_data to split case and control group
##' @param ref_level a character which is the reference group
##' @return a list
##' @examples
##' \dontrun{
##' input = createCountPhe()
##' dds = DEseqObj(input[[1]],input[[2]],ref_level="untreated")
##' }
##' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts
##' @importFrom stats relevel as.formula
##' @export
##'


create_DEseq <- function(count_data,col_data,design_names="condition+type",group_name="condition",ref_level="control"){
  #count_data=input[[1]];col_data=input[[2]];groupNum=1
  if(!is.data.frame(col_data)) stop("col_data is not a data frame")
  if(all(rownames(col_data) == colnames(count_data))){
    print("The id order between gene count file and phenotype file is identical without modification!")
  } else {
    ids <- intersect(colnames(count_data),rownames(col_data))
    count_data <- count_data[, ids]
    col_data <- col_data[ids,,drop=F]
    print(paste0("After modifying, the id order between gene count file and phenotype file is ",all(rownames(col_data) == colnames(count_data))))
    print(paste0("And ",length(ids)," individuals existed in both count_data and col_data are kept for analysis"))
  }


  # gene expression analysis
  testFormula <- paste0("~ ",design_names)
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = as.formula(testFormula) )
  try(`$`(dds,group_name) <- relevel(`$`(dds,group_name), ref = ref_level),silent=TRUE)
  dds <- dds[ rowSums(counts(dds)) > 1 & apply(counts(dds),1,median)>0, ]
  dds <- DESeq(dds)

  return(dds)
}
