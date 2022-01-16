##' @title createCountPhe
##' @description This function would generate a list including a count_data and col_data as input data for testing
##' @param GSE The number of GSE
##' @return a list
##' @examples input = createCountPhe()
##' @import pasilla
##' @importFrom utils read.csv read.table
##' @importFrom stats median
##' @importFrom GEOmeta saveGSE
##' @export
createCountPhe <- function(GSE=NULL){
  UseMethod('createCountPhe')
}

#' @rdname createCountPhe
#' @export
createCountPhe.default <- function(GSE=NULL){
  #library("pasilla")

  # import expression data
  pasCts <- system.file("extdata",
                        "pasilla_gene_counts.tsv",
                        package="pasilla", mustWork=TRUE)
  count_data <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))

  # import phenotype data
  pasAnno <- system.file("extdata",
                         "pasilla_sample_annotation.csv",
                         package="pasilla", mustWork=TRUE)
  col_data <- read.csv(pasAnno, row.names=1)
  col_data <- col_data[,c("condition","type")]
  col_data$condition <- factor(col_data$condition)
  col_data$type <- factor(col_data$type)
  rownames(col_data) <- sub("fb", "", rownames(col_data))

  # check if the order of individuals between expression and phenotype data
  all(rownames(col_data) %in% colnames(count_data))
  all(rownames(col_data) == colnames(count_data))
  count_data <- count_data[, rownames(col_data)]
  all(rownames(col_data) == colnames(count_data))

  # clean expression data by keeping genes with mean count >0, median count >0
  count_data <- count_data[rowMeans(count_data) > 1 & apply(count_data,1,median)>0 ,]

  return(list(count_data,col_data))
}

#' @rdname createCountPhe
#' @export
createCountPhe.character <- function(GSE="GSE24132"){
  #library(GEOmeta);data(phe_test)
  saveGSE(GSE, destdir = "../tmp")
  expr_filename <- dir("../tmp",pattern = paste0("^",GSE,".*GPL.*-matrix.txt$"))
  count_data <- read.table(paste0("../tmp/",expr_filename),header = TRUE,row.names = 1)
  col_data <- GEOmeta::phe_test[GEOmeta::phe_test$GSE==GSE,2:3]
  row.names(col_data) <- col_data[,"GSM"]

  all(rownames(col_data) %in% colnames(count_data))
  all(rownames(col_data) == colnames(count_data))
  count_data <- count_data[, rownames(col_data)]
  all(rownames(col_data) == colnames(count_data))

  return(list(count_data,col_data))
}

