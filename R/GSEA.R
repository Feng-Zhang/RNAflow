# ##' @title GSEA_go
# ##' @description This function would plot volcano based on DESeqDEGres output and save as pdf and png format
# ##' @param change_gene_List a numeric vector of the fold change of genes. It must have name with ENTREZID.
# ##' @param db the database, such as org.Hs.eg.db , org.Mm.eg.db and so on.
# ##' @param type  the type of ids which must be one of idType(OrgDb = db) shows
# ##' @param setReadable A logical value, mapping geneID to gene Symbol. If geneID is Symbol, the value should be FALSE.
# ##' @param ... Aavailable arguments to be passed to gseGO
# ##' @return NULL
# ##' @examples
# ##' \dontrun{
# ##' data(geneList, package="DOSE")
# ##' ego <- GSEAgo(geneList)
# ##' }
# ##' @importFrom clusterProfiler gseGO  setReadable
# ##' @importFrom utils modifyList
# ##' @export
# ##'
#
# GSEA_go <- function(change_gene_List,db=org.Hs.eg.db,type="ENTREZID",setReadable=FALSE,...){
#   #change_gene_List=geneList
#   ego <- NULL
#   change_gene_List <- sort(change_gene_List, decreasing = TRUE)
#   dotargs <- list(...)
#   defargs <- list(geneList=change_gene_List, ont="BP",OrgDb=db,eps=0,keyType = type)
#   ego <- do.call("gseGO",modifyList(defargs,dotargs))
#   if(setReadable) ego <- setReadable(ego, OrgDb = db) #names(x) <- value : 'names' attribute [2] must be the same length as the vector [1]
#   return(ego)
# }
#

##' @title GSEA_kegg
##' @description This function would plot volcano based on DESeqDEGres output and save as pdf and png format
##' @param change_gene_List a numeric vector of the fold change of genes. It must have name with ENTREZID.
##' @param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html', such as hsa for human, mmu for mouse.
##' @param db the annotation database of specific organism
##' @param type  the type of ids which must be one of idType(OrgDb = db) shows
##' @param ... Aavailable arguments to be passed to gseGO
##' @return NULL
##' @examples
##' \dontrun{
##' data(geneList, package="DOSE")
##' ekegg <- GSEAkegg(geneList)
##' }
##' @importFrom clusterProfiler gseKEGG setReadable bitr
##' @importFrom org.Hs.eg.db org.Hs.eg.db
##' @importFrom utils modifyList
##' @export
##'

GSEA_kegg <- function(change_gene_List,organism="hsa",db=org.Hs.eg.db,type="SYMBOL",...){
  ekegg <- NULL
  if(type!="ENTREZID") {
    mappedIds <- bitr(names(change_gene_List), fromType=type, toType="ENTREZID", OrgDb=db)
    change_gene_List <- change_gene_List[mappedIds[,type]]
    names(change_gene_List) <- mappedIds[,"ENTREZID"]
  }
  change_gene_List <- sort(change_gene_List, decreasing = TRUE)
  dotargs <- list(...)
  defargs <- list(change_gene_List,organism=organism,keyType="ncbi-geneid",pvalueCutoff = 0.05)
  ekegg <- do.call("gseKEGG",modifyList(defargs,dotargs))
  ekegg <- setReadable(ekegg, OrgDb = db,keyType = "ENTREZID")
  return(ekegg)
}

