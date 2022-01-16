# ##' @title KEGG
# ##' @description This function would plot volcano based on DESeqDEGres output and save as pdf and png format
# ##' @param ids a character vector of gene id
# ##' @param type  the type of ids which must be one of idType(OrgDb = db) shows
# ##' @param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html', such as hsa for human, mmu for mouse.
# ##' @param db the annotation database of organism, such as org.Hs.eg.db for human
# ##' @param outputDir the path of output
# ##' @param base_name the prefix name of outputs would be saved
# ##' @param ... Aavailable arguments to be passed to enrichKEGG
# ##' @return NULL
# ##' @importFrom stringr str_split_fixed
# ##' @importFrom clusterProfiler bitr enrichKEGG
# ##' @importFrom grDevices pdf jpeg dev.off
# ##' @importFrom graphics barplot
# ##' @importFrom enrichplot dotplot
# ##' @importFrom org.Hs.eg.db org.Hs.eg.db
# ##' @importFrom utils modifyList
# ##' @export
# ##'
#
# KEGG <- function(ids,type="ENTREZID",organism="hsa",db=org.Hs.eg.db,outputDir=".",base_name=NULL,...){
#   dgeKegg <- NULL
#   if(type!="ENTREZID") ids <- bitr(ids, fromType=type, toType="ENTREZID", OrgDb=db)[,2]
#   dotargs <- list(...)
#   defargs <- list(gene = ids, organism=organism,keyType = "ncbi-geneid", pvalueCutoff = 0.05)
#   dgeKegg <- do.call("enrichKEGG",modifyList(defargs,dotargs))
#
#   if(!is.null(base_name)){
#     pdf(file = paste0(outputDir,"/",base_name,"-KEGG-dotplot.pdf"),height=4,width=8)
#     if(!is.null(dgeKegg)) print(dotplot(dgeKegg,showCategory=20))  else plot(0,main="Sorry,No gene can be enriched")
#     dev.off()
#
#     jpeg(filename = paste0(outputDir,"/",base_name,"-KEGG-dotplot.jpg"),height=2000,width=5000,res=300,type="cairo")
#     if(!is.null(dgeKegg)) print(dotplot(dgeKegg,showCategory=20))  else plot(0,main="Sorry,No gene can be enriched")
#     dev.off()
#
#     pdf(file = paste0(outputDir,"/",base_name,"-KEGG-barplot.pdf"),height=4,width=8)
#     if(!is.null(dgeKegg)) print(barplot(dgeKegg, showCategory=20, x = "GeneRatio")) else plot(0,main="Sorry,No gene can be enriched")
#     dev.off()
#
#     jpeg(filename = paste0(outputDir,"/",base_name,"-KEGG-barplot.jpg"),height=2000,width=5000,res=300,type="cairo")
#     if(!is.null(dgeKegg)) print(barplot(dgeKegg, showCategory=20, x = "GeneRatio")) else plot(0,main="Sorry,No gene can be enriched")
#     dev.off()
#   }
#
#   return(dgeKegg)
# }
