# ##' @title GO
# ##' @description This function would plot volcano based on DESeqDEGres output and save as pdf and png format
# ##' @param ids A character vector of gene id
# ##' @param type The type of ids which must be one of idType(OrgDb = db) shows
# ##' @param db The database, such as org.Hs.eg.db , org.Mm.eg.db and so on.
# ##' @param outputDir The path of output
# ##' @param base_name The prefix name of outputs would be saved
# ##' @param ... Aavailable arguments to be passed to enrichGO.
# ##' @return NULL
# ##' @importFrom stringr str_split_fixed
# ##' @importFrom clusterProfiler bitr enrichGO
# ##' @importFrom grDevices jpeg png dev.off
# ##' @importFrom graphics barplot
# ##' @importFrom utils modifyList
# ##' @export
# ##'
#
# GO <- function(ids,type="ENTREZID",db=org.Hs.eg.db,outputDir=".",base_name=NULL,...){
#   if(type!="ENTREZID") ids <- bitr(ids, fromType=type, toType="ENTREZID", OrgDb=db)[,2]
#   dotargs <- list(...)
#   defargs <- list(gene = ids,#括号内必须加逗号，否则是矩阵
#                OrgDb = db,
#                ont = "BP",
#                pAdjustMethod = "BH",
#                pvalueCutoff  = 0.05,  readable = TRUE)
#   ego <- do.call("enrichGO",modifyList(defargs,dotargs)) #enrichGO不识别ENSEMBL
#
#   if(!is.null(base_name)){
#     pdf(file = paste0(outputDir,"/",base_name,"-GO-barplot.pdf"),height=4,width=8)
#     if(nrow(ego)>0) print(barplot(ego, showCategory=20, x = "GeneRatio")) else plot(0,main="Sorry,No gene can be enriched")
#     dev.off()
#     jpeg(filename = paste0(outputDir,"/",base_name,"-GO-barplot.jpg"),height=2000,width=4000,res=300,type="cairo")
#     if(nrow(ego)>0) print(barplot(ego, showCategory=20, x = "GeneRatio")) else plot(0,main="Sorry,No gene can be enriched")
#     dev.off()
#
#
#     pdf(file = paste0(outputDir,"/",base_name,"-GO-dotplot.pdf"),height=4,width=8)
#     if(nrow(ego)>0) print(dotplot(ego,showCategory=20)) else plot(0,main="Sorry,No gene can be enriched")
#     dev.off()
#
#     jpeg(filename = paste0(outputDir,"/",base_name,"-GO-dotplot.jpg"),height=2000,width=4000,res=300,type="cairo")
#     if(nrow(ego)>0) print(dotplot(ego,showCategory=20)) else plot(0,main="Sorry,No gene can be enriched")
#     dev.off()
#   }
#
#   return(ego)
# }




