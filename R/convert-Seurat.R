#' @title Convert Seurat object to SpatialRNA object
#' @description Convert Seurat class to SpatialRNA class as the input data for RCTD analysis using spacexr package.
#' @param seurat A Seurat class object.
#' @param slot Specific assay data to get or set. see [SeuratObject::GetAssayData()] for detail information.
#' @param assay Specific assay to get data from or set data for; defaults to the default assay. see [SeuratObject::GetTissueCoordinates()] for detail information.
#' @param image Name of SpatialImage object to get coordinates for; if NULL, will attempt to select an image automatically.  see [SeuratObject::GetTissueCoordinates()] for detail information.
#' @param nUMI_name Optional, a colname (by pixel barcode) in Seurat object, which refer to total counts or UMI's appearing at each pixel. If not provided, nUMI will be assumed to be the total counts appearing on each pixel.
#' @importFrom SeuratObject  GetAssayData GetTissueCoordinates
#' @importFrom spacexr SpatialRNA
#' @export
Seurat_SpatialRNA <- function(seurat,slot = "counts",assay = "Spatial",image=NULL,nUMI_name="nCount_Spatial"){
  stopifnot("Seurat"==class(seurat))
  counts <- GetAssayData(seurat,slot=slot,assay=assay)
  coords <- GetTissueCoordinates(seurat,image=image)
  coords <- data.frame(x=coords[,"imagecol"],y=coords[,"imagerow"]*(-1),row.names = row.names(coords))
  nUMI <- seurat@meta.data[,nUMI_name]
  names(nUMI) <- colnames(seurat)
  sRNA <- SpatialRNA(coords, counts, nUMI)
  return(sRNA)
}

#' @title Convert Seurat object to Reference object
#' @description convert Seurat class to Reference class as the reference data for RCTD analysis using spacexr package
#' @inheritParams Seurat_SpatialRNA
#' @param cell_type_name A colname in `Seurat@meta.data`, which refer to the group or cell type.
#' @importFrom Seurat GetAssayData
#' @importFrom spacexr Reference
#' @export
Seurat_Reference <- function(seurat,slot = "counts",assay = "RNA",cell_type_name="cell_type",nUMI_name=NULL){
  stopifnot("Seurat"==class(seurat))
  stopifnot(cell_type_name %in% colnames(seurat@meta.data))

  counts <- GetAssayData(seurat,slot=slot,assay=assay)

  cell_types <- seurat@meta.data[,cell_type_name]
  names(cell_types) <- colnames(seurat)
  cell_types <- as.factor(cell_types)

  if(is.null(nUMI_name)){
    nUMI <- NULL
  } else {
    stopifnot(nUMI_name %in% colnames(seurat@meta.data))
    nUMI <- seurat@meta.data[,nUMI_name]
    names(nUMI) <- colnames(seurat)
  }

  reference <- Reference(counts, cell_types=cell_types,nUMI=nUMI)
  return(reference)
}
