#######function to combine two rna data and run SCTv2
#####use output from rna_fq_preprocess pipeline
Run_SCTV2 <- function(RNA_path1,,RNA_path2) {

  rna_path <- RNA_path1 
  seurat_object <- Read10X(data.dir = rna_path)
  seurat_object <- CreateSeuratObject(counts = seurat_object)
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & percent.mt <= 30)
  seurat_object <- seurat_object[rowSums(seurat_object@assays$RNA@counts)>0,]
  
  Chip1 <- seurat_object
  Chip1$chip <- "Chip1"
  
  rna_path <- RNA_path2 
  seurat_object <- Read10X(data.dir = rna_path)
  seurat_object <- CreateSeuratObject(counts = seurat_object)
  
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & percent.mt <= 30)
  seurat_object <- seurat_object[rowSums(seurat_object@assays$RNA@counts)>0,]
  
  Chip2 <- seurat_object
  Chip2$chip <- "Chip2"
  
  object_combined <- merge(Chip1,Chip2,add.cell.ids=c("Chip1","Chip2"))
  object_combined <- DietSeurat(object_combined)
  object_combined <- SCTransform(object_combined, vst.flavor = "v2", verbose = FALSE,vars.to.regress = c("chip")) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
    FindClusters(resolution = 0.6, verbose = FALSE)
  
  return(object_combined)
}
###Run the above function will generate a integrated seurat object with SCTV2 assay. The resulted clusters are used for downstream analysis