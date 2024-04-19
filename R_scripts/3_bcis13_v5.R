#---preload----#####
library(uwot, lib.loc = "/opt/R/4.1.2/lib/R/library")
packageVersion("uwot")
#--uwot version:  0.1.14
library("Seurat")
packageVersion("Seurat")
library(SeuratWrappers)
#---4.1.0
library(dplyr)
library(plyr)
library(tibble)
library(insect)
library(ggplot2)
library(cowplot)
library(tidyr)
library(ComplexHeatmap)
library(scales)
library(copykit)
library(yarrr)
library(future)
library(scquantum)
library(readr)
library(png)
library(uwot)
packageVersion("uwot")
options(max.print = 200)

wd <- c("/volumes/USR2/wangkl/wscDR/wdr_analysis/")
setwd(wd)
source("./scripts/0_wdr_functions.R")
load("./pre_load_data/pre_load_data.rda")
wafer_match_list <- read.csv("./pre_load_data/wafer_match_list.csv", row.names = 1)

chr_name <- unlist(read.table("./pre_load_data/chr_name.txt"), use.names = F)
chr_color <- read.table("./pre_load_data/chr_color.txt")
gene_bin <- read.table("./pre_load_data/hg19_gene_binpos_map.tsv", header = 1)
new_pal = c("#CC0C00B2", "#5C88DAB2", "#84BD00B2", "#FFCD00B2", "#7C878EB2", "#00B5E2B2", "#00AF66B2", "#D2AF81B2",
                      "#FD7446B2", "#46732EB2", "#C1395E", "#E07B42","#D4A2D9", "#8E72D5","#C0EDB9", "#364E4F", "#8EE5EE",
                      "#FFA500", "#458B00", "#CD6090", "#FFAEB9", "#90EE90", "#5f9EA0", "#E6E6FA", "#8B7E66", "#a6cee3","#1f78b4",
                      "#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928",
                      "#3d348b","#f72585","#7678ed","#9b5de5","#2ec4b6", "#fff1d0")
new_pal2 <- new_pal
names(new_pal2) <- paste0("c", 1:length(new_pal2))
superclone_col <- yarrr::piratepal("appletv")
names(superclone_col) <- NULL
hbca_epi_markers <- read.csv("./pre_load_data/HBCA_Epi_markers.csv")

pro_name <- "bcis13t"
pro_name1 <- "bcis13t_chip1"
pro_name2 <- "bcis13t_chip2"
disp_file1 <- c("/volumes/lab/users/wet_lab/instruments/wafergen/Rui_WDR/filter_files/20210913_BCMDCIS13T_Chip1_chip2/133466-2_Chip1_WellList.TXT")
disp_file2 <- c("/volumes/lab/users/wet_lab/instruments/wafergen/Rui_WDR/filter_files/20210913_BCMDCIS13T_Chip1_chip2/133484-2_Chip2_WellList.TXT")

# file.copy(disp_file1, paste0("./metrics/dispense_files/", pro_name1, "_", basename(disp_file1)))
# file.copy(disp_file2, paste0("./metrics/dispense_files/", pro_name2, "_", basename(disp_file2)))

disp1 <- read.delim(paste0("./metrics/dispense_files/", pro_name1, "_", basename(disp_file1)))
disp2 <- read.delim(paste0("./metrics/dispense_files/", pro_name2, "_", basename(disp_file2)))
dim(disp1);dim(disp2);dim(disp3);dim(disp4);dim(disp5);

pro_name_r <- paste0(pro_name, "_rna")
pro_name_d <- paste0(pro_name, "_dna")

###---RNA--------##########
pro_name_r <- paste0(pro_name, "_rna")

myfolder_name_p1 <- "20210912_WDR_bcis13TDFR_ECIS36T-2R"
mysample_name_p1 <- "bcis13T_Chip1_RNA"
raw_path_r1 <- c(paste0("/volumes/USR2/wangkl/wscDR/", myfolder_name_p1, "/data/", 
                        mysample_name_p1, "_exact_hg19_star_smart_2/star_2.7.5/",
                        mysample_name_p1, "_exact_hg19_star_smart_2_Solo.out/GeneFull/filtered_nodup"))
myfolder_name_p2 <- "20221020_WDR_RNAlibs_reseq2_dropDR"
mysample_name_p2 <- "bcis13T_Tube2_cDNA_Chip1_merge"
raw_path_r2 <- c(paste0("/volumes/USR2/wangkl/wscDR/", myfolder_name_p2, "/data/", 
                        mysample_name_p2, "_exact_hg19_star_smart_2/star_2.7.5/",
                        mysample_name_p2, "_exact_hg19_star_smart_2_Solo.out/GeneFull/filtered_nodup"))

myfolder_name_p3 <- "20221005_WDR_RNAlibs_reseq"
mysample_name_p3 <- "bcis13T_Tube2_cDNA_Chip2_merge"
raw_path_r3 <- c(paste0("/volumes/USR2/wangkl/wscDR/", myfolder_name_p3, "/data/", 
                        mysample_name_p3, "_exact_hg19_star_smart_2/star_2.7.5/",
                        mysample_name_p3, "_exact_hg19_star_smart_2_Solo.out/GeneFull/filtered_nodup"))

#-----RNA processing-----
mtx_p1<- Read10X(data.dir = raw_path_r1)
mtx_p1 <- CreateSeuratObject(counts = mtx_p1, min.cells = 1, min.features = 3, project = pro_name1)
names_tmp1 <- mtx_p1@meta.data %>% rownames_to_column()
filt_cells_names_p1 <- disp1 %>% dplyr::filter(!(Sample %in% c("Neg Ctrl", "Pos Ctrl"))) %>%
  dplyr::select(1,2,4) %>% mutate(row_col=paste(Row, Col, sep = "_")) %>% 
  dplyr::select(3:4) %>% left_join(wafer_match_list, by = c("row_col" = "coor")) %>% 
  inner_join(names_tmp1, by = c("RNA_Barcode" = "rowname")) %>% pull(RNA_Barcode)
dim(names_tmp1)
length(filt_cells_names_p1)
mtx_p1 <- subset(mtx_p1, cells = filt_cells_names_p1)

mtx_p2<- Read10X(data.dir = raw_path_r2)
mtx_p2 <- CreateSeuratObject(counts = mtx_p2, min.cells = 1, min.features = 3, project = pro_name2)
names_tmp2 <- mtx_p2@meta.data %>% rownames_to_column()
filt_cells_names_p2 <- disp2 %>% dplyr::filter(!(Sample %in% c("Neg Ctrl", "Pos Ctrl"))) %>% 
  dplyr::select(1,2,4) %>% mutate(row_col=paste(Row, Col, sep = "_")) %>% 
  dplyr::select(3:4) %>% left_join(wafer_match_list, by = c("row_col" = "coor")) %>% 
  inner_join(names_tmp2, by = c("RNA_Barcode" = "rowname")) %>% pull(RNA_Barcode)
dim(names_tmp2)
length(filt_cells_names_p2)
mtx_p2 <- subset(mtx_p2, cells = filt_cells_names_p2)

mtx_p3<- Read10X(data.dir = raw_path_r3)
mtx_p3 <- CreateSeuratObject(counts = mtx_p3, min.cells = 1, min.features = 3, project = pro_name3)
names_tmp3 <- mtx_p3@meta.data %>% rownames_to_column()
filt_cells_names_p3 <- disp3 %>% dplyr::filter(!(Sample %in% c("Neg Ctrl", "Pos Ctrl"))) %>%
  dplyr::select(1,2,4) %>% mutate(row_col=paste(Row, Col, sep = "_")) %>% 
  dplyr::select(3:4) %>% left_join(wafer_match_list, by = c("row_col" = "coor")) %>% 
  inner_join(names_tmp3, by = c("RNA_Barcode" = "rowname")) %>% pull(RNA_Barcode)
dim(names_tmp3)
length(filt_cells_names_p3)
mtx_p3 <- subset(mtx_p3, cells = filt_cells_names_p3)

mtx2 <- merge(mtx_p1, y = c(mtx_p2, mtx_p3), add.cell.ids = c(pro_name1,pro_name2,pro_name3), project = pro_name_r)
mtx2[["percent.mt"]] <- PercentageFeatureSet(mtx2, pattern = "^MT-")

pdf(paste0(c("./figures/"),pro_name_r, c("_singlet_Seurat_QC1.pdf")), height = 5, width = 10)
VlnPlot(mtx2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p1<-ggplot(mtx2@meta.data, aes(x=orig.ident, y=log(nFeature_RNA)))+geom_violin(fill="coral")+ylab("log(Gene count)")+xlab("")+
  geom_jitter(size=0.5)+theme(axis.text.x = element_text(size=10), legend.position = "none")+theme(axis.text.x = element_blank())
p2<-ggplot(mtx2@meta.data, aes(x=orig.ident, y=log(nCount_RNA)))+geom_violin(fill="coral")+ylab("log(Transcript count)")+xlab("")+
  geom_jitter(size=0.5)+theme(axis.text.x = element_text(size=10), legend.position = "none")+theme(axis.text.x = element_blank())
p3<-ggplot(mtx2@meta.data, aes(x=orig.ident, y=percent.mt))+geom_violin(fill="coral")+ylab("percentage of MT")+xlab("")+
  geom_jitter(size=0.5)+theme(axis.text.x = element_text(size=10), legend.position = "none")+theme(axis.text.x = element_blank())
plot_grid(p1, p2, p3, align = "h", ncol = 3)
plot1 <- FeatureScatter(mtx2, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(legend.position = "top")
plot2 <- FeatureScatter(mtx2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ theme(legend.position = "top")
CombinePlots(plots = list(plot1, plot2))
dev.off()

mtx2 <- subset(mtx2, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 20 )
#---add RNA_index info---
mtx2@meta.data <- mtx2@meta.data %>% rownames_to_column() %>% mutate(N7xx_bc = stringr::str_sub(rowname, -16, -9), S5xx_bc = stringr::str_sub(rowname, -8, -1)) %>%
  column_to_rownames()
mtx2 <- SCTransform(mtx2, vars.to.regress = c("percent.mt","nCount_RNA","orig.ident"), verbose = FALSE)
# mtx2 <- NormalizeData(object = mtx2, normalization.method = "LogNormalize", scale.factor = 10000)
# mtx2 <- FindVariableFeatures(mtx2, selection.method = "vst")
all.genes <- rownames(mtx2)

# mtx2 <- RunFastMNN(object.list = SplitObject(mtx2, split.by = "orig.ident"))
# # Note: FastMNN don't use scaled data, but the scaled data will be used for DoHeatmap
# options(future.globals.maxSize= 40960 * 1024^2)
# plan(strategy = "multicore", workers = 20)
# mtx2 <- ScaleData(mtx2, features = all.genes, vars.to.regress = c("percent.mt","nCount_RNA","orig.ident", "S5xx_bc"))#, "percent.rp"))

#---add RNA_index info---
# mtx2@meta.data <- mtx2@meta.data %>% rownames_to_column() %>% mutate(N7xx_bc = stringr::str_sub(rowname, -16, -9), S5xx_bc = stringr::str_sub(rowname, -8, -1)) %>% 
#   column_to_rownames()
# mtx2 <- ScaleData(mtx2, features = all.genes, vars.to.regress = c("orig.ident", "S5xx_bc", "percent.mt","nCount_RNA"))#,"S.Score", "G2M.Score"))
mtx2 <- RunPCA(mtx2, features = all.genes)
# ElbowPlot(mtx2, ndims = 50)

# sel_dim <- 15
# sel_res <- 0.4
# mtx3_1 <- RunUMAP(mtx2, dims = 1:sel_dim, reduction = "mnn")
# mtx3_1_tsu <- FindNeighbors(mtx3_1, dims = 1:sel_dim, reduction = "mnn") %>% FindClusters(resolution = sel_res)

res <- c(0.2,0.4)
dim <- c(5,8,10,15)
pgroup <- c("orig.ident")
plists <- dim_res_finder2(res, dim, mtx2, pgroup)
plist_ru <- plot_grid(plotlist =plists[4:5],ncol = 1)
p_list<-plot_grid(plists[[3]],plist_ru, ncol = 1)
cowplot::ggsave2(paste0(c("./figures/"), pro_name_r, c("_tSNE_test.pdf")), plists[[1]], width = 6*(length(res)+1), height = 5*length(dim))
cowplot::ggsave2(paste0(c("./figures/"), pro_name_r, c("_UMAP_test.pdf")), plists[[2]], width = 6*(length(res)+1), height = 5*length(dim))
cowplot::ggsave2(paste0(c("./figures/"), pro_name_r, c("_top10_test.pdf")), p_list, width = 8*length(dim), height = 30, limitsize = F)

# 
sel_dim <- 8
sel_res <- 0.2
mtx3_1 <- FindNeighbors(mtx2, dims = 1:sel_dim)
mtx3_2 <- FindClusters(mtx3_1, resolution = sel_res)
mtx3_1_tsu <- RunUMAP(mtx3_2, dims = 1:sel_dim)
mtx3_1_tsu <- RunTSNE(mtx3_1_tsu, dims = 1:sel_dim)

p1 <- DimPlot(mtx3_1_tsu, reduction = "tsne", label = T)
p2 <- DimPlot(mtx3_1_tsu, reduction = "tsne", label = T, group.by = "orig.ident")
p3 <- DimPlot(mtx3_1_tsu, reduction = "umap", label = T) #+ scale_color_manual(values = my_color_pal)
p4 <- DimPlot(mtx3_1_tsu, reduction = "umap", label = T, group.by = "orig.ident")
p_list1 <- plot_grid(plotlist = list(p1,p2,p3,p4), ncol = 2)

cowplot::ggsave2(paste0("./figures/", pro_name_r, "_TSNE_UMAP_all_cells.pdf"), p_list1, width = 10, height = 8) 
# cowplot::ggsave2(paste0("./figures/", pro_name_r, "_UMAP_all_cells_by_plates.pdf"), p3, width = 5, height = 4) 

saveRDS(mtx3_1_tsu, paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
mtx3_1_tsu <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))

####---------Gene and UMI number check--
mtx3_1_tsu.markers <- FindAllMarkers(mtx3_1_tsu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mtx3_1_tsu.markers %>% dplyr::group_by(cluster) %>% top_n(2, avg_log2FC)
top10 <- mtx3_1_tsu.markers %>% dplyr::group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(mtx3_1_tsu, features = top10$gene)

typical_marker<-c("COL1A1","COL1A2","APMAP", "ADIPOQ", "PTPRC", "MSR1", "KRT14","KRT15","ESR1", "KRT19", "MMRN1","PECAM1", "VWF",
                  "SYNPO2", "ACTA2","CD2", "CD247","ABCC11", "ACSM1", "RP11-206M11.7", "LINC01088","HPGD","IL18R1","BLK","MS4A1",
                  "BANK1", "BCL11A","FCRL1","ERBB2","GRB7","SERHL2","CENPE", "CDC6", "CIT", "BRCA1","CCNE2")
typical_marker<-c("COL1A1","COL1A2","APMAP", "ADIPOQ", "PTPRC", "MSR1", "KRT14","KRT15","ESR1", "KRT19", "MMRN1", "VWF",
                  "SYNPO2", "ACTA2","CD2", "CD247","ABCC11", "RP11-206M11.7", "LINC01088","HPGD","IL18R1","BLK",
                  "BANK1", "BCL11A","FCRL1","ERBB2","GRB7","SERHL2","CENPE", "CDC6", "CIT", "BRCA1","CCNE2")
n_col <- 6
# dir.create("./fig")
pdf(paste0("./figures/", pro_name, "_typical_marker_feature_plot_2.pdf"), width = 3*n_col, height = 2.8*ceiling(length(typical_marker)/n_col))
FeaturePlot(object = mtx3_1_tsu, feature = typical_marker, cols = c("grey", "blue"), reduction = "umap",ncol = n_col)
dev.off()
pdf(paste0("./figures/", pro_name, "typical_marker_feature_plot_violin_2.pdf"), width = 3*n_col, height = 3*ceiling(length(typical_marker)/n_col))
VlnPlot(object = mtx3_1_tsu, feature = typical_marker, ncol = n_col)
dev.off()


pdf(paste0(c("./figures/"),pro_name_r, c("_singlet_Seurat_top10_Heatmap.pdf")), height = 10, width = 8)
DoHeatmap(mtx3_1_tsu, features = top10$gene)
dev.off()

pdf(paste0(c("./figures/"),pro_name_r, c("_singlet_Seurat_clst_gene_umi_number.pdf")), height = 10, width = 8)
VlnPlot(mtx3_1_tsu, features = c("nCount_RNA","nFeature_RNA"))
dev.off()

###---DNA--------#########
pro_name_d <- paste0(pro_name, "_dna")
raw_path1 <- paste0("./map_seg_output/", pro_name1)
raw_path2 <- paste0("./map_seg_output/", pro_name2)


#-----DNA processing (no filtering applied to all cells)------####
#---t1--
varbin_mtx1 <- readVarbin(raw_path1, remove_Y = TRUE)
varbin_mtx1 # 809
filt_cells_names_tmp1 <- rownames(varbin_mtx1@colData) %>% as.data.frame() %>% mutate(cell = stringr::str_extract(., "c\\d+")) %>% dplyr::rename(cell_name = ".")
filt_cells_names_dis1 <- disp1 %>% dplyr::filter(!(Sample %in% c("Neg Ctrl", "Pos Ctrl"))) %>%  dplyr::select(1,2,4) %>% mutate(row_col=paste(Row, Col, sep = "_")) %>% 
  dplyr::select(3:4) %>% left_join(wafer_match_list, by = c("row_col" = "coor")) %>% mutate(cell2 = tolower(Cell)) %>% 
  inner_join(filt_cells_names_tmp1, by = c("cell2" = "cell")) %>% pull(cell_name)

varbin_mtx_all1 <- varbin_mtx1[,filt_cells_names_dis1]
varbin_mtx_all1@colData # 559

#---t2--
varbin_mtx2 <- readVarbin(raw_path2, remove_Y = TRUE)
varbin_mtx2 # 1287
filt_cells_names_tmp2 <- rownames(varbin_mtx2@colData) %>% as.data.frame() %>% mutate(cell = stringr::str_extract(., "c\\d+")) %>% dplyr::rename(cell_name = ".")
filt_cells_names_dis2 <- disp2 %>% dplyr::filter(!(Sample %in% c("Neg Ctrl", "Pos Ctrl"))) %>%  dplyr::select(1,2,4) %>% mutate(row_col=paste(Row, Col, sep = "_")) %>% 
  dplyr::select(3:4) %>% left_join(wafer_match_list, by = c("row_col" = "coor")) %>% mutate(cell2 = tolower(Cell)) %>% 
  inner_join(filt_cells_names_tmp2, by = c("cell2" = "cell")) %>% pull(cell_name)

varbin_mtx_all2 <- varbin_mtx2[,filt_cells_names_dis2]
varbin_mtx_all2@colData # 838

#------merge no filtered datasets----#####
varbin_mtx_all1@colData$experiment <- pro_name1
varbin_mtx_all2@colData$experiment <- pro_name2
varbin_mtx_all <- cbind(varbin_mtx_all1, varbin_mtx_all2)

varbin_mtx_all_log <- logNorm(varbin_mtx_all, transform = "log2")
saveRDS(varbin_mtx_all_log, file = paste0("./objects/", pro_name_d, c("_without_filtered_all_copykit_rnameta.rds")))
# varbin_mtx_all_log <- readRDS(paste0("./objects/", pro_name_d, c("_without_filtered_all_copykit_rnameta.rds")))

#------Knn smooth-----#######
metadata(varbin_mtx_all_log)$genome <-  "hg19"
varbin_mtx_all_log_knn <- knnSmooth(varbin_mtx_all_log)
varbin_mtx_all_log_knn <- logNorm(varbin_mtx_all_log_knn)
saveRDS(varbin_mtx_all_log_knn, file = paste0("./objects/", pro_name_d, c("_filtered_copykit_rnameta_knn.rds")))
# varbin_mtx_all_log_knn <- readRDS(paste0("./objects/", pro_name_d, c("_filtered_copykit_rnameta_knn.rds")))

#----Find normal cells---#####
varbin_mtx_all_log_knn2 <- findNormalCell(varbin_mtx_all_log_knn)
my_normal <- varbin_mtx_all_log_knn2@colData %>% as.data.frame() %>% dplyr::filter(is_normal)
my_tumor <- varbin_mtx_all_log_knn2@colData %>% as.data.frame() %>% dplyr::filter(!is_normal)
write.csv(my_normal, paste0("./metrics/", pro_name_d, "_knn_non_aneuplod_cells.csv"), row.names = F)
write.csv(my_tumor, paste0("./metrics/", pro_name_d, "_knn_aneuplod_cells.csv"), row.names = F)

# my_normal <- read.csv(paste0("./metrics/", pro_name_d, "_knn_non_aneuplod_cells.csv"))
# my_tumor <- read.csv(paste0("./metrics/", pro_name_d, "_knn_aneuplod_cells.csv"))

my_normal_name <- my_normal %>% pull(sample) %>% janitor::make_clean_names()
my_tumor_name <- my_tumor %>% pull(sample) %>% janitor::make_clean_names()

varbin_mtx_normal <- varbin_mtx_all_log_knn[,my_normal_name]
varbin_mtx_tumor <- varbin_mtx_all_log_knn[,my_tumor_name]

nrow(varbin_mtx_all_log_knn2@colData) # 559
nrow(varbin_mtx_normal@colData)
nrow(varbin_mtx_tumor@colData)

# saveRDS(varbin_mtx_normal, file = paste0("objects/", pro_name_d, c("_knn_normal_copykit_before_rename.rds")))
varbin_mtx_normal <- readRDS(paste0("objects/", pro_name_d, c("_knn_normal_copykit_before_rename.rds")))

# saveRDS(varbin_mtx_tumor, file = paste0("objects/", pro_name_d, c("_knn_tumor_copykit_before_rename.rds")))
varbin_mtx_tumor <- readRDS(paste0("objects/", pro_name_d, c("_knn_tumor_copykit_before_rename.rds")))

###-----Normal cell processing-----#####
varbin_mtx_normal <- readRDS(paste0("objects/", pro_name_d, c("_knn_normal_copykit_before_rename.rds")))

#----------UMAP-----
varbin_mtx_normal_log <- varbin_mtx_normal

set.seed(31)
near_nb <- 20
umap_data <-  data.frame(uwot::umap(log(t(varbin_mtx_normal_log@assays@data$segment_ratios)),
                                    metric = "manhattan", min_dist = 0.1, spread = 3, n_neighbors = near_nb))

umap_data2 <- umap_data %>% dplyr::rename("UMAP_1" = "X1", "UMAP_2" = "X2")
varbin_mtx_normal_log@int_colData$reducedDims <- umap_data2
varbin_mtx_normal_log@colData <- cbind(varbin_mtx_normal_log@colData, umap_data2)

#----------Find clusters----
varbin_mtx_normal_log_temp <- varbin_mtx_normal_log
varbin_mtx_normal_log <- varbin_mtx_normal_log_temp

#---subclone---
set.seed(17)
hdb_data <- dbscan::hdbscan(umap_data2, minPts = 20)
subclones <- paste0("c", as.character(hdb_data$cluster))

varbin_mtx_normal_log@colData <- cbind(varbin_mtx_normal_log@colData, subclones)
num_clst <- length(unique(varbin_mtx_normal_log@colData$subclones))
if(length(new_pal) >= num_clst) {
  new_pal_tmp = new_pal
} else {
  new_pal_tmp = hue_pal()(num_clst)
}
#----remove outlier cluster c0 and cluster that have less 3 cells.
if("c0" %in% names(table(varbin_mtx_normal_log@colData$subclones))){
  varbin_mtx_normal_log@colData$subclones <- factor(varbin_mtx_normal_log@colData$subclones,
                                                    levels = paste0("c", 0:(length(table(subclones))-1)))
  
  p1 <- ggplot(as.data.frame(varbin_mtx_normal_log@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) +
    geom_point(shape = 21, size=2.5, stroke = 0.03) + geom_text(label=subclones, check_overlap = T) +
    scale_fill_manual(values = c("grey",new_pal_tmp)) + theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
  p1
  cowplot::ggsave2(paste0("./figures/", pro_name_d, "_passDNAQC_knn_normal_cells_umap.pdf"), p1, width = 5, height = 4)
  
  clone_num <- table(varbin_mtx_normal_log@colData$subclones)
  kept_clones <- as.data.frame(clone_num) %>% dplyr::filter(Var1 != "c0") %>% dplyr::filter(Freq > 3) 
  
  varbin_mtx_normal_log <- varbin_mtx_normal_log[, subclones %in% kept_clones$Var1]
  varbin_mtx_normal_log@colData
}else{
  varbin_mtx_normal_log@colData$subclones <- factor(varbin_mtx_normal_log@colData$subclones,
                                                    levels = paste0("c", 1:(length(table(subclones)))))
}


p1 <- ggplot(as.data.frame(varbin_mtx_normal_log@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) +
  geom_point(shape = 21, size=2.5, stroke = 0.03) +
  scale_fill_manual(values = new_pal_tmp) + theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_final_filtered_knn_normal_cells_umap.pdf"), p1, width = 5, height = 4)

p2 <- ggplot(as.data.frame(varbin_mtx_normal_log@colData),aes(x = UMAP_1, y = UMAP_2, fill = experiment)) +
  geom_point(shape = 21, size=2.5, stroke = 0.03) +
  # scale_fill_manual(values = tp_col) + theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p2
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_final_filtered_knn_normal_cells_umap_plate.pdf"), p2, width = 5, height = 4)

saveRDS(varbin_mtx_normal_log, file = paste0("objects/", pro_name_d, c("_knn_normal_copykit_subclones.rds")))
# varbin_mtx_normal_log <- readRDS(paste0("objects/", pro_name_d, c("_knn_normal_copykit_subclones.rds")))

#----------Heatmap with clonal status bar----#####
# varbin_mtx_normal_log <- varbin_mtx_normal_log
mtx_srt <- as.data.frame(varbin_mtx_normal_log@colData) %>%
  arrange(factor(subclones, levels = paste0("c", 1:length(unique(varbin_mtx_normal_log@colData$subclones)))))
ht_mtx <- log2(t(varbin_mtx_normal_log@assays@data$segment_ratios))[mtx_srt$sample,]
#----annotation bar--
anno_mtx <- mtx_srt %>% dplyr::select(c("subclones", "experiment"))
rownames(anno_mtx) <- NULL

clst_col <- new_pal_tmp[1:length(unique(anno_mtx$subclones))]
names(clst_col) <- paste0("c", 1:length(unique(anno_mtx$subclones)))

exp_col <- jcolors::jcolors("pal2")[(6-length(unique(anno_mtx$experiment))):5]
names(exp_col) <- c(pro_name1, pro_name2)
#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color,
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue3", "white", "firebrick3"))
ha_row2=rowAnnotation(df = anno_mtx, col = list(subclones=clst_col, experiment=exp_col),
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/", pro_name_d, "_knn_normal_complexHeatmap_DR.pdf"), height = 10, width = 10)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = F, show_row_dend = F,
        row_split = anno_mtx$subclones,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
        # column_title = paste0(pro_name, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff),
        use_raster = T, raster_quality = 5, col = col_vec,
        heatmap_legend_param = list(title = "Log2 (Ratio)",
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row2)
dev.off()


#----------Convert ratio into integer---#####
varbin_mtx_normal_log2 <- varbin_mtx_normal_log
ploidy_pre <- 2
confidence_ratio <- 1
varbin_mtx_normal_log2@colData$ploidy_pre <- ploidy_pre
varbin_mtx_normal_log2@colData$confidence_ratio <- confidence_ratio

varbin_mtx_normal_log3 <- calcInteger(varbin_mtx_normal_log2, assay = "segment_ratios", method = "fixed", 
                                      ploidy_value = colData(varbin_mtx_normal_log2)$ploidy_pre)
varbin_mtx_normal_log3 <- calcConsensus(varbin_mtx_normal_log3, assay = "integer", fun="median", consensus_by = "subclones")

my_sub_clst <- names(varbin_mtx_normal_log3@consensus)
pg_clst <- c(my_sub_clst[apply(varbin_mtx_normal_log3@consensus, 2, max) >=3],
             my_sub_clst[apply(varbin_mtx_normal_log3@consensus, 2, min) <=1])
norm_clst <- setdiff(my_sub_clst, pg_clst)

norm_clst;pg_clst

#----------Merge normal subclones-----######
pg_clst2 <- unique(pg_clst)
my_sub <- varbin_mtx_normal_log2@colData$subclones %>% as.character() %>% as.data.frame()
colnames(my_sub) <- "mysub"
my_sub2 <- my_sub %>% mutate(mysub2 = ifelse(mysub %in% norm_clst, "nc1", plyr::mapvalues(mysub, from = pg_clst2, to = paste0("nc", 2:(length(pg_clst2)+1))))) %>% 
  mutate(mysub2 = factor(mysub2,levels=paste0("nc", 1:(length(pg_clst2)+1))))
varbin_mtx_normal_log2@colData$subclones2 <- my_sub2$mysub2

saveRDS(varbin_mtx_normal_log2, file = paste0("objects/", pro_name_d, c("_knn_normal_copykit_subclones_integer_merged_subclones.rds")))
# varbin_mtx_normal_log2 <- readRDS(paste0("objects/", pro_name_d, c("_knn_normal_copykit_subclones_integer_merged_subclones.rds")))

#----------Heatmap of merged subclones2 with clonal status bar----#####
varbin_mtx_normal_log <- varbin_mtx_normal_log2
mtx_srt <- as.data.frame(varbin_mtx_normal_log@colData) %>%
  arrange(factor(subclones2, levels = paste0("nc", 1:length(unique(varbin_mtx_normal_log@colData$subclones2)))))
ht_mtx <- log2(t(varbin_mtx_normal_log@assays@data$segment_ratios))[mtx_srt$sample,]
#----annotation bar--
anno_mtx <- mtx_srt %>% dplyr::select(c("subclones2", "experiment"))
rownames(anno_mtx) <- NULL

clst_col <- new_pal[1:length(unique(anno_mtx$subclones2))]
names(clst_col) <- paste0("nc", 1:length(unique(anno_mtx$subclones2)))

exp_col <- jcolors::jcolors("pal2")[(6-length(unique(anno_mtx$experiment))):5]
names(exp_col) <- c(pro_name1, pro_name2)
#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color,
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue3", "white", "firebrick3"))
ha_row2=rowAnnotation(df = anno_mtx, col = list(subclones2=clst_col, experiment=exp_col),
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/", pro_name_d, "_knn_normal_complexHeatmap_DR_merged_subclones.pdf"), height = 10, width = 10)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = F, show_row_dend = F,
        row_split = anno_mtx$subclones2,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
        # column_title = paste0(pro_name, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff),
        use_raster = T, raster_quality = 5, col = col_vec,
        heatmap_legend_param = list(title = "Log2 (Ratio)",
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row2)
dev.off()

###-----Tumor cell processing----#####
varbin_mtx_tumor <- readRDS(paste0("objects/", pro_name_d, c("_knn_tumor_copykit_before_rename.rds")))
#----------UMAP-----
varbin_mtx_tumor_log <- varbin_mtx_tumor

set.seed(31)
near_nb <- 20
umap_data <-  data.frame(uwot::umap(log(t(varbin_mtx_tumor_log@assays@data$segment_ratios)), 
                                    metric = "manhattan", min_dist = 0.1, spread = 3, n_neighbors = near_nb))

umap_data2 <- umap_data %>% dplyr::rename("UMAP_1" = "X1", "UMAP_2" = "X2")
varbin_mtx_tumor_log@int_colData$reducedDims <- umap_data2
varbin_mtx_tumor_log@colData <- cbind(varbin_mtx_tumor_log@colData, umap_data2)

#----------Find clusters----
#---subclone---
set.seed(17)
hdb_data <- dbscan::hdbscan(umap_data2, minPts = 20)
subclones <- paste0("c", as.character(hdb_data$cluster))

varbin_mtx_tumor_log@colData <- cbind(varbin_mtx_tumor_log@colData, subclones)
num_clst <- length(unique(varbin_mtx_tumor_log@colData$subclones))
if(length(new_pal) >= num_clst) {
  new_pal_tmp = new_pal
} else {
  new_pal_tmp = hue_pal()(num_clst)
}

#----remove outlier cluster c0 and cluster that have less 3 cells.
if("c0" %in% names(table(varbin_mtx_tumor_log@colData$subclones))){
  varbin_mtx_tumor_log@colData$subclones <- factor(varbin_mtx_tumor_log@colData$subclones,
                                                   levels = paste0("c", 0:(length(table(subclones))-1)))
  
  p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) +
    geom_point(shape = 21, size=2.5, stroke = 0.03) + geom_text(label=subclones, check_overlap = T) +
    scale_fill_manual(values = c("grey",new_pal_tmp)) + theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
  p1
  cowplot::ggsave2(paste0("./figures/", pro_name_d, "_passDNAQC_knn_tumor_cells_umap.pdf"), p1, width = 5, height = 4)
  
  clone_num <- table(varbin_mtx_tumor_log@colData$subclones)
  kept_clones <- as.data.frame(clone_num) %>% dplyr::filter(Var1 != "c0") %>% dplyr::filter(Freq > 3) 
  
  varbin_mtx_tumor_log <- varbin_mtx_tumor_log[, subclones %in% kept_clones$Var1]
  varbin_mtx_tumor_log@colData
}else{
  varbin_mtx_tumor_log@colData$subclones <- factor(varbin_mtx_tumor_log@colData$subclones,
                                                   levels = paste0("c", 1:(length(table(subclones)))))
}


p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) +
  geom_point(shape = 21, size=2.5, stroke = 0.03) +
  scale_fill_manual(values = new_pal_tmp) + theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_final_filtered_knn_tumor_cells_umap.pdf"), p1, width = 5, height = 4)

p2 <- ggplot(as.data.frame(varbin_mtx_tumor_log@colData),aes(x = UMAP_1, y = UMAP_2, fill = experiment)) +
  geom_point(shape = 21, size=2.5, stroke = 0.03) +
  # scale_fill_manual(values = tp_col) + theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p2
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_final_filtered_knn_tumor_cells_umap_plate.pdf"), p2, width = 5, height = 4)

saveRDS(varbin_mtx_tumor_log, file = paste0("objects/", pro_name_d, c("_knn_tumor_copykit_subclones.rds")))
# varbin_mtx_tumor_log <- readRDS(paste0("objects/", pro_name_d, c("_knn_tumor_copykit_subclones.rds")))

#----------Heatmap with clonal status bar----#####
mtx_srt <- as.data.frame(varbin_mtx_tumor_log@colData) %>% mutate(subclones = as.character(subclones)) %>% 
  mutate(subclones = factor(subclones, levels = paste0("c", 1:length(unique(varbin_mtx_tumor_log@colData$subclones))))) %>% 
  arrange(subclones)
ht_mtx <- log2(t(varbin_mtx_tumor_log@assays@data$segment_ratios))[mtx_srt$sample,]
#----annotation bar--
anno_mtx <- mtx_srt %>% dplyr::select(c("subclones", "experiment"))
rownames(anno_mtx) <- NULL
anno_mtx$subclones <- factor(anno_mtx$subclones, levels = paste0("c", 1:length(unique(anno_mtx$subclones))))

clst_col <- new_pal_tmp[1:length(unique(anno_mtx$subclones))]
names(clst_col) <- paste0("c", 1:length(unique(anno_mtx$subclones)))

exp_col <- jcolors::jcolors("pal2")[(6-length(unique(anno_mtx$experiment))):5]
names(exp_col) <- c(pro_name1, pro_name2)
#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color,
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue3", "white", "firebrick3"))
ha_row2=rowAnnotation(df = anno_mtx, col = list(subclones=clst_col, experiment=exp_col),
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/", pro_name_d, "_knn_tumor_complexHeatmap_before_remove_doublets_DR.pdf"), height = 10, width = 10)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = F, show_row_dend = F,
        row_split = anno_mtx$subclones,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
        # column_title = paste0(pro_name, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff),
        use_raster = T, raster_quality = 5, col = col_vec,
        heatmap_legend_param = list(title = "Log2 (Ratio)",
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row2)
dev.off()

#----------Calculate computational ploidy----####
varbin_mtx_tumor_log2 <- varbin_mtx_tumor_log
scq.call <- BiocParallel::bplapply(colnames(assay(varbin_mtx_tumor_log2, "bin_counts")),function(x){
  ## call scquantum, subset the scaffold data if chrY is removed
  a <- ploidy.inference(assay(varbin_mtx_tumor_log2, "bin_counts")[,x], hg19_rg$chr[1:12167], hg19_rg$start[1:12167], hg19_rg$end[1:12167])
  c(x,unlist(a[c(4,11)]))
})
scq_mat <- as_tibble(do.call(rbind, scq.call))
scq_mat2 <- as.data.frame(scq_mat) %>% mutate(ploidy_pre = as.numeric(ploidy), confidence_ratio = as.numeric(confidence_ratio)) %>% 
  dplyr::select("ploidy_pre","confidence_ratio")
varbin_mtx_tumor_log2@colData <- cbind(varbin_mtx_tumor_log2@colData, scq_mat2)

p2 <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% ggplot(aes(x=ploidy_pre, y = confidence_ratio)) + 
  geom_point(size=0.1)  + facet_grid(~subclones)
p2
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_scquantum_ploidy_confidence_by_subclone1.pdf"), p2, 
                 width = length(unique(varbin_mtx_tumor_log2@colData$subclones)), height = 2)

saveRDS(varbin_mtx_tumor_log2, file = paste0("objects/", pro_name_d, c("_tumor_before_remove_doublets_copykit.rds")))
# varbin_mtx_tumor_log2 <- readRDS(paste0("objects/", pro_name_d, c("_tumor_before_remove_doublets_copykit.rds")))

bad_sub <- c("c11","c12")
#----------Remove bad clusters----#######
my_singlets_name <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::filter(!(subclones %in% bad_sub)) %>% pull(sample)
varbin_mtx_tumor_log2 <- varbin_mtx_tumor_log2[,my_singlets_name]

saveRDS(varbin_mtx_tumor_log2, file = paste0("objects/", pro_name_d, c("_tumor_after_remove_doublets_copykit.rds")))
# varbin_mtx_tumor_log2 <- readRDS(paste0("objects/", pro_name_d, c("_tumor_after_remove_doublets_copykit.rds")))

#----------Calculate integer CN of remaining subclones----####
ploidy_sub <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% 
  dplyr::select("subclones", "ploidy_pre") %>% `rownames<-`( NULL ) %>% as_tibble() %>% dplyr::group_by(subclones) %>% dplyr::summarise(sub_ploidy = mean(ploidy_pre))
ploidy_sub
ploidy_sub <- ploidy_sub %>% mutate(sub_ploidy = ifelse(sub_ploidy<2, 2, sub_ploidy))

varbin_mtx_tumor_log3 <- varbin_mtx_tumor_log2
varbin_mtx_tumor_log3@colData$ploidy_subclone <- plyr::mapvalues(varbin_mtx_tumor_log3@colData$subclones, from = ploidy_sub$subclones, to = ploidy_sub$sub_ploidy) %>% 
  as.character() %>% as.numeric()

varbin_mtx_tumor_log3 <- calcInteger(varbin_mtx_tumor_log3, assay = "segment_ratios", method = "fixed", ploidy_value = varbin_mtx_tumor_log3@colData$ploidy_subclone)
varbin_mtx_tumor_log3 <- calcConsensus(varbin_mtx_tumor_log3, assay = "integer", fun="median", consensus_by = "subclones")

ploidy_trunc <- round(2*mean(varbin_mtx_tumor_log3@colData$ploidy_subclone))
eventmat <- getEventMat(varbin_mtx_tumor_log3, bin_adj = 2, ploidy_trunc = ploidy_trunc)

## convert back to segmentation level
popseg_long <- as.data.frame(apply(as.data.frame(t(eventmat %>% dplyr::select(matches("c[0-9]+")))), 1, function(m) {rep.int(m, eventmat$n.bins)}))

attr(popseg_long, "consensus_by") <- "subclones"
attr(popseg_long, "consensus_assay") <- "integer"
copykit::consensus(varbin_mtx_tumor_log3) <- popseg_long
varbin_mtx_tumor_log3 <- runConsensusPhylo(varbin_mtx_tumor_log3)

saveRDS(varbin_mtx_tumor_log3, file = paste0("objects/", pro_name_d, c("_tumor_after_remove_doublets_integerCN_copykit.rds")))
# varbin_mtx_tumor_log3 <- readRDS(paste0("objects/", pro_name_d, c("_tumor_after_remove_doublets_integerCN_copykit.rds")))

#----------Identify subclones with same genetype----#####
#-biggest segment difference of the consensus profiles of two subclones has less than 10 bins
con_tmp <- varbin_mtx_tumor_log3@consensus
con_tmp0 <- con_tmp %>% rownames_to_column()
ref_tem <- varbin_mtx_tumor_log3@rowRanges %>% as.data.frame() %>% rownames_to_column() %>% dplyr::select(1:2)

clst_comp <- data.frame()
for (i in 1:(ncol(con_tmp)-1)) {
  for (j in (i+1):ncol(con_tmp)) {
    diff_loc <- con_tmp0$rowname[con_tmp[,i] != con_tmp[,j]]
    my_sum <- ref_tem %>% dplyr::filter(rowname %in% diff_loc) %>% pull(seqnames) %>% table()
    my_sum2 <- my_sum[my_sum>0]
    if (max(my_sum2) < 10) {
      tmp <- paste(colnames(con_tmp)[i],colnames(con_tmp)[j], sep = ",")
      clst_comp <- rbind(clst_comp, tmp)
    }
  }
}
clst_comp

to_clone1 <- c("c2","c9","c15")
from_clone1 <- c("c3","c10","c16")

subclones2_tmp <- plyr::mapvalues(varbin_mtx_tumor_log3@colData$subclones, 
                                  from = from_clone1, 
                                  to = to_clone1)

subclones2_uni <- unique(subclones2_tmp) %>% sort() %>% as.character()
subclones2_tmp2 <- plyr::mapvalues(subclones2_tmp, from = subclones2_uni, to = paste0("tc", 1:length(subclones2_uni))) %>% as.character() %>% 
  factor(levels = paste0("tc", 1:length(subclones2_uni)))
varbin_mtx_tumor_log2@colData$subclones2 <- subclones2_tmp2

saveRDS(varbin_mtx_tumor_log2, file = paste0("objects/", pro_name_d, c("_tumor_merged_all_subclones_copykit.rds")))
# varbin_mtx_tumor_log2 <- readRDS(paste0("objects/", pro_name_d, c("_tumor_merged_all_subclones_copykit.rds")))

p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log2@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones2)) +
  geom_point(shape = 21, size=2.5, stroke = 0.03) +
  scale_fill_manual(values = new_pal) + theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_final_filtered_knn_tumor_cells_umap_after_remove_doublets.pdf"), p1, width = 5, height = 4)

#----------Heatmap with clonal status bar----#####
mtx_srt <- as.data.frame(varbin_mtx_tumor_log2@colData) %>%
  arrange(factor(subclones2, levels = paste0("tc", 1:length(unique(varbin_mtx_tumor_log2@colData$subclones2)))))
ht_mtx <- log2(t(varbin_mtx_tumor_log2@assays@data$segment_ratios))[mtx_srt$sample,]
#----annotation bar--
anno_mtx <- mtx_srt %>% dplyr::select(c("subclones2", "experiment"))
rownames(anno_mtx) <- NULL

clst_col <- new_pal[1:length(unique(anno_mtx$subclones2))]
names(clst_col) <- paste0("tc", 1:length(unique(anno_mtx$subclones2)))

exp_col <- jcolors::jcolors("pal2")[(6-length(unique(anno_mtx$experiment))):5]
names(exp_col) <- c(pro_name1, pro_name2)
#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color,
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue3", "white", "firebrick3"))
ha_row2=rowAnnotation(df = anno_mtx, col = list(subclones2=clst_col, experiment=exp_col),
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/", pro_name_d, "_knn_tumor_complexHeatmap_DR_after_remove_doublets.pdf"), height = 10, width = 10)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = F, show_row_dend = F,
        row_split = anno_mtx$subclones2,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
        # column_title = paste0(pro_name, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff),
        use_raster = T, raster_quality = 5, col = col_vec,
        heatmap_legend_param = list(title = "Log2 (Ratio)",
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row2)
dev.off()


###-----Merge normal and tumor cells and plot Heatmaps---######
varbin_mtx_all_log <- cbind(varbin_mtx_normal_log2, varbin_mtx_tumor_log2)
saveRDS(varbin_mtx_all_log, file = paste0("objects/", pro_name_d, c("_final_merged_before_rename_copykit.rds")))
# varbin_mtx_all_log <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_before_rename_copykit.rds")))

#---rename subclones---
varbin_mtx_all_log@colData$subclones3 <- plyr::mapvalues(varbin_mtx_all_log@colData$subclones2, 
                                                         from = levels(varbin_mtx_all_log@colData$subclones2), 
                                                         to = paste0("c", 1:length(levels(varbin_mtx_all_log@colData$subclones2))))

ploidy_sub <- varbin_mtx_all_log@colData %>% as.data.frame() %>% 
  dplyr::select("subclones3", "ploidy_pre") %>% `rownames<-`( NULL ) %>% as_tibble() %>% dplyr::group_by(subclones3) %>% 
  dplyr::summarise(sub_ploidy = median(ploidy_pre))
ploidy_sub
ploidy_sub <- ploidy_sub %>% mutate(sub_ploidy = ifelse(sub_ploidy<2, 2, sub_ploidy))

varbin_mtx_all_log@colData$ploidy_subclone <- plyr::mapvalues(varbin_mtx_all_log@colData$subclones3, from = ploidy_sub$subclones3, to = ploidy_sub$sub_ploidy) %>% 
  as.character() %>% as.numeric()

varbin_mtx_all_log2 <- calcInteger(varbin_mtx_all_log, assay = "segment_ratios", method = "fixed", ploidy_value = varbin_mtx_all_log@colData$ploidy_subclone)
varbin_mtx_all_log2 <- calcConsensus(varbin_mtx_all_log2, assay = "integer", fun="median", consensus_by = "subclones3")

ploidy_trunc <- round(2*mean(varbin_mtx_all_log2@colData$ploidy_subclone))
eventmat <- getEventMat(varbin_mtx_all_log2, bin_adj = 2, ploidy_trunc = ploidy_trunc)

## convert back to segmentation level
popseg_long <- as.data.frame(apply(as.data.frame(t(eventmat %>% dplyr::select(matches("c[0-9]+")))), 1, function(m) {rep.int(m, eventmat$n.bins)}))

attr(popseg_long, "consensus_by") <- "subclones3"
attr(popseg_long, "consensus_assay") <- "integer"
copykit::consensus(varbin_mtx_all_log2) <- popseg_long
varbin_mtx_all_log2 <- runConsensusPhylo(varbin_mtx_all_log2)

saveRDS(varbin_mtx_all_log2, file = paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))
# varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))

#----------Integer single cell Heatmap ----#####
# varbin_mtx_all_log2 <- varbin_mtx_all_log2
mtx_srt <- as.data.frame(varbin_mtx_all_log2@colData) %>%
  arrange(factor(subclones3, levels = levels(varbin_mtx_all_log2@colData$subclones3)))
ht_mtx <- t(varbin_mtx_all_log2@assays@data$integer)[mtx_srt$sample,]
#----annotation bar--
anno_mtx <- mtx_srt %>% dplyr::select(c("subclones3", "experiment"))
rownames(anno_mtx) <- NULL

clst_col <- new_pal[1:length(unique(anno_mtx$subclones3))]
names(clst_col) <- levels(varbin_mtx_all_log2@colData$subclones3)

exp_col <- jcolors::jcolors("pal2")[(6-length(unique(anno_mtx$experiment))):5]
names(exp_col) <- c(pro_name1, pro_name2)
#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color,
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

heat_col <- pals::ocean.balance(14)[9:14]
if(max(ht_mtx)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(ht_mtx)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(ht_mtx)-6))
}
names(col_vec) <- 0:max(ht_mtx)

ha_row2=rowAnnotation(df = anno_mtx, col = list(subclones3=clst_col, experiment=exp_col),
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/", pro_name_d, "_merged_integer_complexHeatmap_DR.pdf"), height = 10, width = 10)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = F, show_row_dend = F,
        row_split = anno_mtx$subclones3,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
        # column_title = paste0(pro_name, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff),
        use_raster = T, raster_quality = 5, col = col_vec,
        heatmap_legend_param = list(title = "Integer CNA",
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row2)
dev.off()



#----------Ratio single cell Heatmap----#####
mtx_srt <- as.data.frame(varbin_mtx_all_log2@colData) %>%
  arrange(factor(subclones3, levels = levels(varbin_mtx_all_log2@colData$subclones3)))
ht_mtx <- log2(t(varbin_mtx_all_log2@assays@data$segment_ratios))[mtx_srt$sample,]
#----annotation bar--
anno_mtx <- mtx_srt %>% dplyr::select(c("subclones3", "experiment"))
rownames(anno_mtx) <- NULL

clst_col <- new_pal[1:length(unique(anno_mtx$subclones3))]
names(clst_col) <- levels(varbin_mtx_all_log2@colData$subclones3)

exp_col <- jcolors::jcolors("pal2")[(6-length(unique(anno_mtx$experiment))):5]
names(exp_col) <- c(pro_name1, pro_name2)
#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color,
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue3", "white", "firebrick3"))
ha_row2=rowAnnotation(df = anno_mtx, col = list(subclones3=clst_col, experiment=exp_col),
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/", pro_name_d, "_merged_log2_ratio_complexHeatmap_DR.pdf"), height = 10, width = 10)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = F, show_row_dend = F,
        row_split = anno_mtx$subclones3,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
        # column_title = paste0(pro_name, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff),
        use_raster = T, raster_quality = 5, col = col_vec,
        heatmap_legend_param = list(title = "Log2 (Ratio)",
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row2)
dev.off()


#----------Integer consensus Heatmap----#####
cs_mtx <- t(varbin_mtx_all_log2@consensus)
my_order <- levels(varbin_mtx_all_log2@colData$subclones3)
cs_mtx_order <- cs_mtx[my_order,]
topleft(cs_mtx_order)

clst_col_cs <- new_pal[1:length(unique(varbin_mtx_all_log2@colData$subclones3))]
subclone_cs <- levels(varbin_mtx_all_log2@colData$subclones3)
names(clst_col_cs) <- subclone_cs
ha_row_cs=rowAnnotation(subclones = my_order, col = list(subclones = clst_col_cs), show_annotation_name = F, show_legend = F)

#-----header
b <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/chr_a.txt")
a <- unlist(b, use.names = F)
chr2 <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/chr_chr2.txt")
ha_col=HeatmapAnnotation(foo=anno_text(a, rot = 0, gp = gpar(fontsize =10)), df =chr2, col = list(chr=c("1"="black", "2"="grey")), 
                         show_legend = F, annotation_name_side = "left")

#-----annotated genes---
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
lab45_uniq_ext <- c("CKDN2C","AKT3","CKS1B","RUVBL1","C8orf4","FGFR1","BAG4","MTDH","MYC","PAK1","CDK4","MDM2", "FHIT","NEDD9","FOXK2","CCNE2",
                    "PLA2G10","GRB7","RPS6KB1","PPM1D","CCNE1","YWHAB","ZNF217","AURKA","PTK6","PPP2R2A",
                    "CCND1","NCOA3","TP53","CDKN1A","GADD45A","MDM2","RB1","CDKN2A","CDKN2C","ESR1","SHC1",
                    "PGR","PIK3CA","PTEN","BCL2","GATA3","FGFR4","EGFR","ERBB2","EMSY","BRCA1","FOXA1","BRCA2")
#---from: Chromosome 16 tumor-suppressor genes in breast cancer --- https://doi.org/10.1002/gcc.20318
chr16_tumor_suppressor <- c("E2F4","CTCF","TERF2","FBXL8","LRRC29","CDH5","CDH16","CDH3","CDH11","ATBF1","BBC1","FANCA","GAS11","CPNE7","ZFP276",
                            "CBFA2T3","CYBA","WWOX","TERF21P","SIAH1") 
chr16_tumor_suppressor2 <- c("E2F4","CTCF","TERF2","FBXL8","CDH5","CDH16","CDH3","CDH11","GAS11","CPNE7","ZFP276","CBFA2T3","WWOX","TERF21P") 

gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% c(lab45_uniq_ext, chr16_tumor_suppressor2))
ha_bottom_cs = columnAnnotation(clonal_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
                                                         side = "bottom", labels_gp = gpar(fontsize = 6)))

heat_col <- pals::ocean.balance(14)[9:14]
if(max(cs_mtx)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(cs_mtx)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(cs_mtx)-6))
}
names(col_vec) <- 0:max(cs_mtx)

pdf(paste0("./figures/", pro_name_d, "_merged_consensus_integer_complexHeatmap_DR2.pdf"), height = 4, width = 10)
Heatmap(as.matrix(cs_mtx_order), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        name = "csheatmap", show_row_names = T, row_names_side = "left",
        show_column_names = F, column_title = paste0(pro_name_d, "_consensusHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row_cs, bottom_annotation = ha_bottom_cs)
dev.off()

#----------check small cluster ratio map---#####
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))
small_clst <- c("c2","c3")
chip_dna_path1 <- c("/volumes/USR2/wangkl/wscDR/20221125_WDR_DNAlib_reseq_BCIS13T_Tube1_Chip1/data/waferDR_BCIS13T_tube1_chip1_DNA/waferDR_BCIS13_tube1_chip1-DNA_merged_200/res_200_k")
chip_dna_path2 <- c("/volumes/USR2/wangkl/wscDR/20221126_WDR_DNAlib_reseq_BCIS13T_Tube1_Chip2/data/waferDR_BCIS13T_Tube1_Chip2_DNA/waferDR_BCIS13T_Tube1_Chip2-DNA_200/res_200_k")

for (i in small_clst) {
  cell_mtx <- varbin_mtx_all_log2@colData %>% as.data.frame() %>% dplyr::filter(subclones3 == i) %>% dplyr::select("experiment","subclones3") %>% 
    rownames_to_column() %>% mutate(cell_name = toupper(stringr::str_extract(rowname, "_c\\d+_")))
  
  p_list <- list()
  j=1
  
  p_list <- apply(cell_mtx, 1, FUN = function(x){
    if (stringr::str_detect(x["experiment"], "chip1")) {
      mypath <- chip_dna_path1
    } else if (stringr::str_detect(x["experiment"], "chip2")){
      mypath <- chip_dna_path2
    } else if (stringr::str_detect(x["experiment"], "chip3")){
      mypath <- chip_dna_path3
    }
    mypng <- fs::dir_ls(path = mypath,recurse = T, glob = paste0("*", x["cell_name"], "*.png"))
    img <- readPNG(mypng)
    g <- rasterGrob(img, interpolate=TRUE) 
    p1 <- ggplot()+ annotation_custom(g)
    return(p1)
  })
  
  p_list <- plot_grid(plotlist = p_list, ncol = 6)
  cowplot::ggsave2(paste0("figures/", pro_name_d, "_small_clusters_", i, "_ratio_plots.pdf"), p_list, 
                   width = 6*6, height = nrow(cell_mtx)/6*4, limitsize = F)
}



#----------MEDICC2 Tree----######
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))
## make medicc input matrix
medicc_input <- cbind(varbin_mtx_all_log2@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_all_log2@consensus) %>%
  mutate(diploid=c1) %>%          ## assign diploid cell CN profile to be used as root in the medicc tree
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())

#---remove boarder regions (top and bottom 10 bin of each chr from each subclone)---##
medicc_input <- medicc_input %>%
  group_by(sample_id, chrom) %>%
  mutate(
    CN = case_when(
      row_number() %in% 1:10 ~ CN[11],  # For the first 10 rows, use the CN value of the 11th row
      row_number() %in% (n() - 9):n() ~ CN[n() - 10],  # For the last 10 rows, use the CN value of the 11th row from the bottom
      TRUE ~ CN  # For all other rows, keep the original CN value
    )
  )

dir.create(paste0("metrics/medicc_files/", pro_name_d),recursive = T)
write_tsv(medicc_input, file = paste0("metrics/medicc_files/", pro_name_d, "/",pro_name_d, c("_medicc2_input.tsv")))
# medicc_input <- read_tsv(file = paste0("metrics/medicc_files/", pro_name_d, "/",pro_name_d, c("_medicc2_input.tsv")))

mywd <- getwd()
#----go to terminal and run medicc2---
# conda activate medicc_env
# medicc2 -a CN --total-copy-numbers -j 40 -vv input_path/*_medicc2_input.tsv output_path

tsv_out_path <- paste0(mywd, "/metrics/medicc_files/", pro_name_d, "/")
medic <- ape::read.tree(paste0(tsv_out_path, pro_name_d,"_medicc2_input_final_tree.new"))

my_sub <- factor(levels(varbin_mtx_all_log2@colData$subclones3), levels = levels(varbin_mtx_all_log2@colData$subclones3))
list_samples <- split(my_sub, my_sub)

clst_col <- new_pal[1:length(unique(varbin_mtx_all_log2@colData$subclones3))]
names(clst_col) <- levels(varbin_mtx_all_log2@colData$subclones3)

tree <- ggtree::groupOTU(medic, list_samples)
treeplt <- ggtree::ggtree(ape::ladderize(tree), ladderize = FALSE, size = .2) +
  ggtree::geom_tiplab(size=6, aes(color=group),hjust = -0.6, alpha=1)+
  ggtree::geom_tippoint(size=6, aes(color=group), alpha=1)+
  scale_colour_manual(values = clst_col,breaks = names(clst_col)) +
  geom_text(aes(x=branch, label=plyr::mapvalues(round(branch.length, 0),from = "0",to = ""), vjust=-.5), size = 3) +
  theme(legend.position = "none") + ggtree::geom_rootpoint() 

cowplot::ggsave2(paste0("./figures/", pro_name_d, "_medicc2_tree_merged_subclones.pdf"), treeplt, width = 10, height = 6)

#----------MEDICC2 ME Tree---lumHR----######
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))
## make medicc input matrix
#--other clusters--
other_clst <- c("c2", "c3")
dip_clst <- c("c1")
medicc_input <- cbind(varbin_mtx_all_log2@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_all_log2@consensus) %>%
  mutate(diploid=2) %>%       ## assign diploid cell CN profile to be used as root in the medicc tree
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything()) %>% 
  dplyr::filter(!(sample_id %in% other_clst))

#---remove boarder regions (top and bottom 10 bin of each chr from each subclone)---##
medicc_input <- medicc_input %>%
  group_by(sample_id, chrom) %>%
  mutate(
    CN = case_when(
      row_number() %in% 1:10 ~ CN[11],  # For the first 10 rows, use the CN value of the 11th row
      row_number() %in% (n() - 9):n() ~ CN[n() - 10],  # For the last 10 rows, use the CN value of the 11th row from the bottom
      TRUE ~ CN  # For all other rows, keep the original CN value
    )
  )

dir.create(paste0("metrics/medicc_files_lumhr/", pro_name_d),recursive = T)
write_tsv(medicc_input, file = paste0("metrics/medicc_files_lumhr/", pro_name_d, "/",pro_name_d, c("_lumhr_medicc2_input.tsv")))
# medicc_input <- read_tsv(file = paste0("metrics/medicc_files_lumhr/", pro_name_d, "/",pro_name_d, c("_lumhr_medicc2_input.tsv")))

mywd <- getwd()
#----go to terminal and run medicc2---
# conda activate medicc_env
# medicc2 -a CN --total-copy-numbers -j 40 -vv input_path/*_lumhr_medicc2_input.tsv output_path
# eg: medicc2 -a CN --total-copy-numbers -j 40 -vv ./*_lumhr_medicc2_input.tsv ./
tsv_out_path <- paste0(mywd, "/metrics/medicc_files_lumhr/", pro_name_d, "/")
#---ME tree, need two diploid cell populations to make diploid branch on the most left side.
medic_dis <- read_tsv(paste0(tsv_out_path, pro_name_d, "_lumhr_medicc2_input_pairwise_distances.tsv"))
medic_dis2 <- medic_dis %>% as.data.frame() %>% column_to_rownames(var="sample_id") %>% as.dist()
me_tree <- ape::fastme.bal(medic_dis2)
me_tree <- ape::root.phylo(me_tree, outgroup = "diploid", resolve.root = TRUE)
me_tree <- ape::drop.tip(me_tree, tip = dip_clst)

my_sub <- factor(levels(varbin_mtx_all_log2@colData$subclones3), levels = levels(varbin_mtx_all_log2@colData$subclones3))
my_sub <- my_sub[!(my_sub %in% other_clst)] %>% droplevels()
list_samples <- split(my_sub, my_sub)
clst_col <- new_pal[1:length(unique(varbin_mtx_all_log2@colData$subclones3))]
names(clst_col) <- levels(varbin_mtx_all_log2@colData$subclones3)
clst_col <- clst_col[!(names(clst_col) %in% other_clst)]

tree <- ggtree::groupOTU(me_tree, list_samples)
treeplt <- ggtree::ggtree(ape::ladderize(tree), ladderize = FALSE, size = .2) +
  ggtree::geom_tiplab(size=6, aes(color=group),hjust = -0.6, alpha=1)+
  ggtree::geom_tippoint(size=6, aes(color=group), alpha=1)+
  scale_colour_manual(values = clst_col,breaks = names(clst_col)) +
  geom_text(aes(x=branch, label=plyr::mapvalues(round(branch.length, 0),from = "0",to = ""), vjust=-.5), size = 3) +
  theme(legend.position = "none") + ggtree::geom_rootpoint() 

cowplot::ggsave2(paste0("./figures/", pro_name_d, "_medicc2_ME_tree_merged_subclones_lumhr.pdf"), treeplt, width = 10, height = 6)


####----RNA data integration--------#######
pro_name_r <- paste0(pro_name, "_rna")
rna.cluster.ids<-c("Fibro","LumHR","LumSec","MyoEpi","ECs","Tumor","Tumor2","Tumor3","RNA_miss","Unknown1","Unknown2","Unknown3","Mito", "CAFs")
rna.cluster.col <- c("#D89000","#FF62BC","#00BD5F","#9590FF","#F8766D","#54278f","#BF4000","#3969AC","grey90","grey60","grey70","grey80","#CAB2D6","#fca311")
names(rna.cluster.col) <- rna.cluster.ids

obj_rna <- readRDS(paste0("/volumes/USR1/ruiye/10x/2020ANALYSIS/WDR/WDR_Paper/DCIS_RNA_Object/SCTV2_no_integration/", toupper(pro_name), "_RNA.sct.rds"))
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))
rna_celltype <- c("Tumor","CAFs","ECs")


# typical_marker<-c("COL1A1","COL1A2","APMAP", "PTPRC", "MSR1", "KRT14","KRT15","ESR1", "KRT19", "MMRN1", "VWF",
#                   "SYNPO2", "ACTA2","CD247","ABCC11", "RP11-206M11.7", "LINC01088","HPGD","IL18R1","BLK",
#                   "BANK1", "BCL11A","ERBB2","GRB7","SERHL2","CENPE", "CDC6", "CIT", "BRCA1","CCNE2","RGS5","RGS6")
# 
# n_col <- 6
# pdf(paste0("./figures/", pro_name, "typical_marker_feature_plot_violin_2.pdf"), width = 3*n_col, height = 3*ceiling(length(typical_marker)/n_col))
# VlnPlot(object = obj_rna, feature = typical_marker, ncol = n_col)
# dev.off()

name_dna <- varbin_mtx_all_log2@colData %>% as.data.frame() %>% dplyr::select(c("sample", "experiment")) %>% 
  mutate(newname = toupper(stringr::str_extract(sample, "c\\d+"))) %>% left_join(wafer_match_list, by = c("newname" = "Cell")) %>% 
  mutate(rna_bc2 = paste(toupper(experiment), RNA_Barcode, sep = "_"))

colnames(obj_rna@meta.data) <- ifelse(stringr::str_detect(colnames(obj_rna@meta.data), "RNA"), colnames(obj_rna@meta.data), 
                                      paste0("RNA_",colnames(obj_rna@meta.data)))

DNA_rna_meta <- obj_rna@meta.data %>% rownames_to_column() %>% mutate(rowname = toupper(rowname)) %>% 
  left_join(name_dna, ., by = c("rna_bc2" = "rowname")) %>% 
  dplyr::select(c("nCount_RNA", "nFeature_RNA", "RNA_dispense", "RNA_S5XX", "RNA_percent.mt", "RNA_nCount_SCT", "RNA_nFeature_SCT", 
                  "RNA_seurat_clusters", "RNA_match")) %>% 
  mutate(cell_type = plyr::mapvalues(RNA_seurat_clusters, from = levels(RNA_seurat_clusters), to = rna_celltype))

varbin_mtx_all_log2@colData <- cbind(varbin_mtx_all_log2@colData, DNA_rna_meta)
saveRDS(varbin_mtx_all_log2, file = paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit_RNA_meta.rds")))
# varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit_RNA_meta.rds")))


#------heatmap with clonal status bar----#####
mtx_srt <- as.data.frame(varbin_mtx_all_log2@colData) %>%
  arrange(factor(subclones3, levels = levels(varbin_mtx_all_log2@colData$subclones3)))
ht_mtx <- log2(t(varbin_mtx_all_log2@assays@data$segment_ratios))[mtx_srt$sample,]

#----annotation bar--
anno_mtx <- mtx_srt %>% mutate(rna_clst0 = as.character(cell_type)) %>% mutate(rna_clst = ifelse(is.na(rna_clst0), "RNA_miss",rna_clst0)) %>% 
  dplyr::select(c("subclones3", "rna_clst", "experiment"))
rownames(anno_mtx) <- NULL

clst_col <- new_pal[1:length(unique(anno_mtx$subclones3))]
names(clst_col) <- levels(varbin_mtx_all_log2@colData$subclones3)

rna_col <- rna.cluster.col[names(rna.cluster.col) %in% unique(anno_mtx$rna_clst)]

exp_col <- jcolors::jcolors("pal2")[(6-length(unique(anno_mtx$experiment))):5]
names(exp_col) <- c(pro_name1, pro_name2)

#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue3", "white", "firebrick3"))
ha_row2=rowAnnotation(df = anno_mtx, col = list(subclones3=clst_col, experiment=exp_col, rna_clst=rna_col),
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/", pro_name_d, "_complexHeatmap_combinded_DR.pdf"), height = 10, width = 10)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        row_split = anno_mtx$subclones3,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "Log2 (Ratio)", 
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row2)
dev.off()

#----------Integer consensus Heatmap with cell types bar----#####
cs_mtx <- t(varbin_mtx_all_log2@consensus)
my_order <- levels(varbin_mtx_all_log2@colData$subclones3)
cs_mtx_order <- cs_mtx[my_order,]
topleft(cs_mtx_order)

clst_col_cs <- new_pal[1:length(unique(varbin_mtx_all_log2@colData$subclones3))]
subclone_cs <- levels(varbin_mtx_all_log2@colData$subclones3)
names(clst_col_cs) <- subclone_cs
ha_row_cs=rowAnnotation(subclones3 = my_order, col = list(subclones3 = clst_col_cs), show_annotation_name = F, show_legend = F)

# rna_col <- rna.cluster.col[names(rna.cluster.col) %in% unique(anno_mtx$rna_clst)]
meta_con <- varbin_mtx_all_log2@colData %>% as.data.frame() %>% dplyr::select(subclones3, cell_type) %>% 
  dplyr::group_by(subclones3, cell_type) %>% dplyr::summarise(n=n()) %>% spread(cell_type, n) %>% dplyr::select(-`<NA>`)
meta_con[is.na(meta_con)] <- 0

meta_con_sum <- meta_con %>% column_to_rownames(var = "subclones3") %>% mutate(mysum=apply(., 1, sum))
meta_con_p <- meta_con_sum/meta_con_sum$mysum 
meta_con_p2 <- meta_con_p %>% dplyr::select(-mysum)

rna_col <- rna.cluster.col[names(rna.cluster.col) %in% rna_celltype]
rna_col_bar <- rna_col[colnames(meta_con_p2)]
ha_barplot <- rowAnnotation(bar= anno_barplot(meta_con_p2[my_order,],
                                              gp = gpar(fill = rna_col_bar, 
                                                        col= rna_col_bar)),
                            show_annotation_name = T)

#-----header
b <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/chr_a.txt")
a <- unlist(b, use.names = F)
chr2 <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/chr_chr2.txt")
ha_col=HeatmapAnnotation(foo=anno_text(a, rot = 0, gp = gpar(fontsize =10)), df =chr2, col = list(chr=c("1"="black", "2"="grey")), 
                         show_legend = F, annotation_name_side = "left")

#-----annotated genes---
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
lab45_uniq_ext <- c("CKDN2C","AKT3","CKS1B","RUVBL1","C8orf4","FGFR1","BAG4","MTDH","MYC","PAK1","CDK4","MDM2", "FHIT","NEDD9","FOXK2","CCNE2",
                    "PLA2G10","GRB7","RPS6KB1","PPM1D","CCNE1","YWHAB","ZNF217","AURKA","PTK6","PPP2R2A",
                    "CCND1","NCOA3","TP53","CDKN1A","GADD45A","MDM2","RB1","CDKN2A","CDKN2C","ESR1","SHC1",
                    "PGR","PIK3CA","PTEN","BCL2","GATA3","FGFR4","EGFR","ERBB2","EMSY","BRCA1","FOXA1","BRCA2")
#---from: Chromosome 16 tumor-suppressor genes in breast cancer --- https://doi.org/10.1002/gcc.20318
chr16_tumor_suppressor <- c("E2F4","CTCF","TERF2","FBXL8","LRRC29","CDH5","CDH16","CDH3","CDH11","ATBF1","BBC1","FANCA","GAS11","CPNE7","ZFP276",
                            "CBFA2T3","CYBA","WWOX","TERF21P","SIAH1") 
chr16_tumor_suppressor2 <- c("E2F4","CTCF","TERF2","FBXL8","CDH5","CDH16","CDH3","CDH11","GAS11","CPNE7","ZFP276","CBFA2T3","WWOX","TERF21P") 

gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% c(lab45_uniq_ext, chr16_tumor_suppressor2))
ha_bottom_cs = columnAnnotation(clonal_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
                                                         side = "bottom", labels_gp = gpar(fontsize = 6)))

heat_col <- pals::ocean.balance(14)[9:14]
if(max(cs_mtx)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(cs_mtx)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(cs_mtx)-6))
}
names(col_vec) <- 0:max(cs_mtx)

pdf(paste0("./figures/", pro_name_d, "_merged_consensus_integer_complexHeatmap_combined_DR.pdf"), height = 4, width = 10)
Heatmap(as.matrix(cs_mtx_order), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        name = "csheatmap", show_row_names = T, row_names_side = "left",
        show_column_names = F, column_title = paste0(pro_name_d, "_consensusHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row_cs, bottom_annotation = ha_bottom_cs, right_annotation = ha_barplot)
dev.off()


meta_con2 <- meta_con
meta_con2[meta_con2 <=2] <- 0
meta_con_sum2 <- meta_con2 %>% column_to_rownames(var = "subclones3") %>% mutate(mysum=apply(., 1, sum))
meta_con_p0 <- meta_con_sum2/meta_con_sum2$mysum 
meta_con_p0 <- meta_con_p0 %>% dplyr::select(-mysum)
meta_con_p0[meta_con_p0 == "NaN"] <- 0
rna_col_bar <- rna_col[colnames(meta_con_p0)]

ha_barplot2 <- rowAnnotation(bar= anno_barplot(meta_con_p0[my_order,],
                                               gp = gpar(fill = rna_col_bar, 
                                                         col= "NA")),
                             show_annotation_name = F)

pdf(paste0("./figures/", pro_name_d, "_merged_consensus_integer_complexHeatmap_combined_DR_larger2cell.pdf"), height = 4, width = 10)
Heatmap(as.matrix(cs_mtx_order), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        name = "csheatmap", show_row_names = T, row_names_side = "left",
        show_column_names = F, column_title = paste0(pro_name_d, "_consensusHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row_cs, bottom_annotation = ha_bottom_cs, right_annotation = ha_barplot2)
dev.off()


###---------Integer consensus Heatmap with cell types bar reordered by tree----#####
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit_RNA_meta.rds")))
pro_name_r <- paste0(pro_name, "_rna")
rna.cluster.ids<-c("Fibro","LumHR","LumSec","MyoEpi","ECs","Tumor","Tumor2","Tumor3","RNA_miss","Unknown1","Unknown2","Unknown3","Mito", "CAFs")
rna.cluster.col <- c("#D89000","#FF62BC","#00BD5F","#9590FF","#F8766D","#54278f","#BF4000","#3969AC","grey90","grey60","grey70","grey80","#CAB2D6","#fca311")
names(rna.cluster.col) <- rna.cluster.ids
rna_celltype <- c("Tumor","CAFs","ECs")

cs_mtx <- t(varbin_mtx_all_log2@consensus)
lab <-  (treeplt[["data"]] %>% arrange(y))$label
my_order <- c("c1","c2","c3",rev(lab[grepl("c[0-9]+", lab, perl = T)]))

cs_mtx_order <- cs_mtx[my_order,]
topleft(cs_mtx_order)

clst_col_cs <- new_pal[1:length(unique(varbin_mtx_all_log2@colData$subclones3))]
subclone_cs <- levels(varbin_mtx_all_log2@colData$subclones3)
names(clst_col_cs) <- subclone_cs
ha_row_cs=rowAnnotation(subclones3 = my_order, col = list(subclones3 = clst_col_cs), show_annotation_name = F, show_legend = F)

# rna_col <- rna.cluster.col[names(rna.cluster.col) %in% unique(anno_mtx$rna_clst)]
meta_con <- varbin_mtx_all_log2@colData %>% as.data.frame() %>% dplyr::select(subclones3, cell_type) %>% 
  dplyr::group_by(subclones3, cell_type) %>% dplyr::summarise(n=n()) %>% spread(cell_type, n) %>% dplyr::select(-`<NA>`)
meta_con[is.na(meta_con)] <- 0

meta_con_sum <- meta_con %>% column_to_rownames(var = "subclones3") %>% mutate(mysum=apply(., 1, sum))
meta_con_p <- meta_con_sum/meta_con_sum$mysum 
meta_con_p2 <- meta_con_p %>% dplyr::select(-mysum)

rna_col <- rna.cluster.col[names(rna.cluster.col) %in% rna_celltype]
rna_col_bar <- rna_col[colnames(meta_con_p2)]
ha_barplot <- rowAnnotation(bar= anno_barplot(meta_con_p2[my_order,],
                                              gp = gpar(fill = rna_col_bar, 
                                                        col= rna_col_bar)),
                            show_annotation_name = T)

#-----header
b <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/chr_a.txt")
a <- unlist(b, use.names = F)
chr2 <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/chr_chr2.txt")
ha_col=HeatmapAnnotation(foo=anno_text(a, rot = 0, gp = gpar(fontsize =10)), df =chr2, col = list(chr=c("1"="black", "2"="grey")), 
                         show_legend = F, annotation_name_side = "left")

#-----annotated genes---
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
lab45_uniq_ext <- c("CKDN2C","AKT3","CKS1B","RUVBL1","C8orf4","FGFR1","BAG4","MTDH","MYC","PAK1","CDK4","MDM2", "FHIT","NEDD9","FOXK2","CCNE2",
                    "PLA2G10","GRB7","RPS6KB1","PPM1D","CCNE1","YWHAB","ZNF217","AURKA","PTK6","PPP2R2A",
                    "CCND1","NCOA3","TP53","CDKN1A","GADD45A","MDM2","RB1","CDKN2A","CDKN2C","ESR1","SHC1",
                    "PGR","PIK3CA","PTEN","BCL2","GATA3","FGFR4","EGFR","ERBB2","EMSY","BRCA1","FOXA1","BRCA2")
#---from: Chromosome 16 tumor-suppressor genes in breast cancer --- https://doi.org/10.1002/gcc.20318
chr16_tumor_suppressor <- c("E2F4","CTCF","TERF2","FBXL8","LRRC29","CDH5","CDH16","CDH3","CDH11","ATBF1","BBC1","FANCA","GAS11","CPNE7","ZFP276",
                            "CBFA2T3","CYBA","WWOX","TERF21P","SIAH1") 
chr16_tumor_suppressor2 <- c("E2F4","CTCF","TERF2","FBXL8","CDH5","CDH16","CDH3","CDH11","GAS11","CPNE7","ZFP276","CBFA2T3","WWOX","TERF21P") 

gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% c(lab45_uniq_ext, chr16_tumor_suppressor2))
ha_bottom_cs = columnAnnotation(clonal_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
                                                         side = "bottom", labels_gp = gpar(fontsize = 6)))

heat_col <- pals::ocean.balance(14)[9:14]
if(max(cs_mtx)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(cs_mtx)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(cs_mtx)-6))
}
names(col_vec) <- 0:max(cs_mtx)

meta_con2 <- meta_con
meta_con2[meta_con2 <=2] <- 0
meta_con_sum2 <- meta_con2 %>% column_to_rownames(var = "subclones3") %>% mutate(mysum=apply(., 1, sum))
meta_con_p0 <- meta_con_sum2/meta_con_sum2$mysum 
meta_con_p0 <- meta_con_p0 %>% dplyr::select(-mysum)
meta_con_p0[meta_con_p0 == "NaN"] <- 0
rna_col_bar <- rna_col[colnames(meta_con_p0)]

ha_barplot2 <- rowAnnotation(bar= anno_barplot(meta_con_p0[my_order,],
                                               gp = gpar(fill = rna_col_bar, 
                                                         col= "NA")),
                             show_annotation_name = F)

pdf(paste0("./figures/", pro_name_d, "_merged_consensus_integer_complexHeatmap_combined_DR_larger2cells_tree_ordered.pdf"), height = 4, width = 8)
Heatmap(as.matrix(cs_mtx_order), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        name = "csheatmap", show_row_names = T, row_names_side = "left",
        show_column_names = F, column_title = paste0(pro_name_d, "_consensusHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row_cs, bottom_annotation = ha_bottom_cs, right_annotation = ha_barplot2)
dev.off()

###-----RNA data UMAP----########
obj_rna <- readRDS(paste0("/volumes/USR1/ruiye/10x/2020ANALYSIS/WDR/WDR_Paper/DCIS_RNA_Object/SCTV2_no_integration/", toupper(pro_name), "_RNA.sct.rds"))
obj_rna_newid<-obj_rna
rna_celltype <- c("Tumor","CAFs","ECs")
names(rna_celltype) <- levels(obj_rna_newid)
obj_rna_newid <- RenameIdents(obj_rna_newid, rna_celltype)

rna_col <- rna.cluster.col[names(rna.cluster.col) %in% rna_celltype]
p1 <- DimPlot(obj_rna_newid, reduction = "umap", label = T, pt.size=0.5) + scale_color_manual(values = rna_col) + theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), 
        legend.position = "none", axis.ticks = element_blank(),axis.text = element_blank(), axis.title = element_text(size = 15))#+
p1
# p2 <- LabelClusters(p1, id = "ident", color = unique(ggplot_build(p1)$data[[1]]$colour), size = 5, repel = T,  box.padding = 0, segment.size=0, fontface="bold")
cowplot::ggsave2(paste0("./figures/", pro_name_r, "_RNA_UMAP.pdf"), p1, width = 4, height = 4) 

##---Top 10 markers---#######
obj_rna_newid.markers <- FindAllMarkers(object = obj_rna_newid, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
obj_rna_newid.markers %>% dplyr::group_by(cluster) %>% top_n(2, avg_log2FC)
top10 <- obj_rna_newid.markers %>% dplyr::group_by(cluster) %>% top_n(10, avg_log2FC)

rna_col2 <- rna_col[rna_celltype]
pdf(paste0(c("./figures/"),pro_name_r, c("_Seurat_top10_ComplexHeatmap.pdf")), height = 10, width = 5.5)
DoComHeatmap2(obj_rna_newid, genelist = top10$gene, col_palette=rna_col2, col_title=T, breaks=c(-2,0,2))
dev.off()


#----Plot progenitor cells and tumor cells with module score---########
##---HBCA top 20 markers feature plot---##
#---markers from HBCA Nature paper Supp_Table3--#
obj_rna3 <- obj_rna_newid
typical_basal <- hbca_epi_markers$Basal
typical_lumhr <- hbca_epi_markers$LumHR
typical_lumsec <- hbca_epi_markers$LumSec
typical_marker <- list(typical_basal, typical_lumhr, typical_lumsec)
names(typical_marker) <- c("Basal_ModScore","LumHR_ModScore","LumSec_ModScore")

obj_rna3 <- AddModuleScore(obj_rna3, features = typical_marker, name = "Epi_modscore")
colnames(obj_rna3@meta.data)[grep("Epi_modscore", colnames(obj_rna3@meta.data))] <- names(typical_marker)

#---tumor cells
df_hall <- obj_rna3@meta.data %>% dplyr::filter(subclone3 %in% paste0("c", 4:14)) %>% dplyr::select(c("subclone3",ends_with("ModScore"))) %>% 
  as.data.frame() %>% rownames_to_column("id")
df_hall_long <- df_hall %>% pivot_longer(cols = -c(id, subclone3), names_to = "ModScore", values_to = "value")

p1 <- ggplot(df_hall_long, aes(x=ModScore, y = value, color=ModScore)) + geom_boxplot() + geom_jitter(position = position_jitter(0.1)) + theme_classic() + 
  theme(axis.title.x = element_blank()) + scale_color_manual(values= c("#9590FF","#FF62BC","#00BD5F")) 
p1
cowplot::ggsave2(paste0("./figures/",pro_name_r, "_tumor_epi_modual_score.pdf"), p1, width = 4.5, height = 3)


#---map DNA clusters to RNA UMAP----#####
plist <- list()
j <- 1
my_clones <- levels(obj_rna_newid@meta.data$subclone3)
for (i in my_clones[1:(length(my_clones)-1)]) {
  obj_rna_newid_tmp <- obj_rna_newid
  obj_rna_newid_tmp@meta.data <- obj_rna_newid_tmp@meta.data %>% mutate(subclone4 = ifelse(subclone3 == i, as.character(subclone3), "unsel"))
  
  mut_cell_df <- obj_rna_newid_tmp@meta.data[,c("seurat_clusters","subclone3")] %>% table() %>% as.data.frame() %>% dplyr::filter(Freq >0 & Freq <3)
  mut_cell_chr <-  mut_cell_df %>% dplyr::filter(subclone3 == i) %>% pull(seurat_clusters) %>% as.character()
  obj_rna_newid_tmp@meta.data <- obj_rna_newid_tmp@meta.data %>% mutate(subclone4 = ifelse((subclone4 == i & seurat_clusters %in% mut_cell_chr), "unsel", subclone4))
  
  p1 <- DimPlot(obj_rna_newid_tmp, reduction = "umap", group.by = "subclone4", pt.size = 1.5, order = i) + 
    scale_color_manual(values = c("grey90",new_pal[j])) + ggtitle("")+ theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), 
          axis.ticks = element_blank(),axis.text = element_blank(), axis.title = element_text(size = 15))
  plist[[j]] <- p1
  j <- j +1
}

p_list1 <- plot_grid(plotlist = plist, ncol = 5)
cowplot::ggsave2(paste0("./figures/", pro_name_r, "_DNA_subclone_on_RNA_UMAP.pdf"), p_list1, width = 18, height = 3*(length(my_clones)-1)/5, limitsize = F) 








#----------Integer consensus Heatmap with cell types bar----#####
cs_mtx <- t(varbin_mtx_all_log2@consensus)
my_order <- levels(varbin_mtx_all_log2@colData$subclones3)
cs_mtx_order <- cs_mtx[my_order,]
topleft(cs_mtx_order)

clst_col_cs <- new_pal[1:length(unique(varbin_mtx_all_log2@colData$subclones3))]
subclone_cs <- levels(varbin_mtx_all_log2@colData$subclones3)
names(clst_col_cs) <- subclone_cs
ha_row_cs=rowAnnotation(subclones3 = my_order, col = list(subclones3 = clst_col_cs), show_annotation_name = F, show_legend = F)

# rna_col <- rna.cluster.col[names(rna.cluster.col) %in% unique(anno_mtx$rna_clst)]
meta_con <- varbin_mtx_all_log2@colData %>% as.data.frame() %>% dplyr::select(subclones3, cell_type) %>% 
  dplyr::group_by(subclones3, cell_type) %>% dplyr::summarise(n=n()) %>% spread(cell_type, n) %>% dplyr::select(-`<NA>`)
meta_con[is.na(meta_con)] <- 0

meta_con_sum <- meta_con %>% column_to_rownames(var = "subclones3") %>% mutate(mysum=apply(., 1, sum))
meta_con_p <- meta_con_sum/meta_con_sum$mysum 
meta_con_p2 <- meta_con_p %>% dplyr::select(-mysum)

rna_col <- rna.cluster.col[names(rna.cluster.col) %in% rna_celltype]
rna_col_bar <- rna_col[colnames(meta_con_p2)]
ha_barplot <- rowAnnotation(bar= anno_barplot(meta_con_p2[my_order,],
                                              gp = gpar(fill = rna_col_bar, 
                                                        col= rna_col_bar)),
                            show_annotation_name = T)

#-----header
b <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/chr_a.txt")
a <- unlist(b, use.names = F)
chr2 <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/chr_chr2.txt")
ha_col=HeatmapAnnotation(foo=anno_text(a, rot = 0, gp = gpar(fontsize =10)), df =chr2, col = list(chr=c("1"="black", "2"="grey")), 
                         show_legend = F, annotation_name_side = "left")

#-----annotated genes---
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
lab45_uniq_ext <- c("CKDN2C","AKT3","CKS1B","RUVBL1","C8orf4","FGFR1","BAG4","MTDH","MYC","PAK1","CDK4","MDM2", "FHIT","NEDD9","FOXK2","CCNE2",
                    "PLA2G10","GRB7","RPS6KB1","PPM1D","CCNE1","YWHAB","ZNF217","AURKA","PTK6","PPP2R2A",
                    "CCND1","NCOA3","TP53","CDKN1A","GADD45A","MDM2","RB1","CDKN2A","CDKN2C","ESR1","SHC1",
                    "PGR","PIK3CA","PTEN","BCL2","GATA3","FGFR4","EGFR","ERBB2","EMSY","BRCA1","FOXA1","BRCA2")
#---from: Chromosome 16 tumor-suppressor genes in breast cancer --- https://doi.org/10.1002/gcc.20318
chr16_tumor_suppressor <- c("E2F4","CTCF","TERF2","FBXL8","LRRC29","CDH5","CDH16","CDH3","CDH11","ATBF1","BBC1","FANCA","GAS11","CPNE7","ZFP276",
                            "CBFA2T3","CYBA","WWOX","TERF21P","SIAH1") 
chr16_tumor_suppressor2 <- c("E2F4","CTCF","TERF2","FBXL8","CDH5","CDH16","CDH3","CDH11","GAS11","CPNE7","ZFP276","CBFA2T3","WWOX","TERF21P") 

gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% c(lab45_uniq_ext, chr16_tumor_suppressor2))
ha_bottom_cs = columnAnnotation(clonal_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
                                                         side = "bottom", labels_gp = gpar(fontsize = 6)))

heat_col <- pals::ocean.balance(14)[9:14]
if(max(cs_mtx)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(cs_mtx)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(cs_mtx)-6))
}
names(col_vec) <- 0:max(cs_mtx)

pdf(paste0("./figures/", pro_name_d, "_merged_consensus_integer_complexHeatmap_combined_DR.pdf"), height = 4, width = 10)
Heatmap(as.matrix(cs_mtx_order), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        name = "csheatmap", show_row_names = T, row_names_side = "left",
        show_column_names = F, column_title = paste0(pro_name_d, "_consensusHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row_cs, bottom_annotation = ha_bottom_cs, right_annotation = ha_barplot)
dev.off()




#---calculate basic QC matrics----####
obj_rna <- readRDS(paste0("/volumes/USR1/ruiye/10x/2020ANALYSIS/WDR/WDR_Paper/DCIS_RNA_Object/SCTV2_no_integration/", toupper(pro_name), "_RNA.sct.rds"))
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))

name_dna <- varbin_mtx_all_log2@colData %>% as.data.frame() %>% dplyr::select(c("sample", "experiment")) %>% 
  mutate(newname = toupper(stringr::str_extract(sample, "c\\d+"))) %>% left_join(wafer_match_list, by = c("newname" = "Cell")) %>% 
  mutate(rna_bc2 = paste(toupper(experiment), RNA_Barcode, sep = "_"))
colnames(obj_rna@meta.data) <- ifelse(stringr::str_detect(colnames(obj_rna@meta.data), "RNA"), colnames(obj_rna@meta.data), 
                                      paste0("RNA_",colnames(obj_rna@meta.data)))
obj_rna@meta.data %>% rownames_to_column() %>% mutate(rowname = toupper(rowname)) %>% full_join(name_dna, ., by = c("rna_bc2" = "rowname")) %>% dim()

varbin_mtx_all_log2@colData %>% dim()
varbin_mtx_all_log2@colData %>% as.data.frame() %>% pull(reads_total) %>% mean()
varbin_mtx_all_log2@colData %>% as.data.frame() %>% pull(reads_assigned_bins) %>% mean()
varbin_mtx_all_log2@colData %>% as.data.frame() %>% pull(percentage_duplicates) %>% mean()
varbin_mtx_all_log2@colData %>% as.data.frame() %>% pull(median_bin_count) %>% mean()
obj_rna@meta.data %>% dim()
obj_rna@meta.data$nCount_RNA %>% mean()
obj_rna@meta.data$nFeature_RNA %>% mean()
#----DNA subclone frequency in all cells-----####
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit_RNA_meta.rds")))

celltype_freq <- table(varbin_mtx_all_log2@colData$cell_type, varbin_mtx_all_log2@colData$subclones3) %>% as.data.frame.matrix()
celltype_freq[celltype_freq <=2] <- 0

sum(sapply(celltype_freq, sum))

3/781; 9/781

epi_freq2 <- epi_freq[1:3, 2:7]/1625

freq <- apply(epi_freq2, 2, sum) %>% as.data.frame() %>% rownames_to_column() %>% mutate(freq = `.`*100) %>% 
  mutate(celltype=ifelse(rowname=="c5", "LumHR", "LumSec")) %>% mutate(sample = "ECIS44")


p1 <-  ggplot(freq, aes(x = rowname, y = freq, fill = celltype)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = c("#FF62BC", "#00BD5F")) + 
  theme_classic() +
  xlab("DNA subclones") +
  ylab("Percentage in all epithelial cells") 

cowplot::ggsave2(paste0("./figures/", pro_name_d, "_aneuploidy_in_all_epi_freq.pdf"), p1, width = 4, height = 3)

#---Statistic counts of DNA RNA matching------#########
obj_rna <- readRDS(paste0("/volumes/USR1/ruiye/10x/2020ANALYSIS/WDR/WDR_Paper/DCIS_RNA_Object/SCTV2_no_integration/", toupper(pro_name), "_RNA.sct.rds"))
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))
obj_rna2 <- obj_rna
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
bin_chrpos <- varbin_mtx_all_log2@rowRanges %>% as.data.frame() %>% rownames_to_column() %>% 
  mutate(chr_pos = paste0("chr",seqnames,":",start, "-", end)) %>% dplyr::select(rowname, chr_pos)

subclone_count <- obj_rna2@meta.data %>% dplyr::filter(subclone3 %in% (paste0("c", 2:30))) %>% pull(subclone3) %>% table() 

#---use minimum 25 cell with RNA data as subclone cut-off---
min_cell_num <- 25
min_cell_num <- 50
subclone_count <- subclone_count[subclone_count > min_cell_num] %>% names()
combn_df <- combn(subclone_count, 2) %>% t() %>% as.data.frame() 
combn_df

dna_seg_rna_df <- NULL
for (i in 1:nrow(combn_df)) {
  # i <- 3
  super_clone1 <- combn_df[i,1]
  super_clone2 <- combn_df[i,2]
  super_name1 <- paste(super_clone1, collapse = "_")
  super_name2 <- paste0(super_clone2, collapse = "_")
  
  obj_rna3 <- obj_rna2
  obj_rna3@meta.data <- obj_rna3@meta.data %>% mutate(temp = ifelse(subclone3 %in% super_clone1, "super_c1",
                                                                    ifelse(subclone3 %in% super_clone2, "super_c2", "others")))
  clone1_cells <- obj_rna3@meta.data %>% dplyr::filter(temp == "super_c1") %>% rownames()
  rna_exp_clone1_cells <- obj_rna3@assays$RNA@data[,clone1_cells]
  clone2_cells <- obj_rna3@meta.data %>% dplyr::filter(temp == "super_c2") %>% rownames()
  rna_exp_clone2_cells <- obj_rna3@assays$RNA@data[,clone2_cells]
  
  rna_cons_clone1 <- apply(rna_exp_clone1_cells, 1, mean) %>% as.data.frame() %>% `colnames<-`("mean_scaled_exp_clone1") %>% 
    rownames_to_column() %>% dplyr::filter(mean_scaled_exp_clone1 > 0)
  rna_cons_clone2 <- apply(rna_exp_clone2_cells, 1, mean) %>% as.data.frame() %>% `colnames<-`("mean_scaled_exp_clone2") %>% 
    rownames_to_column() %>% dplyr::filter(mean_scaled_exp_clone2 > 0)
  
  clone1_2_con_exp <- inner_join(rna_cons_clone1, rna_cons_clone2, by = "rowname")
  dna_con <- varbin_mtx_all_log2@consensus %>% rownames_to_column() %>% dplyr::select(rowname,all_of(c(super_clone1, super_clone2)))
  
  dna_rna_merged_cna_exp <- gene_bin %>% mutate(pos = as.character(pos)) %>% left_join(dna_con, ., by = c("rowname" = "pos")) %>% 
    left_join(clone1_2_con_exp, by = c("gene" = "rowname")) %>% left_join(bin_chrpos, by = "rowname")
  
  dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% 
    mutate(mean_scaled_exp_clone1 = ifelse(mean_scaled_exp_clone1>0.05, mean_scaled_exp_clone1, NA), 
           mean_scaled_exp_clone2 = ifelse(mean_scaled_exp_clone2>0.05, mean_scaled_exp_clone2, NA)) %>% 
    mutate(diff = mean_scaled_exp_clone1 - mean_scaled_exp_clone2) %>% mutate(
      normal_clone1 = case_when(
        !is.na(mean_scaled_exp_clone1) & !is.na(mean_scaled_exp_clone2) ~ mean_scaled_exp_clone1/mean_scaled_exp_clone2,
        !is.na(mean_scaled_exp_clone1) & is.na(mean_scaled_exp_clone2) ~ mean_scaled_exp_clone1/0.05,
        is.na(mean_scaled_exp_clone1) & !is.na(mean_scaled_exp_clone2) ~ 0.05/mean_scaled_exp_clone2,
        TRUE ~ NA_real_  # Returns NA for numeric columns
      ),
      normal_clone2 = ifelse(is.na(normal_clone1), NA, 1), 
      cna_diff = get(super_clone1) - get(super_clone2),
      cna_base = 0,
      chr = stringr::str_extract(chr_pos, "chr[^:]+")
    ) %>% dplyr::select(rowname, chr, diff,cna_diff) %>% na.omit(cols="diff")
  
  
  data <- dna_rna_merged_cna_exp_target %>% dplyr::group_by(chr) %>% 
    mutate(group = cumsum(cna_diff != lag(cna_diff, default = first(cna_diff)))) %>% ungroup()
  
  # Grouping by 'chr' and 'group' and then calculating the mean of 'diff'
  grouped_median <- data %>% dplyr::group_by(chr, group) %>%
    dplyr::summarise(median_diff = median(diff), first_cna_diff = first(cna_diff))
  
  write.csv(grouped_median, paste0("./metrics/genotype_seg_vs_sudo_bulk_expr_seg/",pro_name, "_seg_level_genotype_phenotype_statistic_count",super_clone1, "_vs_",super_clone2,".csv"))
  # grouped_median <- read.csv(paste0("./metrics/genotype_seg_vs_sudo_bulk_expr_seg/",pro_name, "_seg_level_genotype_phenotype_statistic_count",super_clone1, "_vs_",super_clone2,".csv"), row.names = 1)
  
  temp_vec <- paste0(super_clone1, "_", super_clone2)
  good_ones <- grouped_median %>% dplyr::filter(first_cna_diff != 0 & abs(median_diff) >=0.03)
  dna_but_no_rna <- grouped_median %>% dplyr::filter(first_cna_diff != 0 & abs(median_diff) < 0.03)
  no_dna_but_rna <- grouped_median %>% dplyr::filter(first_cna_diff == 0 & abs(median_diff) >= 0.03)
  temp_vec <- append(temp_vec, c(nrow(good_ones), nrow(dna_but_no_rna), nrow(no_dna_but_rna)))
  
  dna_seg_rna_df <- rbind(dna_seg_rna_df, temp_vec)
}


colnames(dna_seg_rna_df) <- c("subclone_comb", "matched", "dna_but_no_rna", "no_dna_but_rna")
rownames(dna_seg_rna_df) <- NULL
dna_seg_rna_df <- dna_seg_rna_df %>% as.data.frame()

write.csv(dna_seg_rna_df, paste0("./metrics/",pro_name, "_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_", min_cell_num, "_Cells.csv"))
# dna_seg_rna_df2 <- read.csv(paste0("./metrics/",pro_name, "_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_", min_cell_num, "_Cells.csv"), row.names = 1)

dna_seg_rna_df <- dna_seg_rna_df %>% mutate(matched = as.integer(matched), dna_but_no_rna = as.integer(dna_but_no_rna), 
                                            no_dna_but_rna = as.integer(no_dna_but_rna))
sum_matched <- sum(dna_seg_rna_df$matched)
sum_dna_but_no_rna <- sum(dna_seg_rna_df$dna_but_no_rna)
sum_no_dna_but_rna <- sum(dna_seg_rna_df$no_dna_but_rna)

sum_matched/(sum_matched + sum_dna_but_no_rna)
# 25 cells: 0.7834101
# 50 cells: 0.8387097

#---venn diagram---
library("VennDiagram")
library("RColorBrewer")

area1_value = sum_dna_but_no_rna + sum_matched
area2_value = sum_no_dna_but_rna + sum_matched
if(area1_value > area2_value) {
  category_names <- c("CNA diff", "RNA diff")
  fill_col <- brewer.pal(3, "Set2")[1:2]
} else {
  category_names <- c("RNA diff", "CNA diff")
  fill_col <- brewer.pal(3, "Set2")[2:1]
}

p1 <- draw.pairwise.venn(area1 = area1_value, area2 = area2_value, cross.area = sum_matched, fill = fill_col, 
                         category = category_names, lty = "blank", ext.line.lty = "dashed")
p2 <- plot_grid(grobTree(p1))
p2
cowplot::ggsave2(paste0("./figures/", pro_name, "_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_",min_cell_num,"_Cells_0.03_venn.pdf"), 
                 p2, width = 3, height = 2.5)


#----cut off test---
min_cell_num <- 25
subclone_count <- subclone_count[subclone_count > min_cell_num] %>% names()
combn_df <- combn(subclone_count, 2) %>% t() %>% as.data.frame() 
combn_df

dna_seg_rna_df <- NULL
for (i in 1:nrow(combn_df)) {
  # i <- 3
  super_clone1 <- combn_df[i,1]
  super_clone2 <- combn_df[i,2]
  super_name1 <- paste(super_clone1, collapse = "_")
  super_name2 <- paste0(super_clone2, collapse = "_")
  
  grouped_median <- read.csv(paste0("./metrics/genotype_seg_vs_sudo_bulk_expr_seg/",pro_name, "_seg_level_genotype_phenotype_statistic_count",super_clone1, "_vs_",super_clone2,".csv"), row.names = 1)
  temp_vec <- paste0(super_clone1, "_", super_clone2)
  for (j in c(0.01, 0.02, 0.03, 0.04,0.05)) {
    good_ones <- grouped_median %>% dplyr::filter(first_cna_diff != 0 & abs(median_diff) >= j)
    dna_but_no_rna <- grouped_median %>% dplyr::filter(first_cna_diff != 0 & abs(median_diff) < j)
    no_dna_but_rna <- grouped_median %>% dplyr::filter(first_cna_diff == 0 & abs(median_diff) >= j)
    temp_vec <- append(temp_vec, c(nrow(good_ones), nrow(dna_but_no_rna), nrow(no_dna_but_rna)))
  }
  dna_seg_rna_df <- rbind(dna_seg_rna_df, temp_vec)
}

colnames(dna_seg_rna_df) <- c("subclone_comb", "matched_0.01", "dna_but_no_rna_0.01", "no_dna_but_rna_0.01",
                              "matched_0.02", "dna_but_no_rna_0.02", "no_dna_but_rna_0.02",
                              "matched_0.03", "dna_but_no_rna_0.03", "no_dna_but_rna_0.03",
                              "matched_0.04", "dna_but_no_rna_0.04", "no_dna_but_rna_0.04",
                              "matched_0.05", "dna_but_no_rna_0.05", "no_dna_but_rna_0.05")
rownames(dna_seg_rna_df) <- NULL
dna_seg_rna_df <- dna_seg_rna_df %>% as.data.frame()

write.csv(dna_seg_rna_df, paste0("./metrics/",pro_name, "_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_", min_cell_num, "_Cells_cutoff_test.csv"))
# dna_seg_rna_df <- read.csv(paste0("./metrics/",pro_name, "_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_", min_cell_num, "_Cells_cutoff_test.csv"), row.names = 1)

data_long <- gather(dna_seg_rna_df, key = "category", value = "count", -subclone_comb)
data_long$panel <- sub("^.*\\_", "", data_long$category)

data_long <- data_long %>% mutate(category2 = stringr::str_sub(category, 1, -6), count = as.integer(count))

p1 <- ggplot(data_long, aes(x = subclone_comb, y = count, fill = category2)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~panel, scales = "free_x") +
  theme_cowplot() +
  labs(x = "Subclone Combination", y = "Count", title = paste0(pro_name, " Stacked Bar Plot for cut-off testing")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

cowplot::ggsave2(paste0("./figures/", pro_name, "_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_",min_cell_num,"_Cells_cutoff_test.pdf"), p1, width = 8, height = 5)


#---genetic distance vs DE gene number---#####
obj_rna <- readRDS(paste0("/volumes/USR1/ruiye/10x/2020ANALYSIS/WDR/WDR_Paper/DCIS_RNA_Object/SCTV2_no_integration/", toupper(pro_name), "_RNA.sct.rds"))
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))
obj_rna2 <- obj_rna

subclone_count <- obj_rna2@meta.data %>% dplyr::filter(subclone3 %in% (paste0("c", 2:30))) %>% pull(subclone3) %>% table() 

#---use minimum 25 cell with RNA data as subclone cut-off---
min_cell_num <- 25
subclone_count <- subclone_count[subclone_count > min_cell_num] %>% names()
combn_df <- combn(subclone_count, 2) %>% t() %>% as.data.frame() 
combn_df
de.markers_df <- NULL
for (i in 1:nrow(combn_df)) {
  # i <- 3
  super_clone1 <- combn_df[i,1]
  super_clone2 <- combn_df[i,2]
  super_name1 <- paste(super_clone1, collapse = "_")
  super_name2 <- paste0(super_clone2, collapse = "_")
  
  obj_rna3 <- obj_rna2
  obj_rna3@meta.data <- obj_rna3@meta.data %>% mutate(temp = ifelse(subclone3 %in% super_clone1, "super_c1",
                                                                    ifelse(subclone3 %in% super_clone2, "super_c2", "others")))
  
  de.markers <- FindMarkers(obj_rna3, ident.1 = "super_c1", ident.2 = "super_c2", group.by = 'temp', assay = "RNA", only.pos = F)
  write.csv(de.markers, paste0("./metrics/de_marker_genetic_dis/",pro_name_r, "_default_de_markers_",super_clone1, "_vs_",super_clone2,".csv"))
  # de.markers <- read.csv(paste0("./metrics/de_marker_genetic_dis/",pro_name_r, "_default_de_markers_",super_clone1, "_vs_",super_clone2,".csv"), row.names = 1)
  #--0.263: 1.2X. 0.485: 1.4x. 0.678: 1.6x; 0.848: 1.8x. 1: 2X
  temp_vec <- paste0(super_clone1, "_", super_clone2)
  for (j in c(0.263, 0.485, 0.678, 0.848, 1)) {
    # j <- 0.485
    de_genes_up_down <- de.markers %>% dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC) > j)
    temp_vec <- append(temp_vec, nrow(de_genes_up_down))
  }
  de.markers_df <- rbind(de.markers_df, temp_vec)
}

colnames(de.markers_df) <- c("subclone_comb", "log2fc_0.263", "log2fc_0.485", "log2fc_0.678", "log2fc_0.848", "log2fc_1")
rownames(de.markers_df) <- NULL
de.markers_df <- de.markers_df %>% as.data.frame()

#----DNA subclone genetic distance ---
medic_dis <- read_tsv(paste0("./metrics/medicc_files_lumhr/", pro_name_d, "/", pro_name_d, "_lumhr_medicc2_input_pairwise_distances.tsv"))
medic_dis2 <- medic_dis %>% column_to_rownames(var="sample_id")
medic_dis2 <- medic_dis2[2:(ncol(medic_dis2)-1), 2:(ncol(medic_dis2)-1)] %>% as.matrix()

melted_data <- reshape2::melt(medic_dis2) %>% mutate(subclone_comb = paste0(Var1, "_",Var2)) %>% dplyr::select(value, subclone_comb)
colnames(melted_data) <- c("cna_dist", "subclone_comb")

genetic_phyno_mtx <- left_join(de.markers_df, melted_data, by = c("subclone_comb" = "subclone_comb"))

write.csv(genetic_phyno_mtx, paste0("./metrics/",pro_name, "_genetic_distance_phenotype_diff_subclonesWithMorethan_", min_cell_num, "_Cells.csv"))
# genetic_phyno_mtx <- read.csv(paste0("./metrics/",pro_name_r, "_genetic_distance_phenotype_diff_subclonesWithMorethan50Cells.csv"), row.names = 1)

genetic_phyno_mtx2 <- genetic_phyno_mtx %>% as.data.frame() %>% mutate(log2fc_0.485 = as.integer(log2fc_0.485))
p1 <- ggplot(genetic_phyno_mtx2, aes(x = cna_dist, y = log2fc_0.485)) +
  geom_point() + geom_line() + 
  labs(x = "cna_dist", y = "# of DE genes", title = paste0(pro_name, "_min",min_cell_num,"cells")) +
  theme_cowplot()
p1

cowplot::ggsave2(paste0("./figures/", pro_name, "_genetic_distance_phenotype_diff_",min_cell_num,"_Cells_cutoff.pdf"), p1, width = 3, height = 2.8)

# genetic_phyno_mtx_long <- melt(genetic_phyno_mtx, id.vars = c("subclone_comb", "cna_dist")) %>% as.data.frame() %>% mutate(value = as.integer(value))
# 
# # Create the scatter plot
# ggplot(genetic_phyno_mtx_long, aes(x = cna_dist, y = value, color = variable)) +
#   geom_point() +
#   labs(x = "cna_dist", y = "Value") +
#   theme_minimal()






#------DE gene ratio in clonal or subclonal region of all subclones-----#####
obj_rna <- readRDS(paste0("/volumes/USR1/ruiye/10x/2020ANALYSIS/WDR/WDR_Paper/DCIS_RNA_Object/SCTV2_no_integration/", toupper(pro_name), "_RNA.sct.rds"))
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
obj_rna2 <- obj_rna

subclone_count <- obj_rna2@meta.data %>% dplyr::filter(subclone3 %in% (paste0("c", 2:30))) %>% pull(subclone3) %>% table() 
#---use minimum 25 cell with RNA data as subclone cut-off---
min_cell_num <- 25
subclone_count <- subclone_count[subclone_count > min_cell_num] %>% names()
combn_df <- combn(subclone_count, 2) %>% t() %>% as.data.frame() 
combn_df
de.markers_df <- NULL
temp_vec <- NULL
for (i in 1:nrow(combn_df)) {
  # i <- 3
  super_clone1 <- combn_df[i,1]
  super_clone2 <- combn_df[i,2]
  super_name1 <- paste(super_clone1, collapse = "_")
  super_name2 <- paste0(super_clone2, collapse = "_")
  
  dna_con <- varbin_mtx_all_log2@consensus %>% rownames_to_column() %>% dplyr::select(rowname,all_of(c(super_clone1, super_clone2))) %>% 
    mutate(rowname = as.integer(rowname))
  
  # obj_rna3 <- obj_rna2
  # obj_rna3@meta.data <- obj_rna3@meta.data %>% mutate(temp = ifelse(subclone3 %in% super_clone1, "super_c1",
  #                                                                   ifelse(subclone3 %in% super_clone2, "super_c2", "others")))
  # 
  # de.markers <- FindMarkers(obj_rna3, ident.1 = "super_c1", ident.2 = "super_c2", group.by = 'temp', assay = "RNA", only.pos = F)
  # write.csv(de.markers, paste0("./metrics/de_marker_genetic_dis/",pro_name_r, "_default_de_markers_",super_clone1, "_vs_",super_clone2,".csv"))
  de.markers <- read.csv(paste0("./metrics/de_marker_genetic_dis/",pro_name_r, "_default_de_markers_",super_clone1, "_vs_",super_clone2,".csv"), row.names = 1)
  #--0.263: 1.2X. 0.485: 1.4x. 0.678: 1.6x; 0.848: 1.8x. 1: 2X
  temp_vec <- paste0(super_clone1, "_", super_clone2)
  for (j in c(0.263, 0.485, 0.678, 0.848, 1)) {
    # j <- 0.485
    de_genes_up_down <- de.markers %>% dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC) > j)
    clonal_reg <- de_genes_up_down %>% rownames_to_column() %>% left_join(gene_bin, by = c("rowname" = "gene")) %>% 
      left_join(dna_con, by = c("pos" = "rowname")) %>% 
      dplyr::filter(pos > 0) %>% dplyr::filter(get(super_clone1) == get(super_clone2))
    
    subclonal_reg <- de_genes_up_down %>% rownames_to_column() %>% left_join(gene_bin, by = c("rowname" = "gene")) %>% 
      left_join(dna_con, by = c("pos" = "rowname")) %>% 
      dplyr::filter(pos > 0) %>% dplyr::filter(get(super_clone1) != get(super_clone2))
    temp_vec <- append(temp_vec, c(nrow(clonal_reg), nrow(subclonal_reg)))
  }
  de.markers_df <- rbind(de.markers_df, temp_vec)
}

colnames(de.markers_df) <- c("subclone_comb", "clonal_log2fc_0.263", "subclonal_log2fc_0.263", "clonal_log2fc_0.485","subclonal_log2fc_0.485",
                             "clonal_log2fc_0.678", "subclonal_log2fc_0.678", "clonal_log2fc_0.848", "subclonal_log2fc_0.848", 
                             "clonal_log2fc_1", "subclonal_log2fc_1")
rownames(de.markers_df) <- NULL
de.markers_df <- de.markers_df %>% as.data.frame()

write.csv(de.markers_df, paste0("./metrics/",pro_name, "_DE_gene_in_clonal_subclonal_regions_WithMorethan_", min_cell_num, "_Cells_cutoff_test.csv"))
# de.markers_df <- read.csv(paste0("./metrics/",pro_name, "_DE_gene_in_clonal_subclonal_regions_WithMorethan_", min_cell_num, "_Cells_cutoff_test.csv"),
# row.names = 1)

data_long <- gather(de.markers_df, key = "category", value = "count", -subclone_comb)
data_long2 <- data_long %>% separate(category, into = c("CNA_status","log2fc"), sep = "_", extra = "merge") %>% mutate(count = as.integer(count))

p1 <- ggplot(data_long2, aes(x = subclone_comb, y = count, fill = CNA_status)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~log2fc, scales = "free_x") +
  theme_cowplot() +
  labs(x = "Subclone Combination", y = "Count", title = paste0(pro_name, " Stacked Bar Plot for cut-off testing")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1

cowplot::ggsave2(paste0("./figures/", pro_name, "_DE_gene_in_clonal_subclonal_regions_WithMorethan_",min_cell_num,"_Cells_cutoff_test.pdf"), 
                 p1, width = 8, height = 5)

#---only plot log2fc > 0.485, and with DE gene > 0---
data_long485 <- data_long2 %>% filter(log2fc == "log2fc_0.485" & count >0)

data_long485_sum <- data_long485 %>% dplyr::group_by(CNA_status) %>% dplyr::summarise(sum(count)) 
data_long485_pct <- data_long485_sum[1,2]/(data_long485_sum[1,2] + data_long485_sum[2,2]) 

p1 <- ggplot(data_long485, aes(x = subclone_comb, y = count, fill = CNA_status)) +
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("clonal" = "#424451", "subclonal" = "#F47F20")) + theme_cowplot() +
  labs(x = "Subclone Combination", y = "# of DE genes", 
       title = pro_name,
       subtitle = paste0("clonal_reg_ratio: ",data_long485_sum[1,2], "/(",data_long485_sum[1,2], "+", data_long485_sum[2,2], ") = ", data_long485_pct)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p1

cowplot::ggsave2(paste0("./figures/", pro_name, "_DE_gene_in_clonal_subclonal_regions_WithMorethan_",min_cell_num,"_Cells_log2fc_0.485.pdf"), 
                 p1, width = 4, height = 3)

#---genetic distance vs DE gene number--randomly down-sampled---#####
obj_rna <- readRDS(paste0("/volumes/USR1/ruiye/10x/2020ANALYSIS/WDR/WDR_Paper/DCIS_RNA_Object/SCTV2_no_integration/", toupper(pro_name), "_RNA.sct.rds"))
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))
obj_rna2 <- obj_rna

subclone_count <- obj_rna2@meta.data %>% dplyr::filter(subclone3 %in% (paste0("c", 2:30))) %>% pull(subclone3) %>% table() 

#---use minimum 25 cell with RNA data as subclone cut-off---
min_cell_num <- 25
subsample_min <- min(subclone_count[subclone_count > min_cell_num])
subclone_count <- subclone_count[subclone_count > min_cell_num] %>% names()
combn_df <- combn(subclone_count, 2) %>% t() %>% as.data.frame() 
combn_df

de.markers_df <- NULL
for (i in 1:nrow(combn_df)) {
  # i <- 1
  super_clone1 <- combn_df[i,1]
  super_clone2 <- combn_df[i,2]
  super_name1 <- paste(super_clone1, collapse = "_")
  super_name2 <- paste0(super_clone2, collapse = "_")
  
  obj_rna3 <- obj_rna2
  obj_rna3@meta.data <- obj_rna3@meta.data %>% mutate(temp = ifelse(subclone3 %in% super_clone1, "super_c1",
                                                                    ifelse(subclone3 %in% super_clone2, "super_c2", "others")))
  #---subset cells as the min cells of those subclones===
  set.seed(100)
  super1_sub <- obj_rna3@meta.data %>% rownames_to_column() %>% dplyr::filter(temp == "super_c1") %>%  pull(rowname) %>% sample(size = subsample_min)
  super2_sub <- obj_rna3@meta.data %>% rownames_to_column() %>% dplyr::filter(temp == "super_c2") %>%  pull(rowname) %>% sample(size = subsample_min)
  obj_rna4 <- subset(obj_rna3, cells = c(super1_sub,super2_sub))
  de.markers <- FindMarkers(obj_rna4, ident.1 = "super_c1", ident.2 = "super_c2", group.by = 'temp', assay = "RNA", only.pos = F)
  
  write.csv(de.markers, paste0("./metrics/de_marker_genetic_dis/",pro_name_r, "_default_de_markers_",super_clone1, "_vs_",super_clone2,"_downsampled_",subsample_min,".csv"))
  # de.markers <- read.csv(paste0("./metrics/de_marker_genetic_dis/",pro_name_r, "_default_de_markers_",super_clone1, "_vs_",super_clone2,"_downsampled_",subsample_min,".csv"), row.names = 1)
  #--0.263: 1.2X. 0.485: 1.4x. 0.678: 1.6x; 0.848: 1.8x. 1: 2X
  temp_vec <- paste0(super_clone1, "_", super_clone2)
  for (j in c(0.263, 0.485, 0.678, 0.848, 1)) {
    # j <- 0.485
    de_genes_up_down <- de.markers %>% dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC) > j)
    temp_vec <- append(temp_vec, nrow(de_genes_up_down))
  }
  de.markers_df <- rbind(de.markers_df, temp_vec)
}

colnames(de.markers_df) <- c("subclone_comb", "log2fc_0.263", "log2fc_0.485", "log2fc_0.678", "log2fc_0.848", "log2fc_1")
rownames(de.markers_df) <- NULL
de.markers_df <- de.markers_df %>% as.data.frame()

#----DNA subclone genetic distance ---
medic_dis <- read_tsv(paste0("./metrics/medicc_files_lumhr/", pro_name_d, "/", pro_name_d, "_lumhr_medicc2_input_pairwise_distances.tsv"))
medic_dis2 <- medic_dis %>% column_to_rownames(var="sample_id")
medic_dis2 <- medic_dis2[2:(ncol(medic_dis2)-1), 2:(ncol(medic_dis2)-1)] %>% as.matrix()

melted_data <- reshape2::melt(medic_dis2) %>% mutate(subclone_comb = paste0(Var1, "_",Var2)) %>% dplyr::select(value, subclone_comb)
colnames(melted_data) <- c("cna_dist", "subclone_comb")

genetic_phyno_mtx <- left_join(de.markers_df, melted_data, by = c("subclone_comb" = "subclone_comb"))

write.csv(genetic_phyno_mtx, paste0("./metrics/",pro_name, "_genetic_distance_phenotype_diff_subclonesWithMorethan_", min_cell_num, "_Cells_downsampled_",subsample_min,".csv"))
# genetic_phyno_mtx <- read.csv(paste0("./metrics/",pro_name_r, "_genetic_distance_phenotype_diff_subclonesWithMorethan50Cells.csv"), row.names = 1)

genetic_phyno_mtx2 <- genetic_phyno_mtx %>% as.data.frame() %>% mutate(log2fc_0.485 = as.integer(log2fc_0.485))
p1 <- ggplot(genetic_phyno_mtx2, aes(x = cna_dist, y = log2fc_0.485)) +
  geom_point() + geom_line() + geom_smooth(method = lm) + stat_cor(method = "spearman") + 
  labs(x = "cna_dist", y = "# of DE genes", title = paste0(pro_name, "_min",min_cell_num,"cells_downsampled_",subsample_min)) +
  theme_cowplot()
p1

cowplot::ggsave2(paste0("./figures/", pro_name, "_genetic_distance_phenotype_diff_",min_cell_num,"_Cells_cutoff_downsampled_",subsample_min,".pdf"), p1, width = 3, height = 2.8)







