#---preload----#####
unloadNamespace("copykit")
library(uwot, lib.loc = "/opt/R/4.1.2/lib/R/library")
packageVersion("uwot")
#--uwot version:  0.1.14
library("Seurat")
packageVersion("Seurat")
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
library(forcats)
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
            "#FFA500", "#458B00", "#CD6090", "#FFAEB9", "#90EE90", "#5f9EA0", "#E6E6FA", "#8B7E66")
new_pal2 <- new_pal
names(new_pal2) <- paste0("c", 1:length(new_pal2))
superclone_col <- yarrr::piratepal("appletv")
names(superclone_col) <- NULL

pro_name <- "wdr_nanowell"
pro_name1 <- "wdr_nanowell1"
pro_name2 <- "wdr_nanowell2"
disp_file1 <- c("/volumes/lab/users/wet_lab/instruments/wafergen/Rui_WDR/filter_files/20211025_WDR_231_noNT_test/133444-2_Chip1_WellList.TXT")
disp_file2 <- c("/volumes/lab/users/wet_lab/instruments/wafergen/Rui_WDR/filter_files/20210819_low_dNTPs_high_dCTPs_test/133469-2_chip1_WellList.TXT")
# file.copy(disp_file1, paste0("./metrics/dispense_files/", pro_name1, "_", basename(disp_file1)))
# file.copy(disp_file2, paste0("./metrics/dispense_files/", pro_name2, "_", basename(disp_file2)))
disp1 <- read.delim(paste0("./metrics/dispense_files/", pro_name1, "_", basename(disp_file1)))
disp2 <- read.delim(paste0("./metrics/dispense_files/", pro_name2, "_", basename(disp_file2)))

pro_name_r <- paste0(pro_name, "_rna")
pro_name_d <- paste0(pro_name, "_dna")
###---RNA--------##########
pro_name_r <- paste0(pro_name, "_rna")

myfolder_name_p1 <- "20211029_WDR_MDA231_chip1R_chip2DR-noNT_SpRNA-ECIS44-cDNA2"
mysample_name_p1 <- "WDR_231_Chip1_cDNA_biotin"
raw_path_r1 <- c(paste0("/volumes/USR2/wangkl/wscDR/", myfolder_name_p1, "/data/", 
                        mysample_name_p1, "_exact_hg19_star_smart_2/star_2.7.5/",
                        mysample_name_p1, "_exact_hg19_star_smart_2_Solo.out/GeneFull/filtered_nodup"))
myfolder_name_p2 <- "20210821_WDR_lowdNTP_231DR_ART216orgR_SpRNA-ECIS42-SnubarARC-RNA2"
mysample_name_p2 <- "WDR_low_dNTPs_high_dCTPs_test_cDNA_231"
raw_path_r2 <- c(paste0("/volumes/USR2/wangkl/wscDR/", myfolder_name_p2, "/data/", 
                        mysample_name_p2, "_exact_hg19_star_smart_2/star_2.7.5/",
                        mysample_name_p2, "_exact_hg19_star_smart_2_Solo.out/GeneFull/filtered_nodup"))

#-----RNA processing-----
mtx_p1<- Read10X(data.dir = raw_path_r1)
mtx_p1 <- CreateSeuratObject(counts = mtx_p1, min.cells = 1, min.features = 3, project = pro_name1)
names_tmp1 <- mtx_p1@meta.data %>% rownames_to_column()
filt_cells_names_p1 <- disp1 %>% dplyr::filter(!(Sample %in% c("Neg Ctrl", "Pos Ctrl"))) %>%  dplyr::select(1,2,4) %>% mutate(row_col=paste(Row, Col, sep = "_")) %>% 
  dplyr::select(3:4) %>% left_join(wafer_match_list, by = c("row_col" = "coor")) %>% 
  inner_join(names_tmp1, by = c("RNA_Barcode" = "rowname")) %>% pull(RNA_Barcode)
dim(names_tmp1)
length(filt_cells_names_p1)
mtx_p1 <- subset(mtx_p1, cells = filt_cells_names_p1)

mtx_p2<- Read10X(data.dir = raw_path_r2)
mtx_p2 <- CreateSeuratObject(counts = mtx_p2, min.cells = 1, min.features = 3, project = pro_name2)
names_tmp2 <- mtx_p2@meta.data %>% rownames_to_column()
filt_cells_names_p2 <- disp2 %>% dplyr::filter(!(Sample %in% c("Neg Ctrl", "Pos Ctrl"))) %>%  dplyr::select(1,2,4) %>% mutate(row_col=paste(Row, Col, sep = "_")) %>% 
  dplyr::select(3:4) %>% left_join(wafer_match_list, by = c("row_col" = "coor")) %>% 
  inner_join(names_tmp2, by = c("RNA_Barcode" = "rowname")) %>% pull(RNA_Barcode)
dim(names_tmp2)
length(filt_cells_names_p2)
mtx_p2 <- subset(mtx_p2, cells = filt_cells_names_p2)

mtx2 <- merge(mtx_p1, y = mtx_p2, add.cell.ids = c("nanowell1", "nanowell2"), project = pro_name_r)
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

mtx2 <- subset(mtx2, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 10 )
mtx2 <- NormalizeData(object = mtx2, normalization.method = "LogNormalize", scale.factor = 10000)
mtx2 <- FindVariableFeatures(mtx2, selection.method = "vst")

all.genes <- rownames(mtx2)
library(future)
options(future.globals.maxSize= 40960 * 1024^2)
plan(strategy = "multicore", workers = 20)

mtx2 <- ScaleData(mtx2, features = all.genes, vars.to.regress = c("percent.mt","nCount_RNA","orig.ident"))#,"S.Score", "G2M.Score"))
mtx2 <- RunPCA(mtx2, features = all.genes)
ElbowPlot(mtx2, ndims = 50)

sel_dim <- 20
sel_res <- 0.15
mtx3_1 <- FindNeighbors(mtx2, dims = 1:sel_dim)
mtx3_2 <- FindClusters(mtx3_1, resolution = sel_res)
mtx3_1_tsu <- RunUMAP(mtx3_2, dims = 1:sel_dim)

p2 <- DimPlot(mtx3_1_tsu, reduction = "umap", label = T) #+ scale_color_manual(values = my_color_pal)
p2
p3 <- DimPlot(mtx3_1_tsu, reduction = "umap", label = T, group.by = "orig.ident")
p3
cowplot::ggsave2(paste0("./figures/", pro_name_r, "_UMAP_all_cells.pdf"), p2, width = 5, height = 4) 
cowplot::ggsave2(paste0("./figures/", pro_name_r, "_UMAP_all_cells_by_plates.pdf"), p3, width = 5, height = 4) 

saveRDS(mtx3_1_tsu, paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
mtx3_1_tsu <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))

####---------Gene and UMI number check--
mtx3_1_tsu.markers <- FindAllMarkers(object = mtx3_1_tsu, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
mtx3_1_tsu.markers %>% dplyr::group_by(cluster) %>% top_n(2, avg_log2FC)
top10 <- mtx3_1_tsu.markers %>% dplyr::group_by(cluster) %>% top_n(10, avg_log2FC)

pdf(paste0(c("./figures/"),pro_name_r, c("_singlet_Seurat_top10_Heatmap.pdf")), height = 10, width = 8)
DoHeatmap(mtx3_1_tsu, features = top10$gene)
dev.off()

exp_col <- jcolors::jcolors("pal2")[3:4]
names(exp_col) <- c("wdr_nanowell1", "wdr_nanowell2")
pdf(paste0(c("./figures/"),pro_name_r, c("_singlet_Seurat_top10_ComplexHeatmap.pdf")), height = 10, width = 8)
DoComHeatmapS2(mtx3_1_tsu, genelist = top10$gene, group.by2 = 'orig.ident', group.colors2 = exp_col, breaks=c(-1.5,0,1.5))
dev.off()

pdf(paste0(c("./figures/"),pro_name_r, c("_singlet_Seurat_clst_gene_umi_number.pdf")), height = 10, width = 8)
VlnPlot(mtx3_1_tsu, features = c("nCount_RNA","nFeature_RNA"))
dev.off()


library(future)
options(future.globals.maxSize= 40960 * 1024^2)
plan(strategy = "multicore", workers = 20)
de.markers <- FindMarkers(mtx3_1_tsu, ident.1 = 0, ident.2 = 1, assay = "RNA", only.pos = F, min.pct = 0, logfc.threshold = 0)

sig_de.markers <- de.markers %>% as.data.frame() %>% rownames_to_column() %>% filter(avg_log2FC > 0.485 | avg_log2FC < -0.485) %>% filter(p_val_adj < 0.05)
write.csv(sig_de.markers, paste0("./metrics/", pro_name_r, "_volcano_DE_sig.csv"))

sig_de.markers <- read.csv(paste0("./metrics/", pro_name_r, "_volcano_DE_sig.csv"), row.names = 1)

###Volcano plot wth EnhanceVolcano###
library(EnhancedVolcano)
#--FC: 0.485 == 1.4X fold-change
p1<-EnhancedVolcano(de.markers, lab = rownames(de.markers), x="avg_log2FC", y="p_val_adj", FCcutoff = 0.485, pCutoff = 0.05, titleLabSize = 20,
                      # title = paste0(i, c("_DE_MAST")), labSize = 4.0, pointSize = 2, col = c("black", "grey", "azure4", "red"),
                      xlim = c(-1.5, 1.5),
                      #ylim = c(0, 30),
                      # drawConnectors = T, 
                      widthConnectors = 0.2, colConnectors = "grey30",typeConnectors = "open",
                      endsConnectors = "last", lengthConnectors = unit(0.01, 'npc'),
                      axisLabSize = 20, gridlines.major = FALSE, gridlines.minor = FALSE, border="full",
                      legendLabels = c('NS','Log2FC','adjP < 0.05','Log2FC and adjP'),
                      legendLabSize=15, legendPosition="bottom", colAlpha = 1, xlab = bquote(~Log[2]~"fold change"), ylab = bquote(~-log[10]~adjusted~italic(P)))
p1
p2 <- ggrastr::rasterize(p1, layers='Point', dpi=350)
cowplot::ggsave2(paste0("./figures/",pro_name_r, "_volcano_DE.pdf"), p2, width = 4, height = 6)


#----combine with DNA meta data (only use cells that passed DNA QC)
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
dna_meta_temp <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::select("rna_bc2","superclones","subclones","sample")
rna_meta_temp <- mtx3_1_tsu@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  left_join(dna_meta_temp, by = c("rowname" = "rna_bc2")) %>% column_to_rownames() 
rna_meta_temp2 <- cbind(rna_meta_temp, mtx3_1_tsu@reductions$umap@cell.embeddings)

#---Plot UMAP colored by DNA subclones---
rna_col <- c(ggthemes::calc_pal()(length(levels(varbin_mtx_tumor_log2@colData$seurat_clusters))))

p1 <- ggplot(rna_meta_temp2,aes(x = UMAP_1, y = UMAP_2, fill = seurat_clusters)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = rna_col) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
cowplot::ggsave2(paste0("./figures/", pro_name_r, "_UMAP_all_cells_color_by_rna_clusters.pdf"), p1, width = 5.5, height = 4) 

p2 <- ggplot(rna_meta_temp2,aes(x = UMAP_1, y = UMAP_2, fill = superclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = superclone_col, na.translate =F) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p2
cowplot::ggsave2(paste0("./figures/", pro_name_r, "_UMAP_all_cells_color_by_dna_superclones.pdf"), p2, width = 5, height = 4) 

p3 <- ggplot(rna_meta_temp2,aes(x = UMAP_1, y = UMAP_2, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = new_pal, na.translate =F) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p3
cowplot::ggsave2(paste0("./figures/", pro_name_r, "_UMAP_all_cells_color_by_dna_subclones.pdf"), p3, width = 5, height = 4) 

exp_col <- jcolors::jcolors("pal2")[3:4]
names(exp_col) <- c("wdr_nanowell1", "wdr_nanowell2")
p3 <- ggplot(rna_meta_temp2,aes(x = UMAP_1, y = UMAP_2, fill = orig.ident)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = exp_col, na.translate =F) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p3
cowplot::ggsave2(paste0("./figures/", pro_name_r, "_UMAP_all_cells_color_by_rna_experiments.pdf"), p3, width = 5.5, height = 4) 

#----RNA QC matrics---#######
mtx3_1_tsu <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))

mtx3_1_tsu@meta.data %>% filter(orig.ident == "wdr_nanowell1") %>% pull(nCount_RNA) %>% mean()
mtx3_1_tsu@meta.data %>% filter(orig.ident == "wdr_nanowell2") %>% pull(nCount_RNA) %>% mean()

mtx3_1_tsu@meta.data %>% filter(orig.ident == "wdr_nanowell1") %>% pull(nFeature_RNA) %>% mean()
mtx3_1_tsu@meta.data %>% filter(orig.ident == "wdr_nanowell2") %>% pull(nFeature_RNA) %>% mean()

#----InferCNA using copykat--WDR-MDA231 data---#########
mtx3_1_tsu <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
write.table(mtx3_1_tsu@assays$RNA@counts, file=paste0("./objects/",pro_name_r, "_seurat_rawdata.txt"), quote=F, sep="\t", row.names= T)
###--Run in the server
# devtools::install_github("navinlabcode/copykat") 
setwd("/volumes/USR2/wangkl/wscDR/wdr_analysis/data/copykat/wdr_mda231")
dir.create("/volumes/USR2/wangkl/wscDR/wdr_analysis/data/copykat/wdr_mda231_2")
setwd("/volumes/USR2/wangkl/wscDR/wdr_analysis/data/copykat/wdr_mda231_2")

library(copykat)
copykat_results <- copykat(rawmat=as.matrix(mtx3_1_tsu@assays$RNA@counts), id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, 
                        cell.line = "yes",
                        sam.name=pro_name_r, distance="euclidean", norm.cell.names="",
                        output.seg="FLASE", plot.genes="TRUE", genome="hg20",n.cores=80)

#######-------NEW-------Map cnv (inferred from all cellranger data) to RNA-seq#####
library(gplots)
library(RColorBrewer)
library("ComplexHeatmap")
library("circlize")
library("tidyverse")
library("dendextend")

final <- read.delim("./wdr_nanowell_rna_copykat_CNA_results.txt", header=T)
# clst_ckat <- readRDS("./wdr_nanowell_rna_copykat_clustering_results.rds")
dim(final)

setwd(wd)

###----covert into hg19 genome coordinates---##
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
hg19_position <- varbin_mtx_tumor_log2@rowRanges %>% as.data.frame() %>% dplyr::select("seqnames","start","abspos")
final$chrom <- hg19_position$seqnames
final$chrompos <- hg19_position$start
final$abspos <- hg19_position$abspos

data <- as.matrix(t(final[, 4:ncol(final)]))
chr <- as.integer(final$chrom) %% 2+1
chr2<-as.data.frame(chr)

hc_dist <- amap::Dist(data, method = "euclidean", nbproc = 80)  #----Most time consuming part
saveRDS(hc_dist, file = paste0("/volumes/USR2/wangkl/wscDR/wdr_analysis/data/copykat/wdr_mda231/",pro_name_r, "_hc_dist.rds"))
# hc_dist <- readRDS(file = paste0(pro_name, "_hc_dist.rds"))
hc <- fastcluster::hclust(hc_dist, method = "ward.D2")


#---Add chr # on top of chr_bar--
a <- c()
for (i in 1:23) {
  len <- table(final$chrom)[i]
  if(i<10|i>22){
    a1 <- rep("", times=round(len*0.6))
    a2 <- rep("", times=(len-round(len*0.6)-1))
  }else if (i<21){
    a1 <- rep("", times=round(len*0.8))
    a2 <- rep("", times=(len-round(len*0.8)-1))
  }else{
    a1 <- rep("", times=len)
  }
  if(i<21){
    a <- c(a, a1, i, a2)
  } else if (i==21 | i==22){
    a <- c(a, a1)
  } else {
    a <- c(a, a1, "X", a2)}
}

#######----------------------Plot Complexheatmap----
ha_col=HeatmapAnnotation(foo=anno_text(a, rot = 0, gp = gpar(fontsize =10)), df =chr2, col = list(chr=c("1"="black", "2"="grey")), show_legend = F, 
                         annotation_name_side = "left")

pdf(paste0(c("./figures/"),pro_name_r, c("_copykat_wardD2_complexHeatmap.pdf")), height = 5, width = 7)
Heatmap(data, name = "value", top_annotation = ha_col, left_annotation = ha_row, cluster_rows = hc, show_row_dend = TRUE,
        row_dend_width=unit(25, "mm"), cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE, raster_quality = 5)
dev.off()




#----InferCNA using copykat--10x genomics--MDA231 data---#########
mtx3_1_tsu_10x <- readRDS("/volumes/USR1/ruiye/10x/2020ANALYSIS/WDR/WDR_Paper/cellline_QC_231/wafargen/tenx_231.rds")
# random2000_10x_name <- mtx3_1_tsu_10x@meta.data %>% rownames() %>% sample(2000)
# mtx3_1_tsu_10x_2k <- subset(mtx3_1_tsu_10x, cells = random2000_10x_name)
set.seed(30)
write.table(mtx3_1_tsu_10x@assays$RNA@counts[,sample(ncol(mtx3_1_tsu_10x@assays$RNA@counts), 2000)], 
            file=paste0("./objects/",pro_name_r, "_10xgenomics_random2K_seurat_rawdata.txt"),quote=F, sep="\t", row.names= T)

mda231_10x_2k <- read.table(paste0("./objects/",pro_name_r, "_10xgenomics_random2K_seurat_rawdata.txt"))
# write.table(mtx3_1_tsu_10x_2k@assays$RNA@counts, file=paste0("./objects/",pro_name_r, "_10xgenomics_seurat_rawdata.txt"), quote=F, sep="\t", row.names= T)
###--Run in the server
# devtools::install_github("navinlabcode/copykat") 
wd_10x <- c("/volumes/USR2/wangkl/wscDR/wdr_analysis/data/copykat/wdr_mda231_10xgenomics")
# dir.create(wd_10x)
setwd(wd_10x)
library(copykat)

copykat_results <- copykat(rawmat=mda231_10x_2k, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, 
                           cell.line = "yes",
                           sam.name="mda231_10x", distance="euclidean", norm.cell.names="",
                           output.seg="FLASE", plot.genes="TRUE", genome="hg20",n.cores=80)

#######-------NEW-------Map cnv (inferred from all cellranger data) to RNA-seq#####
library(gplots)
library(RColorBrewer)
library("ComplexHeatmap")
library("circlize")
library("tidyverse")
library("dendextend")

final <- read.delim("./mda231_10x_copykat_CNA_results.txt", header=T)
# clst_ckat <- readRDS("./wdr_nanowell_rna_copykat_clustering_results.rds")
dim(final)

setwd(wd)

###----covert into hg19 genome coordinates---##
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
hg19_position <- varbin_mtx_tumor_log2@rowRanges %>% as.data.frame() %>% dplyr::select("seqnames","start","abspos")
final$chrom <- hg19_position$seqnames
final$chrompos <- hg19_position$start
final$abspos <- hg19_position$abspos

data <- as.matrix(t(final[, 4:ncol(final)]))
chr <- as.integer(final$chrom) %% 2+1
chr2<-as.data.frame(chr)

hc_dist <- amap::Dist(data, method = "euclidean", nbproc = 80)  #----Most time consuming part
saveRDS(hc_dist, file = paste0("/volumes/USR2/wangkl/wscDR/wdr_analysis/data/copykat/wdr_mda231_10xgenomics/mda231_10X_2k_hc_dist.rds"))
# hc_dist <- readRDS(file = paste0(pro_name, "_hc_dist.rds"))
hc <- fastcluster::hclust(hc_dist, method = "ward.D2")


#---Add chr # on top of chr_bar--
a <- c()
for (i in 1:23) {
  len <- table(final$chrom)[i]
  if(i<10|i>22){
    a1 <- rep("", times=round(len*0.6))
    a2 <- rep("", times=(len-round(len*0.6)-1))
  }else if (i<21){
    a1 <- rep("", times=round(len*0.8))
    a2 <- rep("", times=(len-round(len*0.8)-1))
  }else{
    a1 <- rep("", times=len)
  }
  if(i<21){
    a <- c(a, a1, i, a2)
  } else if (i==21 | i==22){
    a <- c(a, a1)
  } else {
    a <- c(a, a1, "X", a2)}
}

#######----------------------Plot Complexheatmap----
ha_col=HeatmapAnnotation(foo=anno_text(a, rot = 0, gp = gpar(fontsize =10)), df =chr2, col = list(chr=c("1"="black", "2"="grey")), show_legend = F, 
                         annotation_name_side = "left")

pdf(paste0(c("./figures/mda231_10xgenomics_random2K_copykat_wardD2_complexHeatmap.pdf")), height = 5, width = 7)
Heatmap(data, name = "value", top_annotation = ha_col, cluster_rows = hc, show_row_dend = TRUE,
        row_dend_width=unit(25, "mm"), cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE, raster_quality = 5)
dev.off()




###---DNA--------#########
pro_name_d <- paste0(pro_name, "_dna")
raw_path1 <- c("/volumes/USR2/wangkl/wscDR/20211029_WDR_MDA231_chip1R_chip2DR-noNT_SpRNA-ECIS44-cDNA2/data/waferDR_MDA231_chip1_DNA_1M/waferDR_MDA231_chip1-DNA-1M_200/res_200_k")
raw_path2 <- c("/volumes/USR2/wangkl/wscDR/20210821_WDR_lowdNTP_231DR_ART216orgR_SpRNA-ECIS42-SnubarARC-RNA2/data/waferDR_231_DNA_1M/waferDR_231_DNA_1M_200/res_200_k")

#-----DNA processing------####
#---nanowell1--
varbin_mtx1 <- readVarbin(raw_path1, remove_Y = TRUE)
varbin_mtx1 # 1655
filter_cells1 = read.table(paste0("./metrics/", pro_name1, "_filtered_bincounts.txt"), header = T)
filt_cells_names1 <- colnames(filter_cells1)[4:ncol(filter_cells1)] %>% janitor::make_clean_names()

filt_cells_names_tmp1 <- filt_cells_names1 %>% as.data.frame() %>% mutate(cell = stringr::str_extract(., "c\\d+")) %>% dplyr::rename(cell_name = ".")
filt_cells_names_dis1 <- disp1 %>% dplyr::filter(!(Sample %in% c("Neg Ctrl", "Pos Ctrl"))) %>%  dplyr::select(1,2,4) %>% mutate(row_col=paste(Row, Col, sep = "_")) %>% 
  dplyr::select(3:4) %>% left_join(wafer_match_list, by = c("row_col" = "coor")) %>% mutate(cell2 = tolower(Cell)) %>% 
  inner_join(filt_cells_names_tmp1, by = c("cell2" = "cell")) %>% pull(cell_name)
length(filt_cells_names1) # 1576
length(filt_cells_names_dis1) # 1569

varbin_mtx_tumor1 <- varbin_mtx1[,filt_cells_names_dis1]
varbin_mtx_tumor1@colData

#---nanowell2--
varbin_mtx2 <- readVarbin(raw_path2, remove_Y = TRUE)
varbin_mtx2 # 1458
filter_cells2 = read.table(paste0("./metrics/", pro_name2, "_filtered_bincounts.txt"), header = T)
filt_cells_names2 <- colnames(filter_cells2)[4:ncol(filter_cells2)] %>% janitor::make_clean_names()

filt_cells_names_tmp2 <- filt_cells_names2 %>% as.data.frame() %>% mutate(cell = stringr::str_extract(., "c\\d+")) %>% dplyr::rename(cell_name = ".")
filt_cells_names_dis2 <- disp2 %>% dplyr::filter(!(Sample %in% c("Neg Ctrl", "Pos Ctrl"))) %>%  dplyr::select(1,2,4) %>% mutate(row_col=paste(Row, Col, sep = "_")) %>% 
  dplyr::select(3:4) %>% left_join(wafer_match_list, by = c("row_col" = "coor")) %>% mutate(cell2 = tolower(Cell)) %>% 
  inner_join(filt_cells_names_tmp2, by = c("cell2" = "cell")) %>% pull(cell_name)
length(filt_cells_names2) # 1367
length(filt_cells_names_dis2) # 1367

varbin_mtx_tumor2 <- varbin_mtx2[,filt_cells_names_dis2]
varbin_mtx_tumor2@colData

#------merge 2 nanowell data sets----
varbin_mtx_tumor1@colData$experiment <- "nanowell1"
varbin_mtx_tumor2@colData$experiment <- "nanowell2"

varbin_mtx_tumor <- cbind(varbin_mtx_tumor1, varbin_mtx_tumor2)
#---rename cell names---
newname <- rownames(varbin_mtx_tumor@colData)
exper <- varbin_mtx_tumor@colData$experiment
rownames(varbin_mtx_tumor@colData) <- paste(rownames(varbin_mtx_tumor@colData), varbin_mtx_tumor@colData$experiment, sep = "_")
varbin_mtx_tumor@colData$sample <- paste(rownames(varbin_mtx_tumor@colData), varbin_mtx_tumor@colData$experiment, sep = "_")

saveRDS(varbin_mtx_tumor, file = paste0("objects/", pro_name_d, c("_filtered_copykit_before_rename.rds")))
varbin_mtx_tumor <- readRDS(paste0("objects/", pro_name_d, c("_filtered_copykit_before_rename.rds")))

varbin_mtx_tumor2 <- varbin_mtx_tumor
varbin_mtx_tumor2 <- countBreakpoints(varbin_mtx_tumor2)

#---add meta info--
name_dna <- varbin_mtx_tumor@colData %>% as.data.frame() %>% dplyr::select(c("sample", "experiment")) %>% 
  mutate(newname = toupper(stringr::str_extract(sample, "c\\d+"))) %>% left_join(wafer_match_list, by = c("newname" = "Cell")) %>% 
  mutate(rna_bc2 = paste(experiment, RNA_Barcode, sep = "_"))

rna_file <- mtx3_1_tsu
DNA_rna_meta <- rna_file@meta.data %>% rownames_to_column() %>% left_join(name_dna, ., by = c("rna_bc2" = "rowname")) %>% dplyr::select(-1)

varbin_mtx_tumor@colData <- cbind(varbin_mtx_tumor@colData, DNA_rna_meta)
varbin_mtx_tumor_log <- logNorm(varbin_mtx_tumor, transform = "log2")
saveRDS(varbin_mtx_tumor_log, file = paste0("./objects/", pro_name_d, c("_filtered_copykit_rnameta.rds")))
varbin_mtx_tumor_log <- readRDS(paste0("./objects/", pro_name_d, c("_filtered_copykit_rnameta.rds")))

#----UMAP-----
# varbin_mtx_tumor_log2 <- varbin_mtx_tumor_log
set.seed(31)
near_nb <- 15
umap_data <-  data.frame(uwot::umap(log(t(varbin_mtx_tumor_log@assays@data$segment_ratios)), 
                                    metric = "manhattan", min_dist = 0.1, spread = 1, n_neighbors = near_nb))

umap_data2 <- umap_data %>% dplyr::rename("UMAP_1" = "X1", "UMAP_2" = "X2")
varbin_mtx_tumor_log@int_colData$reducedDims <- umap_data2
varbin_mtx_tumor_log@colData <- cbind(varbin_mtx_tumor_log@colData, umap_data2)

#---find clusters----
#--superclone---
g_major <-bluster::makeSNNGraph(umap_data2, k = 70)
superclones <- as.factor(paste0("s", igraph::membership(igraph::components(g_major))))
varbin_mtx_tumor_log@colData$superclones <- superclones

ggplot(as.data.frame(varbin_mtx_tumor_log@colData),aes(x = UMAP_1, y = UMAP_2, fill = superclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = superclone_col) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())

#---subclone---
set.seed(17)
hdb_data <- dbscan::hdbscan(umap_data2, minPts = nrow(umap_data2)*0.015)
subclones <- paste0("c", as.character(hdb_data$cluster))

varbin_mtx_tumor_log@colData <- cbind(varbin_mtx_tumor_log@colData, subclones)
varbin_mtx_tumor_log@colData$subclones <- factor(varbin_mtx_tumor_log@colData$subclones, 
                                                 levels = paste0("c", 0:(length(table(subclones))-1)))

p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = c("grey",new_pal)) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_passDNAQC_cells_umap.pdf"), p1, width = 5, height = 4)

#----remove outlier cluster c0 and cluster that have less 3 cells.
clone_num <- table(varbin_mtx_tumor_log@colData$subclones)
kept_clones <- as.data.frame(clone_num) %>% dplyr::filter(Var1 != "c0") %>% dplyr::filter(Freq > 3) 

varbin_mtx_tumor_log2 <- varbin_mtx_tumor_log[, subclones %in% kept_clones$Var1]
varbin_mtx_tumor_log2@colData

p0 <- ggplot(as.data.frame(varbin_mtx_tumor_log@colData),aes(x = UMAP_1, y = UMAP_2, fill = superclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = superclone_col) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p0
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_final_filtered_cells_umap_superclones.pdf"), p1, width = 5, height = 4)

p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log2@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = new_pal) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_final_filtered_cells_umap.pdf"), p1, width = 5, height = 4)

exp_col <- jcolors::jcolors("pal2")[3:4]
names(exp_col) <- c("nanowell1", "nanowell2")
p2 <- ggplot(shuf(as.data.frame(varbin_mtx_tumor_log2@colData)),aes(x = UMAP_1, y = UMAP_2, fill = experiment)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = exp_col) + theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p2
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_final_filtered_cells_umap_plate.pdf"), p2, width = 5, height = 4)

rna_col <- c(ggthemes::calc_pal()(length(levels(varbin_mtx_tumor_log2@colData$seurat_clusters))),"grey10")
p3 <- ggplot(as.data.frame(varbin_mtx_tumor_log2@colData),aes(x = UMAP_1, y = UMAP_2, fill = seurat_clusters)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = rna_col, na.translate =F) + theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p3
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_final_filtered_cells_umap_subclone.pdf"), p3, width = 5, height = 3.5)

p4 <- plot_grid(plotlist=list(p0,p1,p2,p3), ncol=4, align='h')
p4
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_final_filtered_cells_umap_merge.pdf"), p4, width = 16, height = 2.8)

p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log2@colData)) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = superclones), size=8) + scale_color_manual(values = superclone_col) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill = subclones), shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = new_pal) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_final_filtered_cells_umap_subandsuper_clones.pdf"), p1, width = 5, height = 4)

p2 <- ggplot(as.data.frame(varbin_mtx_tumor_log2@colData)) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = superclones), size=8) + scale_color_manual(values = superclone_col) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill = seurat_clusters), shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = rna_col, na.translate =F) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p2
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_final_filtered_cells_umap_subandsuper_clones_RNA_clusters.pdf"), p2, width = 5, height = 4)

saveRDS(varbin_mtx_tumor_log2, file = paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
#------heatmap with clonal status bar----#####
mtx_srt <- as.data.frame(varbin_mtx_tumor_log2@colData) %>%
  arrange(factor(subclones, levels = paste0("c", 1:length(unique(varbin_mtx_tumor_log2@colData$subclones)))))
library(dplyr)
shuf_name <- mtx_srt %>% mutate(sample2 = paste(sample, experiment, sep = "_")) %>% group_by(subclones) %>% 
  group_modify(~ mutate(., sample2 = sample(sample2))) %>%
  ungroup() %>% pull(sample2)
rownames(mtx_srt) <-paste(mtx_srt$sample, mtx_srt$experiment, sep = "_")
mtx_srt <- mtx_srt[shuf_name,]
ht_mtx <- log2(t(varbin_mtx_tumor_log2@assays@data$segment_ratios))[mtx_srt$sample,]

#----annotation bar--
anno_mtx <- mtx_srt %>% mutate(rna_clst = ifelse(is.na(seurat_clusters), "RNA_miss",seurat_clusters)) %>% 
  dplyr::select(c("superclones","subclones", "rna_clst", "experiment"))
rownames(anno_mtx) <- NULL

sclst_col <- superclone_col[1:length(unique(anno_mtx$superclones))]
names(sclst_col) <- paste0("s", 1:length(unique(anno_mtx$superclones)))
clst_col <- new_pal[1:length(unique(anno_mtx$subclones))]
names(clst_col) <- paste0("c", 1:length(unique(anno_mtx$subclones)))

rna_col <- c(ggthemes::calc_pal()(length(unique(anno_mtx$rna_clst))-1),"grey50")
names(rna_col) <- sort(unique(anno_mtx$rna_clst))
exp_col <- c("#9970ab", "#66c2a5")
names(exp_col) <- c("nanowell1", "nanowell2")
#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue3", "white", "firebrick3"))
ha_row2=rowAnnotation(df = anno_mtx, col = list(superclones=sclst_col, subclones=clst_col, rna_clst=rna_col, experiment=exp_col), 
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

#----clonal status--
low_cutoff <- -0.15
up_cutoff <-  0.15
cell_cutoff <- 0.88
neu_cutoff <- 0.90
my_sample <- clonality_log_trinary_neu(log_ratio_df=ht_mtx, lower_cutoff = low_cutoff, upper_cutoff = up_cutoff, 
                                         cell_pct = cell_cutoff, neu_pct = neu_cutoff)
table(my_sample)
# cs <- as.data.frame(my_sample)
sig_de.markers <- read.csv(paste0("./metrics/", pro_name_r, "_volcano_DE_sig.csv"), row.names = 1)
gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% sig_de.markers$rowname)

gene_bin_sel_sample <- as.data.frame(my_sample) %>% rownames_to_column() %>% dplyr::mutate(rowname = as.numeric(rowname)) %>% 
  dplyr::right_join(gene_bin_sel, by = c("rowname" = "pos")) %>% 
  dplyr::mutate(my_color = plyr::mapvalues(my_sample, from = c("cCNA", "sCNA","neu"), to = c("#424451","#F47F20","grey92"))) 

my_col = structure(c("#424451","grey92","#F47F20"), names = c("cCNA","neu","sCNA"))
my_sample_df <- as.data.frame(t(my_sample))

cs <- as.data.frame(my_sample)
ha_bottom = columnAnnotation(df = cs, col = list(my_sample=c("cCNA"="#424451", "neu" = "grey92","sCNA"="#F47F20")), 
                             clonal_state = anno_mark(at=gene_bin_sel_sample$rowname, labels = gene_bin_sel_sample$gene, 
                                                      side = "bottom", labels_gp = gpar(fontsize = 14, col = gene_bin_sel_sample$my_color)))

pdf(paste0("./figures/", pro_name_d, "_complexHeatmap_DR.pdf"), height = 10, width = 10)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE,
        row_split = anno_mtx$subclones,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
        column_title = paste0(pro_name, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff),
        use_raster = T, raster_quality = 5, col = col_vec,
        heatmap_legend_param = list(title = "Log2 (Ratio)",
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row2, bottom_annotation = ha_bottom)
dev.off()

#----Plot gene freq with and without dosage effect---
df <- gene_bin_sel_sample
p1 <- ggplot(df, aes(x=my_sample, fill=my_sample)) + geom_bar() + theme_classic() + 
  scale_fill_manual(values=c("cCNA"="#7f7f7f", "sCNA"="grey70")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  labs(x="my_sample", y="gene count")
p1
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_barplot_gene_freq_dosage_effect.pdf"), p1, width = 2.5, height = 2.5)

#----------Calculate computational ploidy----####
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

saveRDS(varbin_mtx_tumor_log3, file = paste0("objects/", pro_name_d, c("_tumor_integerCN_copykit.rds")))
# varbin_mtx_tumor_log3 <- readRDS(paste0("objects/", pro_name_d, c("_tumor_integerCN_copykit.rds")))

#----------Integer consensus Heatmap----#####
cs_mtx <- t(varbin_mtx_tumor_log3@consensus)
my_order <- unique(varbin_mtx_tumor_log3@colData$subclones) %>% as.character() %>% sort()
cs_mtx_order <- cs_mtx[my_order,]
topleft(cs_mtx_order)

clst_col_cs <- new_pal[1:length(unique(varbin_mtx_tumor_log3@colData$subclones))]
subclone_cs <- unique(varbin_mtx_tumor_log3@colData$subclones)
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

pdf(paste0("./figures/", pro_name_d, "_merged_consensus_integer_complexHeatmap.pdf"), height = 6, width = 10)
Heatmap(as.matrix(cs_mtx_order), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        name = "csheatmap", show_row_names = T, row_names_side = "left",
        show_column_names = F, column_title = paste0(pro_name_d, "_consensusHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row_cs, bottom_annotation = ha_bottom_cs)
dev.off()


#----------MEDICC2 Tree----######
varbin_mtx_tumor_log3 <- readRDS(paste0("objects/", pro_name_d, c("_tumor_integerCN_copykit.rds")))
## make medicc input matrix
medicc_input <- cbind(varbin_mtx_tumor_log3@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_tumor_log3@consensus) %>%
  mutate(diploid=2) %>%          ## assign diploid cell CN profile to be used as root in the medicc tree
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
#---ME tree, need two diploid cell populations to make diploid branch on the most left side.
medic_dis <- read_tsv(paste0(tsv_out_path, pro_name_d, "_medicc2_input_pairwise_distances.tsv"))
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
  theme(legend.position = "none") + ggtree::geom_rootpoint() + ggtree::theme_tree2()
treeplt

cowplot::ggsave2(paste0("./figures/", pro_name_d, "_medicc2_ME_tree_merged_subclones_lumhr.pdf"), treeplt, width = 8, height = 4.5)



#---DE analysis of subclones 1 versus 7 ---#####
obj_rna <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
dna_meta_temp <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::select("rna_bc2","superclones","subclones","sample")

obj_rna2 <- obj_rna
obj_rna2@meta.data <- obj_rna2@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  left_join(dna_meta_temp, by = c("rowname" = "rna_bc2")) %>% column_to_rownames() 

clone_num1 <- 1
clone_num2 <- 7
super_clone1 <- paste0("c", clone_num1)
super_clone2 <- paste0("c", clone_num2)
super_name1 <- paste(super_clone1, collapse = "_")
super_name2 <- paste0(super_clone2, collapse = "_")
obj_rna2@meta.data <- obj_rna2@meta.data %>% mutate(temp = ifelse(subclones %in% super_clone1, "super_c1", 
                                                                  ifelse(subclones %in% super_clone2, "super_c2", "others")))

de.markers <- FindMarkers(obj_rna2, ident.1 = "super_c1", ident.2 = "super_c2", group.by = 'temp', assay = "RNA", 
                          only.pos = F, min.pct = 0, logfc.threshold = 0)
de.markers
###Volcano plot wth EnhanceVolcano###
library(EnhancedVolcano)
#--FC: 0.485 == 1.4X fold-change
p1<-EnhancedVolcano(de.markers, lab = rownames(de.markers), x="avg_log2FC", y="p_val_adj", FCcutoff = 0.485, pCutoff = 0.05, titleLabSize = 20,
                    # title = paste0(i, c("_DE_MAST")), 
                    labSize = 3.0, 
                    # pointSize = 2, col = c("black", "grey", "azure4", "red"),
                    xlim = c(-1.5, 1.5),
                    ylim = c(0, 35),
                    # drawConnectors = T, 
                    widthConnectors = 0.2, colConnectors = "grey30",typeConnectors = "open",
                    endsConnectors = "last", lengthConnectors = unit(0.01, 'npc'),
                    axisLabSize = 20, gridlines.major = FALSE, gridlines.minor = FALSE, border="full",
                    legendLabels = c('NS','Log2FC','adjP < 0.05','Log2FC and adjP'),
                    legendLabSize=15, legendPosition="bottom", colAlpha = 1, xlab = bquote(~Log[2]~"fold change"), ylab = bquote(~-log[10]~adjusted~italic(P)))
p1
p2 <- ggrastr::rasterize(p1, layers='Point', dpi=350)
cowplot::ggsave2(paste0("./figures/",pro_name_r, "_",super_name1, "_vs_", super_name2, "_volcano_DE.pdf"), p2, width = 4, height = 6)

#------heatmap with dosage effect status bar----#####
#---define copy number dosage effect: Due to the bad RNA-seq data quality, probably just define yes or no.
#---For yes: all subclones have same copy number, otherwise no.
varbin_mtx_tumor_log3 <- readRDS(paste0("objects/", pro_name_d, c("_tumor_integerCN_copykit.rds")))
cs_mtx <- t(varbin_mtx_tumor_log3@consensus)
my_order <- c(super_clone1, super_clone2)

cs_mtx <- cs_mtx %>% as.data.frame() %>% dplyr::filter(rownames(cs_mtx) %in% my_order)
cs_mtx_order <- cs_mtx[my_order,]
topleft(cs_mtx_order)

#----annotation bar--
anno_mtx <- rownames(cs_mtx_order) %>% as.data.frame() %>% dplyr::rename(subclones=".") %>% 
  mutate(superclones = ifelse(subclones %in% super_clone1, "superclone1", "superclone2")) %>% 
  dplyr::select(c("subclones", "superclones"))
rownames(anno_mtx) <- NULL

clst_col <- c(new_pal[clone_num1], new_pal[clone_num2])
names(clst_col) <- anno_mtx$subclones

sclst_col <- superclone_col[1:length(unique(anno_mtx$superclones))]
names(sclst_col) <- unique(anno_mtx$superclones) %>% sort()

#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")
#-----annotated genes---
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
de_genes_up <- de.markers %>% dplyr::filter(p_val_adj<0.05 & avg_log2FC > 0.485)
de_genes_down <- de.markers %>% dplyr::filter(p_val_adj<0.05 & avg_log2FC < -0.485)

dos_state <- apply(cs_mtx_order, 2, function(row) {ifelse(length(unique(row)) == 1, "cCNA", "sCNA")})
dos_state_df <- dos_state %>% as.data.frame() %>% dplyr::rename(dos_state = ".")

gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% rownames(rbind(de_genes_up, de_genes_down))) %>% 
  mutate(col=ifelse(dos_state[pos] == "sCNA", "#F47F20", "#424451"))

ha_bottom_cs = columnAnnotation(df = dos_state_df, col = list(dos_state=c("sCNA"="#F47F20", "cCNA"="#424451")), 
                                exp_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
                                                      side = "bottom", labels_gp = gpar(fontsize = 10, col = gene_bin_sel$col)))

# ha_bottom_cs = columnAnnotation(exp_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
#                                                       side = "bottom", labels_gp = gpar(fontsize = 14, col = gene_bin_sel$col)))

heat_col <- pals::ocean.balance(14)[9:14]
if(max(cs_mtx)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(cs_mtx)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(cs_mtx)-6))
}
names(col_vec) <- 0:max(cs_mtx)

ha_row_cs=rowAnnotation(df = anno_mtx, col = list(subclones=clst_col, superclones=sclst_col),
                        show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/", pro_name_d, "_complexHeatmap_tumor_superclone_",super_name1, "_vs_", super_name2, "_with_DE_dosage_annotations.pdf"), height = 2.5, width = 10)
Heatmap(as.matrix(cs_mtx_order), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        name = "csheatmap", show_row_names = T, row_names_side = "left",
        show_column_names = F, column_title = paste0(pro_name_d, "_consensusHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row_cs, bottom_annotation = ha_bottom_cs)
dev.off()


#----Plot gene freq with and without dosage effect---
gene_bin_sel_sample <- gene_bin_sel %>% dplyr::mutate(my_dos = plyr::mapvalues(col, from = c("#424451","#F47F20"), to = c("cCNA", "sCNA"))) 

df <- gene_bin_sel_sample
p1 <- ggplot(df, aes(x=my_dos, fill=my_dos)) + geom_bar() + theme_classic() + 
  scale_fill_manual(values=c("sCNA"="#F47F20", "cCNA"="#424451")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  labs(x="my_dos", y="gene count")
p1
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_barplot_gene_freq_dosage_effect_", super_name1, "_vs_", super_name2, ".pdf"), p1, width = 2.5, height = 2.5)

#---DE analysis of subclones 2 versus 7 ---#####
obj_rna <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
dna_meta_temp <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::select("rna_bc2","superclones","subclones","sample")
  
obj_rna2 <- obj_rna
obj_rna2@meta.data <- obj_rna2@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  left_join(dna_meta_temp, by = c("rowname" = "rna_bc2")) %>% column_to_rownames() 

clone_num1 <- 2
clone_num2 <- 7
super_clone1 <- paste0("c", clone_num1)
super_clone2 <- paste0("c", clone_num2)
super_name1 <- paste(super_clone1, collapse = "_")
super_name2 <- paste0(super_clone2, collapse = "_")
obj_rna2@meta.data <- obj_rna2@meta.data %>% mutate(temp = ifelse(subclones %in% super_clone1, "super_c1", 
                                                                    ifelse(subclones %in% super_clone2, "super_c2", "others")))
  
de.markers <- FindMarkers(obj_rna2, ident.1 = "super_c1", ident.2 = "super_c2", group.by = 'temp', assay = "RNA", 
                            only.pos = F, min.pct = 0, logfc.threshold = 0)
de.markers
###Volcano plot wth EnhanceVolcano###
library(EnhancedVolcano)
#--FC: 0.485 == 1.4X fold-change
p1<-EnhancedVolcano(de.markers, lab = rownames(de.markers), x="avg_log2FC", y="p_val_adj", FCcutoff = 0.485, pCutoff = 0.05, titleLabSize = 20,
                    # title = paste0(i, c("_DE_MAST")), 
                    labSize = 3.0, 
                    # pointSize = 2, col = c("black", "grey", "azure4", "red"),
                    # xlim = c(-1.5, 1.5),
                    ylim = c(0, 6),
                    # drawConnectors = T, 
                    widthConnectors = 0.2, colConnectors = "grey30",typeConnectors = "open",
                    endsConnectors = "last", lengthConnectors = unit(0.01, 'npc'),
                    axisLabSize = 20, gridlines.major = FALSE, gridlines.minor = FALSE, border="full",
                    legendLabels = c('NS','Log2FC','adjP < 0.05','Log2FC and adjP'),
                    legendLabSize=15, legendPosition="bottom", colAlpha = 1, xlab = bquote(~Log[2]~"fold change"), ylab = bquote(~-log[10]~adjusted~italic(P)))
p1
p2 <- ggrastr::rasterize(p1, layers='Point', dpi=350)
cowplot::ggsave2(paste0("./figures/",pro_name_r, "_",super_name1, "_vs_", super_name2, "_volcano_DE.pdf"), p2, width = 4, height = 6)
  
#------heatmap with dosage effect status bar----#####
#---define copy number dosage effect: Due to the bad RNA-seq data quality, probably just define yes or no.
#---For yes: all subclones have same copy number, otherwise no.
cs_mtx <- t(varbin_mtx_tumor_log3@consensus)
my_order <- c(super_clone1, super_clone2)

cs_mtx <- cs_mtx %>% as.data.frame() %>% dplyr::filter(rownames(cs_mtx) %in% my_order)
cs_mtx_order <- cs_mtx[my_order,]
topleft(cs_mtx_order)

#----annotation bar--
anno_mtx <- rownames(cs_mtx_order) %>% as.data.frame() %>% dplyr::rename(subclones=".") %>% 
  mutate(superclones = ifelse(subclones %in% super_clone1, "superclone1", "superclone2")) %>% 
  dplyr::select(c("subclones", "superclones"))
rownames(anno_mtx) <- NULL

clst_col <- c(new_pal[clone_num1], new_pal[clone_num2])
names(clst_col) <- anno_mtx$subclones

sclst_col <- superclone_col[1:length(unique(anno_mtx$superclones))]
names(sclst_col) <- unique(anno_mtx$superclones) %>% sort()

#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")
#-----annotated genes---
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
de_genes_up <- de.markers %>% dplyr::filter(p_val_adj<0.05 & avg_log2FC > 0.485)
de_genes_down <- de.markers %>% dplyr::filter(p_val_adj<0.05 & avg_log2FC < -0.485)

dos_state <- apply(cs_mtx_order, 2, function(row) {ifelse(length(unique(row)) == 1, "cCNA", "sCNA")})
dos_state_df <- dos_state %>% as.data.frame() %>% dplyr::rename(dos_state = ".")

gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% rownames(rbind(de_genes_up, de_genes_down))) %>% 
  mutate(col=ifelse(dos_state[pos] == "sCNA", "#F47F20", "#424451"))

ha_bottom_cs = columnAnnotation(df = dos_state_df, col = list(dos_state=c("sCNA"="#F47F20", "cCNA"="#424451")), 
                                exp_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
                                                      side = "bottom", labels_gp = gpar(fontsize = 10, col = gene_bin_sel$col)))

# ha_bottom_cs = columnAnnotation(exp_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
#                                                       side = "bottom", labels_gp = gpar(fontsize = 14, col = gene_bin_sel$col)))

heat_col <- pals::ocean.balance(14)[9:14]
if(max(cs_mtx)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(cs_mtx)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(cs_mtx)-6))
}
names(col_vec) <- 0:max(cs_mtx)

ha_row_cs=rowAnnotation(df = anno_mtx, col = list(subclones=clst_col, superclones=sclst_col),
                        show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/", pro_name_d, "_complexHeatmap_tumor_superclone_",super_name1, "_vs_", super_name2, "_with_DE_dosage_annotations.pdf"), height = 4, width = 10)
Heatmap(as.matrix(cs_mtx_order), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        name = "csheatmap", show_row_names = T, row_names_side = "left",
        show_column_names = F, column_title = paste0(pro_name_d, "_consensusHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row_cs, bottom_annotation = ha_bottom_cs)
dev.off()


#----Plot gene freq with and without dosage effect---
gene_bin_sel_sample <- gene_bin_sel %>% dplyr::mutate(my_dos = plyr::mapvalues(col, from = c("#424451","#F47F20"), to = c("sCNA", "cCNA"))) 

df <- gene_bin_sel_sample
p1 <- ggplot(df, aes(x=my_dos, fill=my_dos)) + geom_bar() + theme_classic() + 
  scale_fill_manual(values=c("cCNA"="#F47F20", "sCNA"="#424451")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  labs(x="my_dos", y="gene count")
p1
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_barplot_gene_freq_dosage_effect_c7_vs_c11.pdf"), p1, width = 2.5, height = 2.5)



  
#---DE analysis of subclones 4 versus 7 ---#####
obj_rna <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
dna_meta_temp <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::select("rna_bc2","superclones","subclones","sample")

obj_rna2 <- obj_rna
obj_rna2@meta.data <- obj_rna2@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  left_join(dna_meta_temp, by = c("rowname" = "rna_bc2")) %>% column_to_rownames() 

clone_num1 <- 4
clone_num2 <- 7
super_clone1 <- paste0("c", clone_num1)
super_clone2 <- paste0("c", clone_num2)
super_name1 <- paste(super_clone1, collapse = "_")
super_name2 <- paste0(super_clone2, collapse = "_")
obj_rna2@meta.data <- obj_rna2@meta.data %>% mutate(temp = ifelse(subclones %in% super_clone1, "super_c1", 
                                                                  ifelse(subclones %in% super_clone2, "super_c2", "others")))

de.markers <- FindMarkers(obj_rna2, ident.1 = "super_c1", ident.2 = "super_c2", group.by = 'temp', assay = "RNA", 
                          only.pos = F, min.pct = 0, logfc.threshold = 0)
de.markers
###Volcano plot wth EnhanceVolcano###
library(EnhancedVolcano)
#--FC: 0.485 == 1.4X fold-change
p1<-EnhancedVolcano(de.markers, lab = rownames(de.markers), x="avg_log2FC", y="p_val_adj", FCcutoff = 0.485, pCutoff = 0.05, titleLabSize = 20,
                    # title = paste0(i, c("_DE_MAST")), 
                    labSize = 3.0, 
                    # pointSize = 2, col = c("black", "grey", "azure4", "red"),
                    xlim = c(-1.5, 1.5),
                    ylim = c(0, 13),
                    # drawConnectors = T, 
                    widthConnectors = 0.2, colConnectors = "grey30",typeConnectors = "open",
                    endsConnectors = "last", lengthConnectors = unit(0.01, 'npc'),
                    axisLabSize = 20, gridlines.major = FALSE, gridlines.minor = FALSE, border="full",
                    legendLabels = c('NS','Log2FC','adjP < 0.05','Log2FC and adjP'),
                    legendLabSize=15, legendPosition="bottom", colAlpha = 1, xlab = bquote(~Log[2]~"fold change"), ylab = bquote(~-log[10]~adjusted~italic(P)))
p1
p2 <- ggrastr::rasterize(p1, layers='Point', dpi=350)
cowplot::ggsave2(paste0("./figures/",pro_name_r, "_",super_name1, "_vs_", super_name2, "_volcano_DE.pdf"), p2, width = 4, height = 6)

#------heatmap with dosage effect status bar----#####
#---define copy number dosage effect: Due to the bad RNA-seq data quality, probably just define yes or no.
#---For yes: all subclones have same copy number, otherwise no.
cs_mtx <- t(varbin_mtx_tumor_log3@consensus)
my_order <- c(super_clone1, super_clone2)

cs_mtx <- cs_mtx %>% as.data.frame() %>% dplyr::filter(rownames(cs_mtx) %in% my_order)
cs_mtx_order <- cs_mtx[my_order,]
topleft(cs_mtx_order)

#----annotation bar--
anno_mtx <- rownames(cs_mtx_order) %>% as.data.frame() %>% dplyr::rename(subclones=".") %>% 
  mutate(superclones = ifelse(subclones %in% super_clone1, "superclone1", "superclone2")) %>% 
  dplyr::select(c("subclones", "superclones"))
rownames(anno_mtx) <- NULL

clst_col <- c(new_pal[clone_num1], new_pal[clone_num2])
names(clst_col) <- anno_mtx$subclones

sclst_col <- superclone_col[1:length(unique(anno_mtx$superclones))]
names(sclst_col) <- unique(anno_mtx$superclones) %>% sort()

#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")
#-----annotated genes---
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
de_genes_up <- de.markers %>% dplyr::filter(p_val_adj<0.05 & avg_log2FC > 0.485)
de_genes_down <- de.markers %>% dplyr::filter(p_val_adj<0.05 & avg_log2FC < -0.485)

dos_state <- apply(cs_mtx_order, 2, function(row) {ifelse(length(unique(row)) == 1, "cCNA", "sCNA")})
dos_state_df <- dos_state %>% as.data.frame() %>% dplyr::rename(dos_state = ".")

gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% rownames(rbind(de_genes_up, de_genes_down))) %>% 
  mutate(col=ifelse(dos_state[pos] == "sCNA", "#F47F20", "#424451"))

ha_bottom_cs = columnAnnotation(df = dos_state_df, col = list(dos_state=c("sCNA"="#F47F20", "cCNA"="#424451")), 
                                exp_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
                                                      side = "bottom", labels_gp = gpar(fontsize = 10, col = gene_bin_sel$col)))

# ha_bottom_cs = columnAnnotation(exp_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
#                                                       side = "bottom", labels_gp = gpar(fontsize = 14, col = gene_bin_sel$col)))

heat_col <- pals::ocean.balance(14)[9:14]
if(max(cs_mtx)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(cs_mtx)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(cs_mtx)-6))
}
names(col_vec) <- 0:max(cs_mtx)

ha_row_cs=rowAnnotation(df = anno_mtx, col = list(subclones=clst_col, superclones=sclst_col),
                        show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/", pro_name_d, "_complexHeatmap_tumor_superclone_",super_name1, "_vs_", super_name2, "_with_DE_dosage_annotations.pdf"), height = 4, width = 10)
Heatmap(as.matrix(cs_mtx_order), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        name = "csheatmap", show_row_names = T, row_names_side = "left",
        show_column_names = F, column_title = paste0(pro_name_d, "_consensusHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row_cs, bottom_annotation = ha_bottom_cs)
dev.off()


#----Plot gene freq with and without dosage effect---
gene_bin_sel_sample <- gene_bin_sel %>% dplyr::mutate(my_dos = plyr::mapvalues(col, from = c("#424451","#F47F20"), to = c("sCNA", "cCNA"))) 

df <- gene_bin_sel_sample
p1 <- ggplot(df, aes(x=my_dos, fill=my_dos)) + geom_bar() + theme_classic() + 
  scale_fill_manual(values=c("cCNA"="#F47F20", "sCNA"="#424451")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  labs(x="my_dos", y="gene count")
p1
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_barplot_gene_freq_dosage_effect_c4_vs_c7.pdf"), p1, width = 2.5, height = 2.5)




#---Make DNA-RNA co-line plots (c4 vs c7)---#########
obj_rna <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
dna_meta_temp <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::select("rna_bc2","superclones","subclones","sample")

obj_rna2 <- obj_rna
obj_rna2@meta.data <- obj_rna2@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  left_join(dna_meta_temp, by = c("rowname" = "rna_bc2")) %>% column_to_rownames() 

clone_num1 <- 4
clone_num2 <- 7
super_clone1 <- paste0("c", clone_num1)
super_clone2 <- paste0("c", clone_num2)
super_name1 <- paste(super_clone1, collapse = "_")
super_name2 <- paste0(super_clone2, collapse = "_")
obj_rna2@meta.data <- obj_rna2@meta.data %>% mutate(temp = ifelse(subclones %in% super_clone1, "super_c1", 
                                                                  ifelse(subclones %in% super_clone2, "super_c2", "others")))

varbin_mtx_tumor_log3 <- readRDS(paste0("objects/", pro_name_d, c("_tumor_integerCN_copykit.rds")))
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
bin_chrpos <- varbin_mtx_tumor_log3@rowRanges %>% as.data.frame() %>% rownames_to_column() %>% 
  mutate(chr_pos = paste0("chr",seqnames,":",start, "-", end)) %>% dplyr::select(rowname, chr_pos)

clone1_cells <- obj_rna2@meta.data %>% dplyr::filter(temp == "super_c1") %>% rownames()
# rna_exp_clone1_cells <- obj_rna2@assays$RNA@scale.data[,clone1_cells]
rna_exp_clone1_cells <- obj_rna2@assays$RNA@data[,clone1_cells]
clone2_cells <- obj_rna2@meta.data %>% dplyr::filter(temp == "super_c2") %>% rownames()
# rna_exp_clone2_cells <- obj_rna2@assays$RNA@scale.data[,clone2_cells]
rna_exp_clone2_cells <- obj_rna2@assays$RNA@data[,clone2_cells]
# obj_clone1_2 <- subset(obj_rna2, cells = c(clone1_cells, clone2_cells))
# all.genes <- rownames(obj_clone1_2)
# obj_clone1_2_scaled <- ScaleData(obj_clone1_2, features = all.genes)
# rna_exp_clone1_cells <- obj_clone1_2_scaled@assays$RNA@scale.data[,clone1_cells]
# rna_exp_clone2_cells <- obj_clone1_2_scaled@assays$RNA@scale.data[,clone2_cells]


rna_cons_clone1 <- apply(rna_exp_clone1_cells, 1, mean) %>% as.data.frame() %>% `colnames<-`("med_scaled_exp_clone1") %>% 
  rownames_to_column() %>% dplyr::filter(med_scaled_exp_clone1 > 0)
rna_cons_clone2 <- apply(rna_exp_clone2_cells, 1, mean) %>% as.data.frame() %>% `colnames<-`("med_scaled_exp_clone2") %>% 
  rownames_to_column() %>% dplyr::filter(med_scaled_exp_clone2 > 0)

clone1_2_con_exp <- inner_join(rna_cons_clone1, rna_cons_clone2, by = "rowname")

dna_con <- varbin_mtx_tumor_log3@consensus %>% rownames_to_column() %>% dplyr::select(rowname,all_of(c(super_clone1, super_clone2)))

varbin_mtx_tumor_log3@consensus

dna_rna_merged_cna_exp <- gene_bin %>% mutate(pos = as.character(pos)) %>% left_join(dna_con, ., by = c("rowname" = "pos")) %>% 
  left_join(clone1_2_con_exp, by = c("gene" = "rowname")) %>% left_join(bin_chrpos, by = "rowname")

# dna_rna_merged_cna_exp %>% dplyr::filter(gene == "VAMP8")
# dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% dplyr::filter(rowname > 1339 & rowname < 1349)
# dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% dplyr::filter(rowname > 1200 & rowname < 1500)
dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% dplyr::filter(stringr::str_detect(chr_pos, "chr3:"))



plot_s <- list()
j<-1
#---checked chr20-chrX, due to the chr size and gene expression number, it's can't do the smoonth. 
for (i in paste0("chr",c(1:16),":")) {
  # i <- "chr16"
  dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% dplyr::filter(stringr::str_detect(chr_pos, i))
  dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp_target %>% mutate(med_scaled_exp_clone1 = ifelse(med_scaled_exp_clone1>0.05, med_scaled_exp_clone1, NA), 
                                                                            med_scaled_exp_clone2 = ifelse(med_scaled_exp_clone2>0.05, med_scaled_exp_clone2, NA))
  dna_rna_merged_cna_exp_target_2 <- dna_rna_merged_cna_exp_target %>% mutate(diff = med_scaled_exp_clone1 - med_scaled_exp_clone2)
  
  dna_rna_merged_cna_exp_target_3 <- dna_rna_merged_cna_exp_target_2 %>% mutate(
    normal_clone1 = case_when(
      !is.na(med_scaled_exp_clone1) & !is.na(med_scaled_exp_clone2) ~ med_scaled_exp_clone1/med_scaled_exp_clone2,
      !is.na(med_scaled_exp_clone1) & is.na(med_scaled_exp_clone2) ~ med_scaled_exp_clone1/0.05,
      is.na(med_scaled_exp_clone1) & !is.na(med_scaled_exp_clone2) ~ 0.05/med_scaled_exp_clone2,
      TRUE ~ NA_real_  # Returns NA for numeric columns
    ),
    normal_clone2 = ifelse(is.na(normal_clone1), NA, 1), 
    cna_diff = get(super_clone1) - get(super_clone2),
    cna_base = 0
  )
  
  dna_vamp8_mtx <- dna_rna_merged_cna_exp_target_3 %>% dplyr::select(rowname, diff,cna_diff)
  dna_vamp8_mtx$rowname_numeric <- as.numeric(as.character(dna_vamp8_mtx$rowname))
  transformation_factor <- ifelse(max(dna_vamp8_mtx$cna_diff) >0, 
                                  2*max(dna_vamp8_mtx$cna_diff, na.rm = TRUE) / max(dna_vamp8_mtx$diff, na.rm = TRUE), 
                                  2*(-min(dna_vamp8_mtx$cna_diff, na.rm = TRUE) / max(dna_vamp8_mtx$diff, na.rm = TRUE)))
  
  p1 <- ggplot(dna_vamp8_mtx, aes(x = rowname_numeric)) +
    geom_point(aes(y = diff)) + geom_smooth(aes(y = diff), method = "gam", na.rm = T, se = F, color = "red") + 
    geom_line(aes(y = cna_diff / transformation_factor), color = "blue") +
    scale_y_continuous("Diff of log(gene expression)",
                       sec.axis = sec_axis(~ . * transformation_factor, name = "CNA Diff")) +
    labs(x = "Genomic bins", title = i) + theme_cowplot()
  
  p1
  plot_s[[j]]<-p1
  j<-j+1
}

p_list<-plot_grid(plotlist = plot_s, ncol = 5)

cowplot::ggsave2(paste0("./figures/", pro_name, "_dna_rna_line_plot_",super_clone1, "_vs_",super_clone2, "_chr1_18.pdf"), p_list, width = 25, height = 16)



#---Make DNA-RNA co-line plots (c1 vs c7)---#########
obj_rna <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
dna_meta_temp <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::select("rna_bc2","superclones","subclones","sample")

obj_rna2 <- obj_rna
obj_rna2@meta.data <- obj_rna2@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  left_join(dna_meta_temp, by = c("rowname" = "rna_bc2")) %>% column_to_rownames() 

clone_num1 <- 1
clone_num2 <- 7
super_clone1 <- paste0("c", clone_num1)
super_clone2 <- paste0("c", clone_num2)
super_name1 <- paste(super_clone1, collapse = "_")
super_name2 <- paste0(super_clone2, collapse = "_")
obj_rna2@meta.data <- obj_rna2@meta.data %>% mutate(temp = ifelse(subclones %in% super_clone1, "super_c1", 
                                                                  ifelse(subclones %in% super_clone2, "super_c2", "others")))

varbin_mtx_tumor_log3 <- readRDS(paste0("objects/", pro_name_d, c("_tumor_integerCN_copykit.rds")))
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
bin_chrpos <- varbin_mtx_tumor_log3@rowRanges %>% as.data.frame() %>% rownames_to_column() %>% 
  mutate(chr_pos = paste0("chr",seqnames,":",start, "-", end)) %>% dplyr::select(rowname, chr_pos)

clone1_cells <- obj_rna2@meta.data %>% dplyr::filter(temp == "super_c1") %>% rownames()
# rna_exp_clone1_cells <- obj_rna2@assays$RNA@scale.data[,clone1_cells]
rna_exp_clone1_cells <- obj_rna2@assays$RNA@data[,clone1_cells]
clone2_cells <- obj_rna2@meta.data %>% dplyr::filter(temp == "super_c2") %>% rownames()
# rna_exp_clone2_cells <- obj_rna2@assays$RNA@scale.data[,clone2_cells]
rna_exp_clone2_cells <- obj_rna2@assays$RNA@data[,clone2_cells]
# obj_clone1_2 <- subset(obj_rna2, cells = c(clone1_cells, clone2_cells))
# all.genes <- rownames(obj_clone1_2)
# obj_clone1_2_scaled <- ScaleData(obj_clone1_2, features = all.genes)
# rna_exp_clone1_cells <- obj_clone1_2_scaled@assays$RNA@scale.data[,clone1_cells]
# rna_exp_clone2_cells <- obj_clone1_2_scaled@assays$RNA@scale.data[,clone2_cells]


rna_cons_clone1 <- apply(rna_exp_clone1_cells, 1, mean) %>% as.data.frame() %>% `colnames<-`("med_scaled_exp_clone1") %>% 
  rownames_to_column() %>% dplyr::filter(med_scaled_exp_clone1 > 0)
rna_cons_clone2 <- apply(rna_exp_clone2_cells, 1, mean) %>% as.data.frame() %>% `colnames<-`("med_scaled_exp_clone2") %>% 
  rownames_to_column() %>% dplyr::filter(med_scaled_exp_clone2 > 0)

clone1_2_con_exp <- inner_join(rna_cons_clone1, rna_cons_clone2, by = "rowname")

dna_con <- varbin_mtx_tumor_log3@consensus %>% rownames_to_column() %>% dplyr::select(rowname,all_of(c(super_clone1, super_clone2)))

varbin_mtx_tumor_log3@consensus

dna_rna_merged_cna_exp <- gene_bin %>% mutate(pos = as.character(pos)) %>% left_join(dna_con, ., by = c("rowname" = "pos")) %>% 
  left_join(clone1_2_con_exp, by = c("gene" = "rowname")) %>% left_join(bin_chrpos, by = "rowname")

# dna_rna_merged_cna_exp %>% dplyr::filter(gene == "VAMP8")
# dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% dplyr::filter(rowname > 1339 & rowname < 1349)
# dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% dplyr::filter(rowname > 1200 & rowname < 1500)
# dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% dplyr::filter(stringr::str_detect(chr_pos, "chr3:"))



plot_s <- list()
j<-1
#---checked chr20-chrX, due to the chr size and gene expression number, it's can't do the smoonth. 
for (i in paste0("chr",c(1:19),":")) {
  # i <- "chr1:"
  dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% dplyr::filter(stringr::str_detect(chr_pos, i))
  dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp_target %>% mutate(med_scaled_exp_clone1 = ifelse(med_scaled_exp_clone1>0.05, med_scaled_exp_clone1, NA), 
                                                                            med_scaled_exp_clone2 = ifelse(med_scaled_exp_clone2>0.05, med_scaled_exp_clone2, NA))
  dna_rna_merged_cna_exp_target_2 <- dna_rna_merged_cna_exp_target %>% mutate(diff = med_scaled_exp_clone1 - med_scaled_exp_clone2)
  
  dna_rna_merged_cna_exp_target_3 <- dna_rna_merged_cna_exp_target_2 %>% mutate(
    normal_clone1 = case_when(
      !is.na(med_scaled_exp_clone1) & !is.na(med_scaled_exp_clone2) ~ med_scaled_exp_clone1/med_scaled_exp_clone2,
      !is.na(med_scaled_exp_clone1) & is.na(med_scaled_exp_clone2) ~ med_scaled_exp_clone1/0.05,
      is.na(med_scaled_exp_clone1) & !is.na(med_scaled_exp_clone2) ~ 0.05/med_scaled_exp_clone2,
      TRUE ~ NA_real_  # Returns NA for numeric columns
    ),
    normal_clone2 = ifelse(is.na(normal_clone1), NA, 1), 
    cna_diff = get(super_clone1) - get(super_clone2),
    cna_base = 0
  )
  
  dna_vamp8_mtx <- dna_rna_merged_cna_exp_target_3 %>% dplyr::select(rowname, diff,cna_diff)
  dna_vamp8_mtx$rowname_numeric <- as.numeric(as.character(dna_vamp8_mtx$rowname))
  transformation_factor <- ifelse(max(dna_vamp8_mtx$cna_diff) >0, 
                                  2*max(dna_vamp8_mtx$cna_diff, na.rm = TRUE) / max(dna_vamp8_mtx$diff, na.rm = TRUE), 
                                  2*(-min(dna_vamp8_mtx$cna_diff, na.rm = TRUE) / max(dna_vamp8_mtx$diff, na.rm = TRUE)))
  
  p1 <- ggplot(dna_vamp8_mtx, aes(x = rowname_numeric)) +
    geom_point(aes(y = diff)) + geom_smooth(aes(y = diff), method = "gam", na.rm = T, se = F, color = "red") + 
    geom_line(aes(y = cna_diff / transformation_factor), color = "blue") +
    scale_y_continuous("Diff of log(gene expression)",
                       sec.axis = sec_axis(~ . * transformation_factor, name = "CNA Diff")) +
    labs(x = "Genomic bins", title = i) + theme_cowplot()
  
  p1
  plot_s[[j]]<-p1
  j<-j+1
}

p_list<-plot_grid(plotlist = plot_s, ncol = 5)

cowplot::ggsave2(paste0("./figures/", pro_name, "_dna_rna_line_plot_chr1_19.pdf"), p_list, width = 25, height = 16)

#----staticstic counts of DNA RNA matching---
dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% 
  mutate(med_scaled_exp_clone1 = ifelse(med_scaled_exp_clone1>0.05, med_scaled_exp_clone1, NA), 
         med_scaled_exp_clone2 = ifelse(med_scaled_exp_clone2>0.05, med_scaled_exp_clone2, NA)) %>% 
  mutate(diff = med_scaled_exp_clone1 - med_scaled_exp_clone2) %>% mutate(
    normal_clone1 = case_when(
      !is.na(med_scaled_exp_clone1) & !is.na(med_scaled_exp_clone2) ~ med_scaled_exp_clone1/med_scaled_exp_clone2,
      !is.na(med_scaled_exp_clone1) & is.na(med_scaled_exp_clone2) ~ med_scaled_exp_clone1/0.05,
      is.na(med_scaled_exp_clone1) & !is.na(med_scaled_exp_clone2) ~ 0.05/med_scaled_exp_clone2,
      TRUE ~ NA_real_  # Returns NA for numeric columns
    ),
    normal_clone2 = ifelse(is.na(normal_clone1), NA, 1), 
    cna_diff = get(super_clone1) - get(super_clone2),
    cna_base = 0,
    chr = stringr::str_extract(chr_pos, "chr[^:]+")
  ) %>% dplyr::select(rowname, chr, diff,cna_diff) %>% na.omit(cols="diff")


data <- dna_rna_merged_cna_exp_target %>% 
  dplyr::group_by(chr) %>% 
  dplyr::mutate(group = cumsum(cna_diff != lag(cna_diff, default = dplyr::first(cna_diff)))) %>%
  ungroup()

# Grouping by 'chr' and 'group' and then calculating the mean of 'diff'
grouped_median <- data %>%
  dplyr::group_by(chr, group) %>%
  dplyr::summarise(median_diff = median(diff), first_cna_diff = dplyr::first(cna_diff))

grouped_median %>% filter(first_cna_diff != 0)
grouped_median %>% filter(first_cna_diff == 0 & abs(median_diff) > 0.02) %>% view()

# p2 <- ggplot(rna_vamp8_mtx_2, aes(x = as.factor(rowname), y = diff)) +
#   geom_point() + theme_bw() +
#   labs(x = "genomic bins", y = "Diff gene exp")
# p2
# 
# 
# 
# dna_vamp8_mtx <- dna_rna_merged_cna_exp_target_3 %>% dplyr::select(rowname, cna_diff, cna_c7, chr_pos)
# dna_vamp8_mtx_long <- dna_vamp8_mtx %>% pivot_longer(cols = c(cna_diff, cna_c7),
#                                                      names_to = "variable", values_to = "value")
# 
# rna_vamp8_mtx <- dna_rna_merged_cna_exp_target_3 %>% dplyr::select(rowname, normal_clone1, normal_clone2, chr_pos)
# rna_vamp8_mtx_long <- rna_vamp8_mtx %>% pivot_longer(cols = c(normal_clone1, normal_clone2),
#                                                      names_to = "variable", values_to = "value")
# 
# 
# #---Another way is DNA one line, RNA one line, both in the same plot. 
# 
# # rna_vamp8_mtx_2 <- dna_rna_merged_cna_exp_target_2 %>% dplyr::select(rowname, diff, chr_pos)
# # rna_vamp8_mtx_long <- rna_vamp8_mtx_2 %>% pivot_longer(cols = c(med_scaled_exp_clone1, med_scaled_exp_clone2),
# #                                                      names_to = "variable", values_to = "value")
# # 
# # p2 <- ggplot(rna_vamp8_mtx_2, aes(x = as.factor(rowname), y = diff)) +
# #   geom_point() + theme_bw() +
# #   labs(x = "genomic bins", y = "Diff gene exp") 
# # p2
# 
# #=============++
# dna_vamp8_mtx <- dna_rna_merged_cna_exp_target %>% dplyr::select(rowname, all_of(super_clone1), all_of(super_clone2), chr_pos)
# dna_vamp8_mtx_long <- dna_vamp8_mtx %>% pivot_longer(cols = c(super_clone1, super_clone2),
#                names_to = "variable", values_to = "value")
# # rna_vamp8_mtx <- dna_rna_merged_cna_exp_target %>% dplyr::select(rowname, med_scaled_exp_clone1, med_scaled_exp_clone2, chr_pos)
# # rna_vamp8_mtx_long <- rna_vamp8_mtx %>% pivot_longer(cols = c(med_scaled_exp_clone1, med_scaled_exp_clone2),
# #                                                      names_to = "variable", values_to = "value")
# #==============+++
# p1 <- ggplot(dna_vamp8_mtx_long, aes(x = as.factor(rowname), y = value, color = variable, group = variable)) +
#   geom_line() + theme_bw() +
#   labs(y = "DNA copy number") +
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 1)) +
#   scale_color_manual(values = new_pal2[c(super_clone1, super_clone2)]) +
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank())
# 
# p1
# p2 <- ggplot(rna_vamp8_mtx_long, aes(x = as.factor(rowname), y = value, color = variable, group = variable)) +
#   geom_point() + geom_smooth(method = "loess", na.rm = T, se = F) + 
#   theme_bw() +
#   labs(x = "genomic bins", y = "normalized gene exp") +
#   scale_color_manual(values = c("#5C88DAB2","#00AF66B2"))
# p2
# p3 <- plot_grid(plotlist=list(p1,p2), ncol=1, align='v')
# p3
# cowplot::ggsave2(paste0("./figures/", pro_name, "_dna_rna_VAMP8.pdf"), p3, width = 5, height = 4)
# 
# 
# dna_rna_merged_cna_exp %>% dplyr::filter(gene == "PTS")
# dna_rna_merged_cna_exp %>% dplyr::filter(gene == "KIRREL3")
# 
# dna_rna_merged_cna_exp_PTS <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% dplyr::filter(rowname > 7900 & rowname < 8100)
# 
# dna_PTS_mtx <- dna_rna_merged_cna_exp_PTS %>% dplyr::select(rowname, c2, c7, chr_pos)
# dna_PTS_mtx_long <- dna_PTS_mtx %>% pivot_longer(cols = c(c2, c7),
#                                                      names_to = "variable", values_to = "value")
# rna_PTS_mtx <- dna_rna_merged_cna_exp_PTS %>% dplyr::select(rowname, med_scaled_exp_clone1, med_scaled_exp_clone2, chr_pos)
# rna_PTS_mtx_long <- rna_PTS_mtx %>% pivot_longer(cols = c(med_scaled_exp_clone1, med_scaled_exp_clone2),
#                                                      names_to = "variable", values_to = "value")
# 
# p1 <- ggplot(dna_PTS_mtx_long, aes(x = as.factor(rowname), y = value, color = variable, group = variable)) +
#   geom_line() + theme_bw() +
#   labs(y = "DNA copy number") +
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 1)) +
#   scale_color_manual(values = new_pal2[c(super_clone1, super_clone2)]) +
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank())
# 
# p2 <- ggplot(rna_PTS_mtx_long, aes(x = as.factor(rowname), y = value, color = variable, group = variable)) +
#   geom_point() + theme_bw() +
#   labs(x = "genomic bins", y = "Scaled gene exp") +
#   scale_color_manual(values = c("#5C88DAB2","#00AF66B2"))
# 
# p3 <- plot_grid(plotlist=list(p1,p2), ncol=1, align='v')
# p3
# cowplot::ggsave2(paste0("./figures/", pro_name, "_dna_rna_PTS.pdf"), p3, width = 8, height = 4)




#---Make DNA-RNA co-line plots (c2 vs c4)---#########
obj_rna <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
dna_meta_temp <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::select("rna_bc2","superclones","subclones","sample")

obj_rna2 <- obj_rna
obj_rna2@meta.data <- obj_rna2@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  left_join(dna_meta_temp, by = c("rowname" = "rna_bc2")) %>% column_to_rownames() 

clone_num1 <- 2
clone_num2 <- 4
super_clone1 <- paste0("c", clone_num1)
super_clone2 <- paste0("c", clone_num2)
super_name1 <- paste(super_clone1, collapse = "_")
super_name2 <- paste0(super_clone2, collapse = "_")
obj_rna2@meta.data <- obj_rna2@meta.data %>% mutate(temp = ifelse(subclones %in% super_clone1, "super_c1", 
                                                                  ifelse(subclones %in% super_clone2, "super_c2", "others")))

varbin_mtx_tumor_log3 <- readRDS(paste0("objects/", pro_name_d, c("_tumor_integerCN_copykit.rds")))
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
bin_chrpos <- varbin_mtx_tumor_log3@rowRanges %>% as.data.frame() %>% rownames_to_column() %>% 
  mutate(chr_pos = paste0("chr",seqnames,":",start, "-", end)) %>% dplyr::select(rowname, chr_pos)

clone1_cells <- obj_rna2@meta.data %>% dplyr::filter(temp == "super_c1") %>% rownames()
# rna_exp_clone1_cells <- obj_rna2@assays$RNA@scale.data[,clone1_cells]
rna_exp_clone1_cells <- obj_rna2@assays$RNA@data[,clone1_cells]
clone2_cells <- obj_rna2@meta.data %>% dplyr::filter(temp == "super_c2") %>% rownames()
# rna_exp_clone2_cells <- obj_rna2@assays$RNA@scale.data[,clone2_cells]
rna_exp_clone2_cells <- obj_rna2@assays$RNA@data[,clone2_cells]
# obj_clone1_2 <- subset(obj_rna2, cells = c(clone1_cells, clone2_cells))
# all.genes <- rownames(obj_clone1_2)
# obj_clone1_2_scaled <- ScaleData(obj_clone1_2, features = all.genes)
# rna_exp_clone1_cells <- obj_clone1_2_scaled@assays$RNA@scale.data[,clone1_cells]
# rna_exp_clone2_cells <- obj_clone1_2_scaled@assays$RNA@scale.data[,clone2_cells]


rna_cons_clone1 <- apply(rna_exp_clone1_cells, 1, mean) %>% as.data.frame() %>% `colnames<-`("med_scaled_exp_clone1") %>% 
  rownames_to_column() %>% dplyr::filter(med_scaled_exp_clone1 > 0)
rna_cons_clone2 <- apply(rna_exp_clone2_cells, 1, mean) %>% as.data.frame() %>% `colnames<-`("med_scaled_exp_clone2") %>% 
  rownames_to_column() %>% dplyr::filter(med_scaled_exp_clone2 > 0)

clone1_2_con_exp <- inner_join(rna_cons_clone1, rna_cons_clone2, by = "rowname")
dna_con <- varbin_mtx_tumor_log3@consensus %>% rownames_to_column() %>% dplyr::select(rowname,all_of(c(super_clone1, super_clone2)))
dna_rna_merged_cna_exp <- gene_bin %>% mutate(pos = as.character(pos)) %>% left_join(dna_con, ., by = c("rowname" = "pos")) %>% 
  left_join(clone1_2_con_exp, by = c("gene" = "rowname")) %>% left_join(bin_chrpos, by = "rowname")


plot_s <- list()
j<-1
#---checked chr20-chrX, due to the chr size and gene expression number, it's can't do the smoonth. 
for (i in paste0("chr",c(1:15),":")) {
  #i <- "chr15:"
  dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% dplyr::filter(stringr::str_detect(chr_pos, i))
  dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp_target %>% mutate(med_scaled_exp_clone1 = ifelse(med_scaled_exp_clone1>0.05, med_scaled_exp_clone1, NA), 
                                                                            med_scaled_exp_clone2 = ifelse(med_scaled_exp_clone2>0.05, med_scaled_exp_clone2, NA))
  dna_rna_merged_cna_exp_target_2 <- dna_rna_merged_cna_exp_target %>% mutate(diff = med_scaled_exp_clone1 - med_scaled_exp_clone2)
  
  dna_rna_merged_cna_exp_target_3 <- dna_rna_merged_cna_exp_target_2 %>% mutate(
    normal_clone1 = case_when(
      !is.na(med_scaled_exp_clone1) & !is.na(med_scaled_exp_clone2) ~ med_scaled_exp_clone1/med_scaled_exp_clone2,
      !is.na(med_scaled_exp_clone1) & is.na(med_scaled_exp_clone2) ~ med_scaled_exp_clone1/0.05,
      is.na(med_scaled_exp_clone1) & !is.na(med_scaled_exp_clone2) ~ 0.05/med_scaled_exp_clone2,
      TRUE ~ NA_real_  # Returns NA for numeric columns
    ),
    normal_clone2 = ifelse(is.na(normal_clone1), NA, 1), 
    cna_diff = get(super_clone1) - get(super_clone2),
    cna_base = 0
  )
  
  dna_vamp8_mtx <- dna_rna_merged_cna_exp_target_3 %>% dplyr::select(rowname, diff,cna_diff)
  dna_vamp8_mtx$rowname_numeric <- as.numeric(as.character(dna_vamp8_mtx$rowname))
  transformation_factor <- ifelse(max(dna_vamp8_mtx$cna_diff) >0, 
                                  2*max(dna_vamp8_mtx$cna_diff, na.rm = TRUE) / max(dna_vamp8_mtx$diff, na.rm = TRUE), 
                                  2*(-min(dna_vamp8_mtx$cna_diff, na.rm = TRUE) / max(dna_vamp8_mtx$diff, na.rm = TRUE)))
  
  p1 <- ggplot(dna_vamp8_mtx, aes(x = rowname_numeric)) +
    geom_point(aes(y = diff)) + geom_smooth(aes(y = diff), method = "gam", na.rm = T, se = F, color = "red") + 
    geom_line(aes(y = cna_diff / transformation_factor), color = "blue") +
    scale_y_continuous("Diff of log(gene expression)",
                       sec.axis = sec_axis(~ . * transformation_factor, name = "CNA Diff")) +
    labs(x = "Genomic bins", title = i) + theme_cowplot()
  
  p1
  plot_s[[j]]<-p1
  j<-j+1
}

p_list<-plot_grid(plotlist = plot_s, ncol = 5)

cowplot::ggsave2(paste0("./figures/", pro_name, "_dna_rna_line_plot_",super_clone1, "_vs_",super_clone2, "_chr1_16.pdf"), p_list, width = 25, height = 12)

#----staticstic counts of DNA RNA matching---
dna_rna_merged_cna_exp_target <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% 
  mutate(med_scaled_exp_clone1 = ifelse(med_scaled_exp_clone1>0.05, med_scaled_exp_clone1, NA), 
         med_scaled_exp_clone2 = ifelse(med_scaled_exp_clone2>0.05, med_scaled_exp_clone2, NA)) %>% 
  mutate(diff = med_scaled_exp_clone1 - med_scaled_exp_clone2) %>% mutate(
    normal_clone1 = case_when(
      !is.na(med_scaled_exp_clone1) & !is.na(med_scaled_exp_clone2) ~ med_scaled_exp_clone1/med_scaled_exp_clone2,
      !is.na(med_scaled_exp_clone1) & is.na(med_scaled_exp_clone2) ~ med_scaled_exp_clone1/0.05,
      is.na(med_scaled_exp_clone1) & !is.na(med_scaled_exp_clone2) ~ 0.05/med_scaled_exp_clone2,
      TRUE ~ NA_real_  # Returns NA for numeric columns
    ),
    normal_clone2 = ifelse(is.na(normal_clone1), NA, 1), 
    cna_diff = get(super_clone1) - get(super_clone2),
    cna_base = 0,
    chr = stringr::str_extract(chr_pos, "chr[^:]+")
  ) %>% dplyr::select(rowname, chr, diff,cna_diff) %>% na.omit(cols="diff")


data <- dna_rna_merged_cna_exp_target %>% 
  dplyr::group_by(chr) %>% 
  mutate(group = cumsum(cna_diff != lag(cna_diff, default = first(cna_diff)))) %>%
  ungroup()

# Grouping by 'chr' and 'group' and then calculating the mean of 'diff'
grouped_median <- data %>%
  dplyr::group_by(chr, group) %>%
  dplyr::summarise(median_diff = median(diff), first_cna_diff = first(cna_diff))

grouped_median %>% filter(first_cna_diff != 0)
grouped_median %>% filter(first_cna_diff == 0 & abs(median_diff) > 0.02) %>% view()

# p2 <- ggplot(rna_vamp8_mtx_2, aes(x = as.factor(rowname), y = diff)) +
#   geom_point() + theme_bw() +
#   labs(x = "genomic bins", y = "Diff gene exp")
# p2
# 
# 
# 
# dna_vamp8_mtx <- dna_rna_merged_cna_exp_target_3 %>% dplyr::select(rowname, cna_diff, cna_c7, chr_pos)
# dna_vamp8_mtx_long <- dna_vamp8_mtx %>% pivot_longer(cols = c(cna_diff, cna_c7),
#                                                      names_to = "variable", values_to = "value")
# 
# rna_vamp8_mtx <- dna_rna_merged_cna_exp_target_3 %>% dplyr::select(rowname, normal_clone1, normal_clone2, chr_pos)
# rna_vamp8_mtx_long <- rna_vamp8_mtx %>% pivot_longer(cols = c(normal_clone1, normal_clone2),
#                                                      names_to = "variable", values_to = "value")
# 
# 
# #---Another way is DNA one line, RNA one line, both in the same plot. 
# 
# # rna_vamp8_mtx_2 <- dna_rna_merged_cna_exp_target_2 %>% dplyr::select(rowname, diff, chr_pos)
# # rna_vamp8_mtx_long <- rna_vamp8_mtx_2 %>% pivot_longer(cols = c(med_scaled_exp_clone1, med_scaled_exp_clone2),
# #                                                      names_to = "variable", values_to = "value")
# # 
# # p2 <- ggplot(rna_vamp8_mtx_2, aes(x = as.factor(rowname), y = diff)) +
# #   geom_point() + theme_bw() +
# #   labs(x = "genomic bins", y = "Diff gene exp") 
# # p2
# 
# #=============++
# dna_vamp8_mtx <- dna_rna_merged_cna_exp_target %>% dplyr::select(rowname, all_of(super_clone1), all_of(super_clone2), chr_pos)
# dna_vamp8_mtx_long <- dna_vamp8_mtx %>% pivot_longer(cols = c(super_clone1, super_clone2),
#                names_to = "variable", values_to = "value")
# # rna_vamp8_mtx <- dna_rna_merged_cna_exp_target %>% dplyr::select(rowname, med_scaled_exp_clone1, med_scaled_exp_clone2, chr_pos)
# # rna_vamp8_mtx_long <- rna_vamp8_mtx %>% pivot_longer(cols = c(med_scaled_exp_clone1, med_scaled_exp_clone2),
# #                                                      names_to = "variable", values_to = "value")
# #==============+++
# p1 <- ggplot(dna_vamp8_mtx_long, aes(x = as.factor(rowname), y = value, color = variable, group = variable)) +
#   geom_line() + theme_bw() +
#   labs(y = "DNA copy number") +
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 1)) +
#   scale_color_manual(values = new_pal2[c(super_clone1, super_clone2)]) +
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank())
# 
# p1
# p2 <- ggplot(rna_vamp8_mtx_long, aes(x = as.factor(rowname), y = value, color = variable, group = variable)) +
#   geom_point() + geom_smooth(method = "loess", na.rm = T, se = F) + 
#   theme_bw() +
#   labs(x = "genomic bins", y = "normalized gene exp") +
#   scale_color_manual(values = c("#5C88DAB2","#00AF66B2"))
# p2
# p3 <- plot_grid(plotlist=list(p1,p2), ncol=1, align='v')
# p3
# cowplot::ggsave2(paste0("./figures/", pro_name, "_dna_rna_VAMP8.pdf"), p3, width = 5, height = 4)
# 
# 
# dna_rna_merged_cna_exp %>% dplyr::filter(gene == "PTS")
# dna_rna_merged_cna_exp %>% dplyr::filter(gene == "KIRREL3")
# 
# dna_rna_merged_cna_exp_PTS <- dna_rna_merged_cna_exp %>% mutate(rowname = as.integer(rowname)) %>% dplyr::filter(rowname > 7900 & rowname < 8100)
# 
# dna_PTS_mtx <- dna_rna_merged_cna_exp_PTS %>% dplyr::select(rowname, c2, c7, chr_pos)
# dna_PTS_mtx_long <- dna_PTS_mtx %>% pivot_longer(cols = c(c2, c7),
#                                                      names_to = "variable", values_to = "value")
# rna_PTS_mtx <- dna_rna_merged_cna_exp_PTS %>% dplyr::select(rowname, med_scaled_exp_clone1, med_scaled_exp_clone2, chr_pos)
# rna_PTS_mtx_long <- rna_PTS_mtx %>% pivot_longer(cols = c(med_scaled_exp_clone1, med_scaled_exp_clone2),
#                                                      names_to = "variable", values_to = "value")
# 
# p1 <- ggplot(dna_PTS_mtx_long, aes(x = as.factor(rowname), y = value, color = variable, group = variable)) +
#   geom_line() + theme_bw() +
#   labs(y = "DNA copy number") +
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 1)) +
#   scale_color_manual(values = new_pal2[c(super_clone1, super_clone2)]) +
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank())
# 
# p2 <- ggplot(rna_PTS_mtx_long, aes(x = as.factor(rowname), y = value, color = variable, group = variable)) +
#   geom_point() + theme_bw() +
#   labs(x = "genomic bins", y = "Scaled gene exp") +
#   scale_color_manual(values = c("#5C88DAB2","#00AF66B2"))
# 
# p3 <- plot_grid(plotlist=list(p1,p2), ncol=1, align='v')
# p3
# cowplot::ggsave2(paste0("./figures/", pro_name, "_dna_rna_PTS.pdf"), p3, width = 8, height = 4)




#----clonal status---testing------#######
for (i in c(0.25,0.3,0.35)) {

low_cutoff <- -i
up_cutoff <-  i
cell_cutoff <- 0.85
neu_cutoff <- 0.90
my_sample <- clonality_log_trinary_neu(log_ratio_df=ht_mtx, lower_cutoff = low_cutoff, upper_cutoff = up_cutoff, 
                                       cell_pct = cell_cutoff, neu_pct = neu_cutoff)
table(my_sample)
# cs <- as.data.frame(my_sample)
sig_de.markers <- read.csv(paste0("./metrics/", pro_name_r, "_volcano_DE_sig.csv"), row.names = 1)
gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% sig_de.markers$rowname)

gene_bin_sel_sample <- as.data.frame(my_sample) %>% rownames_to_column() %>% dplyr::mutate(rowname = as.numeric(rowname)) %>% 
  dplyr::right_join(gene_bin_sel, by = c("rowname" = "pos")) %>% 
  dplyr::mutate(my_color = plyr::mapvalues(my_sample, from = c("sCNA", "cCNA","neu"), to = c("grey70","#7f7f7f","grey92"))) 

my_col = structure(c("grey70","grey92","#7f7f7f"), names = c("sCNA","neu","cCNA"))
my_sample_df <- as.data.frame(t(my_sample))

cs <- as.data.frame(my_sample)
ha_bottom = columnAnnotation(df = cs, col = list(my_sample=c("sCNA"="grey70", "neu" = "grey92","cCNA"="#7f7f7f")), 
                             clonal_state = anno_mark(at=gene_bin_sel_sample$rowname, labels = gene_bin_sel_sample$gene, 
                                                      side = "bottom", labels_gp = gpar(fontsize = 14, col = gene_bin_sel_sample$my_color)))

pdf(paste0("./figures/", pro_name_d, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff, "_complexHeatmap_DR.pdf"), 
    height = 10, width = 10)
print(Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        row_split = anno_mtx$subclones,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
        column_title = paste0(pro_name, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "Log2 (Ratio)", 
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row2, bottom_annotation = ha_bottom))
dev.off()
}
for (i in c(0.1,0.15)) {
  
  low_cutoff <- -i
  up_cutoff <-  i
  cell_cutoff <- 0.85
  neu_cutoff <- 0.90
  my_sample <- clonality_log_trinary_neu(log_ratio_df=ht_mtx, lower_cutoff = low_cutoff, upper_cutoff = up_cutoff, 
                                         cell_pct = cell_cutoff, neu_pct = neu_cutoff)
  table(my_sample)
  # cs <- as.data.frame(my_sample)
  sig_de.markers <- read.csv(paste0("./metrics/", pro_name_r, "_volcano_DE_sig.csv"), row.names = 1)
  gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% sig_de.markers$rowname)
  
  gene_bin_sel_sample <- as.data.frame(my_sample) %>% rownames_to_column() %>% dplyr::mutate(rowname = as.numeric(rowname)) %>% 
    dplyr::right_join(gene_bin_sel, by = c("rowname" = "pos")) %>% 
    dplyr::mutate(my_color = plyr::mapvalues(my_sample, from = c("sCNA", "cCNA","neu"), to = c("grey70","#7f7f7f","grey92"))) 
  
  my_col = structure(c("grey70","grey92","#7f7f7f"), names = c("sCNA","neu","cCNA"))
  my_sample_df <- as.data.frame(t(my_sample))
  
  cs <- as.data.frame(my_sample)
  ha_bottom = columnAnnotation(df = cs, col = list(my_sample=c("sCNA"="grey70", "neu" = "grey92","cCNA"="#7f7f7f")), 
                               clonal_state = anno_mark(at=gene_bin_sel_sample$rowname, labels = gene_bin_sel_sample$gene, 
                                                        side = "bottom", labels_gp = gpar(fontsize = 14, col = gene_bin_sel_sample$my_color)))
  
  pdf(paste0("./figures/", pro_name_d, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff, "_complexHeatmap_DR.pdf"), 
      height = 10, width = 10)
  print(Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
                row_split = anno_mtx$subclones,
                name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
                column_title = paste0(pro_name, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff),
                use_raster = T, raster_quality = 5, col = col_vec, 
                heatmap_legend_param = list(title = "Log2 (Ratio)", 
                                            title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
                top_annotation = ha_col, left_annotation = ha_row2, bottom_annotation = ha_bottom))
  dev.off()
}
for (i in c(0.88,0.92, 0.95)) {
  low_cutoff <- -0.15
  up_cutoff <-  0.15
  cell_cutoff <- i
  neu_cutoff <- 0.90
  my_sample <- clonality_log_trinary_neu(log_ratio_df=ht_mtx, lower_cutoff = low_cutoff, upper_cutoff = up_cutoff, 
                                         cell_pct = cell_cutoff, neu_pct = neu_cutoff)
  table(my_sample)
  # cs <- as.data.frame(my_sample)
  sig_de.markers <- read.csv(paste0("./metrics/", pro_name_r, "_volcano_DE_sig.csv"), row.names = 1)
  gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% sig_de.markers$rowname)
  
  gene_bin_sel_sample <- as.data.frame(my_sample) %>% rownames_to_column() %>% dplyr::mutate(rowname = as.numeric(rowname)) %>% 
    dplyr::right_join(gene_bin_sel, by = c("rowname" = "pos")) %>% 
    dplyr::mutate(my_color = plyr::mapvalues(my_sample, from = c("sCNA", "cCNA","neu"), to = c("grey70","#7f7f7f","grey92"))) 
  
  my_col = structure(c("grey70","grey92","#7f7f7f"), names = c("sCNA","neu","cCNA"))
  my_sample_df <- as.data.frame(t(my_sample))
  
  cs <- as.data.frame(my_sample)
  ha_bottom = columnAnnotation(df = cs, col = list(my_sample=c("sCNA"="grey70", "neu" = "grey92","cCNA"="#7f7f7f")), 
                               clonal_state = anno_mark(at=gene_bin_sel_sample$rowname, labels = gene_bin_sel_sample$gene, 
                                                        side = "bottom", labels_gp = gpar(fontsize = 14, col = gene_bin_sel_sample$my_color)))
  
  pdf(paste0("./figures/", pro_name_d, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff, "_complexHeatmap_DR.pdf"), 
      height = 10, width = 10)
  print(Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
                row_split = anno_mtx$subclones,
                name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
                column_title = paste0(pro_name, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff),
                use_raster = T, raster_quality = 5, col = col_vec, 
                heatmap_legend_param = list(title = "Log2 (Ratio)", 
                                            title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
                top_annotation = ha_col, left_annotation = ha_row2, bottom_annotation = ha_bottom))
  dev.off()
}
# pdf(paste0("./figures/", pro_name_d, "_complexHeatmap_DR.pdf"), height = 10, width = 10)
# Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
#         row_split = anno_mtx$subclones,
#         name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
#         column_title = paste0(pro_name, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff),
#         use_raster = T, raster_quality = 5, col = col_vec, 
#         heatmap_legend_param = list(title = "Log2 (Ratio)", 
#                                     title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
#         top_annotation = ha_col, left_annotation = ha_row2, bottom_annotation = ha_bottom)
# dev.off()

#----calculate overdispersion---######
bin_count <- varbin_mtx_tumor_log2@assays@data$bin_counts
bin_count_overdisp <- purrr::map_dfr(bin_count, overdispersion) %>% t() %>% as.data.frame()
bin_count_overdisp2 <- bin_count_overdisp %>% rownames_to_column() %>% dplyr::rename(over_disp = V1) %>% dplyr::select("over_disp")
varbin_mtx_tumor_log2@colData <- cbind(varbin_mtx_tumor_log2@colData, bin_count_overdisp2) 

saveRDS(varbin_mtx_tumor_log2, file = paste0("./objects/", pro_name_d, c("_final_filtered_overdisp_copykit.rds")))
# varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_overdisp_copykit.rds")))

#---calculate basic QC matrics----####
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
varbin_mtx_tumor_log2@colData$experiment %>% table()

varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::filter(experiment == "nanowell1") %>% pull(reads_total) %>% mean()
varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::filter(experiment == "nanowell2") %>% pull(reads_total) %>% mean()

varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::filter(experiment == "nanowell1") %>% pull(reads_assigned_bins) %>% mean()
varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::filter(experiment == "nanowell2") %>% pull(reads_assigned_bins) %>% mean()

varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::filter(experiment == "nanowell1") %>% pull(percentage_duplicates) %>% mean()
varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::filter(experiment == "nanowell2") %>% pull(percentage_duplicates) %>% mean()

varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::filter(experiment == "nanowell1") %>% pull(median_bin_count) %>% mean()
varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::filter(experiment == "nanowell2") %>% pull(median_bin_count) %>% mean()

obj_rna <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
dna_meta_temp <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::select("rna_bc2","superclones","subclones","sample","experiment")

full_list <- obj_rna@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  full_join(dna_meta_temp, by = c("rowname" = "rna_bc2")) %>% column_to_rownames() 


#---correlation scDNA profiles to copykat inferred CNA profiles----####
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
final <- read.delim("./data/copykat/wdr_mda231/wdr_nanowell_rna_copykat_CNA_results.txt", header=T)

wdr_dna <- log2(varbin_mtx_tumor_log2@assays@data$segment_ratios)
copycat_dna <- as.data.frame(final[, 4:ncol(final)])

seg_wdr_dna <- apply(wdr_dna, 1, median)
copycat_dna <- apply(copycat_dna, 1, median)
cor(seg_wdr_dna, copycat_dna)


#---Statistic counts of DNA RNA matching------#########
obj_rna <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
varbin_mtx_tumor_log3 <- readRDS(paste0("objects/", pro_name_d, c("_tumor_integerCN_copykit.rds")))
dna_meta_temp <- varbin_mtx_tumor_log3@colData %>% as.data.frame() %>% dplyr::select("rna_bc2","superclones","subclones","sample")

gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
bin_chrpos <- varbin_mtx_tumor_log3@rowRanges %>% as.data.frame() %>% rownames_to_column() %>% 
  mutate(chr_pos = paste0("chr",seqnames,":",start, "-", end)) %>% dplyr::select(rowname, chr_pos)

obj_rna2 <- obj_rna
obj_rna2@meta.data <- obj_rna2@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  left_join(dna_meta_temp, by = c("rowname" = "rna_bc2")) %>% column_to_rownames() 
subclone_count <- obj_rna2@meta.data$subclones %>% table() 

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
  obj_rna3@meta.data <- obj_rna3@meta.data %>% mutate(temp = ifelse(subclones %in% super_clone1, "super_c1",
                                                                    ifelse(subclones %in% super_clone2, "super_c2", "others")))
  clone1_cells <- obj_rna3@meta.data %>% dplyr::filter(temp == "super_c1") %>% rownames()
  rna_exp_clone1_cells <- obj_rna3@assays$RNA@data[,clone1_cells]
  clone2_cells <- obj_rna3@meta.data %>% dplyr::filter(temp == "super_c2") %>% rownames()
  rna_exp_clone2_cells <- obj_rna3@assays$RNA@data[,clone2_cells]
  
  rna_cons_clone1 <- apply(rna_exp_clone1_cells, 1, mean) %>% as.data.frame() %>% `colnames<-`("mean_scaled_exp_clone1") %>% 
    rownames_to_column() %>% dplyr::filter(mean_scaled_exp_clone1 > 0)
  rna_cons_clone2 <- apply(rna_exp_clone2_cells, 1, mean) %>% as.data.frame() %>% `colnames<-`("mean_scaled_exp_clone2") %>% 
    rownames_to_column() %>% dplyr::filter(mean_scaled_exp_clone2 > 0)
  
  clone1_2_con_exp <- inner_join(rna_cons_clone1, rna_cons_clone2, by = "rowname")
  dna_con <- varbin_mtx_tumor_log3@consensus %>% rownames_to_column() %>% dplyr::select(rowname,all_of(c(super_clone1, super_clone2)))
  
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
# dna_seg_rna_df <- read.csv(paste0("./metrics/",pro_name, "_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_", min_cell_num, "_Cells.csv"), row.names = 1)

dna_seg_rna_df <- dna_seg_rna_df %>% mutate(matched = as.integer(matched), dna_but_no_rna = as.integer(dna_but_no_rna), 
                                            no_dna_but_rna = as.integer(no_dna_but_rna))
sum_matched <- sum(dna_seg_rna_df$matched)
sum_dna_but_no_rna <- sum(dna_seg_rna_df$dna_but_no_rna)
sum_no_dna_but_rna <- sum(dna_seg_rna_df$no_dna_but_rna)

sum_matched/(sum_matched + sum_dna_but_no_rna)
# 25 cells: 0.7319277
# 50 cells: 0.7075472

#---venn diagram---
#---venn diagram---
library("VennDiagram")
library("RColorBrewer")

area1_value = sum_dna_but_no_rna + sum_matched
area2_value = sum_no_dna_but_rna + sum_matched
if(area1_value > area2_value) {
  category_names <- c("CNA diff", "RNA diff")
} else {
  category_names <- c("RNA diff", "CNA diff")
}

p1 <- draw.pairwise.venn(area1 = area1_value, area2 = area2_value, cross.area = sum_matched, fill = brewer.pal(3, "Set2")[2:1], 
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
obj_rna <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
dna_meta_temp <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::select("rna_bc2","superclones","subclones","sample")

obj_rna2 <- obj_rna
obj_rna2@meta.data <- obj_rna2@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  left_join(dna_meta_temp, by = c("rowname" = "rna_bc2")) %>% column_to_rownames() 
subclone_count <- obj_rna2@meta.data$subclones %>% table() 

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
  obj_rna3@meta.data <- obj_rna3@meta.data %>% mutate(temp = ifelse(subclones %in% super_clone1, "super_c1",
                                                                    ifelse(subclones %in% super_clone2, "super_c2", "others")))
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
medic_dis <- read_tsv(paste0("./metrics/medicc_files/", pro_name_d, "/", pro_name_d, "_medicc2_input_pairwise_distances.tsv"))
medic_dis2 <- medic_dis %>% column_to_rownames(var="sample_id")
medic_dis2 <- medic_dis2[1:(ncol(medic_dis2)-1), 1:(ncol(medic_dis2)-1)] %>% as.matrix()

melted_data <- reshape2::melt(medic_dis2) %>% mutate(subclone_comb = paste0(Var1, "_",Var2)) %>% dplyr::select(value, subclone_comb)
colnames(melted_data) <- c("cna_dist", "subclone_comb")

genetic_phyno_mtx <- left_join(de.markers_df, melted_data, by = c("subclone_comb" = "subclone_comb"))

write.csv(genetic_phyno_mtx, paste0("./metrics/",pro_name, "_genetic_distance_phenotype_diff_subclonesWithMorethan_", min_cell_num, "_Cells.csv"))
# genetic_phyno_mtx <- read.csv(paste0("./metrics/",pro_name_r, "_genetic_distance_phenotype_diff_subclonesWithMorethan50Cells.csv"), row.names = 1)

genetic_phyno_mtx2 <- genetic_phyno_mtx %>% as.data.frame() %>% mutate(log2fc_0.485 = as.integer(log2fc_0.485))
p1 <- ggplot(genetic_phyno_mtx2, aes(x = cna_dist, y = log2fc_0.485)) +
  geom_point() + geom_line() + sm_statCorr() +
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
obj_rna <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
varbin_mtx_tumor_log3 <- readRDS(paste0("objects/", pro_name_d, c("_tumor_integerCN_copykit.rds")))
dna_meta_temp <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::select("rna_bc2","superclones","subclones","sample")
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)

obj_rna2 <- obj_rna
obj_rna2@meta.data <- obj_rna2@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  left_join(dna_meta_temp, by = c("rowname" = "rna_bc2")) %>% column_to_rownames() 

subclone_count <- obj_rna2@meta.data$subclones %>% table() 
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
  
  dna_con <- varbin_mtx_tumor_log3@consensus %>% rownames_to_column() %>% dplyr::select(rowname,all_of(c(super_clone1, super_clone2))) %>% 
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
obj_rna <- readRDS(paste0("./objects/",pro_name_r, "_seurat_obj.rds"))
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_DNA_filtered_copykit.rds")))
dna_meta_temp <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% dplyr::select("rna_bc2","superclones","subclones","sample")

obj_rna2 <- obj_rna
obj_rna2@meta.data <- obj_rna2@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  left_join(dna_meta_temp, by = c("rowname" = "rna_bc2")) %>% column_to_rownames() 
subclone_count <- obj_rna2@meta.data$subclones %>% table() 

#---use minimum 25 cell with RNA data as subclone cut-off---
min_cell_num <- 25
subsample_min <- min(subclone_count[subclone_count > min_cell_num])
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
  obj_rna3@meta.data <- obj_rna3@meta.data %>% mutate(temp = ifelse(subclones %in% super_clone1, "super_c1",
                                                                    ifelse(subclones %in% super_clone2, "super_c2", "others")))
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
medic_dis <- read_tsv(paste0("./metrics/medicc_files/", pro_name_d, "/", pro_name_d, "_medicc2_input_pairwise_distances.tsv"))
medic_dis2 <- medic_dis %>% column_to_rownames(var="sample_id")
medic_dis2 <- medic_dis2[1:(ncol(medic_dis2)-1), 1:(ncol(medic_dis2)-1)] %>% as.matrix()

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

# genetic_phyno_mtx_long <- melt(genetic_phyno_mtx, id.vars = c("subclone_comb", "cna_dist")) %>% as.data.frame() %>% mutate(value = as.integer(value))
# 
# # Create the scatter plot
# ggplot(genetic_phyno_mtx_long, aes(x = cna_dist, y = value, color = variable)) +
#   geom_point() +
#   labs(x = "cna_dist", y = "Value") +
#   theme_minimal()




