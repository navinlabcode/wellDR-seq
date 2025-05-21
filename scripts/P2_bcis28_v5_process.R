#---preload----#####
library("Seurat")
library(SeuratWrappers)
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
options(max.print = 200)

wd <- c("./wscDR/")
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

pro_name <- "bcis28t2chips"
pro_name1 <- "bcis28t_chip1"
pro_name2 <- "bcis28t_chip2"


pro_name_r <- paste0(pro_name, "_rna")
pro_name_d <- paste0(pro_name, "_dna")


###---DNA--------#########
pro_name_d <- paste0(pro_name, "_dna")
raw_path1 <- paste0("./map_seg_output/", pro_name1)
raw_path2 <- paste0("./map_seg_output/", pro_name2)


#-----DNA processing (no filtering applied to all cells)------####
#---t1--
varbin_mtx1 <- readVarbin(raw_path1, remove_Y = TRUE)
varbin_mtx1 
varbin_mtx_all1 <- varbin_mtx1
varbin_mtx_all1@colData 

#---t2--
varbin_mtx2 <- readVarbin(raw_path2, remove_Y = TRUE)
varbin_mtx2 
varbin_mtx_all2 <- varbin_mtx2[,filt_cells_names_dis2]
varbin_mtx_all2@colData 



#------merge no filtered datasets----#####
varbin_mtx_all1@colData$experiment <- pro_name1
varbin_mtx_all2@colData$experiment <- pro_name2
# varbin_mtx_all3@colData$experiment <- pro_name3
# varbin_mtx_all4@colData$experiment <- pro_name4

varbin_mtx_all <- cbind(varbin_mtx_all1, varbin_mtx_all2)#,varbin_mtx_all3)#, varbin_mtx_all4)

varbin_mtx_all_log <- logNorm(varbin_mtx_all, transform = "log2")
saveRDS(varbin_mtx_all_log, file = paste0("./objects/", pro_name_d, c("_without_filtered_all_copykit_rnameta.rds")))
# varbin_mtx_all_log <- readRDS(paste0("./objects/", pro_name_d, c("_without_filtered_all_copykit_rnameta.rds")))

#------Knn smooth-----#######
metadata(varbin_mtx_all_log)$genome <-  "hg19"
varbin_mtx_all_log_knn <- knnSmooth2(varbin_mtx_all_log)
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

saveRDS(varbin_mtx_normal, file = paste0("objects/", pro_name_d, c("_knn_normal_copykit_before_rename.rds")))
# varbin_mtx_normal <- readRDS(paste0("objects/", pro_name_d, c("_knn_normal_copykit_before_rename.rds")))

saveRDS(varbin_mtx_tumor, file = paste0("objects/", pro_name_d, c("_knn_tumor_copykit_before_rename.rds")))
# varbin_mtx_tumor <- readRDS(paste0("objects/", pro_name_d, c("_knn_tumor_copykit_before_rename.rds")))

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
names(exp_col) <- c(pro_name1, pro_name2)#, pro_name3)#, pro_name4)
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
names(exp_col) <- c(pro_name1, pro_name2)#, pro_name3)#, pro_name4)
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
hdb_data <- dbscan::hdbscan(umap_data2, minPts = 10)
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
  cowplot::ggsave2(paste0("./figures/", pro_name_d, "_passDNAQC_knn_tumor_cells_umap.pdf"), p1, width = 5, height = 4)
  p1
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
names(exp_col) <- c(pro_name1, pro_name2)#, pro_name3)#, pro_name4)
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

#----------Remove bad clusters----#######
bad_sub <- c("c6","c8")
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

# clone1 <- "c7"
# clone2 <- "c8"
# con_tmp0 <- con_tmp %>% rownames_to_column()
# diff_loc <- con_tmp0$rowname[con_tmp[,clone1] != con_tmp[,clone2]]
# varbin_mtx_tumor_log3@rowRanges %>% as.data.frame() %>% rownames_to_column() %>% dplyr::filter(rowname %in% diff_loc) %>%
#      left_join(con_tmp0[,c("rowname",clone1,clone2)], by = c("rowname" = "rowname"))

to_clone1 <- "c3"
from_clone1 <- paste0("c",4:5)

subclones2_tmp <- plyr::mapvalues(varbin_mtx_tumor_log3@colData$subclones, 
                                  from = from_clone1, 
                                  to = rep(to_clone1, times = length(from_clone1)))

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
names(exp_col) <- c(pro_name1, pro_name2)#, pro_name3)#, pro_name4)
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

#----------Calculate computational ploidy----##
varbin_mtx_all_log_tmp <- varbin_mtx_all_log
scq.call <- BiocParallel::bplapply(colnames(assay(varbin_mtx_all_log_tmp, "bin_counts")),function(x){
  ## call scquantum, subset the scaffold data if chrY is removed
  a <- ploidy.inference(assay(varbin_mtx_all_log_tmp, "bin_counts")[,x], hg19_rg$chr[1:12167], hg19_rg$start[1:12167], hg19_rg$end[1:12167])
  c(x,unlist(a[c(4,11)]))
})
scq_mat <- as_tibble(do.call(rbind, scq.call))
scq_mat2 <- as.data.frame(scq_mat) %>% mutate(ploidy_pre2 = as.numeric(ploidy), confidence_ratio2 = as.numeric(confidence_ratio)) %>% 
  dplyr::select("ploidy_pre2","confidence_ratio2")
varbin_mtx_all_log_tmp@colData <- cbind(varbin_mtx_all_log_tmp@colData, scq_mat2)

p2 <- varbin_mtx_all_log_tmp@colData %>% as.data.frame() %>% ggplot(aes(x=ploidy_pre2, y = confidence_ratio2)) + 
  geom_point(size=0.1)  + facet_grid(~subclones2)
p2
cowplot::ggsave2(paste0("./figures/", pro_name_d, "_scquantum_ploidy_confidence_by_subclone2_all_cells.pdf"), p2, 
                 width = length(unique(varbin_mtx_all_log_tmp@colData$subclones)), height = 2)

saveRDS(varbin_mtx_all_log_tmp, file = paste0("objects/", pro_name_d, c("_final_merged_before_rename_copykit_withdoublets.rds")))
# varbin_mtx_all_log_tmp <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_before_rename_copykit_withdoublets.rds")))

bad_sub <- c("nc2")
#----------Remove bad clusters----##
my_singlets_name <- varbin_mtx_all_log_tmp@colData %>% as.data.frame() %>% dplyr::filter(!(subclones2 %in% bad_sub)) %>% pull(sample)
varbin_mtx_all_log <- varbin_mtx_all_log[,my_singlets_name]

saveRDS(varbin_mtx_all_log, file = paste0("objects/", pro_name_d, c("_final_merged_before_rename_copykit_remove_doublets.rds")))
# varbin_mtx_all_log <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_before_rename_copykit_remove_doublets.rds")))

#---rename subclones---
varbin_mtx_all_log@colData$subclones3 <- plyr::mapvalues(varbin_mtx_all_log@colData$subclones2, 
                                                         from = varbin_mtx_all_log@colData$subclones2 %>% unique() %>% sort(), 
                                                         to = paste0("c", 1:length(unique(varbin_mtx_all_log@colData$subclones2)))) %>% 
  factor(levels = paste0("c", 1:length(unique(varbin_mtx_all_log@colData$subclones2))))

ploidy_sub <- varbin_mtx_all_log@colData %>% as.data.frame() %>% 
  dplyr::select("subclones3", "ploidy_pre") %>% `rownames<-`( NULL ) %>% as_tibble() %>% dplyr::group_by(subclones3) %>% 
  dplyr::summarise(sub_ploidy = median(ploidy_pre))
ploidy_sub
ploidy_sub <- ploidy_sub %>% mutate(sub_ploidy = ifelse(sub_ploidy<2, 2, sub_ploidy))

varbin_mtx_all_log@colData$ploidy_subclone <- plyr::mapvalues(varbin_mtx_all_log@colData$subclones3, from = ploidy_sub$subclones3, to = ploidy_sub$sub_ploidy) %>% 
  as.character() %>% as.numeric()

varbin_mtx_all_log2 <- calcInteger(varbin_mtx_all_log, assay = "segment_ratios", method = "fixed", ploidy_value = varbin_mtx_all_log@colData$ploidy_subclone)
varbin_mtx_all_log2 <- calcConsensus(varbin_mtx_all_log2, assay = "integer", fun="median", consensus_by = "subclones3")

saveRDS(varbin_mtx_all_log2, file = paste0("objects/", pro_name_d, c("_merged_integerCN_copykit_new.rds")))

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
names(exp_col) <- c(pro_name1, pro_name2)#, pro_name3)#, pro_name4)
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
names(exp_col) <- c(pro_name1, pro_name2)#, pro_name3)#, pro_name4)
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
b <- read.table("./pre_load_data/chr_a.txt")
a <- unlist(b, use.names = F)
chr2 <- read.table("./pre_load_data/chr_chr2.txt")
ha_col=HeatmapAnnotation(foo=anno_text(a, rot = 0, gp = gpar(fontsize =10)), df =chr2, col = list(chr=c("1"="black", "2"="grey")), 
                         show_legend = F, annotation_name_side = "left")

#-----annotated genes---
gene_bin <- read.table("./pre_load_data/hg19_gene_binpos_map.tsv", header = 1)

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
        top_annotation = ha_col, left_annotation = ha_row_cs)
dev.off()



#----------MEDICC2 Tree----######
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))
## make medicc input matrix
medicc_input <- cbind(varbin_mtx_all_log2@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_all_log2@consensus) %>%
  mutate(diploid=c1) %>%          ## assign diploid cell CN profile to be used as root in the medicc tree
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())

#---remove boarder regions (top and bottom 10 bin of each chr from each subclone)---#######
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
medicc_input <- read_tsv(file = paste0("metrics/medicc_files/", pro_name_d, "/",pro_name_d, c("_medicc2_input.tsv")))

# 
# 
# write_tsv(medicc_input2, file = paste0("metrics/medicc_files/", pro_name_d, "/",pro_name_d, c("_medicc2_no_boarder_input.tsv")))
# medicc_input2 <- read_tsv(file = paste0("metrics/medicc_files/", pro_name_d, "/",pro_name_d, c("_medicc2_no_boarder_input.tsv")))

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

cowplot::ggsave2(paste0("./figures/", pro_name_d, "_medicc2_tree_merged_subclones.pdf"), treeplt, width = 10, height = 4)


#----------MEDICC2 ME Tree---lumHR----######
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))
## make medicc input matrix
#--other clusters--
other_clst <- c("c2","c3")
dip_clst <- c("c1")
medicc_input <- cbind(varbin_mtx_all_log2@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_all_log2@consensus) %>%
  mutate(diploid=2) %>%       ## assign diploid cell CN profile to be used as root in the medicc tree
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything()) %>% 
  dplyr::filter(!(sample_id %in% other_clst))

#---remove boarder regions (top and bottom 10 bin of each chr from each subclone)---#######
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

cowplot::ggsave2(paste0("./figures/", pro_name_d, "_medicc2_ME_tree_merged_subclones_lumhr.pdf"), treeplt, width = 6, height = 2.5)

####----RNA data integration--------#######
pro_name_r <- paste0(pro_name, "_rna")
rna.cluster.ids<-c("Fibro","LumHR","LumSec","MyoEpi","ECs","Tumor","Tumor2","Tumor3","RNA_miss","Unknown1","Unknown2","Unknown3","Mito")
rna.cluster.col <- c("#D89000","#FF62BC","#00BD5F","#9590FF","#F8766D","#54278f","#BF4000","#3969AC","grey90","grey60","grey70","grey80","#CAB2D6")
names(rna.cluster.col) <- rna.cluster.ids

obj_rna <- readRDS("./objects/BCIS28T_RNA_2chip.sct.rds")

varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit.rds")))
rna_celltype <- c("LumHR","MyoEpi","Fibro","LumSec")

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
varbin_mtx_all_log2 <- readRDS(paste0("objects/", pro_name_d, c("_final_merged_integerCN_copykit_RNA_meta.rds")))


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
names(exp_col) <- c(pro_name1, pro_name2)#, pro_name3, pro_name4)

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
b <- read.table("./pre_load_data/chr_a.txt")
a <- unlist(b, use.names = F)
chr2 <- read.table("./pre_load_data/chr_chr2.txt")
ha_col=HeatmapAnnotation(foo=anno_text(a, rot = 0, gp = gpar(fontsize =10)), df =chr2, col = list(chr=c("1"="black", "2"="grey")), 
                         show_legend = F, annotation_name_side = "left")

#-----annotated genes---
gene_bin <- read.table("./pre_load_data/hg19_gene_binpos_map.tsv", header = 1)


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
        top_annotation = ha_col, left_annotation = ha_row_cs, right_annotation = ha_barplot)
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

pdf(paste0("./figures/", pro_name_d, "_merged_consensus_integer_complexHeatmap_combined_DR_larger2cells.pdf"), height = 4, width = 10)
Heatmap(as.matrix(cs_mtx_order), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        name = "csheatmap", show_row_names = T, row_names_side = "left",
        show_column_names = F, column_title = paste0(pro_name_d, "_consensusHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row_cs, bottom_annotation = ha_bottom_cs, right_annotation = ha_barplot2)
dev.off()

###-----RNA data UMAP----########
obj_rna_newid<-obj_rna
rna_celltype <- c("LumHR","MyoEpi","Fibro","LumSec")
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







