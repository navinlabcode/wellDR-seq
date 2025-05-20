#---preload----#####
options(bitmapType="cairo")
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
library(scquantum)
library(readr)
library(png)
library(beeswarm)
library(useful)
library(parallel)
options(max.print = 200)
source("./0_wdr_functions.R")
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
hbca_epi_markers <- read.csv("./pre_load_data/HBCA_Epi_markers.csv")
new_pal2 <- new_pal
names(new_pal2) <- paste0("c", 1:length(new_pal2))

#---The correlation between the copy number segment difference relative to diploid and the normalized mean expression difference(Figure 5I)-----####
gene_bin <- read.table("./pre_load_data/hg19_gene_binpos_map.tsv", header = 1) 
names(gene_bin) <- c("gene", "bins")

process_sample <- function(sample_path, sample_name) {
  tryCatch({
    cat(sprintf("Processing sample: %s\n", sample_name))
    
    # Load RNA object for the given sample
    obj_rna <- readRDS(sample_path)
    # Subclones' names
    subclones <- as.vector(obj_rna@meta.data$subclone3)
    # Gene expression matrix for each cell
    object_bulk <- as.data.frame(t(obj_rna@assays$RNA@data))
    
    # Add subclone index
    object_bulk$subclone <- subclones
    
    # Replace zeros with NA for mean and median calculations
    all_genes <- colnames(object_bulk)
    # mean expression of each subclone
    # object_bulk_remove_zero <- object_bulk %>%
    #       dplyr::mutate(across(where(is.numeric), ~replace(., . == 0, NA)))
    object_bulk_mean <- object_bulk %>%
          dplyr:: group_by(subclone) %>%
          dplyr:: summarise(across(everything(),   ~mean(., na.rm = TRUE)))
    object_bulk_mean <- as.data.frame(t(object_bulk_mean))
    # Set column names and remove first row (which contains subclone names)
    #object_bulk_mean <- object_bulk_mean %>% column_to_rownames("subclones")
    colnames(object_bulk_mean) <- object_bulk_mean[1, ]
    object_bulk_mean <- object_bulk_mean[-1, ]
    object_bulk_mean$gene <- rownames(object_bulk_mean)
    
    
    # Get rid of unassigned subclone cells
    object_bulk_mean <- object_bulk_mean %>%
      select(-matches("^(unassigned|unsigned)$"))
    
    cat(sprintf("Start processing dna\n"))
    #sample_name <- "ecis44t"
    pro_name_d <- paste0(sample_name, "_dna")
    varbin_mtx_all_log2 <- readRDS(paste0("./objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))
    gene_bin_sel <- gene_bin %>% 
      dplyr::filter(gene %in% all_genes )
                         
    raw_mtx_cnv <- as.data.frame(varbin_mtx_all_log2@consensus)

    cutoff <- c(5, 10, 15, 20)
    clone_name <- vector("character", length = 0)
    cat(sprintf("Start processing cnv_cutoff\n"))
    for (j in 1:4){
      #i<-4
      cut <- cutoff[j]
      col <- as.data.frame(varbin_mtx_all_log2@colData) %>%
        select(sample,subclones3,nCount_RNA) %>%
        dplyr::group_by(subclones3) %>%
        dplyr::summarise(cell_count = sum(!is.na(nCount_RNA))) %>%
        filter(cell_count >= cut)
      #assign(paste0("clone_name_with_cutoff_",cut),col$subclones3)
      clone_name <- col$subclones3

    
      ## compared with diploid c1
      cnv_changes_c1 <- raw_mtx_cnv
      # store new colname for cnv-c1
      new_col <- vector("character", length = 0)
      # store new colname for expression of each clone
      x <- vector("character", length = 0)
      # cnv minus diploid
      for (i in 2:(length(clone_name))) {
          # generate new colname
          #i<- 8
          new_col_name <- paste0("cnv_", clone_name[i], "_c1")
          new_col <- append(new_col, new_col_name)
          x <- append(x,paste0(clone_name[i],"_exp"))
          # calculate each subclones' cnv minus c1 and store them in new cols
          cnv_changes_c1[new_col_name] <- raw_mtx_cnv %>% select(clone_name[i]) - raw_mtx_cnv$c1
        
      }
      # add bins information
      cnv_changes_c1$bins <- as.numeric(gsub("V","",rownames(cnv_changes_c1)))
      # find genes for each bin
      names(gene_bin_sel) <- c("gene", "bins")
      cnv_genes_mtx_c1 <- left_join(cnv_changes_c1, gene_bin_sel, by = "bins")
      # find expression data for each gene
      cnv_genes_mean_exp_mtx_c1 <- left_join(cnv_genes_mtx_c1,object_bulk_mean, by = "gene" )
      names(cnv_genes_mean_exp_mtx_c1) <- gsub(".y$", "_exp", names(cnv_genes_mean_exp_mtx_c1))
      names(cnv_genes_mean_exp_mtx_c1) <- gsub(".x$", "_cnv", names(cnv_genes_mean_exp_mtx_c1))

      
      
      #MEAN
      # calculate mean expression of each bin
      exp_columns <- grepl("exp$",colnames(cnv_genes_mean_exp_mtx_c1))
      cnv_genes_mean_exp_mtx_c1[exp_columns] <- lapply(cnv_genes_mean_exp_mtx_c1[exp_columns], as.numeric)
      cnv_genes_mean_exp_mtx_c1 <- cnv_genes_mean_exp_mtx_c1 %>%
                      dplyr::mutate(across(ends_with("_exp"), ~ ifelse(. > 0.05, ., NA)))
      gene_per_bin <- cnv_genes_mean_exp_mtx_c1 %>%
                      dplyr::group_by(bins) %>%
                      dplyr::summarise(count=n()) 
      median(gene_per_bin$count)
      # calculate expression minus c1 for each clone
      for (i in 1:(length(clone_name)-1)) {
        new_col_name <- paste0("sub_", new_col[i])
        x_col <- x[i] 
        cnv_genes_mean_exp_mtx_c1[[new_col_name]] <- cnv_genes_mean_exp_mtx_c1[[x_col]] - cnv_genes_mean_exp_mtx_c1$c1_exp
        
      }
      
      for (i in 1:(length(clone_name)-1)) {
        group <- 1
        group_name <- paste0("group_",new_col[i])
        cnv <- new_col[i]
        for (j in 1:(nrow(cnv_genes_mean_exp_mtx_c1) - 1)) {
          cnv_genes_mean_exp_mtx_c1[[group_name]][j] <- group
          
          if (cnv_genes_mean_exp_mtx_c1[[cnv]][j+1] == cnv_genes_mean_exp_mtx_c1[[cnv]][j]) {
            cnv_genes_mean_exp_mtx_c1[[group_name]][j+1] <- group
          } else {
            group <- group + 1
            cnv_genes_mean_exp_mtx_c1[[group_name]][j+1] <- group
          }
        }
      }
      sample_cnv_genes_mean_exp_mtx_c1_bin <-  cnv_genes_mean_exp_mtx_c1
      
      # separate to different mtxes
      sample_mean_mtx_name <- vector("character", length = 0)
      for (i in 1:(length(clone_name)-1)) {
            group <- paste0("group_",new_col[i])
            exp_sub_col <- paste0("sub_", new_col[i])
            cnv_col <- new_col[i]
            name <- paste0("cnv_genes_exp_mean_mtx_bin_",new_col[i])
            mtx <- data.frame(sample_cnv_genes_mean_exp_mtx_c1_bin[group],
                                   sample_cnv_genes_mean_exp_mtx_c1_bin[cnv_col],sample_cnv_genes_mean_exp_mtx_c1_bin[exp_sub_col],
                                   bins = sample_cnv_genes_mean_exp_mtx_c1_bin$bins)
            mtx_group <- mtx %>%
                dplyr::group_by(.data[[group]]) %>%
                dplyr::summarise(across(.cols = starts_with("sub_"), ~median(., na.rm = TRUE), .names = "{.col}"))
            mtx <- left_join(mtx,mtx_group,by =group)
            colnames(mtx) <- c("group","cnv","sub_bin","bins","sub_group")
            assign(name,mtx)
            sample_mean_mtx_name <- append(sample_mean_mtx_name,name)
      }
          
  
      #add all clones to one mtx for mean
      cnv_mtx <- data.frame()
      for(i in 1:(length(clone_name)-1)){
        #i <-2 
        mtx <- get(sample_mean_mtx_name[i])
        colnames(mtx) <- c("group","cnv","sub_bin","bins","sub_group")
        filter_seg_size <- mtx %>%
          dplyr::group_by(bins) %>%
          dplyr::summarise(count=n(),group = max(group)) %>%
          dplyr::group_by(group) %>%
          dplyr::summarise(count=n()) %>%
          filter(count >10)
        mtx <- mtx %>%
          filter(group %in% filter_seg_size$group)
        mtx_tmp <- mtx %>%
          select(-bins,-sub_bin) %>%
          distinct() %>%
          select(-group)
        cnv_mtx <- rbind(cnv_mtx, mtx_tmp)
      }
      cnv_mtx$cnv <- factor(cnv_mtx$cnv , sort(unique(cnv_mtx$cnv)))
      
      cnv_mean_mtx_sample <- cnv_mtx
      cnv_mean_mtx_sample$sample <- sample_name
      write.csv(cnv_mean_mtx_sample,paste0("./cnv_exp/",sample_name,"_cnv_new_mtx_with_cutoff_", cut, ".csv"),row.names = F)

      
      
  }
    cat(sprintf("Finished processing sample: %s\n", sample_name))
  }, error = function(e) {
    cat(sprintf("Error processing sample: %s - %s\n", sample_name, e$message))
  })
}

# Define paths and names for all samples
samples <- list(
  list(path = "objects/bcis66t_rna_final_merged_integerCN.rds", name = "bcis66t"),
  list(path = "objects/bcis106t_rna_final_merged_integerCN.rds", name = "bcis106t"),
  list(path = "objects/bcis74t_rna_final_merged_integerCN.rds", name = "bcis74t"),
  list(path = "objects/bcis51t_rna_final_merged_integerCN.rds", name = "bcis51t"),
  list(path = "objects/bcis70t_rna_final_merged_integerCN.rds", name = "bcis70t"),
  list(path = "objects/ecis57t_rna_final_merged_integerCN.rds", name = "ecis57t"),
  list(path = "objects/BCIS13T_Chip1_2_RNA.rds", name = "bcis13t"),
  list(path = "objects/BCIS28T_Chip1_2_RNA.rds", name = "bcis28t"),
  list(path = "objects/ECIS48T_Chip1_2_RNA.rds", name = "ecis48t"),
  list(path = "objects/ECIS25T_Chip1_2_3_RNA.rds", name = "ecis25t"),
  list(path = "objects/ECIS36T_Chip1_2_RNA.rds", name = "ecis36t"),
  list(path = "objects/ECIS44T_Chip1_2_3_4_5_RNA.rds", name = "ecis44t")
)

lapply(samples, function(sample) {
  process_sample(sample$path, sample$name)
})


# plot of segment level
## read files
new_samples <- c("bcis51t","bcis66t","bcis106t","bcis70t","bcis74t","ecis57t","bcis13t", "bcis28t", "ecis25t", "ecis36t", "ecis44t", "ecis48t")
cutoffs <- c(5, 10, 15, 20)

for (sample in new_samples) {
  for (cutoff in cutoffs) {
    file_name <- paste0("./cnv_exp/",sample,"_cnv_new_mtx_with_cutoff_", cutoff, ".csv")
    if (file.exists(file_name)) {  
      assign(paste0("cnv_new_mtx_", sample, "_with_cutoff_", cutoff),read.csv(file_name))
    } else {
      cat(paste0("File ", file_name, " does not exist.\n"))
    }
  }
}

## with_cutoff_20
data_list <- list(
  bcis13t = cnv_new_mtx_bcis13t_with_cutoff_20,
  bcis28t = cnv_new_mtx_bcis28t_with_cutoff_20,
  ecis25t = cnv_new_mtx_ecis25t_with_cutoff_20,
  ecis36t = cnv_new_mtx_ecis36t_with_cutoff_20,
  ecis44t = cnv_new_mtx_ecis44t_with_cutoff_20,
  bcis51t = cnv_new_mtx_bcis51t_with_cutoff_20,
  bcis66t = cnv_new_mtx_bcis66t_with_cutoff_20,
  bcis106t = cnv_new_mtx_bcis106t_with_cutoff_20,
  bcis70t = cnv_new_mtx_bcis70t_with_cutoff_20,
  bcis74t = cnv_new_mtx_bcis74t_with_cutoff_20,
  ecis57t = cnv_new_mtx_ecis57t_with_cutoff_20
)

processed_data_list <- lapply(names(data_list), function(sample_name) {
  data <- data_list[[sample_name]]
  processed_data <- data %>%
    dplyr::group_by(cnv) %>%
    dplyr::summarise(sub_group = mean(sub_group, na.rm = TRUE)) %>%
    mutate(sample = sample_name)
  
  return(processed_data)
})

cnv_mtx_new_mean_scale_with_cutoff_20 <- bind_rows(processed_data_list)
cnv_mtx_new_mean_scale_with_cutoff_20$cnv <- factor(cnv_mtx_new_mean_scale_with_cutoff_20$cnv,sort(unique(cnv_mtx_new_mean_scale_with_cutoff_20$cnv)))

new_clone_color <- c("#004529","#088247","#7CC767",
"#223D6C","#5D90BA","#9ecae1",
"#8c510a","#bf812d","#D8D155",
"#b30000","#d73027","#fc9272",
"#7A142C","#E0367A","#df65b0",
"#4a1486","#6a51a3","#9e9ac8",
"#252525","#737373")

new_samples <- c("bcis51t","bcis66t","bcis106t","bcis70t","bcis74t","ecis57t","bcis13t", "bcis28t", "ecis25t", "ecis36t", "ecis44t", "ecis48t")
pnames<-        c("P11","P9","P10","P8","P7","P12","P5","P2","P1","P6","P4","P3")
pname_mtx <- data.frame(pnames,sample = new_samples)
pname_mtx$pnames <- factor(pname_mtx$pnames, levels = paste0("P",1:12))
sample_platte <- new_clone_color[1:length(pname_mtx$pnames )]
names(sample_platte) <- pnames

cnv_mtx_new_mean_scale_with_cutoff_20 <- left_join(cnv_mtx_new_mean_scale_with_cutoff_20, pname_mtx, by = "sample") 


library(ggpubr)

library(mgcv)
model <- gam(sub_group ~ s(cnv), data = cnv_mtx_new_mean_scale_with_cutoff_20, method = "REML")
predicted_values <- predict(model)
actual_values <- cnv_mtx_new_mean_scale_with_cutoff_20$sub_group

cor_test <- cor.test(predicted_values, actual_values)
print(cor_test)
# 	Pearson's product-moment correlation
# 
# data:  predicted_values and actual_values
# t = 20.146, df = 59, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8925911 0.9602603
# sample estimates:
#       cor 
# 0.9343869 


p4 <- ggplot(cnv_mtx_new_mean_scale_with_cutoff_20) +
  geom_point( aes(x = cnv, y = sub_group,colour = pnames)) + 
  geom_line(aes(x = cnv, y = sub_group,colour = pnames)) +  
  scale_colour_manual(values = sample_platte)+
  scale_x_continuous(breaks = c(-1:22))+
  geom_smooth(aes(x = cnv, y = sub_group), method = "gam", na.rm = T, se = T,color = "black")  +
  theme_bw() +
  theme(
    panel.grid = element_blank(),          
    axis.line = element_line(color = "black")  
  ) 
p4
cowplot::ggsave2(paste0("./figures/cnv_exp_with_cutoffs_20.pdf"), p4, width = 5, height = 3)




#---Statistic counts of DNA RNA matching------#########
process_sample <- function(sample_path, sample_name) {
  tryCatch({
    cat(sprintf("Processing sample: %s\n", sample_name))
    obj_rna <- readRDS(sample_path)
    pro_name <- sample_name
    pro_name_d <- paste0(sample_name,"_dna")
    varbin_mtx_all_log2 <- readRDS(paste0("./objects/",pro_name_d,"_merged_integerCN_copykit_RNA_meta_new.rds"))
    obj_rna2 <- obj_rna
    gene_bin <- read.table("./pre_load_data/hg19_gene_binpos_map.tsv", header = 1)
    bin_chrpos <- varbin_mtx_all_log2@rowRanges %>% as.data.frame() %>% rownames_to_column() %>% 
      mutate(chr_pos = paste0("chr",seqnames,":",start, "-", end)) %>% dplyr::select(rowname, chr_pos)
    
    subclone_count <- obj_rna2@meta.data %>% dplyr::filter(subclone3 %in% (paste0("c", 2:30))) %>% pull(subclone3) %>% table() 
    
    #---use minimum 20 cell with RNA data as subclone cut-off---
    min_cell_num <- 20
    subclone_count <- subclone_count[subclone_count > min_cell_num] %>% names()
    combn_df <- combn(subclone_count, 2) %>% t() %>% as.data.frame() 
    combn_df
    
    dna_seg_rna_df <- NULL
    for (i in 1:nrow(combn_df)) {
      #i <- 1
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
      dna_con$rowname <- gsub("V","",dna_con$rowname)
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
      dna_rna_merged_cna_exp_target$cna_diff
      
      data <- dna_rna_merged_cna_exp_target %>% 
        dplyr::group_by(chr) %>% 
        mutate(group = cumsum(cna_diff != lag(cna_diff, default = dplyr::first(cna_diff)))) %>% ungroup()
      
      # Grouping by 'chr' and 'group' and then calculating the mean of 'diff'
      grouped_median <- data %>% dplyr::group_by(chr, group) %>%
        dplyr::summarise(median_diff = median(diff), first_cna_diff = dplyr::first(cna_diff))
      
      write.csv(grouped_median, 
                paste0("./metrics/genotype_seg_vs_sudo_bulk_expr_seg/",pro_name, "_seg_level_genotype_phenotype_statistic_count",super_clone1,"_vs_",super_clone2,".csv"))
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
    #dna_seg_rna_df2 <- read.csv(paste0("./metrics/",pro_name, "_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_", min_cell_num, "_Cells.csv"), row.names = 1)
    
    dna_seg_rna_df <- dna_seg_rna_df %>% mutate(matched = as.integer(matched), dna_but_no_rna = as.integer(dna_but_no_rna), 
                                                no_dna_but_rna = as.integer(no_dna_but_rna))
    sum_matched <- sum(dna_seg_rna_df$matched)
    sum_dna_but_no_rna <- sum(dna_seg_rna_df$dna_but_no_rna)
    sum_no_dna_but_rna <- sum(dna_seg_rna_df$no_dna_but_rna)
    
    sum_matched/(sum_matched + sum_dna_but_no_rna)
    
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
    cowplot::ggsave2(paste0("./figures/", pro_name, "_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_",min_cell_num,"_Cells_0.03_venn.pdf"), 
                     p2, width = 3, height = 2.5)
    
    
    
    
    cat(sprintf("Finished processing sample: %s\n", sample_name))
  }, error = function(e) {
    cat(sprintf("Error processing sample: %s - %s\n", sample_name, e$message))
  })
}

# Define paths and names for all samples
samples <- list(
  list(path = "objects/bcis66t_rna_final_merged_integerCN.rds", name = "bcis66t"),
  list(path = "objects/bcis106t_rna_final_merged_integerCN.rds", name = "bcis106t"),
  list(path = "objects/bcis74t_rna_final_merged_integerCN.rds", name = "bcis74t"),
  list(path = "objects/bcis51t_rna_final_merged_integerCN.rds", name = "bcis51t"),
  list(path = "objects/bcis70t_rna_final_merged_integerCN.rds", name = "bcis70t"),
  list(path = "objects/ecis57t_rna_final_merged_integerCN.rds", name = "ecis57t"),
  list(path = "objects/BCIS13T_Chip1_2_RNA.rds", name = "bcis13t"),
  list(path = "objects/BCIS28T_Chip1_2_RNA.rds", name = "bcis28t"),
  list(path = "objects/ECIS48T_Chip1_2_RNA.rds", name = "ecis48t"),
  list(path = "objects/ECIS25T_Chip1_2_3_RNA.rds", name = "ecis25t"),
  list(path = "objects/ECIS36T_Chip1_2_RNA.rds", name = "ecis36t"),
  list(path = "objects/ECIS44T_Chip1_2_3_4_5_RNA.rds", name = "ecis44t")
)

lapply(samples, function(sample) {
  process_sample(sample$path, sample$name)
})



process_sample <- function(sample_path, sample_name) {
  tryCatch({
    cat(sprintf("Processing sample: %s\n", sample_name))
    obj_rna <- readRDS(sample_path)
    pro_name <- sample_name
    pro_name_d <- paste0(sample_name,"_dna")
    varbin_mtx_all_log2 <- readRDS(paste0("./objects/",pro_name_d,"_merged_integerCN_copykit_RNA_meta_new.rds"))
    obj_rna2 <- obj_rna
    gene_bin <- read.table("./pre_load_data/hg19_gene_binpos_map.tsv", header = 1)
    bin_chrpos <- varbin_mtx_all_log2@rowRanges %>% as.data.frame() %>% rownames_to_column() %>% 
      mutate(chr_pos = paste0("chr",seqnames,":",start, "-", end)) %>% dplyr::select(rowname, chr_pos)
    #obj_rna2@meta.data$subclone3 <- obj_rna2@meta.data$subclones
    subclone_count <- obj_rna2@meta.data %>% dplyr::filter(subclone3 %in% (paste0("c", 2:30))) %>% pull(subclone3) %>% table()
    #---use minimum 20 cell with RNA data as subclone cut-off---
    min_cell_num <- 20
    subclone_count <- subclone_count[subclone_count > min_cell_num] %>% names()
    combn_df <- combn(subclone_count, 2) %>% t() %>% as.data.frame() 
    combn_df
    
    dna_seg_rna_df <- NULL
    for (i in 1:nrow(combn_df)) {
      #i <- 1
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
      dna_con$rowname <- gsub("V","",dna_con$rowname)
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
      dna_rna_merged_cna_exp_target$cna_diff
      
      data <- dna_rna_merged_cna_exp_target %>% 
        dplyr::group_by(chr) %>% 
        mutate(group = cumsum(cna_diff != lag(cna_diff, default = dplyr::first(cna_diff)))) %>% ungroup()
      
      # Grouping by 'chr' and 'group' and then calculating the mean of 'diff'
      grouped_median <- data %>% dplyr::group_by(chr, group) %>%
        dplyr::summarise(median_diff = median(diff), first_cna_diff = dplyr::first(cna_diff))
      
      write.csv(grouped_median, 
                paste0("./metrics/genotype_seg_vs_sudo_bulk_expr_seg/",pro_name, "_seg_level_genotype_phenotype_statistic_count",super_clone1,"_vs_",super_clone2,".csv"))
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
    #dna_seg_rna_df2 <- read.csv(paste0("./metrics/",pro_name, "_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_", min_cell_num, "_Cells.csv"), row.names = 1)
    
    dna_seg_rna_df <- dna_seg_rna_df %>% mutate(matched = as.integer(matched), dna_but_no_rna = as.integer(dna_but_no_rna), 
                                                no_dna_but_rna = as.integer(no_dna_but_rna))
    sum_matched <- sum(dna_seg_rna_df$matched)
    sum_dna_but_no_rna <- sum(dna_seg_rna_df$dna_but_no_rna)
    sum_no_dna_but_rna <- sum(dna_seg_rna_df$no_dna_but_rna)
    
    sum_matched/(sum_matched + sum_dna_but_no_rna)

    
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
    cowplot::ggsave2(paste0("./figures/", pro_name, "_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_",min_cell_num,"_Cells_0.03_venn.pdf"), 
                     p2, width = 3, height = 2.5)
    
    
    cat(sprintf("Finished processing sample: %s\n", sample_name))
  }, error = function(e) {
    cat(sprintf("Error processing sample: %s - %s\n", sample_name, e$message))
  })
}

# Define paths and names for all samples
samples <- list(
  list(path = "objects/bcis66t_rna_final_merged_integerCN.rds", name = "bcis66t"),
  list(path = "objects/bcis106t_rna_final_merged_integerCN.rds", name = "bcis106t"),
  list(path = "objects/bcis74t_rna_final_merged_integerCN.rds", name = "bcis74t"),
  list(path = "objects/bcis51t_rna_final_merged_integerCN.rds", name = "bcis51t"),
  list(path = "objects/bcis70t_rna_final_merged_integerCN.rds", name = "bcis70t"),
  list(path = "objects/ecis57t_rna_final_merged_integerCN.rds", name = "ecis57t"),
  list(path = "objects/BCIS13T_Chip1_2_RNA.rds", name = "bcis13t"),
  list(path = "objects/BCIS28T_Chip1_2_RNA.rds", name = "bcis28t"),
  list(path = "objects/ECIS48T_Chip1_2_RNA.rds", name = "ecis48t"),
  list(path = "objects/ECIS25T_Chip1_2_3_RNA.rds", name = "ecis25t"),
  list(path = "objects/ECIS36T_Chip1_2_RNA.rds", name = "ecis36t"),
  list(path = "objects/ECIS44T_Chip1_2_3_4_5_RNA.rds", name = "ecis44t"),
  list(path = "objects/WDR_231_RNA_downsampled_filtered.rds", name = "wdr_nanowell")
)

lapply(samples, function(sample) {
  process_sample(sample$path, sample$name)
})



#--- Stacked bar plot(Figure 5H)------#####
my_palette<- c("#8DA0CB","#FC8D62","#66C2A5")
min_cell_num <- 20
new_samples <- c("bcis51t","bcis66t","bcis106t","bcis70t","bcis74t","ecis57t","bcis13t", "bcis28t", "ecis25t", "ecis36t","wdr_nanowell")
sum_matched<- vector()
sum_dna_but_no_rna<- vector()
sum_no_dna_but_rna<- vector()
for (pro_name in new_samples){
  dna_seg_rna_df <- read.csv(paste0("./metrics/",pro_name, "_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_", min_cell_num, "_Cells.csv"), row.names = 1)
  dna_seg_rna_df <- dna_seg_rna_df %>% mutate(matched = as.integer(matched), dna_but_no_rna = as.integer(dna_but_no_rna), 
                                              no_dna_but_rna = as.integer(no_dna_but_rna))
  sum_matched <- rbind(sum_matched,sum(dna_seg_rna_df$matched))
  sum_dna_but_no_rna <- rbind(sum_dna_but_no_rna,sum(dna_seg_rna_df$dna_but_no_rna))
  sum_no_dna_but_rna <- rbind(sum_no_dna_but_rna,sum(dna_seg_rna_df$no_dna_but_rna))
  
  
}
stacker_bar_chart <- data.frame(new_samples,matched=sum_matched,CNA_diff =sum_dna_but_no_rna,RNA_diff=sum_no_dna_but_rna)
stacker_long <- stacker_bar_chart %>%
  pivot_longer(cols = -new_samples, names_to = "sum", values_to = "value")

stacker_percentage <- stacker_long %>%
  dplyr::group_by(new_samples) %>%
  dplyr::mutate(Total = sum(value), Proportion = value / Total) %>%
  dplyr::ungroup()


stacker_sorted <- stacker_percentage %>%
  dplyr::group_by(new_samples) %>%
  dplyr::mutate(matched_proportion = sum(Proportion[sum == "matched"])) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(matched_proportion))
stacker_sorted$sum <- factor(stacker_sorted$sum, levels = c("matched","CNA_diff","RNA_diff"))

new_samples <- c("bcis51t","bcis66t","bcis106t","bcis70t","bcis74t","ecis57t","bcis13t", "bcis28t", "ecis25t", "ecis36t", "ecis44t", "ecis48t","wdr_nanowell")
pnames<-        c("P11","P9","P10","P8","P7","P12","P5","P2","P1","P6","P4","P3","wdr_nanowell")
pname_mtx <- data.frame(pnames,sample = new_samples)
stacker_sorted <- left_join(stacker_sorted, pname_mtx ,by = c("new_samples"= "sample") )

p1 <- ggplot(stacker_sorted, aes(x = reorder(pnames, matched_proportion), y = log2(value+1), fill = sum)) +
  geom_bar(stat = "identity", position = "stack")+
  labs(
    x = "Samples",
    y = "log2"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line()
  )+scale_fill_manual(values=my_palette)

p2 <- ggplot(stacker_sorted, aes(x = reorder(pnames, matched_proportion), y = Proportion, fill = sum)) +
  geom_bar(stat = "identity", position = "stack") +  
  labs(
    x = "Samples",
    y = "Proportion"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line()
  ) +scale_fill_manual(values=my_palette)

p3 <- plot_grid(plotlist=list(p1,p2), ncol=1, align='v', rel_heights = c(1,2))
p3
cowplot::ggsave2(paste0("./figures/all_samples_genotype_seg_vs_sudo_bulk_expr_seg_Statistic_WithMorethan_",min_cell_num,"_Cells_0.03_venn.pdf"), 
                 p3, width = 4, height = 4)





#---find cancer cells(Figure 5A)----######
#load data
cs_mtx_all <- read.csv(paste0("./metrics/consensus_of_all_samples.csv"),row.names = 1)
clone_seg <- read.csv(paste0("./metrics/all_subclones_seg_info.csv"),row.names = 1)
mtx_srt <- read.csv(paste0("./metrics/cells_info_of_all_samples.csv"),row.names = 1)
rna.cluster.ids<-c("Fibro","LumHR","LumSec","MyoEpi","Peri","Tumor","VasEndo","Tcell","RNA_miss","Unknown1")
rna.cluster.col<-c("#D89000","#FF62BC","#00BD5F","#9590FF","#F8766D","#54278f","#BF4000", "#3969AC","gray90","#fff1d0")
names(rna.cluster.col) <- rna.cluster.ids
rui <- read.csv("./unique_subclone.txt",header = F) #need Rui
rui$V1 <- tolower(rui$V1)
rui_clones <- rui$V1[1:57]


set.seed(17)
clones2_col <- sample(unique(c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"), 
                               RColorBrewer::brewer.pal(12, "Paired"),rcartocolor::carto_pal(12, "Bold"),
                               rcartocolor::carto_pal(12, "Pastel"), rcartocolor::carto_pal(12, "Safe"),
                               viridis::viridis(22, option = "D"))), size = length(rui_clones))
names(clones2_col)<- unique(anno_mtx$new_clone_2)

mtx_srt_cancer <- mtx_srt %>% filter(new_clone %in% rui_clones )
mtx_srt_cancer_clone <- mtx_srt %>% filter(new_clone %in% rui_clones ) %>% dplyr::group_by(new_clone) %>% dplyr::summarise(sample = max(sample))
cs_srt_cancer <- left_join(mtx_srt_cancer,cs_mtx_all, by = c("new_clone" = "rowname")) %>% 
  select(-c(2:5)) %>% column_to_rownames(var ="sample")
cs_srt_cancer_clone <- left_join(mtx_srt_cancer_clone,cs_mtx_all, by = c("new_clone" = "rowname")) %>% column_to_rownames(var ="sample") %>% 
  select(-1)
anno_mtx <- mtx_srt_cancer %>% mutate(rna_clst0 = as.character(cell_type)) %>% 
  mutate(rna_clst = ifelse(is.na(rna_clst0), "RNA_miss",rna_clst0)) %>% 
  dplyr::select(c("new_clone", "rna_clst"))
anno_mtx$rna_clst <- gsub("ECs","VasEndo",anno_mtx$rna_clst)
anno_mtx$rna_clst <- gsub("CAFs","Fibro",anno_mtx$rna_clst)
anno_mtx$rna_clst <- gsub("2","",anno_mtx$rna_clst)
anno_mtx$rna_clst <- gsub("3","",anno_mtx$rna_clst)
anno_mtx$rna_clst <- gsub("4","",anno_mtx$rna_clst)
rownames(anno_mtx) <- NULL

rna_col <- rna.cluster.col[names(rna.cluster.col) %in% unique(anno_mtx$rna_clst)]
mtx_srt_cancer_clone <- mtx_srt %>% filter(new_clone %in% rui_clones ) %>% dplyr::group_by(new_clone) %>% 
  dplyr::summarise(sample = max(sample))

dist_matrix <- dist(cs_srt_cancer_clone, method = "euclidean")
hc <- hclust(dist_matrix, method = "complete")  
sample_order <- hc$order  
cluster_4 <- sample_order[1:3]
cluster_5 <- sample_order[4:21]
cluster_2 <- sample_order[22:25]
cluster_1 <- sample_order[26:31]
cluster_3 <- sample_order[32:40]
cluster_6 <- sample_order[41:length(sample_order)]
new_order <- c(cluster_1,cluster_2,cluster_3,cluster_4,cluster_5,cluster_6)
mtx_srt_cancer_clone <- mtx_srt_cancer_clone[new_order,]

mtx_srt_cancer <- mtx_srt %>% filter(new_clone %in% rui_clones )
mtx_srt_cancer$new_clone <- factor(mtx_srt_cancer$new_clone , levels = mtx_srt_cancer_clone$new_clone)
mtx_srt_cancer <- mtx_srt_cancer[order(mtx_srt_cancer$new_clone), ]
cs_srt_cancer <- left_join(mtx_srt_cancer,cs_mtx_all, by = c("new_clone" = "rowname")) %>% 
  select(-c(2:5)) %>% column_to_rownames(var ="sample")
#add patient info
mtx_srt_cancer_new <- mtx_srt_cancer %>% mutate(sample_name = gsub("_.*","",new_clone)) %>% left_join(pname_mtx, by = c("sample_name"= "sample")) %>% 
  mutate(new_clone_2 = paste0(pnames,"_",subclones3))

anno_mtx <- mtx_srt_cancer_new %>% dplyr::select(c("new_clone_2","pnames"))
rownames(anno_mtx) <- NULL
cs_srt_cancer <- data.frame(lapply(cs_srt_cancer, function(x) ifelse(x > 10, 10, x)))


#---heatmap of cancer_subclone_consensus
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color,
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

heat_col <- pals::ocean.balance(14)[9:14]
if(max(cs_srt_cancer)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(cs_srt_cancer)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(cs_srt_cancer)-6))
}
names(col_vec) <- 0:max(cs_srt_cancer)

ha_row2=rowAnnotation(df = anno_mtx, col = list(pnames = sample_platte,new_clone_2 =clones2_col),
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/cancer_subclone_consensus.pdf"), height = 20, width = 10)  
Heatmap(as.matrix(cs_srt_cancer), cluster_columns = FALSE, border = TRUE, cluster_rows = F, show_row_dend = F,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(cs_srt_cancer), " single cells"),
        use_raster = T, raster_quality = 5, col = col_vec,
        heatmap_legend_param = list(title = "Integer CNA",
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row2)
dev.off()
