#---find sporadic somatic CNA events---####

#load data
pro_name<- "ecis57t"
pro_name_d <-paste0(pro_name,"_dna")
varbin_mtx_ecis57t <- readRDS(paste0("./objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))

pro_name<- "bcis74t"
pro_name_d <-paste0(pro_name,"_dna")
varbin_mtx_bcis74t <- readRDS(paste0("objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))

pro_name<- "bcis51t"
pro_name_d <-paste0(pro_name,"_dna")
varbin_mtx_bcis51t <- readRDS(paste0("objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))

pro_name<- "bcis106t"
pro_name_d <-paste0(pro_name,"_dna")
varbin_mtx_bcis106t <- readRDS(paste0("objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))

pro_name<- "bcis70t"
pro_name_d <-paste0(pro_name,"_dna")
varbin_mtx_bcis70t <- readRDS(paste0("objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))

pro_name<- "bcis13t"
pro_name_d <-paste0(pro_name,"_dna")
varbin_mtx_bcis13t <- readRDS(paste0("objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))

pro_name<- "bcis28t"
pro_name_d <-paste0(pro_name,"_dna")
varbin_mtx_bcis28t <- readRDS(paste0("objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))

pro_name<- "ecis44t"
pro_name_d <-paste0(pro_name,"_dna")
varbin_mtx_ecis44t <- readRDS(paste0("objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))

pro_name<- "ecis48t"
pro_name_d <-paste0(pro_name,"_dna")
varbin_mtx_ecis48t <- readRDS(paste0("objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))

pro_name<- "bcis66t"
pro_name_d <-paste0(pro_name,"_dna")
varbin_mtx_bcis66t <- readRDS(paste0("objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))

pro_name<- "ecis25t"
pro_name_d <-paste0(pro_name,"_dna")
varbin_mtx_ecis25t <- readRDS(paste0("objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))

pro_name<- "ecis36t"
pro_name_d <-paste0(pro_name,"_dna")
varbin_mtx_ecis36t <- readRDS(paste0("objects/", pro_name_d, c("_merged_integerCN_copykit_RNA_meta_new.rds")))

# rna color palette
rna.cluster.ids<-c("Fibro","LumHR","LumSec","MyoEpi","Peri","Tumor","VasEndo","Tcell","RNA_miss")
rna.cluster.col<-c("#D89000","#FF62BC","#00BD5F","#9590FF","#F8766D","#54278f","#BF4000", "#3969AC","gray90")
names(rna.cluster.col) <- rna.cluster.ids

new_clone_color <- c("#004529","#088247","#7CC767",
                     "#223D6C","#5D90BA","#9ecae1",
                     "#8c510a","#bf812d","#D8D155",
                     "#b30000","#d73027","#fc9272",
                     "#7A142C","#E0367A","#df65b0",
                     "#4a1486","#6a51a3","#9e9ac8",
                     "#252525","#737373")

samples <- c("bcis51t","bcis66t","bcis106t","bcis70t","bcis74t","ecis57t","bcis13t", "bcis28t", "ecis25t", "ecis36t", "ecis44t", "ecis48t")
#prepare small subclone from each sample
cs_mtx_ecis57t <- as.data.frame(varbin_mtx_ecis57t@consensus)
ht_mtx_ecis57t <- as.data.frame(t(varbin_mtx_ecis57t@assays@data$integer))%>% rownames_to_column()
colnames(cs_mtx_ecis57t) <- paste0("ecis57t_",colnames(cs_mtx_ecis57t))
mtx_srt_ecis57t <- as.data.frame(varbin_mtx_ecis57t@colData)%>%
  select(sample,experiment,subclones3,cell_type)%>%
  mutate(new_clone = paste0("ecis57t_",subclones3))

cs_mtx_bcis74t <- as.data.frame(varbin_mtx_bcis74t@consensus)
ht_mtx_bcis74t <- as.data.frame(t(varbin_mtx_bcis74t@assays@data$integer))%>% rownames_to_column()
colnames(cs_mtx_bcis74t) <- paste0("bcis74t_",colnames(cs_mtx_bcis74t))
mtx_srt_bcis74t <- as.data.frame(varbin_mtx_bcis74t@colData)%>%
  select(sample,experiment,subclones3,cell_type)%>%
  mutate(new_clone = paste0("bcis74t_",subclones3))


cs_mtx_bcis51t <- as.data.frame(varbin_mtx_bcis51t@consensus)
ht_mtx_bcis51t <- as.data.frame(t(varbin_mtx_bcis51t@assays@data$integer))%>% rownames_to_column()
colnames(cs_mtx_bcis51t) <- paste0("bcis51t_",colnames(cs_mtx_bcis51t))
mtx_srt_bcis51t <- as.data.frame(varbin_mtx_bcis51t@colData)%>%
  select(sample,experiment,subclones3,cell_type)%>%
  mutate(new_clone = paste0("bcis51t_",subclones3))

cs_mtx_bcis106t <- as.data.frame(varbin_mtx_bcis106t@consensus)
ht_mtx_bcis106t <- as.data.frame(t(varbin_mtx_bcis106t@assays@data$integer))%>% rownames_to_column()
colnames(cs_mtx_bcis106t) <- paste0("bcis106t_",colnames(cs_mtx_bcis106t))
mtx_srt_bcis106t <- as.data.frame(varbin_mtx_bcis106t@colData)%>%
  select(sample,experiment,subclones3,cell_type)%>%
  mutate(new_clone = paste0("bcis106t_",subclones3))

cs_mtx_bcis70t <- as.data.frame(varbin_mtx_bcis70t@consensus)
ht_mtx_bcis70t <- as.data.frame(t(varbin_mtx_bcis70t@assays@data$integer))%>% rownames_to_column()
colnames(cs_mtx_bcis70t) <- paste0("bcis70t_",colnames(cs_mtx_bcis70t))
mtx_srt_bcis70t <- as.data.frame(varbin_mtx_bcis70t@colData)%>%
  select(sample,experiment,subclones3,cell_type)%>%
  mutate(new_clone = paste0("bcis70t_",subclones3))

cs_mtx_bcis28t <- as.data.frame(varbin_mtx_bcis28t@consensus)
ht_mtx_bcis28t <- as.data.frame(t(varbin_mtx_bcis28t@assays@data$integer))%>% rownames_to_column()
colnames(cs_mtx_bcis28t) <- paste0("bcis28t_",colnames(cs_mtx_bcis28t))
mtx_srt_bcis28t <- as.data.frame(varbin_mtx_bcis28t@colData)%>%
  select(sample,experiment,subclones3,cell_type)%>%
  mutate(new_clone = paste0("bcis28t_",subclones3))

cs_mtx_bcis13t <- as.data.frame(varbin_mtx_bcis13t@consensus)
ht_mtx_bcis13t <- as.data.frame(t(varbin_mtx_bcis13t@assays@data$integer))%>% rownames_to_column()
colnames(cs_mtx_bcis13t) <- paste0("bcis13t_",colnames(cs_mtx_bcis13t))
mtx_srt_bcis13t <- as.data.frame(varbin_mtx_bcis13t@colData)%>%
  select(sample,experiment,subclones3,cell_type)%>%
  mutate(new_clone = paste0("bcis13t_",subclones3))

cs_mtx_ecis44t <- as.data.frame(varbin_mtx_ecis44t@consensus)
ht_mtx_ecis44t <- as.data.frame(t(varbin_mtx_ecis44t@assays@data$integer))%>% rownames_to_column()
colnames(cs_mtx_ecis44t) <- paste0("ecis44t_",colnames(cs_mtx_ecis44t))
mtx_srt_ecis44t <- as.data.frame(varbin_mtx_ecis44t@colData)%>%
  select(sample,experiment,subclones3,cell_type)%>%
  mutate(new_clone = paste0("ecis44t_",subclones3))

cs_mtx_ecis48t <- as.data.frame(varbin_mtx_ecis48t@consensus)
ht_mtx_ecis48t <- as.data.frame(t(varbin_mtx_ecis48t@assays@data$integer))%>% rownames_to_column()
colnames(cs_mtx_ecis48t) <- paste0("ecis48t_",colnames(cs_mtx_ecis48t))
mtx_srt_ecis48t <- as.data.frame(varbin_mtx_ecis48t@colData)%>%
  select(sample,experiment,subclones3,cell_type)%>%
  mutate(new_clone = paste0("ecis48t_",subclones3))

cs_mtx_bcis66t <- as.data.frame(varbin_mtx_bcis66t@consensus)
ht_mtx_bcis66t <- as.data.frame(t(varbin_mtx_bcis66t@assays@data$integer))%>% rownames_to_column()
colnames(cs_mtx_bcis66t) <- paste0("bcis66t_",colnames(cs_mtx_bcis66t))
mtx_srt_bcis66t <- as.data.frame(varbin_mtx_bcis66t@colData)%>%
  select(sample,experiment,subclones3,cell_type)%>%
  mutate(new_clone = paste0("bcis66t_",subclones3))

cs_mtx_ecis25t <- as.data.frame(varbin_mtx_ecis25t@consensus)
ht_mtx_ecis25t <- as.data.frame(t(varbin_mtx_ecis25t@assays@data$integer))%>% rownames_to_column()
colnames(cs_mtx_ecis25t) <- paste0("ecis25t_",colnames(cs_mtx_ecis25t))
mtx_srt_ecis25t <- as.data.frame(varbin_mtx_ecis25t@colData)%>%
  select(sample,experiment,subclones3,cell_type)%>%
  mutate(new_clone = paste0("ecis25t_",subclones3))

cs_mtx_ecis36t <- as.data.frame(varbin_mtx_ecis36t@consensus)
ht_mtx_ecis36t <- as.data.frame(t(varbin_mtx_ecis36t@assays@data$integer))%>% rownames_to_column()
colnames(cs_mtx_ecis36t) <- paste0("ecis36t_",colnames(cs_mtx_ecis36t))
mtx_srt_ecis36t <- as.data.frame(varbin_mtx_ecis36t@colData)%>%
  select(sample,experiment,subclones3,cell_type)%>%
  mutate(new_clone = paste0("ecis36t_",subclones3))



cs_mtx_ecis57t <-  cbind(varbin_mtx_ecis57t@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_ecis57t@consensus) %>%
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything()) %>% mutate(sample_id = paste0("ecis57t_",sample_id)) 
cs_mtx_bcis74t <-  cbind(varbin_mtx_bcis74t@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_bcis74t@consensus) %>%
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())%>% mutate(sample_id = paste0("bcis74t_",sample_id)) 
cs_mtx_bcis51t <-  cbind(varbin_mtx_bcis51t@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_bcis51t@consensus) %>%
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())%>% mutate(sample_id = paste0("bcis51t_",sample_id)) 
cs_mtx_bcis106t <-  cbind(varbin_mtx_bcis106t@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_bcis106t@consensus) %>%
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())%>% mutate(sample_id = paste0("bcis106t_",sample_id)) 
cs_mtx_bcis13t <-  cbind(varbin_mtx_bcis13t@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_bcis13t@consensus) %>%
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())%>% mutate(sample_id = paste0("bcis13t_",sample_id)) 
cs_mtx_bcis70t <-  cbind(varbin_mtx_bcis70t@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_bcis70t@consensus) %>%
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())%>% mutate(sample_id = paste0("bcis70t_",sample_id)) 
cs_mtx_bcis28t <-  cbind(varbin_mtx_bcis28t@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_bcis28t@consensus) %>%
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())%>% mutate(sample_id = paste0("bcis28t_",sample_id)) 
cs_mtx_ecis44t <-  cbind(varbin_mtx_ecis44t@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_ecis44t@consensus) %>%
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())%>% mutate(sample_id = paste0("ecis44t_",sample_id)) 
cs_mtx_ecis48t <-  cbind(varbin_mtx_ecis48t@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_ecis48t@consensus) %>%
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())%>% mutate(sample_id = paste0("ecis48t_",sample_id)) 
cs_mtx_bcis66t <-  cbind(varbin_mtx_bcis66t@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_bcis66t@consensus) %>%
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())%>% mutate(sample_id = paste0("bcis66t_",sample_id)) 
cs_mtx_ecis25t <-  cbind(varbin_mtx_ecis25t@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_ecis25t@consensus) %>%
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())%>% mutate(sample_id = paste0("ecis25t_",sample_id)) 
cs_mtx_ecis36t <-  cbind(varbin_mtx_ecis36t@rowRanges %>% dplyr::as_tibble() %>% dplyr::select(seqnames, start, end), varbin_mtx_ecis36t@consensus) %>%
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())%>% mutate(sample_id = paste0("ecis36t_",sample_id)) 

#consensus for each subclone
cs_mtx <- rbind(cs_mtx_ecis57t,cs_mtx_bcis74t,cs_mtx_bcis51t,cs_mtx_bcis106t,cs_mtx_bcis70t,cs_mtx_bcis13t,cs_mtx_bcis28t,cs_mtx_ecis44t,cs_mtx_ecis48t,cs_mtx_bcis66t,cs_mtx_ecis25t,cs_mtx_ecis36t) #%>% t() %>% as.data.frame() %>% rownames_to_column()

cs_mtx <- cs_mtx %>%
  dplyr::group_by(sample_id,chrom) %>%
  dplyr::mutate(
    CN = case_when(
      row_number() %in% 1:10 ~ CN[11],  # For the first 10 rows, use the CN value of the 11th row
      row_number() %in% (n() - 9):n() ~ CN[n() - 10],  # For the last 10 rows, use the CN value of the 11th row from the bottom
      TRUE ~ CN  # For all other rows, keep the original CN value
    )
  )

cs_mtx_new <- cs_mtx %>%
  spread(key = sample_id, value = CN) %>% select(-start, -end)

cs_mtx_new <-cs_mtx_new[,-1]
ht_mtx <- rbind(ht_mtx_ecis57t,ht_mtx_bcis74t,ht_mtx_bcis51t,ht_mtx_bcis106t,ht_mtx_bcis70t,ht_mtx_bcis13t,ht_mtx_bcis28t,ht_mtx_ecis44t,ht_mtx_ecis48t,ht_mtx_bcis66t,ht_mtx_ecis25t,ht_mtx_ecis36t)
mtx_srt <- rbind(mtx_srt_ecis57t,mtx_srt_bcis74t,mtx_srt_bcis51t,mtx_srt_bcis106t,mtx_srt_bcis70t,mtx_srt_bcis13t,mtx_srt_bcis28t,mtx_srt_ecis44t,mtx_srt_ecis48t,mtx_srt_bcis66t,mtx_srt_ecis25t,mtx_srt_ecis36t)

write.csv(mtx_srt, paste0("./metrics/cells_info_of_all_samples.csv"))
#mtx_srt <- read.csv(paste0("./metrics/cells_info_of_all_samples.csv"),row.names = 1)

cs_mtx_all <- cs_mtx_new %>% t() %>% as.data.frame() %>% rownames_to_column()
cs_mtx_all
write.csv(cs_mtx_all, paste0("./metrics/consensus_of_all_samples.csv"))
#cs_mtx_all <- read.csv(paste0("./metrics/consensus_of_all_samples.csv"),row.names = 1)

clone_seg_num <- vector()
clone_seg_name <- vector()
clone_seg_binsize <- vector()
for(i in colnames(cs_mtx_new)){
  group <- 1
  group_name <- paste0("group_",i)
  for (j in 1:(nrow(cs_mtx_new) - 1)) {
    #j <- 1
    cs_mtx_new[[group_name]][j] <- group
    if (cs_mtx_new[[i]][j+1] == cs_mtx_new[[i]][j]) {
      cs_mtx_new[[group_name]][j+1] <- group
    } else {
      group <- group + 1
      cs_mtx_new[[group_name]][j+1] <- group
    }
  }
  tmp <- cs_mtx_new %>% filter(cs_mtx_new[[i]] !=2) 
  clone_seg_num <- append(clone_seg_num,length(unique(tmp[[group_name]])))
  clone_seg_binsize <-append(clone_seg_binsize,length(tmp[[group_name]]))
  clone_seg_name <- append(clone_seg_name,i)
}
clone_seg <- data.frame(clone_seg_name,clone_seg_num,clone_seg_binsize)
write.csv(clone_seg, paste0("./metrics/all_subclones_seg_info.csv"))
clone_seg <- read.csv(paste0("./metrics/all_subclones_seg_info.csv"),row.names = 1)

# filter sporadic clusters
cell_num <- 40
cutoff <- 5
clone_cutoff <- mtx_srt %>% dplyr::group_by(new_clone) %>% dplyr::summarise(count = n()) %>% filter(count < cell_num)
clone_seg_cutoff <- clone_seg %>% filter(clone_seg_num <=cutoff & clone_seg_num>0 & clone_seg_binsize <= 1000)
clone_cutoff <- clone_cutoff %>% filter(new_clone %in% clone_seg_cutoff$clone_seg_name & !new_clone %in% c("ecis44t_c8","bcis74t_c13"))
mtx_srt_filter <- mtx_srt %>% filter(new_clone %in% clone_cutoff$new_clone)
mtx_srt_new <- left_join(mtx_srt_filter,cs_mtx_all, by = c("new_clone" = "rowname")) %>% select(-c(2:5)) %>% column_to_rownames(var ="sample")

clone_cutoff
write.csv(clone_cutoff, paste0("./metrics/names_of_small_event_samples.csv"))
#clone_cutoff <- read.csv(paste0("./metrics/names_of_small_event_samples.csv"),row.names = 1)

clone_seg_cutoff_eventnum <- clone_seg_cutoff %>% filter(clone_seg_name %in%  clone_cutoff$new_clone) 
clone_seg_cutoff_eventnum <- sum(clone_seg_cutoff_eventnum$clone_seg_num)

mtx_srt_filter$cell_type <- gsub("ECs","VasEndo",mtx_srt_filter$cell_type)
mtx_srt_filter$cell_type <- gsub("CAFs","Fibro",mtx_srt_filter$cell_type)
mtx_srt_filter$cell_type <- gsub("2","",mtx_srt_filter$cell_type)
mtx_srt_filter$cell_type <- gsub("3","",mtx_srt_filter$cell_type)

write.csv(mtx_srt_filter, paste0("./metrics/cells_info_of_small_event_samples.csv"))
mtx_srt_filter <- read.csv(paste0("./metrics/cells_info_of_small_event_samples.csv"),row.names = 1)



# sporadic_cna_subclones_with_cell_num_40_and_seg_num_5_subclone_consensus_hcluster_heatmap
anno_mtx <- mtx_srt_filter %>% mutate(rna_clst0 = as.character(cell_type)) %>% 
  mutate(rna_clst = ifelse(is.na(rna_clst0), "RNA_miss",rna_clst0)) %>% 
  dplyr::select(c("new_clone", "rna_clst"))

rownames(anno_mtx) <- NULL
rna_col <- rna.cluster.col[names(rna.cluster.col) %in% unique(anno_mtx$rna_clst)]
clst_col <- new_pal[1:length(unique(anno_mtx$new_clone))]
names(clst_col) <- unique(anno_mtx$new_clone)

ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color,
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")
heat_col <- pals::ocean.balance(14)[9:14]
if(max(mtx_srt_new)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(mtx_srt_new)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(mtx_srt_new)-6))
}
names(col_vec) <- 0:max(mtx_srt_new)
ha_row2=rowAnnotation(df = anno_mtx, col = list(new_clone=clst_col, rna_clst=rna_col),
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/small_subclones_with_cell_num_40_and_seg_num_5_subclone_consensus_hcluster.pdf"), height = 10, width = 10)  
Heatmap(as.matrix(mtx_srt_new), cluster_columns = FALSE, border = TRUE, cluster_rows = TRUE, show_row_dend = F,
        #row_split = anno_mtx$new_clone,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(mtx_srt_new), " single cells"),
        column_title = paste0("with_cell_num_",cell_num, "_with_seg_num_", cutoff),
        use_raster = T, raster_quality = 5, col = col_vec,
        heatmap_legend_param = list(title = "Integer CNA",
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row2)
dev.off()


#---seperate epi and non-epi---#####
# load data 
cs_mtx_all <- read.csv(paste0("./metrics/consensus_of_all_samples.csv"),row.names = 1)
clone_seg <- read.csv(paste0("./metrics/all_subclones_seg_info.csv"),row.names = 1)
mtx_srt <- read.csv(paste0("./metrics/cells_info_of_all_samples.csv"),row.names = 1)
clone_cutoff <- read.csv(paste0("./metrics/names_of_small_event_samples.csv"),row.names = 1)
mtx_srt_filter <- read.csv(paste0("./metrics/cells_info_of_small_event_samples.csv"),row.names = 1)
rna.cluster.ids<-c("Fibro","LumHR","LumSec","MyoEpi","Peri","Tumor","VasEndo","Tcell","RNA_miss","Unknown1")
rna.cluster.col<-c("#D89000","#FF62BC","#00BD5F","#9590FF","#F8766D","#54278f","#BF4000", "#3969AC","gray90","#fff1d0")
names(rna.cluster.col) <- rna.cluster.ids



# small_subclones_with_cell_num_40_and_seg_num_5_subclone_consensus
mtx_srt_filter_new <- mtx_srt_filter %>% mutate(if_epi = case_when(cell_type %in% c("LumSec", "LumHR", "MyoEpi") ~ "epithelial_cell",
                                                                   is.na(cell_type) ~ "RNA_miss",
                                                                   cell_type == "Tumor" ~ "Tumor",
                                                                   TRUE ~ "non_epithelial_cell"))
anno_mtx <- mtx_srt_filter_new %>% mutate(rna_clst0 = as.character(cell_type)) %>% 
  mutate(rna_clst = ifelse(is.na(rna_clst0), "RNA_miss",rna_clst0)) %>% 
  dplyr::select(c("new_clone", "rna_clst","if_epi"))

rownames(anno_mtx) <- NULL


rna_col <- rna.cluster.col[names(rna.cluster.col) %in% unique(anno_mtx$rna_clst)]
clst_col <- new_pal[1:length(unique(anno_mtx$new_clone))]
names(clst_col) <- unique(anno_mtx$new_clone)


ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color,
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

heat_col <- pals::ocean.balance(14)[9:14]
if(max(mtx_srt_new)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(mtx_srt_new)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(mtx_srt_new)-6))
}
names(col_vec) <- 0:max(mtx_srt_new)
platte_2 <- c("#879e82","#c39f72","grey90","#54278f")
names(platte_2) <- c("epithelial_cell","non_epithelial_cell","RNA_miss","Tumor")

ha_row2=rowAnnotation(df = anno_mtx, col = list(new_clone=clst_col, rna_clst=rna_col,if_epi =platte_2),
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/small_subclones_with_cell_num_40_and_seg_num_5_subclone_consensus.pdf"), height = 10, width = 10)  
Heatmap(as.matrix(mtx_srt_new), cluster_columns = FALSE, border = TRUE, cluster_rows = F, show_row_dend = F,
        row_split = anno_mtx$new_clone,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(mtx_srt_filter_new), " single cells"),
        column_title = paste0("with_cell_num_",cell_num, "_with_seg_num_", cutoff),
        use_raster = T, raster_quality = 5, col = col_vec,
        heatmap_legend_param = list(title = "Integer CNA",
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row2)
dev.off()


#---Heatmap of single cell integer copy number profiles for non-cancer epithelial cells and stromal cells---#
epi_clones <- c("ecis44t_c2","ecis44t_c3","ecis44t_c4","ecis44t_c5","ecis44t_c6","ecis44t_c7","bcis28t_c2","ecis36t_c15",
                "ecis48t_c2","ecis48t_c3","ecis48t_c6","ecis48t_c7","bcis70t_c2")
mtx_srt_filter_new <- mtx_srt_filter %>% mutate(if_epi = case_when(cell_type %in% c("LumSec", "LumHR", "MyoEpi") ~ "epithelial_cell",
                                                                   is.na(cell_type) ~ "RNA_miss",
                                                                   cell_type == "Tumor" ~ "Tumor",
                                                                   TRUE ~ "non_epithelial_cell")) %>% filter(if_epi != "Tumor") %>%
  mutate(epi_clone = ifelse(new_clone %in% epi_clones,"epi","non_epi"))


# match patient names to sample names
new_samples <- c("bcis51t","bcis66t","bcis106t","bcis70t","bcis74t","ecis57t","bcis13t", "bcis28t", "ecis25t", "ecis36t", "ecis44t", "ecis48t")
pnames<-        c("P11","P9","P10","P8","P7","P12","P5","P2","P1","P6","P4","P3")
pname_mtx <- data.frame(pnames,sample = new_samples)
pname_mtx$pnames <- factor(pname_mtx$pnames, levels = paste0("P",1:12))


mtx_srt_filter_new_if_epi <- mtx_srt_filter_new %>% mutate(cell_type_new =paste0(new_clone,"_",cell_type),sample_name = gsub("_.*","",new_clone)) %>%
  left_join(pname_mtx, by = c("sample_name"="sample")) %>% mutate(new_clone_names = paste0(pnames,"_",subclones3))
mtx_srt_filter_new_group <- mtx_srt_filter_new %>% dplyr::group_by(new_clone,cell_type) %>% dplyr::summarise(count = n()) %>% 
  filter(count >= 3) %>% mutate(cell_type_new =paste0(new_clone,"_",cell_type))
mtx_srt_filter_new_if_epi <- mtx_srt_filter_new_if_epi %>% filter(cell_type_new %in% mtx_srt_filter_new_group$cell_type_new) %>% 
  filter(!new_clone %in% c("bcis51t_c2","ecis48t_c2","ecis48t_c4","ecis57t_c3","bcis28t_c3"))



test <- left_join(mtx_srt_filter_new_if_epi,cs_mtx_all, by = c("new_clone" = "rowname")) %>% 
  select(-c(2:11)) %>% column_to_rownames(var ="sample")
anno_mtx <- mtx_srt_filter_new_if_epi %>% mutate(rna_clst0 = as.character(cell_type)) %>% 
  mutate(rna_clst = ifelse(is.na(rna_clst0), "RNA_miss",rna_clst0)) %>% 
  dplyr::select(c("new_clone_names", "rna_clst","epi_clone","if_epi"))
rownames(anno_mtx) <- NULL

rna_col <- rna.cluster.col[names(rna.cluster.col) %in% unique(anno_mtx$rna_clst)]
clst_col <- new_pal[1:length(unique(anno_mtx$new_clone_names))]
names(clst_col) <- unique(anno_mtx$new_clone_names)

anno_mtx$epi_clone <- factor(anno_mtx$epi_clone, levels = )

platte_1 <- c("#008585","#e5c185")
names(platte_1) <- c("epi","non_epi")
platte_2 <- c("#879e82","#c39f72","grey90")
names(platte_2) <- c("epithelial_cell","non_epithelial_cell","RNA_miss")
#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color,
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

heat_col <- pals::ocean.balance(14)[9:14]
if(max(test)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(test)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(test)-6))
}
names(col_vec) <- 0:max(test)

ha_row2=rowAnnotation(df = anno_mtx, col = list(new_clone_names=clst_col, rna_clst=rna_col, epi_clone=platte_1,if_epi =platte_2),
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))


pdf(paste0("./figures/small_subclones_with_cell_num_40_and_seg_num_5_subclone_consensus_if_epi.pdf"), height = 10, width = 10)  
Heatmap(as.matrix(test), cluster_columns = FALSE, border = TRUE, cluster_rows = T, show_row_dend = F,
        row_split = anno_mtx$epi_clone,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(test), " single cells"),
        #column_title = "non epi cells",
        use_raster = T, raster_quality = 5, col = col_vec,
        heatmap_legend_param = list(title = "Integer CNA",
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row2)
dev.off()

#---plot figure4G---#####
dist_matrix <- dist(test, method = "euclidean")  
hc <- hclust(dist_matrix, method = "complete")  
sample_order <- hc$order  

mtx_srt_filter_new_if_epi$new_clone <- factor(mtx_srt_filter_new_if_epi$new_clone, 
                                              levels = c("ecis44t_c6","ecis44t_c7","ecis48t_c7","ecis36t_c15","bcis28t_c2","ecis44t_c3","ecis44t_c2","ecis48t_c3","bcis70t_c2","ecis44t_c4","ecis48t_c6","ecis44t_c5","bcis13t_c3","bcis106t_c2","ecis57t_c2","bcis13t_c2"))


sorted_indices <- unlist(
  lapply(levels(mtx_srt_filter_new_if_epi$new_clone), function(level) {
    group_samples <- which(mtx_srt_filter_new_if_epi$new_clone == level)  
    ordered_group <- group_samples[order(match(group_samples, sample_order))]  
    return(ordered_group)
  })
)

test_reordered <- test[sorted_indices, ]

mtx_srt_filter_new_if_epi_test <- mtx_srt_filter_new_if_epi[sorted_indices,]
anno_mtx <- mtx_srt_filter_new_if_epi_test %>% mutate(rna_clst0 = as.character(cell_type)) %>% 
  mutate(rna_clst = ifelse(is.na(rna_clst0), "RNA_miss",rna_clst0)) %>% 
  dplyr::select(c("new_clone_names", "rna_clst","epi_clone","if_epi"))
rownames(anno_mtx) <- NULL

rna_col <- rna.cluster.col[names(rna.cluster.col) %in% unique(anno_mtx$rna_clst)]
clst_col <- new_pal[1:length(unique(anno_mtx$new_clone_names))]
names(clst_col) <- unique(anno_mtx$new_clone_names)

anno_mtx$epi_clone <- factor(anno_mtx$epi_clone, levels = )

platte_1 <- c("#008585","#e5c185")
names(platte_1) <- c("epi","non_epi")
platte_2 <- c("#879e82","#c39f72","grey90")
names(platte_2) <- c("epithelial_cell","non_epithelial_cell","RNA_miss")
#-----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color,
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

heat_col <- pals::ocean.balance(14)[9:14]
if(max(test)<=8){
  col_vec <- c("#181C43","#3787BA","#F2F2F2FF", heat_col[1:(max(test)-2)])
}else{
  col_vec <- c("#181C43","#3787BA","white", heat_col, rep("#3C0912", times= max(test)-6))
}
names(col_vec) <- 0:max(test)

ha_row2=rowAnnotation(df = anno_mtx, col = list(new_clone_names=clst_col, rna_clst=rna_col, epi_clone=platte_1,if_epi =platte_2),
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/small_subclones_with_cell_num_40_and_seg_num_5_subclone_consensus_if_epi_by_order.pdf"), height = 10, width = 10)  
Heatmap(as.matrix(test_reordered), cluster_columns = FALSE, border = TRUE, cluster_rows = F, show_row_dend = F,
        row_split = anno_mtx$epi_clone,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(test), " single cells"),
        #column_title = "non epi cells",
        use_raster = T, raster_quality = 5, col = col_vec,
        heatmap_legend_param = list(title = "Integer CNA",
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row2)
dev.off()


#---freq_of_small_clones_among_their_samples(plot figure4I)---#####
mtx_srt_cell_num <-  mtx_srt %>% dplyr::group_by(new_clone) %>% dplyr::summarise(count = n()) %>% mutate(sample = gsub("_.*", "",new_clone))%>% left_join(pname_mtx, by = c("sample"="sample")) %>% mutate(panme= paste0(pnames,"_",gsub("^[^_]*_", "",new_clone)) )
sample_sum <- mtx_srt_cell_num %>% dplyr::group_by(sample) %>% dplyr::summarise(count = sum(count)) 
mtx_cell_all <- left_join(mtx_srt_cell_num,sample_sum, by ="sample") %>% filter(new_clone %in% mtx_srt_filter_new_if_epi$new_clone) %>% mutate(freq = count.x/count.y) 
mtx_cell_all <-arrange(mtx_cell_all, desc(freq))
mtx_cell_all$panme <- factor(mtx_cell_all$panme,levels = mtx_cell_all$panme )
p1<- ggplot(mtx_cell_all, aes(x = panme, y = freq,fill = panme)) +
  geom_bar(stat = "identity") + scale_fill_manual(values =clst_col) + theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1))

p1
cowplot::ggsave2(paste0("./figures/freq_of_small_clones_among_their_samples.pdf"), p1, width = 8, height = 5)





#---pie plot for cell types(figure 4H)-------####
cell_type_mtx <- mtx_srt_filter_new_if_epi %>% filter(cell_type != is.na(cell_type)) %>% dplyr::group_by(cell_type) %>% 
  dplyr::summarise(count = n())
cell_type_mtx$percentage <- round(cell_type_mtx$count / sum(cell_type_mtx$count) * 100, 1)
cell_type_mtx$label <- paste0(cell_type_mtx$percentage, "%")

p1 <- ggplot(cell_type_mtx, aes(x = "", y = count, fill = cell_type)) +  
  geom_bar(stat = "identity", width = 1) +               
  coord_polar(theta = "y") +                             
  theme_void() +                                         
  labs(fill = "cell_type")  +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))+
  scale_fill_manual(values =rna_col)
p1

cowplot::ggsave2(paste0("./figures/small_events_cell_type_pie_plot.pdf"), p1, width = 6, height = 4)

event_mtx <- data.frame(
  event = c("chr1_amplification","chrx_deletion","chrx_amplification","others"),
  count = c(5,5,2,22)
)
event_mtx$percentage <- round(event_mtx$count / sum(event_mtx$count) * 100, 1)
event_mtx$label <- paste0(event_mtx$percentage, "%")
p1 <- ggplot(event_mtx, aes(x = "", y = count, fill = event)) +  
  geom_bar(stat = "identity", width = 1) +               
  coord_polar(theta = "y") +                             
  theme_void() +                                         
  labs(fill = "event")  +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))
p1
cowplot::ggsave2(paste0("./figures/sp_events_pie_plot.pdf"), p1, width = 6, height = 4)


chr1_amp <- c("bcis28t_c2","ecis36t_c15","ecis44t_c6","ecis44t_c7","ecis48t_c7")
chrx_amp <- c("bcis13t_c2")
chrx_del <- c("bcis106t_c2","bcis13t_c3")

sp_event_subclones_mtx <- mtx_srt_filter_new_if_epi %>% filter(cell_type != is.na(cell_type) & new_clone %in% chrx_del) %>% 
  dplyr::group_by(cell_type) %>% dplyr::summarise(count = n())
sp_event_subclones_mtx$percentage <- round(sp_event_subclones_mtx$count / sum(sp_event_subclones_mtx$count) * 100, 1)
sp_event_subclones_mtx$label <- paste0(sp_event_subclones_mtx$percentage, "%")

p1 <- ggplot(sp_event_subclones_mtx, aes(x = "", y = percentage, fill = cell_type)) +  
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "cell_type_in_special_event_chrx_del") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  scale_fill_manual(values = rna.cluster.col) + theme_bw() + 
  theme(panel.grid = element_blank())
p1
cowplot::ggsave2(paste0("./figures/chrx_del_cell_type_stack_plot.pdf"), p1, width = 4, height = 4)



sp_event_subclones_mtx <- mtx_srt_filter_new_if_epi %>% filter(cell_type != is.na(cell_type) & new_clone %in% chr1_amp) %>% 
  dplyr::group_by(cell_type) %>% dplyr::summarise(count = n())
sp_event_subclones_mtx$percentage <- round(sp_event_subclones_mtx$count / sum(sp_event_subclones_mtx$count) * 100, 1)
sp_event_subclones_mtx$label <- paste0(sp_event_subclones_mtx$percentage, "%")

p1 <- ggplot(sp_event_subclones_mtx, aes(x = "", y = count, fill = cell_type)) +  
  geom_bar(stat = "identity", width = 1) +               
  coord_polar(theta = "y") +                             
  theme_void() +                                         
  labs(fill = "cell_type",
       title = "cell_type_in_special_event_chr1_amp")  +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))+
  scale_fill_manual(values =rna_col)
p1

cowplot::ggsave2(paste0("./figures/chr1_amp_cell_type_stack_plot.pdf"), p1, width = 4, height = 4)

sp_event_subclones_mtx <- mtx_srt_filter_new_if_epi %>% filter(cell_type != is.na(cell_type) & new_clone %in% chrx_amp) %>% 
  dplyr::group_by(cell_type) %>% dplyr::summarise(count = n())
sp_event_subclones_mtx$percentage <- round(sp_event_subclones_mtx$count / sum(sp_event_subclones_mtx$count) * 100, 1)
sp_event_subclones_mtx$label <- paste0(sp_event_subclones_mtx$percentage, "%")

p1 <- ggplot(sp_event_subclones_mtx, aes(x = "", y = count, fill = cell_type)) +  
  geom_bar(stat = "identity", width = 1) +               
  coord_polar(theta = "y") +                             
  theme_void() +                                         
  labs(fill = "cell_type",
       title = "cell_type_in_special_event_chrx_amp")  +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))+
  scale_fill_manual(values =rna_col)
p1

cowplot::ggsave2(paste0("./figures/chrx_amp_cell_type_stack_plot.pdf"), p1, width = 4, height = 4)




# epi and non epi
cell_type_mtx <- mtx_srt_filter_new_if_epi %>% filter(epi_clone == "epi" & cell_type != is.na(cell_type)) %>% 
  dplyr::group_by(cell_type) %>%   dplyr::summarise(count = n())
cell_type_mtx$percentage <- round(cell_type_mtx$count / sum(cell_type_mtx$count) * 100, 1)
cell_type_mtx$label <- paste0(cell_type_mtx$percentage, "%")

p1 <- ggplot(cell_type_mtx, aes(x = "", y = percentage, fill = cell_type)) +  
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "epi") +  
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  scale_fill_manual(values = rna.cluster.col) + theme_bw() + 
  theme(panel.grid = element_blank())
p1
cowplot::ggsave2(paste0("./figures/small_events_cell_type_stack_plot_epi.pdf"), p1, width = 4, height = 4)


cell_type_mtx <- mtx_srt_filter_new_if_epi %>% filter(epi_clone == "non_epi" & cell_type != is.na(cell_type)) %>% 
  dplyr::group_by(cell_type) %>%   dplyr::summarise(count = n())
cell_type_mtx$percentage <- round(cell_type_mtx$count / sum(cell_type_mtx$count) * 100, 1)
cell_type_mtx$label <- paste0(cell_type_mtx$percentage, "%")

p1 <- ggplot(cell_type_mtx, aes(x = "", y = percentage, fill = cell_type)) +  
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "non_epi") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  scale_fill_manual(values = rna.cluster.col) + theme_bw() + 
  theme(panel.grid = element_blank())
p1
cowplot::ggsave2(paste0("./figures/small_events_cell_type_stack_plot_non_epi.pdf"), p1, width = 4, height = 4)







