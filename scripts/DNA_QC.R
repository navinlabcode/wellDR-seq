library(copykit)
library(tidyverse)
library(stringr)
library(dplyr)
library(tidyr)
library(cowplot)
library(useful)
library(tibble)
library(SummarizedExperiment)
library(vctrs)
library(ggplot2)
require(RColorBrewer)
library(Homo.sapiens)
library(ggpubr)


source("/wellDR-seq-main/R_scripts/0_wdr_functions.R")
bin_coords <- read.csv("/wellDR-seq-main/pre_load_data/bin_coords.csv")

######preprocess all DNA data
#---wellDR1----####
#-folder contains ./final_result and ./metrics
met_path <- c("./map_seg_output/wellDR1")
met_path2 <- paste0(met_path, "/final_result")
met_path3 <- paste0(met_path, "/metrics")
system(paste0("mkdir -p ", met_path2, "; mkdir -p ", met_path3))
system(paste0("cp ",met_path0, "/final_result/uber.* ", met_path2, "/; cp ", met_path0, "/metrics/all_stat_metrics.txt ",met_path3,"/"))

wellDR1_metrics <- read.table(paste0(met_path, "/metrics/all_stat_metrics.txt"), sep = "\t", header = T) %>% 
  dplyr::rename(sample = "Sample.Name") %>% mutate(dups_percentage = DupsRemoved/TotalReads) 

wellDR1 <- readVarbin(met_path, remove_Y = TRUE, clean_names = FALSE)
wellDR1 <- filterCells(wellDR1, resolution = 0.8)
wellDR1 <- findAneuploidCells(wellDR1)

wellDR1_metadata <- as.data.frame(colData(wellDR1))
wellDR1_metadata_metrics <- left_join(wellDR1_metadata, wellDR1_metrics)
write.table(wellDR1_metadata_metrics,"./metrics/wellDR1_metadata.metrics_newnormal.txt", sep = "\t",row.names = F)


wellDR1 <- wellDR1[,SummarizedExperiment::colData(wellDR1)$filtered == "kept"]
bin_coords2 <- bin_coords %>% dplyr::filter(chrom != 24)
write.table(cbind(bin_coords2,assay(wellDR1,'bin_counts')), 
            "./metrics/wellDR1_filtered_bincounts.txt",row.names = F)


#---wellDR2----####
#-folder contains ./final_result and ./metrics
met_path <- c("./map_seg_output/wellDR2")
met_path2 <- paste0(met_path, "/final_result")
met_path3 <- paste0(met_path, "/metrics")
system(paste0("mkdir -p ", met_path2, "; mkdir -p ", met_path3))
system(paste0("cp ",met_path0, "/final_result/uber.* ", met_path2, "/; cp ", met_path0, "/metrics/all_stat_metrics.txt ",met_path3,"/"))

wellDR2_metrics <- read.table(paste0(met_path, "/metrics/all_stat_metrics.txt"), sep = "\t", header = T) %>% 
  dplyr::rename(sample = "Sample.Name") %>% mutate(dups_percentage = DupsRemoved/TotalReads) 

wellDR2 <- readVarbin(met_path, remove_Y = TRUE, clean_names = FALSE)
wellDR2 <- filterCells(wellDR2, resolution = 0.8)
wellDR2 <- findAneuploidCells(wellDR2)

wellDR2_metadata <- as.data.frame(colData(wellDR2))
wellDR2_metadata_metrics <- left_join(wellDR2_metadata, wellDR2_metrics)
write.table(wellDR2_metadata_metrics,"./metrics/wellDR2_metadata.metrics_newnormal.txt", sep = "\t",row.names = F)


wellDR2 <- wellDR2[,SummarizedExperiment::colData(wellDR2)$filtered == "kept"]
bin_coords2 <- bin_coords %>% dplyr::filter(chrom != 24)
write.table(cbind(bin_coords2,assay(wellDR2,'bin_counts')), 
            "./metrics/wellDR2_filtered_bincounts.txt",row.names = F)

#---wellDR3----####
#-folder contains ./final_result and ./metrics
met_path <- c("./map_seg_output/wellDR3")
met_path2 <- paste0(met_path, "/final_result")
met_path3 <- paste0(met_path, "/metrics")
system(paste0("mkdir -p ", met_path2, "; mkdir -p ", met_path3))
system(paste0("cp ",met_path0, "/final_result/uber.* ", met_path2, "/; cp ", met_path0, "/metrics/all_stat_metrics.txt ",met_path3,"/"))

wellDR3_metrics <- read.table(paste0(met_path, "/metrics/all_stat_metrics.txt"), sep = "\t", header = T) %>% 
  dplyr::rename(sample = "Sample.Name") %>% mutate(dups_percentage = DupsRemoved/TotalReads) 

wellDR3 <- readVarbin(met_path, remove_Y = TRUE, clean_names = FALSE)
wellDR3 <- filterCells(wellDR3, resolution = 0.8)
wellDR3 <- findAneuploidCells(wellDR3)

wellDR3_metadata <- as.data.frame(colData(wellDR3))
wellDR3_metadata_metrics <- left_join(wellDR3_metadata, wellDR3_metrics)
write.table(wellDR3_metadata_metrics,"./metrics/wellDR3_metadata.metrics_newnormal.txt", sep = "\t",row.names = F)


wellDR3 <- wellDR3[,SummarizedExperiment::colData(wellDR3)$filtered == "kept"]
bin_coords2 <- bin_coords %>% dplyr::filter(chrom != 24)
write.table(cbind(bin_coords2,assay(wellDR3,'bin_counts')), 
            "./metrics/wellDR3_filtered_bincounts.txt",row.names = F)


######---- Arc-well-mda231 -----########
#-folder contains ./final_result and ./metrics
met_path <- c("./map_seg_output/arc_well")
met_path2 <- paste0(met_path, "/final_result")
met_path3 <- paste0(met_path, "/metrics")
system(paste0("mkdir -p ", met_path2, "; mkdir -p ", met_path3))
system(paste0("cp ",met_path0, "/final_result/uber.* ", met_path2, "/; cp ", met_path0, "/metrics/all_stat_metrics.txt ",met_path3,"/"))

arc_well_metrics <- read.table(paste0(met_path, "/metrics/all_stat_metrics.txt"), sep = "\t", header = T) %>% 
  dplyr::rename(sample = "Sample.Name") %>% mutate(dups_percentage = DupsRemoved/TotalReads) 

arc_well <- readVarbin(met_path, remove_Y = TRUE, clean_names = FALSE)
arc_well <- filterCells(arc_well, resolution = 0.8)
arc_well <- findAneuploidCells(arc_well)

arc_well_metadata <- as.data.frame(colData(arc_well))
arc_well_metadata_metrics <- left_join(arc_well_metadata, arc_well_metrics)
write.table(arc_well_metadata_metrics,"./metrics/arc_well_metadata.metrics_newnormal.txt", sep = "\t",row.names = F)

arc_well <- arc_well[,SummarizedExperiment::colData(arc_well)$filtered == "kept"]
bin_coords2 <- bin_coords %>% dplyr::filter(chrom != 24)
write.table(cbind(bin_coords2,assay(arc_well,'bin_counts')), 
            "./metrics/arc_well_filtered_bincounts.txt",row.names = F)


######---- 10X CNA -----########
#-folder contains ./final_result and ./metrics
met_path <- c("./map_seg_output/tenx_CNV")
met_path2 <- paste0(met_path, "/final_result")
met_path3 <- paste0(met_path, "/metrics")
system(paste0("mkdir -p ", met_path2, "; mkdir -p ", met_path3))
system(paste0("cp ",met_path0, "/final_result/uber.* ", met_path2, "/; cp ", met_path0, "/metrics/all_stat_metrics.txt ",met_path3,"/"))

tenx_CNV_metrics <- read.table(paste0(met_path, "/metrics/all_stat_metrics.txt"), sep = "\t", header = T) %>% 
  dplyr::rename(sample = "Sample.Name") %>% mutate(dups_percentage = DupsRemoved/TotalReads) 

tenx_CNV <- readVarbin2(met_path, remove_Y = TRUE, clean_names = FALSE)
tenx_CNV <- filterCells(tenx_CNV, resolution = 0.8)
tenx_CNV <- findAneuploidCells(tenx_CNV)

tenx_CNV_metadata <- as.data.frame(colData(tenx_CNV))
tenx_CNV_metadata_metrics <- left_join(tenx_CNV_metadata, tenx_CNV_metrics)
write.table(tenx_CNV_metadata_metrics,"./metrics/tenx_CNV_metadata.metrics_newnormal.txt", sep = "\t",row.names = F)


tenx_CNV <- tenx_CNV[,SummarizedExperiment::colData(tenx_CNV)$filtered == "kept"]
bin_coords2 <- bin_coords %>% dplyr::filter(chrom != 24)
write.table(cbind(bin_coords2,assay(tenx_CNV,'bin_counts')), 
            "./metrics/tenx_CNV_filtered_bincounts.txt",row.names = F)


######---- DLP+ SE 50bp -----########
#-folder contains ./final_result and ./metrics

met_path <- c("./map_seg_output/dlp_plus")
met_path2 <- paste0(met_path, "/final_result")
met_path3 <- paste0(met_path, "/metrics")
system(paste0("mkdir -p ", met_path2, "; mkdir -p ", met_path3))
system(paste0("cp ",met_path0, "/final_result/uber.* ", met_path2, "/; cp ", met_path0, "/metrics/all_stat_metrics.txt ",met_path3,"/"))

dlp_plus_metrics <- read.table(paste0(met_path, "/metrics/all_stat_metrics.txt"), sep = "\t", header = T) %>% 
  dplyr::rename(sample = "Sample.Name") %>% mutate(dups_percentage = DupsRemoved/TotalReads) 

dlp_plus <- readVarbin2(met_path, remove_Y = TRUE, clean_names = FALSE)
dlp_plus <- filterCells(dlp_plus, resolution = 0.8)
dlp_plus <- findAneuploidCells(dlp_plus)

dlp_plus_metadata <- as.data.frame(colData(dlp_plus))
dlp_plus_metadata_metrics <- left_join(dlp_plus_metadata, dlp_plus_metrics)
write.table(dlp_plus_metadata_metrics,"./metrics/dlp_plus_metadata.metrics_newnormal.txt", sep = "\t",row.names = F)


dlp_plus <- dlp_plus[,SummarizedExperiment::colData(dlp_plus)$filtered == "kept"]
bin_coords2 <- bin_coords %>% dplyr::filter(chrom != 24)
write.table(cbind(bin_coords2,assay(dlp_plus,'bin_counts')), 
            "./metrics/dlp_plus_filtered_bincounts.txt",row.names = F)



#---DOP-PCR-SE 50bp----####
met_path <- c("./map_seg_output/DOP_PCR")
met_path2 <- paste0(met_path, "/final_result")
met_path3 <- paste0(met_path, "/metrics")
system(paste0("mkdir -p ", met_path2, "; mkdir -p ", met_path3))
system(paste0("cp ",met_path0, "/final_result/uber.* ", met_path2, "/; cp ", met_path0, "/metrics/all_stat_metrics.txt ",met_path3,"/"))

DOP_PCR_metrics <- read.table(paste0(met_path, "/metrics/all_stat_metrics.txt"), sep = "\t", header = T) %>% 
  dplyr::rename(sample = "Sample.Name") %>% mutate(dups_percentage = DupsRemoved/TotalReads) 

DOP_PCR <- readVarbin(met_path, remove_Y = TRUE, clean_names = FALSE)
DOP_PCR <- filterCells(DOP_PCR, resolution = 0.8)
DOP_PCR <- findAneuploidCells(DOP_PCR)

DOP_PCR_metadata <- as.data.frame(colData(DOP_PCR))
DOP_PCR_metadata_metrics <- left_join(DOP_PCR_metadata, DOP_PCR_metrics)
write.table(DOP_PCR_metadata_metrics,"./metrics/DOP_PCR_metadata.metrics_newnormal.txt", sep = "\t",row.names = F)


DOP_PCR <- DOP_PCR[,SummarizedExperiment::colData(DOP_PCR)$filtered == "kept"]
bin_coords2 <- bin_coords %>% dplyr::filter(chrom != 24)
write.table(cbind(bin_coords2,assay(DOP_PCR,'bin_counts')), 
            "./metrics/DOP_PCR_filtered_bincounts.txt",row.names = F)


######---- DEFND-seq SE 50bp BJ Fibroblast-----########
#-folder contains ./final_result and ./metrics

#####due to data qulity is so bad, we skip the filter cell step for DEFND-seq, otherwise only 6 cells are left
met_path <- c("./map_seg_output/DEFND")
met_path2 <- paste0(met_path, "/final_result")
met_path3 <- paste0(met_path, "/metrics")
system(paste0("mkdir -p ", met_path2, "; mkdir -p ", met_path3))
system(paste0("cp ",met_path0, "/final_result/uber.* ", met_path2, "/; cp ", met_path0, "/metrics/all_stat_metrics.txt ",met_path3,"/"))

DEFND_metrics <- read.table(paste0(met_path, "/metrics/all_stat_metrics.txt"), sep = "\t", header = T) %>% 
  dplyr::rename(sample = "Sample.Name") %>% mutate(dups_percentage = DupsRemoved/TotalReads) 

DEFND <- readVarbin2(met_path, remove_Y = TRUE, clean_names = FALSE)
DEFND <- filterCells(DEFND, resolution = 0.8)
DEFND <- findAneuploidCells(DEFND)

DEFND_metadata <- as.data.frame(colData(DEFND))
DEFND_metadata_metrics <- left_join(DEFND_metadata, DEFND_metrics)
write.table(DEFND_metadata_metrics,"./metrics/DEFND_metadata.metrics_newnormal.txt", sep = "\t",row.names = F)


#DEFND <- DEFND[,SummarizedExperiment::colData(DEFND)$filtered == "kept"]
bin_coords2 <- bin_coords %>% dplyr::filter(chrom != 24)
write.table(cbind(bin_coords2,assay(DEFND,'bin_counts')), 
            "./metrics/DEFND_filtered_bincounts.txt",row.names = F)

##########Calculate overdispersion metrics
pro_name <- c("wellDR1","wellDR2","wellDR3","arc_well","tenx_CNV","dlp_plus","DOP_PCR","DEFND")
tech_name <- c("WellDR1","WellDR2","wellDR3","Arc_well","10X CNV","DLP+","DOP-PCR","DEFND-seq")

all_meta <- data.frame()
for (i in 1:length(pro_name)) {
  print(paste0("Now is runnig: ", pro_name[i]))
  
  filter_cells = read.table(paste0("./metrics/", pro_name[i], "_filtered_bincounts.txt"), header = T, check.names = F)
  meta_file <- read.table(paste0("./metrics/", pro_name[i], "_metadata.metrics_newnormal.txt"), header = T, check.names = F) 
  
  #-----select filtered cells
  varbin_mtx <- readVarbin(paste0("./map_seg_output/", pro_name[i]), remove_Y = TRUE, clean_names = F)
  filt_cells_names <- colnames(filter_cells)[4:ncol(filter_cells)] 
  varbin_mtx_filter <- varbin_mtx[,filt_cells_names]
  merged_meta <- as.data.frame(varbin_mtx_filter@colData) 
  
  #---add meta data
  rownames(meta_file) <- meta_file$sample
  meta_file2 <- meta_file %>% mutate(dups_percentage = DupsRemoved/TotalReads) %>% mutate(sample_name = pro_name[i]) %>% 
    rownames_to_column() %>% 
    dplyr::select(c("rowname", "filter_corr_value","ReadsKept","MedianBinCount","dups_percentage", "sample_name","filtered","is_aneuploid"))
  merged_meta2 <- merged_meta %>% left_join(meta_file2, by = c("sample" = "rowname"))
  varbin_mtx_filter@colData <- cbind(varbin_mtx_filter@colData, 
                                     merged_meta2[,c("filter_corr_value","ReadsKept","MedianBinCount","dups_percentage", "sample_name","filtered","is_aneuploid")])
  
  #----calculate over dispersion---
  bin_count <- varbin_mtx_filter@assays@data$bin_counts
  bin_count_overdisp <- map_dfr(bin_count, overdispersion) %>% t() %>% as.data.frame()
  bin_count_overdisp2 <- bin_count_overdisp %>% rownames_to_column() %>% dplyr::rename(over_disp = V1) %>% dplyr::select("over_disp")
  varbin_mtx_filter@colData <- cbind(varbin_mtx_filter@colData, bin_count_overdisp2) 
  
  write.table(varbin_mtx_filter@colData , paste0("./metrics/", pro_name[i],"_basicQC_meta.txt"), sep = "\t", quote = F, row.names = F)
  all_meta <- rbind(all_meta, as.data.frame(varbin_mtx_filter@colData))
}

all_meta$tech <- plyr::mapvalues(all_meta$sample_name, from = pro_name, to = tech_name)

######downsample to 120 cells per each method
set.seed(1212)
all_meta_s120 <- all_meta %>% group_by(sample_name) %>% sample_n(120) %>% ungroup()

#----tech comparing QC and Figure 2a-----
all_overdisp_s120_comp <- all_meta_s120 %>% mutate(tech = as.factor(tech)) %>% 
  mutate(tech = fct_relevel(tech, tech_name))

pdf("DNA_QC_overdispersion.pdf",height = 6,width = 6,useDingbats = F)
ggplot(all_overdisp_s120_comp) + 
  ggbeeswarm::geom_quasirandom(aes(x = tech, y = over_disp, fill = tech), shape = 21, dodge.width = .8) +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),legend.position = "none", 
        strip.background = element_rect(fill = "white"))  +
  ylab("Overdispersion") + xlab("") 
dev.off()

#########calculate breadth of coverage
#####add corresponding bam file path as the same order as pro_name
bam_path <- c()

for (i in 1: length(pro_name)) {
  my_pro_name <- pro_name[i]
  my_bam_path <- bam_path[i]
  bin_f <- read.table(paste0("./metrics/", my_pro_name, "_filtered_bincounts.txt"), header = T)
  output_path <- paste0("./data/QC_downsampled_reads/", my_pro_name, "/")
  #----Find bam with over 500K reads---
  system(paste0("mkdir -p ./data/QC_downsampled_reads/", my_pro_name, "/data"))
  link_bam_files(bin_f, my_bam_path, output_path, my_pro_name, target_reads = 750000)
  #---link bam with over 500k reads---
  system(paste0("cd ./data/QC_downsampled_reads/", my_pro_name, "/data; cat ../", my_pro_name, " | xargs -I % ln -s % ."))
}

#------downsample bam files to 500k reads and calculate the breadth of coverage---
#----run the snakemake pipeline in terminal--
# source ~/.bashrc
# source ~/.bash_profile
# conda activate snakemake
# snakemake --snakefile ../ds_breadth.smk --cores 10

#----merge breadth of coverage data of all QC samples--
dir.create("./metrics/QC_downsample_reads")

all_cov_breadth <- data.frame()
for (i in 1:length(pro_name)) {
  print(paste0("Now is runnig: ", pro_name[i]))
  my_cov_path <- paste0("./data/QC_downsample_reads/", pro_name[i], "/covfile/")
  my_cov <- calc_coverage(path = my_cov_path) %>% mutate(sample = pro_name[i])
  write.table(my_cov, paste0("./metrics/QC_downsample_reads/",pro_name[i], "_coverage_breadth.txt"), sep = "\t", quote = F, row.names = F)
  all_cov_breadth <- rbind(all_cov_breadth, my_cov)
}

all_cov_breadth$tech <- plyr::mapvalues(all_cov_breadth$sample, from = pro_name, to = tech_name)

#----sampling cells ----
# sampling cells to make group sizes more similar
set.seed(1212)
all_cov_s120 <- all_cov_breadth %>% group_by(sample) %>% sample_n(120) %>% ungroup()

#---plotting breadth of coverage----#####
all_cov_s120_comp <- all_cov_s120 %>% mutate(tech = as.factor(tech)) %>% mutate(tech = fct_relevel(tech, tech_name))

pdf("DNA_QC_breadth_of_coverage.pdf",height = 6,width = 6,useDingbats = F)
ggplot(all_cov_s120_comp) + ggbeeswarm::geom_quasirandom(aes(x = tech, y = breadth, fill = tech), shape = 21, dodge.width = .8) +
  theme_cowplot()  + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), 
        legend.position = "none", strip.background = element_rect(fill = "white"))  +
  # ggtitle("Index of Dispersion") +
  ylab("Breadth of coverage") +
  xlab("") 
dev.off()

##############Construct lorenz curve
library(pracma)  # for the trapz() function
tech_name <- c("WellDR1",   "WellDR2",   "WellDR3",   "Arc_well",  "10X CNV" ,"DLP+",     "DOP-PCR",   "DEFND-seq")
#########first merge bam files from the 120 downsampled single cell bam files from each method. Then calculate whole genome coverage using tools like bedtools or others
###make sure the order of path to the coverage file is the same as tech_name
coverage_file_LIST <-list()
lorenz_data_list <- list()
gini_index_result <- c()
for(i in 1: length(coverage_file_LIST)){
  coverage_file <- coverage_file_LIST[[i]]
  coverage_data <- read.table(coverage_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(coverage_data) <- c("chromosome", "depth", "bases_at_depth", "chrom_size", "fraction")
  genome_data <- coverage_data[coverage_data$chromosome == "genome", ]
  head(genome_data)
  df <- genome_data %>% arrange(depth)
  # Calculate totals
  total_bases <- sum(df$bases_at_depth)
  total_reads <- sum(df$depth * df$bases_at_depth)
  
  # Compute cumulative fractions for bases and reads
  df <- df %>%
    dplyr::mutate(cum_bases = cumsum(as.numeric(bases_at_depth)) / total_bases,
                  cum_reads = cumsum(as.numeric(depth * bases_at_depth)) / total_reads)
  
  # Prepend the (0, 0) point to the cumulative distributions, remove depth=0 row
  x <- c(0, df$cum_bases)
  y <- c(0, df$cum_reads)
  
  # Calculate area under the Lorenz curve using the trapezoidal rule
  area_under_curve <- trapz(x, y)
  # Calculate the Gini index: Gini = 1 - 2 * (area under Lorenz curve)
  gini_index <- 1 - 2 * area_under_curve
  print(paste("Gini index:", gini_index))
  
  # Create a dataframe for plotting
  lorenz_data <- data.frame(cum_bases = x, cum_reads = y)
  lorenz_data$tech <- tech_name[i]
  
  gini_index_result[i] <- gini_index
  names(gini_index_result)[i] <- tech_name[i]
  
  lorenz_data_list[[i]] <- lorenz_data
  names(lorenz_data_list)[i] <- tech_name[i]
}

lorenz_data <- do.call("rbind",lorenz_data_list)
lorenz_data$tech <- factor(lorenz_data$tech,levels =tech_name )

# Plot the Lorenz curve along with the line of equality
pdf("genome_coverage.pdf",height=5,width = 5,useDingbats = F)
ggplot(lorenz_data, aes(x = cum_bases, y = cum_reads,color=tech)) +
  geom_line(size = 0.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
               linetype = "dashed", color = "black")+
  labs(x = "Cumulative fraction of bases", 
       y = "Cumulative fraction of reads",
       title = "Lorenz Curve of Genome Coverage") +
  theme_minimal()+theme(aspect.ratio = 1)
dev.off()

##############compare wellDR-seq MDA231 scDNA-seq data vs bulk MDA231 isogenic clones (78 in total) DNA-seq data
MDA231_Bulk <- read_rds("MDA231_Bulk.rds")
MDA231_merge <- read_rds("wellDR_MDA231_DNA.rds")
MDA231_combiend <- cbind(MDA231_merge,MDA231_Bulk)

#######make bulk consensus line
MDA231_combiend@colData$group <- ifelse(MDA231_combiend@colData$sample %in% MDA231_merge@colData$sample,"SC","Bulk")
table(MDA231_combiend@colData$group)

MDA231_combiend <- calcConsensus(MDA231_combiend, assay = "segment_ratios", fun="mean", consensus_by = "group")

chr_ranges <-
  as.data.frame(SummarizedExperiment::rowRanges(MDA231_combiend))
chr_lengths <- rle(as.numeric(chr_ranges$seqnames))$lengths

# obtaining first and last row of each chr
chr_ranges_start <- chr_ranges %>%
  dplyr::group_by(seqnames) %>%
  dplyr::arrange(seqnames, start) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::ungroup()

chr_ranges_end <- chr_ranges %>%
  dplyr::group_by(seqnames) %>%
  dplyr::arrange(seqnames, start) %>%
  dplyr::filter(dplyr::row_number() == dplyr::n()) %>%
  dplyr::ungroup()

# Creating data frame object for chromosome rectangles shadows
chrom_rects <- data.frame(
  chr = chr_ranges_start$seqnames,
  xstart = as.numeric(chr_ranges_start$abspos),
  xend = as.numeric(chr_ranges_end$abspos)
)
xbreaks <- rowMeans(chrom_rects %>%
                      dplyr::select(
                        xstart,
                        xend
                      ))

if (nrow(chrom_rects) == 24) {
  chrom_rects$colors <- rep(
    c("white", "gray"),
    length(chr_lengths) / 2
  )
} else {
  chrom_rects$colors <- c(rep(
    c("white", "gray"),
    (length(chr_lengths) / 2)
  ), "white")
}

# Creating the geom_rect object
ggchr_back <-
  list(geom_rect(
    data = chrom_rects,
    aes(
      xmin = xstart,
      xmax = xend,
      ymin = -Inf,
      ymax = Inf,
      fill = colors
    ),
    alpha = .2
  ), scale_fill_identity())

sec_breaks <- c(0, 0.5e9, 1e9, 1.5e9, 2e9, 2.5e9, 3e9)
sec_labels <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)


# theme
ggaes <- list(
  scale_x_continuous(
    breaks = xbreaks,
    labels = gsub("chr", "", chrom_rects$chr),
    position = "top",
    expand = c(0, 0),
    sec.axis = sec_axis(
      ~.,
      breaks = sec_breaks,
      labels = sec_labels,
      name = "genome position (Gb)"
    )
  ),
  theme_classic(),
  theme(
    axis.text.x = element_text(
      angle = 0,
      vjust = .5,
      size = 15
    ),
    axis.text.y = element_text(size = 15),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1.3
    ),
    legend.position = "right",
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 15)
  )
)
####################
# obtaining and wrangling data
####################


con <- consensus(MDA231_combiend)

con_l <- con %>%
  dplyr::mutate(abspos = chr_ranges$abspos) %>%
  tidyr::gather(
    key = "group",
    value = "segment_ratio",
    -abspos
  )


pdf("pseudobulk_all_sc_all_bulk_line_plot.pdf",height = 4,width = 10,useDingbats = F)
ggplot(con_l) +
  ggchr_back +
  ggaes +
  geom_line(aes(abspos, segment_ratio, color = group),
            size = 1
  ) +
  labs(
    x = "",
    y = "consensus segment ratio"
  )+scale_color_manual(values = YZRY_color_set2[c(2,3)])
dev.off()

#######calculate correlation of CN profiles between pseudobulk wellDR-seq scDNA-seq data and bulk WGS data 
correlation <- cor.test(MDA231_combiend@consensus$Bulk,MDA231_combiend@consensus$SC,method = "pearson")

