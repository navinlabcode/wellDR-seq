#####-----------Process-------######
library(dplyr)
library(tidyr)
library(cowplot)
library(useful)
library(tibble)
# library(copykit)
library(tidyverse)
library(ggplot2)
require(RColorBrewer)
library(Homo.sapiens)
library(ggpubr)
options(max.print = 200)

setwd("/volumes/USR2/wangkl/wscDR/wdr_analysis/")
source("./scripts/0_wdr_functions.R")
load("./pre_load_data/pre_load_data.rda")
gene_bin <- read.table("./pre_load_data/hg19_gene_binpos_map.tsv", header = 1)
new_pal = c("#CC0C00B2", "#5C88DAB2", "#84BD00B2", "#FFCD00B2", "#7C878EB2", "#00B5E2B2", "#00AF66B2", "#D2AF81B2", 
            "#FD7446B2", "#46732EB2", "#C1395E", "#E07B42","#D4A2D9", "#8E72D5","#C0EDB9", "#364E4F", "#8EE5EE",
            "#FFA500", "#458B00", "#CD6090", "#FFAEB9", "#90EE90", "#5f9EA0", "#E6E6FA", "#8B7E66")
new_pal2 <- new_pal
names(new_pal2) <- paste0("c", 1:length(new_pal2))
tp_col <- c("#EA3291", "#96C942")
colorpal <- c("#66C5CC","#8BE0A4","#F6CF71","#F89C74","#DCB0F2","#87C55F","#9EB9F3","#FE88B1","#C9DB74","#B497E7")

pro_name <- c("wdr_nanowell1","wdr_nanowell2","arc_well", "tenx", "dlp_plus", "dop_merge")
tech_name <- c( "WDR1", "WDR2", "Arc-well", "TenX", "DLP+", "DOP-PCR")
raw_path0 <- paste0("./map_seg_output/", pro_name)

####----calculate over dispersion for cell line and QC data------
#-----create copykit and meta columns of filtered cells---
all_meta <- data.frame()
for (i in 1:length(pro_name)) {
print(paste0("Now is runnig: ", pro_name[i]))

filter_cells = read.table(paste0("./metrics/", pro_name[i], "_filtered_bincounts.txt"), header = T, check.names = F)
meta_file <- read.table(paste0("./metrics/", pro_name[i], "_metadata.metrics_newnormal.txt"), header = T, check.names = F) 

#-----select filtered cells
varbin_mtx <- readVarbinCNA(paste0("./map_seg_output/", pro_name[i]), remove_Y = TRUE, clean_names = F)
filt_cells_names <- colnames(filter_cells)[4:ncol(filter_cells)] 
varbin_mtx_filter <- varbin_mtx[,filt_cells_names]
merged_meta <- as.data.frame(varbin_mtx_filter@colData) 

#---add meta data
rownames(meta_file) <- meta_file$sample
meta_file2 <- meta_file %>% mutate(dups_percentage = DupsRemoved/TotalReads) %>% mutate(sample_name = pro_name[i]) %>% 
  rownames_to_column() %>% 
  dplyr::select(c("rowname", "filter_corr_value","ReadsKept","MedianBinCount","dups_percentage", "sample_name"))
merged_meta2 <- merged_meta %>% left_join(meta_file2, by = c("sample" = "rowname"))
varbin_mtx_filter@colData <- cbind(varbin_mtx_filter@colData, 
                                   merged_meta2[,c("filter_corr_value","ReadsKept","MedianBinCount","dups_percentage", "sample_name")])

#----calculate over dispersion---
bin_count <- varbin_mtx_filter@assays@data$bin_counts
bin_count_overdisp <- map_dfr(bin_count, overdispersion) %>% t() %>% as.data.frame()
bin_count_overdisp2 <- bin_count_overdisp %>% rownames_to_column() %>% dplyr::rename(over_disp = V1) %>% dplyr::select("over_disp")
varbin_mtx_filter@colData <- cbind(varbin_mtx_filter@colData, bin_count_overdisp2) 

write.table(varbin_mtx_filter@colData , paste0("./metrics/", pro_name[i],"_basicQC_meta.txt"), sep = "\t", quote = F, row.names = F)
all_meta <- rbind(all_meta, as.data.frame(varbin_mtx_filter@colData))
}
write.table(all_meta, file = "./metrics/QCsamples_merged_basicQC_meta.txt", sep = "\t", quote = F, row.names = F)

all_meta$tech <- plyr::mapvalues(all_meta$sample_name, from = pro_name, to = tech_name)
write.table(all_meta, file = "./metrics/QCsamples_merged_basicQC_meta_with_tech.txt", sep = "\t", quote = F, row.names = F)
# all_meta <- read.table(file = "./metrics/QCsamples_merged_basicQC_meta_with_tech.txt", sep = "\t", header = T)

#----sampling cells ----
# sampling cells to make group sizes more similar
set.seed(31)
all_meta_s120 <- all_meta %>% group_by(sample_name) %>% sample_n(120) %>% ungroup()

write.table(all_meta_s120, "./metrics/sampled_120cells_QCsamples_merged_basicQC_meta_with_tech.txt", sep = "\t", quote = F, row.names = F)
all_meta_s120 <- read.table(file = "./metrics/sampled_120cells_QCsamples_merged_basicQC_meta_with_tech.txt", sep = "\t", header = T)

#----tech comparing QC and Figure 2a-----
all_overdisp_s120_comp <- all_meta_s120 %>% mutate(tech = as.factor(tech)) %>% 
  mutate(tech = fct_relevel(tech, tech_name))

# all_overdisp_s120_comp %>% dplyr::filter(tech == "Arc-well-fresh-315A") %>% pull(over_disp) %>% median()
# all_overdisp_s120_comp %>% dplyr::filter(tech == "Arc-well-fresh-315A") %>% pull(over_disp) %>% mad()
# all_overdisp_s120_comp %>% dplyr::filter(tech == "Arc-well-MDAMB231") %>% pull(over_disp) %>% median()
# all_overdisp_s120_comp %>% dplyr::filter(tech == "Arc-well-MDAMB231") %>% pull(over_disp) %>% mad()
# all_overdisp_s120_comp %>% dplyr::filter(tech == "ACT-MDA231") %>% pull(over_disp) %>% median()
# all_overdisp_s120_comp %>% dplyr::filter(tech == "ACT-MDA231") %>% pull(over_disp) %>% mad()


p1 <- ggboxplot(all_overdisp_s120_comp, x = "tech", y = "over_disp",color = "tech", palette = "colorpal")+
  stat_compare_means(ref.group = "WDR1")     
p2 <- ggboxplot(all_overdisp_s120_comp, x = "tech", y = "over_disp",color = "tech", palette = "colorpal")+
  stat_compare_means(ref.group = "WDR2")

p3 <- plot_grid(plotlist=list(p1,p2), ncol=1, align='v')
p3
cowplot::ggsave2("./figures/all_tech_comp_overdispersion_wilcox_text.pdf", p3, width = 10, height = 10, limitsize = F)

p1 <- ggplot(all_overdisp_s120_comp) + 
  ggbeeswarm::geom_quasirandom(aes(x = tech, y = over_disp, fill = tech), shape = 21, dodge.width = .8) +
  theme_cowplot() + scale_fill_manual(values = colorpal) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),legend.position = "none", 
        strip.background = element_rect(fill = "white"))  +
  ylab("Overdispersion") + xlab("") 
p1

cowplot::ggsave2("./figures/tech_comp_over_dispersion.pdf", p1, width = 4, height = 4)






####----calculate breadth of coverage for cell line and QC data (750k reads)------
##---Results are very good!!! More sampled reads, the results are better!!!!
pro_name <- c("wdr_nanowell1","wdr_nanowell2","arc_well", "tenx", "dlp_plus", "dop_merge")
#-----find bam files with over 750k reads---####
bam_path <- c("/volumes/USR2/wangkl/wscDR/20211029_WDR_MDA231_chip1R_chip2DR-noNT_SpRNA-ECIS44-cDNA2/data/waferDR_MDA231_chip1_DNA_1M/waferDR_MDA231_chip1-DNA-1M_200/res_200_k/sort/",
              "/volumes/USR2/wangkl/wscDR/20210821_WDR_lowdNTP_231DR_ART216orgR_SpRNA-ECIS42-SnubarARC-RNA2/data/waferDR_231_DNA_1M/waferDR_231_DNA_1M_200/res_200_k/sort/",
              "/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/20210206_MDA231/data/MDA231-9x_200/res_200_k/sort/",
              "/volumes/seq/projects/CNA_projects/10X_CNA/Breast/TN7_TN17/TN17_CBS_output/CBS_output/output/sort/",
              "/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/dop_pcr_reprocess/dop_pcr_merged_200/res_200_k/sort/",
              "/volumes/seq/external_data/laks_2019/dlp_plus/varbin_200kb/output/sort/")

for (i in 1:6) {
  my_pro_name <- pro_name[i]
  my_bam_path <- bam_path[i]
  bin_f <- read.table(paste0("./metrics/", my_pro_name, "_filtered_bincounts.txt"), header = T)
  output_path <- paste0("./data/QC_downsample_750k/", my_pro_name, "/")
  #----Find bam with over 500K reads---
  system(paste0("mkdir -p ./data/QC_downsample_750k/", my_pro_name, "/data"))
  link_bam_files(bin_f, my_bam_path, output_path, my_pro_name, target_reads = 750000)
  #---link bam with over 500k reads---
  system(paste0("cd ./data/QC_downsample_750k/", my_pro_name, "/data; cat ../", my_pro_name, " | xargs -I % ln -s % ."))
}

#------downsample bam files to 750k reads and calculate the breadth of coverage---
#----run the snakemake pipeline in terminal--
# source ~/.bashrc
# source ~/.bash_profile
# conda activate snakemake
# snakemake --snakefile ./scripts/snakemake_files/ds_breadth.smk --cores 20

#----merge breadth of coverage data of all QC samples--
dir.create("./metrics/breadth_cov_750k")
all_cov_breadth <- data.frame()
for (i in 1:length(pro_name)) {
  print(paste0("Now is runnig: ", pro_name[i]))
  my_cov_path <- paste0("./data/QC_downsample_750k/", pro_name[i], "/covfile/")
  my_cov <- calc_coverage(path = my_cov_path) %>% mutate(sample = pro_name[i])
  write.table(my_cov, paste0("./metrics/breadth_cov_750k/",pro_name[i], "_coverage_breadth.txt"), sep = "\t", quote = F, row.names = F)
  all_cov_breadth <- rbind(all_cov_breadth, my_cov)
}

write.table(all_cov_breadth, "./metrics/QCsamples_merged_coverage_breadth_750k.txt", sep = "\t", quote = F, row.names = F)
all_cov_breadth <- read.table(file = "./metrics/QCsamples_merged_coverage_breadth_750k.txt", sep = "\t", header = T)

all_cov_breadth$tech <- plyr::mapvalues(all_cov_breadth$sample, from = pro_name, to = tech_name)
write.table(all_cov_breadth, file = "./metrics/QCsamples_merged_coverage_breadth_with_tech_750k.txt", sep = "\t", quote = F, row.names = F)
# all_cov_breadth <- read.table(file = "./metrics/QCsamples_merged_coverage_breadth_with_tech.txt", sep = "\t", header = T)

#----sampling cells ----
# sampling cells to make group sizes more similar
set.seed(35)
all_cov_s120 <- all_cov_breadth %>% group_by(sample) %>% sample_n(120) %>% ungroup()
write.table(all_cov_s120, "./metrics/sampled_120cells_QCsamples_merged_coverage_breadth_with_tech.txt", 
            sep = "\t", quote = F, row.names = F)
all_cov_s120 <- read.table("./metrics/sampled_120cells_QCsamples_merged_coverage_breadth_with_tech.txt", sep = "\t", header = T)
#---plotting figure 2b----#####
all_cov_s120_comp <- all_cov_s120 %>% mutate(tech = as.factor(tech)) %>% mutate(tech = fct_relevel(tech, tech_name))

p1 <- ggplot(all_cov_s120_comp) + ggbeeswarm::geom_quasirandom(aes(x = tech, y = breadth, fill = tech), shape = 21, dodge.width = .8) +
  theme_cowplot() + scale_fill_manual(values = colorpal) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), 
        legend.position = "none", strip.background = element_rect(fill = "white"))  +
  # ggtitle("Index of Dispersion") +
  ylab("Breadth of coverage") +
  xlab("") 
p1
cowplot::ggsave2("./figures/tech_comp_breadth_cov_750k.pdf", p1, width = 4, height = 4)

all_cov_s120_comp %>% dplyr::filter(tech == "WDR1") %>% pull(breadth) %>% median()
all_cov_s120_comp %>% dplyr::filter(tech == "WDR1") %>% pull(breadth) %>% mad()
all_cov_s120_comp %>% dplyr::filter(tech == "WDR2") %>% pull(breadth) %>% median()
all_cov_s120_comp %>% dplyr::filter(tech == "WDR2") %>% pull(breadth) %>% mad()

ggboxplot(all_cov_s120_comp, x = "tech", y = "breadth",color = "tech", palette = "colorpal")+
  stat_compare_means(ref.group = "WDR1")     

ggboxplot(all_cov_s120_comp, x = "tech", y = "breadth",color = "tech", palette = "colorpal")+
  stat_compare_means(ref.group = "WDR2")

plist <- list()
j=1
for (i in 1:200) {
  set.seed(i)
  all_cov_s80 <- all_cov_breadth %>% group_by(sample) %>% sample_n(120) %>% ungroup()
  all_cov_s80_comp <- all_cov_s80 %>% mutate(tech = as.factor(tech)) %>% mutate(tech = fct_relevel(tech, tech_name))
  
  p1 <- ggplot(all_cov_s80_comp) + ggbeeswarm::geom_quasirandom(aes(x = tech, y = breadth, fill = tech), shape = 21, dodge.width = .8) +
    theme_cowplot() + scale_fill_manual(values = colorpal) +
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
          legend.position = "none", strip.background = element_rect(fill = "white"))  +
    # ggtitle("Index of Dispersion") +
    ylab("Breadth of coverage") +
    xlab(i)
  plist[[j]] <- p1
  j <- j + 1
}

p_list1<-plot_grid(plotlist = plist, ncol = 20)
cowplot::ggsave2("./figures/tech_comp_breadth_cov_test.pdf", p_list1, width = 4*20, height = 4*10, limitsize = F)





#---plot overdispersion and breadth of coverage together----#####
p1 <- ggplot(all_overdisp_s120_comp) + 
  ggbeeswarm::geom_quasirandom(aes(x = tech, y = over_disp, fill = tech), shape = 21, dodge.width = .8) +
  theme_cowplot() + scale_fill_manual(values = colorpal) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),legend.position = "none", 
        strip.background = element_rect(fill = "white"))  +
  ylab("Overdispersion") + xlab("") 

p2 <- ggplot(all_cov_s120_comp) + ggbeeswarm::geom_quasirandom(aes(x = tech, y = breadth, fill = tech), shape = 21, dodge.width = .8) +
  theme_cowplot() + scale_fill_manual(values = colorpal) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), 
        legend.position = "none", strip.background = element_rect(fill = "white"))  +
  ylab("Breadth of coverage") + xlab("") #+ scale_y_continuous(labels = scales::scientific_format())
p2
p3 <- plot_grid(plotlist=list(p1,p2), ncol=1, align='v')
p3
cowplot::ggsave2("./figures/all_tech_comp_overdispersion_breadthcov.pdf", p3, width = 4, height = 7.5)


