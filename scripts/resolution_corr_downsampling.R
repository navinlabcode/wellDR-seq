library(dplyr)
library(tidyr)
library(data.table)
library(purrr)
library(stringr)
library(ggbeeswarm)
library(ggplot2)
library(tidyverse)
# function 
corr_cells <- function(original,
                       downsampled) {
  
  ###
  # Loops across all cells matching by their names and calculates spearman
  # correlation
  
  
  # original: uber.sample.seg.txt matrix of segment ratios from original files
  # downsampled: uber.sample.seg.txt matrix of segment ratios from downsampled files
  ###
  
  # keeping only cells that were downsampled
  original <- original[,names(downsampled)]
  
  # sanity check
  stopifnot(identical(names(downsampled), names(original)))
  
  # cell names vector
  cell_names <- names(downsampled)[-c(1:3)]
  
  # running correlation
  cor_list <- BiocParallel::bplapply(cell_names, function(x)
    cor(downsampled[,x], original[,x], method = 'spearman')
  )
  
  # binding to a vector
  cor_vector <- do.call(c, cor_list)
  names(cor_vector) <- cell_names
  cor_vector
  
}

create_corr_df <- function(cor_1M,
                           cor_750k,
                           cor_500k,
                           cor_250k,
                           cor_125k,
                           cor_75k,
                           cor_50k) {
  
  ###
  # returns a data frame of all correlations
  ###
  cat(
    length(cor_1M),
    length(cor_750k),
    length(cor_500k),
    length(cor_250k),
    length(cor_125k),
    length(cor_75k),
    length(cor_50k)
  )
  
  data.frame(cell = c(
    names(cor_1M),
    names(cor_750k),
    names(cor_500k),
    names(cor_250k),
    names(cor_125k),
    names(cor_75k),
    names(cor_50k)
  ),
  correlation = c(cor_1M,
                  cor_750k,
                  cor_500k,
                  cor_250k,
                  cor_125k,
                  cor_75k,
                  cor_50k),
  n_reads = c(rep('1M', length(cor_1M)),
              rep('750k', length(cor_750k)),
              rep('500k', length(cor_500k)),
              rep('250k', length(cor_250k)),
              rep('125k', length(cor_125k)),
              rep('75k', length(cor_75k)),
              rep('50k', length(cor_50k))
  )) %>% 
    mutate(n_reads = fct_relevel(n_reads, n_reads_levels))
  
}

# load segmentation data from cnv pipeline
MDA231_chip1_200k_seg <- read.delim("./uber.MDA231_chip1_200k.seg.txt") 
MDA231_chip1_100k_seg <- read.delim("./uber.MDA231_chip1_100k.seg.txt")
MDA231_chip1_50k_seg <- read.delim("./uber.outdir_chip1_sorted_bam_50k.seg.txt")


MDA231_chip1_200k_seg$chrompos <- MDA231_chip1_200k_seg$chrompos+1
MDA231_chip1_200k_seg$abspos <- MDA231_chip1_200k_seg$abspos+1

MDA231_chip1_50k_seg$chrompos <- MDA231_chip1_50k_seg$chrompos+1
MDA231_chip1_50k_seg$abspos <- MDA231_chip1_50k_seg$abspos+1

#compare 50k and 100k
mtx_compare_50k_100k <- MDA231_chip1_50k_seg %>%
  left_join(MDA231_chip1_100k_seg, by = "abspos") %>%
  select(-chrom.y,-chrompos.y)
mtx_compare_50k_100k <- as.data.table(mtx_compare_50k_100k)
#find all col end up with y
cols_to_fill_50k_100k <- names(mtx_compare_50k_100k) %>% 
  stringr::str_subset("\\.y$") 
#fill NA
mtx_compare_50k_100k[, (cols_to_fill_50k_100k) := lapply(.SD, nafill, type = "locf"), .SDcols = cols_to_fill_50k_100k]
mtx_compare_50k_100k[, (cols_to_fill_50k_100k) := lapply(.SD, function(x) {
  if (is.na(x[1]) & length(x) > 1) {
    replace(x, 1, x[2])
  } else {
    x
  }
}), .SDcols = cols_to_fill_50k_100k]
cols_50k <- names(mtx_compare_50k_100k) %>% stringr::str_subset("\\.x$")
cols_100k <- names(mtx_compare_50k_100k) %>% stringr::str_subset("\\.y$")
#calculate correlation by cell pair 
prefixes <- gsub("\\.x$", "", cols_50k)
matched_pairs <- setNames(map2(cols_50k, paste0(prefixes, ".y"), ~ .x), prefixes)

corr_test_50k_100k <- map2_df(
  matched_pairs,
  map(matched_pairs, ~ mtx_compare_50k_100k[[.x]]),
  ~ {
    col_y <- str_replace(.x, "\\.x$", ".y")
    result <- tryCatch(
      cor.test(mtx_compare_50k_100k[[.x]], mtx_compare_50k_100k[[col_y]], method = "spearman"),
      error = function(e) list(p.value = NA, estimate = NA, statistic = NA)
    )
    tibble(
      prefix = str_remove(.x, "\\.x$"),
      method = "Spearman",
      rho = result$estimate,
      p_value = result$p.value,
      statistic = result$statistic
    )
  }
)

corr_test_50k_100k<-corr_test_50k_100k[-1,]
corr_test_50k_100k<-corr_test_50k_100k[-1,]



#compare 50k and 200k
mtx_compare_50k_200k <- MDA231_chip1_50k_seg %>%
  left_join(MDA231_chip1_200k_seg, by = "abspos") %>%
  select(-chrom.y,-chrompos.y)
mtx_compare_50k_200k <- as.data.table(mtx_compare_50k_200k)
cols_to_fill_50k_200k <- names(mtx_compare_50k_200k) %>% 
  stringr::str_subset("\\.y$") 
mtx_compare_50k_200k[, (cols_to_fill_50k_200k) := lapply(.SD, nafill, type = "locf"), .SDcols = cols_to_fill_50k_200k]
mtx_compare_50k_200k[, (cols_to_fill_50k_200k) := lapply(.SD, function(x) {
  if (is.na(x[1]) & length(x) > 1) {
    replace(x, 1, x[2])
  } else {
    x
  }
}), .SDcols = cols_to_fill_50k_200k]
cols_50k <- names(mtx_compare_50k_200k) %>% stringr::str_subset("\\.x$")
cols_200k <- names(mtx_compare_50k_200k) %>% stringr::str_subset("\\.y$")

prefixes <- stringr::str_remove(cols_50k, "\\.x$")
matched_pairs <- setNames(map2(cols_50k, paste0(prefixes, ".y"), ~ .x), prefixes)

corr_test_50k_200k <- map2_df(
  matched_pairs,
  map(matched_pairs, ~ mtx_compare_50k_200k[[.x]]),
  ~ {
    col_y <- str_replace(.x, "\\.x$", ".y")
    result <- tryCatch(
      cor.test(mtx_compare_50k_200k[[.x]], mtx_compare_50k_200k[[col_y]], method = "spearman"),
      error = function(e) list(p.value = NA, estimate = NA, statistic = NA)
    )
    tibble(
      prefix = str_remove(.x, "\\.x$"),
      method = "Spearman",
      rho = result$estimate,
      p_value = result$p.value,
      statistic = result$statistic
    )
  }
)
corr_test_50k_200k<-corr_test_50k_200k[-1,]
corr_test_50k_200k<-corr_test_50k_200k[-1,]

#compare 100k and 200k
mtx_compare_100k_200k <- MDA231_chip1_100k_seg %>%
  left_join(MDA231_chip1_200k_seg, by = "abspos") %>%
  select(-chrom.y,-chrompos.y)
mtx_compare_100k_200k <- as.data.table(mtx_compare_100k_200k)
cols_to_fill_100k_200k <- names(mtx_compare_100k_200k) %>% 
  stringr::str_subset("\\.y$") 
mtx_compare_100k_200k[, (cols_to_fill_100k_200k) := lapply(.SD, nafill, type = "locf"), .SDcols = cols_to_fill_100k_200k]
mtx_compare_100k_200k[, (cols_to_fill_100k_200k) := lapply(.SD, function(x) {
  if (is.na(x[1]) & length(x) > 1) {
    replace(x, 1, x[2])
  } else {
    x
  }
}), .SDcols = cols_to_fill_100k_200k]
cols_100k <- names(mtx_compare_100k_200k) %>% stringr::str_subset("\\.x$")
cols_200k <- names(mtx_compare_100k_200k) %>% stringr::str_subset("\\.y$")

prefixes <- stringr::str_remove(cols_100k, "\\.x$")
matched_pairs <- setNames(map2(cols_100k, paste0(prefixes, ".y"), ~ .x), prefixes)

corr_test_100k_200k <- map2_df(
  matched_pairs,
  map(matched_pairs, ~ mtx_compare_100k_200k[[.x]]),
  ~ {
    col_y <- str_replace(.x, "\\.x$", ".y")
    result <- tryCatch(
      cor.test(mtx_compare_100k_200k[[.x]], mtx_compare_100k_200k[[col_y]], method = "spearman"),
      error = function(e) list(p.value = NA, estimate = NA, statistic = NA)
    )
    tibble(
      prefix = str_remove(.x, "\\.x$"),
      method = "Spearman",
      rho = result$estimate,
      p_value = result$p.value,
      statistic = result$statistic
    )
  }
)
corr_test_100k_200k<-corr_test_100k_200k[-1,]
corr_test_100k_200k<-corr_test_100k_200k[-1,]


plot_data <- rbind(data.frame(rho = corr_test_50k_100k$rho, test = "50k_100k"),
                   data.frame(rho = corr_test_50k_200k$rho, test = "50k_200k"),
                   data.frame(rho = corr_test_100k_200k$rho, test = "100k_200k"))

p1 <- ggplot(plot_data,aes(x=test,y=as.numeric(rho),fill = test,color = test))+
  geom_point(size=0.1)+
  geom_jitter(size=0.1)+
  ylim(0,1) + theme_bw() +theme(panel.grid = element_blank(),
                                axis.text.x = element_text(angle = 45))
p1

mean(corr_test_50k_100k$rho)
sd(corr_test_50k_100k$rho)
mean(corr_test_50k_200k$rho)
sd(corr_test_50k_200k$rho)
mean(corr_test_100k_200k$rho)
sd(corr_test_100k_200k$rho)





n_reads_levels <- c('1M', '750k', '500k', '250k', '125k', '75k', '50k')

wafer_mda231p_original_seg <- read.table("/sibcb1/wangkailelab1/pengjinyu/mda231/uber.MDA231_chip1_200k.seg.txt", header = T)
wafer_mda231p_1M_seg <- read.table("./downsample_corr/wdr_mda231_chip1/wdr_mda231_chip1_1M/final_result/uber.wdr_mda231_chip1_1M.seg.txt", header = T)
names(wafer_mda231p_original_seg) <- str_remove(names(wafer_mda231p_original_seg), '.sorted')
wafer_mda231p_750k_seg <- read.table("./downsample_corr/wdr_mda231_chip1/wdr_mda231_chip1_750k/final_result/uber.wdr_mda231_chip1_750k.seg.txt", header = T)
wafer_mda231p_500k_seg <- read.table("./downsample_corr/wdr_mda231_chip1/wdr_mda231_chip1_500k/final_result/uber.wdr_mda231_chip1_500k.seg.txt", header = T)
wafer_mda231p_250k_seg <- read.table("./downsample_corr/wdr_mda231_chip1/wdr_mda231_chip1_250k/final_result/uber.wdr_mda231_chip1_250k.seg.txt", header = T)
wafer_mda231p_125k_seg <- read.table("./downsample_corr/wdr_mda231_chip1/wdr_mda231_chip1_125k/final_result/uber.wdr_mda231_chip1_125k.seg.txt", header = T)
wafer_mda231p_75k_seg <- read.table("./downsample_corr/wdr_mda231_chip1/wdr_mda231_chip1_75k/final_result/uber.wdr_mda231_chip1_75k.seg.txt", header = T)
wafer_mda231p_50k_seg <- read.table("./downsample_corr/wdr_mda231_chip1/wdr_mda231_chip1_50k/final_result/uber.wdr_mda231_chip1_50k.seg.txt", header = T)


wafer_mda231p_cor_1M <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_1M_seg)
wafer_mda231p_cor_750k <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_750k_seg)
wafer_mda231p_cor_500k <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_500k_seg)
wafer_mda231p_cor_250k <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_250k_seg)
wafer_mda231p_cor_125k <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_125k_seg)
wafer_mda231p_cor_75k <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_75k_seg)
wafer_mda231p_cor_50k <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_50k_seg)

wafer_mda231p_corr_df <- create_corr_df(wafer_mda231p_cor_1M,
                                        wafer_mda231p_cor_750k,
                                        wafer_mda231p_cor_500k,
                                        wafer_mda231p_cor_250k,
                                        wafer_mda231p_cor_125k,
                                        wafer_mda231p_cor_75k,
                                        wafer_mda231p_cor_50k
)


p1 <- ggplot(wafer_mda231p_corr_df, aes(x = n_reads, y = correlation)) + 
  ylim(c(0, 1)) + 
  geom_jitter(aes(color = n_reads),position = position_jitter(seed = 1, width = 0.2), size = 0.01) +
  geom_violin(alpha=0, size = 0.2) + 
  stat_summary(fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75)) + 
  cowplot::theme_cowplot() +
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p1

wafer_mda231p_corr_df %>% dplyr::group_by(n_reads) %>% dplyr::summarise(cor=mean(correlation),sd=sd(correlation))
summary_df <- wafer_mda231p_corr_df %>%
  dplyr::group_by(n_reads) %>%
  dplyr::summarise(
    mean_cor = mean(correlation),
    sd_cor = sd(correlation),
    ymin = mean_cor - sd_cor,
    ymax = mean_cor + sd_cor
  ) %>%
  mutate(label = sprintf("%.2f±%.3f", mean_cor, sd_cor))

p2 <- p1 + 
  geom_text(data = summary_df,
            aes(x = n_reads, y = 1, label = label),
            vjust = -0.5, size = 2, color = "red")
p2

cowplot::ggsave2("./figures/downsample_mda231p_spearman_cor_pvln.pdf", p2, width = 6, height = 5)
