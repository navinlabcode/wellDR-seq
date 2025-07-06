#########Compare median genes per cell across different chemistries
###########RNA performance comparisons of wellDR-seq vs Takara and 10X data
#######load downsampled (~26K reads per cell) scRNA-seq data


#####Takeara Bio 3'DE MDA231 data
takara_231 <- subset(takara_231, subset = nFeature_RNA > 500 & nFeature_RNA<10000 & percent.mt < 10)
takara_231 <- NormalizeData(takara_231)
takara_231$experiment <- "takaraDE"

############10X MDA231 data
tenx_231 <- subset(tenx_231, subset = nFeature_RNA > 500 & nFeature_RNA<10000 & percent.mt < 10)
tenx_231 <- NormalizeData(tenx_231)
tenx_231$experiment <- "tenx"

##########wellDR-seq scRNA-seq data
WDR_Wafar_231_RNA_chip1 <- subset(WDR_Wafar_231_RNA_chip1, subset = nFeature_RNA > 500 & nFeature_RNA<10000 & percent.mt < 10)
WDR_Wafar_231_RNA_chip1$experiment <- "WellDR1"

WDR_Wafar_231_RNA_chip2 <- subset(WDR_Wafar_231_RNA_chip2, subset = nFeature_RNA > 500 & nFeature_RNA<10000 & percent.mt < 10)
WDR_Wafar_231_RNA_chip2$experiment <- "WellDR2"

WDR_Wafar_231_RNA_chip3 <- subset(WDR_Wafar_231_RNA_chip3, subset = nFeature_RNA > 500 & nFeature_RNA<10000 & percent.mt < 10)
WDR_Wafar_231_RNA_chip3$experiment <- "WellDR3"

######DEFND-seq data
DEFEND_RNA <- subset(DEFEND_RNA, subset = nFeature_RNA > 500 & nFeature_RNA<10000 & percent.mt < 10)
DEFEND_RNA <- NormalizeData(DEFEND_RNA)
DEFEND_RNA$experiment <- "DEFEND-seq"

#######merge data across platforms and perfrom QC comparisons
merged_231_acroos_platform_RNA <- merge(WDR_Wafar_231_RNA_chip1,c(WDR_Wafar_231_RNA_chip2,WDR_Wafar_231_RNA_chip3,DEFEND_RNA,takara_231,tenx_231),add.cell.ids=c("WellDR1","WellDR2","WellDR3","DEFEND-seq","Takara","Tenx"))
merged_231_acroos_platform_RNA$experiment <- factor(merged_231_acroos_platform_RNA$experiment,levels = c("WellDR1","WellDR2","WellDR3","DEFEND-seq","takaraDE","tenx"))
Idents(merged_231_acroos_platform_RNA) <- "experiment"


####downsample to same cell numbers with 1000 permutations
# Load necessary libraries
library(dplyr)


# Function to perform downsampling, Wilcoxon test, and calculate mean group values
downsample_and_test <- function(df) {
  # Downsample each group to 800 rows
  sampled_df <- df %>%
    group_by(experiment) %>%
    sample_n(500) %>%
    ungroup()
  
  # Perform pairwise Wilcoxon rank sum test with Bonferroni adjustment
  result <- pairwise.wilcox.test(sampled_df$nFeature_RNA, sampled_df$experiment, p.adjust.method = "bonferroni")
  
  # Calculate mean group values
  mean_values <- sampled_df %>%
    group_by(experiment) %>%
    summarise(mean_value = mean(nFeature_RNA)) %>%
    ungroup()
  
  return(list(p_values = result$p.value, means = mean_values))
}

# Number of iterations
n_iterations <- 1000

# Initialize lists to store p-value matrices and mean values
p_value_list <- vector("list", n_iterations)
mean_value_list <- vector("list", n_iterations)

# Perform downsampling and testing 1000 times
set.seed(1212)  # Set seed for reproducibility
for (i in 1:n_iterations) {
  res <- downsample_and_test(data_comparison)
  p_value_list[[i]] <- res$p_values
  mean_value_list[[i]] <- res$means
}

# Convert list of p-value matrices to a 3D array
p_value_array <- simplify2array(p_value_list)

# Calculate the median p-values across iterations
final_p_value_matrix <- apply(p_value_array, c(1, 2), median, na.rm = TRUE)

# Convert list of mean values to a dataframe
mean_value_df <- bind_rows(mean_value_list, .id = "iteration")

# Calculate the mean of mean values for each group
final_mean_values <- mean_value_df %>%
  group_by(experiment) %>%
  summarise(mean_of_means = mean(mean_value, na.rm = TRUE))


mean_value_df$experiment <- factor(mean_value_df$experiment,levels = c("WellDR1","WellDR2","WellDR3","DEFEND-seq","takaraDE","tenx"))

pdf("RNA_QC_gene_count.pdf",height = 6,width = 6,useDingbats = F)
ggplot(mean_value_df, aes(y = mean_value, x = experiment,color=experiment)) + geom_boxplot(position = position_dodge2(0.9),outlier.shape=NA)+labs(x=NULL,y="Mean genes per cell from each permutation")+geom_point(position=position_jitterdodge(jitter.width=1.5),size=0.01)+theme_classic2()+theme(axis.text.x = element_text(size = 12, face = "bold",angle = 90,hjust = 1,vjust = 0.5),axis.text.y = element_text(size = 16,face = "bold"),axis.title=element_text(size=16,face = "bold"))+NoLegend()+ylim(0,3500)
dev.off()

###########Calculate correlations between wellDR-seq and Takera Bio and 10X data
mda231_combined <- merge(WDR_Wafar_231_RNA_chip1,c(WDR_Wafar_231_RNA_chip2,WDR_Wafar_231_RNA_chip3,takara_231,tenx_231))

#####construct pseudobulk with deseq2 
exp_pseudobulk <- AggregateExpression(mda231_combined,assays = "RNA",slot = "counts",return.seurat = T,group.by = "experiment")

exp_pseudobulk <- data.frame(exp_pseudobulk@assays$RNA@counts)
meta_bulk <- data.frame(colnames(exp_pseudobulk))
rownames(meta_bulk) <- meta_bulk$colnames.exp_pseudobulk.
library(DESeq2)
exp_dds_bulk <- DESeqDataSetFromMatrix(countData = exp_pseudobulk,
                                       colData = meta_bulk,
                                       design = ~ 1)
exp_dds_bulk <- vst(exp_dds_bulk, blind=TRUE)
exp_dds_bulk <- assay(exp_dds_bulk)

exp_dds_bulk <- data.frame(exp_dds_bulk)

library(corrplot)
M<-cor(exp_dds_bulk,method = "pearson")
pdf("Gene_expression_correlation_WDR_10X_Takara2.pdf",height = 6,width = 6,useDingbats = F)
corrplot(M, method="color",type="lower",order="hclust",addCoef.col = "white",tl.col="black",diag=F,col  = colorRampPalette(c("red","white","#00BFFF"))(100))
dev.off()

