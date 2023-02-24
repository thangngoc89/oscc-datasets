# Set working directory
setwd("~/Library/CloudStorage/OneDrive-UMP/NCKH_Nhatnam")
options(stringsAsFactors = FALSE)

#CPU core
options(mc.cores = parallel::detectCores())

# Load library for RColorBrewer
library(RColorBrewer)
# Load library for pheatmap
library(pheatmap)
# Load library for tidyverse
library(tidyverse)
# Load rhdf5
library(rhdf5)
# Load tximport
library(tximport)
# Including gene names into transcript-level analysis
library(biomaRt)
# Load DESeq2
library(DESeq2)
# other pakage
library(ggplot2)
# ConsensusClusterPlus
library(ConsensusClusterPlus)


# Set color
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heat_colors <- brewer.pal(n = 6, name = "YlOrRd")


# Read metadata
metadata <- read.table("./11.Data/R_output/230216/metadata_withclusters.csv", sep= ",",header = TRUE)

#Load data
txi.kallisto <- readRDS("./11.Data/R_output/230216/txi.kallisto.RDS")

# Import with cluster information
metadata$cluster2 <- factor(metadata$cluster2, levels = c("1","2"))
metadata$cluster3 <- factor(metadata$cluster3, levels = c("1","2","3"))
metadata$cluster4 <- factor(metadata$cluster4, levels = c("1","2","3","4"))
metadata$cluster5 <- factor(metadata$cluster5, levels = c("1","2","3","4","5"))

# 3 clusters
ddsMat <- DESeqDataSetFromTximport(txi.kallisto, metadata, ~ cluster2)

# Pre-filtering the dataset
nrow(ddsMat)
ddsMat <- ddsMat[rowSums(counts(ddsMat)) > 1, ] # keep genes with total counts >1
nrow(ddsMat)

# Determine the size factors to use for normalization
dds1 <- estimateSizeFactors(ddsMat)
sizeFactors(dds1)

# Differential expression analysis
dds2 <- DESeq(dds1)
saveRDS(dds2, "./11.Data/R_output/230216/dds2.RDS")
#dds2 <- readRDS("./11.Data/R_output/230216/dds2.RDS")

resultsNames(dds2)
res <- results(dds2)
mcols(res, use.names = TRUE)
summary(res)

#cluster3_2_vs_1
# lower the false discovery rate threshold
res2.05 <- results(dds2, alpha = 0.05, name = "cluster2_2_vs_1")
summary(res2.05)

# Shrink the log2 fold change estimates to be more accurate
dds2_res <- lfcShrink(dds2, coef=  "cluster2_2_vs_1",type="apeglm"
                      ,res = res2.05)
summary(dds2_res)

pdf("./11.Data/R_output/230216/fig_4_MA_plot_cluster2_2_vs_1.pdf", useDingbats = F, height = 6, width =6)
plotMA(dds2_res, ylim=c(-14,12), colSig = "cadetblue")
dev.off()

# Subset the Output3 to only return the significant genes with p-adjusted values less than 0.05
dds2_res_all <- data.frame(dds2_res)
dds2_res_all$Gene <- rownames(dds2_res_all)
dds2_res_sig <- subset(dds2_res_all, padj < 0.01)
dds2_res_sig <- dplyr :: arrange(dds2_res_sig, desc(dds2_res_sig$log2FoldChange))
write.csv(dds2_res_sig,"./11.Data/R_output/230216/cluster2_2_vs_1.csv")

#cluster3_3_vs_1
# lower the false discovery rate threshold
#res3.05 <- results(dds2, alpha = 0.05, name = "cluster3_3_vs_1")
#summary(res3.05)

# Shrink the log2 fold change estimates to be more accurate
#dds3_res <- lfcShrink(dds2, coef=  "cluster3_3_vs_1",type="normal"
#                      ,res = res)
#summary(dds3_res)

#pdf("./11.Data/R_output/230216/fig_4_MA_plot_cluster3_3_vs_1.pdf", useDingbats = F, height = 6, width =6)
#plotMA(dds3_res, ylim=c(-8,8), colSig = "cadetblue")
#dev.off()

# Subset the Output3 to only return the significant genes with p-adjusted values less than 0.05
#dds3_res_all <- data.frame(dds3_res)
#dds3_res_all$Gene <- rownames(dds3_res_all)
#dds3_res_sig <- subset(dds3_res_all, padj < 0.1)
#dds3_res_sig <- dplyr :: arrange(dds3_res_sig, desc(dds3_res_sig$log2FoldChange))
#write.csv(dds3_res_sig,"./11.Data/R_output/230216/cluster3_3_vs_1")


# Group by clusters
metadata2 <- dplyr :: arrange(metadata, metadata$cluster2)
df <- as.data.frame(metadata2$cluster2)
rownames(df) <- metadata2$SampleID
colnames(df) <- "Group"


Gene <- unique(dds2_res_sig$Gene
#               , dds3_res_sig$Gene
               )

# Subset normalized counts to significant genes
sig_norm_counts_all <- counts(dds2,normalized=T)[Gene, ]
sig_norm_counts_all <- sig_norm_counts_all[,metadata2$X]


# Plot heatmap
pdf("./11.Data/R_output/230216/fig_5_heatmap.pdf", useDingbats = F, height = 6, width = 4)
pheatmap(sig_norm_counts_all, 
         color = heat_colors,
         #cutree_cols=2,
         cluster_rows = T, cluster_cols=F,
         #clustering_distance_cols="manhattan",
         #clustering_method = "single",
         show_colnames=F,
         show_rownames =F,
         annotation_col = df, 
         fontsize_col= 7,
         scale = "row", main= "Significant different genes")
dev.off()

#plotCounts


dds2_res_sig2 <- subset(dds2_res_sig, baseMean > 5)
dds2_res_sig2 <- subset(dds2_res_sig, log2FoldChange > 3 | log2FoldChange < -3)
dds2_res_sig2 <- dplyr :: arrange(dds2_res_sig2, desc(dds2_res_sig2$log2FoldChange))
write.csv(dds2_res_sig2,"./11.Data/R_output/230216/cluster2_2_vs_1strictly.csv")

# Subset normalized counts to significant genes
sig_norm_counts_all2 <- counts(dds2,normalized=T)[dds2_res_sig2$Gene, ]
sig_norm_counts_all2 <- sig_norm_counts_all2[,metadata2$X]


# Plot heatmap
pdf("./11.Data/R_output/230216/fig_5_heatmap_highest_genes.pdf", useDingbats = F, height = 5, width = 4)
pheatmap(sig_norm_counts_all2, 
         color = heat_colors,
         #cutree_cols=2,
         cluster_rows = T, cluster_cols=F,
         #clustering_distance_cols="manhattan",
         #clustering_method = "single",
         show_colnames=F,
         show_rownames =F,
         annotation_col = df, 
         fontsize_col= 7,
         scale = "row", main= "Significant different genes")
dev.off()


# Plot count
genelist <- as.data.frame(rownames(dds1))
gene <- dds2_res_sig$Gene
gene <- c("TNFSF12", "TNFRSF25", "TNFRSF12A", "WNT5A", "CTHRC1", "SFRP2")

plotCounts(dds2, gene=gene, intgroup="cluster2", normalized = TRUE, transform = F)

pdf("./11.Data/R_output/230216/fig_6_plotCounts_gene.pdf", useDingbats = F, height = 4, width = 7)

plotCounts <- list()
plot <- list()
for (i in 1:length(gene)) {
  plotCounts[[i]] <- plotCounts(dds2, gene= gene[[i]], intgroup =c("cluster2","SampleID"), returnData=TRUE)
  plot[[i]] <- ggplot(plotCounts[[i]], aes(x=cluster2, y=log(count), color =SampleID)) +
    geom_point(size =4) +
    labs(title = gene[[i]]) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
}
print(plot)


dev.off()








