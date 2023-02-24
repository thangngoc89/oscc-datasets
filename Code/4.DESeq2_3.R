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

#Load data
dds2 <- readRDS("./11.Data/R_output/230216/dds2.RDS")

# Plot count
gene <- c("TNFSF12", "TNFRSF25", "TNFRSF12A", "WNT5A", "CTHRC1", "SFRP2")

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








