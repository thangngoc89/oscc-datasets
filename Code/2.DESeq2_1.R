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

# Check biomart
listMarts()
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl"
#                         ,host = 'jul2022.archive.ensembl.org' # verison 107
                         ) 

# Check databset
listDatasets(mart)
# Check attributes
listAttributes(mart)
# Check filters
listFilters(mart)

tx2gene <- biomaRt::getBM(attributes = c("ensembl_transcript_id_version", "external_gene_name"), mart = mart)

tx2gene <- dplyr::rename(tx2gene, TXNAME = ensembl_transcript_id_version,
                          GENEID = external_gene_name)

#Save and load
#saveRDS(tx2gene,"./11.Data/R_output/tx2gene.RDS")
tx2gene <- readRDS("./11.Data/R_output/tx2gene.RDS")
genesymbol <- as.data.frame(tx2gene$GENEID)
genesymbol$`tx2gene$GENEID` <- ifelse(genesymbol$`tx2gene$GENEID` == '', NA, genesymbol$`tx2gene$GENEID`)
genesymbol <- na.omit(genesymbol)
genesymbol <- unique(genesymbol)
genesymbol <- dplyr :: arrange(genesymbol, genesymbol$`tx2gene$GENEID`)
write.csv(genesymbol,"./11.Data/R_output/230216/genesymbol.csv")

# Read metadata
metadata <- read.table("./11.Data/Metadata/Metadata.txt", sep= "\t",header = TRUE)

# Load Kallisto files
name_files <- metadata$SampleID
kallisto_files <- file.path("./12.Dataforweb/2. Cơ sở dữ liệu/3. Processed files", name_files, "abundance.h5")
names(kallisto_files) <- name_files
all(file.exists(kallisto_files))

# Replace empty GeneID info
tx2gene$GENEID <- ifelse(tx2gene$GENEID == '', NA, tx2gene$GENEID)
tx2gene <- tx2gene[complete.cases(tx2gene), ]

# Import Kallisto
txi.kallisto <- tximport(kallisto_files, type = "kallisto", tx2gene = tx2gene, varReduce= T, txOut = F)
head(txi.kallisto$counts)
names(txi.kallisto)

#Save and load
saveRDS(txi.kallisto,"./11.Data/R_output/230216/txi.kallisto.RDS")
txi.kallisto <- readRDS("./11.Data/R_output/230216/txi.kallisto.RDS")


###
# Import with SampleID
ddsMat <- DESeqDataSetFromTximport(txi.kallisto, metadata, ~ SampleID)

# Pre-filtering the dataset
nrow(ddsMat)
ddsMat <- ddsMat[rowSums(counts(ddsMat)) > 1, ] # keep genes with total counts >1
nrow(ddsMat)

# Determine the size factors to use for normalization
dds1 <- estimateSizeFactors(ddsMat)
sizeFactors(dds1)

# Transform
vsd <- varianceStabilizingTransformation(ddsMat)

# PCA plot
pdf("./11.Data/R_output/230216/fig_1_PCA_pre.pdf", useDingbats = F, height = 3, width =4)
plotPCA(vsd, "SampleID")
plotPCA(vsd, "CaseID")
plotPCA(vsd, "PatientID")
plotPCA(vsd, "seq_date")
dev.off()

assay(vsd) <- limma :: removeBatchEffect(assay(vsd), dds1$seq_date)

pdf("./11.Data/R_output/230216/fig_1_PCA_pos.pdf", useDingbats = F, height = 3, width =4)
plotPCA(vsd, "SampleID")
plotPCA(vsd, "CaseID")
plotPCA(vsd, "PatientID")
plotPCA(vsd, "seq_date")
dev.off()

# Sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)

# Cluster
plotCluster <- list()
pdf("./11.Data/R_output/230216/fig_2_sample_distance.pdf", useDingbats = F, height = 5, width = 5.5)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, show_colnames=T,
         show_rownames =T)
for (i in 2:5) {
  plotCluster[[i]] <- pheatmap(t(sampleDistMatrix), kmeans_k=i,
                               #clustering_distance_rows = sampleDists,
                               #clustering_distance_cols = sampleDists,
                               col = colors, show_colnames=T,
                               show_rownames =T)
}
print(plotCluster[[i]])
dev.off()


###
# ConsensusClusterPlus
# Prepare data 
# Get data matrix

d <- counts(dds1,normalized=T)
write.table(d,"./11.Data/R_output/GSEA/Expression/norm_counts_forGSEA.txt", sep="\t", quote=FALSE, row.names=T)


d <- counts(dds1,normalized=T)
# Top 5,000 most variable genes
mads <- apply(d,1,mad)
d <- d[rev(order(mads))[1:5000],]
# Center to median
#d <- sweep(d,1, apply(d,1,median,na.rm=T))


# ConsensusClusterPlus
title <- tempdir()
results <- ConsensusClusterPlus(d,maxK=5,reps=1000,pItem=0.8,pFeature=1,
                                 title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf")

icl <- calcICL(results,title=title,plot="pdf")
icl[["clusterConsensus"]]


pdf("./11.Data/R_output/230216/fig_3_ConsensusCluster.pdf", useDingbats = F, height = 5, width = 5)
ConsensusClusterPlus(d,maxK=5,reps=1000,pItem=0.8,pFeature=1,
                     title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279)
icl <- calcICL(results,title=title)

dev.off()

results

# Add cluster to metadata
metadata$cluster2 <- results[[2]][["consensusClass"]]
metadata$cluster3 <- results[[3]][["consensusClass"]]
metadata$cluster4 <- results[[4]][["consensusClass"]]
metadata$cluster5 <- results[[5]][["consensusClass"]]

write.csv(metadata,"./11.Data/R_output/230216/metadata_withclusters.csv")









