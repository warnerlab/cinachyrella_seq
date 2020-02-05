library(reshape)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(plyr)
library(readr)

setwd("~/Desktop/test_2/")

counts <- read.csv(file="gene_counts_annotation.csv", header=T, sep=";")
counts <-counts[-c(1,3)]

counts_summed = ddply(counts, "accession_number", summarise,
                      X1C = sum(X1C),
                      X1O1 = sum(X1O1),
                      X1D1 = sum(X1D1),
                      X1OD1 = sum(X1OD1), 
                      X1O24 = sum(X1O24),
                      X1D24 = sum(X1D24),
                      X1OD24 = sum(X1OD24),
                      X2C = sum(X2C),
                      X2O1 = sum(X2O1),
                      X2D1 = sum(X2D1),
                      X2OD1 = sum(X2OD1), 
                      X2O24 = sum(X2O24),
                      X2D24 = sum(X2D24),
                      X2OD24 = sum(X2OD24),
                      X3C = sum(X3C),
                      X3O1 = sum(X3O1),
                      X3D1 = sum(X3D1),
                      X3OD1 = sum(X3OD1), 
                      X3O24 = sum(X3O24),
                      X3D24 = sum(X3D24),
                      X3OD24 = sum(X3OD24))

write.table(counts_summed, file="counts_summed.txt", sep="\t", quote=F)

rownames(counts_summed) = counts_summed[,(1)]
counts_summed <- counts_summed[,-(1)]

m <- data.matrix(counts_summed)

meta <- as.data.frame(read.table(file="meta.txt", sep="\t", header=T))

rownames(meta)=colnames(m)
colData <- meta
countData <- m
rownames(countData) <- row.names(countData)

all(names(data) %in% rownames(meta))
all(names(data) == rownames(meta))

dds <- DESeqDataSetFromMatrix(countData=m, colData=meta, design =~Condition)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType=c("local"))
plotDispEsts(dds, main="Per-gene dispersion estimates")

normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)

vsdds <- varianceStabilizingTransformation(dds)

write.table(assay(vsdds), file="variance_stabilized.txt", sep="\t", quote=F)

z <- assay(vsdds)
z <- as.data.frame(z)

#PCA plot creation
#code the conditions
meta$Condition
treatment<- meta$Condition
#code the replicates
experiment <- factor(c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3))
transpose <- t(z) #this puts samples as rows, genes as columns 
transpose_df <- as.data.frame(transpose)
pca.data <- prcomp(transpose_df) #this is the pca function that I like to use. For this data set, I have 9772 genes. I have to make this explicit here, because I added the IDs in the line above (so my data frame has 9773 “samples”)
scores = as.data.frame(pca.data$x) 
summary(pca.data) #this will give you the proportion of variance explained by PC1, PC2, etc. usually we only plot PC1 and PC2 since these explain most of the variance.


p <- ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores), colour=treatment, shape=experiment)) + scale_shape_manual(values=c(16,17,18)) + geom_point(size=6) + scale_fill_hue(l=40) + coord_fixed(ratio=1, xlim=c(-250, 450), ylim=c(-180, 250)) 
PCA <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
PCA
ggsave(PCA, filename = "PCA.pdf", height=8, width=10)

# Other figures

res <- DESeq(dds, test=c("Wald"), fitType=c("local"))
resultsNames(res)
rdds <- results(res)
summary(rdds)



write.table(results(res, contrast=c("Condition", "oil_1", "control")), file="control_vs_O1_test.txt", sep="\t", quote=F)
write.table(results(res, contrast=c("Condition", "oil_24", "control")), file="control_vs_O24_test.txt", sep="\t", quote=F)
write.table(results(res, contrast=c("Condition", "dispersant_1", "control")), file="control_vs_D1_test.txt", sep="\t", quote=F)
write.table(results(res, contrast=c("Condition", "dispersant_24", "control")), file="control_vs_D24_test.txt", sep="\t", quote=F)
write.table(results(res, contrast=c("Condition", "oil_dispersant_1", "control")), file="control_vs_OD1_test.txt", sep="\t", quote=F)
write.table(results(res, contrast=c("Condition", "oil_dispersant_24", "control")), file="control_vs_OD24_test.txt", sep="\t", quote=F)
write.table(results(res, contrast=c("Condition", "oil_24", "oil_1")), file="O1_vs_O24_test.txt", sep="\t", quote=F)
write.table(results(res, contrast=c("Condition", "dispersant_24", "dispersant_1")), file="D1_vs_D24_test.txt", sep="\t", quote=F)
write.table(results(res, contrast=c("Condition", "oil_dispersant_24", "oil_dispersant_1")), file="OD1_vs_OD24_test.txt", sep="\t", quote=F)



C_O1 <- read.table("control_vs_O1_test.txt", header=T, sep="\t")
C_D1 <- read.table("control_vs_D1_test.txt", header=T, sep="\t")
C_OD1 <- read.table("control_vs_OD1_test.txt", header=T, sep="\t")
C_O24 <- read.table("control_vs_O24_test.txt", header=T, sep="\t")
C_D24 <- read.table("control_vs_D24_test.txt", header=T, sep="\t")
C_OD24 <- read.table("control_vs_OD24_test.txt", header=T, sep="\t")
O1_O24 <- read.table("O1_vs_O24_test.txt", header=T, sep="\t")
D1_D24 <- read.table("D1_vs_D24_test.txt", header=T, sep="\t")
OD1_OD24 <- read.table("OD1_vs_OD24_test.txt", header=T, sep="\t")

# Volcano plot creation

C_O1$sig <- as.factor(abs(C_O1$log2FoldChange) > 2 & C_O1$padj < 0.05) 

#this didn't work
#test_table <- filter(C_O1, sig == "TRUE")

#you can also use square bracket subsetting:
test_table <- C_O1[which(C_O1$sig=='TRUE'),]

#use dim to show the both dimensions
dim(test_table)

volc_C_O1 <- ggplot(C_O1) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=sig)) +
  geom_vline(xintercept=-2) +
  geom_vline(xintercept=2) +
  geom_hline(yintercept=1.3) +
  ggtitle("Control vs Oil 1h") +
  xlab("Log2 Fold Change") + 
  ylab("-log10(padj)") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


a <- volc_C_O1 + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=30), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),panel.border = element_rect(colour = "black", fill=NA, size=0.8), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none" )
a
ggsave(a, filename = "Control_vs_O1.png", height=8, width=10)  


C_D1$sig1 <- as.factor(abs(C_D1$log2FoldChange) > 2 & C_D1$padj < 0.05)  

test_table <- filter(C_D1, sig1 == "TRUE")
length(test_table)

volc_C_D1 <- ggplot(C_D1) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=sig1)) +
  geom_vline(xintercept=-2) +
  geom_vline(xintercept=2) +
  geom_hline(yintercept=1.3) +
  ggtitle("Control vs Dispersant 1h") +
  xlab("Log2 Fold Change") + 
  ylab("-log10(padj)") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


b <- volc_C_D1 + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=30), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),panel.border = element_rect(colour = "black", fill=NA, size=0.8), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")
b
ggsave(b, filename = "Control_vs_D1.png", height=8, width=10)  


C_OD1$sig2 <- as.factor(abs(C_OD1$log2FoldChange) > 2 & C_OD1$padj < 0.05)  

test_table <- filter(C_OD1, sig2 == "TRUE")
length(test_table)

volc_C_OD1 <- ggplot(C_OD1) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=sig2)) +
  geom_vline(xintercept=-2) +
  geom_vline(xintercept=2) +
  geom_hline(yintercept=1.3) +
  ggtitle("Control vs Oil:Dispersant 1h") +
  xlab("Log2 Fold Change") + 
  ylab("-log10(padj)") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


c <- volc_C_OD1 + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=30), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),panel.border = element_rect(colour = "black", fill=NA, size=0.8), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")
c
ggsave(c, filename = "Control_vs_OD1.png", height=8, width=10)  


C_O24$sig3 <- as.factor(abs(C_O24$log2FoldChange) > 2 & C_O24$padj < 0.05)  

test_table <- filter(C_O24, sig3 == "TRUE")
length(test_table)

volc_C_O24 <- ggplot(C_O24) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=sig3)) +
  geom_vline(xintercept=-2) +
  geom_vline(xintercept=2) +
  geom_hline(yintercept=1.3) +
  ggtitle("Control vs Oil 24h") +
  xlab("Log2 Fold Change") + 
  ylab("-log10(padj)") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


d <- volc_C_O24 + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=30), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),panel.border = element_rect(colour = "black", fill=NA, size=0.8), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")
d
ggsave(d, filename = "Control_vs_O24.png", height=8, width=10)  


C_D24$sig4 <- as.factor(abs(C_D24$log2FoldChange) > 2 & C_D24$padj < 0.05)  

test_table <- filter(C_D24, sig4 == "TRUE")
length(test_table)

volc_C_D24 <- ggplot(C_D24) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=sig4)) +
  geom_vline(xintercept=-2) +
  geom_vline(xintercept=2) +
  geom_hline(yintercept=1.3) +
  ggtitle("Dontrol vs Dispersant 24h") +
  xlab("Log2 Fold Change") + 
  ylab("-log10(padj)") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


e <- volc_C_D24 + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=30), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),panel.border = element_rect(colour = "black", fill=NA, size=0.8), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")
e
ggsave(e, filename = "Control_vs_D24.png", height=8, width=10)  


C_OD24$sig5 <- as.factor(abs(C_OD24$log2FoldChange) > 2 & C_OD24$padj < 0.05)  

test_table <- filter(C_OD24, sig5 == "TRUE")
length(test_table)

volc_C_OD24 <- ggplot(C_OD24) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=sig5)) +
  geom_vline(xintercept=-2) +
  geom_vline(xintercept=2) +
  geom_hline(yintercept=1.3) +
  ggtitle("control vs Oil:Dispersant 24h") +
  xlab("Log2 Fold Change") + 
  ylab("-log10(padj)") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


f <- volc_C_OD24 + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=30), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),panel.border = element_rect(colour = "black", fill=NA, size=0.8), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")
f
ggsave(f, filename = "Control_vs_OD24.png", height=8, width=10)  


O1_O24$sig6 <- as.factor(abs(O1_O24$log2FoldChange) > 2 & O1_O24$padj < 0.05)

test_table <- filter(O1_O24, sig6 == "TRUE")
length(test_table)

volc_O1_O24 <- ggplot(O1_O24) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=sig6)) +
  geom_vline(xintercept=-2) +
  geom_vline(xintercept=2) +
  geom_hline(yintercept=1.3) +
  ggtitle("Oil 1h vs Oil 24h") +
  xlab("Log2 Fold Change") + 
  ylab("-log10(padj)") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


g <- volc_O1_O24 + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=30), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),panel.border = element_rect(colour = "black", fill=NA, size=0.8), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")
g
ggsave(g, filename = "O1_vs_O24.png", height=8, width=10)  


D1_D24$sig7 <- as.factor(abs(D1_D24$log2FoldChange) > 2 & D1_D24$padj < 0.05)  

test_table <- filter(D1_D24, sig7 == "TRUE")
length(test_table)

volc_D1_D24 <- ggplot(D1_D24) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=sig7)) +
  geom_vline(xintercept=-2) +
  geom_vline(xintercept=2) +
  geom_hline(yintercept=1.3) +
  ggtitle("Dispersant 1h vs Dispersant 24h") +
  xlab("Log2 Fold Change") + 
  ylab("-log10(padj)") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


h <- volc_D1_D24 + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=30), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),panel.border = element_rect(colour = "black", fill=NA, size=0.8), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")
h
ggsave(h, filename = "D1_vs_D24.png", height=8, width=10)  


OD1_OD24$sig8 <- as.factor(abs(OD1_OD24$log2FoldChange) > 2 & OD1_OD24$padj < 0.05)

test_table <- filter(OD1_OD24, sig8 == "TRUE")
length(test_table)

volc_OD1_OD24 <- ggplot(OD1_OD24) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=sig8)) +
  geom_vline(xintercept=-2) +
  geom_vline(xintercept=2) +
  geom_hline(yintercept=1.3) +
  ggtitle("Oil:Dispersant 1h vs Oil:Dispersant 24h") +
  xlab("Log2 Fold Change") + 
  ylab("-log10(padj)") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


i <- volc_OD1_OD24 + theme_bw() +theme(plot.title=element_text(hjust=0.5, size=30), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),panel.border = element_rect(colour = "black", fill=NA, size=0.8), panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")
i
ggsave(i, filename = "OD1_vs_OD24.png", height=8, width=10)  





#heatmap creation
threshold <- subset(rdds, abs(log2FoldChange) > 2 & padj < 0.05) 
norm <- normalized_counts[rownames(threshold),]       

annotation <- data.frame(sampletype=meta[,'Condition'], 
                         row.names=rownames(meta))       

heat_colors <- brewer.pal(6, "YlOrRd")
#if other colors wanted can set own color palette using the colors() function 
# colors on the palette available at 
# https://www.r-bloggers.com/creating-color-palettes-in-r/
# my_cols <- colors()[c(256,154,134)]

pheatmap(rdds, color = heat_colors, cluster_rows = T, show_rownames=F,
         annotation= annotation, border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)







