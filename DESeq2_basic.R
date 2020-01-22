
#There were some issues with your file. The first is that it had DOS like new line coding which R doesn't like.
#THe second issues was there were a bunch of '#N/A' s. this is problematice because both # and / are not parsed by a lof of algorithms.
#I ran this to convert those to NA
# cat gene_counts_annotation.txt | tr -d '\r' | sed 's/#N\/A/NA/g' > gene_counts_annotation_clean.txt
#I also had to remove the descriptions line because there are too many special characters

#read in counts where each gene is a row and each column is a sample
counts <- read.table(file="gene_counts_annotation_clean_nodescriptions.txt", sep="\t", header=T, quote = "")

#We'll set aside the query = descriptions for now. We can add them back later:
#if you're new to R this is called bracket subsetting
annotations <- counts[1:2]

library("plyr")

#this will set aside the no blast transcripts
noblast<- counts[which(is.na(counts$accession_number)),]
write.table(noblast, file="RawCounts_NoBlast.txt", sep="\t", quote =F, row.names=F)

#remove those same no blast hits from the counts for now (we can add them back)
counts <- counts[-c(which(is.na(counts$accession_number))),]

#remove the TRINITY column and description so that ddply will work properly
counts <-counts[-c(1)]
#Sum counts by top blast hit
#add your column names on the lines below

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


#some formating
rownames(counts_summed) = counts_summed[,(1)]
counts_summed <- counts_summed[,-(1)]
#coerce to matrix
m <- data.matrix(counts_summed)

#drop low counts
#change the numbers to reflect your cutoff. Below is 5 counts in more than 12 samples
m <- m[rowSums(m > 5) >=5,]
#convert nas to 0
m[is.na(m)] <- 0
#round counts
m <- apply(m, 1:2, round)

#get the meta data for DEseq in
meta <- as.data.frame(read.table(file="meta.txt", sep="\t", header=T))

rownames(meta)=colnames(m)
colData <- meta
countData <- m
rownames(countData) <- row.names(countData)
library(DESeq2)

#sanity checks 
all(rownames(colData) %in% colnames(countData))
all(rownames(colData) == colnames(countData))

#choose a design that suites you
dds <-  DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~Condition)

## Estimate library size and dispersion
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds, main="DESeq: Per-gene dispersion estimates")

counts.dds <- counts(dds, normalized = TRUE)

## Principal components biplot on variance stabilized data, color-coded by condition-librarytype
vsdds <- varianceStabilizingTransformation(dds)

write.table(assay(vsdds), file="variance_stabilized.txt", sep="\t", quote=F)

z <- assay(vsdds)
z <- as.data.frame(z)

#PCA
#code the conditions
meta$Condition
group <- meta$Condition
#code the replicates
grouprep <- factor(c(21,21,21,21,21,21,21,22,22,22,22,22,22,22,23,23,23,23,23,23,23))
transpose <- t(z) #this puts samples as rows, genes as columns 
transpose_df <- as.data.frame(transpose)
pca.data <- prcomp(transpose_df) #this is the pca function that I like to use. For this data set, I have 9772 genes. I have to make this explicit here, because I added the IDs in the line above (so my data frame has 9773 “samples”)
scores = as.data.frame(pca.data$x) 
summary(pca.data) #this will give you the proportion of variance explained by PC1, PC2, etc. usually we only plot PC1 and PC2 since these explain most of the variance.


p <- ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores), colour=factor(group), shape=factor(grouprep))) + geom_point(size=6) + scale_fill_hue(l=40) + coord_fixed(ratio=1, xlim=c(-200, 300), ylim=c(-200, 200)) 
p
ggsave(p, filename = "PCA.pdf", height=8, width=12)