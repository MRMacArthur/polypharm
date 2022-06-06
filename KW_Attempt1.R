#Load packages
library(ggplot2)
library(org.Mm.eg.db)
library(limma)
library(edgeR)
library(clusterProfiler)
library(eulerr)
library(RColorBrewer)
library(Glimma)
library(gplots)

# Set ggplot theme style (more setups)
myTheme <- theme(panel.background = element_blank(),
                 axis.line = element_line(color = "black"),
                 text = element_text(color = 'black', size = 12),
                 legend.key=element_blank())

# Load count data
rawData <- read.csv("ausCounts.csv")
geneVec <- rawData[,1]
rownames(rawData) <- rawData[,1]
rawData <- rawData[-1]

# Make group assignments
colnames(rawData) <- sub(".fq.gz.subread.BAM", "", colnames(rawData))
txGroup <- factor(c("Control", "Deprescribed", "Control", "Polypharmacy", 
                    "Control", "Polypharmacy", "Control", "Deprescribed", 
                    "Control", "Deprescribed", "Control", "Polypharmacy", 
                    "Polypharmacy", "Deprescribed", "Polypharmacy", 
                    "Polypharmacy"),
                  levels = c("Control", "Polypharmacy", "Deprescribed"))
color.tx <- c("darkgreen","darkred","cornflowerblue")[txGroup]

# Make DGE object
x <- DGEList(counts = rawData, genes = geneVec)

# Add gene symbol based on ENSEMBL
x$genes$symbol <- mapIds(org.Mm.eg.db, rownames(x),
                         keytype = "ENSEMBL", column = "SYMBOL")
# Add ENTREZ ID based on ENSEMBL
x$genes$entrez <- mapIds(org.Mm.eg.db, rownames(x),
                         keytype = "ENSEMBL", column = "ENTREZID")
x$samples$txGroup <- txGroup

# Filtering lowly transcribed and untranscribed genes
## Sum read counts from all samples, and determine whether it is above (good) or below (bad) the set threshold (i.e., keep gene if 5 or more samples have >5 cpm)
keepS <- rowSums(cpm(x) > 5) >= 5
## See how many genes pass the filter (true)
table(keepS)
## Remove those that fail (in this case, from 55k to 11k) - i.e., remove all untranscribed/lowly transcribed genes
x_5 <- x[keepS, , keep.lib.sizes = F]
dim(x_5)

# Library size
par(mfrow=c(1,1))
barplot(x_5$samples$lib.size/1e06, names = colnames(x_5), col = color.tx, las = 2,ylim=c(0,30))
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")
legend("topright", fill = 2:4, legend=c("Control", "Polypharmacy", "Deprescribed"),  
       col=color.tx, cex=0.8)

# Obtain cpm and log-cpm
xcpm <- cpm(x)
xlcpm <- cpm(x, log=TRUE)

# Density plot 
## Color setup 
nsamples <- ncol(x)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
col = brewer.pal(nsamples, "Dark2")

## For unfiltered
par(mfrow=c(1,2))
plot(density(xlcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main= "A. Raw data", xlab = "Log-cpm")
for (i in 2:nsamples){
  den <- density(xlcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
## For filtered 
xlcpm_5 <- cpm(x_5, log=TRUE)
plot(density(xlcpm_5[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main= "B. Filtered data using >5 cpm cutoff", xlab = "Log-cpm")
for (i in 2:nsamples){
  den <- density(xlcpm_5[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

# TMM normalization 
x_norm_5 = calcNormFactors(x_5, method = "TMM")
x_norm_5$samples$norm.factors

# Boxplot
## Before normalization
par(mfrow=c(1,2))
boxplot(xlcpm_5, las=2, col=color.tx, main="")
abline(h=median(xlcpm_5),col="blue")
title(main="A. Unnormalized data", ylab="Log-cpm")

## After normalization 
xlcpm_norm_5 = cpm(x_norm_5, log=TRUE)
boxplot(xlcpm_norm_5, las=2, col=color.tx, main="")
abline(h=median(xlcpm_norm_5),col="blue")
title(main="B. Normalized data using >5 cpm cutoff", ylab="Log-cpm")

## MDS plot
# Dims 1 and 2
MDS_5 <- plotMDS(xlcpm_norm_5, labels=txGroup, col=color.tx)
title(main="A. Filtered by >5 cpm cutoff")

# Dims 3 and 4
plotMDS(xlcpm_norm_5, labels=txGroup, col=color.tx, dim=c(3,4))
title(main="C. Filtered by >5 cpm cutoff (dim 3/4)")

#Interactive MDS plot
labels <- paste(colnames(xlcpm_norm_5))
glMDSPlot(xlcpm_norm_5, labels=labels, groups=txGroup, col=color.tx, folder="Interactive data viz", launch=T)

# Heatmap
var_genes = apply(xlcpm_norm_5, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
highly_variable_lcpm <- xlcpm_norm_5[select_var,]
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycol <- colorpanel(1000, "blue", "white", "red")
par(mar = c(4.1, 4.4, 4.1, 1.9))
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", 
          main="Top 500 most variable genes across samples",
          ColSideColors=color.tx,scale="row", distfun = '')
legend(x = 1, y = 200, col = color.tx, 
       legend = txGroup, bty = "n", 
       title = "Treatment")

# DGE analysis using limma-voom pipeline
## Design matrix to set up groups
model <- model.matrix(~0 + txGroup)
colnames(model) <- gsub("txGroup", "", colnames(model))
model

## Define pairwise comparisons
groupComps <- makeContrasts("Control vs Polypharmacy" = Control-Polypharmacy,
                            "Control vs Deprescribed" = Control-Deprescribed,
                            "Polypharmacy vs Deprescribed" = Polypharmacy-Deprescribed,
                            levels = model)

## voom plot
v_5 <- voom(x_norm_5, model, plot = T)

## Linear modelling and SA plot
vfit_5 <- lmFit(v_5, model)
vfit_5 <- contrasts.fit(vfit_5, contrasts = groupComps)
efit_5 <- eBayes(vfit_5, trend=F)
plotSA(efit_5)
title(main="SA plot")

## Number of DE genes
summary(decideTests(efit_5))

# Volcano plot
## Pairwise comparison
ConPolyTop <- topTable(efit_5, coef = 1, sort.by = "P", n = Inf)
ConPolyTop$sig <- F; ConPolyTop[ConPolyTop$P.Value < 0.05, ]$sig <- T

PolyDepTop <- topTable(efit_5, coef = 2, sort.by = "P", n = Inf)
PolyDepTop$sig <- F; PolyDepTop[PolyDepTop$P.Value < 0.05, ]$sig <- T

ConDepTop <- topTable(efit_5, coef = 3, sort.by = "P", n = Inf)
ConDepTop$sig <- F; ConDepTop[ConDepTop$P.Value < 0.05, ]$sig <- T

ggplot(ConPolyTop, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.7) + labs(x = "Log2 fold change", 
                                 y = "-log10 adjusted p value",
                                 color = "Significant") +
  myTheme

ggplot(PolyDepTop, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.7) + labs(x = "Log2 fold change", 
                                 y = "-log10 adjusted p value",
                                 color = "Significant") +
  myTheme

ggplot(ConDepTop, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.7) + labs(x = "Log2 fold change", 
                                 y = "-log10 adjusted p value",
                                 color = "Significant") +
  myTheme

table(ConPolyTop$sig)
table(PolyDepTop$sig)
table(ConDepTop$sig)

# Taking a look at the most DE gene from each pairwise comparison
## First gene from the topTable (ConDepTop) is ENSMUSG00000093930 (Hmgcs1):
par(mfrow=c(1,3))
stripchart(v_5$E["ENSMUSG00000093930",]~txGroup,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,method="jitter",
           ylab="Normalized log2 transcript level",main="Hmgcs1")

## First gene from the topTable (ConPolyTop) is ENSMUSG00000031355 (Arhgap6):
stripchart(v_5$E["ENSMUSG00000031355",]~txGroup,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,method="jitter",
           ylab="Normalized log2 transcript level",main="Arhgap6")

## First gene from the topTable (PolyDepTop) is ENSMUSG00000024807 (Syvn1):
stripchart(v_5$E["ENSMUSG00000024807",]~txGroup,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,method="jitter",
           ylab="Normalized log2 transcript level",main="Syvn1")

# Interactive volcano (for ctrl vs poly)
tfit_5 = treat(vfit_5)
dT = decideTests(tfit_5, method = NULL, trend = F, adjust.method = NULL, p.value = 0.05)
glXYPlot(x=efit_5$coefficients[,1], y=-log(efit_5$p.value[,1]),
         xlab="logFC", ylab="-log(p-value)", main="Ctrl v Poly", 
         counts=v_5$E, groups=txGroup, status=dT[,1], 
         side.main="ENTREZID",folder="Interactive data viz")
