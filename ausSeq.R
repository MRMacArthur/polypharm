library(ggplot2)
library(org.Mm.eg.db)
library(limma)
library(edgeR)
library(clusterProfiler)
library(eulerr)

myTheme <- theme(panel.background = element_blank(),
                 axis.line = element_line(color = "black"),
                 text = element_text(color = 'black', size = 12),
                 legend.key=element_blank())

# Load count data
rawData <- read.csv("ausCounts.csv")
geneVec <- rawData[,1]
rownames(rawData) <- rawData[,1]
rawData <- rawData[-1]

colnames(rawData)

# Make group assignments
txGroup <- factor(c("Con", "DBIdep", "Con", "DBI", "Con", 
             "DBI", "Con", "DBIdep", "Con", "DBIdep", 
             "Con", "DBI", "DBI", "DBIdep", "DBI", 
             "DBI"),
             levels = c("Con", "DBI", "DBIdep"))

# Make DGE object
xS <- DGEList(counts = rawData, genes = geneVec)
xS$genes$symbol <- mapIds(org.Mm.eg.db, rownames(xS),
                          keytype = "ENSEMBL", column = "SYMBOL")
xS$genes$entrez <- mapIds(org.Mm.eg.db, rownames(xS),
                          keytype = "ENSEMBL", column = "ENTREZID")

# Filter low expression genes
keepS <- rowSums(cpm(xS) > 5) >= 5
table(keepS)
xS <- xS[keepS, , keep.lib.sizes = F]
xS <- calcNormFactors(xS)

#QC Plots
barplot(xS$samples$lib.size, names = colnames(xS), las = 2)
title("Library sizes")
logcounts <- cpm(xS,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

scatPlot <- plotMDS(xS)
scatMat <- data.frame(cbind(as.numeric(scatPlot$x), as.numeric(scatPlot$y)))
colnames(scatMat) <- c("x", "y")

scatMat$txGroup <- txGroup
ggplot(scatMat, aes(x = x, y = y, color = txGroup)) +
  geom_point(size = 3) + 
  stat_ellipse(level = 0.8) +
  myTheme + labs(x = "Dim1", y = "Dim2", color = "Tx")

#generate models
modelMatInt <- model.matrix(~ 0 + txGroup)
colnames(modelMatInt) <- gsub("txGroup", "", colnames(modelMatInt))

groupComps <- makeContrasts(Con_DBI = Con-DBI,
                            Con_DBIdep = Con-DBIdep,
                            DBI_DBIdep = DBI-DBIdep,
                            levels = modelMatInt)
v <- voom(xS, modelMatInt, plot = T)
vfit <- lmFit(v, modelMatInt)
vfit <- contrasts.fit(vfit, contrasts = groupComps)
efit <- eBayes(vfit)
plotSA(efit)

summary(decideTests(efit))

ConDBITop <- topTable(efit, coef = 1, sort.by = "P", n = Inf)
ConDBITop$sig <- F; ConDBITop[ConDBITop$P.Value < 0.05, ]$sig <- T
DBI_DBIdepTop <- topTable(efit, coef = 2, sort.by = "P", n = Inf)
DBI_DBIdepTop$sig <- F; DBI_DBIdepTop[DBI_DBIdepTop$P.Value < 0.05, ]$sig <- T
ConDBIdepTop <- topTable(efit, coef = 3, sort.by = "P", n = Inf)
ConDBIdepTop$sig <- F; ConDBIdepTop[ConDBIdepTop$P.Value < 0.05, ]$sig <- T

ggplot(ConDBITop, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.7) + labs(x = "Log2 fold change", 
                                 y = "-log10 adjusted p value",
                                 color = "Significant") +
  myTheme

ggplot(DBI_DBIdepTop, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.7) + labs(x = "Log2 fold change", 
                                 y = "-log10 adjusted p value",
                                 color = "Significant") +
  myTheme

ggplot(ConDBIdepTop, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.7) + labs(x = "Log2 fold change", 
                                 y = "-log10 adjusted p value",
                                 color = "Significant") +
  myTheme

table(ConDBITop$sig)
table(DBI_DBIdepTop$sig)
table(ConDBIdepTop$sig)

ConDBISig <- ConDBITop[ConDBITop$P.Value < 0.05,]$entrez
ConDBISigUp <- ConDBITop[ConDBITop$P.Value < 0.05 & ConDBITop$logFC > 0,]$entrez
ConDBISigDown <- ConDBITop[ConDBITop$P.Value < 0.05 & ConDBITop$logFC < 0,]$entrez
ConDBISigClust <- ConDBITop[ConDBITop$logFC < -1.2,]$entrez

DBI_DBIdepSig <- DBI_DBIdepTop[DBI_DBIdepTop$P.Value < 0.05,]$entrez
DBI_DBIdepSigUp <- DBI_DBIdepTop[DBI_DBIdepTop$P.Value < 0.05 & DBI_DBIdepTop$logFC > 0,]$entrez
DBI_DBIdepSigDown <- DBI_DBIdepTop[DBI_DBIdepTop$P.Value < 0.05 & DBI_DBIdepTop$logFC < 0,]$entrez

ConDBIdepSig <- ConDBIdepTop[ConDBIdepTop$P.Value < 0.05,]$entrez
ConDBIdepSigUp <- ConDBIdepTop[ConDBIdepTop$P.Value < 0.05 & ConDBIdepTop$logFC > 0,]$entrez
ConDBIdepSigDown <- ConDBIdepTop[ConDBIdepTop$P.Value < 0.05 & ConDBIdepTop$logFC < 0,]$entrez

DBIup_DBIdepDown <- intersect(ConDBISigUp, DBI_DBIdepSigUp)
DBIup_depUp <- intersect(ConDBISigUp, DBI_DBIdepSigDown)


ConDBIPath <- enrichGO(ConDBISig,
                       universe = xS$genes$entrez,
                       OrgDb = org.Mm.eg.db,
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = T,
                       pool = T)

ConDBIPathUp <- enrichGO(ConDBISigUp,
                    universe = xS$genes$entrez,
                    OrgDb = org.Mm.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = T,
                    pool = T)

ConDBIPathDown <- enrichGO(ConDBISigDown,
                       universe = xS$genes$entrez,
                       OrgDb = org.Mm.eg.db,
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,
                       readable = T,
                       pool = T)

ConDBIPathClust <- enrichGO(ConDBISigClust,
                         universe = xS$genes$entrez,
                         OrgDb = org.Mm.eg.db,
                         ont = "ALL",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = T,
                         pool = T)

DBI_DBIdepPath <- enrichGO(DBI_DBIdepSig,
                       universe = xS$genes$entrez,
                       OrgDb = org.Mm.eg.db,
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = T,
                       pool = T)

DBI_DBIdepPathUp <- enrichGO(DBI_DBIdepSigUp,
                    universe = xS$genes$entrez,
                    OrgDb = org.Mm.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = T,
                    pool = T)

DBI_DBIdepPathDown <- enrichGO(DBI_DBIdepSigDown,
                       universe = xS$genes$entrez,
                       OrgDb = org.Mm.eg.db,
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = T,
                       pool = T)

ConDBIdepPath <- enrichGO(ConDBIdepSig,
                       universe = xS$genes$entrez,
                       OrgDb = org.Mm.eg.db,
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = T,
                       pool = T)

ConDBIdepPathUp <- enrichGO(ConDBIdepSigUp,
                     universe = xS$genes$entrez,
                     OrgDb = org.Mm.eg.db,
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05,
                     readable = T,
                     pool = T)

ConDBIdepPathDown <- enrichGO(ConDBIdepSigDown,
                     universe = xS$genes$entrez,
                     OrgDb = org.Mm.eg.db,
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05,
                     readable = T,
                     pool = T)



DBIup_DepDown <- enrichGO(DBIup_DBIdepDown,
                       universe = xS$genes$entrez,
                       OrgDb = org.Mm.eg.db,
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = T,
                       pool = T)

dotplot(DBIup_DepDown)

dotplot(ConDBIPathUp, showCategory = 20)
dotplot(ConDBIPathDown)
dotplot(ConDBIPathClust, showCategory  = 20)

dotplot(DBI_DBIdepPathUp)
dotplot(DBI_DBIdepPathDown, showCategory = 20)

dotplot(ConDBIdepPath)
dotplot(ConDBIdepPathUp)
dotplot(ConDBIdepPathDown)

cpmX <- cpm(xS,log=F)
cpmX <- data.frame(t(cpmX))
cpmX$txGroup <- txGroup

ggplot(cpmX, aes(x = txGroup, y = ENSMUSG00000026812)) +
  geom_boxplot() + 
  labs(x = "Tx", y = "Mvd (cpm)") + myTheme

#WGCNA
library(WGCNA)

cpmOut <- cpm(xS, log = F, normalized.lib.sizes = T)
cpmT <- t(cpmOut)
rownames(cpmT) <- sub("_1.*", "", rownames(cpmT))

metaDat <- read.csv("metadata.csv")
rownames(metaDat) <- metaDat$ID
metaDat <- metaDat[,3:42]

cpmT <- cpmT[rownames(cpmT) %in% rownames(metaDat),]
metaDat <- metaDat[rownames(metaDat) %in% rownames(cpmT), ]

sampleTree <- hclust(dist(cpmT), method = "average")
traitColors <- numbers2colors(metaDat, signed = F)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(metaDat),
                    main = "")

options(stringsAsFactors = F)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(cpmT, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5);par(mfrow = c(1,2));cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");# this line corresponds to using an R^2 cut-off of habline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", 
     type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(cpmT, power = 7,
                       TOMType = "unsigned", 
                       minModuleSize = 30,
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ausTom",
                       verbose = 3)

mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",dendroLabels = FALSE, 
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                    main = "")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
#save(MEs, moduleLabels, moduleColors, geneTree,
#     file = "sarcoModules.RData")

nGenes <- ncol(cpmT)
nSamples <- nrow(cpmT)

MEsO <- moduleEigengenes(cpmT, moduleColors)$eigengenes
MEs <- orderMEs(MEsO)

moduleTraitCor <- cor(MEs, metaDat, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nSamples)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitP, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(metaDat),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = T,
               cex.text = 0.5,
               zlim = c(-1,1))

FI = as.data.frame(metaDat$FI)
names(FI) = "FI"
modNamesFI <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(cpmT, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("MM", modNamesFI, sep="")
names(MMPvalue) <- paste("p.MM", modNamesFI, sep="")
geneTraitSigFI <- as.data.frame(cor(cpmT, FI, use = "p"))
FIPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSigFI), nSamples))
names(geneTraitSigFI) <- paste("FI.", names(FI), sep="")
names(FIPvalue) = paste("p.FI", names(FI), sep="")

module = "cyan"
column = match(module, modNamesFI)
moduleGenes <- moduleColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSigFI[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", names(FI)),
                   main = "", cex.main=1.2, cex.lab=1.2, cex.axis=1.2)

cyanGenes <- rownames(geneModuleMembership[moduleColors == "cyan",])
salmonGenes <- rownames(geneModuleMembership[moduleColors == "salmon",])
magentaGenes <- rownames(geneModuleMembership[moduleColors == "magenta",])
blackGenes <- rownames(geneModuleMembership[moduleColors == "black",])
lightcyanGenes <- rownames(geneModuleMembership[moduleGenes,])
grey60Genes <- rownames(geneModuleMembership[moduleGenes,])
purpleGenes <- rownames(geneModuleMembership[moduleGenes,])
yellowGenes <- rownames(geneModuleMembership[moduleGenes,])



cyanEnt <- mapIds(org.Mm.eg.db, cyanGenes,
                     keytype = "ENSEMBL", column = "ENTREZID")

salmonEnt <- mapIds(org.Mm.eg.db, salmonGenes,
                  keytype = "ENSEMBL", column = "ENTREZID")
salmonSym <- mapIds(org.Mm.eg.db, salmonGenes,
                    keytype = "ENSEMBL", column = "SYMBOL")

magentaEnt <- mapIds(org.Mm.eg.db, magentaGenes,
                     keytype = "ENSEMBL", column = "ENTREZID")
magentaSym <- mapIds(org.Mm.eg.db, magentaGenes,
                     keytype = "ENSEMBL", column = "SYMBOL")

blackEnt <- mapIds(org.Mm.eg.db, blackGenes,
                   keytype = "ENSEMBL", column = "ENTREZID")
blackSym <- mapIds(org.Mm.eg.db, blackGenes,
                   keytype = "ENSEMBL", column = "SYMBOL")

yellowEnt <- mapIds(org.Mm.eg.db, yellowGenes,
                    keytype = "ENSEMBL", column = "ENTREZID")

purpleEnt <- mapIds(org.Mm.eg.db, purpleGenes,
                    keytype = "ENSEMBL", column = "ENTREZID")
purpleSym <- mapIds(org.Mm.eg.db, purpleGenes,
                    keytype = "ENSEMBL", column = "SYMBOL")

grey60Ent <- mapIds(org.Mm.eg.db, grey60Genes,
                    keytype = "ENSEMBL", column = "ENTREZID")
grey60Sym <- mapIds(org.Mm.eg.db, grey60Genes,
                    keytype = "ENSEMBL", column = "SYMBOL")

lightcyanEnt <- mapIds(org.Mm.eg.db, lightcyanGenes,
                       keytype = "ENSEMBL", column = "ENTREZID")
lightcyanSym <- mapIds(org.Mm.eg.db, lightcyanGenes,
                       keytype = "ENSEMBL", column = "SYMBOL")

allGenes <- mapIds(org.Mm.eg.db, colnames(cpmT),
                   keytype = "ENSEMBL", column = "ENTREZID")

modEnrich <- enrichGO(gene = cyanEnt[!is.na(cyanEnt)],
                      universe = ConDBIdepTop$entrez,
                      OrgDb = org.Mm.eg.db,
                      ont = "All",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      pvalueCutoff = 0.05,
                      readable = T)

dotplot(modEnrich)
cnetplot(modEnrich)

allDat <- merge(cpmT, metaDat, by = "row.names")

ggplot(allDat, aes(x = allDat$ENSMUSG00000021048, y = Rate)) +
  geom_point() + labs(x = "Ncor1 (normalized counts)", y = "Gastroc weight") +
  myTheme + geom_smooth(method = "lm", se = F)

