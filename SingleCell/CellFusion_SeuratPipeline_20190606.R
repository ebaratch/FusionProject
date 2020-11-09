## Cell Fusion Analysis Pipeline

## run on cluster using: qsub -I -l mem=128gb
## module load R/3.5.2
## R

library(Matrix)
library(dplyr)
library(readr)
library(rdetools)
library(data.table)
library(ggplot2)
library(iterators)
library(Seurat)
sessionInfo()

# Set working directory 
#     RDS files were saved at each step along the way, but aren't included here
#     Diversity outputs at the end of the pipeline are included in the diversity-results folder

####################################################################
##### Step 1: Load all the data
####################################################################
MDA.data <- Read10X(data.dir = "./MDA231/outs/filtered_gene_bc_matrices/MegCustomGenome/")
MDA <- CreateSeuratObject(counts = MDA.data, project = "MDA231")
SUM.data <- Read10X(data.dir = "./SUM159/outs/filtered_gene_bc_matrices/MegCustomGenome/")
SUM <- CreateSeuratObject(counts = SUM.data, project = "SUM159")
P2.data <- Read10X(data.dir = "./P2/outs/filtered_gene_bc_matrices/MegCustomGenome/")
P2 <- CreateSeuratObject(counts = P2.data, project = "P2")
P10.data <- Read10X(data.dir = "./P10/outs/filtered_gene_bc_matrices/MegCustomGenome/")
P10 <- CreateSeuratObject(counts = P10.data, project = "P10")

####################################################################
#### Step 1b: Quality Control/Subset the Data Individually 
####          Done on the per sample  basis
####################################################################

# ============== MDA-MB-231 =====================
mito.genes <- grep(pattern = "^MT-", x = rownames(x=MDA), value = TRUE)
percent.mitoMDA <- Matrix::colSums(x=GetAssayData(object=MDA, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=MDA, slot='counts'))
MDA[['percent.mito']] <- percent.mitoMDA

## Descriptive stats:
summary(MDA@meta.data$nFeature_RNA) 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#991    2921    3450    3560    4151    6746 
mean(MDA@meta.data$nFeature_RNA, na.rm=TRUE) # 3559.824
sd(MDA@meta.data$nFeature_RNA, na.rm=TRUE) # 943.8025
## Determine bounds
nFeatUpperMDA <- mean(MDA@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(MDA@meta.data$nFeature_RNA, na.rm=TRUE) # 5447.429
nFeatLowerMDA <- mean(MDA@meta.data$nFeature_RNA, na.rm=TRUE) - 2*sd(MDA@meta.data$nFeature_RNA, na.rm=TRUE) # 1672.219
perMitoUpperMDA <- mean(MDA@meta.data$percent.mito, na.rm=TRUE) + 2*sd(MDA@meta.data$percent.mito, na.rm=TRUE) # 0.1122103

pdf("./MDAonly_VlnPlot_nFeature-nCount-perMito_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=MDA, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
print(vPlot)
dev.off()
pdf("./MDAonly_VlnPlot_nFeature-with-cutoffs_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=MDA, features = c("nFeature_RNA"))+geom_hline(yintercept=nFeatUpperMDA, linetype="dashed", color="red", size=1)+geom_hline(yintercept=nFeatLowerMDA, linetype="dashed", color="red", size=1)
print(vPlot)
dev.off()
pdf("./MDAonly_VlnPlot_perMito-with-cutoff_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=MDA, features = c("percent.mito"))+geom_hline(yintercept=perMitoUpperMDA, linetype="dashed", color="red", size=1)
print(vPlot)
dev.off()
pdf("./MDAonly_ScatterPlot_RNAcount-vs-mito_MegCustom_2019-06-06.pdf ")
fsPlot <- FeatureScatter(object=MDA, feature1 = "nCount_RNA", feature2 = "percent.mito")
print(fsPlot)
dev.off()
pdf("./MDAonly_ScatterPlot_nGenes-vs-nRNA_MegCustom_2019-06-06.pdf")
fsPlot <- FeatureScatter(object=MDA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(fsPlot)
dev.off()

# We filter out cells that have unique feature counts over 2*SD+Mean or less than Mean-2*SD and mitochondrial content over 2*SD+Mean
MDA <- subset(x=MDA, subset = nFeature_RNA > nFeatLowerMDA & nFeature_RNA < nFeatUpperMDA & percent.mito < perMitoUpperMDA)
MDA <- RenameCells(object=MDA, add.cell.id="MDA231")

# ================== SUM159 ==========================
percent.mitoSUM <- Matrix::colSums(x=GetAssayData(object=SUM, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=SUM, slot='counts'))
SUM[['percent.mito']] <- percent.mitoSUM

## Descriptive stats:
summary(SUM@meta.data$nFeature_RNA) 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1093    2654    3062    3103    3501    6341 
mean(SUM@meta.data$nFeature_RNA, na.rm=TRUE) # 3102.927
sd(SUM@meta.data$nFeature_RNA, na.rm=TRUE) # 725.0492
## Determine bounds
nFeatUpperSUM <- mean(SUM@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(SUM@meta.data$nFeature_RNA, na.rm=TRUE) # 4553.026
nFeatLowerSUM <- mean(SUM@meta.data$nFeature_RNA, na.rm=TRUE) - 2*sd(SUM@meta.data$nFeature_RNA, na.rm=TRUE) # 1652.829
perMitoUpperSUM <- mean(SUM@meta.data$percent.mito, na.rm=TRUE) + 2*sd(SUM@meta.data$percent.mito, na.rm=TRUE) # 0.05529778

pdf("./SUMonly_VlnPlot_nFeature-nCount-perMito_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=SUM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
print(vPlot)
dev.off()
pdf("./SUMonly_VlnPlot_nFeature-with-cutoffs_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=SUM, features = c("nFeature_RNA"))+
  geom_hline(yintercept=nFeatUpperSUM, linetype="dashed", color="red", size=1)+
  geom_hline(yintercept=nFeatLowerSUM, linetype="dashed", color="red", size=1)
print(vPlot)
dev.off()
pdf("./SUMonly_VlnPlot_perMito-with-cutoff_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=SUM, features=c("percent.mito"))+
  geom_hline(yintercept=perMitoUpperSUM, linetype="dashed", color="red", size=1)
print(vPlot)
dev.off()
pdf("./SUMonly_ScatterPlot_RNAcount-vs-mito_MegCustom_2019-06-06.pdf ")
fsPlot <- FeatureScatter(object=SUM, feature1="nCount_RNA", feature2="percent.mito")
print(fsPlot)
dev.off()
pdf("./SUMonly_ScatterPlot_nGenes-vs-nRNA_MegCustom_2019-06-06.pdf")
fsPlot <- FeatureScatter(object=SUM, feature1="nCount_RNA", feature2="nFeature_RNA")
print(fsPlot)
dev.off()

# We filter out cells that have unique feature counts over 2*SD+Mean or less than Mean-2*SD and mitochondrial content over 2*SD+Mean
SUM <- subset(x=SUM, subset = nFeature_RNA > nFeatLowerSUM & nFeature_RNA < nFeatUpperSUM & percent.mito < perMitoUpperSUM)
SUM <- RenameCells(object=SUM, add.cell.id="SUM159")

# ================= Passage 2 =======================
percent.mitoP2 <- Matrix::colSums(x=GetAssayData(object=P2, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=P2, slot='counts'))
P2[['percent.mito']] <- percent.mitoP2

## Descriptive stats:
summary(P2@meta.data$nFeature_RNA) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 884    2463    2861    2987    3417    7552
mean(P2@meta.data$nFeature_RNA, na.rm=TRUE) # 2987.317
sd(P2@meta.data$nFeature_RNA, na.rm=TRUE) # 783.0067
## Determine bounds
nFeatUpperP2 <- mean(P2@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(P2@meta.data$nFeature_RNA, na.rm=TRUE) # 4553.33
nFeatLowerP2 <- mean(P2@meta.data$nFeature_RNA, na.rm=TRUE) - 2*sd(P2@meta.data$nFeature_RNA, na.rm=TRUE) # 1421.304
perMitoUpperP2 <- mean(P2@meta.data$percent.mito, na.rm=TRUE) + 2*sd(P2@meta.data$percent.mito, na.rm=TRUE) # 0.06463388

pdf("./P2only_VlnPlot_nFeature-nCount-perMito_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=P2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
print(vPlot)
dev.off()
pdf("./P2only_VlnPlot_nFeature-with-cutoffs_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=P2, features = c("nFeature_RNA"))+
  geom_hline(yintercept=nFeatUpperP2, linetype="dashed", color="red", size=1)+
  geom_hline(yintercept=nFeatLowerP2, linetype="dashed", color="red", size=1)
print(vPlot)
dev.off()
pdf("./P2only_VlnPlot_perMito-with-cutoff_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=P2, features=c("percent.mito"))+
  geom_hline(yintercept=perMitoUpperP2, linetype="dashed", color="red", size=1)
print(vPlot)
dev.off()
pdf("./P2only_ScatterPlot_RNAcount-vs-mito_MegCustom_2019-06-06.pdf ")
fsPlot <- FeatureScatter(object=P2, feature1="nCount_RNA", feature2="percent.mito")
print(fsPlot)
dev.off()
pdf("./P2only_ScatterPlot_nGenes-vs-nRNA_MegCustom_2019-06-06.pdf")
fsPlot <- FeatureScatter(object=P2, feature1="nCount_RNA", feature2="nFeature_RNA")
print(fsPlot)
dev.off()

# We filter out cells that have unique feature counts over 2*SD+Mean or less than Mean-2*SD and mitochondrial content over 2*SD+Mean
P2 <- subset(x=P2, subset = nFeature_RNA > nFeatLowerP2 & nFeature_RNA < nFeatUpperP2 & percent.mito < perMitoUpperP2)
P2 <- RenameCells(object=P2, add.cell.id="P2")

# ===================== Passage 10 ==========================
percent.mitoP10 <- Matrix::colSums(x=GetAssayData(object=P10, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=P10, slot='counts'))
P10[['percent.mito']] <- percent.mitoP10

## Descriptive stats:
summary(P10@meta.data$nFeature_RNA) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1239    2556    3132    3196    3739    8111 
mean(P10@meta.data$nFeature_RNA, na.rm=TRUE) # 3196.279
sd(P10@meta.data$nFeature_RNA, na.rm=TRUE) # 901.6977
## Determine bounds
nFeatUpperP10 <- mean(P10@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(P10@meta.data$nFeature_RNA, na.rm=TRUE) # 4999.674
nFeatLowerP10 <- mean(P10@meta.data$nFeature_RNA, na.rm=TRUE) - 2*sd(P10@meta.data$nFeature_RNA, na.rm=TRUE) # 1392.883
perMitoUpperP10 <- mean(P10@meta.data$percent.mito, na.rm=TRUE) + 2*sd(P10@meta.data$percent.mito, na.rm=TRUE) # 0.08535412

pdf("./P10only_VlnPlot_nFeature-nCount-perMito_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=P10, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
print(vPlot)
dev.off()
pdf("./P10only_VlnPlot_nFeature-with-cutoffs_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=P10, features = c("nFeature_RNA"))+
  geom_hline(yintercept=nFeatUpperP10, linetype="dashed", color="red", size=1)+
  geom_hline(yintercept=nFeatLowerP10, linetype="dashed", color="red", size=1)
print(vPlot)
dev.off()
pdf("./P10only_VlnPlot_perMito-with-cutoff_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=P10, features=c("percent.mito"))+
  geom_hline(yintercept=perMitoUpperP10, linetype="dashed", color="red", size=1)
print(vPlot)
dev.off()
pdf("./P10only_ScatterPlot_RNAcount-vs-mito_MegCustom_2019-06-06.pdf ")
fsPlot <- FeatureScatter(object=P10, feature1="nCount_RNA", feature2="percent.mito")
print(fsPlot)
dev.off()
pdf("./P10only_ScatterPlot_nGenes-vs-nRNA_MegCustom_2019-06-06.pdf")
fsPlot <- FeatureScatter(object=P10, feature1="nCount_RNA", feature2="nFeature_RNA")
print(fsPlot)
dev.off()

# We filter out cells that have unique feature counts over 2*SD+Mean or less than Mean-2*SD and mitochondrial content over 2*SD+Mean
P10 <- subset(x=P10, subset = nFeature_RNA > nFeatLowerP10 & nFeature_RNA < nFeatUpperP10 & percent.mito < perMitoUpperP10)
P10 <- RenameCells(object=P10, add.cell.id="P10")

####################################################################
#### Step 2: Combine Data
####################################################################

CellFusion1.combined <- merge(x=MDA, y=SUM, add.cell.id1 = "M", add.cell.id2 = "S", merge.data=TRUE, project = "MDA+SUM")
CellFusion2.combined <- merge(x = CellFusion1.combined, y = P2,  add.cell.id2 = "P2", merge.data=TRUE, project = "MDA+SUM")
CellFusion.combined <- merge(x = CellFusion2.combined, y = P10,  add.cell.id2 = "P10", merge.data=TRUE, project = "MDA+SUM")

pdf("./CombinedCellFusion_VlnPlot_nFeature-nCount-perMito_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=CellFusion.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
print(vPlot)
dev.off()

pdf("./CombinedCellFusion_VlnPlot_nFeature-nCount-perMito_distribOnly_MegCustom_2019-06-06.pdf")
vPlot <- VlnPlot(object=CellFusion.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol=3, pt.size=0)
print(vPlot)
dev.off()

saveRDS(CellFusion.combined, file = "./CombinedMDA231+SUM159_individSubset_MegCustom_SeuratObj_20190606.rds")
# CellFusion.combined <- readRDS("./CombinedMDA231+SUM159_individSubset_SeuratObj_20190405.rds")

####################################################################
#### Step 3: Normalize
####################################################################

CellFusion.combined <- NormalizeData(object = CellFusion.combined , normalization.method = "LogNormalize", scale.factor = 1e4)

####################################################################
##### Step 4: Detection of variable features across the single cells
####################################################################

CellFusion.combined <- FindVariableFeatures(object=CellFusion.combined, mean.cutoff=c(0.0125, 3), dispersion.cutoff=c(0.5, Inf))
length(x=VariableFeatures(object=CellFusion.combined))  
## number of variable features detected after filtering, normalization: 2000 (default max)

## genes with high cell-to-cell variation
top10 <- head(x=VariableFeatures(object=CellFusion.combined), 10)
top25 <- head(x=VariableFeatures(object=CellFusion.combined), 25)
top100 <- head(x=VariableFeatures(object=CellFusion.combined), 100)
pdf("./CellFusion_DispersionPlot_VariableGenes_Top25_MegCustom_2019-06-06.pdf", height=8, width=11)
plot1 <- VariableFeaturePlot(object=CellFusion.combined)
plot2 <- LabelPoints(plot=plot1, points=top25, repel=TRUE)
cPlot <- CombinePlots(plots = list(plot1, plot2))
print(cPlot)
dev.off() 

pdf("./CellFusion_DispersionPlot_VariableGenes_Top25_LabeledOnly_MegCustom_2019-06-06.pdf", height=8, width=10)
print(plot2)
dev.off()

write.csv(top100, "./CellFusion_Top100_VariableGenes_2019-06-06.csv")

####################################################################
##### Step 5: Scaling the data and removing unwanted sources of variation 
####################################################################

CellFusion.combined <- ScaleData(object=CellFusion.combined, features=rownames(x=CellFusion.combined), vars.to.regress = "nCount_RNA") #c("nCount_RNA", "percent.mito"))
saveRDS(CellFusion.combined, file = "./CombinedMDA231+SUM159_individSubset_SeuratObj_noNA_MegCustom_thrStep5_20190606.rds")
# CellFusion.combined <- readRDS("./CombinedMDA231+SUM159_individSubset_SeuratObj_noNA_thrStep5_20190415.rds")

####################################################################
##### Step 6: Perform linear dimensional reduction 
####################################################################

CellFusion.combined <- RunPCA(object = CellFusion.combined , features = VariableFeatures(object = CellFusion.combined), verbose = TRUE, npcs=200)
saveRDS(CellFusion.combined, file = "./CombinedMDA231+SUM159_individSubset_SeuratObj_noNA_MegCustom_thrStep6_20190606.rds")

####################################################################
##### Step 7: Cluster the Cells
####################################################################

CellFusion.combined <- FindNeighbors(object = CellFusion.combined, dims=1:50, reduction="pca")
CellFusion.combined <- FindClusters(object = CellFusion.combined, resolution=0.6)

# CLUSTERING RESULTS
# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
# Number of nodes: 10059
# Number of edges: 413785
#
# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
#  [----|----|----|----|----|----|----|----|----|----|
#     **************************************************|
#   Maximum modularity in 10 random starts: 0.8931
#  Number of communities: 12
#  Elapsed time: 13 seconds

saveRDS(CellFusion.combined, file = "./CombinedMDA231+SUM159_individSubset_SeuratObj_noNA_MegCustom_thrStep7_20190606.rds")

####################################################################
##### Step 8: Run Non-linear dimensional reduction: UMAP 
####################################################################
### picked up here on laptop, ran into issues on the cluster with UMAP

CellFusion.combined <- readRDS("./CombinedMDA231+SUM159_individSubset_SeuratObj_noNA_MegCustom_thrStep7_20190606.rds")

CellFusion.combined <- RunUMAP(object = CellFusion.combined, reduction = "pca")
pdf("./CombinedCellFusion_UMAP_res0.6_12clusters_byCluster_MegCustom_2019-06-07.pdf")
uPlot <- DimPlot(CellFusion.combined, reduction.use = "umap")
print(uPlot)
dev.off()
pdf("./CombinedCellFusion_UMAP_res0.6_12clusters_byType_MegCustom_2019-06-07.pdf")
uPlot <- DimPlot(CellFusion.combined, reduction.use = "umap", group.by='orig.ident')
print(uPlot)
dev.off()

pdf("./CombinedCellFusion_UMAP_res0.6_FeaturePlot-of-Custom-Genes_2019-06-07.pdf")
features <- c("EGFP", "mCherry")
fPlot <- FeaturePlot(object = CellFusion.combined, features = features)
print(fPlot)
dev.off()

pdf("./CombinedCellFusion_RidgePlot-of-Custom-Genes_2019-06-07.pdf")
rPlot <- RidgePlot(object = CellFusion.combined, features = features, ncol=2)
print(rPlot)
dev.off()

pdf("./CombinedCellFusion_ViolinPlot-of-Custom-Genes_2019-06-07.pdf")
vPlot <- VlnPlot(object = CellFusion.combined, features = features)
print(vPlot)
dev.off()

pdf("./CombinedCellFusion_DotPlot-of-Custom-Genes_2019-06-07.pdf")
dPlot <- DotPlot(object = CellFusion.combined, features = features, group.by='orig.ident') + RotatedAxis()
print(dPlot)
dev.off()

pdf("./CombinedCellFusion_Heatmap-of-Custom-Genes_2019-06-07.pdf")
dhmPlot <- DoHeatmap(object = subset(x = CellFusion.combined, downsample = 250), features = features, size = 3)
print(dhmPlot)
dev.off()

saveRDS(CellFusion.combined, file = "./CombinedMDA231+SUM159_individSubset_SeuratObj_noNA_MegCustom_thrStep8_20190607.rds")

####################################################################
##### Step 9: Diversity Score
####################################################################

file1 <- "./outputDataNames_noNA_MegCustom_20190607.csv";
file2 <- "./outputData_noNA_MegCustom_20190607.csv"

cat(CellFusion.combined@meta.data$orig.ident, file=file1, sep=",\n")
cat(CellFusion.combined@meta.data$RNA_snn_res.0.6, file=file2, sep=",\n")

mydat1 <- read.csv(file2)
mydat2 <- read.csv(file1)
fulldat <- cbind(mydat2[1],mydat1[1])
fulltab <- as.data.table(fulldat)
# name table columns
names(fulltab)[1] <- paste("UMI")
names(fulltab)[2] <- paste("cluster")
# group data based on clusters
group_by(fulltab, cluster)
# create a table counting unqiue UMIs/cells per cluster
tabPerClus <- fulltab %>% group_by(cluster) %>% count()
type <- sub("\\_.*","",fulltab$UMI)
fulltab <- cbind(fulltab, type)

tmp <- matrix(ncol=1, nrow=dim(tabPerClus)[1])
for(i in 1:dim(tabPerClus)[1]){
  tmp[i] <- tabPerClus$n[i]/dim(fulltab)[1]
}
tmp <- as.data.frame(tmp)
names(tmp)[1] <- paste("freq")
tabPerClusFreq <- cbind(as.data.frame(tabPerClus), tmp)

calcqD <- function(dat, q){
  diversity <- 0.0;
  if(q == 1.0){
    for(row in 1:dim(dat)[1]){
      diversity <- diversity + log((dat$freq[row])^(dat$freq[row]))
    }
    diversity <- abs(diversity)
  }else{
    for(row in 1:dim(dat)[1]){
      diversity <- diversity + (dat$freq[row])^q
    }
    diversity <- diversity^(1/(1-q))
  }
}

qDp01 <- calcqD(tabPerClusFreq,0.01)
qDp1 <- calcqD(tabPerClusFreq,0.1)
qD1 <- calcqD(tabPerClusFreq,1.000001)
qD2 <- calcqD(tabPerClusFreq,2.0)
qD10 <- calcqD(tabPerClusFreq,10.0)
qD100 <- calcqD(tabPerClusFreq,100.0)

clus <- {}
cbind(clus,qDp01,qDp1,qD1,qD2,qD10,qD100)

# creates a table containing the number of cells per cluster per condition
tabPerClusType <- na.omit(fulltab) %>% group_by(cluster, type) %>% count()
# create a table containing the number of unique per condition (lung1norm, etc.)
tabPerType <- na.omit(fulltab) %>% group_by(type) %>% count()

# calculate the frequency of each cells out of total number of cells per condition in each cluster
# these are depended on total number of cells in the specific condition (tabPerType)
# the fraction of cells in cluster X out of all cells in condition Y
tmp <- matrix(ncol=1, nrow=dim(tabPerClusType)[1])
for(i in 1:dim(tabPerClusType)[1]){
  if(tabPerClusType$type[i]==tabPerType$type[1]) {
    ## MDA231
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[1]
  } else if (tabPerClusType$type[i]==tabPerType$type[2]) {
    ## P2
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[2]
  } else if (tabPerClusType$type[i]==tabPerType$type[3]) {
    ## P10
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[3]
  } else if (tabPerClusType$type[i]==tabPerType$type[4]) {
    ## SUM159
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[4]
  } else {   
    tmp[i] <- 0.0
  }
}
tmp <- as.data.frame(tmp)
names(tmp)[1] <- paste("freq")
tabPerClusType <- cbind(as.data.frame(tabPerClusType), tmp)

# filtering the overall table to have a smaller table with just 1 condition in it
mdaClus <- filter(tabPerClusType,type=="MDA231")
sumClus <- filter(tabPerClusType,type=="SUM159")
p2Clus <- filter(tabPerClusType,type=="P2")
p10Clus <- filter(tabPerClusType,type=="P10")

# function to calculate diversity score for a given q
calcqD <- function(dataset, q){
  diversity <- 0.0;
  for(row in 1:dim(dataset)[1]){
    diversity <- diversity + (dataset$freq[row])^q
  }
  diversity <- diversity^(1/(1-q))
}

# calculates diversity scores based on cluster frequencies for a range of q
qRange <- logspace(-2, 2, n = 1000)
qDmda <- calcqD(mdaClus,qRange)
qDmdatab <- as.data.frame(cbind(qRange,qDmda))
qDsum <- calcqD(sumClus,qRange)
qDsumtab <- as.data.frame(cbind(qRange,qDsum))
qDp2 <- calcqD(p2Clus,qRange)
qDp2tab <- as.data.frame(cbind(qRange,qDp2))
qDp10 <- calcqD(p10Clus,qRange)
qDp10tab <- as.data.frame(cbind(qRange,qDp10))

# saved the qD tables
write.csv(qDmdatab, file="Parental-MDAMD231-qD-table_Custom_2019-06-07.csv")
write.csv(qDsumtab, file="Parental-SUM159-qD-table_Custom_2019-06-07.csv")
write.csv(qDp2tab, file="Passage2-qD-table_Custom_2019-06-07.csv")
write.csv(qDp10tab, file="Passage10-qD-table_Custom_2019-06-07.csv")

# plots all data together
pdf("./CellFusion_DiversityCurve_12clus_Custom_2019-06-07.pdf")
myPlot <- ggplot(data=qDmdatab, aes(x=qRange, y=qDmda), color="blue3") +
  geom_smooth(size=4) + 
  scale_x_log10() +
  geom_smooth(data=qDsumtab, aes(x=qRange, y=qDsum),size=4, color="red3") +
  geom_smooth(data=qDp2tab, aes(x=qRange, y=qDp2),size=4, color="orchid3") +
  geom_smooth(data=qDp10tab, aes(x=qRange, y=qDp10),size=4, color="purple4", ylim=c(0,12)) +
  labs(title="Diversity Score",x="q", y = "qD") +
  theme(legend.position="right") + 
  theme_light()
print(myPlot)
dev.off()

## Legend 
## MDA231 = Blue; SUM159= Red; P2 = Light Purple; P10 = Dark purple

scale_color_manual(labels=c("MDA231","SUM159","Early Fusion Passage","Late Fusion Passage"), 
                   breaks = c("MDA231","SUM159","Early Fusion Passage","Late Fusion Passage"),
                   values = c("blue3","red3","orchid3","purple4"),
                   guide = guide_legend(override.aes=aes(fill=NA)))

## Plotting individual Q's
lowQ <- c(qDmdatab[1,2], qDsumtab[1,2], qDp2tab[1,2], qDp10tab[1,2])
pdf("./CellFusion_DiversityCurve_12clus_lowQ_Custom_2019-06-07.pdf")
bPlot <- barplot(lowQ, ylab="diversity", ylim=c(0,12))
print(bPlot)
dev.off()

highQ <- c(qDmdatab[1000,2], qDsumtab[1000,2], qDp2tab[1000,2], qDp10tab[1000,2])
pdf("./CellFusion_DiversityCurve_12clus_highQ_Custom_2019-06-07.pdf")
bPlot <- barplot(highQ, ylab="diversity", ylim=c(0,12))
print(bPlot)
dev.off()

qOf1 <- c(mean(qDmdatab[500:501,2]), mean(qDsumtab[500:501,2]),mean(qDp2tab[500:501,2]),mean(qDp10tab[500:501,2])) 
pdf("./CellFusion_DiversityCurve_12clus_qOf1_Custom_2019-06-07.pdf")
bPlot <- barplot(qOf1, ylab="diversity", ylim=c(0,12))
print(bPlot)
dev.off()

####################################################################
##### CHECKING COMPOSITION OF  CELLS 
####################################################################

fulltab %>% group_by(type) %>% count()
# A tibble: 4 x 2
# Groups:   type [4]
#  type       n
#  <chr>  <int>
#  1 MDA231  2718
#  2 P10     2023
#  3 P2      2588
#  4 SUM159  2729

fulltab %>% group_by(cluster) %>% count()
# A tibble: 12 x 2
# Groups:   cluster [12]
#  cluster     n
#  <int> <int>
#  1       1  2357
#  2       2  2074
#  3       3  1109
#  4       4   845
#  5       5   833
#  6       6   594
#  7       7   414
#  8       8   409
#  9       9   390
#  10      10   389
#  11      11   364
#  12      12   280

## Per cluster per type breakdown: 
A <- fulltab %>% group_by(cluster) %>% count(type)
write.csv(A, file="CellBreakdown_PerClusterPerType_Custom_2019-06-07.csv")


####################################################################
##### OPTIONAL: JACK STRAWING ANALYSIS
####################################################################

# NOTE: This process can take a long time for big datasets, comment out for expediency.
# More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time

CellFusion.combined <- readRDS("./CombinedMDA231+SUM159_individSubset_SeuratObj_noNA_MegCustom_thrStep8_20190607.rds")
CellFusion.combined <- JackStraw(object=CellFusion.combined, num.replicate=200, dims=50)
CellFusion.combined <- ScoreJackStraw(object=CellFusion.combined, dims=1:50)

pdf("./CellFusion_JackStrawPlot_200x_50dim_MegCustom_2019-07-03.pdf")
jPlot <- JackStrawPlot(object = CellFusion.combined, dims = 1:50)
print(jPlot)
dev.off()
pdf("./CellFusion_JackStrawPlot_200x_50dim_MegCustom_2019-07-03_NO-LEGEND.pdf")
jPlot <- JackStrawPlot(object = CellFusion.combined, dims = 1:50)+theme(legend.position="none")
print(jPlot)
dev.off()
pdf("./CellFusion_JackStraw_ElbowPlot_200x_50dim_MegCustom_2019-07-03.pdf")
ePlot <- ElbowPlot(object = CellFusion.combined, ndims=50)
print(ePlot)
dev.off()

saveRDS(CellFusion.combined, file="./CombinedMDA231+SUM159_individSubset_SeuratObj_noNA_MegCustom_withJackStraw_thrStep8_20190703.rds")
