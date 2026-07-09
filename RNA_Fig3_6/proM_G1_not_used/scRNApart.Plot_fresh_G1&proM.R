---
  #title: "1.run-Seurat_Fig1"
  #date: "2023-02-23"
  #R version: 3.6.3
---
##First part: Seurat QC  
##install.package ----
#install.packages("devtools")
#library(devtools)
#devtools::install_version(package = 'Seurat', version = package_version('3.2.3'))

#install "BiomaRT" package
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("openssl")
#BiocManager::install("biomaRt")

#BiocManager::install("reshape2")
#BiocManager::install("plyr")
#BiocManager::install("mime")
##load.package ----
library(dplyr) #version = package_version('1.0.10')
library(Seurat) #version = package_version('3.2.3')
library(patchwork)  #version = package_version('1.1.2')
library(ggplot2)  #version = package_version('3.4.0')
library(biomaRt) #version = package_version('2.42.1')
library(stringr)  #version = package_version('1.5.0')
##load.dataset and gerenate matrix ----
gene_name <- read.csv("/Users/hxr/Desktop/WJH/fresh_proM_G1//fresh_w_300.gene_names.txt",sep="\t",header = T)
freshCell <- readRDS("/Users/hxr/Desktop/WJH/fresh_proM_G1//fresh_w_300.dgecounts.rds")
barcodes <- read.csv("/Users/hxr/Desktop/WJH/fresh_proM_G1//reads_for_zUMIs.samples_fresh_300.txt",sep="\t")

##choose to use freshCell$readcount or freshCell$umicount to creat seurat objects----
#readscounts
reads <-as.matrix(freshCell$readcount$inex$all)
rownames(reads) <- gene_name$gene_name[match(rownames(freshCell$readcount$inex$all),gene_name$gene_id)] ## change gene id as gene name
#umicounts
reads <-as.matrix(freshCell$umicount$inex$all)
rownames(reads) <- gene_name$gene_name[match(rownames(freshCell$umicount$inex$all),gene_name$gene_id)] ## change gene id as gene name

## transfer as seurat objects
inex_reads <- Seurat::CreateSeuratObject(counts=reads,project="freshCell-inex-reads")

##check the objects ----
head(inex_reads@assays$RNA@counts[1:10,1:2])
rownames(inex_reads@assays$RNA@counts[1:10,1:2])

## calculate percentage contamination ----
inex_reads <-  PercentageFeatureSet(inex_reads,'^MT-',col.name="percent_mito")
inex_reads <-  PercentageFeatureSet(inex_reads,'^RP[SL]',col.name="percent_ribo")
ribo_genes <- rownames(inex_reads)[grep("^RP[SL]",rownames(inex_reads))]
features <- c("percent_ribo")
#"nCount_RNA","nFeature_RNA"
VlnPlot(inex_reads,group.by='orig.ident',features=features,pt.size=0.1,ncol=2)
FeatureScatter(inex_reads,"nCount_RNA","nFeature_RNA",group.by = "orig.ident",pt.size=0.5)
FeatureScatter(inex_reads,"nFeature_RNA","percent_ribo")
FeatureScatter(inex_reads,"nFeature_RNA","percent_mito")

#only select cells with > a features and genes in > b cells ----
selecteda <- WhichCells(inex_reads,expression=nFeature_RNA>800)
selectedb <- rownames(inex_reads)[Matrix::rowSums(inex_reads)>2] #readcounts>3; umicounts>2

inex_reads.filt <- subset(inex_reads,features=selectedb,cells=selecteda)
dim(inex_reads)
dim(inex_reads.filt)

#add original.ID to filtered objects----
rownames(barcodes) <- barcodes[,4]
inex_reads.filt$original.ID <- NA
inex_reads.filt$original.ID[rownames(barcodes)] <- as.character(barcodes[,5])
#add original.batch to filtered objects----
inex_reads.filt$original.batch <- NA
inex_reads.filt$original.batch[rownames(barcodes)] <- as.character(barcodes[,6])

# see the top 10 genes with mapped reads ----
par(mar=c(4,10,2,1))
c <- inex_reads.filt@assays$RNA@counts
c <- Matrix::t(Matrix::t(c)/Matrix::colSums(c)) * 100
most_expressed <- order(apply(c,1,median), decreasing=T)[20:1]
boxplot(as.matrix(t(as.matrix(c[most_expressed,]))),cex=0.1, las=1, xlab="% total count per cell",col=(scales::hue_pal())(20)[20:1],horizontal=TRUE)

# cell cycle gene ----
inex_reads.filt <- NormalizeData(inex_reads.filt)
inex_reads.filt <- CellCycleScoring(object=inex_reads.filt,g2m.features = cc.genes$g2m.genes,s.features=cc.genes$s.genes)
VlnPlot(inex_reads.filt,features=c('S.Score',"G2M.Score"),group.by="orig.ident",ncol=4,pt.size=0.1)

## Nomalize+feature gene+runPCA ----
inex_reads.filt <- NormalizeData(inex_reads.filt,normalization.method ="LogNormalize",scale.factor=10000)
inex_reads.filt <- FindVariableFeatures(inex_reads.filt, selection.method="vst",nfeatures=2000)
top10 <- head(VariableFeatures(inex_reads.filt),10)
p1 <- VariableFeaturePlot(inex_reads.filt)
p2 <- LabelPoints(plot=p1,points=top10,repel=TRUE)
p2

all.genes <- rownames(inex_reads.filt)
inex_reads.filt <- ScaleData(inex_reads.filt, features=all.genes,verbose = FALSE)
inex_reads.filt <- RunPCA(inex_reads.filt,npcs=20)
print(inex_reads.filt[["pca"]],dims=1:5,nfeatures=5)
#VizDimLoadings(inex_reads.filt,dims=1:2,reduction="pca")
DimPlot(inex_reads.filt,reduction="pca")
DimHeatmap(inex_reads.filt,dims=1:5,balanced=TRUE)

##determine the dimensionality----
inex.filt <- JackStraw(inex_reads.filt,num.replicate=100)
ElbowPlot(inex_reads.filt)

##reduce dimention and plot----
#inex_reads.filt <- FindNeighbors(inex_reads.filt, dims = 1:8) #reads_used
inex_reads.filt <- FindNeighbors(inex_reads.filt, dims = 1:5) #umi_used
inex_reads.filt <- FindClusters(inex_reads.filt, resolution = 0.5)
inex_reads.filt <- RunUMAP(inex_reads.filt, reduction = "pca", dims = 1:5,min.dist = 0.3)
p1 <- DimPlot(inex_reads.filt,reduction = "umap",label=TRUE,group.by="ident")
p1
p2 <- DimPlot(inex_reads.filt,reduction = "umap",label=TRUE,group.by="original.ID")
p2
p3 <- DimPlot(inex_reads.filt,reduction = "umap",label=TRUE,group.by="original.batch")
p3

## Second part: get rid of batch effect by using Harmony
##install.package----
#git clone https://github.com/immunogenomics/harmony.git
#In "RunHarmony.R"; remove "reduction.key <- Seurat::Key(reduction.save, quite = TRUE") 
#add "reduction.key = reduction.save"
#library(devtools)
#install_local('~/R/x86_64-pc-linux-gnu-library/3.6/harmony-master/harmony/',force=T)

##load.package----
library(harmony) #version = package_version(‘0.1.1’)  


##Get the separate batch seurat objects ----
fresh_reads.batch1 <- subset(inex_reads.filt, subset = original.batch == "1")
fresh_reads.batch2 <- subset(inex_reads.filt, subset = original.batch == "2")

# Merge raw samples
fresh_merge <- merge(x = fresh_reads.batch1,
                     y = fresh_reads.batch2,
                     merge.data = TRUE)

# Perform log-normalization and feature selection, as well as SCT normalization on global object
fresh_merge <- fresh_merge %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>%
  SCTransform(vars.to.regress = c("percent_mito"))

# Calculate PCs using variable features determined by SCTransform (3000 by default)
fresh_merge <- RunPCA(fresh_merge, assay = "SCT", npcs = 50) 

fresh_merge_harmonized <- RunHarmony(fresh_merge, 
                                group.by.vars = c("original.batch"),
                                reduction = "pca", assay.use = "SCT", reduction.save = "harmony",
                                theta = 10) #umicounts used 7/readscounts used 2
ElbowPlot(fresh_merge_harmonized)
#fresh_merge_harmonized <- RunUMAP(fresh_merge_harmonized, reduction = "harmony", assay = "SCT", dims = 1:8) #reads_acount
fresh_merge_harmonized <- RunUMAP(fresh_merge_harmonized, reduction = "harmony", assay = "SCT", dims = 1:4) #umi_account

fresh_merge_harmonized <- FindNeighbors(object = fresh_merge_harmonized, reduction = "harmony")
fresh_merge_harmonized <- FindClusters(fresh_merge_harmonized, resolution =  0.3)
DimPlot(fresh_merge_harmonized, group.by = c("original.batch", "ident", "original.ID"), ncol = 3)

#save dot position in Umap
fresh_merge_harmonized_umap = fresh_merge_harmonized@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(clusters = fresh_merge_harmonized@meta.data$seurat_clusters) %>% cbind(originalbatch = fresh_merge_harmonized$original.batch) %>% cbind(originalbatch = fresh_merge_harmonized$original.ID)
fresh_merge_harmonized = fresh_merge_harmonized@reductions$harmony@cell.embeddings %>% as.data.frame() %>% cbind(clusters = fresh_merge_harmonized@meta.data$seurat_clusters) %>% cbind(originalbatch = fresh_merge_harmonized$original.batch) %>% cbind(originalbatch = fresh_merge_harmonized$original.ID)
write.table(fresh_merge_harmonized_umap, file = "/Users/hxr/Desktop/WJH/fresh_proM_G1/fresh_G1&proM_umap_3.csv",quote=F) 
write.table(fresh_merge_harmonized, file = "/Users/hxr/Desktop/WJH/fresh_proM_G1/fresh_G1&proM_harmony_3.csv",quote=F) 



##plot the feature gene----
S_gene <-c("ORC1","CCNE1","CCNE2","MCM6")
G2M_gene <-c("WEE1","CDK1","CCN1","NUSAP1","AURKA","CCNA2","CCNB2")

#FeaturePlot(fresh_merge_harmonized,features = G2M_gene)
VlnPlot(fresh_merge_harmonized,features=G2M_gene)
VlnPlot(fresh_merge_harmonized,features=S_gene)
#check the absolute numbers
VlnPlot(fresh_merge_harmonized, features = G2M_gene, slot = "counts", log = TRUE) 
VlnPlot(fresh_merge_harmonized, features = S_gene, slot = "counts", log = TRUE)
#check the mito/ribo and nCount/nFeature
VlnPlot(fresh_merge_harmonized, features = c("percent_mito",'percent_ribo'),slot = "counts", log = TRUE)
VlnPlot(fresh_merge_harmonized, features = c("nCount_RNA", "nFeature_RNA"))

##try use densUMAP to plot UMAP----
#Sys.setenv(RETICULATE_PYTHON = "/Users/hxr/miniconda3/envs/r-reticulate-python38/bin/python")
library(reticulate)
#py_config()
Fresh_merge_dens <-Seurat::RunUMAP(fresh_merge_harmonized, dims = 1:4, umap.method = "umap-learn", densmap=TRUE)
Fresh_merge_dens <- FindNeighbors(object = Fresh_merge_dens, reduction = "harmony")
Fresh_merge_dens <- FindClusters(Fresh_merge_dens, resolution =  0.3)
DimPlot(Fresh_merge_dens, group.by = c("original.batch", "ident", "original.ID"), ncol = 3)

