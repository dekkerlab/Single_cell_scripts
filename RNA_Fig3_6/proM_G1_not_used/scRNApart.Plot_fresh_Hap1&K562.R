---
  #title: "1.run-Seurat_Fig1"
  #date: "2023-02-23"
  #R version: 3.6.3
---
  
  ##load.package ----
library(dplyr) #version = package_version('1.0.10')
library(Seurat) #version = package_version('3.2.3')
library(patchwork)  #version = package_version('1.1.2')
library(ggplot2)  #version = package_version('3.4.0')
library(biomaRt) #version = package_version('2.42.1')
library(stringr)  #version = package_version('1.5.0')
##load.dataset and gerenate matrix ----
Fresh_2cells <- readRDS("/Users/hxr/Desktop/WJH/fresh_proM_G1//Fresh_Hap1_K562.dgecounts.rds")
barcodes <- read.csv("/Users/hxr/Desktop/WJH/fresh_proM_G1//Fresh_reads_for_zUMIs.samples.txt",sep="\t")
gene_name <- read.csv("/Users/hxr/Desktop/WJH/fresh_proM_G1//Fresh_Hap1_K562.gene_names.txt",sep="\t",header = T)

##choose to use Fresh_2cells$readcount or Fresh_2cells$umicount to creat seurat objects (I used umicount)----
reads <-as.matrix(Fresh_2cells$umicount$inex$all)
rownames(reads) <- gene_name$gene_name[match(rownames(Fresh_2cells$umicount$inex$all),gene_name$gene_id)] ## change gene id as gene name

## transfer as seurat objects
inex_reads <- Seurat::CreateSeuratObject(counts=reads,project="Fresh_2cells-reads")

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
selecteda <- WhichCells(inex_reads,expression=(nFeature_RNA < 10000)) #readcounts>980
#selecteda <- WhichCells(inex_reads,expression=nFeature_RNA>700) #umicounts>700
selectedb <- rownames(inex_reads)[Matrix::rowSums(inex_reads)>2] #readcounts_exon>2; umicounts_inex>1

inex_reads.filt <- subset(inex_reads,features=selectedb,cells=selecteda)
dim(inex_reads)
dim(inex_reads.filt)

#add original.ID to filtered objects----
rownames(barcodes) <- barcodes[,4]
inex_reads.filt$original.ID <- NA
inex_reads.filt$original.ID[rownames(barcodes)] <- as.character(barcodes[,6])
#add original.batch to filtered objects----
inex_reads.filt$original.batch <- NA
inex_reads.filt$original.batch[rownames(barcodes)] <- as.character(barcodes[,5])

# see the top 10 genes with mapped reads ----
par(mar=c(4,10,2,1))
c <- inex_reads.filt@assays$RNA@counts
c <- Matrix::t(Matrix::t(c)/Matrix::colSums(c)) * 100
most_expressed <- order(apply(c,1,median), decreasing=T)[20:1]
boxplot(as.matrix(t(as.matrix(c[most_expressed,]))),cex=0.1, las=1, xlab="% total count per cell",col=(scales::hue_pal())(20)[20:1],horizontal=TRUE)

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
DimPlot(inex_reads.filt,reduction="pca",group.by="original.ID")
DimPlot(inex_reads.filt,reduction="pca",group.by="ident")
DimPlot(inex_reads.filt,reduction="pca",group.by="original.batch")
DimHeatmap(inex_reads.filt,dims=1:5,balanced=TRUE)

##determine the dimensionality----
inex.filt <- JackStraw(inex_reads.filt,num.replicate=100)
ElbowPlot(inex_reads.filt)

##reduce dimention and plot----
#inex_reads.filt <- FindNeighbors(inex_reads.filt, dims = 1:20) #reads_used
inex_reads.filt <- FindNeighbors(inex_reads.filt, dims = 1:7) #umi_used
inex_reads.filt <- FindClusters(inex_reads.filt, resolution = 0.5)
inex_reads.filt <- RunUMAP(inex_reads.filt, reduction = "pca", dims = 1:10,min.dist = 0.3)
p1 <- DimPlot(inex_reads.filt,reduction = "umap",label=TRUE,group.by="ident")
p1
p2 <- DimPlot(inex_reads.filt,reduction = "umap",label=TRUE,group.by="original.ID")
p2
p3 <- DimPlot(inex_reads.filt,reduction = "umap",label=TRUE,group.by="original.batch")
p3

## Second part: get rid of batch effect by using Harmony
##install.package----
#git clone https://github.com/immunogenomics/harmony.git
#In "RunHarmony.R"; remove "reduction.key <- Seurat::Key(reduction.save, quite = TRUE" 
#add "reduction.key = reduction.save"
#library(devtools)
#install_local('~/R/x86_64-pc-linux-gnu-library/3.6/harmony-master/harmony/',force=T)




##load.package----
library(harmony) #version = package_version(‘0.1.1’)


##Get the separate batch seurat objects ----
fresh_2cells_umi.batch1 <- subset(inex_reads.filt, subset = original.batch == "1")
fresh_2cells_umi.batch2 <- subset(inex_reads.filt, subset = original.batch == "2")

# Merge raw samples
fresh_2cells_list <- list(fresh_2cells_umi.batch1,fresh_2cells_umi.batch2)
fresh_2cells_merge <- Reduce(function(x, y) merge(x, y, merge.data = TRUE), fresh_2cells_list)

# Perform log-normalization and feature selection, as well as SCT normalization on global object
fresh_2cells_merge <- fresh_2cells_merge %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>%
  SCTransform(vars.to.regress = c("percent_mito"))

# Calculate PCs using variable features determined by SCTransform (3000 by default)
fresh_2cells_merge <- RunPCA(fresh_2cells_merge, assay = "SCT", npcs = 50) 

fresh_2cells_merge_harmonized <- RunHarmony(fresh_2cells_merge, 
                                        group.by.vars = c("original.batch"),
                                        reduction = "pca", assay.use = "SCT", reduction.save = "harmony",
                                        theta = 6) #umicounts used 6/readscounts_exon_all used 6
ElbowPlot(fresh_2cells_merge_harmonized)
fresh_2cells_merge_harmonized <- RunUMAP(fresh_2cells_merge_harmonized, reduction = "harmony", assay = "SCT", dims = 1:8) #reads_acount
fresh_2cells_merge_harmonized <- RunUMAP(fresh_2cells_merge_harmonized, reduction = "harmony", assay = "SCT", dims = 1:5) #umi_account;nfeature > 980,dim=5  

fresh_2cells_merge_harmonized <- FindNeighbors(object = fresh_2cells_merge_harmonized, reduction = "harmony")
fresh_2cells_merge_harmonized <- FindClusters(fresh_2cells_merge_harmonized, resolution =  0.2) #if nfeature > 980; umi used 0.2
DimPlot(fresh_2cells_merge_harmonized, group.by = c("original.batch", "ident", "original.ID"), ncol = 3)
DimPlot(fresh_2cells_merge_harmonized,reduction="pca", group.by = c("original.batch", "ident", "original.ID"), ncol = 3)

#save dot position in Umap
fresh_2cells_merge_harmonized_umap = fresh_2cells_merge_harmonized@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(clusters = fresh_2cells_merge_harmonized@meta.data$seurat_clusters) %>% cbind(originalbatch = fresh_2cells_merge_harmonized$original.batch) %>% cbind(originalID = fresh_2cells_merge_harmonized$original.ID)
#fresh_2cells_merge_harmonized = fresh_2cells_merge_harmonized@reductions$harmony@cell.embeddings %>% as.data.frame() %>% cbind(clusters = fresh_2cells_merge_harmonized@meta.data$seurat_clusters) %>% cbind(originalbatch = fresh_2cells_merge_harmonized$original.batch) %>% cbind(originalbatch = fresh_2cells_merge_harmonized$original.ID)
write.table(fresh_2cells_merge_harmonized_umap, file = "/Users/hxr/Desktop/WJH/fresh_proM_G1/fresh_Hap1&K562_umap.csv",quote=F) 
#write.table(fresh_2cells_merge_harmonized, file = "/Users/hxr/Desktop/WJH/fresh_proM_G1/fresh_Hap1&K562_harmony.csv",quote=F) 

#get cluster info of each cells and save into a file
cluster_identities <- Idents(object = fresh_2cells_merge_harmonized)
ident_labels <- as.matrix(unlist(cluster_identities))
write.table(ident_labels, file = "/Users/hxr/Desktop/WJH/fresh_proM_G1/cluster_ident.txt", sep = "\t", row.names = TRUE, col.names = F)


##plot the feature gene---- 
S_gene <-c("ORC1","CCNE1","CCNE2","MCM6")
G2M_gene <-c("WEE1","CDK1","CCN1","NUSAP1","AURKA","CCNA2","CCNB2")

#FeaturePlot(fresh_merge_harmonized,features = G2M_gene)
VlnPlot(fixlysis_merge_harmonized,features=G2M_gene)
VlnPlot(fixlysis_merge_harmonized,features=S_gene)
FeaturePlot(fixlysis_merge_harmonized,features=G2M_gene)
FeaturePlot(fixlysis_merge_harmonized,features=S_gene)
#check the absolute numbers
VlnPlot(fixlysis_merge_harmonized, features = G2M_gene, slot = "counts", log = TRUE) 
VlnPlot(fixlysis_merge_harmonized, features = S_gene, slot = "counts", log = TRUE)
#check the mito/ribo and nCount/nFeature
VlnPlot(fixlysis_merge_harmonized, features = c("percent_mito",'percent_ribo'),slot = "counts", log = TRUE)
VlnPlot(fixlysis_merge_harmonized, features = c("nCount_RNA", "nFeature_RNA"))



##transfer to Densmap by umap-learn implement into Seurat (discard)---- 
## (Not suitable for data with batch effect)
## Have to install Seurat_develop version and activate reticulate package:
#install.packages('remotes')
#remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')
#library(Seurat)
#codes:
#Fresh_2cells_dens <-Seurat::RunUMAP(fresh_2cells_merge_harmonized, dims = 1:5, umap.method = "umap-learn", densmap=TRUE)
#Fresh_2cells_dens <- FindNeighbors(object = Fresh_2cells_dens, reduction = "harmony")
#Fresh_2cells_dens <- FindClusters(Fresh_2cells_dens, resolution =  0.3)
#DimPlot(Fresh_2cells_dens, group.by = c("original.batch", "ident", "original.ID"), ncol = 3)

#set the enviroment and python into specific path
#Sys.setenv(RETICULATE_PYTHON = "/Users/hxr/miniconda3/envs/r-reticulate-python38/bin/python")
#library(reticulate)
#py_config()

##transfer to Densmap do harmony first then Umap plot (discard)---- 
#BiocManager::install("densvis")
#library(densvis)
#library(tensorflow)
#Fresh_2cells_dens <- densmap(fresh_2cells_merge_harmonized, reduction.use = "harmony", verbose = FALSE)
#fresh_2cells_harmony <- RunHarmony(fresh_2cells_merge,assay.use = "SCT", group.by.vars = "original.batch")
#fresh_2cells_harmony <- RunPCA(fresh_2cells_harmony, features = VariableFeatures(fresh_2cells_harmony))


#export raw matrix out and switch to python to plot densmap---- 
#I choosed orginal matrix, assay = 'SCT', slot = 'counts' finally 
#(assay = 'SCT' means, SCT normalized and log-transformed expressed level; slot = 'counts' means not undergone further downstream processing or dimensionality reduction; comparing with slot = 'data', which undergone normalization, scaling, and dimensionality reduction )
reduced_data <- GetAssayData(fresh_2cells_harmony, assay = "SCT", slot = "counts") #use original matrix finally
#reduced_data <- fresh_2cells_harmony@reductions$harmony@cell.embeddings #use harmonied matrix
write.table(reduced_data, file = "/Users/hxr/Desktop/WJH/fresh_proM_G1/SCT_matrix.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


#get batch1 normalized dataset (play for fun)
#cannot separated into two clusters if cell<96 in Densmap
#fresh_2cells_merge_harmonized.batch1 <- subset(fresh_2cells_merge_harmonized, subset = original.batch == "1")
#reduced_data <- GetAssayData(fresh_2cells_merge_harmonized.batch1, assay = "RNA", slot = "counts") #use original matrix
#write.table(reduced_data, file = "/Users/hxr/Desktop/WJH/fresh_proM_G1/SCT_matrix_batch1_counts.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

