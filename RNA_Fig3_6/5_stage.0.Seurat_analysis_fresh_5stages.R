---
  #title: "1.run-Seurat"
  #author: "Yu Fu (QC) & Xiangru (batch effect remove, imputation, pseudotime,variances)" 
  #date: "2023-02-23"
  #R version: 3.6.3
  ---
  
  ##load.package ----
library(dplyr) #version = package_version('1.0.10')
library(Seurat) #version = package_version('3.2.3') later: 4.0.1
library(patchwork)  #version = package_version('1.1.2')
library(ggplot2)  #version = package_version('3.4.0')
library(biomaRt) #version = package_version('2.42.1')
library(stringr)  #version = package_version('1.5.0')

##load sup H1-mHep dataset and gerenate matrix ----
H1_100k <- readRDS("/Users/hxr/Desktop/SSC_scripts/4DN/4DN_5stage_fresh_sup/4DN_H1.dgecounts.rds")
gene_name <- read.csv("/Users/hxr/Desktop/SSC_scripts/4DN/4DN_5stage_fresh_sup/4DN_H1.gene_names.txt",sep="\t",header = T)
reads <-as.matrix(H1_100k$readcount$inex$all)
rownames(reads) <- gene_name$gene_name[match(rownames(H1_100k$readcount$inex$all),gene_name$gene_id)] ## change gene id as gene name
H1_sup <- Seurat::CreateSeuratObject(counts=reads,project="H1_sup")

DE_100k <- readRDS("/Users/hxr/Desktop/SSC_scripts/4DN/4DN_5stage_fresh_sup/4DN_DE.dgecounts.rds")
reads <-as.matrix(DE_100k$readcount$inex$all)
rownames(reads) <- gene_name$gene_name[match(rownames(DE_100k$readcount$inex$all),gene_name$gene_id)] ## change gene id as gene name
DE_sup <- Seurat::CreateSeuratObject(counts=reads,project="DE_sup")

HB_100k <- readRDS("/Users/hxr/Desktop/SSC_scripts/4DN/4DN_5stage_fresh_sup/4DN_HB.dgecounts.rds")
reads <-as.matrix(HB_100k$readcount$inex$all)
rownames(reads) <- gene_name$gene_name[match(rownames(HB_100k$readcount$inex$all),gene_name$gene_id)] ## change gene id as gene name
HB_sup <- Seurat::CreateSeuratObject(counts=reads,project="HB_sup")

iHep_100k <- readRDS("/Users/hxr/Desktop/SSC_scripts/4DN/4DN_5stage_fresh_sup/4DN_iHep.dgecounts.rds")
reads <-as.matrix(iHep_100k$readcount$inex$all)
rownames(reads) <- gene_name$gene_name[match(rownames(iHep_100k$readcount$inex$all),gene_name$gene_id)] ## change gene id as gene name
iHep_sup <- Seurat::CreateSeuratObject(counts=reads,project="iHep_sup")

mHep_100k <- readRDS("/Users/hxr/Desktop/SSC_scripts/4DN/4DN_5stage_fresh_sup/4DN_mHep.dgecounts.rds")
reads <-as.matrix(mHep_100k$readcount$inex$all)
rownames(reads) <- gene_name$gene_name[match(rownames(mHep_100k$readcount$inex$all),gene_name$gene_id)] ## change gene id as gene name
mHep_sup <- Seurat::CreateSeuratObject(counts=reads,project="mHep_sup")

#Merge raw samples & load barcodes----
differentiation_sup_list <- list(H1_sup,DE_sup,HB_sup,iHep_sup,mHep_sup)
differentiation_sup <- Reduce(function(x, y) merge(x, y, merge.data = TRUE), differentiation_sup_list)
barcodes <- read.csv("/Users/hxr/Desktop/SSC_scripts/4DN/4DN_5stage_fresh_sup/5stage_sup_reads_for_zUMIs.samples.txt",sep="\t")

##check the objects ----
head(differentiation_sup@assays$RNA@counts[1:10,1:10])
rownames(differentiation_sup@assays$RNA@counts[1:10,1:10])
colnames(differentiation_sup@assays$RNA@counts[1:10,1:10])
## save the matrix
reads <-as.matrix(differentiation_sup@assays$RNA@counts)
write.table(reads, file = "/Users/hxr/Desktop/SSC_scripts/4DN/4DN_fresh_5stages_reads_sup_matrix.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#add original.ID to objects----
rownames(barcodes) <- barcodes[,4]
differentiation_sup$original.ID <- NA
differentiation_sup$original.ID[rownames(barcodes)] <- as.character(barcodes[,5])
#add original.batch to objects
differentiation_sup$original.batch <- NA
differentiation_sup$original.batch[rownames(barcodes)] <- as.character(barcodes[,6])
#add original.plate to objects
differentiation_sup$original.plate <- NA
differentiation_sup$original.plate[rownames(barcodes)] <- as.character(barcodes[,7])

#change the colnames of differentiation_sup as sample ID ----
matching_indices <- match(colnames(differentiation_sup@assays$RNA@counts),barcodes$BC)
colnames(differentiation_sup@assays$RNA@counts)[!is.na(matching_indices)] <- barcodes$sample[matching_indices[!is.na(matching_indices)]] ## change cell as sample id
colnames(differentiation_sup@assays$RNA@data)[!is.na(matching_indices)] <- barcodes$sample[matching_indices[!is.na(matching_indices)]] ## change cell as sample id
colnames(differentiation_sup@assays$RNA@counts[1:10,1:10])
colnames(differentiation_sup@assays$RNA@data[1:10,1:10])


##load old dataset and gerenate matrix ----
differentiation_100k <- readRDS("/Users/hxr/Desktop/SSC_scripts/4DN/4DN_5stages.dgecounts.rds")
barcodes <- read.csv("/Users/hxr/Desktop/SSC_scripts/4DN/4DN_fresh_5stages_reads_for_zUMIs.samples.txt",sep="\t")
gene_name <- read.csv("/Users/hxr/Desktop/SSC_scripts/4DN/4DN_5stages.gene_names.txt",sep="\t",header = T)
##choose to use differentiation_100k$readcount to creat seurat objects (I used readcount)
reads <-as.matrix(differentiation_100k$readcount$inex$all)
rownames(reads) <- gene_name$gene_name[match(rownames(differentiation_100k$readcount$inex$all),gene_name$gene_id)] ## change gene id as gene name
colnames(reads) <- barcodes$sample[match(colnames(differentiation_100k$readcount$inex$all),barcodes$BC)] ## change cell as sample id

## transfer as seurat objects
differentiation_old <- Seurat::CreateSeuratObject(counts=reads,project="differentiation_old")

##check the objects 
head(differentiation_old@assays$RNA@counts[1:10,1:10])
rownames(differentiation_old@assays$RNA@counts[1:10,1:10])
colnames(differentiation_old@assays$RNA@counts[1:10,1:10])

#merge two batch of matrix together----
count1 <- as.matrix(differentiation_old@assays$RNA@counts)
count2 <- as.matrix(differentiation_sup@assays$RNA@counts)
all_genes <- union(rownames(count1),rownames(count2))
all_barcodes <- union(colnames(count1),colnames(count2))
count1_full <- Matrix::Matrix(0,nrow = length(all_genes),ncol = length(all_barcodes),sparse = TRUE)
rownames(count1_full) <- all_genes
colnames(count1_full) <- all_barcodes
count2_full <-count1_full
count1_full[rownames(count1),colnames(count1)] <- count1
count2_full[rownames(count2),colnames(count2)] <- count2
merge_counts <- count1_full + count2_full
differentiation_all <- CreateSeuratObject(counts = merge_counts)
## save the matrix
write.table(merge_counts, file = "/Users/hxr/Desktop/SSC_scripts/4DN/4DN_fresh_5stages_reads_all_matrix.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#prepare barcode infomations----
barcodes_old <- read.csv("/Users/hxr/Desktop/SSC_scripts/4DN/4DN_fresh_5stages_reads_for_zUMIs.samples.txt",sep="\t")
rownames(barcodes_old) <- barcodes_old[,3]
barcodes_new <- read.csv("/Users/hxr/Desktop/SSC_scripts/4DN/4DN_5stage_fresh_sup/5stage_sup_reads_for_zUMIs.samples.txt",sep="\t")
rownames(barcodes_new) <- barcodes_new[,3]
new_rownames_to_add <- setdiff(rownames(barcodes_new),rownames(barcodes_old))
barcodes_new_selected <- barcodes_new[new_rownames_to_add,,drop = FALSE]
colnames(barcodes_new_selected) <- colnames(barcodes_old)
barcodes_all <- rbind(barcodes_old,barcodes_new_selected)
## save the barcode
write.table(barcodes_all, file = "/Users/hxr/Desktop/SSC_scripts/4DN/4DN_fresh_5stages_barcode_all.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
#add original.ID----
differentiation_all$original.ID <- NA
differentiation_all$original.ID[rownames(barcodes_all)] <- as.character(barcodes_all[,5])
differentiation_all$original.ID<-factor(differentiation_all$original.ID, levels = c('H1','DE','HB','iHep','mHep'))
#add original.batch to differentiation_all objects
differentiation_all$original.batch <- NA
differentiation_all$original.batch[rownames(barcodes_all)] <- as.character(barcodes_all[,6])
#add original.plate to differentiation_all objects
differentiation_all$original.plate <- NA
differentiation_all$original.plate[rownames(barcodes_all)] <- as.character(barcodes_all[,7])

## calculate percentage contamination ----
differentiation_all <-  PercentageFeatureSet(differentiation_all,'^MT-',col.name="percent_mito")
differentiation_all <-  PercentageFeatureSet(differentiation_all,'^RP[SL]',col.name="percent_ribo")
ribo_genes <- rownames(differentiation_all)[grep("^RP[SL]",rownames(differentiation_all))]
features <- c("percent_ribo","percent_mito")
#"nCount_RNA","nFeature_RNA"
VlnPlot(differentiation_all,group.by='orig.ident',features=features,pt.size=0.1,ncol=2)
VlnPlot(differentiation_all,group.by='orig.ident',features=c("nCount_RNA","nFeature_RNA"),pt.size=0.1,ncol=2)
FeatureScatter(differentiation_all,"nCount_RNA","nFeature_RNA",group.by="original.ID",pt.size=0.5)
FeatureScatter(differentiation_all,"nFeature_RNA","percent_ribo")
FeatureScatter(differentiation_all,"nFeature_RNA","percent_mito")
#save the "nCountRNA" VS "nFeatureRNA" 
scatter_data <- data.frame(nCount_RNA = FetchData(differentiation_all,"nCount_RNA")$nCount_RNA, nFeature_RNA = FetchData(differentiation_all,"nFeature_RNA")$nFeature_RNA, original_ID = differentiation_all@meta.data$original.ID)
write.table(scatter_data, file = "/Users/hxr/Desktop/5stage_scRNA/4DN_5stage_ReadsVSgenecounts_nofilter.csv",quote=F) 

#remove ribo gene
differentiation_all_ribo <- subset(differentiation_all, features = rownames(differentiation_all)[!(rownames(differentiation_all) %in% ribo_genes)])
dim(differentiation_all)
dim(differentiation_all_ribo)
#only select cells with > a features and genes in > b cells ----
#selecteda <- WhichCells(differentiation_all,expression=(nFeature_RNA < 10000)) #readcounts>980
selecteda <- WhichCells(differentiation_all_ribo,expression=nFeature_RNA>200) #umicounts>700
selectedb <- rownames(differentiation_all_ribo)[Matrix::rowSums(differentiation_all_ribo)>1] #readcounts_exon>2; umicounts_inex>1
differentiation_all.filt1 <- subset(differentiation_all_ribo,features=selectedb,cells=selecteda)
dim(differentiation_all)
dim(differentiation_all.filt1)
VlnPlot(differentiation_all.filt1,group.by='orig.ident',features=c("nCount_RNA","nFeature_RNA"),pt.size=0.1,ncol=2)
VlnPlot(differentiation_all.filt1, features = c("percent_mito",'percent_ribo'),slot = "counts", log = TRUE)

#set cells with < 5000 features ---
selectedc <- WhichCells(differentiation_all.filt1,expression=nFeature_RNA<11000) #genecounts<3000 remove the tail in violin plot
differentiation_all.filt2 <- subset(differentiation_all.filt1,cells=selectedc)
dim(differentiation_all.filt2)
VlnPlot(differentiation_all.filt2,group.by='orig.ident',features=c("nCount_RNA","nFeature_RNA"),pt.size=0.1,ncol=2)
selectedd <- WhichCells(differentiation_all.filt2,expression=nCount_RNA<2000000) #genecounts<3000 remove the tail in violin plot
differentiation_all.filt3 <- subset(differentiation_all.filt2,cells=selectedd)
dim(differentiation_all.filt3)
VlnPlot(differentiation_all.filt3,group.by='orig.ident',features=c("nCount_RNA","nFeature_RNA"),pt.size=0.1,ncol=2)
VlnPlot(differentiation_all.filt3, features = c("percent_mito",'percent_ribo'),slot = "counts", log = TRUE)

#set cells with < 30% ribo (ES stage has many ribo) ---
selectedc <- WhichCells(differentiation_all.filt3,expression=percent_ribo<35) #percent_ribo <30% remove the tail in violin plot
differentiation_all.filt4 <- subset(differentiation_all.filt3,cells=selectedc)
dim(differentiation_all.filt4)
VlnPlot(differentiation_all.filt4,group.by='orig.ident',features=c("nCount_RNA","nFeature_RNA"),pt.size=0.1,ncol=2)
VlnPlot(differentiation_all.filt4, features = c("percent_mito",'percent_ribo'),slot = "counts", log = TRUE)
selectedd <- WhichCells(differentiation_all.filt4,expression=percent_ribo>0.0001) #percent_mito >1% remove the tail in violin plot
differentiation_all.filt <- subset(differentiation_all.filt4,cells=selectedd)
dim(differentiation_all.filt)
VlnPlot(differentiation_all.filt,group.by='orig.ident',features=c("nCount_RNA","nFeature_RNA"),pt.size=0.1,ncol=2)
VlnPlot(differentiation_all.filt, group.by='orig.ident',features = c("percent_mito",'percent_ribo'),slot = "counts", log = TRUE)
VlnPlot(differentiation_all.filt, group.by='original.ID',features=c("nCount_RNA","nFeature_RNA"),slot = "counts", log = TRUE)

#save the violinplot into pdf
#pdf("/Users/hxr/Desktop/5stage_scRNA/4DN_5stage_ReadsVSfeaturegene.pdf", width = 8, height = 6)  # Adjust size as needed
#VlnPlot(differentiation_all.filt, group.by = 'original.ID', 
#        features = c("nCount_RNA", "nFeature_RNA"), 
#        slot = "counts", log = TRUE)
#dev.off()

# see the top 10 genes with mapped reads ----
par(mar=c(4,10,2,1))
c <- differentiation_all.filt@assays$RNA@counts
c <- Matrix::t(Matrix::t(c)/Matrix::colSums(c)) * 100
most_expressed <- order(apply(c,1,median), decreasing=T)[20:1]
boxplot(as.matrix(t(as.matrix(c[most_expressed,]))),cex=0.1, las=1, xlab="% total count per cell",col=(scales::hue_pal())(20)[20:1],horizontal=TRUE)

## Nomalize+feature gene+runPCA ----
differentiation_all.filt <- NormalizeData(differentiation_all.filt,normalization.method ="LogNormalize",scale.factor=10000)
differentiation_all.filt <- FindVariableFeatures(differentiation_all.filt, selection.method="vst",nfeatures=5000)
top10 <- head(VariableFeatures(differentiation_all.filt),10)
p1 <- VariableFeaturePlot(differentiation_all.filt)
p2 <- LabelPoints(plot=p1,points=top10,repel=TRUE)
p2
#manually remove 'LINC00486' the outlier genes
differentiation_all.filt <- subset(differentiation_all.filt, features = setdiff(rownames(differentiation_all.filt),'LINC00486'))
p1 <- VariableFeaturePlot(differentiation_all.filt)
p2 <- LabelPoints(plot=p1,points=top10,repel=TRUE)
p2

#run PCA
all.genes <- rownames(differentiation_all.filt)
differentiation_all.filt <- ScaleData(differentiation_all.filt, features=all.genes,verbose = FALSE)
differentiation_all.filt <- RunPCA(differentiation_all.filt,npcs=40)
print(differentiation_all.filt[["pca"]],dims=1:5,nfeatures=5)
#VizDimLoadings(differentiation_all.filt,dims=1:2,reduction="pca")
DimPlot(differentiation_all.filt,reduction="pca",group.by="original.ID")
DimPlot(differentiation_all.filt,reduction="pca",group.by="ident")
FeaturePlot(differentiation_all.filt,reduction="pca",features='nFeature_RNA')
FeaturePlot(differentiation_all.filt,reduction="pca",features='nCount_RNA')
DimPlot(differentiation_all.filt,reduction="pca",group.by="original.batch")
DimPlot(differentiation_all.filt,reduction="pca",group.by="original.plate")

#save the umap dot informations into csv----
differentiation_pca= differentiation_all.filt@reductions$pca@cell.embeddings %>% as.data.frame() %>% cbind(clusters = differentiation_all.filt@active.ident) %>% cbind(originalbatch = differentiation_all.filt$original.batch) %>% cbind(originalID = differentiation_all.filt$original.ID)%>%cbind(Genecount = differentiation_all.filt$nFeature_RNA)%>%cbind(Readscount = differentiation_all.filt$nCount_RNA)
write.table(differentiation_pca, file = "/Users/hxr/Desktop/5stage_scRNA/4DN_5stage_pca.csv",quote=F) 


##determine the dimensionality----
differentiation_all.filt <- JackStraw(differentiation_all.filt,num.replicate=100)
ElbowPlot(differentiation_all.filt)
plot(differentiation_all.filt@reductions$pca@stdev)

##reduce dimention and plot----
differentiation_all.filt <- FindNeighbors(differentiation_all.filt, dims = 1:20) #reads_used
differentiation_all.filt <- FindClusters(differentiation_all.filt, resolution = 0.25)
differentiation_all.filt <- RunUMAP(differentiation_all.filt, reduction = "pca", dims = 1:20,min.dist = 0.3)
p1 <- DimPlot(differentiation_all.filt,reduction = "umap",label=TRUE,group.by="ident")
p1
p2 <- DimPlot(differentiation_all.filt,reduction = "umap",label=TRUE,group.by="original.ID")
p2
p3 <- DimPlot(differentiation_all.filt,reduction = "umap",label=TRUE,group.by="original.batch")
p3
DimPlot(differentiation_all.filt,reduction = "umap",label=TRUE,group.by="original.plate")

#save the umap dot informations into csv----
p1_umap = differentiation_all.filt@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(clusters = differentiation_all.filt@active.ident) %>% cbind(originalbatch = differentiation_all.filt$original.batch) %>% cbind(originalID = differentiation_all.filt$original.ID)
write.table(p1_umap, file = "/Users/hxr/Desktop/5stage_scRNA/4DN_5stage_umap.csv",quote=F) 


##plot the feature gene---- 
H1_gene <-c("OCT4","NANOG","SOX2",'ZFP42','KLF4','ESRRB','LIN28','NR5A2')
DE_gene <-c("SOX17",'FOXA2','GATA4','LEFTY1','GATA6','CXCR4','EOMES')
mHep_gene <-c("ALB","HNF4A",'AFP','CK8','TBX','FOXA2','TTR','SERPINA1')

#FeaturePlot(features = H1_gene)
VlnPlot(differentiation_all.filt,features=H1_gene)
VlnPlot(differentiation_all.filt,features=DE_gene)
VlnPlot(differentiation_all.filt,features=mHep_gene)
p1 <-FeaturePlot(differentiation_all.filt,features = H1_gene[2],cols  = c("lightgrey" ,"#2C6FD8"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#2C6FD8', high = 'blue')
p2 <-FeaturePlot(differentiation_all.filt,features = H1_gene[3],cols  = c("lightgrey" ,"#2C6FD8"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#2C6FD8', high = 'blue')
p3 <-FeaturePlot(differentiation_all.filt,features = H1_gene[4],cols  = c("lightgrey" ,"#2C6FD8"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = '#2C6FD8', high = 'blue')
p4 <-FeaturePlot(differentiation_all.filt,features = H1_gene[5],cols  = c("lightgrey" ,"#2C6FD8"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = '#2C6FD8', high = 'blue')
p5 <-FeaturePlot(differentiation_all.filt,features = H1_gene[8],cols  = c("lightgrey" ,"#2C6FD8"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = '#2C6FD8', high = 'blue')
p1+p2
p3+p4+p5
#FeaturePlot(features = DE_gene)
p6 <- FeaturePlot(differentiation_all.filt,features = DE_gene[1],cols  = c("lightgrey" ,"#199F62"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#199F62', high = 'green')
p7 <- FeaturePlot(differentiation_all.filt,features = DE_gene[2],cols  = c("lightgrey" ,"#199F62"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#199F62', high = 'green')
p7_2 <- FeaturePlot(differentiation_all.filt,features = DE_gene[3],cols  = c("lightgrey" ,"#199F62"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#199F62', high = 'green')
p7_4 <- FeaturePlot(differentiation_all.filt,features = DE_gene[5],cols  = c("lightgrey" ,"#199F62"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#199F62', high = 'green')
p7_5 <- FeaturePlot(differentiation_all.filt,features = DE_gene[6],cols  = c("lightgrey" ,"#199F62"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#199F62', high = 'green')
p7_6 <- FeaturePlot(differentiation_all.filt,features = DE_gene[7],cols  = c("lightgrey" ,"#199F62"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#199F62', high = 'green')
p7_7 <- FeaturePlot(differentiation_all.filt,features = DE_gene[8],cols  = c("lightgrey" ,"#199F62"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#199F62', high = 'green')
p6+p7+p7_2
p7_4+p7_5+p7_6+p7_7
#FeaturePlot(features = mHep_gene)
p8 <- FeaturePlot(differentiation_all.filt,features = mHep_gene[1],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p9 <- FeaturePlot(differentiation_all.filt,features = mHep_gene[2],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p10 <- FeaturePlot(differentiation_all.filt,features = mHep_gene[3],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p11 <- FeaturePlot(differentiation_all.filt,features = mHep_gene[6],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p12 <- FeaturePlot(differentiation_all.filt,features = mHep_gene[7],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p13 <- FeaturePlot(differentiation_all.filt,features = mHep_gene[8],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p8
p9+p10
p11+p12+p13

#plot DEpHBn or DEnHBp genes
DEpHBn <- c("LEFTY1","LEFT2","CER1")
DEnHBp <- c("AFP","FGA","FGG","FGB",'APOA1','APOA2','APOA4','APOA5','APOC3')

p <-FeaturePlot(differentiation_all.filt,features = DEpHBn,cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.3)
p+scale_color_gradient(trans = 'log', low = 'orange', high = 'yellow')

p1 <- FeaturePlot(differentiation_all.filt,features = DEnHBp[1],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p2 <- FeaturePlot(differentiation_all.filt,features = DEnHBp[2],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p3 <- FeaturePlot(differentiation_all.filt,features = DEnHBp[3],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p4 <- FeaturePlot(differentiation_all.filt,features = DEnHBp[4],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p5 <- FeaturePlot(differentiation_all.filt,features = DEnHBp[5],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p6 <- FeaturePlot(differentiation_all.filt,features = DEnHBp[6],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p7 <- FeaturePlot(differentiation_all.filt,features = DEnHBp[7],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p8 <- FeaturePlot(differentiation_all.filt,features = DEnHBp[8],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p9 <- FeaturePlot(differentiation_all.filt,features = DEnHBp[9],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')

p1+p2+p3
p4+p5+p6
p7+p8+p9

#save the feature gene information into umap----
all_features <- unique(c(H1_gene,DE_gene, mHep_gene, DEpHBn, DEnHBp))
expr_values <- FetchData(differentiation_all.filt, vars = all_features)
p1_umap <- cbind(p1_umap, expr_values)
write.table(p1_umap, file = "/Users/hxr/Desktop/5stage_scRNA/4DN_5stage_umap_featuregene.csv",quote=F) 


##Select cells from DE to mHep---- 
DE_mHep_obj <- subset(differentiation_all.filt, subset = original.ID %in% c('DE','HB','iHep','mHep'))
DE_mHep_obj <- NormalizeData(DE_mHep_obj,normalization.method ="LogNormalize",scale.factor=10000)
DE_mHep_obj <- FindVariableFeatures(DE_mHep_obj, selection.method="vst",nfeatures=5000)
top10 <- head(VariableFeatures(DE_mHep_obj),10)
p1 <- VariableFeaturePlot(DE_mHep_obj)
p2 <- LabelPoints(plot=p1,points=top10,repel=TRUE)
p2

all.genes <- rownames(DE_mHep_obj)
DE_mHep_obj <- ScaleData(DE_mHep_obj, features=all.genes,verbose = FALSE)
DE_mHep_obj <- RunPCA(DE_mHep_obj,npcs=60)
print(DE_mHep_obj[["pca"]],dims=1:5,nfeatures=5)
#VizDimLoadings(DE_mHep_obj,dims=1:2,reduction="pca")
DimPlot(DE_mHep_obj,reduction="pca",group.by="original.ID")
DimPlot(DE_mHep_obj,reduction="pca",group.by="ident")
FeaturePlot(DE_mHep_obj,reduction="pca",features='nFeature_RNA')
FeaturePlot(DE_mHep_obj,reduction="pca",features='nCount_RNA')
DimPlot(DE_mHep_obj,reduction="pca",group.by="original.batch")
DimPlot(DE_mHep_obj,reduction="pca",group.by="original.plate")
DimHeatmap(DE_mHep_obj,dims=1:5,balanced=TRUE)

##determine the dimensionality----
DE_mHep_obj <- JackStraw(DE_mHep_obj,num.replicate=100)
ElbowPlot(DE_mHep_obj)
plot(DE_mHep_obj@reductions$pca@stdev)

##reduce dimention and plot----
DE_mHep_obj <- FindNeighbors(DE_mHep_obj, dims = 1:39) #reads_used
DE_mHep_obj <- FindClusters(DE_mHep_obj, resolution = 0.25)
DE_mHep_obj <- RunUMAP(DE_mHep_obj, reduction = "pca", dims = 1:39,min.dist = 0.3)
p1 <- DimPlot(DE_mHep_obj,reduction = "umap",label=TRUE,group.by="ident")
p1
p2 <- DimPlot(DE_mHep_obj,reduction = "umap",label=TRUE,group.by="original.ID")
p2
p3 <- DimPlot(DE_mHep_obj,reduction = "umap",label=TRUE,group.by="original.batch")
p3
DimPlot(DE_mHep_obj,reduction = "umap",label=TRUE,group.by="original.plate")

#save Umap dot position into csv
DEtomHep_umap = DE_mHep_obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(clusters = DE_mHep_obj@active.ident) %>% cbind(originalbatch = DE_mHep_obj$original.batch) %>% cbind(originalID = DE_mHep_obj$original.ID)
write.table(DEtomHep_umap, file = "/Users/hxr/Desktop/5stage_scRNA/4DN_4stage_umap.csv",quote=F) 

##plot the feature gene---- 
#FeaturePlot(features = DE_gene)
p6 <- FeaturePlot(DE_mHep_obj,features = DE_gene[1],cols  = c("lightgrey" ,"#199F62"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#199F62', high = 'green')
p7 <- FeaturePlot(DE_mHep_obj,features = DE_gene[2],cols  = c("lightgrey" ,"#199F62"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#199F62', high = 'green')
p6+p7
#FeaturePlot(features = mHep_gene)
p8 <- FeaturePlot(DE_mHep_obj,features = mHep_gene[1],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p9 <- FeaturePlot(DE_mHep_obj,features = mHep_gene[2],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p10 <- FeaturePlot(DE_mHep_obj,features = mHep_gene[3],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p11 <- FeaturePlot(DE_mHep_obj,features = mHep_gene[6],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p8
p9+p10

#plot DEpHBn or DEnHBp genes
#DE+HB-
p <-FeaturePlot(DE_mHep_obj,features = DEpHBn,cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.3)
p+scale_color_gradient(trans = 'log', low = 'orange', high = 'yellow')
#DE-HB+
p1 <- FeaturePlot(DE_mHep_obj,features = DEnHBp[1],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p2 <- FeaturePlot(DE_mHep_obj,features = DEnHBp[2],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p3 <- FeaturePlot(DE_mHep_obj,features = DEnHBp[3],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p4 <- FeaturePlot(DE_mHep_obj,features = DEnHBp[4],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p5 <- FeaturePlot(DE_mHep_obj,features = DEnHBp[5],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p6 <- FeaturePlot(DE_mHep_obj,features = DEnHBp[6],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p7 <- FeaturePlot(DE_mHep_obj,features = DEnHBp[7],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p8 <- FeaturePlot(DE_mHep_obj,features = DEnHBp[8],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p9 <- FeaturePlot(DE_mHep_obj,features = DEnHBp[9],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')

p1+p2+p3
p4+p5+p6
p7+p8+p9

#save these feature genes information into umap----
all_features <- unique(c(DE_gene, mHep_gene, DEpHBn, DEnHBp))
expr_values <- FetchData(DE_mHep_obj, vars = all_features)
DEtomHep_umap <- cbind(DEtomHep_umap, expr_values)
write.table(DEtomHep_umap, file = "/Users/hxr/Desktop/5stage_scRNA/4DN_4stage_umap_featuregene.csv", quote = FALSE, sep = ",", row.names = TRUE)

##trajectory from HB to mHep
##Select cells from HB to mHep (not used)---- 
HB_mHep_obj <- subset(differentiation_all.filt, subset = original.ID %in% c('HB','iHep','mHep'))
HB_mHep_obj <- NormalizeData(HB_mHep_obj,normalization.method ="LogNormalize",scale.factor=10000)
HB_mHep_obj <- FindVariableFeatures(HB_mHep_obj, selection.method="vst",nfeatures=5000)
top10 <- head(VariableFeatures(HB_mHep_obj),10)
p1 <- VariableFeaturePlot(HB_mHep_obj)
p2 <- LabelPoints(plot=p1,points=top10,repel=TRUE)
p2

all.genes <- rownames(HB_mHep_obj)
HB_mHep_obj <- ScaleData(HB_mHep_obj, features=all.genes,verbose = FALSE)
HB_mHep_obj <- RunPCA(HB_mHep_obj,npcs=60)
print(HB_mHep_obj[["pca"]],dims=1:5,nfeatures=5)
#VizDimLoadings(DE_mHep_obj,dims=1:2,reduction="pca")
DimPlot(HB_mHep_obj,reduction="pca",group.by="original.ID")
DimPlot(HB_mHep_obj,reduction="pca",group.by="ident")
FeaturePlot(HB_mHep_obj,reduction="pca",features='nFeature_RNA')
FeaturePlot(HB_mHep_obj,reduction="pca",features='nCount_RNA')
DimPlot(HB_mHep_obj,reduction="pca",group.by="original.batch")
DimPlot(HB_mHep_obj,reduction="pca",group.by="original.plate")

##determine the dimensionality----
HB_mHep_obj <- JackStraw(HB_mHep_obj,num.replicate=100)
ElbowPlot(HB_mHep_obj)
plot(HB_mHep_obj@reductions$pca@stdev)

##reduce dimention and plot----
HB_mHep_obj <- FindNeighbors(HB_mHep_obj, dims = 1:40) #reads_used
HB_mHep_obj <- FindClusters(HB_mHep_obj, resolution = 0.25)
HB_mHep_obj <- RunUMAP(HB_mHep_obj, reduction = "pca", dims = 1:40,min.dist = 0.3)
p1 <- DimPlot(HB_mHep_obj,reduction = "umap",label=TRUE,group.by="ident")
p1
p2 <- DimPlot(HB_mHep_obj,reduction = "umap",label=TRUE,group.by="original.ID")
p2
p3 <- DimPlot(HB_mHep_obj,reduction = "umap",label=TRUE,group.by="original.batch")
p3
DimPlot(HB_mHep_obj,reduction = "umap",label=TRUE,group.by="original.plate")

##plot the feature gene---- 
#FeaturePlot(features = mHep_gene)
p8 <- FeaturePlot(HB_mHep_obj,features = mHep_gene[1],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p9 <- FeaturePlot(HB_mHep_obj,features = mHep_gene[2],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p10 <- FeaturePlot(HB_mHep_obj,features = mHep_gene[3],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p11 <- FeaturePlot(HB_mHep_obj,features = mHep_gene[6],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.1)+scale_color_gradient(trans = 'log', low = '#E97D30', high = 'yellow')
p8
p9+p10

#plot DEpHBn or DEnHBp genes
#DE-HB+
p1 <- FeaturePlot(HB_mHep_obj,features = DEnHBp[1],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p2 <- FeaturePlot(HB_mHep_obj,features = DEnHBp[2],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p3 <- FeaturePlot(HB_mHep_obj,features = DEnHBp[3],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p4 <- FeaturePlot(HB_mHep_obj,features = DEnHBp[4],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p5 <- FeaturePlot(HB_mHep_obj,features = DEnHBp[5],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p6 <- FeaturePlot(HB_mHep_obj,features = DEnHBp[6],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p7 <- FeaturePlot(HB_mHep_obj,features = DEnHBp[7],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p8 <- FeaturePlot(HB_mHep_obj,features = DEnHBp[8],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')
p9 <- FeaturePlot(HB_mHep_obj,features = DEnHBp[9],cols  = c("lightgrey" ,"#E97D30"),pt.size = 0.05)+scale_color_gradient(trans = 'log', low = 'green', high = '#199F62')

p1+p2+p3
p4+p5+p6
p7+p8+p9



### variance analysis below:

#Get whole SCT residual matrix-----
#differentiation_all.filt.SCT <- SCTransform(differentiation_all.filt, verbose = FALSE) #get top 3000 variances genes
differentiation_all.filt.SCT <- SCTransform(
  differentiation_all.filt,
  return.only.var.genes = FALSE,
  verbose = FALSE
)
residuals_matrix <- GetAssayData(differentiation_all.filt.SCT, assay = 'SCT', slot = 'scale.data')
write.csv(as.matrix(residuals_matrix), file = "/Users/hxr/Desktop/5stage_scRNA/SCT_residuals_matrix.csv")

# Transpose it: now rows = cells, columns = genes
residuals_df <- as.data.frame(t(as.matrix(residuals_matrix)))

# Extract stage info from column names
residuals_df$cell <- rownames(residuals_df)
residuals_df$stage <- sub("-.*", "", residuals_df$cell)  # Extract stage from "H1-xxx", "DE-xxx", etc.

#plot
library(dplyr)
residual_summary <- residuals_df %>%
  rowwise() %>%
  mutate(residual_mean = mean(c_across(where(is.numeric)))) %>%
  ungroup()

ggplot(residual_summary, aes(x = stage, y = residual_mean)) +
  geom_violin(fill = "orchid", color = "black", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean Pearson Residuals per Cell by Stage",
    x = "Stage",
    y = "Mean Pearson Residual"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Split the Seurat object by 'original.ID' and check the variance separatly----
seurat_list <- SplitObject(differentiation_all.filt, split.by = "original.ID")


#Get dispersion.scaled (z-scored log dispersion) from whole matrix-----
Vst_variable <- FindVariableFeatures(differentiation_all.filt, selection.method = "vst", nfeatures = Inf, verbose = FALSE)
Vst_wholematrix <- Vst_variable@assays$RNA@meta.features[, c("vst.mean", "vst.variance", "vst.variance.expected", "vst.variance.standardized")]
Vst_wholematrix$gene <- rownames(Vst_wholematrix)
write.csv(Vst_wholematrix, "/Users/hxr/Desktop/5stage_scRNA/vst_dispersion_wholematrix.csv", row.names = FALSE)


#get SCT matrix per five stage----
#Run SCTransform per stage and extract residuals----
residuals_list <- list()  # <- initialize here
for (stage in names(seurat_list)) {
  message("Running SCTransform on stage: ", stage)
  obj <- seurat_list[[stage]]
  obj <- SCTransform(obj, verbose = FALSE)
  # Get residuals matrix (genes x cells)
  residuals <- GetAssayData(obj, assay = "SCT", slot = "scale.data")
  # Rename columns to include stage (if needed)
  colnames(residuals) <- paste0(stage, "-", colnames(residuals))
  # Store
  residuals_list[[stage]] <- residuals
}

# Combine all residuals (genes x all cells)
combined_residuals <- do.call(cbind, residuals_list)
# Save the full matrix
write.csv(as.matrix(combined_residuals), file = "/Users/hxr/Desktop/5stage_scRNA/SCT_residuals_5stages_separately.csv")

# Transpose to cells x genes
residuals_df <- as.data.frame(t(as.matrix(combined_residuals)))
residuals_df$cell <- rownames(residuals_df)
# Extract stage from cell name (e.g., "H1-cellX" → "H1")
residuals_df$stage <- sub("-.*", "", residuals_df$cell)
# Calculate per-cell residual average
library(dplyr)
residual_summary <- residuals_df %>%
  rowwise() %>%
  mutate(residual_mean = mean(c_across(where(is.numeric)))) %>%
  ungroup()

# Plot
library(ggplot2)
ggplot(residual_summary, aes(x = stage, y = residual_mean)) +
  geom_violin(fill = "lightblue", color = "black", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean Pearson Residual per Cell by Stage (SCTransform)",
    x = "Stage",
    y = "Mean Pearson Residual"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#Second methods: Get variable features from differentiation_all.filt by FindVariableFeatures per stages
# 1. Get dispersion.scaled (z-scored log dispersion) for all genes per stage----
# ------------------- H1 -------------------
H1_obj <- seurat_list[["H1"]]
DefaultAssay(H1_obj) <- "RNA"
H1_obj <- NormalizeData(H1_obj, verbose = FALSE)
H1_obj <- FindVariableFeatures(H1_obj, selection.method = "vst", nfeatures = Inf, verbose = FALSE)
vst_H1 <- H1_obj@assays$RNA@meta.features[, c("vst.mean", "vst.variance", "vst.variance.expected", "vst.variance.standardized")]
vst_H1$gene <- rownames(vst_H1)
vst_H1$stage <- "H1"

# ------------------- HB -------------------
HB_obj <- seurat_list[["HB"]]
DefaultAssay(HB_obj) <- "RNA"
HB_obj <- NormalizeData(HB_obj, verbose = FALSE)
HB_obj <- FindVariableFeatures(HB_obj, selection.method = "vst", nfeatures = Inf, verbose = FALSE)
vst_HB <- HB_obj@assays$RNA@meta.features[, c("vst.mean", "vst.variance", "vst.variance.expected", "vst.variance.standardized")]
vst_HB$gene <- rownames(vst_HB)
vst_HB$stage <- "HB"

# ------------------- DE -------------------
DE_obj <- seurat_list[["DE"]]
DefaultAssay(DE_obj) <- "RNA"
DE_obj <- NormalizeData(DE_obj, verbose = FALSE)
DE_obj <- FindVariableFeatures(DE_obj, selection.method = "vst", nfeatures = Inf, verbose = FALSE)
vst_DE <- DE_obj@assays$RNA@meta.features[, c("vst.mean", "vst.variance", "vst.variance.expected", "vst.variance.standardized")]
vst_DE$gene <- rownames(vst_DE)
vst_DE$stage <- "DE"

# ------------------- mHep -------------------
mHep_obj <- seurat_list[["mHep"]]
DefaultAssay(mHep_obj) <- "RNA"
mHep_obj <- NormalizeData(mHep_obj, verbose = FALSE)
mHep_obj <- FindVariableFeatures(mHep_obj, selection.method = "vst", nfeatures = Inf, verbose = FALSE)
vst_mHep <- mHep_obj@assays$RNA@meta.features[, c("vst.mean", "vst.variance", "vst.variance.expected", "vst.variance.standardized")]
vst_mHep$gene <- rownames(vst_mHep)
vst_mHep$stage <- "mHep"

# ------------------- iHep -------------------
iHep_obj <- seurat_list[["iHep"]]
DefaultAssay(iHep_obj) <- "RNA"
iHep_obj <- NormalizeData(iHep_obj, verbose = FALSE)
iHep_obj <- FindVariableFeatures(iHep_obj, selection.method = "vst", nfeatures = Inf, verbose = FALSE)
vst_iHep <- iHep_obj@assays$RNA@meta.features[, c("vst.mean", "vst.variance", "vst.variance.expected", "vst.variance.standardized")]
vst_iHep$gene <- rownames(vst_iHep)
vst_iHep$stage <- "iHep"

# ------------------- Combine all -------------------
vst_all_stages <- bind_rows(vst_H1, vst_DE, vst_HB, vst_iHep, vst_mHep)

#save
write.csv(vst_all_stages, "/Users/hxr/Desktop/5stage_scRNA/vst_dispersion_by_stage.csv", row.names = FALSE)

#plot
library(tidyr)
# Convert to long format for plotting multiple features
vst_long <- vst_all_stages %>%
  pivot_longer(
    cols = c(vst.mean, vst.variance, vst.variance.expected, vst.variance.standardized),
    names_to = "metric",
    values_to = "value"
  )

ggplot(vst_long, aes(x = stage, y = value)) +
  geom_violin(fill = "lightblue", color = "black", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Gene Variability Metrics Across Stages",
    x = "Stage",
    y = "Value"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

