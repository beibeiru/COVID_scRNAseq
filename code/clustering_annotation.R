library(Seurat)

outputPath <- "integration_oneGEO/"
dir.create(outputPath)
dir.create(paste0(outputPath,"marker"))


sc.integrated <- readRDS(file = "integration_oneGEO.rds")

sc.integrated <- FindNeighbors(sc.integrated, dims = 1:50)
sc.integrated <- FindClusters(sc.integrated, resolution = 1)

clusterRes <- data.frame(
	cellID=names(sc.integrated @ active.ident),
	clusterID=sc.integrated @ active.ident,
	stringsAsFactors=FALSE
	)
write.table(clusterRes,paste0(outputPath,"cluster_results.tsv"),quote=F,sep="\t",row.names=F)


for(reduc in c("umap","tsne"))
{
	g <- DimPlot(sc.integrated, reduction = reduc, group.by = "seurat_clusters")
	ggplot2::ggsave(paste0(outputPath,reduc,"/corrected_cluster.jpg"), g, width = 18, height = 16, dpi=200, units = "cm")
}


DefaultAssay(sc.integrated) <- "RNA"

markers.list <- list(
	level1=c("AGER","SFTPC","SCGB3A2","TPPP3","KRT5","CD68","FCN1","CD1C","TPSB2","CD3D","IGHG4","MS4A1","VWF","DCN"),
	level2.1=c("FCN1","SPP1","FABP4"),
	level2.2=c("KLRF1","KLRC1","KLRD1","NKG7","CD8A","CD8B","CD4","CCR7","IL7R","CTLA4","FOXP3","IL2RA","TYMS","MKI67")
	)

for(i in 1:length(markers.list))
{
	g <- FeaturePlot(sc.integrated, markers.list[[i]], reduction = "tsne",ncol=5)
	ggplot2::ggsave(paste0(outputPath,"marker/",names(markers.list)[i],"_feature_tsne.jpg"), g, width = 75, height = 45, dpi=300, units = "cm")
	
	g <- FeaturePlot(sc.integrated, markers.list[[i]], reduction = "umap",ncol=5)
	ggplot2::ggsave(paste0(outputPath,"marker/",names(markers.list)[i],"_feature_umap.jpg"), g, width = 75, height = 45, dpi=300, units = "cm")
	
	g <- VlnPlot(sc.integrated, markers.list[[i]], ncol=5)
	ggplot2::ggsave(paste0(outputPath,"marker/",names(markers.list)[i],"_vln.jpg"), g, width = 75, height = 45, dpi=300, units = "cm")

	g <- DotPlot(sc.integrated, features=rev(markers.list[[i]]), cols = c("white", "red"), dot.scale = 8)
	ggplot2::ggsave(paste0(outputPath,"marker/",names(markers.list)[i],"_dot.jpg"), g, width = 35, height = 40, dpi=300, units = "cm")
}

if(outputPath == "integration_oneGEO/"){
cluster_annotation_vector <- c(
"0","Myeloid","Macrophage",
"1","Myeloid","Monocyte",
"2","Myeloid","Macrophage",
"3","Myeloid","Macrophage",
"4","Myeloid","Macrophage",
"5","Myeloid","Macrophage",
"6","T","T",
"7","Myeloid","Macrophage",
"8","Myeloid","Macrophage",
"9","Myeloid","Macrophage",
"10","Myeloid","Macrophage",
"11","Myeloid","Macrophage",
"12","Myeloid","Monocyte",
"13","T","T",
"14","T","T",
"15","Myeloid","Macrophage",
"16","Epithelial","Epithelial",
"17","T","T",
"18","Myeloid","Macrophage",
"19","B","Plasma",
"20","Myeloid","Macrophage",
"21","Epithelial","Epithelial",
"22","Myeloid","Dendritic",
"23","B","B",
"24","Myeloid","Macrophage",
"25","Myeloid","Mast"
)
}else{
cluster_annotation_vector <- c(
"0","Myeloid","Macrophage",
"1","Epithelial","Epithelial",
"2","Myeloid","Macrophage",
"3","Myeloid","Macrophage",
"4","Epithelial","Epithelial",
"5","Myeloid","Macrophage",
"6","Myeloid","Monocyte",
"7","Epithelial","Epithelial",
"8","Myeloid","Macrophage",
"9","Myeloid","Macrophage",
"10","T","T",
"11","Myeloid","Macrophage",
"12","T","T",
"13","Mensenchymal","Endothelial",
"14","Epithelial","Epithelial",
"15","Myeloid","Macrophage",
"16","T","T",
"17","Epithelial","Epithelial",
"18","Epithelial","Epithelial",
"19","B","Plasma",
"20","Mensenchymal","Fibroblast",
"21","Myeloid","Dendritic",
"22","Epithelial","Epithelial",
"23","Myeloid","Macrophage",
"24","Epithelial","Epithelial",
"25","Myeloid","Mast",
"26","Epithelial","Epithelial",
"27","B","B",
"28","Epithelial","Epithelial"
}

cluster_annotation <- matrix(
	cluster_annotation_vector,
	ncol=3,
	byrow=TRUE)
colnames(cluster_annotation) <- c("clusterID","cellType","cellTypeSub")

write.table(cluster_annotation,paste0(outputPath,"cluster_annotation.txt"),quote=F,sep="\t",row.names=F)







cluster_results <- as.matrix(read.table(paste0(outputPath,"cluster_results.tsv"),as.is=T,sep="\t",header=T))

Tclusters <- cluster_annotation[cluster_annotation[,"cellType"]=="T","clusterID"]
Tcells <- cluster_results[cluster_results[,"clusterID"]%in%Tclusters,"cellID"]

DefaultAssay(sc.integrated) <- "integrated"
sc.integrated.T <- sc.integrated[,Tcells]

sc.integrated.T <- RunPCA(sc.integrated.T)
sc.integrated.T <- RunTSNE(sc.integrated.T, dims = 1:50, check_duplicates = FALSE)
sc.integrated.T <- RunUMAP(sc.integrated.T, dims = 1:50)

sc.integrated.T <- FindNeighbors(sc.integrated.T, dims = 1:50)
sc.integrated.T <- FindClusters(sc.integrated.T,resolution=1)

clusterRes <- data.frame(
	cellID=names(sc.integrated.T @ active.ident),
	clusterID=sc.integrated.T @ active.ident,
	stringsAsFactors=FALSE
	)
write.table(clusterRes,paste0(outputPath,"cluster_T_results.tsv"),quote=F,sep="\t",row.names=F)


for(reduc in c("umap","tsne"))
{
	g <- DimPlot(sc.integrated.T, reduction = reduc, group.by = "seurat_clusters",label=TRUE)
	ggplot2::ggsave(paste0(outputPath,reduc,"/corrected_cluster_T.jpg"), g, width = 18, height = 16, dpi=200, units = "cm")
}




DefaultAssay(sc.integrated.T) <- "RNA"

markers.list <- list(
	level2.2sub=c("KLRF1","KLRC1","KLRD1","NKG7","CD8A","CD8B","CD4","CCR7","IL7R","CTLA4","FOXP3","IL2RA","TYMS","MKI67")
	)

for(i in 1:length(markers.list))
{
	g <- FeaturePlot(sc.integrated.T, markers.list[[i]], reduction = "tsne",ncol=5)
	ggplot2::ggsave(paste0(outputPath,"marker/",names(markers.list)[i],"_feature_tsne.jpg"), g, width = 75, height = 45, dpi=300, units = "cm")
	
	g <- FeaturePlot(sc.integrated.T, markers.list[[i]], reduction = "umap",ncol=5)
	ggplot2::ggsave(paste0(outputPath,"marker/",names(markers.list)[i],"_feature_umap.jpg"), g, width = 75, height = 45, dpi=300, units = "cm")
	
	g <- VlnPlot(sc.integrated.T, markers.list[[i]], ncol=5)
	ggplot2::ggsave(paste0(outputPath,"marker/",names(markers.list)[i],"_vln.jpg"), g, width = 75, height = 45, dpi=300, units = "cm")

	g <- DotPlot(sc.integrated.T, features=rev(markers.list[[i]]), cols = c("white", "red"), dot.scale = 8)
	ggplot2::ggsave(paste0(outputPath,"marker/",names(markers.list)[i],"_dot.jpg"), g, width = 35, height = 40, dpi=300, units = "cm")
}

if(outputPath == "integration_oneGEO/"){
cluster_annotation_T_vector <- c(
"0","CD8",
"1","NK",
"2","CD4",
"3","CD8",
"4","NK",
"5","CD4",
"6","Treg",
"7","Proliferating",
"8","Proliferating",
"9","CD8",
"10","CD8",
"11","CD8"
)
}else{
cluster_annotation_T_vector <- c(
"0","CD4",
"1","CD8",
"2","Proliferating",
"3","NK",
"4","CD8",
"5","CD8",
"6","NK",
"7","Proliferating",
"8","Proliferating",
"9","Treg",
"10","Proliferating"
)
}

cluster_annotation_T <- matrix(
	cluster_annotation_T_vector,
	ncol=2,
	byrow=TRUE)
colnames(cluster_annotation_T) <- c("TclusterID","cellTypeSub")

write.table(cluster_annotation_T,paste0(outputPath,"cluster_T_annotation.txt"),quote=F,sep="\t",row.names=F)

