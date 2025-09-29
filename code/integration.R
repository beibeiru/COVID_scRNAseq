library(Seurat)

outputPath <- "integration_oneGEO/"
dir.create(outputPath)
dir.create(paste0(outputPath,"tsne"))
dir.create(paste0(outputPath,"umap"))

ids <- c(
	"GSM4339769_C141",
	"GSM4339770_C142",
	"GSM4339772_C144",
	"GSM4339771_C143",
	"GSM4339773_C145",
	"GSM4339774_C146",
	"GSM3489182_Donor_01",
	"GSM3489185_Donor_02",
	"GSM3489187_Donor_03",
	"GSM3489189_Donor_04",
	"GSM3489191_Donor_05",
	"GSM3489193_Donor_06",
	"GSM3489195_Donor_07",
	"GSM3489197_Donor_08"
	)

batches <- c("M1","M2","M3","S1","S2","S3","HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8")

ids <- c(
	"GSM4339769_C141",
	"GSM4339770_C142",
	"GSM4339772_C144",
	"GSM4339771_C143",
	"GSM4339773_C145",
	"GSM4339774_C146",
	"GSM4475048_C51",
	"GSM4475049_C52",
	"GSM4475050_C100",
	"GSM4475051_C148",
	"GSM4475052_C149",
	"GSM4475053_C152"	
	)

batches <- c("M1","M2","M3","S1","S2","S3","H1","H2","H3","S4","S5","S6")

statMat <- data.frame()

for(i in 1:length(ids))
{
	if(grepl("Donor",ids[i]))
	{
		sc.data <- Read10X_h5(paste0("GSE122960_RAW/",ids[i],"_filtered_gene_bc_matrices_h5.h5"), use.names = TRUE, unique.features = TRUE)
	}else{
		sc.data <- Read10X_h5(paste0("GSE145926_RAW/",ids[i],"_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
	}
	colnames(sc.data) <- paste0(batches[i],"_",colnames(sc.data))
	
	sc <- CreateSeuratObject(counts = sc.data, project = batches[i])

	sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
	#head(sc@meta.data, 5)
	
	sc <- subset(sc, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & nCount_RNA >= 1000 & percent.mt <= 10)
	
	Matrix::writeMM(sc@assays$RNA@counts,file=paste0("data/",batches[i],"_counts.mtx"))
	writeLines(sc@assays$RNA@counts@Dimnames[[1]],paste0("data/",batches[i],"_genes.tsv"))
	writeLines(sc@assays$RNA@counts@Dimnames[[2]],paste0("data/",batches[i],"_cells.tsv"))
	
	statMat[batches[i],1] <- dim(sc)[1]
	statMat[batches[i],2] <- dim(sc)[2]
	
	if(!exists("sc.combined"))
	{
		sc.combined <- sc
	}else{
		sc.combined <- merge(sc.combined, y = sc, project = batches[i])
	}
	
	if(!exists("sc.list"))
	{
		sc.list <- list()
	}
	sc.list[[batches[i]]] <- sc
}

# sc.combined

sc.combined <- NormalizeData(sc.combined, normalization.method = "LogNormalize", scale.factor = 10000)
sc.combined <- FindVariableFeatures(sc.combined, selection.method = "vst", nfeatures = 2000)

sc.combined <- ScaleData(sc.combined)

sc.combined <- RunPCA(sc.combined)
sc.combined <- RunTSNE(sc.combined, dims = 1:50, check_duplicates = FALSE)
sc.combined <- RunUMAP(sc.combined, dims = 1:50)

for(reduc in c("umap","tsne"))
{
	g <- DimPlot(sc.combined, reduction = reduc, group.by = "orig.ident")
	ggplot2::ggsave(paste0(outputPath,reduc,"/uncorrected.jpg"), g, width = 18, height = 16, dpi=200, units = "cm")
}

# sc.integrated

for (i in 1:length(sc.list)) {
    sc.list[[i]] <- NormalizeData(sc.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
    sc.list[[i]] <- FindVariableFeatures(sc.list[[i]], selection.method = "vst", nfeatures = 2000)
}

sc.anchors <- FindIntegrationAnchors(object.list = sc.list, dims = 1:50)
sc.integrated <- IntegrateData(anchorset = sc.anchors, dims = 1:50)

sc.integrated <- ScaleData(sc.integrated)

sc.integrated <- RunPCA(sc.integrated)
sc.integrated <- RunTSNE(sc.integrated, dims = 1:50, check_duplicates = FALSE)
sc.integrated <- RunUMAP(sc.integrated, dims = 1:50)

for(reduc in c("umap","tsne"))
{
	g <- DimPlot(sc.integrated, reduction = reduc, group.by = "orig.ident")
	ggplot2::ggsave(paste0(outputPath,reduc,"/corrected_batch.jpg"), g, width = 18, height = 16, dpi=200, units = "cm")
}

saveRDS(sc.integrated, file = paste0(gsub("/","",outputPath),".rds"))
