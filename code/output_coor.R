library(Seurat)


outputPath <- "integration_oneGEO/"

sc.integrated <- readRDS(file = "integration_oneGEO.rds")

write.table(sc.integrated@reductions$umap@cell.embeddings,paste0(outputPath,"cell_coordinates_umap.tsv"),quote=F,sep="\t")
write.table(sc.integrated@reductions$tsne@cell.embeddings,paste0(outputPath,"cell_coordinates_tsne.tsv"),quote=F,sep="\t")


outputPath <- "integration_twoGEO/"

sc.integrated <- readRDS(file = "integration_twoGEO.rds")

write.table(sc.integrated@reductions$umap@cell.embeddings,paste0(outputPath,"cell_coordinates_umap.tsv"),quote=F,sep="\t")
write.table(sc.integrated@reductions$tsne@cell.embeddings,paste0(outputPath,"cell_coordinates_tsne.tsv"),quote=F,sep="\t")
