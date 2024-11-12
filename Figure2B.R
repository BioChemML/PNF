s1 <- readRDS("Kidney_IGF_PNF_raw.rds")

s1 <- subset(s1, subset = nFeature_RNA > 50 & nFeature_RNA < 1100)

s1 <- NormalizeData(s1, normalization.method = "LogNormalize", scale.factor = 10000)

VariableFeatures(s1) <- split(row.names(s1@meta.data), s1@meta.data$stim) %>% lapply(function(cells_use) {
  s1[,cells_use] %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    VariableFeatures()
}) %>% unlist %>% unique

s1 <- s1 %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(features = VariableFeatures(s1), npcs = 30, verbose = FALSE)

s1 <- s1 %>% 
  RunHarmony("stim", plot_convergence = FALSE, nclust = 50, theta = 3, sigma = 0.2, max_iter = 10, early_stop = F)

cluster_annotations <- c(
  "1" = "Proximal Tubule", "3" = "Proximal Tubule", "4" = "Proximal Tubule",
  "5" = "Proximal Tubule", "25" = "Proximal Tubule", "47" = "Proximal Tubule",
  "51" = "Proximal Tubule", "21" = "Proximal Tubule", "49" = "Proximal Tubule",
  "2" = "Distal Convoluted Tubule", "19" = "Distal Convoluted Tubule", 
  "22" = "Distal Convoluted Tubule", "27" = "Distal Convoluted Tubule",
  "6" = "Capilary Endothelial Cell", "9" = "Capilary Endothelial Cell",
  "18" = "Capilary Endothelial Cell", "8" = "Fibroblast_1", "11" = "Fibroblast_2",
  "15" = "Principal Cell", "42" = "Principal Cell", "48" = "Principal Cell",
  "39" = "Podocyte", "46" = "Podocyte", "16" = "Connecting Tubule", 
  "23" = "Connecting Tubule", "38" = "Connecting Tubule",
  "7" = "Thick Ascending Limb", "24" = "Thick Ascending Limb", 
  "29" = "Thick Ascending Limb",
  "33" = "KC_1", "34" = "KC_1", 
  "35" = "KC_1", "26" = "KC_1", 
  "37" = "KC_1", "50" = "KC_1", 
  "40" = "KC_1", "52" = "KC_1", 
  "13" = "KC_1", "14" = "KC_1", 
  "28" = "KC_1", "30" = "KC_1", 
  "32" = "KC_1", "43" = "KC_1", 
  "41" = "KC_1", 
  "17" = "Macrophage/Immune cells", "31" = "Vascular Smooth Muscle Cell", 
  "45" = "Parietal Epithelial Cell", "44" = "Macula Densa", 
  "20" = "KC_2", "36" = "KC_2", 
  "10" = "KC_2", "12" = "KC_2"
)

seurat_obj@meta.data$cell_type <- NA

for (cluster in names(cluster_annotations)) {
  cell_type <- cluster_annotations[[cluster]]
  s1@meta.data$cell_type[s1@meta.data$seurat_clusters == as.numeric(cluster)] <- cell_type
}

s1 <- SetIdent(s1, value = "seurat_clusters")

s1 <- RenameIdents(s1, cluster_annotations)

s1 <- FindNeighbors(s1,reduction = "harmony", dims = 1:30)

s1 <- FindClusters(s1, resolution = 2.5)

merge_s <- RunUMAP(merge_s,reduction = "harmony", min.dist = 0.1, dims = 1:20)

UMAP_s1_1 <- DimPlot(s1, reduction = "umap", label = TRUE, raster=FALSE)