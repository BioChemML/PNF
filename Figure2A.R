s1 <- readRDS("Kidney_IGF_PNF_annotated.rds")

s1 <- SetIdent(s1, value = "cell_type")

annotated_markers <- FindAllMarkers(
  s1, 
  only.pos = FALSE,            
  min.pct = 0.1,         
)

selected_cell_types <- c(
  "Capilary Endothelial Cell", "Connecting Tubule", "Distal Convoluted Tubule",
  "Macula Densa", "Parietal Epithelial Cell", "Podocyte", "Principal Cell",
  "Proximal Tubule", "Thick Ascending Limb", "Vascular Smooth Muscle Cell",
  "Macrophage/Immune cells", "Fibroblast_1", "Fibroblast_2", "KC_1", "KC_2"
)

top5 <- annotated_markers %>%
  filter(cluster %in% selected_cell_types) %>%
  mutate(cluster = factor(cluster, levels = selected_cell_types)) %>%
  arrange(cluster) %>% 
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()

s1_selected <- subset(s1, idents = selected_cell_types)

Idents(s1_selected) <- factor(Idents(s1_selected), levels = selected_cell_types)

s1_selected_meta <- s1_selected@meta.data %>%
  mutate(cell = rownames(.))

set.seed(42)

sampled_cells <- s1_selected_meta %>%
  group_by(cell_type = Idents(s1_selected)) %>%
  sample_n(100, replace = FALSE) %>%
  pull(cell)

s1_sampled <- subset(s1_selected, cells = sampled_cells)
s1_sampled <- ScaleData(s1_sampled, features = top5$gene)

custom_colors <- c("#CCDFE9FF", "#E9ECEBFF", "#EAE0CDFF", "#DAC9A5FF", "#CAB17EFF", 
                   "#CB8E62", "#BD7848", "#C1654F", "#B14C33", "#913E2A", 
                   "#803726", "#713122", "#6B2E20")

heatmap <- DoHeatmap(s1_sampled, features = top5$gene, slot = 'data') + 
  scale_fill_gradient(colors = custom_colors)