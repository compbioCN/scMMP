library(scRNAtoolVis)
library(cowplot)
library(rstatix)
library(ggpubr)
library(ggbeeswarm)
library(msigdbr)
library(GSVA)
library(tidytext)
library(plyr)
library(dplyr)
library(tibble)
library(reshape2)
library(viridis)
library(ggrepel)
options(reticulate.conda_binary = "/home/Shared/miniconda3/bin/conda")
library(SCP)
library(Seurat)

#Fig2A-----------
colors <- c("#F5BC6E","#D5E7AC","#C7AED5","#F6C6CB","#8AB6D6")
dp <- cbind(epi_merge@meta.data, epi_merge@reductions$umap@cell.embeddings)
colnames(dp)[colnames(dp) %in% c("umap_1", "umap_2")] <- c("UMAP_1", "UMAP_2")
centroid <- ddply(dp, .(cluster), summarise, label_x = median(UMAP_1), label_y = median(UMAP_2))
ggplot(dp, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = colors) +
  geom_label_repel(
    data = centroid, 
    aes(x = label_x, y = label_y, label = cluster),
    color = "black",       
    fill = NA,             
    label.size = NA,       
    size = 4,
    show.legend = FALSE    
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  ggtitle("") +
  theme(
    plot.title = element_blank(),
    text = element_text(family = "Helvetica",size = 14),
    axis.title.x = element_text(family = "Helvetica",size = 16),
    axis.title.y = element_text(family = "Helvetica",size = 16),
    legend.position = "none"
  )
ggsave(filename = "Fig/Fig2A.pdf",width = 5,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica")

genes <- c("COL17A1","ASS1","CCL2","KRT1","TOP2A","TGM3")
for(gene in genes){
  umap_data <- epi_merge@reductions$umap@cell.embeddings
  gene_expression <- FetchData(epi_merge, vars = c(gene))
  plot_data <- cbind(umap_data, gene_expression)
  colnames(plot_data) <- c("UMAP_1", "UMAP_2", "gene")
  plot_data$gene <- ifelse(plot_data$gene > 2, 2, plot_data$gene)
  ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = gene)) +
    geom_point(size = 0.3)+
    scale_color_gradientn(colors = c('grey90', 'Violet', 'Dark Violet'), limits = c(0, 2) )+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    theme(
      text = element_text(family = "Helvetica",size = 14),
      axis.title.x = element_text(family = "Helvetica",size = 16),
      axis.title.y = element_text(family = "Helvetica",size = 16),
      legend.position = "none"
    )+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.title = element_blank(),  
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'transparent'), 
          plot.background=element_rect(fill="transparent"),color = NA)+ 
    ggtitle(gene)+
    theme(plot.title = element_text(size = 60, face = "bold", hjust = 0.5,family = "Helvetica"))
  ggsave(filename = paste0("Fig/Fig2A_",gene,".pdf"),width = 5,height = 5.9,dpi = 600,device = cairo_pdf, family = "Helvetica")
}
#Fig2B-------------
markers <- c(
  "IFGBP3","ASS1","COL17A1","SLC2A1","MYC","CXCL14","KRT19",
  "KRT1","KRT10","SPRR1B","CLDN4",
  "TOP2A","MKI67",
  "FLG","IVL","TGM3",
  "CCL2","IRF1"
)
markers <- markers[18:1]
DotPlot(epi_merge,features = markers,scale = TRUE)+RotatedAxis()+
  scale_color_gradientn(colours = colorRampPalette(c("#4080B4", "#FCFCB2", "#E52F21"))(100))+
  theme(
    axis.text.x = element_text(family = "Helvetica"),
    axis.text.y = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  )+
  coord_flip()
ggsave(filename = "Fig/Fig2B.pdf",width = 4.5,height = 7.5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig2C-------------
Idents(epi_merge) <- "cluster" 
expr <- AverageExpression(epi_merge, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] 
expr <- as.matrix(expr)
x <- msigdbr_collections()
human_GO = msigdbr(species = "Homo sapiens", #物种
                   category = "C5",
                   subcategory = "BP") %>% 
  dplyr::select(gs_name,gene_symbol)#选择gene symbol或者ID
human_KEGG_Set = human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
human_GO_Set = human_GO %>% split(x = .$gene_symbol, f = .$gs_name)
GENESET_LIST <- GENESET_LIST %>% split(x = .$gene_symbol, f = .$gs_name)
params <- GSVA::gsvaParam(exprData = expr,
                          geneSets = GENESET_LIST,
                          kcdf = "Gaussian")
go_cluster <- GSVA::gsva(params)
rownames(go_cluster) <- go_cluster$...1
go_cluster <- go_cluster[,-1]
go_cluster <- go_cluster[,c(1,3,4,2,5)]
go_cluster_long <- melt(as.matrix(go_cluster), varnames = c("ID", "Cluster"), value.name = "gsva_Score")
ggplot(go_cluster_long, aes(x = Cluster, y = ID, fill = gsva_Score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(colors = colorRampPalette(c("navy", "white", "firebrick3"))(50)) +
  labs(x = "", y = "", fill = "Score") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(
      family = "Helvetica",
      size = 14,
      margin = margin(l = 10, r = 0),  
      hjust = 0  
    ),
    axis.ticks.y = element_line(),  
    axis.ticks.length.y = unit(0.2, "cm"),  
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      family = "Helvetica",
      size = 14
    ),
    panel.grid = element_blank(),
    legend.text = element_text(family = "Helvetica", size = 14),
    legend.title = element_text(family = "Helvetica", size = 14),
    axis.text.y.right = element_text(
      family = "Helvetica",
      size = 14,
      hjust = 0
    ),
    axis.line.y.right = element_line(), 
    axis.ticks.y.right = element_line()  
  ) +
  scale_y_discrete(position = "right")
ggsave(filename = "Fig//Fig2C.pdf",width = 9.9,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig2D-------------
colors <- c("#726BAE",  "#E5948E","#60A897", "#EA945A", "#78C4D4")
dp <- cbind(epi_merge@meta.data, epi_merge@reductions$umap@cell.embeddings)
colnames(dp)[colnames(dp) %in% c("umap_1", "umap_2")] <- c("UMAP_1", "UMAP_2")
centroid <- ddply(dp, .(Sub_celltype), summarise, label_x = median(UMAP_1), label_y = median(UMAP_2))
ggplot(dp, aes(x = UMAP_1, y = UMAP_2, color = Sub_celltype)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = colors) +
  geom_label_repel(
    data = centroid, 
    aes(x = label_x, y = label_y, label = Sub_celltype),
    color = "black",       
    fill = NA,            
    label.size = NA,       
    size = 4,
    show.legend = FALSE    
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  ggtitle("") +
  theme(
    plot.title = element_blank(),
    text = element_text(family = "Helvetica",size = 14),
    axis.title.x = element_text(family = "Helvetica",size = 16),
    axis.title.y = element_text(family = "Helvetica",size = 16),
    legend.position = "none"
  )
ggsave(filename = "Fig/Fig2D.pdf",width = 5,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig2E-----------
colors <- c("#726BAE",  "#E5948E","#60A897", "#EA945A", "#78C4D4")
sample.prop <- as.data.frame(prop.table(table(epi_merge$Sub_celltype,epi_merge$Patient_merge)))
colnames(sample.prop) <- c("celltype","sample","proportion")
ggplot(sample.prop,aes(x=sample,y=proportion,fill=celltype))+
  geom_bar(stat = "identity",position = "fill")+
  theme_bw()+
  ylab("Proportion")+
  xlab("")+
  guides(fill=guide_legend(title = NULL))+
  scale_fill_manual(values = colors)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15,family = "Helvetica"),  
    axis.text.y = element_text(size = 15,family = "Helvetica"),  
    axis.title.x = element_text(size = 15,family = "Helvetica"),  
    axis.title.y = element_text(size = 15,family = "Helvetica"),  
    axis.ticks.length = unit(0.3, "cm"),
    plot.title = element_text(hjust = 0.5, size = 15,family = "Helvetica"),  
    legend.title = element_text(size = 15,family = "Helvetica"),  
    legend.text = element_text(size = 15,family = "Helvetica") 
  )
ggsave(filename = "Fig/Fig2E.pdf",width = 4.5,height = 5.5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig2F----------
gene_df <- FetchData(
  object = basal,
  vars = c("WNT4")
)
metadata_column <- basal@meta.data[, c("Sub_celltype","Patient_merge","sample_merge"), drop = FALSE]
combined_data <- cbind(metadata_column,gene_df)
colnames(combined_data) <- c("Sub_celltype", "Patient_merge","sample_merge", "gene")
stat_test <- combined_data %>%
  wilcox_test(gene ~ Patient_merge) %>%
  adjust_pvalue(method = "BH") %>%  
  add_significance("p.adj") %>%
  add_xy_position(x = "Patient_merge", dodge = 0.8)  

ggviolin(
  combined_data,
  x = "Patient_merge",
  y = "gene",
  color = "Patient_merge",
  palette = c("#E34256", "#4A89B0", "#CBC3E3"),
  add = c("boxplot"),
  add.params = list(alpha = 0.8,  
                    size = 0.5,   
                    jitter = 0.2  
  ),shape = 21
) +
  stat_compare_means(
    comparisons = list(
      c("MMPTissue", "PVTissue"),  
      c("PVTissue", "HealthyTissue"),  
      c("MMPTissue", "HealthyTissue")  
    ),
    method = "wilcox.test", 
    label = "p.signif",      
    step.increase = 0.1,     
    tip.length = 0.01,family = "Helvetica",size =5
  ) +theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
    axis.line = element_line(colour = "black", linewidth = 0.5),            
    panel.grid.major = element_blank(),                                     
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(family = "Helvetica",angle = 45,hjust = 1,size = 15),
    axis.text.y = element_text(family = "Helvetica",size = 15),
    axis.title.y = element_text(family = "Helvetica",size = 15),
    legend.text = element_text(family = "Helvetica",size = 15),
    legend.title = element_text(family = "Helvetica",size = 15),
    text = element_text(family = "Helvetica")
  )+ 
  xlab("")+
  ylab("The expression levels of WNT4")+
  ylim(-0.25,3.1)+
  labs(color = "Group")
ggsave(filename = "Fig/Fig2F.pdf",width = 5,height = 5.5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig2G-----------
ego_all$Description <- factor(ego_all$Description, levels = ego_all$Description[order(ego_all$FoldEnrichment, decreasing = FALSE)])
ego_all <- ego_all[order(ego_all$Count,decreasing = T),]
ego_all$ID <- rep(1:ceiling(nrow(ego_all) / 3), each = 3, length.out = nrow(ego_all))
custom_labels <- c("Basal" = "Basal","Granular"="Granular","Spinous"="Spinous", "Proliferating" = "Proliferating", "CCL2+_Epi" = "CCL2+_Epi" )
ego_all_sorted <- ego_all[order(ego_all$Count,decreasing = T),]
ego_all_sorted <- ego_all %>%
  group_by(ID) %>%              
  arrange(desc(Count), .by_group = TRUE) %>%  
  ungroup()    
ego_all <- ego_all %>%
  mutate(Description_ordered = fct_reorder2(Description, ID, -Count))
ego_all$ID <- factor(ego_all$ID,levels = celltype_order)
ggplot(ego_all_sorted, aes(x = reorder_within(Description, Count, ID), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_gradient(low = "#276D9F", high = "#67C9F2" ) + 
  theme_minimal() +
  scale_x_reordered() +
  labs(x = "GO Term", y = "Gene Count", title = "GO Enrichment Analysis") +
  facet_wrap(~  ID, scales = "free_y",ncol = 1,labeller = labeller(ID = custom_labels),strip.position = "right") +  
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold",family = "Helvetica"),
    axis.text.x = element_text(size = 14, face = "bold",family = "Helvetica"),
    axis.title = element_text(size = 14, face = "bold",family = "Helvetica"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5,family = "Helvetica"),
    legend.title = element_text(size = 14,family = "Helvetica"),
    legend.text = element_text(size = 14,family = "Helvetica"),
    strip.text = element_text(size = 14, face = "bold", hjust = 0.5,family = "Helvetica"),
    strip.background = element_rect(color = "black", size = 1)
  )
ggsave(filename = "Fig/Fig2G.pdf",width = 10,height = 8,dpi = 600,device = cairo_pdf, family = "Helvetica")






# 覆盖上皮细胞的detailed_cell_type
common_cells <- colnames(cells_subset_filtered1)[colnames(cells_subset_filtered1) %in% colnames(seurat_metadata)]
epi_indices <- match(common_cells, colnames(cells_subset_filtered1))
seurat_indices <- match(common_cells, colnames(seurat_metadata))
seurat_metadata$detailed_subtype[seurat_indices] <- cells_subset_filtered1$detailed_cell_type[epi_indices]


custom_palette = {
  "Adipocytes":"#023FA5",
  "B_cells":"#7D87B9",
  "CCL18+CD163+_TAMs":"#BEC1D4",
  "CLDN2+_Epi": "#D6BCC0",
  "COL17A1+_Epi": "#BB7784",
  "Endothelial":"#8E063B",
  "FOSB+FOS+_Epi": "#4A6FE3",
  "LOR+_Epi": "#8595E1",
  "RGS1+_TAMs":"#A5D8FF",
  "MMP1+MMP10+_CAFs":"#B5BBE3",
  #"MMP1+_CAFs":"#B5BBE3",
  "MUC4+Epi": "#E6AFB9",
  "Myocytes":"#E07B91",
  "PLA2G2A+_CAFs":"#D33F6A",
  #"PLA2G2A+_CAFs":"#D33F6A",    
  "Plasma":"#11C638",
  "SERPINB4+_Epi": "#8DD593",
  "CXCL8+CAFs":"#C6DEC7",
  #"SFRP4+_CAFs":"#C6DEC7",
  "SOX2+_Epi": "#EAD3C6",
  "SOX14+_Epi": "#F0B98D",    
  "SPP1+_TAMs":"#EF9708",
  "ACTA2+_CAFs":"#9CDED6",
  #"TDO2+_CAFs":"#9CDED6",
  "SERPINB12+_CAFs":"#FFD166",
  "TGM3+_Epi": "#0FCFC0",
  "T_cells":"#D5EAE7",
  "VSMCs":"#F3E1EB",
  "others":"#F8C4E1" 
}

TLS_palette = {
  "B cells":"#60A897",
  "Plasmacytoid dendritic cells":"#023FA5",
  "T cells":"#E05F48",
  "MRC1⁺/CD163⁺ macrophages":"#FFD166"
}

