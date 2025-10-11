library(Seurat)
library(cowplot)
library(rstatix)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(ggbeeswarm)
library(CellChat)
library(tidytext)
library(plyr)
library(dplyr)
library(tibble)
library(reshape2)
library(viridis)
library(ggrepel)
options(reticulate.conda_binary = "/home/Shared/miniconda3/bin/conda")
library(SCP)

#Fig3A--------
colors <- c("#9ECAE1","#C7AED5")
MMP$sub_celltype <- factor(MMP$sub_celltype,levels = c("CD4+ Tnaive","CD4+ Trm","CD4+ Tem","Treg","Th17","CD8+ Tnaive","CD8+ Tem","CD8+ Tex","MAIT","NKT","NK"))
Idents(MMP) <- MMP$sub_celltype
dp <- cbind(MMP@meta.data, MMP@reductions$umap@cell.embeddings)
colnames(dp)[colnames(dp) %in% c("umap_1", "umap_2")] <- c("UMAP_1", "UMAP_2")
centroid <- ddply(dp, .(Sample), summarise, label_x = median(UMAP_1), label_y = median(UMAP_2))
axis_length <- 2  
ggplot(dp, aes(x = UMAP_1, y = UMAP_2, color = Sample)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = colors) +
  geom_label_repel(
    data = centroid, 
    aes(x = label_x, y = label_y, label = Sample),
    color = "black",       
    fill = NA,             
    label.size = NA,       
    size = 4.5,
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
  theme(
    plot.title = element_blank(),
    text = element_text(family = "Helvetica",size = 16),
    axis.title.x = element_text(family = "Helvetica",size = 16),
    axis.title.y = element_text(family = "Helvetica",size = 16),
    legend.position = "none"
  )+
  ggtitle("")+
  theme(plot.title = element_blank(),
        text = element_text(family = "Helvetica"),
        legend.position = "none")
ggsave(filename = "Fig/Fig3A",width = 5.3,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig3B--------
colors <- c("#C7AED5","#A5AA99","#60A897","#EA945A","#ACD48A","#F7DBF0","#E5948E","#276D9F","#8AB6D6","#86A667","#E7C5DB")
MMP$sub_celltype <- factor(MMP$sub_celltype,levels = c("CD4+ Tnaive","CD4+ Trm","CD4+ Tem","Treg","Th17","CD8+ Tnaive","CD8+ Tem","CD8+ Tex","MAIT","NKT","NK"))
dp <- cbind(MMP@meta.data, MMP@reductions$umap@cell.embeddings)
colnames(dp)[colnames(dp) %in% c("umap_1", "umap_2")] <- c("UMAP_1", "UMAP_2")
centroid <- ddply(dp, .(sub_celltype), summarise, label_x = median(UMAP_1), label_y = median(UMAP_2))
axis_length <- 2  
ggplot(dp, aes(x = UMAP_1, y = UMAP_2, color = sub_celltype)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = colors) +
  geom_label_repel(
    data = centroid, 
    aes(x = label_x, y = label_y, label = sub_celltype),
    color = "black",        
    fill = NA,             
    label.size = NA,       
    size = 4.5,
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
  theme(
    plot.title = element_blank(),
    text = element_text(family = "Helvetica",size = 16),
    axis.title.x = element_text(family = "Helvetica",size = 16),
    axis.title.y = element_text(family = "Helvetica",size = 16),
    legend.position = "none"
  )+
  ggtitle("")+
  theme(plot.title = element_blank(),
        text = element_text(family = "Helvetica"),
        legend.position = "none")
ggsave(filename = "Fig/Fig3B.pdf",width = 5.5,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig3C--------
colors <- c("#C7AED5","#A5AA99","#60A897","#EA945A","#ACD48A","#F7DBF0","#E5948E","#276D9F","#8AB6D6","#86A667","#E7C5DB")
sample.prop <- as.data.frame(prop.table(table(MMP$sub_celltype,MMP$Sample)))
colnames(sample.prop) <- c("celltype","sample","proportion")
sample.prop$sample <- factor(sample.prop$sample,levels = c("MMPTissue","HealthyTissue"))
ggplot(sample.prop,aes(x=sample,y=proportion,fill=celltype))+
  geom_bar(stat = "identity",position = "fill")+
  theme_bw()+
  ylab("Proportion")+
  xlab("")+
  guides(fill=guide_legend(title = NULL))+
  scale_fill_manual(values = colors)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15,family = "Helvetica"),  # x轴标签字体大小
    axis.text.y = element_text(size = 15,family = "Helvetica"),  # y轴标签字体大小
    axis.title.x = element_text(size = 15,family = "Helvetica"),  # x轴标题字体大小
    axis.title.y = element_text(size = 15,family = "Helvetica"),  # y轴标题字体大小
    axis.ticks.length = unit(0.3, "cm"),
    plot.title = element_text(hjust = 0.5, size = 15,family = "Helvetica"),  # 标题字体大小
    legend.title = element_text(size = 15,family = "Helvetica"),  # 图例标题字体大小
    legend.text = element_text(size = 15,family = "Helvetica")  # 图例文本字体大小
  )
ggsave(filename = "Fig/Fig3C.pdf",width = 4,height = 4.8,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig3D--------
markers <- c(
  "CD3E","CD4","CCR7","LEF1",
  "CD69","PTGER4","IL7R","CXCR6","MYADM",
  "ANXA1","ANXA2","RGS1","CD44",
  "FOXP3","IL2RA","IKZF2",
  "RORA","RORC","IL17A",
  "GZMK","CD8A","CD8B","KLF2",
  "CXCR3",
  "ITGAE",
  "BHLHE40","S1PR1","PRF1",
  "ZBTB16","NCR3","TRDC",
  "GNLY","NKG7","KLRC1"
)
markers <- markers[35:1]
DotPlot(MMP,features = markers,scale = TRUE)+RotatedAxis()+
  scale_color_gradientn(colours = colorRampPalette(c("#4080B4", "#FCFCB2", "#E52F21"))(100))+
  coord_flip(clip = "off")+
  theme(
    axis.text.x = element_text(family = "Helvetica"),
    axis.text.y = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica"),
    axis.text.y.right = element_text(family = "Helvetica"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave(filename = "Fig/Fig3D.pdf",width = 5.3,height = 9.5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig3E-----------

anti <- read.csv("~/MMP/analysis/7.0.0_MMP/li_SOUPX/MsigDB/GO0050900.csv",header = F)
gene <- as.data.frame(anti[,2])
MMPvsHea <- AddModuleScore(MMP,features = gene, name = 'acti_score') 
actiscore <- MMPvsHea@meta.data[,c("Sample","sub_celltype","acti_score1")]
actiscore$sub_celltype <-as.factor(actiscore$sub_celltype)
ggboxplot(actiscore, 
          x = "Sample", 
          y = "acti_score1", 
          fill = "Sample",  
          color = "Sample", 
          facet.by = "sub_celltype", 
          xlab = "", 
          ylab = "Leukocyte Migration", 
          ylim = c(-0.1, 0.4), 
          ncol = 6,
          palette = "npc",  
          outlier.size = 0.8,      
          size = 0.5 ,               
          outlier.color = "black" 
) +
  scale_fill_manual(values = rep("white", length(unique(actiscore$Sample)))) + # 填充设为白色
  scale_color_manual(values = colors) + 
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1,family = "Helvetica",size = 12), 
    axis.text.y = element_text(family = "Helvetica",size = 12),
    axis.title.y = element_text(family = "Helvetica",size = 14),
    strip.text = element_text(
      family = "Helvetica",  
      size = 14            
    ),
    legend.position = "none"
  ) +

  geom_signif(
    comparisons = list(c("MMPTissue", "HealthyTissue")), 
    map_signif_level = TRUE, 
    textsize = 4, 
    size = 0.5,family = "Helvetica"
  )
ggsave(filename = "Fig/Fig3E.pdf",width = 10,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig3F-------
TH17_MMP <- subset(MMP,sub_celltype == "Th17")
gene_df <- FetchData(
  object = TH17_MMP,
  vars = c("PDCD1")
)
TH17_MMP[["RNA3"]] <- as(object = TH17_MMP[["RNA"]], Class = "Assay")  
DefaultAssay(TH17_MMP) <- 'RNA3'
TH17_MMP[['RNA']] <- NULL
TH17_MMP <- RenameAssays(TH17_MMP, RNA3  = 'RNA')

SCP::FeatureStatPlot( srt = TH17_MMP,
                      group.by="Sample",trend_ptsize = 20,pt.size = 20,sig_labelsize = 7,
                      stat.by= c("PDCD1"), stat_size = 15, add_box = TRUE, comparisons = list( c("MMPTissue", "HealthyTissue")),
                      bg.by='Sample')+
  xlab("")+
  ylab("The expression levels of PDCD1")+
  ggtitle("")+
  theme(
    axis.text.x = element_text(family = "Helvetica",size = 14),
    axis.text.y = element_text(family = "Helvetica",size = 14),
    axis.title.y = element_text(family = "Helvetica",size = 14),
    legend.text = element_text(family = "Helvetica",size = 14),
    legend.title = element_text(family = "Helvetica",size = 15),
    legend.position = "none",
    plot.title.position = "plot",  # 确保标题不占用布局空间
    plot.title = element_text(size = 0)
  )
ggsave(filename = "Fig/Fig3F.pdf",width = 2.5,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig3G-------
cellchat@idents <- factor(cellchat@idents,levels = c(levels(MMP)[11:1],"Basal"))
source_use <- "Basal"
target_use <- levels(cellchat@idents)[-12]
netVisual_bubble(cellchat, sources.use = source_use, targets.use = target_use, signaling = c("MIF","CD99","CLEC","LAMININ","APP","CXCL"),comparison = c(1, 2),color.text = c("skyblue","pink"),
                 max.dataset = 2, title.name = "Increased signaling in TIL", angle.x = 45, remove.isolate = T,font.size = 14,font.size.title = 16)+
  scale_colour_gradientn(colors = colorRampPalette(c("#8AB6D6", "white", "#E05F48"))(100))+
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 12.5,family = "Helvetica"),
    legend.text = element_text(size = 14,family = "Helvetica"),   # 进一步确保图例文本大小
    legend.title = element_text(size = 14,family = "Helvetica")    # 进一步确保图例标题大小
  )
ggsave(filename = "./Fig/Fig3G.pdf",width = 11.5,height = 5.8,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig3H-------
avg_tf_expr <- tf_data %>%
  group_by(Sample, sub_celltype) %>%
  summarise_at(vars(everything()), mean, na.rm = TRUE)
avg_tf_expr_long <- avg_tf_expr %>%
  pivot_longer(cols = -c(Sample, sub_celltype), names_to = "TF", values_to = "Mean_Expression")
avg_tf_expr_wide <- avg_tf_expr_long %>%
  unite(Sample_sub_celltype, Sample, sub_celltype, sep = "_") %>%  
  pivot_wider(names_from = Sample_sub_celltype, values_from = Mean_Expression)
heatmap_data <- as.matrix(avg_tf_expr_wide[, -1])  
rownames(heatmap_data) <- avg_tf_expr_wide$TF
heatmap_data <- as.data.frame(heatmap_data)
colnames(heatmap_data) <- gsub("Tissue","",colnames(heatmap_data))
heatmap_data <- heatmap_data[,c(1,12,2,13,3,14,4,15,5,16,6,17,7,18,8,19,9,20,10,21,11,22)]
heatmap_data <- heatmap_data[,c(12:22,1:11)]
annotation_col <- data.frame(Sample = gsub("_.*", "", colnames(heatmap_data)),
                             Celltype = gsub(".*_", "", colnames(heatmap_data)))
rownames(annotation_col) <- colnames(heatmap_data)
gene_ids_T <- c("FOSL2","PKNOX1","SP3","RARA","BCL11B") 
gene_ids_Th17 <- c("STAT5A","NFATC2","NFKB1","SMAD4","RXRB")
gene_ids_leukocyte_activation <- c("MYB","XBP1","STAT6","JUNB","TBX21")
gene_ids_myeloid_cell_differentiation <- c("KLF10","ZBTB7A","CEBPA","ETS1","HOXA9")
TF_groups <- data.frame(
  Group =c(rep("T cell differentiation",5),rep("Th17 cell differentiation",5),rep("Leukocyte activation involved in immune response",5),rep("Myeloid cell differentiation",5)),
  TF = c(gene_ids_T,gene_ids_Th17,gene_ids_leukocyte_activation,gene_ids_myeloid_cell_differentiation)
)
rownames(TF_groups) <- paste0(TF_groups$TF,"(+)")
sample_groups <- TF_groups$Group
custom_dist_row <- as.dist(outer(sample_groups, sample_groups, `!=`)) 
cluster_rows <- cutree(p$tree_row, k = 4)
block_anno <- data.frame(
  Block = TF_groups$Group
)
pdf(file = "./Fig/heatmap_TF_T_Cells.pdf",width = 7,height = 6)
pheatmap::pheatmap(heatmap_data[paste0(TF_groups$TF,"(+)"),],angle_col = 90,
                   scale = "row",  
                   cluster_rows = T,treeheight_row = 0,breaks = my_breaks <- seq(-2, 2, length.out = 100) ,
                   cluster_cols = T,treeheight_col = 0,
                   border_color = "gray60",cutree_cols = 2,cutree_rows = 4,
                   color = colorRampPalette(c("#8AB6D6", "white", "#E05F48"))(100),
                   show_rownames = TRUE, labels_row = c(paste0(TF_groups$TF,"(+)")) ,annotation_names_row = TRUE, 
                   show_colnames = TRUE, 
                   annotation_col = annotation_col,
                   clustering_distance_cols = custom_dist_col,
                   clustering_distance_rows = custom_dist_row,
                   annotation_colors = list(Sample = c("Healthy" = "#726BAE", "MMP" = "#F5BC6E"),
                                            Celltype = c("CD4+ Tnaive" = "#C7AED5","CD4+ Trm" = "#A5AA99","CD4+ Tem" = "#60A897",
                                                         "Treg" ="#EA945A","Th17" = "#ACD48A","CD8+ Tnaive" = "#F7DBF0",
                                                         "CD8+ Tem" = "#E5948E","CD8+ Tex"="#276D9F","MAIT"="#8AB6D6",
                                                         "NKT" = "#86A667","NK" = "#E7C5DB",
                                            )
                   )
)
dev.off()

#Fig3I-J-----------------
options(reticulate.conda_binary = "/home/Shared/miniconda3/bin/conda")
sce_Th17 <- subset(sce,sub_celltype == "Th17")
sce_Th17$Sample <- factor(sce_Th17$Sample,levels = c("MMPTissue","HealthyTissue"))
for(TF in c("STAT6(+)","NFKB1")){
  expr_data <- FetchData(sce_Th17, vars = c(TF, "Sample"))
  # 计算Wilcoxon检验
  p_value <- wilcox.test(expr_data[expr_data$Sample == levels(expr_data$Sample)[1], "STAT6(+)"],
                         expr_data[expr_data$Sample == levels(expr_data$Sample)[2], "STAT6(+)"])$p.value
  
  # 格式化p值符号
  p_symbol <- ifelse(p_value < 0.001, "***",
                     ifelse(p_value < 0.01, "**",
                            ifelse(p_value < 0.05, "*", "ns")))
  RidgePlot(sce_Th17,features=TF,group.by = "Sample",cols = c("#F5BC6E","#726BAE"))+
    ylab("")+
    ggtitle("")+
    xlab(paste0("The activity score of ",TF))+
    annotate("text", 
             x = Inf, y = Inf, 
             label = paste0("p = 1e-04"),
             hjust = 1.1, vjust = 1.5,  
             size = 5, color = "black") +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 14,family = "Helvetica"),
      axis.text.y = element_blank(),
      axis.title.x = element_text(
        hjust = 0.5,
        vjust = 1,
        margin = margin(t = 10),
        size = 15,
        color = "black"
      ),
      panel.border = element_rect(
        colour = "gray30",
        linetype = "solid", 
        size = 0.8          
      )
    )
  ggsave(filename = paste0("Fig/Fig",TF,".pdf"),width = 5.5,height = 4,dpi = 600,device = cairo_pdf, family = "Helvetica")
}


