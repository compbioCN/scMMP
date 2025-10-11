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


#Fig5A-----------
colors <- c("#C7AED5","#A5AA99","#60A897","#EA945A","#ACD48A",
            "#F7DBF0","#E5948E","#276D9F","#8AB6D6","#86A667",
            "#E7C5DB","#78C4D4","#CCD7DD","#D5E7AC")
sce$Celltype <- factor(sce$Celltype,levels = c("CD4+ Tnaive","CD4+ Trm","CD4+ Tcm","CD4+ Tem","CD4+ Temra",
                                               "CD8+ Tnaive","CD8+ Tem","CD8+ Temra","CD8+ Tex",
                                               "Th17","Treg","Tprolif","γδ T","MAIT"))
levels(sce$Celltype)[levels(sce$Celltype) == "γδ T"] <- "γδ_T"
dp <- cbind(sce@meta.data, sce@reductions$umap@cell.embeddings)
colnames(dp)[colnames(dp) %in% c("umap_1", "umap_2")] <- c("UMAP_1", "UMAP_2")
centroid <- ddply(dp, .(Celltype), summarise, label_x = median(UMAP_1), label_y = median(UMAP_2))
axis_length <- 2  
ggplot(dp, aes(x = UMAP_1, y = UMAP_2, color = Celltype)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = colors) +
  geom_label_repel(
    data = centroid, 
    aes(x = label_x, y = label_y, label = Celltype),
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
ggsave(filename = "Fig/Fig5A.pdf",width = 5.5,height = 5,dpi = 600,device = cairo_pdf,family="Helvetica") #保存特殊字符γδ_T

#Fig5B--------
colors <- c("#260F99FF","#422CB2FF","#6551CCFF","#8F7EE5FF","#BFB2FFFF")

dp <- cbind(sce@meta.data, sce@reductions$umap@cell.embeddings)
colnames(dp)[colnames(dp) %in% c("umap_1", "umap_2")] <- c("UMAP_1", "UMAP_2")
centroid <- ddply(dp, .(cloneSize), summarise, label_x = median(UMAP_1), label_y = median(UMAP_2))
axis_length <- 2 
ggplot(dp, aes(x = UMAP_1, y = UMAP_2, color = cloneSize)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = colors) +
  geom_label_repel(
    data = centroid, 
    aes(x = label_x, y = label_y, label = cloneSize),
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
ggsave(filename = "Fig/Fig5B.pdf",width = 5.5,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica") 

#Fig5D-------
sce_sub <- subset(sce,Sample == "MMPBlood")
clonalOccupy(sce_sub, 
             x.axis = "Celltype", 
             proportion = TRUE, 
             label = FALSE) +
  scale_fill_manual(values = colors) +  
  scale_color_manual(values = colors) + 
  theme(
    text = element_text(family = "Helvetica", size = 14),
    title = element_text(size = 15),
    plot.title = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    plot.margin = margin(t = 20, r = 10, b = 20, l = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) # 外边框
  )
ggsave(filename = "Fig/Fig5D.pdf",width = 7.5,height = 4,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig5E---------
sce_sub <- subset(sce,Sample == "HealthyBlood")
clonalOccupy(sce_sub, 
             x.axis = "Celltype", 
             proportion = TRUE, 
             label = FALSE) +
  scale_fill_manual(values = colors) +  
  scale_color_manual(values = colors) + 
  theme(
    text = element_text(family = "Helvetica", size = 14),
    title = element_text(size = 15),
    plot.title = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    plot.margin = margin(t = 20, r = 10, b = 20, l = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) # 外边框
  )
ggsave(filename = "Fig/Fig5E.pdf",width = 7.5,height = 4,dpi = 600,device = cairo_pdf, family = "Helvetica")


