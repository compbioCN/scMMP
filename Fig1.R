library(scRNAtoolVis)
library(cowplot)
library(rstatix)
library(ggpubr)
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
#Fig1C-----------
markers <- c(
  "GZMK","CD8A","CD8B","CD3D","CD3E","IL7R","CCR7",
  "MS4A1","CD19","CD37","CD79A",
  "MZB1","XBP1","CD138","SDC1",
  "LYZ","CD14","CD68","CD163",
  "KRT14","KRT5","JUP","SOX2",
  "DCN","COL1A1","CFD","MME",
  "TAGLN","MYL9","FOXC2","ACTA2",
  "PECAM1","VWF","CDH5","ENG",
  "ACTA1","MYL1","MYH2",
  "DCT","PMEL","TYRP1"
)
Idents(MMP_all) <- "Celltype"
DotPlot(MMP_all,features = markers,scale = TRUE)+RotatedAxis()+
  scale_color_gradientn(colours = colorRampPalette(c("#4080B4", "#FCFCB2", "#E52F21"))(100))+
  theme(
    axis.text.x = element_text(family = "Helvetica"),
    axis.text.y = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica",size = 13),
    legend.title = element_text(family = "Helvetica",size = 13)
  )
ggsave(filename = "Fig/Fig1C.pdf",width = 18,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica")
#Fig1D----------
colors <- c("#00A08A","#F2AD00")
#提取seurat object中的umap坐标及metadata
dp <- cbind(MMP_all@meta.data, MMP_all@reductions$umap@cell.embeddings)
colnames(dp)[colnames(dp) %in% c("umap_1", "umap_2")] <- c("UMAP_1", "UMAP_2")
centroid <- ddply(dp, .(Celltype_label), summarise, label_x = median(UMAP_1), label_y = median(UMAP_2))
#绘制UMAP
ggplot(dp, aes(x = UMAP_1, y = UMAP_2, color = Celltype_label)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = colors) +
  geom_label_repel(
    data = centroid, 
    aes(x = label_x, y = label_y, label = Celltype_label),
    color = "black",       
    fill = NA,            
    label.size = NA,       
    size = 5.5,
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
    axis.title.x = element_text(family = "Helvetica",size = 20),
    axis.title.y = element_text(family = "Helvetica",size = 20),
    legend.position = "none"
  )
ggsave(filename = "Fig/Fig1D.pdf",width = 5,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig1E------------
colors <- c("#C7AED5","#A5AA99","#60A897","#EA945A","#ACD48A","#F7DBF0","#E5948E","#276D9F","#8AB6D6","#86A667")
#提取seurat object中的umap坐标及metadata
dp <- cbind(MMP_all@meta.data, MMP_all@reductions$umap@cell.embeddings)
colnames(dp)[colnames(dp) %in% c("umap_1", "umap_2")] <- c("UMAP_1", "UMAP_2")
centroid <- ddply(dp, .(Celltype), summarise, label_x = median(UMAP_1), label_y = median(UMAP_2))
#绘制UMAP
ggplot(dp, aes(x = UMAP_1, y = UMAP_2, color = Celltype)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = colors) +
  geom_label_repel(
    data = centroid, 
    aes(x = label_x, y = label_y, label = Celltype),
    color = "black",        
    fill = NA,            
    label.size = NA,       
    size = 5.5,
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
    axis.title.x = element_text(family = "Helvetica",size = 20),
    axis.title.y = element_text(family = "Helvetica",size = 20),
    legend.position = "none"
  )
ggsave(filename = "Fig/Fig1E.pdf",width = 5,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica")
#Fig1F-----------
colors <- c("#C7AED5","#A5AA99","#60A897","#EA945A")
immune <- subset(MMP_all,Celltype_label == "Immune")
immune$Celltype <- droplevels(immune$Celltype)
sample.prop <- as.data.frame(prop.table(table(immune$Celltype)))
colnames(sample.prop) <- c("celltype","proportion")
sample.prop$lab <- paste0(sample.prop$celltype,"(",round(sample.prop$proportion,2) * 100,"%)")
ggplot(sample.prop, aes(x = "", y = proportion, fill = lab)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = color) +
  coord_polar(theta = "y", start = 0) +
  theme_void() +  
  ggtitle("Immune") + 
  theme(
    plot.title = element_text(
      size = 30,
      family = "Helvetica",
      hjust = 0.5,  
      margin = margin(b = -20)  
    ),
    legend.text = element_text(size = 20, family = "Helvetica"),
    legend.title = element_blank(),
    legend.position = "bottom", 
    legend.box.spacing = unit(-1.3, "cm")  
  )
ggsave("Fig/Fig1F.pdf",width = 7,height = 6,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig1G---------
colors <- c("#ACD48A","#F7DBF0","#E5948E","#276D9F","#8AB6D6","#86A667")
nonimmune <- subset(MMP_all,Celltype_label == "Non-immune")
nonimmune$Celltype <- droplevels(nonimmune$Celltype)
sample.prop <- as.data.frame(prop.table(table(nonimmune$Celltype)))
colnames(sample.prop) <- c("celltype","proportion")
sample.prop$lab <- paste0(sample.prop$celltype,"(",round(sample.prop$proportion,2) * 100,"%)")
ggplot(sample.prop, aes(x = "", y = proportion, fill = lab)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = color) +
  coord_polar(theta = "y", start = 0) +
  theme_void() +  
  ggtitle("Non-immune") +  
  theme(
    plot.title = element_text(
      size = 30,
      family = "Helvetica",
      hjust = 0.5,  
      margin = margin(b = -20) 
    ),
    legend.text = element_text(size = 20, family = "Helvetica"),
    legend.title = element_blank(),
    legend.position = "bottom",  
    legend.box.spacing = unit(-1.3, "cm") 
  )

ggsave("Fig/Fig1G.pdf",width = 8,height = 6,dpi = 600,device = cairo_pdf, family = "Helvetica")
#Fig1H------------
colors <- c("#C7AED5","#A5AA99","#60A897","#EA945A","#ACD48A","#F7DBF0","#E5948E","#276D9F","#8AB6D6","#86A667")
DEGs <- FindAllMarkers(MMP_all,group.by = "Celltype",logfc.threshold = 0.25,min.pct = 0.25)
mymarkers <- DEGs$gene
celltype_order <- c("T/NK","B","Plasma","Myeloid","Epithelial","Fibroblasts","Myofibroblasts","Endothelial","Myocytes","Melanocytes")
DEGs$cluster <- factor(DEGs$cluster,levels = celltype_order)
caPairs_up <- data.frame(
  x = c("T/NK","B","Plasma","Myeloid","Epithelial",'Fibroblasts',"Myofibroblasts","Endothelial","Myocytes","Melanocytes"),
  label=c("168 genes","187 genes","89 genes","206 genes","304 genes","263 genes","219 genes","211 genes","137 genes","160 genes")
)
jjVolcano(DEGs,
           log2FC.cutoff = 2.5, 
           pSize = 0.5, 
           base_size = 25,
           aesCol = c('indianred1','gray70'), 
           topGeneN = 0, 
           celltypeSize = 5.5,
           legend.position = c(0.93,0.95),
           tile.col = colors,myMarkers = mymarkers)+
  #geom_text(data = caPairs, aes(x=x, label=label), y=0,  size=5.5,colour = "black")+
  geom_text(data = caPairs_up, aes(x=x, label=label), y=16, vjust=2, size=5.5,colour = "black")+

  theme(text = element_text(family = "Helvetica",size = 20),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
        axis.text.y = element_text(family = "Helvetica",size = 18),
        axis.title.x = element_text(family = "Helvetica",size = 18),  
        axis.title.y = element_text(family = "Helvetica",size = 18),
        legend.text = element_text(size = 16,family = "Helvetica"),
        legend.position = "right",
        legend.box.spacing = unit(-0.2, "cm")
  )+
  xlab("")
ggsave(filename = "Fig/Fig1H.pdf",width = 18,height = 5,dpi=600,device = cairo_pdf, family = "Helvetica")




