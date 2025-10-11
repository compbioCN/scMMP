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

#Fig4A---------
colors <- c("#6C6BA6","#E5948E")
dp <- cbind(MMP_Myeloid@meta.data, MMP_Myeloid@reductions$umap@cell.embeddings)
colnames(dp)[colnames(dp) %in% c("umap_1", "umap_2")] <- c("UMAP_1", "UMAP_2")
centroid <- ddply(dp, .(Sample), summarise, label_x = median(UMAP_1), label_y = median(UMAP_2))
ggplot(dp, aes(x = UMAP_1, y = UMAP_2, color = Sample)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = colors) +
  geom_label_repel(
    data = centroid, 
    aes(x = label_x, y = label_y, label = Sample),
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
  theme(
    plot.title = element_blank(),
    text = element_text(family = "Helvetica",size = 14),
    axis.title.x = element_text(family = "Helvetica",size = 16),
    axis.title.y = element_text(family = "Helvetica",size = 16),
    legend.position = "none"
  )+
  ggtitle("") +
  xlim(-8, 10) +
  ylim(-15, 10)
ggsave(filename = "Fig/Fig4A.pdf",width = 5,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig4B---------
colors <- c("#D5E7AC","#67C9F2","#E5948E","#F4D4D1")
MMP_Myeloid$sub_celltype <- droplevels(MMP_Myeloid$sub_celltype)
celltype_order <- c("Monocytes","Macrophages","cDCs","Mitotic")
MMP_Myeloid$Celltype <- factor(MMP_Myeloid$Celltype,levels = celltype_order)
dp <- cbind(MMP_Myeloid@meta.data, MMP_Myeloid@reductions$umap@cell.embeddings)
colnames(dp)[colnames(dp) %in% c("umap_1", "umap_2")] <- c("UMAP_1", "UMAP_2")
centroid <- ddply(dp, .(sub_celltype), summarise, label_x = median(UMAP_1), label_y = median(UMAP_2))
ggplot(dp, aes(x = UMAP_1, y = UMAP_2, color = sub_celltype)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = celltype_color) +
  geom_label_repel(
    data = centroid, 
    aes(x = label_x, y = label_y, label = sub_celltype),
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
  theme(
    plot.title = element_blank(),
    text = element_text(family = "Helvetica",size = 14),
    axis.title.x = element_text(family = "Helvetica",size = 16),
    axis.title.y = element_text(family = "Helvetica",size = 16),
    legend.position = "none"
  )+
  ggtitle("") +
  xlim(-8, 10) +
  ylim(-15, 10)
ggsave(filename = "Fig/Fig4B.pdf",width = 5,height = 5,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig4C----------
sample.prop <- as.data.frame(prop.table(table(MMP_Myeloid$sub_celltype,MMP_Myeloid$Sample)))
colnames(sample.prop) <- c("celltype","sample","proportion")
ggplot(sample.prop,aes(x=sample,y=proportion,fill=celltype))+
  geom_bar(stat = "identity",position = "fill")+
  theme_bw()+
  ylab("Proportion")+
  xlab("")+
  guides(fill=guide_legend(title = NULL))+
  scale_fill_manual(values = celltype_color)+
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
ggsave(filename = "Fig/Fig4C.pdf",width = 4.5,height = 6,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig4D-----------
markers <- c(
  #"CD14","CD68","CCR2","ITGAX","FCGR3A",
  "FCN1","AQP9","LILRA5","IL1RN",
  "CCL2","MRC1","PLTP","MSR1","MAF","NET1",
  "CLEC10A","CD1C","CD1E","FCGR2B",
  "MKI67","PCNA","UBE2C"
)
DotPlot(MMP_Myeloid,features = markers,scale = TRUE,)+RotatedAxis()+
  scale_color_gradientn(colours = colorRampPalette(c("#4080B4", "#FCFCB2", "#E52F21"))(100))+
  theme(
    axis.text.x = element_text(family = "Helvetica"),
    axis.text.y = element_text(family = "Helvetica"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.text = element_text(family = "Helvetica",size = 10),
    legend.title = element_text(family = "Helvetica",size = 10)
  )
ggsave(filename = "Fig/Fig4D.pdf",width = 12,height = 3.6,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig4E----------
sample.prop <- as.data.frame(prop.table(table(MMP_Myeloid$sub_celltype)))
colnames(sample.prop) <- c("celltype","proportion")
sample.prop$lab <- paste0(sample.prop$celltype,"(",round(sample.prop$proportion,2) * 100,"%)")
sample.prop$lab <- factor(sample.prop$lab,levels = c("Monocytes(6%)","Macrophages(32%)","cDCs(60%)","Mitotic(2%)"))
ggplot(sample.prop, aes(x = "", y = proportion, fill = lab)) + 
  geom_bar(stat = "identity", width = 1)+ 
  scale_fill_manual(values = celltype_color)+
  coord_polar(theta = "y", start = 0)+
  theme_void() +  
  theme(
    legend.text = element_text(size = 15,family = "Helvetica"),  
    legend.title = element_text(size = 18,face = "bold",family = "Helvetica")
  )+
  labs(fill = "Celltype")
ggsave("Fig/Fig4D.pdf",width = 6,height = 6,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig4F-----------
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,signaling = c("TNF","ICAM","CXCL","ITGB2","MIF","MHC-II","IL1"),color.use = c("lightblue","pink"))+
  theme(axis.text.y = element_text(color = "black",family = "Helvetica"),
        axis.text.x = element_text(family = "Helvetica"),
        axis.title.x = element_text(family = "Helvetica"),
        legend.text = element_text(family = "Helvetica")
  )
ggsave(filename = "Fig/Fig4F.pdf", width = 4, height = 4,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig4G----------
source_use <- c("Fibroblasts")
target_use <- c("cDCs","Monocytes","Macrophages","Mitotic")
netVisual_bubble(cellchat, sources.use = target_use, targets.use = source_use,  comparison = c(1, 2), angle.x = 45)+
  theme(
    axis.text.x = element_text(family = "Helvetica",size = 13),
    axis.text.y = element_text(family = "Helvetica",size = 13),
    legend.text = element_text(family = "Helvetica",size = 10),
    legend.title = element_text(family = "Helvetica",size = 13)
  )
ggsave("Fig/Fig4G.pdf",  width = 6.5, height = 8,dpi = 600,device = cairo_pdf, family = "Helvetica")

#Fig4H-J----------
pathways = c("TNF","ICAM","MIF")
for(i in pathways){
  pathways.show <- i
  weight.max <- getMaxWeight(cci_list, slot.name = c("netP"), attribute = pathways.show) 
  pdf(paste0("Fig/Fig4_",i,".pdf"), width = 8, height = 5.2, family = "Helvetica")
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cci_list)) {
    netVisual_aggregate1(cci_list[[i]], signaling = pathways.show, layout = "circle", 
                         edge.weight.max = weight.max[1], edge.width.max = 10, 
                         signaling.name = paste(pathways.show, names(cci_list)[i]))
  }
  dev.off()
}




