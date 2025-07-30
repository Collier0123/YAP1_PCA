library(FactoMineR)
library(ggplot2)
library(factoextra)
library(RColorBrewer)
library(dplyr)
# Data preparation
Raw_yap_quant$Group <- as.factor(Raw_yap_quant$Group)
Raw_yap_quant$Position <- as.factor(Raw_yap_quant$Position)
Raw_yap_quant$Inj.Status <- as.factor(Raw_yap_quant$Inj.Status)
Raw_yap_quant$Condition <- interaction(Raw_yap_quant$Position, Raw_yap_quant$Inj.Status)
Raw_yap_quant$Condition <- factor(Raw_yap_quant$Condition)
Raw_yap_quant$GroupCondition <- interaction(Raw_yap_quant$Group, Raw_yap_quant$Position, Raw_yap_quant$Inj.Status)
Raw_yap_quant$GroupCondition <- factor(Raw_yap_quant$GroupCondition, 
                                       levels = c("siNTC.Outer.Non-inj", "siNTC.Outer.Inj", 
                                                  "siNTC.Inner.Non-inj", "siNTC.Inner.Inj",
                                                  "siTead4.Outer.Non-inj", "siTead4.Outer.Inj",
                                                  "siTead4.Inner.Non-inj", "siTead4.Inner.Inj"))
# Perform PCA
pca_data_yap <- Raw_yap_quant[, c("Area.Cyto","Area.Nuc","TCCF.Nuc","TCCF.Cyto","Av.N.C")]
pca_result_yap <- PCA(pca_data_yap, ncp = 5, graph = FALSE)

#pre-plot analysis - data overview
library(corrplot)
pca_result_yap$eig
fviz_eig(pca_result_yap, addlabels = TRUE, ylim = c(0, 50))
fviz_pca_var(pca_result_yap, repel = TRUE)
fviz_contrib(pca_result_yap,"var")
fviz_contrib(pca_result_yap,"var", axes = 2)
fviz_contrib(pca_result_yap,"var", axes = 3)
fviz_contrib(pca_result_yap, "var", axes = 1:2)
fviz_contrib(pca_result_yap, "var", axes = 1:3)
corrplot(pca_result_yap$var$coord)

#Enhanced feature scree-plot
scree_plot <- fviz_eig(pca_result_yap, addlabels = FALSE, barfill = "black", barcolor = "black") +
  labs(
    title = "Scree Plot",         
    x = "Dimensions",             
    y = "Percentage of explained 
    varience (%)"  
  ) +
  theme_minimal() +                         
  theme(
    plot.title = element_text(size = 45,face = "bold", hjust = 0.5, margin = margin(r=60,b =50)),
    plot.subtitle = element_text(size=40, hjust =0.5, margin = margin(b=50, r=60)),
    axis.title = element_text(size = 40, hjust = 0.5),                           
    axis.text = element_text(size = 33), 
    panel.grid.major = element_line(color = "gray"), 
    panel.grid.minor = element_blank(),                              
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 20))
  )
print(scree_plot)

#Values for enhanced feature graphs
#dim vaues
eig_values <- pca_result_yap$eig
dim1_var <- format(round(eig_values[1, 2], 1), nsmall = 1)
dim2_var <- format(round(eig_values[2, 2], 1), nsmall = 1)
dim3_var <- format(round(eig_values[3, 2], 1), nsmall = 1)
dim4_var <- format(round(eig_values[4, 2], 1), nsmall = 1)

#Row labels - dynamic
variable_name_mapping <- c(
  "TCCF.Nuc" = "Nuclear YAP1 
TCCF",
  "TCCF.Cyto" = "Cytoplasmic YAP1 
TCCF",
  "Av.N.C" = "Nuclear:Cytoplasmic
YAP1 Ratio",
  "Area.Cyto" = "Cytoplasmic area 
(µm^2)",
  "Area.Nuc" = "Nuclear area 
(µm^2)"
)
yap_variable_names <- rownames(pca_result_yap$var$coord)
yap_labels_row <- variable_name_mapping[yap_variable_names]


#Column lables - dynamic
yap_correlation_matrix <- pca_result_yap$var$coord
yap_labels_column <- paste0("Dim", 1:ncol(yap_correlation_matrix))

#Matrix mapping
rownames(yap_correlation_matrix) <- yap_labels_row
colnames(yap_correlation_matrix) <- yap_labels_column


#Enhanced corr plot (CHECK THAT TITLES MATCH BASIC!)
mwb_palette <- colorRampPalette(c("magenta","white", "black"))


corrplot(yap_correlation_matrix, 
         col = mwb_palette(200),          
         method = "circle",               
         tl.col = "black",                 
         tl.offset = 0.9,
         tl.cex = 1.6,                     
         tl.srt = 0,                   
         cl.cex = 1,
         cl.pos = "r")   
#Contribution plots
#Dim1
fviz_contrib_plot1 <- fviz_contrib(pca_result_yap, "var", axes = 1)
fviz_contrib_data1 <- fviz_contrib_plot1$data
fviz_contrib_data1$variable <- variable_name_mapping[rownames(fviz_contrib_data1)]
fviz_contrib_data1$variable <- reorder(fviz_contrib_data1$variable, -fviz_contrib_data1$contrib)
ggplot(fviz_contrib_data1, aes(x = variable, y = contrib)) +
  geom_bar(stat = "identity", color="black", fill = "black", width = 0.3) +
  geom_hline(yintercept = mean(fviz_contrib_data1$contrib), linetype = "dashed", color = "red") +
  labs(
    title = "Variable Contributions to Dimension 1",
    x = NULL,
    y = "Contribution (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 26,face = "bold", hjust = 0.5, margin = margin(r=60, b=20)),
    plot.subtitle = element_text(size=26, hjust =0.5, margin = margin(b=50, r=60)),
    axis.title = element_text(size = 22, hjust = 0.5),                           
    axis.text = element_text(size = 20), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20, color = "black"),
    panel.grid.major = element_line(color = "gray"), 
    panel.grid.minor = element_blank(),                              
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 20))
  )
#Compare with basic
fviz_contrib(pca_result_yap,"var")

#Dim2
fviz_contrib_plot2 <- fviz_contrib(pca_result_yap, "var", axes = 2)
fviz_contrib_data2 <- fviz_contrib_plot2$data
fviz_contrib_data2$variable <- variable_name_mapping[rownames(fviz_contrib_data2)]
fviz_contrib_data2$variable <- reorder(fviz_contrib_data2$variable, -fviz_contrib_data2$contrib)
ggplot(fviz_contrib_data2, aes(x = variable, y = contrib)) +
  geom_bar(stat = "identity", color="black", fill = "black", width = 0.3) +
  geom_hline(yintercept = mean(fviz_contrib_data2$contrib), linetype = "dashed", color = "red") +
  labs(
    title = "Variable Contributions to Dimension 2",
    x = NULL,
    y = "Contribution (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 26,face = "bold", hjust = 0.5, margin = margin(r=60, b=20)),
    plot.subtitle = element_text(size=26, hjust =0.5, margin = margin(b=50, r=60)),
    axis.title = element_text(size = 22, hjust = 0.5),                           
    axis.text = element_text(size = 20), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20, color = "black"),
    panel.grid.major = element_line(color = "gray"), 
    panel.grid.minor = element_blank(),                              
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 20))
  )
#Compare with basic
fviz_contrib(pca_result_yap,"var", axes = 2)

#Dim3
fviz_contrib_plot3 <- fviz_contrib(pca_result_yap, "var", axes = 3)
fviz_contrib_data3 <- fviz_contrib_plot3$data
fviz_contrib_data3$variable <- variable_name_mapping[rownames(fviz_contrib_data3)]
fviz_contrib_data3$variable <- reorder(fviz_contrib_data3$variable, -fviz_contrib_data3$contrib)
ggplot(fviz_contrib_data3, aes(x = variable, y = contrib)) +
  geom_bar(stat = "identity", color="black", fill = "black", width = 0.3) +
  geom_hline(yintercept = mean(fviz_contrib_data3$contrib), linetype = "dashed", color = "red") +
  labs(
    title = "Variable Contributions to Dimension 3",
    x = NULL,
    y = "Contribution (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 26,face = "bold", hjust = 0.5, margin = margin(r=60, b=20)),
    plot.subtitle = element_text(size=26, hjust =0.5, margin = margin(b=50, r=60)),
    axis.title = element_text(size = 22, hjust = 0.5),                           
    axis.text = element_text(size = 20), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20, color = "black"),
    panel.grid.major = element_line(color = "gray"), 
    panel.grid.minor = element_blank(),                              
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 20))
  )
#Compare with basic
fviz_contrib(pca_result_yap, "var", axes = 3)


#SCATTER PLOT - final
#Create levels
ordered_levels_norm <- c(sort(unique(Raw_yap_quant$GroupCondition[grep("siNTC", Raw_yap_quant$GroupCondition)])),
                         sort(unique(Raw_yap_quant$GroupCondition[grep("siTead4", Raw_yap_quant$GroupCondition)])))
Raw_yap_quant$GroupCondition <- factor(Raw_yap_quant$GroupCondition, levels = ordered_levels_norm)

# Define the specific color scheme
condition_colours <- c(
  "siNTC.Outer.Non-inj" = "#B2BEB5",   
  "siNTC.Outer.Inj" = "black",       
  "siNTC.Inner.Non-inj" = "#8fd9fb",   
  "siNTC.Inner.Inj" = "#00AFFF",       
  "siTead4.Outer.Non-inj" = "#FF0000", 
  "siTead4.Outer.Inj" = "#880808",     
  "siTead4.Inner.Non-inj" = "#ffb343", 
  "siTead4.Inner.Inj" = "#fe8205"      
)

#Legend labels
group_condition_labels <- c(
  "siNTC.Outer.Non-inj" = "siNTC Non-inj. Outer",
  "siNTC.Outer.Inj" = "siNTC Inj. Outer",
  "siNTC.Inner.Non-inj" = "siNTC Non-inj. Inner",
  "siNTC.Inner.Inj" = "siNTC Inj. Inner",
  "siTead4.Outer.Non-inj" = "siTead4 Non-inj. Outer",
  "siTead4.Outer.Inj" = "siTead4 Inj. Outer",
  "siTead4.Inner.Non-inj" = "siTead4 Non-inj. Inner",
  "siTead4.Inner.Inj" = "siTead4 Inj. Inner"
)

#2D 4-axis plot values
library(ggplot2)
ind_coords_Raw_yap_quant <- as.data.frame(pca_result_yap$ind$coord)
plot_data_Raw_yap_quant <- data.frame(
  Dim1 = ind_coords_Raw_yap_quant[, 1],
  Dim2 = ind_coords_Raw_yap_quant[, 2],
  Dim3 = ind_coords_Raw_yap_quant[, 3],
  Dim4 = ind_coords_Raw_yap_quant[, 4],
  GroupCondition = Raw_yap_quant$GroupCondition
)

centroid_coords <- as.data.frame(pca_result_yap$ind$coord) %>%
  dplyr::mutate(GroupCondition = Raw_yap_quant$GroupCondition) %>%
  dplyr::group_by(GroupCondition) %>%
  dplyr::summarize(
    Dim.1 = mean(Dim.1),
    Dim.2 = mean(Dim.2)
  )

#Basic plot data - NEED THIS TO AQCUIRE CORRECT ELLIPSES, make sure to use correct habillage
fviz_pca_plot <- fviz_pca_ind(pca_result_yap,
                              habillage = Raw_yap_quant$GroupCondition,  
                              palette = condition_colours,
                              addEllipses = TRUE,
                              ellipse.level = 0.95,
                              ellipse.alpha = 0,
                              geom = "point",
                              repel = FALSE
) +
  scale_shape_manual(values = rep(19, length(unique(Raw_yap_quant$GroupCondition)))) +
  guides(shape = "none")

#2Dims
fviz_pca_plot +
  geom_hline(yintercept = 0, linetype = "solid", color = "#C0C0C0") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "#C0C0C0") +
  geom_point(data = centroid_coords, aes(x = Dim.1, y = Dim.2),color = "black", fill = condition_colours, 
             size = 3, shape = 24, stroke = 0.5) +
  scale_color_manual(values = condition_colours, name = "Cell Profile",
                     labels = group_condition_labels[levels(Raw_yap_quant$GroupCondition)]
  ) +
  labs(
    x = paste0("Dim1 (", dim1_var, "%)"),
    y = paste0("Dim2 (", dim2_var, "%)"),
    title = "PCA"
  ) +
  guides(
    color = guide_legend(title = "Cell Profile"),
    size = guide_legend(title = paste0("Dim3 (", dim3_var, "%)")),
    fill = "none"  
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, margin = margin(b=20, l=100)),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    plot.margin = margin(0, 0, 0, 0))

