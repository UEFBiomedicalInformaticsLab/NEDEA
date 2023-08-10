# Define global options for this script 
data_balance_method <- "none"
model <- "rf"
feature_type <- "CombinedDisAdr2Gene"

plot_data <- feature_matrix
plot_data$Class <- substr(row.names(plot_data), 1, 3)


preProcess <- caret::preProcess(feature_matrix, method = c("zv", "center", "scale") )
feature_matrix <- predict(object = preProcess, newdata = feature_matrix)

# Plot PCA scatter plot
pca_scatter <- ggplot() +
  geom_point(data = plot_data, 
             aes(x = plot_data[,1], y = plot_data[,2], color = Class), 
             size = 0.5,
             shape = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(margin = margin(1,1,1,1)),
        text = element_text(size = 8), 
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
  ggtitle(disease) +
  xlab(colnames(plot_data)[1]) + 
  ylab(colnames(plot_data)[2])
