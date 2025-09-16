library(tidyverse)


array <- as.matrix(read.csv("data/nanopore_sequencing/array_library/array_Lib.csv", header = FALSE))

array <- as.data.frame(array)

array <-  distinct(array)

colnames(array) <- paste0("Position_", 1:10)

array_long <- array %>% 
  pivot_longer(cols = 1:10, names_to = "Position", values_to = "gRNA") 

array_long$Position <- factor(str_extract(array_long$Position, "\\d+"), levels = paste0(1:10))

gRNA_coverage <- array_long %>% 
  dplyr::filter(gRNA %in% c(1:10)) %>% 
  group_by(gRNA) %>% 
  summarise(count = n())

position_coverage <- array_long %>% 
  dplyr::filter(gRNA %in% c(1:10)) %>% 
  group_by(Position) %>% 
  summarise(count = n())

gRNA_count <- array_long %>% 
  dplyr::filter(gRNA %in% c(1:10)) %>% 
  group_by(Position, gRNA) %>% 
  summarise(count = n())

# Barplot of gRNA coverage

gRNA_coverage_plot <- ggplot(gRNA_coverage, aes(x = factor(gRNA), y = count)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(y = "Total Count", x = "", title = "") +
  geom_text(aes(label = count), hjust = 1.5, size = 5, color = "black") +
  theme_minimal() +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
gRNA_coverage_plot
# ggsave(file="gRNA_coverage_plot.svg", plot=gRNA_coverage_plot, width=10, height=8)

position_coverage_plot <- ggplot(position_coverage, aes(x = factor(Position), y = count)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(x = "Position", y = "Coverage", title = "Coverage per position") +
  geom_text(aes(label = count), vjust = -0.5, size = 5, color = "black") +
  theme_minimal() +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18),
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
position_coverage_plot
# ggsave(file="position_coverage_plot.svg", plot=position_coverage_plot, width=10, height=8)

heatmap <- ggplot(gRNA_count, aes(x = Position, y = factor(gRNA), fill = count)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = c("white", "lightblue", "darkblue"),
    limits = c(0, 70)  # Set legend scale range
  ) +
  labs(x = "Position within the array", y = "gRNA") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18),
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.direction = "horizontal"
  )
heatmap
# ggsave(file="array_library.svg", plot=heatmap, width=10, height=8)

final_plot <- heatmap + gRNA_coverage_plot + 
  plot_layout(widths = c(3, 1), guides = "collect") & 
  theme(legend.position = "top", legend.direction = "horizontal")
final_plot

# ggsave(file="array_library_final.svg", plot=final_plot, width=16, height=12)
