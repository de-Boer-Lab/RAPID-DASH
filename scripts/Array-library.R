library(tidyverse)


array <- as.matrix(read.csv("../Plasmidsaurus/array_Lib_regPCA.csv", header = FALSE))

array <- as.data.frame(array)

array <-  distinct(array)

colnames(array) <- paste0("Position_", 1:10)

array_long <- array %>% 
  pivot_longer(cols = 1:10, names_to = "Position", values_to = "gRNA") 

array_long$Position <- factor(str_extract(array_long$Position, "\\d+"), levels = paste0(1:10))


gRNA_count <- array_long %>% 
  dplyr::filter(gRNA %in% c(1:10)) %>% 
  group_by(Position, gRNA) %>% 
  summarise(count = n())


heatmap <- ggplot(gRNA_count, aes(x = Position, y = factor(gRNA), fill = count)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("white", "lightblue", "darkblue")) +
  labs(x = "Position within the array", y = "gRNA") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.ticks = element_blank(),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

ggsave(file="array_library.svg", plot=heatmap, width=10, height=8)

