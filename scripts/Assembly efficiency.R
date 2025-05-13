library(tidyverse)

#Bulk Plasmid Sequencing

gRNA_dist <-  read_tsv("data/nanopore_sequencing/bulk_plasmid_sequencing/Bulk_gRNAarrays.tsv", col_names = T)

gRNA_dist_long <-  gRNA_dist %>% 
  pivot_longer(cols = - `No of gRNAs in the array`, names_to = "Assembly", values_to = "Count")

gRNA_dist_wide <- gRNA_dist_long %>%
  pivot_wider(id_cols = `Assembly`, 
              names_from = `No of gRNAs in the array`, names_prefix = "Count_",
              values_from = Count)

percent_count <- gRNA_dist_wide %>% 
  transmute(across(c(Count_0:`Count_Not in any  range`), ~ .x*100 /Count_Total)) 

percent_count <- percent_count %>% 
  pivot_longer(c(Count_0:`Count_Not in any  range`), names_to = "gRNA_count", values_to = "Percentage") %>% 
  group_by(gRNA_count) %>% 
  mutate(percent_mean = mean(Percentage))

percent_count <- percent_count %>% 
  mutate(gRNA_count = str_extract(gRNA_count, "[>0-9]+"))


percent_count <- percent_count %>% 
  mutate(gRNA_count = case_when(is.na(gRNA_count) ~ "Not in any range",
                                TRUE  ~ gRNA_count ))

p <- percent_count %>% 
  mutate(gRNA_count  = factor(gRNA_count, levels = c(0:10, ">10", "Not in any range"))) %>% 
  filter(gRNA_count != "Not in any range") %>% 
  ggplot(aes(gRNA_count, Percentage)) +
  geom_bar(stat = "summary", fill = "#97B0B6", width = 0.6) +
  geom_point(color = "#525352", size = 2, alpha = 0.8) +
  labs(
    title = "Bulk Plasmid Sequencing",
    y = "Percentage of sequenced arrays", 
    x = "Number of gRNAs in an array"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.ticks = element_blank(),
    panel.grid.minor= element_blank(),
    panel.grid.major.x =  element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# ggsave(file="BulkPlasmidSeq.jpg", plot=p, width=10, height=8)

#Clonal Plasmid sequencing

gRNA_dist <-  read_tsv("data/nanopore_sequencing/clonal_plasmid_sequencing/Clonal_gRNAarrays.tsv", col_names = T)

gRNA_dist_long <-  clonal_dist %>% 
  pivot_longer(cols = - `No of gRNAs in the array`, names_to = "Assembly", values_to = "Count")

gRNA_dist_wide <- gRNA_dist_long %>%
  pivot_wider(id_cols = `Assembly`, 
              names_from = `No of gRNAs in the array`, names_prefix = "Count_",
              values_from = Count)

percent_count <- gRNA_dist_wide %>% 
  transmute(across(c(Count_0:`Count_Not in any  range`), ~ .x*100 /Count_Total)) 

percent_count <- percent_count %>% 
  pivot_longer(c(Count_0:`Count_Not in any  range`), names_to = "gRNA_count", values_to = "Percentage")

percent_count <- percent_count %>% 
  mutate(gRNA_count = str_extract(gRNA_count, "[>0-9]+"))


percent_count <- percent_count %>% 
  mutate(gRNA_count = case_when(is.na(gRNA_count) ~ "Not in any range",
                                TRUE  ~ gRNA_count ))

p <- percent_count %>% 
  mutate(gRNA_count  = factor(gRNA_count, levels = c(0:10, ">10", "Not in any range"))) %>% 
  filter(gRNA_count != "Not in any range") %>% 
  ggplot(aes(gRNA_count, Percentage)) +
  geom_bar(stat = "summary", fill = "#97B0B6", width = 0.6) +
  # geom_point(color = "#525352", size = 2, alpha = 0.8) +
  labs(
    title = "Clonal Plasmid Sequencing",
    y = "Percentage of sequenced arrays", 
    x = "Number of gRNAs in an array"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.ticks = element_blank(),
    panel.grid.minor= element_blank(),
    panel.grid.major.x =  element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )
p

ggsave(file="ClonalPlasmidSeq.jpg", plot=p, width=10, height=8)
