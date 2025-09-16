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

rounded_percent_count <- as.data.frame(round(t(percent_count), 2))

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
  geom_text(aes(y = percent_mean, label = round(percent_mean, 1)), 
            vjust = -2.75, size = 4) +
  labs(
    y = "Percentage of arrays", 
    x = "Number of gRNAs"
  ) +
  lims(y = c(0, 100)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.ticks = element_blank(),
    panel.grid.minor= element_blank(),
    panel.grid.major.x =  element_blank(),
  )
p

# ggsave(file="BulkPlasmidSeq.svg", plot=p, width=10, height=8)

#Clonal Plasmid sequencing

clonal_dist <-  read_tsv("data/nanopore_sequencing/clonal_plasmid_sequencing/Clonal_gRNAarrays.tsv", col_names = T)

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
  dplyr::filter(gRNA_count != "Not in any range") %>% 
  ggplot(aes(gRNA_count, Percentage)) +
  geom_bar(stat = "summary", fill = "#97B0B6", width = 0.6) +
  # geom_point(color = "#525352", size = 2, alpha = 0.8) +
  labs(
    y = "Percentage of arrays", 
    x = "Number of gRNAs"
  ) +
  lims(y = c(0, 100)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.ticks = element_blank(),
    panel.grid.minor= element_blank(),
    panel.grid.major.x =  element_blank(),
  )
p

# ggsave(file="ClonalPlasmidSeq.svg", plot=p, width=10, height=8)

## Mutation Rate Estimation


mutations <- read_tsv("data/nanopore_sequencing/clonal_plasmid_sequencing/MutationRates.txt", col_names = T)

mutations <- mutations |> 
  mutate(Clone = as.factor(Clone)) |> 
  rename(
    `PCA Primer` = PCA_handle,
    `Spacer Oligo` = Spacer_oligo
  )

n <- mutations |> 
  ggplot(aes(x = Total)) +
  stat_count(fill = "#97B0B6") +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)), # use computed count from stat_count
    vjust = -0.5,
    size = 5
  ) +
  labs(
    x = "Total Mutations per array",
    y = "Number of full length assemblies"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.ticks = element_blank(),
    panel.grid.minor= element_blank(),
    panel.grid.major.x =  element_blank(),
  )

#Plot histogram of total mutations using stat_count and add total as geom_text

n <- mutations |> 
  ggplot(aes(x = Total)) +
  stat_count(fill = "#97B0B6") +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)), # use computed count from stat_count
    vjust = -0.5,
    size = 5
  ) +
  labs(
    x = "Total Mutations per array",
    y = "Number of full length assemblies"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.ticks = element_blank(),
    panel.grid.minor= element_blank(),
    panel.grid.major.x =  element_blank(),
  )
n
ggsave(file="MutationsTotal.svg", plot=n, width=12, height=10)

mutations <- tibble("PCA primer" = 6, "gRNA" = 21)


mutations_longer <- mutations |> 
  pivot_longer(cols = c(`PCA primer`, `gRNA`), 
               names_to = "Category", 
               values_to = "Mutations") |> 
  mutate(Category = factor(Category, levels = c("PCA primer", "gRNA")))

m <- mutations_longer |> 
  ggplot(aes(x = Category, y = (Mutations/280)*100, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = round((Mutations/280)*100,2)), # use computed count from stat_count
    vjust = -0.5,
    size = 5
  ) +
  labs(
    x = "Mutated Region",
    y = "Percentage of gRNA units"
  ) +
  scale_fill_manual(
    values = c(
      "PCA primer"   = "purple",
      "gRNA" = "darkgreen"))+
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.ticks = element_blank(),
    panel.grid.minor= element_blank(),
    panel.grid.major.x =  element_blank(),
    legend.position = "none"
  )


m
ggsave(file="MutationsByPositionbygRNA.svg", plot=m, width=8, height=10)
