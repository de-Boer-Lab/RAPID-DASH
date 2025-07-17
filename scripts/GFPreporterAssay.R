# Load required libraries for data manipulation, flow cytometry, and plotting
library(tidyverse)   # Data manipulation and visualization
library(flowCore)    # Reading and processing flow cytometry data
library(ggcyto)      # Flow cytometry specific plotting functions

# List all .fcs files in the specified directory and get full file paths
file_paths <- list.files("data/flowcytometry/GFPreporter/Singlets/", pattern="fcs", full.names = TRUE)

# Read each FCS file, extract data and add a sample identifier column
data_list <- lapply(file_paths, function(file) {
  # Read the FCS file, do not truncate the range
  flow_data <- read.FCS(file, truncate_max_range = FALSE)
  
  # Extract the sample name from the file name
  sample_name <- basename(file)
  
  # Convert the flow cytometry expression data to a data frame
  data_frame <- as.data.frame(exprs(flow_data))
  
  # Extract a sample identifier using a regular expression and add it as a new column
  data_frame$Sample <- str_extract(sample_name, "(?<=_).*(?=_[^_]+$)")
  
  return(data_frame)
})

# Combine the list of data frames into one data frame for all samples
all_data <- do.call(rbind, data_list)

all_data <- all_data |> 
  mutate(sample_name = str_extract(Sample, "(?<=-).*?(?=-)"),
         Rep = case_when(
           str_detect(Sample, "-.*?-.*?[AD]") ~ 1,  # Matches wells with A or D -> replicate 1
           str_detect(Sample, "-.*?-.*?[BE]") ~ 2,  # Matches wells with B or E -> replicate 2
           str_detect(Sample, "-.*?-.*?[CF]") ~ 3   # Matches wells with C or F -> replicate 3
         ))

# Group data by sample name and replicate, then count the number of events (rows)
all_data |> 
  group_by(sample_name, Rep) |> 
  summarise(n = n())

# Create a bar plot of the number of events recorded for each sample and replicate
all_data |>
  group_by(sample_name, Rep) |> 
  summarise(n = n()) |> 
  ggplot(aes(x = factor(sample_name), y = n, fill = factor(Rep))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "No of events recorded",
    x = "Sample",
    y = "Events recorded"
  ) +
  theme_minimal()

# For each sample and replicate, compute the range of the 'FL1-A' channel values
all_data |> 
  group_by(sample_name, Rep) |> 
  summarise(range = range(`FL1-A`))

# Calculate the percentage of GFP positive cells (FL1-A > 4000) for each sample and replicate
percent_data <- all_data |>
  group_by(sample_name, Rep) |>
  summarize(percentage = mean(`FL1-A` > 4000) * 100, .groups = "drop") |> 
  mutate(Sample = paste(sample_name, "_", Rep, sep = ""))

# Compute the mean GFP+ percentage per sample across replicates, and assign a Category based on the sample name
percent_data_mean <- percent_data |> 
  group_by(sample_name) |> 
  summarise(mean_GFP = mean(percentage)) |> 
  mutate(
    Category = case_when(
      str_detect(sample_name, "^G\\d+$|Array") ~ "gRNA Array",        # Samples starting with "G" and a number or containing "Array"
      TRUE ~ "Individual gRNA Vector"                                  # All other samples
    )
  )

# stat_results <- compare_means(percentage ~ sample_name, data = percent_data,
#                               method = "t.test", ref.group = ".all.",
#                               p.adjust.method = "BH")
# 
# stat_results <- stat_results |>
#   select(group2, p.adj, p.signif) |>
#   rename(sample_name = group2) |>
#   mutate(p_label = if_else(
#     sample_name == "G0",
#     paste0("p = ", p.adj, " (", p.signif, ")"),  # No rounding for G0
#     paste0("p = ", round(p.adj, 4), " (", p.signif, ")")  # Rounded for others
#   ))
# 
# stat_results
# 
# percent_data_mean <- percent_data_mean |>
#   mutate(y.position = mean_GFP + 5) |>    # adjust vertically  
#   left_join(stat_results, by = "sample_name")  # Join with statistical results

# Create a refined bar plot with custom ordering and labeling:
GFP_reporter_plot <- ggplot(data = percent_data_mean, 
    aes(x = factor(sample_name, levels = c(paste0("G", 0:10), "GFPonly")))) +
    # Plot the mean GFP+ percentages as bars, colored by Category
    geom_bar(aes(y = mean_GFP), fill = "lightgreen", 
             stat = "identity", width = 0.6) +
    # Overlay individual replicate points
    geom_point(data = percent_data, aes(y = percentage), 
               color = "#525352", size = 2, alpha = 0.8) +
    labs(y = "Percentage of GFP+ve cells", 
      x = NULL
    ) +
    # geom_text(aes(x = sample_name, y = y.position, label = p_label),
    #         inherit.aes = FALSE, size = 6) +
    coord_flip() +  # Flip the coordinate system for horizontal bars
    # Customize the x-axis labels to display more informative names for each sample
    scale_x_discrete(labels = c("GFPonly" = "GFP targeting gRNA",
                                "G10" = "G10", "G9" = "G9", "G8" = "G8",
                                "G7" = "G7", "G6" = "G6", "G5" = "G5",
                                "G4" = "G4", "G3" = "G3", "G2" = "G2",
                                "G1" = "G1", "G0" = "NTC")) +
    # Enhance the theme for clarity and visual appeal
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.ticks = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )

GFP_reporter_plot


#save the plot to a file
ggsave(file = "GFP_reporter_plot.svg", plot = GFP_reporter_plot,
       width = 15, height = 12)


get_sig_label <- function(p) {
  if (is.na(p)) return("")
  else if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}


#statistical analysis using ANOVA and post-hoc Tukey's HSD test
anova_result <- aov(percentage ~ sample_name, data = percent_data)
postHOC <- tukey_hsd(anova_result, p.adjust.method = "BH")

pairwise_results <- postHOC |> 
  mutate(group1 = factor(group1, levels = c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "GFPonly")),
         group2 = factor(group2, levels = c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "GFPonly"))) 

all_groups <- union(pairwise_results$group1, pairwise_results$group2)
all_combos <- expand.grid(group1 = all_groups, group2 = all_groups)

pairwise_results_complete <- all_combos %>%
  left_join(pairwise_results, by = c("group1", "group2"))


p_log_matrix <- pairwise_results_complete %>%
  select(group1, group2, p.adj) |> 
  mutate(p.adj = -log10(p.adj)) |> 
  pivot_wider(names_from = group2, values_from = p.adj) %>%
  column_to_rownames("group1") %>%
  as.matrix()

# TO make the matrix symmetrical
for (i in 1:nrow(p_log_matrix)) {
  for (j in 1:ncol(p_log_matrix)) {
    if (is.na(p_log_matrix[i, j]) && !is.na(p_log_matrix[j, i])) {
      p_log_matrix[i, j] <- p_log_matrix[j, i]
    } else if (!is.na(p_log_matrix[i, j]) && is.na(p_log_matrix[j, i])) {
      p_log_matrix[j, i] <- p_log_matrix[i, j]
    }
  }
}
diag(p_log_matrix) <- NA  # Optional

group_levels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "GFPonly")
p_log_matrix <- p_log_matrix[group_levels, group_levels]
p_log_matrix

pairwise_long <- as.data.frame(p_log_matrix) %>%
  rownames_to_column(var = "Group1") %>%
  pivot_longer(-Group1, names_to = "Group2", values_to = "p_log")

group_levels <- rownames(p_log_matrix)
pairwise_long <- pairwise_long %>%
  mutate(
    group1 = factor(Group1, levels = group_levels),
    group2 = factor(Group2, levels = group_levels)
  )

#Create significance labels 
pairwise_long <- pairwise_long |> 
  mutate(sig_label = case_when(
    is.na(p_log) ~ "",
    p_log > -log10(0.001) ~ "***",
    p_log > -log10(0.01)  ~ "**",
    p_log > -log10(0.05)  ~ "*",
    TRUE ~ ""
  )) |> 
  mutate(
    label_text = ifelse(!is.na(p_log), paste0(round(p_log, 1), sig_label), ""),
    sig_flag = ifelse(sig_label != "", TRUE, FALSE)
  )

pairwise_long

pairwise_comparisons_heatmap <- ggplot(pairwise_long, aes(x = Group2, y = Group1, fill = p_log)) +
  geom_tile(color = "white", linewidth = 0.5) +  # base heatmap
  geom_text(aes(label = label_text), size = 5) +  # asterisks
  scale_fill_gradient(
    low = "white", high = "steelblue",
    name = expression(-log[10](p[adj])),
    na.value = "grey90"
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file = "pairwise_comparisons_heatmap.svg",plot = pairwise_comparisons_heatmap, width = 14, height = 12)
