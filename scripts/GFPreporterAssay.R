# Load required libraries for data manipulation, flow cytometry, and plotting
library(tidyverse)   # Data manipulation and visualization
library(flowCore)    # Reading and processing flow cytometry data
library(ggcyto)      # Flow cytometry specific plotting functions
library(ggridges)    # Ridge plots (if needed for later visualization)

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

all_data <- all_data %>% 
  mutate(sample_name = str_extract(Sample, "(?<=-).*?(?=-)"),
         Rep = case_when(
           str_detect(Sample, "-.*?-.*?[AD]") ~ 1,  # Matches wells with A or D -> replicate 1
           str_detect(Sample, "-.*?-.*?[BE]") ~ 2,  # Matches wells with B or E -> replicate 2
           str_detect(Sample, "-.*?-.*?[CF]") ~ 3   # Matches wells with C or F -> replicate 3
         ))

# Group data by sample name and replicate, then count the number of events (rows)
all_data %>% 
  group_by(sample_name, Rep) %>% 
  summarise(n = n())

# Create a bar plot of the number of events recorded for each sample and replicate
all_data %>%
  group_by(sample_name, Rep) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = factor(sample_name), y = n, fill = factor(Rep))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "No of events recorded",
    x = "Sample",
    y = "Events recorded"
  ) +
  theme_minimal()

# For each sample and replicate, compute the range of the 'FL1-A' channel values
all_data %>% 
  group_by(sample_name, Rep) %>% 
  summarise(range = range(`FL1-A`))

# Calculate the percentage of GFP positive cells (FL1-A > 4000) for each sample and replicate
percent_data <- all_data %>%
  group_by(sample_name, Rep) %>%
  summarize(percentage = mean(`FL1-A` > 4000) * 100, .groups = "drop") %>% 
  mutate(Sample = paste(sample_name, "_", Rep, sep = ""))

# Compute the mean GFP+ percentage per sample across replicates, and assign a Category based on the sample name
percent_data_mean <- percent_data %>% 
  group_by(sample_name) %>% 
  summarise(mean_GFP = mean(percentage)) %>% 
  mutate(
    Category = case_when(
      str_detect(sample_name, "^G\\d+$|Array") ~ "gRNA Array",        # Samples starting with "G" and a number or containing "Array"
      TRUE ~ "Individual gRNA Vector"                                  # All other samples
    )
  )

# Display the summarized mean percentages per sample
percent_data_mean

# Create a refined bar plot with custom ordering and labeling:
ggplot(data = percent_data_mean, 
       aes(x = factor(sample_name, levels = c(paste0("G", 0:10), "GFPonly", "Ind_gRNA_GFP")))) +
  # Plot the mean GFP+ percentages as bars, colored by Category
  geom_bar(aes(y = mean_GFP, fill = Category), 
           stat = "identity", width = 0.6) +
  # Overlay individual replicate points
  geom_point(data = percent_data, aes(y = percentage), 
             color = "#525352", size = 2, alpha = 0.8) +
  # Manually define the fill colors and labels for the categories
  scale_fill_manual(values = c("lightgreen", "grey"),
                    labels = c("gRNA Array", "Individual gRNA Vector")) +
  labs(
    title = "GFP activation",
    y = "Percentage of GFP+ve cells", 
    x = NULL
  ) +
  coord_flip() +  # Flip the coordinate system for horizontal bars
  # Customize the x-axis labels to display more informative names for each sample
  scale_x_discrete(labels = c("Ind_gRNA_GFP" = "Pooled gRNA vectors",
                              "GFPonly" = "GFP targeting gRNA",
                              "G10" = "G10", "G9" = "G9", "G8" = "G8",
                              "G7" = "G7", "G6" = "G6", "G5" = "G5",
                              "G4" = "G4", "G3" = "G3", "G2" = "G2",
                              "G1" = "G1", "G0" = "G0")) +
  # Enhance the theme for clarity and visual appeal
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
