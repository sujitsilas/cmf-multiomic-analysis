# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(car)
library(stringr)
library(ggforce)

# Define a function for plotting
create_plot <- function(data, type) {
  ggplot(data, aes(x = Sample, y = Expression, fill = Concentration)) + 
    geom_point() +  # add points for individual values
    geom_bar(stat = "summary", fun = "mean", position = "dodge") +  # bar plot for mean
    geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.8), width = 0.25) +  # error bars
    ggforce::facet_row(~Gene, scales = "free_y") +  # facet by gene
    stat_compare_means(method = "t.test") +  # statistical significance
    ylab("TPM")+
    theme_minimal() +  # minimal theme for clarity
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15),  # x-axis text size
      axis.text.y = element_text(angle = 0, hjust = 1, size = 15),  # y-axis text size
      axis.title.y = element_text(angle = 90, size = 15),  # y-axis title size
      legend.title = element_blank(),  # remove legend title
      strip.text = element_text(size = 18),  # increase size of facet titles
      strip.background = element_blank()  # removes facet background
    )
}

# Function to process RNA data
process_rna_data <- function(file_path, genes_of_interest) {
  rna_matrix <- read.csv(file_path, header = TRUE)
  rownames(rna_matrix) <- make.names(rna_matrix$X, unique = TRUE)
  rna_matrix$X <- NULL
  
  # Filter out genes with total expression <= 10 across all samples
  rna_matrix_filtered <- rna_matrix[rowSums(rna_matrix) > 10, ]
  
  # Check which genes are available in the filtered matrix
  data <- rna_matrix_filtered[rownames(rna_matrix_filtered) %in% genes_of_interest, ] %>% 
    t() %>% 
    as.data.frame()
  
  data$Type <- str_split(rownames(data), pattern = "0uM|1uM|0nM", simplify = TRUE)[, 1]
  data$Concentration <- str_split(rownames(data), pattern = "ESCs|EpiLCs|_REP1|_REP2|_REP3", simplify = TRUE)[, 2]
  
  data_long <- pivot_longer(data, cols = genes_of_interest, values_to = "Expression", names_to = "Gene")
  data_long$Concentration <- ifelse(data_long$Concentration == "0nM", "0uM", data_long$Concentration)
  data_long$Sample <- paste0(data_long$Type, " ", data_long$Concentration)
  
  list(
    ESCs = filter(data_long, Type == "ESCs"),
    EpiLCs = filter(data_long, Type == "EpiLCs")
  )
}


# Process and plot RNA data
genes_of_interest_rna <- c("Dnmt3b", "Gapdh")
rna_data <- process_rna_data('RNA_CMF.csv', genes_of_interest_rna)

create_plot(rna_data$ESCs, "ESCs")
create_plot(rna_data$EpiLCs, "EpiLCs")




# Process and plot metabolomics data
# Define a function for plotting
create_plot_m <- function(data, type) {
  ggplot(data, aes(x = Sample, y = Expression, fill = Concentration)) + 
    geom_point() +  # add points for individual values
    geom_bar(stat = "summary", fun = "mean", position = "dodge") +  # bar plot for mean
    geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.8), width = 0.25) +  # error bars
    ggforce::facet_row(~Gene, scales = "free_y") +  # facet by gene
    stat_compare_means(method = "t.test") +  # statistical significance
    theme_minimal() +  # minimal theme for clarity
    ylab("Relative amount")+
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15),  # x-axis text size
      axis.text.y = element_text(angle = 0, hjust = 1, size = 15),  # y-axis text size
      axis.title.y = element_text(angle = 90, size = 15),  # y-axis title size
      legend.title = element_blank(),  # remove legend title
      strip.text = element_text(size = 18),  # increase size of facet titles
      strip.background = element_blank()  # removes facet background
    )
}

# Function to process RNA data
process_mtb_data <- function(file_path, genes_of_interest_metabo) {
  mtb_matrix <- read.csv(file_path, header = TRUE, row.names = 1)

  # Check which genes are available in the filtered matrix
  data <- mtb_matrix[rownames(mtb_matrix) %in% genes_of_interest_metabo, ] %>% t() %>% as.data.frame()
  
  data$Type <- str_split(rownames(data), pattern = "0uM|1uM|0nM", simplify = TRUE)[, 1]
  data$Concentration <- str_split(rownames(data), pattern = "ESCs|EpiLCs|_REP1|_REP2|_REP3", simplify = TRUE)[, 2]
  
  data_long <- pivot_longer(data, cols = genes_of_interest_metabo, values_to = "Expression", names_to = "Gene")
  data_long$Concentration <- ifelse(data_long$Concentration == "0nM", "0uM", data_long$Concentration)
  data_long$Sample <- paste0(data_long$Type, " ", data_long$Concentration)
  
  list(
    ESCs = filter(data_long, Type == "ESCs"),
    EpiLCs = filter(data_long, Type == "EpiLCs")
  )
}


genes_of_interest_metabo <- c("Malate")
metabo_data <- process_mtb_data(file_path = 'METABO_CMF.csv', genes_of_interest_metabo)

create_plot_m(metabo_data$ESCs, "ESCs")
create_plot_m(metabo_data$EpiLCs, "EpiLCs")
