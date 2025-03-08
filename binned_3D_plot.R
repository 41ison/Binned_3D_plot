# plot the 3D diagram for intensity, retention time, and m/z values
# download this script and place it in the same directory as the DIA-NN report file

library(tidyverse) # for data wrangling
library(here) # for file path management
library(arrow) # to read the report.parquet file

# load the DIA-NN report file
diann_report <- arrow::read_parquet("report.parquet")

{
# Function to extract binned data from DIA-NN report for 3D plot
binning_data <- function(data, num_RT_bins, num_Mz_bins) {
  z_matrix <- data %>%
    dplyr::select(RT, Precursor.Mz, Precursor.Normalised) %>%
    dplyr::mutate(
      RT_bin = cut(RT, breaks = num_RT_bins, labels = FALSE),  # Convert RT into bins
      Mz_bin = cut(Precursor.Mz, breaks = num_Mz_bins, labels = FALSE)  # Convert m/z values into bins
    ) %>%
    group_by(RT_bin, Mz_bin) %>%
    summarize(Precursor.Normalised = mean(Precursor.Normalised, na.rm = TRUE), .groups = "drop") %>% 
    pivot_wider(names_from = Mz_bin, values_from = Precursor.Normalised, values_fill = 0) %>%
    dplyr::select(-RT_bin) %>%
    as.matrix()
  
  # Extract the vector containing the unique values for x axis (retention time)
  x_vals <- data %>%
    dplyr::select(RT, Precursor.Mz, Precursor.Normalised) %>%
    dplyr::mutate(
      RT_bin = cut(RT, breaks = num_RT_bins, labels = FALSE),  # Convert RT into bins
      Mz_bin = cut(Precursor.Mz, breaks = num_Mz_bins, labels = FALSE)  # Convert m/z values into bins
    ) %>%
    group_by(RT_bin, Mz_bin) %>%
    summarize(Precursor.Normalised = mean(Precursor.Normalised, na.rm = TRUE), .groups = "drop") %>% 
    pull(RT_bin) %>%
    unique()
  
  # Extract the vector containing the unique values for y axis (precursor m/z)
  y_vals <- data %>% 
    dplyr::select(RT, Precursor.Mz, Precursor.Normalised) %>%
    dplyr::mutate(
      RT_bin = cut(RT, breaks = num_RT_bins, labels = FALSE),  # Convert RT into bins
      Mz_bin = cut(Precursor.Mz, breaks = num_Mz_bins, labels = FALSE)  # Convert m/z values into bins
    ) %>%
    group_by(RT_bin, Mz_bin) %>%
    summarize(Precursor.Normalised = mean(Precursor.Normalised, na.rm = TRUE), .groups = "drop") %>% 
    pull(Mz_bin) %>%
    unique()
  
  # this function will return a list containing the z_matrix, x_vals, and y_vals that will be used to plot the 3D diagram using the graphics::persp() function
  return(list(z_matrix = z_matrix, x_vals = x_vals, y_vals = y_vals))
}

}

# filter the sample data for plotting
sample_binned <- binning_data(diann_report %>% 
                                dplyr::filter(Run == "HeLa_01"), 50, 50)

# save the 3D plot as a png file
png("plot/3D_plot.png", width = 10, height = 10, units = "in", res = 300)
persp(
  x = sample_binned$x_vals, 
  y = sample_binned$y_vals, 
  z = sample_binned$z_matrix, 
  col = "forestgreen", 
  theta = 60, phi = 15, expand = 0.4,
  main = "LC-MS/MS dimensions",
  xlab = "Retention time (min)",
  ylab = "Precursor m/z",
  zlab = "Precursor intensity",
  scale = TRUE,
  ticktype = "simple",
  box = TRUE, axes = TRUE, nticks = 4,
  border = NULL,
  shade = 0.1
)
dev.off()