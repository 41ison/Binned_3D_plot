# Binned 3D plot for LC-MS/MS data

The following function will extract binned retention time (min), precursor m/z, and intensity values ​​from the DIA-NN report and return a list containing the information needed to create a 3D plot.
The reason I binned the data is because I couldn't find a workaround to plot many values ​​using R libraries like `graphics` or `GA`. If you know a way to do it, please let me know.

you will need the tidyverse library appended.
❗I'm assuming you know how to import report.parquet into your R session and do all the necessary data wrangling. If you have trouble doing this step, you can check out the repository for an in-depth analysis of the report.parquet file [here](https://github.com/41ison/QC4DIANN).

```r
library(tidyverse)
```

# Function to extract binned data from DIA-NN report for 3D plot
The function binning_data(data, num_RT_bins, num_Mz_bins) has 3 arguments:
1. data: your DIA-NN report file
2. num_RT_bins: the number of bins used to group retention time units
3. num_Mz_bins: the number of bins used to group precursor m/z values

```r
binning_data <- function(data, num_RT_bins, num_Mz_bins) {
z_matrix <- data %>%
    dplyr::select(RT, Precursor.Mz, Precursor.Normalised) %>%
    dplyr::mutate(
      RT_bin = cut(RT, breaks = num_RT_bins, labels = FALSE),  # Convert RT into bins
      Mz_bin = cut(Precursor.Mz, breaks = num_Mz_bins, labels = FALSE)  # Convert Mz into bins
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
      Mz_bin = cut(Precursor.Mz, breaks = num_Mz_bins, labels = FALSE)  # Convert Mz into bins
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

# This will return a list containing the z_matrix, x_vals, and y_vals that will be used to plot the 3D diagram using the graphics::persp() function
  return(list(z_matrix = z_matrix, x_vals = x_vals, y_vals = y_vals))
}
```

Once the function is loaded, all you need to do is import the DIA-NN report and select the bins that work for you.
Here, you need to find the optimal binning for your data. You can try starting with 50, then increasing by +10 until the function stops working.

```r
sample_binned <- binning_data(diann_report %>%
                                  dplyr::filter(Run == "HeLa_01"), 80, 80)
```

The 3D plot can be extracted using the 3 items in the sample_binned list.

```r
png("3D_plot.png", width = 10, height = 10, units = "in", res = 300)
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
```
<p align="center">
<img src="https://github.com/41ison/Binned_3D_plot/blob/main/3D_plot.png" width="500">
</p>
