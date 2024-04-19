library(flowClust)
library(flowCore)
library(tcltk)

# Choose FCS files using a file dialog window
fcs_files <- tk_choose.files(multi = TRUE)

# Choose the directory to save the PDF plots
output_directory <- tk_choose.dir()

# Create the output directory if it doesn't exist
if (!file.exists(output_directory)) {
  dir.create(output_directory)
}

# Empty list to store results for each FCS file
results_list <- list()


# Fix threshold for unstained samples (CD14 threshold) - A10, B10 and C10 ONLY Ctrls here
target_patterns <- c('A10_A10', 'B10_B10', 'C10_C10')

# Function to check if any pattern is in a path
contains_any_pattern <- function(path, patterns) {
  # Check each pattern; return TRUE if any pattern matches
  any(sapply(patterns, function(pattern) {
    grepl(pattern, path, fixed = TRUE)
  }))
}
matches <- sapply(file_paths, contains_any_pattern, patterns = target_patterns)

# Filter the paths
unstained_files <- file_paths[matches]
# Calculate median of CD14 for unstained samples
cd14_thresholds <- c()
for (file in unstained_files){
  unstained_file <- read.FCS(file, truncate_max_range = FALSE)
  cd14_channel_threshold <- quantile(exprs(unstained_file)[, "V677/20-A"], c(0.95))
  cd14_channel_threshold <- cd14_channel_threshold[[1]]
  cd14_thresholds <- c(cd14_thresholds, cd14_channel_threshold)
}
cd14_channel_threshold <- median(cd14_thresholds)

failed_wells <- c()
# Loop over each selected FCS file
for (file_path in fcs_files) {
  
  #Skip unstained files - here with the prefix "A", "B" or "C"
  to_skip <- c("_B", "_A", "_C")
  matches <- sapply(to_skip, function(substring) {
    grepl(substring, file_path)
  })
  if(any(matches)){
    next()
  }
  
  # Load your FCS file
  your_data <- read.FCS(file_path)
  
  # Stage 1: Perform clustering analysis on FSC channels to gate out singlets, number of fixed clusters from one to two, 200 EM iterations
  res1 <- flowClust(
    your_data,
    varNames = c("FSC-A", "FSC-H"),
    K = 1:2,
    B = 200
  )
  
  # Calculate percentiles for adjusting x-axis limits
  percentiles <- quantile(exprs(your_data)[, "FSC-A"], c(0.01, 1))
  
  # Set the x-axis limits using the percentiles
  xlim <- c(percentiles[1], percentiles[2])
  
  # Plot Stage 1 results
  plot_file_name <- file.path(output_directory, paste0("Stage1_", basename(file_path), ".pdf"))
  pdf(plot_file_name)
  plot(res1[[1]],
       data = your_data,
       level = 0.8,
       z.cutoff = 0,
       xlim = xlim)
  dev.off()
  
  # Stage 2: Perform clustering on SSC and CD14 channel to identify CD14+ cell population, number of fixed clusters from one to five, 2000 EM iterations
  your_data2 <- your_data[your_data %in% res1, ]
  
  # Set threshold for gating on V677/20-A CD14 channel
  threshold <- cd14_channel_threshold  # Adjust threshold value as needed
  
  # Create a logical vector indicating whether the values in V677/20-A channel exceed the threshold
  gate_indices <- exprs(your_data2)[, "V677/20-A"] > threshold
  
  # Apply the gating by subsetting the data
  your_data2_gated <- your_data2[gate_indices, ]
  
  # Perform clustering on the gated data
  res2 <- flowClust(
    your_data2_gated,
    varNames = c("V677/20-A", "SSC-A"),  # Adjust the channel names accordingly
    K = 1:3,
    B = 2000
  )
  
  # Store results in the results_list
  results_list[[file_path]] <- list(
    res1 = res1,
    res2 = res2
  )
  
  # Show Stage 2 results
  print(criterion(res2, "BIC"))
  print(summary(res2[[3]]))
  print(ruleOutliers(res2[[3]]))
  
  
  #insert well_name here, strsplit of file_path? For example: well_name <- strsplit(strsplit(strsplit(file_path, '//')[[1]][2], '\\.')[[1]][1], '_')[[1]][3] 
  if (any(is.na(res2[[3]]@sigma) | is.infinite(res2[[3]]@sigma))) {
    print('failed!') 
    failed_wells <- c(failed_wells, well_name)
    next()
  }
  
  # Plot Stage 2 clustering scatter plot with adjusted x-axis range
  plot_file_name <- file.path(output_directory, paste0("Stage2_", basename(file_path), ".pdf"))
  pdf(plot_file_name)
  plot(res2[[3]],
       data = your_data2_gated,
       level = 0.8,
       z.cutoff = 0.8,
       xlim = c(min(exprs(your_data2_gated)[, "V677/20-A"]), max(exprs(your_data2_gated)[, "V677/20-A"])))
  dev.off()
}
