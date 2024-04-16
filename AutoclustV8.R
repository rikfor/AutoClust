##########Automated Cell clustering##########
#Aman Mebrahtu 



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

# Loop over each selected FCS file
for (file_path in fcs_files) {
  # Load your FCS file
  your_data <- read.FCS(file_path)
  
  # Stage 1: Perform clustering analysis on FSC channels to gate out singlets, nummber of fixed clusters from one to two, 200 EM iterations
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
  
  # Stage 2: Perform clustering on SSC and CD14 channel to identify CD14+ cell popluation, number of fixed clusters from one to five, 2000 EM iterations
  your_data2 <- your_data[your_data %in% res1, ]
  res2 <- flowClust(
    your_data2,
    varNames = c("V677/20-A", "SSC-A"),  # Adjust the channel names accordingly
    K = 1:5,
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
  
  # Plot Stage 2 scatter with adjusted x-axis range
  plot_file_name <- file.path(output_directory, paste0("Stage2_", basename(file_path), ".pdf"))
  pdf(plot_file_name)
  plot(res2[[3]],
       data = your_data2,
       level = 0.8,
       z.cutoff = 0.8,
       xlim = c(min(exprs(your_data2)[, "V677/20-A"]), max(exprs(your_data2)[, "V677/20-A"])))
  dev.off()
}
