#https://github.com/RGLab/flowClust
# Article FlowClust : https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-145
# Article FlowCore: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-106

#The proposed approach is based on multivariate $t$ mixture models with the Box-Cox transformation.  
#This approach generalizes Gaussian mixture models by modeling outliers using $t$ distributions and allowing for clusters taking non-ellipsoidal 
#shapes upon proper data transformation.  Parameter estimation is carried out using an Expectation-Maximization (EM) algorithm which simultaneously handles 
# outlier identification and transformation selection.

#This **flowClust** package consists of a core function to implement the aforementioned clustering methodology.  
# Its source code is built in C for optimal utilization of system resources. 
#Currently, **flowClust** provides two options for estimating the number of clusters when it is unknown, namely, the Bayesian Information Criterion (BIC) 
# and the Integrated Completed Likelihood (ICL).

# **flowClust** is built in a way such that it is highly integrated with **flowCore**,the core package for flow cytometry that provides data structures and basic manipulation of flow cytometry data.  
# Please read https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-106  for details about actual implementation.


#####################
####### TODO ########
##### Fix so we automatically change channel names (maybe here: https://rdrr.io/bioc/flowWorkspace/man/cytoframe-labels.html) - no e.g. CD14 pERK are named weird

# [1] "B530/30-A"
# [1] "ILKBa-A488"
# [1] "-------------"
# [1] "V677/20-A"
# [1] "CD14-BV650"
# [1] "-------------"
# [1] "YG610/20-A"
# [1] "pERK-PE-Dazzle"
# [1] "-------------"
# [1] "YG780/60-A"
# [1] "p38-PE-Cy7"

# [1] "-------------"

library(flowClust)
library(flowCore)
library(Seurat)
library(dplyr)
library(tidyverse)

fcs_directory <- "/Users/rikardforlin/Forskning/Data/Hamilton/96Well-VbottompFlow_stimuli LPS_240406/"
# List all CSV files in the directory
file_paths <- list.files(path = fcs_directory, pattern = "\\.fcs$", full.names = TRUE)
files_to_look_at <- c(#'Specimen_001_G3_G03.fcs', 'Specimen_001_G1_G01.fcs', #'Specimen_001_E7_E07.fcs', 'Specimen_001_D1_D01.fcs', 'Specimen_001_A10_A10.fcs','Specimen_001_F2_F02.fcs', 'Specimen_001_B10_B10.fcs',
                      # 'Specimen_001_C10_C10.fcs',
                      'Specimen_001_F4_F04.fcs')
#flowfile <- read.FCS('/Users/rikardforlin/Forskning/Data/Hamilton/96Well-VbottompFlow_stimuli LPS_240406/Specimen_001_B6_B06.fcs', truncate_max_range = FALSE)
median_MFI_CD14_overall <- list()
median_MFI_PERK_overall <- list()
median_MFI_P38_overall <- list()
median_MFI_ILKBA_overall <- list()
#for (fcs_file in file_paths){
#for (i in 1:20){
for (fcs_file_tolookat in files_to_look_at){
  fcs_file_tolookat <- 'Specimen_001_F4_F04.fcs'
  well_name <- strsplit(strsplit(fcs_file_tolookat, '\\.')[[1]][1], '_')[[1]][3]
  
  #fcs_file <- '/Users/rikardforlin/Forskning/Data/Hamilton/96Well-VbottompFlow_stimuli LPS_240406/Specimen_001_E7_E07.fcs'#file_paths[[1]] #1
  fcs_file <- paste0(fcs_directory, fcs_file_tolookat)
  
  print(fcs_file)
  flowfile <- read.FCS(fcs_file, truncate_max_range = FALSE)
  
  #IDENTIFY OUTLIERS:
  #The following code performs an analysis with one cluster using the two scattering parameters:
  # The main purpose of performing an analysis with one cluster here is to identify outliers, which will be removed from subsequent analysis.
  res1 <-
    flowClust(
      flowfile,
      varNames = c("FSC-A", "SSC-A"),
      K = 1, #number of clusters
      B = 100, #maximum number of EM iterations
      min = c(0, 0), #Filtering, minimum values for varNames[1] and varNames[2]. 
      max = c(100000, 1000000)#Filtering, maximum values for varNames[1] and varNames[2]. 
    )
  #ruleOutliers(res1) <- list(level = 0.90)
  # Next, we would like to proceed with an analysis using the two fluorescence parameters on cells selected from the first stage.  
  # The following code performs the analysis with the number of clusters being fixed from one to six in turn
  
  flowfile2 <- flowfile[flowfile %in% res1, ]
  res2 <-
    flowClust(
      flowfile2,
      varNames = c("V677/20-A", "SSC-A"), #CD14
      K = 3:5,
      B = 100,
      min = c(0, 0), #Filtering, minimum values for varNames[1] and varNames[2]. 
      max = c(100000, 1000000)#Filtering, maximum values for varNames[1] and varNames[2]. 
    )

  #?flowClust
  # Assuming res is your list of models
  #logLikes <- numeric(length(res2))  # Create a numeric vector to store log likelihoods
  # Loop through each model in the list
  #for(i in seq_along(res2)) {
  #  logLikes[i] <- res2[[i]]@logLike  # Extract log likelihood
  #}
  #criterion(res2, 'BIC')
  #criterion(res2, 'ICL')
  #logLikes
  # Find the index of the model with the highest log likelihood
  #highest_LL <- which.max(logLikes)
  #res2 <- res2[[highest_LL]]
  #?plot
  # plot(rr[[i]],
  #     data = flowfile2,
  # #level = 0.8,
  # #z.cutoff = 0,
  # )
  
  expression_data <- as.data.frame(flowfile2@exprs)
  
  rr <- res2
  median_MFI_CD14_clusterlevel <- list()
  median_MFI_PERK_clusterlevel <- list()
  median_MFI_P38_clusterlevel <- list()
  median_MFI_ILKBA_clusterlevel <- list()
  #length(rr)
  
  for (i in 1:length(rr)){
    res2 <- rr[[i]]
    cluster_labels <- res2@label
    #Fill outliers (NAN values) with 0
    cluster_labels <- lapply(cluster_labels, function(x) { ifelse(is.na(x), 0, x) })
    #Add cluster belonging to expression data
    expression_data$Cluster <- list(cluster_labels)
    unique_clusters <- unique(cluster_labels)
    number_of_clusters <- length(unique_clusters) -1 #-1 to remove outliers as a cluster
    
    savestr <- paste0("/Users/rikardforlin/Forskning/Data/Hamilton/automated_cluster_figures/", well_name, '_', number_of_clusters, 'clusters', '.png')
  
    png(savestr)
    plot(res2,
       data = flowfile2,
       #level = 0.8,
       #z.cutoff = 0,
       main = paste0(well_name, '_', number_of_clusters, 'Clusters')
        )
    dev.off()
    
    # Initialize a list to store mean MFI results for each cluster
    median_MFI_CD14 <- list()
    median_MFI_PERK <- list()
    median_MFI_P38 <- list()
    median_MFI_ILKBA <- list()
    
    for(cluster in unique_clusters) {
      if (cluster == 0){
        next()
      }
      # Subset cells belonging to the current cluster
      cluster_cells <- expression_data[cluster_labels == cluster, ]
      
      # Calculate mean MFI for each marker in the subset
      median_MFI <- median(cluster_cells$`V677/20-A`) #CD14 #`YG610/20-A`) #pERK
      median_MFIPERK <- median(cluster_cells$`YG610/20-A`) #pERK
      median_MFIP38 <- median(cluster_cells$`YG780/60-A`) #p38
      median_MFIILKBA <- median(cluster_cells$`B530/30-A`) #p38
      
      # Store the result in the list
      median_MFI_CD14[[paste("Cluster", cluster)]] <- median_MFI
      median_MFI_PERK[[paste("Cluster", cluster)]] <- median_MFIPERK
      median_MFI_P38[[paste("Cluster", cluster)]] <- median_MFIP38
      median_MFI_ILKBA[[paste("Cluster", cluster)]] <- median_MFIILKBA
    }
    median_MFI_CD14_clusterlevel[[paste0(i+2, "_Cluster")]] <- median_MFI_CD14
    median_MFI_PERK_clusterlevel[[paste0(i+2, "_Cluster")]] <- median_MFI_PERK
    median_MFI_P38_clusterlevel[[paste0(i+2, "_Cluster")]] <- median_MFI_P38
    median_MFI_ILKBA_clusterlevel[[paste0(i+2, "_Cluster")]] <- median_MFI_ILKBA
    
    }
  
  median_MFI_CD14_overall[[strsplit(strsplit(fcs_file, "/")[[1]][8], ".fcs")[[1]][1]]] <- median_MFI_CD14_clusterlevel
  median_MFI_PERK_overall[[strsplit(strsplit(fcs_file, "/")[[1]][8], ".fcs")[[1]][1]]] <- median_MFI_PERK_clusterlevel
  median_MFI_P38_overall[[strsplit(strsplit(fcs_file, "/")[[1]][8], ".fcs")[[1]][1]]] <- median_MFI_P38_clusterlevel
  median_MFI_ILKBA_overall[[strsplit(strsplit(fcs_file, "/")[[1]][8], ".fcs")[[1]][1]]] <- median_MFI_ILKBA_clusterlevel
  
}



# # Initialize an empty list to store individual rows
# rows_list <- list()
# # Iterate over the list to create rows
# index <- 1
# for (group_name in names(median_MFI_CD14_overall)) {
#   for (item_name in names(median_MFI_CD14_overall[[group_name]])) {
#     # Current item data
#     item_data <- median_MFI_CD14_overall[[group_name]][[item_name]]
#     
#     # Create a row as a named vector (or list) and include group and item names
#     row_data <- c(Group = group_name, Item = item_name, item_data)
#     
#     # Add to the list of rows
#     rows_list[[index]] <- row_data
#     index <- index + 1
#   }
# }
# 
# # Convert the list of rows into a dataframe with bind_rows
# df_cd14 <- bind_rows(lapply(rows_list, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
# 
# # If needed, convert character columns that should be numeric
# numeric_cols <- sapply(df_cd14, function(x) all(!is.na(suppressWarnings(as.numeric(x)))))
# df_cd14[numeric_cols] <- lapply(df_cd14[numeric_cols], function(x) as.numeric(as.character(x)))
# 
# write.csv2(df_cd14, '/Users/rikardforlin/Forskning/Data/Hamilton/flowclust_tmp/CD14_temp.csv')
# 
# 
# 
# # Initialize an empty list to store individual rows
# rows_list <- list()
# # Iterate over the list to create rows
# index <- 1
# for (group_name in names(median_MFI_PERK_overall)) {
#   for (item_name in names(median_MFI_PERK_overall[[group_name]])) {
#     # Current item data
#     item_data <- median_MFI_PERK_overall[[group_name]][[item_name]]
#     
#     # Create a row as a named vector (or list) and include group and item names
#     row_data <- c(Group = group_name, Item = item_name, item_data)
#     
#     # Add to the list of rows
#     rows_list[[index]] <- row_data
#     index <- index + 1
#   }
# }
# 
# # Convert the list of rows into a dataframe with bind_rows
# df_perk <- bind_rows(lapply(rows_list, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
# 
# # If needed, convert character columns that should be numeric
# numeric_cols <- sapply(df_perk, function(x) all(!is.na(suppressWarnings(as.numeric(x)))))
# df_perk[numeric_cols] <- lapply(df_perk[numeric_cols], function(x) as.numeric(as.character(x)))
# 
# write.csv2(df_perk, '/Users/rikardforlin/Forskning/Data/Hamilton/flowclust_tmp/pERK_temp.csv')
# 
# 
# df <- merge(df_perk, df_cd14)
# df
# 



# Example nested_list modification for demonstration
# nested_list <- list(
#   "file1" = list(
#     "2_clusters" = list(Marker1 = 100, Marker2 = 200, Marker3 = 300),
#     "3_clusters" = list(Marker1 = 150, Marker2 = 250, Marker3 = 350, Marker4 = 450)
#   ),
#   "file2" = list(
#     "2_clusters" = list(Marker1 = 110, Marker2 = 210, Marker3 = 310),
#     "3_clusters" = list(Marker1 = 160, Marker2 = 260, Marker3 = 360, Marker4 = 460, Marker5 = 560)
#   )
# )




# Flatten the nested list into a dataframe while handling variable numbers of markers
df_p38 <- map_df(names(median_MFI_P38_overall), function(file_name) {
  map_df(names(median_MFI_P38_overall[[file_name]]), function(cluster_name) {
    markers <- median_MFI_P38_overall[[file_name]][[cluster_name]]
    tibble(
      FileName = file_name,
      Cluster = cluster_name,
      # Dynamically create columns for each marker
      !!!markers
    )
  })
})
write.csv2(df_p38, '/Users/rikardforlin/Forskning/Data/Hamilton/flowclust_tmp/dfp38_temp.csv')


df_cd14 <- map_df(names(median_MFI_CD14_overall), function(file_name) {
  map_df(names(median_MFI_CD14_overall[[file_name]]), function(cluster_name) {
    markers <- median_MFI_CD14_overall[[file_name]][[cluster_name]]
    tibble(
      FileName = file_name,
      Cluster = cluster_name,
      # Dynamically create columns for each marker
      !!!markers
    )
  })
})
write.csv2(df_cd14, '/Users/rikardforlin/Forskning/Data/Hamilton/flowclust_tmp/dfcd14_temp.csv')


df_perk <- map_df(names(median_MFI_PERK_overall), function(file_name) {
  map_df(names(median_MFI_PERK_overall[[file_name]]), function(cluster_name) {
    markers <- median_MFI_PERK_overall[[file_name]][[cluster_name]]
    tibble(
      FileName = file_name,
      Cluster = cluster_name,
      # Dynamically create columns for each marker
      !!!markers
    )
  })
})
write.csv2(df_perk, '/Users/rikardforlin/Forskning/Data/Hamilton/flowclust_tmp/dfperk_temp.csv')


df_ilkba <- map_df(names(median_MFI_ILKBA_overall), function(file_name) {
  map_df(names(median_MFI_ILKBA_overall[[file_name]]), function(cluster_name) {
    markers <- median_MFI_ILKBA_overall[[file_name]][[cluster_name]]
    tibble(
      FileName = file_name,
      Cluster = cluster_name,
      # Dynamically create columns for each marker
      !!!markers
    )
  })
})
write.csv2(df_ilkba, '/Users/rikardforlin/Forskning/Data/Hamilton/flowclust_tmp/dfilkba_temp.csv')


