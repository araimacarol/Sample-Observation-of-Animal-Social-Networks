setwd("C:/Users/rcappaw/Desktop/R/Workflow/igraphEpi-New")

#| echo: false
# L/oad required libraries
options("install.lock"=FALSE)
#library(lightgbm)
library(bonsai)
#library(treesnip)
library(kernelshap)
library(future)
library(future.apply)
library(shapviz)
library(MASS)
library(factoextra)
library(igraph)
library(dplyr)
library(tidyr)
library(modelr)
library(janitor)
library(tidymodels)
library(ggplot2)
library(cowplot)
library(furrr)
library(themis)



# ##-----species----##
# #aves-wildbird-network
# #insecta-ant-trophallaxis-colony1
# #insecta-beetle-group-c
# #insecta-beetle-group-t
# #mammalia-baboon-grooming
# #mammalia-baboon-association
# #mammalia-primate-association
# #mammalia-raccoon-proximity
# #mammalia-voles-bhp-trapping
# #mammalia-voles-plj-trapping
# #mammalia-voles-kcs-trapping
# #mammalia-voles-rob-trapping
# #reptilia-tortoise-network-pv
# #reptilia-tortoise-network-bsv
# #reptilia-tortoise-network-fi
# #reptilia-tortoise-network-cs
# #reptilia-tortoise-network-lm
# #Chimp
# #insecta-beetle
# #Badger
# #aves-wildbird-network
# #insecta-ant-colony1
# #aves-weaver-social



##########################################################################
# Comparing features and Boxplots 
##########################################################################

CompareFeats_and_combinedPlots=function(data="Racoon_Proximity.csv",
                      x_label="Days"
                     
){
 
  #######################################################
  #----compare features----### 
  ########################################################
    df=read.csv(data,header=T, sep=",", fill=T)
  colnames(df) <- gsub("\\.", " ", colnames(df))
  df <- na.omit(df)
  df_numeric <- as_tibble(df[, sapply(df, is.numeric)])
  
  cols_to_normalize <- c("Mean eccentricity",
                         "Mean path length", "Graph energy", 
                         "Modularity", "Diameter", "Betweenness centrality", 
                         "Transitivity", 
                         "Spectral radius", "Eigen centrality", 
                         "Degree centrality", "Mean degree", 
                         "Min cut", "Fiedler value", "Normalized Fiedler", 
                         "Closeness centrality", 
                         "Degree assortativity")
  
  # Apply scaling operation using mutate() and across()
  df_scaled <- df %>%
    mutate(across(cols_to_normalize, ~ ./Nodes))
  # Step 2: Replace original columns with scaled ones
  df[, cols_to_normalize] <- df_scaled[, cols_to_normalize]
  
  ###----Complete networks----###
  comp_taxa_df=head(df, 1)
  ###----Sample networks----###
  sample_taxa_df=df[-1,]
  
  sample_taxa_df=sample_taxa_df %>%
    mutate(Time_series = 1:nrow( sample_taxa_df))

  ###----Select columns needed for plots
  new_column_names=c('Normalized Fiedler','Fiedler value','Transitivity',
                     'Spectral radius','Modularity',
                     'Eigen centrality','Degree centrality',
                     'Degree assortativity','Time_series')

  sample_taxa_df1=sample_taxa_df%>%
    dplyr::select(new_column_names)

  ###----Transform data to long format
  sample_taxa_df_long <- sample_taxa_df1 %>%
    gather(key = "Graph_Features", value = "Value", -Time_series)

  sample_taxa_df_long$Graph_Features=factor(sample_taxa_df_long$Graph_Features)

  ###----Get the levels of PredictedClass
  class_levels <- levels(factor(sample_taxa_df$Predicted_Class))

  max_value <- max(sample_taxa_df_long$Value, na.rm = TRUE)

  # Set the y-axis position for the points just above the maximum value
  y_position <- max_value + 0.02

  ##############################################################################
  ###----Create the combined plot with dual axes
  # distinct_colors <- viridis_pal(option = "E")(9)
  # Create a vector of distinct colors
  distinct_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#A65628", "#F781BF","turquoise","yellow3")
  # Box and whisker plot for num_of_nodes
  p1=ggplot(sample_taxa_df_long,
            aes(x = Graph_Features, y = Value, color = Graph_Features))+
    geom_boxplot() +
    labs(#title = "Box and Whisker Plot",
      x = "Graph_Features",
      y = "Values")+
    scale_color_manual(values = distinct_colors) +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 16),
          text = element_text(size = 16),
          axis.text= element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())


  p2=ggplot() +
  geom_point(data = sample_taxa_df_long, aes(x = Time_series, y = Value, color = Graph_Features))+
  geom_line(data = sample_taxa_df_long, aes(x = Time_series, y = Value, color = Graph_Features), size = 1) +
  scale_color_manual(values = distinct_colors) +
    geom_point(data = sample_taxa_df,
               aes(x = Time_series, y = y_position, shape = Predicted_Class),
               size = 4) +
    labs(
      x = x_label,
      y = "Values") +
    #scale_color_discrete(name = "Graph Features") +
    scale_shape_manual(name = "Predicted_Class", values = setNames(seq_along(class_levels), class_levels)) +
    theme_classic()+
    theme(axis.text.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          text = element_text(size = 16),
          axis.text= element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())


  # Combine the plots side by side
  combined_plot <- plot_grid(p2, p1,  ncol = 2)#,labels = c("A)", "B)"),)
  res=list()
  res$Combined_Plots=combined_plot 
  res$Scaled_Features=df
  return(res)
}


x=CompareFeats_and_combinedPlots(data="Racoon_Proximity.csv",
                                        x_label="Days")

##############################################################################
# Distance calculation of closest feature
#############################################################################
df=x$Scaled_Features
y_new=df[,-c(1,2,3,4)]

# Calculate Euclidean distance between each row and the first row
distances <- apply(y_new, 1, function(row) sqrt(sum((row - y_new[1, ])^2)))

# Create a data frame with distances and corresponding row numbers
distances_df <- data.frame(RowNumber = 1:nrow(y_new), Distance = distances)

# Remove the first row as it corresponds to the complete network
distances_df <- distances_df[-1, ]

# Find the index of the row with the least distance
min_distance_index <- which.min(distances_df$Distance)

# Extract the row with the least distance
row_with_least_distance <- df[min_distance_index, ]
# 
# # Extract the last row representing the complete network
complete_network<- df[1,]

res=rbind(row_with_least_distance,complete_network)

#Calculate absolute differences between row_with_least_distance and complete_network_row
differences <- abs(row_with_least_distance[,-c(1,2,3,4)] - 
                     complete_network[,-c(1,2,3,4)])

# Find the column with the smallest absolute difference
closest_column <- names(differences)[which.min(differences)]
###################################################################
# END
###################################################################

