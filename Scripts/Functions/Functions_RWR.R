
func_statisticsRWR <- function(visiting_prob_matrix) {
  # Initialize empty vectors to store IQR and number of outliers for each column
  iqr_values <- numeric(ncol(visiting_prob_matrix))
  num_outliers <- numeric(ncol(visiting_prob_matrix))
  
  # Loop through each column to compute IQR and number of outliers
  for (i in 1:ncol(visiting_prob_matrix)) {
    column_data <- visiting_prob_matrix[, i]
    
    # Calculate the quartiles
    lower_quartile <- quantile(column_data, 0.25)
    upper_quartile <- quantile(column_data, 0.75)
    
    # Calculate IQR
    iqr_values[i] <- IQR(column_data)
    
    # Calculate outliers: those outside 1.5 * IQR from the quartiles
    lower_bound <- lower_quartile - 1.5 * iqr_values[i]
    upper_bound <- upper_quartile + 1.5 * iqr_values[i]
    outliers <- column_data[column_data < lower_bound | column_data > upper_bound]
    
    # Store the number of outliers
    num_outliers[i] <- length(outliers)
  }
  
  # Print IQR values
  #print(paste("IQR values for each column: ", iqr_values))
  
  # Print number of outliers
  #print(paste("Number of outliers for each column: ", num_outliers))
  
  return(list(IQRs = iqr_values, OUTs = num_outliers))
}




# Function to identify the threshold for RWR probabilities using elbow method
# Find index of maximum distance
func_RWR_threshold <- function(probabilities){

  # Sort the input probabilities in descending order
  sorted_probabilities <- sort(probabilities, decreasing = TRUE)

  # Calculate distances to line connecting first and last points
  n <- length(sorted_probabilities)
  line_endpoint1 <- c(1, unname(sorted_probabilities[1]))
  line_endpoint2 <- c(n, unname(sorted_probabilities[n]))
  line_slope <- (line_endpoint2[2] - line_endpoint1[2]) / (line_endpoint2[1] - line_endpoint1[1])
  line_intercept <- line_endpoint1[2] - line_slope * line_endpoint1[1]

  # Calculate perpendicular distances from each point to the line
  distances <- abs((line_slope * (1:n) - sorted_probabilities + line_intercept) / sqrt(line_slope^2 + 1))

  # Find index of maximum distance, this is your elbow point
  elbow_point <- unname(which.max(distances))

  # To get the actual value at the elbow point
  elbow_value <- sorted_probabilities[elbow_point]

  # Count the number of points before the elbow point
  points_before_elbow <- elbow_point - 1

  return(list(ELB = unname(elbow_value), CNT = points_before_elbow))
}
