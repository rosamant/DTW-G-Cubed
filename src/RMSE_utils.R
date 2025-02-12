find_closest_by_x <- function(obs_x, predicted_data) {
  diff_x <- abs(predicted_data$U1464_Depth - obs_x)
  
  closest_index <- which.min(diff_x)
  
  return(predicted_data[closest_index, "U1482_Depth"])
}

calculate_rmse_for_line <- function(predicted_data, observed_data) {
  predicted_depths <- vector("numeric", nrow(observed_data))
  
  for (i in 1:nrow(observed_data)) {
    obs_x <- observed_data$U1464_Depth[i]
    predicted_depths[i] <- find_closest_by_x(obs_x, predicted_data)
  }
  
  observed_depths <- observed_data$U1482_Depth
  rmse <- sqrt(mean((observed_depths - predicted_depths)^2))
  return(rmse)
}