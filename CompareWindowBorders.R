# Create logical matrix
logical_matrix <- compare.window

# Function to find border coordinates
find_border_coords <- function(mat) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  # Initialize a matrix to store border coordinates
  border_coords <- matrix(nrow = 0, ncol = 2)
  
  # Check horizontal borders
  for (i in 1:nr) {
    for (j in 1:(nc-1)) {
      if (mat[i, j] != mat[i, j + 1]) {
        border_coords <- rbind(border_coords, c(i, j))
        border_coords <- rbind(border_coords, c(i, j + 1))
      }
    }
  }
  
  # Check vertical borders
  for (j in 1:nc) {
    for (i in 1:(nr-1)) {
      if (mat[i, j] != mat[i + 1, j]) {
        border_coords <- rbind(border_coords, c(i, j))
        border_coords <- rbind(border_coords, c(i + 1, j))
      }
    }
  }
  
  return(border_coords)
}

# Find the border coordinates
border_coords <- find_border_coords(logical_matrix)
