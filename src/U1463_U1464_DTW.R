install.packages(setdiff(c("dtw", "DescTools", "astrochron"), rownames(installed.packages())))
source("src/win.f.R")
source("src/RMSE_utils.R")

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Picard1 and U1464 datasets

Picard1 <- read.csv("data/PICARD 1.csv", header=TRUE, stringsAsFactors=FALSE)
Picard1=Picard1[c(83:7750),] # Eocene-Miocene Unconformity
head(Picard1)
plot(Picard1, type="l", xlim = c(150, 1300), ylim = c(0, 50))

U1464 <- read.csv("data/U1464-HSGR.csv", header=TRUE, stringsAsFactors=FALSE)

#Biostratigraphic Data
U1463_U1464_depth <- read.csv("data/U1463-U1464_Depth.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting attenuated signal
M1 = Gmean(U1464[c(1:551),2])
M2 = Gmean(U1464[c(535:900),2])
SD1 = Gsd(U1464[c(1:551),2])
SD2 = Gsd(U1464[c(535:900),2])
U1464[c(1:553),2]=(U1464[c(1:553),2]+(M2-M1))*(SD1/SD2)

head(U1464)
plot(U1464, type="l", xlim = c(0, 800), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Picard1_interpolated <- astrochron::linterp(Picard1, dt = 0.2, genplot = F)
U1464_interpolated <- astrochron::linterp(U1464, dt = 0.2, genplot = F)

# Scaling the data
Pmean = DescTools::Gmean(Picard1_interpolated$GR)
Pstd = DescTools::Gsd(Picard1_interpolated$GR)
Picard1_scaled = (Picard1_interpolated$GR - Pmean)/Pstd
Picard1_rescaled = data.frame(Picard1_interpolated$DEPT, Picard1_scaled)

Umean = DescTools::Gmean(U1464_interpolated$HSGR)
Ustd = DescTools::Gsd(U1464_interpolated$HSGR)
U1464_scaled = (U1464_interpolated$HSGR - Umean)/Ustd
U1464_rescaled = data.frame(U1464_interpolated$DEPTH_WMSF, U1464_scaled)

# Resampling the data using moving window statistics
Picard1_scaled = mwStats(Picard1_rescaled, cols = 2, win=3, ends = T)
Picard1_standardized = data.frame(Picard1_scaled$Center_win, Picard1_scaled$Average)

U1464_scaled = mwStats(U1464_rescaled, cols = 2, win=3, ends = T)
U1464_standardized = data.frame(U1464_scaled$Center_win, U1464_scaled$Average)

# Plotting the rescaled and resampled data
plot(Picard1_standardized, type="l", xlim = c(150, 1300), ylim = c(-20, 20), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
plot(U1464_standardized, type="l", xlim = c(0, 800), ylim = c(-20, 20), xlab = "U1464 Resampled Depth", ylab = "Normalized GR")

#### DTW with step pattern asymmetric but no window ####

# Perform dtw
system.time(al_U1464_p1_ap <- dtw(U1464_standardized$U1464_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetric, open.begin = F, open.end = T))
plot(al_U1464_p1_ap, "threeway")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth_ap = astrochron::tune(U1464_standardized, cbind(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(U1464_on_Picard1_depth_ap, col = "red")

#### DTW with step pattern asymmetricP05 but no window ####

# Perform dtw
system.time(al_U1464_p1_ap05 <- dtw(U1464_standardized$U1464_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP05, open.begin = F, open.end = T))
plot(al_U1464_p1_ap05, "threeway")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth_ap05 = astrochron::tune(U1464_standardized, cbind(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap05$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap05$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(U1464_on_Picard1_depth_ap05, col = "red")

#### DTW with step pattern asymmetricP1 but no window ####

# Perform dtw
system.time(al_U1464_p1_ap1 <- dtw::dtw(U1464_standardized$U1464_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1, open.begin = F, open.end = T))
plot(al_U1464_p1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth_ap1 = astrochron::tune(U1464_standardized, cbind(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap1$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(U1464_on_Picard1_depth_ap1, col = "red")

#### DTW with step pattern asymmetricP2 but no window ####

# Perform dtw
system.time(al_U1464_p1_ap2 <- dtw::dtw(U1464_standardized$U1464_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP2, open.begin = F, open.end = T))
plot(al_U1464_p1_ap2, "threeway")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth_ap2 = astrochron::tune(U1464_standardized, cbind(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap2$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap2$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(U1464_on_Picard1_depth_ap2, col = "red")

#### DTW with step pattern asymmetricP1 but Sakoe Chiba window ####

# Perform dtw
system.time(al_U1464_p1_ap3 <- dtw::dtw(U1464_standardized$U1464_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1, window.type = "sakoechiba", window.size = 2000, open.begin = F, open.end = T))
plot(al_U1464_p1_ap3, "threeway")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth_ap3 = astrochron::tune(U1464_standardized, cbind(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap1$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(U1464_on_Picard1_depth_ap3, col = "red")

#### DTW with stratigraphy-optimized step pattern asymmetricP1.1 but no knowledge-based window ####

# Perform dtw
system.time(al_U1464_p1_ap11 <- dtw::dtw(U1464_standardized$U1464_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = T))
plot(al_U1464_p1_ap11, "threeway")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth_ap11 = astrochron::tune(U1464_standardized, cbind(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap11$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap11$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(U1464_on_Picard1_depth_ap11, col = "red")

# DTW Distance
al_U1464_p1_ap11$normalizedDistance
al_U1464_p1_ap11$distance


#### DTW with stratigraphy-optimized step pattern asymmetricP1.1 and knowledge-based window ####

# create matrix for the knowledge-based window

compare.window <- matrix(data=TRUE,nrow=nrow(U1464_standardized),ncol=nrow(Picard1_standardized))
image(x=Picard1_standardized[,1],y=U1464_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- DescTools::Closest(190, Picard1_standardized[,1],which=TRUE)
base_1_y <- DescTools::Closest(50, U1464_standardized[,1],which=TRUE)

# Depth values for second datum
base_2_x <- DescTools::Closest(270, Picard1_standardized[,1],which=TRUE)
base_2_y <- DescTools::Closest(120, U1464_standardized[,1],which=TRUE)

# Depth values for third datum
base_3_x <- DescTools::Closest(310, Picard1_standardized[,1],which=TRUE)
base_3_y <- DescTools::Closest(184, U1464_standardized[,1],which=TRUE)

# Depth values for fourth datum
base_4_x <- DescTools::Closest(390, Picard1_standardized[,1],which=TRUE)
base_4_y <- DescTools::Closest(275, U1464_standardized[,1],which=TRUE)

# Depth values for fifth datum
base_5_x <- DescTools::Closest(1010, Picard1_standardized[,1],which=TRUE)
base_5_y <- DescTools::Closest(700, U1464_standardized[,1],which=TRUE)

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(U1464_standardized), ncol = nrow(Picard1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+200):nrow(U1464_standardized),1:(base_1_x-200)] <- 0
compare.window[1:(base_1_y-200),(base_1_x+200):ncol(compare.window)] <- 0

compare.window[(base_2_y+400):nrow(U1464_standardized),1:(base_2_x-400)] <- 0
compare.window[1:(base_2_y-200),(base_2_x+200):ncol(compare.window)] <- 0

compare.window[(base_3_y+300):nrow(U1464_standardized),1:(base_3_x-300)] <- 0
compare.window[1:(base_3_y-400),(base_3_x+400):ncol(compare.window)] <- 0

compare.window[(base_4_y+300):nrow(U1464_standardized),1:(base_4_x-300)] <- 0
compare.window[1:(base_4_y-500),(base_4_x+500):ncol(compare.window)] <- 0

compare.window[(base_5_y+100):nrow(U1464_standardized),1:(base_5_x-100)] <- 0
compare.window[1:(base_5_y-100),(base_5_x+100):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Picard1_standardized[,1],y=U1464_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Picard1_standardized[,1],y=U1464_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a knowledge-based window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with knowledge-based window
system.time(al_U1464_p1_ap12 <- dtw::dtw(U1464_standardized$U1464_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_U1464_p1_ap12, type = "threeway")

# DTW Distance measure
al_U1464_p1_ap12$normalizedDistance
al_U1464_p1_ap12$distance

# Dtw knowledge-based window plot

image(y = Picard1_standardized[,1], x = U1464_standardized[,1], z = compare.window, useRaster = T)
lines(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap12$index1], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap12$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth_ap12 = astrochron::tune(U1464_standardized, cbind(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap12$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap12$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (U1464)")
lines(U1464_on_Picard1_depth_ap12, col = "red")

# Retrieve the corresponding depth values for step pattern asymmetric
Picard1_DTW_Depth_ap <- Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap$index2]
U1464_DTW_Depth_ap <- U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap$index1]

# Create a data frame to show the depths of Picard 1 and U1464 side by side
Picard1_U1464_Depth_ap <- data.frame(Picard1_Depth = Picard1_DTW_Depth_ap,  U1464_Depth = U1464_DTW_Depth_ap)

plot(Picard1_U1464_Depth_ap, type = "l")

# Retrieve the corresponding depth values for step pattern asymmetricP05
Picard1_DTW_Depth_ap05 <- Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap05$index2]
U1464_DTW_Depth_ap05 <- U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap05$index1]

# Create a data frame to show the depths of Picard 1 and U1464 side by side
Picard1_U1464_Depth_ap05 <- data.frame(Picard1_Depth = Picard1_DTW_Depth_ap05,  U1464_Depth = U1464_DTW_Depth_ap05)

plot(Picard1_U1464_Depth_ap05, type = "l")

# Retrieve the corresponding depth values for step pattern asymmetricP1
Picard1_DTW_Depth_ap1 <- Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap1$index2]
U1464_DTW_Depth_ap1 <- U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap1$index1]

# Create a data frame to show the depths of Picard 1 and U1464 side by side
Picard1_U1464_Depth_ap1 <- data.frame(Picard1_Depth = Picard1_DTW_Depth_ap1,  U1464_Depth = U1464_DTW_Depth_ap1)

plot(Picard1_U1464_Depth_ap1, type = "l")

# Retrieve the corresponding depth values for step pattern asymmetricP02
Picard1_DTW_Depth_ap2 <- Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap2$index2]
U1464_DTW_Depth_ap2 <- U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap2$index1]

# Create a data frame to show the depths of Picard 1 and U1464 side by side
Picard1_U1464_Depth_ap2 <- data.frame(Picard1_Depth = Picard1_DTW_Depth_ap2,  U1464_Depth = U1464_DTW_Depth_ap2)

plot(Picard1_U1464_Depth_ap2, type = "l")

# Retrieve the corresponding depth values for step pattern asymmetricP1 and Sakoe Chiba window
Picard1_DTW_Depth_ap3 <- Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap3$index2]
U1464_DTW_Depth_ap3 <- U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap3$index1]

# Create a data frame to show the depths of Picard 1 and U1464 side by side
Picard1_U1464_Depth_ap3 <- data.frame(Picard1_Depth = Picard1_DTW_Depth_ap3,  U1464_Depth = U1464_DTW_Depth_ap3)

plot(Picard1_U1464_Depth_ap3, type = "l")

# Retrieve the corresponding depth values for step pattern asymmetricP1.1 and no knowledge-based window
Picard1_DTW_Depth_ap11 <- Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap11$index2]
U1464_DTW_Depth_ap11 <- U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap11$index1]

# Create a data frame to show the depths of Picard 1 and U1464 side by side
Picard1_U1464_Depth_ap11 <- data.frame(Picard1_Depth = Picard1_DTW_Depth_ap11,  U1464_Depth = U1464_DTW_Depth_ap11)

plot(Picard1_U1464_Depth_ap11, type = "l")

# Retrieve the corresponding depth values for step pattern asymmetricP1.1 and knowledge-based window
Picard1_DTW_Depth_ap12 <- Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap12$index2]
U1464_DTW_Depth_ap12 <- U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap12$index1]

# Create a data frame to show the depths of Picard 1 and U1464 side by side
Picard1_U1464_Depth_ap12 <- data.frame(Picard1_Depth = Picard1_DTW_Depth_ap12,  U1464_Depth = U1464_DTW_Depth_ap12)

plot(Picard1_U1464_Depth_ap12, type = "l")


############################
# Calculate RMSE
###########################


find_closest_by_x <- function(obs_x, predicted_data) {
  diff_x <- abs(predicted_data$Picard1_Depth - obs_x)
  
  closest_index <- which.min(diff_x)
  
  return(predicted_data[closest_index, "U1464_Depth"])
}

calculate_rmse_for_line <- function(predicted_data, observed_data) {
  predicted_depths <- vector("numeric", nrow(observed_data))
  
  for (i in 1:nrow(observed_data)) {
    obs_x <- observed_data$Picard1_Depth[i]
    predicted_depths[i] <- find_closest_by_x(obs_x, predicted_data)
  }
  
  observed_depths <- observed_data$U1464_Depth
  rmse <- sqrt(mean((observed_depths - predicted_depths)^2))
  return(rmse)
}

predicted_lines <- list(
  Picard1_U1464_Depth_ap,
  Picard1_U1464_Depth_ap05,
  Picard1_U1464_Depth_ap1,
  Picard1_U1464_Depth_ap2,
  Picard1_U1464_Depth_ap3,
  Picard1_U1464_Depth_ap11,
  Picard1_U1464_Depth_ap12
)
rmse_values <- numeric(length(predicted_lines))

for (i in 1:length(predicted_lines)) {
  rmse_values[i] <- calculate_rmse_for_line(predicted_lines[[i]], U1463_U1464_depth)
}

best_line_index <- which.min(rmse_values)
best_rmse <- rmse_values[best_line_index]

print(paste("Best predicted line: Line", best_line_index))
print(paste("RMSE for the best predicted line:", best_rmse))
print(rmse_values)

############################
# Make Figure
###########################

setwd("Figures/")

pdf(file = "Figure 12.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 12.png", width = 5000, height = 6000, res = 600)

par(mar=c(5,5,1,1))
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), widths=c(1,1,1), heights=c(1.4,1,1.25))

plot(U1463_U1464_depth[,1], U1463_U1464_depth[,2], type ="n", ylim = c(800,0), xlim = c(150, 1050), xaxs = "i", yaxs = "i", axes = F, xlab = "", ylab = "")
#points(U1463_U1464_depth[,1], U1463_U1464_depth[,2], pch = 15, cex = 2, bg = "black")

segments(333,51.33,333,59.46, lwd = 2)
segments(333,51.33,333,60.84, lwd = 2)
segments(391,126.3,391,136.87, lwd = 2)
segments(404,97.42,404,126.3, lwd = 2)
segments(440,126.3,440,163, lwd = 2)
segments(426,173.23,426,189.65, lwd = 2)
segments(463,189.65,463,228.35, lwd = 2)
segments(493,228.35,493,281.39, lwd = 2)
segments(545,293.15,545,312.75, lwd = 2)
segments(1010,579.91,1010,707.57, lwd = 2)


segments(328.24,55.43,338.05,55.43, lwd = 2)
segments(328.24,55.63,338.05,55.63, lwd = 2)
segments(384.89,131.6,397.26,131.6, lwd = 2)
segments(397.26,111.92,411.63,111.92, lwd = 2)
segments(434.61,144,447.75,144, lwd = 2)
segments(397.11,181.43,454.11,181.43, lwd = 2)
segments(447.7,209,480.71,209, lwd = 2)
segments(454.11,255,532.22,255, lwd = 2)
segments(532.22,302.95,558.44,302.95, lwd = 2)
segments(1010,643.74,1010,643.71, lwd = 2)

axis(1, at = c(1050, seq(150,1050,100)), cex.axis = 1.25)
axis(2, at = c(750, seq(0,750,150)), cex.axis = 1.25)
lines(Picard1_U1464_Depth_ap, type = "l", col = "brown", lwd = 2)
lines(Picard1_U1464_Depth_ap05, type = "l", col = "cyan", lwd = 2)
lines(Picard1_U1464_Depth_ap1, type = "l", col = "blue", lwd = 4)
lines(Picard1_U1464_Depth_ap2, type = "l", col = "violet", lwd = 2)
lines(Picard1_U1464_Depth_ap3, type = "l", col = "red", lwd = 2)
lines(Picard1_U1464_Depth_ap11, type = "l", col = "darkorange", lwd = 2)
lines(Picard1_U1464_Depth_ap12, type = "l", col = "green", lwd = 2)
mtext("Picard-1 Depth (m)", side = 1, line = 2.5, cex = 0.9)
mtext("U1464 Depth (m)", side = 2, line = 2.5, cex = 0.9)
text(175,750,"(a)", cex = 1.5)

par(mar=c(0,5,0,1))

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "Natural Gamma Radiation (gAPI)", cex.lab = 1.25)
lines(U1464_on_Picard1_depth_ap1, col = "blue")
text(500,-18,"(b) Correlation with asymmetricP1 step pattern and no window", cex = 1.5)
axis(2, at = c(-20,0,20), cex.axis = 1.25, las = 1)
legend(x = max(Picard1_standardized$Picard1_scaled.Center_win)-250, y = max(Picard1_standardized$Picard1_scaled.Average)+8, legend = c("Picard-1", "U1464"), col = c("black", "blue"), lty = 1, lwd = 3, cex = 1.5, bty = "n")
arrows(x0=640, y0 =12 ,x1 = 590, y1 = 8, length = 0.10, angle = 30, code = 2, lwd = 1.25)
arrows(x0=710, y0 =12 ,x1 = 710, y1 = -4, length = 0.10, angle = 30, code = 2, lwd = 1.25)
text(x = 670, y = 14, "Over-stretching", cex = 1.25)

par(mar=c(5,5,1,1))

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, yaxt = "n", 
     xlab = "Picard-1 Depth (meters)", ylab = "Natural Gamma Radiation (gAPI)", cex.lab = 1.25)
lines(U1464_on_Picard1_depth_ap11, col = "darkorange")
axis(1, at = c(150, seq(300,1300,100)), cex.axis = 1.25, las = 1)
axis(2, at = c(-20,0,20), cex.axis = 1.25, las = 1)
text(540,-18,"(c) Correlation with stratigraphy-optimized step pattern and no window", cex = 1.5)
legend(x = max(Picard1_standardized$Picard1_scaled.Center_win)-250, y = max(Picard1_standardized$Picard1_scaled.Average)+8, legend = c("Picard-1", "U1464"), col = c("black", "darkorange"), lty = 1, lwd = 3, cex = 1.5, bty = "n")

dev.off()



############################
# Make Supplementary Figure
###########################

setwd("Figures/")

png(filename = "Figure 14.png", width = 6000, height = 8000, res = 600)

par(mar=c(0,5,1,1))
layout(matrix(c(1,2,3,4,5), 5, 1, byrow = TRUE), widths=c(1,1,1,1,1), heights=c(1,1,1,1,1.3))

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1464_on_Picard1_depth_ap, col = "brown")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
legend(x = max(Picard1$DEPT)-450, y = max(Picard1$GR)+2, legend = c("Picard-1", "U1464"), col = c("black", "orange"), lty = 1, lwd = 3, cex = 1.70, bty = "n")
text(270,15,"(a) asymmetric", cex = 1.5)

par(mar=c(0,5,0,1))

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim =  c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1464_on_Picard1_depth_ap05, col = "cyan")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(280,15,"(b) Asymmetric P05", cex = 1.5)

par(mar=c(0,5,0,1))

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1464_on_Picard1_depth_ap2, col = "violet")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(280,15,"(c) Asymmetric P2", cex = 1.5)

par(mar=c(0,5,0,1))

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1464_on_Picard1_depth_ap3, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(380,15,"(d) Asymmetric P1 + SakoeChiba window", cex = 1.5)

par(mar=c(5,5,0,1))

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1464_on_Picard1_depth_ap12, col = "green")
axis(1, at = c(150, seq(300,1300,100)), cex.axis = 1.5, las = 1)
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(480,15,"(e) Stratigraphy-optimized step pattern + Knowledge-based window", cex = 1.5)

dev.off()
