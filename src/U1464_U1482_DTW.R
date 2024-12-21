install.packages(setdiff(c("dtw", "DescTools", "astrochron"), rownames(installed.packages())))
source("src/win.f.R")
source("src/RMSE_utils.R")

rm(list = ls())

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import U1464 and U1482 datasets

U1464 <- read.csv("data/U1464-HSGR.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting attenuated signal
M1 = Gmean(U1464[c(1:551),2])
M2 = Gmean(U1464[c(535:900),2])
SD1 = Gsd(U1464[c(1:551),2])
SD2 = Gsd(U1464[c(535:900),2])
U1464[c(1:553),2]=(U1464[c(1:553),2]+(M2-M1))*(SD1/SD2)

U1464_1 = U1464[c(1:2953),]
head(U1464_1)
plot(U1464_1, type="l", xlim = c(0, 500), ylim = c(0, 60))

U1482 <- read.csv("data/U1482.csv", header=TRUE, stringsAsFactors=FALSE)
head(U1482)
plot(U1482, type="l", xlim = c(0, 550), ylim = c(0, 40))

# Biostratigraphic Data
U1464_U1482_depth <- read.csv("data/U1464-U1482_Depth.csv", header=TRUE, stringsAsFactors=FALSE)

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets

U1464_1_interpolated <- astrochron::linterp(U1464_1, dt = 0.2, genplot = F)
U1482_interpolated <- astrochron::linterp(U1482, dt = 0.2, genplot = F)

# Scaling the data
Smean = DescTools::Gmean(U1464_1_interpolated$HSGR)
Sstd = DescTools::Gsd(U1464_1_interpolated$HSGR)
U1464_1_scaled = (U1464_1_interpolated$HSGR - Smean)/Sstd
U1464_1_rescaled = data.frame(U1464_1_interpolated$DEPTH_WMSF, U1464_1_scaled)

Umean = DescTools::Gmean(U1482_interpolated$GR)
Ustd = DescTools::Gsd(U1482_interpolated$GR)
U1482_scaled = (U1482_interpolated$GR - Umean)/Ustd
U1482_rescaled = data.frame(U1482_interpolated$DEPT, U1482_scaled)

# Resampling the data using moving window statistics
U1464_1_scaled = astrochron::wStats(U1464_1_rescaled, cols = 2, win=3, ends = T)
U1464_1_standardized = data.frame(U1464_1_scaled$Center_win, U1464_1_scaled$Average)

U1482_scaled = astrochron::mwStats(U1482_rescaled, cols = 2, win=3, ends = T)
U1482_standardized = data.frame(U1482_scaled$Center_win, U1482_scaled$Average)

# Plotting the rescaled and resampled data
plot(U1464_1_standardized, type="l", xlim = c(0, 500), ylim = c(-20, 20), xlab = "South Galapagos1 Resampled Depth", ylab = "Normalized GR")
plot(U1482_standardized, type="l", xlim = c(0, 550), ylim = c(-10, 15), xlab = "U1482 Resampled Depth", ylab = "Normalized GR")

#### DTW with step pattern asymmetric but no window ####

# Perform dtw
system.time(al_U1482_U1464_ap <- dtw(U1482_standardized$U1482_scaled.Average, U1464_1_standardized$U1464_1_scaled.Average, keep.internals = T, step.pattern = asymmetric, open.begin = F, open.end = T))
plot(al_U1482_U1464_ap, "threeway")

# Tuning the standardized data on reference depth scale
U1482_on_U1464_depth_ap = astrochron::tune(U1482_standardized, cbind(U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap$index1s], U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 500), xlab = "U1464 Resampled Depth", ylab = "Normalized GR")
lines(U1482_on_U1464_depth_ap, col = "red")

#### DTW with step pattern asymmetricP05 but no window ####

# Perform dtw
system.time(al_U1482_U1464_ap05 <- dtw(U1482_standardized$U1482_scaled.Average, U1464_1_standardized$U1464_1_scaled.Average, keep.internals = T, step.pattern = asymmetricP05, open.begin = F, open.end = T))
plot(al_U1482_U1464_ap05, "threeway")

# Tuning the standardized data on reference depth scale
U1482_on_U1464_depth_ap05 = astrochron::tune(U1482_standardized, cbind(U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap05$index1s], U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap05$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 500), xlab = "U1464 Resampled Depth", ylab = "Normalized GR")
lines(U1482_on_U1464_depth_ap05, col = "red")

#### DTW with step pattern asymmetricP1 but no window ####

# Perform dtw
system.time(al_U1482_U1464_ap1 <- dtw::dtw(U1482_standardized$U1482_scaled.Average, U1464_1_standardized$U1464_1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1, open.begin = F, open.end = T))
plot(al_U1482_U1464_ap1, "threeway")

# Tuning the standardized data on reference depth scale
U1482_on_U1464_depth_ap1 = astrochron::tune(U1482_standardized, cbind(U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap1$index1s], U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap1$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 500), xlab = "U1464 Resampled Depth", ylab = "Normalized GR")
lines(U1482_on_U1464_depth_ap1, col = "red")

#### DTW with step pattern asymmetricP2 but no window ####

# Perform dtw
system.time(al_U1482_U1464_ap2 <- dtw::dtw(U1482_standardized$U1482_scaled.Average, U1464_1_standardized$U1464_1_scaled.Average, keep.internals = T, step.pattern = asymmetricP2, open.begin = F, open.end = T))
plot(al_U1482_U1464_ap2, "threeway")

# Tuning the standardized data on reference depth scale
U1482_on_U1464_depth_ap2 = astrochron::tune(U1482_standardized, cbind(U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap2$index1s], U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap2$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 500), xlab = "U1464 Resampled Depth", ylab = "Normalized GR")
lines(U1482_on_U1464_depth_ap2, col = "red")

#### DTW with step pattern asymmetricP1 but Sakoe Chiba window ####

# Perform dtw
system.time(al_U1482_U1464_ap3 <- dtw::dtw(U1482_standardized$U1482_scaled.Average, U1464_1_standardized$U1464_1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1, window.type = "sakoechiba", window.size = 2000, open.begin = F, open.end = T))
plot(al_U1482_U1464_ap3, "threeway")

# Tuning the standardized data on reference depth scale
U1482_on_U1464_depth_ap3 = astrochron::tune(U1482_standardized, cbind(U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap1$index1s], U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap1$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 500), xlab = "U1464 Resampled Depth", ylab = "Normalized GR")
lines(U1482_on_U1464_depth_ap3, col = "red")

#### DTW with stratigraphy-optimized step pattern asymmetricP1.1 but no knowledge-based window ####

# Perform dtw
system.time(al_U1482_U1464_ap11 <- dtw::dtw(U1482_standardized$U1482_scaled.Average, U1464_1_standardized$U1464_1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = T))
plot(al_U1482_U1464_ap11, "threeway")

# Tuning the standardized data on reference depth scale
U1482_on_U1464_depth_ap11 = astrochron::tune(U1482_standardized, cbind(U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap11$index1s], U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap11$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 500), xlab = "U1464 Resampled Depth", ylab = "Normalized GR")
lines(U1482_on_U1464_depth_ap11, col = "red")

# DTW Distance
al_U1482_U1464_ap11$normalizedDistance
al_U1482_U1464_ap11$distance


#### DTW with stratigraphy-optimized step pattern asymmetricP1.1 and knowledge-based window ####

# create matrix for the knowledge-based window

compare.window <- matrix(data=TRUE,nrow=nrow(U1482_standardized),ncol=nrow(U1464_1_standardized))
image(x=U1464_1_standardized[,1],y=U1482_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- DescTools::Closest(52, U1464_1_standardized[,1],which=TRUE)
base_1_y <- DescTools::Closest(100, U1482_standardized[,1],which=TRUE)

# Depth values for second datum
base_2_x <- DescTools::Closest(128, U1464_1_standardized[,1],which=TRUE)
base_2_y <- DescTools::Closest(145, U1482_standardized[,1],which=TRUE)

# Depth values for third datum
base_3_x <- DescTools::Closest(190, U1464_1_standardized[,1],which=TRUE)
base_3_y <- DescTools::Closest(174, U1482_standardized[,1],which=TRUE)

# Depth values for fourth datum
base_4_x <- DescTools::Closest(282, U1464_1_standardized[,1],which=TRUE)
base_4_y <- DescTools::Closest(236, U1482_standardized[,1],which=TRUE)

# Depth values for fifth datum
base_5_x <- DescTools::Closest(313, U1464_1_standardized[,1],which=TRUE)
base_5_y <- DescTools::Closest(290, U1482_standardized[,1],which=TRUE)

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(U1482_standardized), ncol = nrow(U1464_1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+50):nrow(U1482_standardized),1:(base_1_x-50)] <- 0
compare.window[1:(base_1_y-50),(base_1_x+50):ncol(compare.window)] <- 0

compare.window[(base_2_y+50):nrow(U1482_standardized),1:(base_2_x-50)] <- 0
compare.window[1:(base_2_y-50),(base_2_x+50):ncol(compare.window)] <- 0

compare.window[(base_3_y+50):nrow(U1482_standardized),1:(base_3_x-50)] <- 0
compare.window[1:(base_3_y-50),(base_3_x+50):ncol(compare.window)] <- 0
 
compare.window[(base_4_y+50):nrow(U1482_standardized),1:(base_4_x-50)] <- 0
compare.window[1:(base_4_y-50),(base_4_x+50):ncol(compare.window)] <- 0

compare.window[(base_5_y+100):nrow(U1482_standardized),1:(base_5_x-100)] <- 0
compare.window[1:(base_5_y-100),(base_5_x+100):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=U1464_1_standardized[,1],y=U1482_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=U1464_1_standardized[,1],y=U1482_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Perform dtw with knowledge-based window
system.time(al_U1482_U1464_ap12 <- dtw::dtw(U1482_standardized$U1482_scaled.Average, U1464_1_standardized$U1464_1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_U1482_U1464_ap12, type = "threeway")

# DTW Distance measure
al_U1482_U1464_ap12$normalizedDistance
al_U1482_U1464_ap12$distance

image(y = U1464_1_standardized[,1], x = U1482_standardized[,1], z = compare.window, useRaster = T)
lines(U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap12$index1], U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap12$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
U1482_on_U1464_depth_ap12 = astrochron::tune(U1482_standardized, cbind(U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap12$index1s], U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap12$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 500), xlab = "U1464 Resampled Depth", ylab = "Normalized GR (U1482)")
lines(U1482_on_U1464_depth_ap12, col = "red")

# Retrieve the corresponding depth values for step pattern asymmetric
U1464_DTW_Depth_ap <- U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap$index2]
U1482_DTW_Depth_ap <- U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap$index1]

# Create a data frame to show the depths of U1464 and U1482 side by side
U1464_U1482_Depth_ap <- data.frame(U1464_Depth = U1464_DTW_Depth_ap,  U1482_Depth = U1482_DTW_Depth_ap)

plot(U1464_U1482_Depth_ap, type = "l")

# Retrieve the corresponding depth values for step pattern asymmetricP05
U1464_DTW_Depth_ap05 <- U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap05$index2]
U1482_DTW_Depth_ap05 <- U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap05$index1]

# Create a data frame to show the depths of U1464 and U1482 side by side
U1464_U1482_Depth_ap05 <- data.frame(U1464_Depth = U1464_DTW_Depth_ap05,  U1482_Depth = U1482_DTW_Depth_ap05)

plot(U1464_U1482_Depth_ap05, type = "l")

# Retrieve the corresponding depth values for step pattern asymmetricP1
U1464_DTW_Depth_ap1 <- U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap1$index2]
U1482_DTW_Depth_ap1 <- U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap1$index1]

# Create a data frame to show the depths of U1464 and U1482 side by side
U1464_U1482_Depth_ap1 <- data.frame(U1464_Depth = U1464_DTW_Depth_ap1,  U1482_Depth = U1482_DTW_Depth_ap1)

plot(U1464_U1482_Depth_ap1, type = "l")

# Retrieve the corresponding depth values for step pattern asymmetricP02
U1464_DTW_Depth_ap2 <- U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap2$index2]
U1482_DTW_Depth_ap2 <- U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap2$index1]

# Create a data frame to show the depths of U1464 and U1482 side by side
U1464_U1482_Depth_ap2 <- data.frame(U1464_Depth = U1464_DTW_Depth_ap2,  U1482_Depth = U1482_DTW_Depth_ap2)

plot(U1464_U1482_Depth_ap2, type = "l")

# Retrieve the corresponding depth values for step pattern asymmetricP1 and Sakoe Chiba window
U1464_DTW_Depth_ap3 <- U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap3$index2]
U1482_DTW_Depth_ap3 <- U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap3$index1]

# Create a data frame to show the depths of U1464 and U1482 side by side
U1464_U1482_Depth_ap3 <- data.frame(U1464_Depth = U1464_DTW_Depth_ap3,  U1482_Depth = U1482_DTW_Depth_ap3)

plot(U1464_U1482_Depth_ap3, type = "l")

# Retrieve the corresponding depth values for step pattern asymmetricP1.1 and no knowledge-based window
U1464_DTW_Depth_ap11 <- U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap11$index2]
U1482_DTW_Depth_ap11 <- U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap11$index1]

# Create a data frame to show the depths of U1464 and U1482 side by side
U1464_U1482_Depth_ap11 <- data.frame(U1464_Depth = U1464_DTW_Depth_ap11,  U1482_Depth = U1482_DTW_Depth_ap11)

plot(U1464_U1482_Depth_ap11, type = "l")

# Retrieve the corresponding depth values for step pattern asymmetricP1.1 and knowledge-based window
U1464_DTW_Depth_ap12 <- U1464_1_standardized$U1464_1_scaled.Center_win[al_U1482_U1464_ap12$index2]
U1482_DTW_Depth_ap12 <- U1482_standardized$U1482_scaled.Center_win[al_U1482_U1464_ap12$index1]

# Create a data frame to show the depths of U1464 and U1482 side by side
U1464_U1482_Depth_ap12 <- data.frame(U1464_Depth = U1464_DTW_Depth_ap12,  U1482_Depth = U1482_DTW_Depth_ap12)

plot(U1464_U1482_Depth_ap12, type = "l")


############################
# Calculate RMSE
###########################

predicted_lines1 <- list(
  U1464_U1482_Depth_ap,
  U1464_U1482_Depth_ap05,
  U1464_U1482_Depth_ap1,
  U1464_U1482_Depth_ap2,
  U1464_U1482_Depth_ap3,
  U1464_U1482_Depth_ap11,
  U1464_U1482_Depth_ap12
)

rmse_values1 <- numeric(length(predicted_lines))

for (i in 1:length(predicted_lines1)) {
  rmse_values1[i] <- calculate_rmse_for_line(predicted_lines1[[i]], U1464_U1482_depth)
}

best_line_index1 <- which.min(rmse_values1)
best_rmse1 <- rmse_values1[best_line_index1]

print(paste("Best predicted line: Line", best_line_index1))
print(paste("RMSE for the best predicted line:", best_rmse1))
print(rmse_values1)


############################
# Make Figure
###########################

setwd("Figures/")

pdf(file = "Figure 13.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 13.png", width = 5000, height = 6000, res = 600)

par(mar=c(5,5,1,1))
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), widths=c(1,1,1), heights=c(1.4,1,1.25))

plot(U1464_U1482_depth[,1], U1464_U1482_depth[,2], type ="n", ylim = c(300,0), xlim = c(0, 400), xaxs = "i", yaxs = "i", axes = F, xlab = "", ylab = "")
points(U1464_U1482_depth[,1], U1464_U1482_depth[,2], pch = 15, cex = 2, bg = "black")

axis(1, at = c(400, seq(0,400,50)), cex.axis = 1.25)
axis(2, at = c(300, seq(0,300,50)), cex.axis = 1.25)
lines(U1464_U1482_Depth_ap, type = "l", col = "brown", lwd = 2)
lines(U1464_U1482_Depth_ap05, type = "l", col = "cyan", lwd = 2)
lines(U1464_U1482_Depth_ap1, type = "l", col = "blue", lwd = 5)
lines(U1464_U1482_Depth_ap2, type = "l", col = "violet", lwd = 2)
lines(U1464_U1482_Depth_ap3, type = "l", col = "red", lwd = 2)
lines(U1464_U1482_Depth_ap11, type = "l", col = "darkorange", lwd = 2)
lines(U1464_U1482_Depth_ap12, type = "l", col = "green", lwd = 2)
mtext("U1464 depth (m)", side = 1, line = 2.5, cex = 0.9)
mtext("U1482 Depth (m)", side = 2, line = 2.5, cex = 0.9)
text(10,290,"(a)", cex = 1.5)

par(mar=c(0,5,0,1))

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 450), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "U1464 Depth (meters)", ylab = "Natural Gamma Radiation (gAPI)", cex.lab = 1.25)
lines(U1482_on_U1464_depth_ap11, col = "darkorange")
axis(2, at = c(-20,0,20), cex.axis = 1.25, las = 1)
text(150,-18,"(b) Correlation with stratigraphy-optimized step pattern and no window", cex = 1.5)
legend(x = max(U1464_1_standardized$U1464_1_scaled.Center_win)-100, y = max(U1464_1_standardized$U1464_1_scaled.Average)+7, legend = c("U1464", "U1482"), col = c("black", "darkorange"), lty = 1, lwd = 3, cex = 1.5, bty = "n")
arrows(x0=120, y0 =12 ,x1 = 155, y1 = 5, length = 0.10, angle = 30, code = 2, lwd = 1.25)
text(x = 100, y = 14, "Over-compression", cex = 1.25)
arrows(x0=300, y0 =-2 ,x1 = 320, y1 = 6, length = 0.10, angle = 30, code = 1, lwd = 1.25)
text(x = 340, y = 8, "Over-stretching", cex = 1.25)

par(mar=c(5,5,1,1))

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 450), axes = FALSE, yaxt = "n", 
     xlab = "U1464 Depth (meters)", ylab = "Natural Gamma Radiation (gAPI)", cex.lab = 1.25)
lines(U1482_on_U1464_depth_ap12, col = "green")
axis(1, at = c(450, seq(0,450,50)), cex.axis = 1.25, las = 1)
axis(2, at = c(-20,0,20), cex.axis = 1.25, las = 1)
text(185,-18,"(c) Correlation with stratigraphy-optimized step pattern and knowledge-based window", cex = 1.5)
legend(x = max(U1464_1_standardized$U1464_1_scaled.Center_win)-100, y = max(U1464_1_standardized$U1464_1_scaled.Average)+7, legend = c("U1464", "U1482"), col = c("black", "green"), lty = 1, lwd = 3, cex = 1.5, bty = "n")

dev.off()


############################
# Make Supplementary Figure
###########################

setwd("Figures/")

png(filename = "Figure 15.png", width = 6000, height = 8000, res = 600)

par(mar=c(0,5,1,1))
layout(matrix(c(1,2,3,4,5), 5, 1, byrow = TRUE), widths=c(1,1,1,1,1), heights=c(1,1,1,1,1.3))

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 450), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1482_on_U1464_depth_ap, col = "brown")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(40,15,"(a) asymmetric", cex = 1.5)

par(mar=c(0,5,0,1))

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim =  c(0, 450), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1482_on_U1464_depth_ap05, col = "cyan")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(50,15,"(b) Asymmetric P05", cex = 1.5)

par(mar=c(0,5,0,1))

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 450), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1482_on_U1464_depth_ap1, col = "blue")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(50,15,"(c) Asymmetric P1", cex = 1.5)

par(mar=c(0,5,0,1))

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 450), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1482_on_U1464_depth_ap2, col = "violet")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(50,15,"(d) Asymmetric P2", cex = 1.5)

par(mar=c(5,5,0,1))

plot(U1464_1_standardized, type = "l", ylim = c(-20, 20), xlim = c(0, 450), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "U1464 Depth (meters)", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1482_on_U1464_depth_ap3, col = "red")
axis(1, at = c(450, seq(0,450,50)), cex.axis = 1.5, las = 1)
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(75,15,"(e) Asymmetric P1 + SakoeChiba window", cex = 1.5)

dev.off()
