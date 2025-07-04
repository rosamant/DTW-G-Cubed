install.packages(setdiff(c("dtw", "DescTools", "astrochron"), rownames(installed.packages())))
source("src/win.f.R")

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Picard1 and Angel2 datasets

Picard1 <- read.csv("data/PICARD 1.csv", header=TRUE, stringsAsFactors=FALSE)
Picard1=Picard1[c(83:7750),] # Eocene-Miocene Unconformity
head(Picard1)
plot(Picard1, type="l", xlim = c(150, 1300), ylim = c(0, 50))

Angel2 <- read.csv("data/Angel2.csv", header=TRUE, stringsAsFactors=FALSE)
Angel2=Angel2[c(1:8419),] # Oligocene-Miocene
head(Angel2)
plot(Angel2, type="l", xlim = c(100, 1400), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Picard1_interpolated <- astrochron::linterp(Picard1, dt = 0.2, genplot = F)
Angel2_interpolated <- astrochron::linterp(Angel2, dt = 0.2, genplot = F)

# Scaling the data
Pmean = DescTools::Gmean(Picard1_interpolated$GR)
Pstd = DescTools::Gsd(Picard1_interpolated$GR)
Picard1_scaled = (Picard1_interpolated$GR - Pmean)/Pstd
Picard1_rescaled = data.frame(Picard1_interpolated$DEPT, Picard1_scaled)

Amean = DescTools::Gmean(Angel2_interpolated$GR)
Astd = DescTools::Gsd(Angel2_interpolated$GR)
Angel2_scaled = (Angel2_interpolated$GR - Amean)/Astd
Angel2_rescaled = data.frame(Angel2_interpolated$DEPT, Angel2_scaled)

# Resampling the data using moving window statistics
Picard1_scaled = mwStats(Picard1_rescaled, cols = 2, win=3, ends = T)
Picard1_standardized = data.frame(Picard1_scaled$Center_win, Picard1_scaled$Average)

Angel2_scaled = mwStats(Angel2_rescaled, cols = 2, win=3, ends = T)
Angel2_standardized = data.frame(Angel2_scaled$Center_win, Angel2_scaled$Average)

# Plotting the rescaled and resampled data
plot(Picard1_standardized, type="l", xlim = c(150, 1300), ylim = c(-20, 20), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
plot(Angel2_standardized, type="l", xlim = c(100, 1400), ylim = c(-20, 20), xlab = "Angel2 Resampled Depth", ylab = "Normalized GR")

#### DTW with step pattern asymmetricP1 and no knowledge-based step pattern or window  ####

# Perform dtw
system.time(al_a1_p1_p1 <- dtw::dtw(Angel2_standardized$Angel2_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1, open.begin = F, open.end = T))
plot(al_a1_p1_p1, "threeway")

# Tuning the standardized data on reference depth scale
Angel2_on_Picard1_depth = astrochron::tune(Angel2_standardized, cbind(Angel2_standardized$Angel2_scaled.Center_win[al_a1_p1_p1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_a1_p1_p1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (metres)", ylab = "Gamma Ray (gAPI)", cex.lab = 1.25)
lines(Angel2_on_Picard1_depth, col = "red")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(-20,-10,0,10,20), cex.axis = 1.25, las = 1)
legend(x = max(Picard1_standardized$Picard1_scaled.Center_win)-320, y = max(Picard1_standardized$Picard1_scaled.Average)+8, legend = c("Picard-1", "Angel-2"), col = c("black", "red"), lty = 1, lwd = 3, cex = 1.5, bty = "n")

# DTW Distance 
al_a1_p1_p1$normalizedDistance
al_a1_p1_p1$distance

#### DTW with stratigraphy-optimized step pattern asymmetricP1.1 and knowledge-based window ####

# create matrix for the knowledge-based window

compare.window <- matrix(data=TRUE,nrow=nrow(Angel2_standardized),ncol=nrow(Picard1_standardized))
image(x=Picard1_standardized[,1],y=Angel2_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- DescTools::Closest(260, Picard1_standardized[,1],which=TRUE)
base_1_y <- DescTools::Closest(225, Angel2_standardized[,1],which=TRUE)

# Depth values for second datum
base_2_x <- DescTools::Closest(370, Picard1_standardized[,1],which=TRUE)
base_2_y <- DescTools::Closest(350, Angel2_standardized[,1],which=TRUE)

# Depth values for third datum
base_3_x <- DescTools::Closest(430, Picard1_standardized[,1],which=TRUE)
base_3_y <- DescTools::Closest(430, Angel2_standardized[,1],which=TRUE)

# Depth values for fourth datum
base_4_x <- DescTools::Closest(545, Picard1_standardized[,1],which=TRUE)
base_4_y <- DescTools::Closest(630, Angel2_standardized[,1],which=TRUE)

# Depth values for fifth datum
base_5_x <- DescTools::Closest(1010, Picard1_standardized[,1],which=TRUE)
base_5_y <- DescTools::Closest(975, Angel2_standardized[,1],which=TRUE)

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Angel2_standardized), ncol = nrow(Picard1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+200):nrow(Angel2_standardized),1:(base_1_x-200)] <- 0
compare.window[1:(base_1_y-200),(base_1_x+200):ncol(compare.window)] <- 0

compare.window[(base_2_y+200):nrow(Angel2_standardized),1:(base_2_x-200)] <- 0
compare.window[1:(base_2_y-200),(base_2_x+200):ncol(compare.window)] <- 0

compare.window[(base_3_y+400):nrow(Angel2_standardized),1:(base_3_x-400)] <- 0
compare.window[1:(base_3_y-400),(base_3_x+400):ncol(compare.window)] <- 0

compare.window[(base_4_y+400):nrow(Angel2_standardized),1:(base_4_x-400)] <- 0
compare.window[1:(base_4_y-400),(base_4_x+400):ncol(compare.window)] <- 0

compare.window[(base_5_y+400):nrow(Angel2_standardized),1:(base_5_x-400)] <- 0
compare.window[1:(base_5_y-400),(base_5_x+400):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Picard1_standardized[,1],y=Angel2_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Picard1_standardized[,1],y=Angel2_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Perform dtw with knowledge-based window
system.time(al_a1_p1_ap1 <- dtw::dtw(Angel2_standardized$Angel2_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_a1_p1_ap1, type = "threeway")

# DTW Distance measure
al_a1_p1_ap1$normalizedDistance
al_a1_p1_ap1$distance

# Dtw Density plot

plot(al_a1_p1_ap1, type = "density", xaxt = "n", yaxt = "n", xlab = "Angel-2", ylab = "Picard-1", cex.lab = 1.25)
axis(1, at = c(776,2026,3276,4526,5776), labels = c(250,500,750,1000,1250), cex.axis = 1.25)
axis(2, at = c(251,1251,2251,3251,4251,5251), labels = c(200,400,600,800,1000,1200), cex.axis = 1.25)
points(border_coords, cex = 0.3)

# Dtw knowledge-based window plot

image(y = Picard1_standardized[,1], x = Angel2_standardized[,1], z = compare.window, useRaster = T)
lines(Angel2_standardized$Angel2_scaled.Center_win[al_a1_p1_ap1$index1], Picard1_standardized$Picard1_scaled.Center_win[al_a1_p1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Angel2_on_Picard1_depth_cw = tune(Angel2_standardized, cbind(Angel2_standardized$Angel2_scaled.Center_win[al_a1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_a1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (metres)", ylab = "Gamma Ray (gAPI)", cex.lab = 1.25)
lines(Angel2_on_Picard1_depth_cw, col = "red")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(-20,-10,0,10,20), cex.axis = 1.25, las = 1)
legend(x = max(Picard1_standardized$Picard1_scaled.Center_win)-320, y = max(Picard1_standardized$Picard1_scaled.Average)+8, legend = c("Picard-1", "Angel-2"), col = c("black", "red"), lty = 1, lwd = 3, cex = 1.5, bty = "n")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Angel2_originalGR_on_Picard1_depth = data.frame(Angel2_on_Picard1_depth_cw$X1, Angel2_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (metres)", ylab = "Gamma Ray (gAPI)", cex.lab = 1.25)
lines(Angel2_originalGR_on_Picard1_depth, col = "red")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.25, las = 1)
legend(x = max(Picard1_originalGR$Picard1_standardized.Picard1_scaled.Center_win)-340, y = max(Picard1_originalGR$Picard1_interpolated.GR)+6, legend = c("Picard-1", "Angel-2"), col = c("black", "red"), lty = 1, lwd = 3, cex = 1.5, bty = "n")
