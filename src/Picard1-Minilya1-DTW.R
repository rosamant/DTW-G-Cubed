install.packages(setdiff(c("dtw", "DescTools", "astrochron"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Picard1 and Minilya1 datasets

Picard1 <- read.csv("data/PICARD 1.csv", header=TRUE, stringsAsFactors=FALSE)
Picard1=Picard1[c(83:7750),] # Eocene-Miocene Unconformity
head(Picard1)
plot(Picard1, type="l", xlim = c(150, 1300), ylim = c(0, 50))

Minilya1 <- read.csv("data/Minilya_1.csv", header=TRUE, stringsAsFactors=FALSE)
Minilya1=Minilya1[c(1:5613),] # Eocene-Miocene Unconformity
head(Minilya1)
plot(Minilya1, type="l", xlim = c(150, 1100), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Picard1_interpolated <- linterp(Picard1, dt = 0.2, genplot = F)
Minilya1_interpolated <- linterp(Minilya1, dt = 0.2, genplot = F)

# Scaling the data
Pmean = Gmean(Picard1_interpolated$GR)
Pstd = Gsd(Picard1_interpolated$GR)
Picard1_scaled = (Picard1_interpolated$GR - Pmean)/Pstd
Picard1_rescaled = data.frame(Picard1_interpolated$DEPT, Picard1_scaled)

Mmean = Gmean(Minilya1_interpolated$GR)
Mstd = Gsd(Minilya1_interpolated$GR)
Minilya1_scaled = (Minilya1_interpolated$GR - Mmean)/Mstd
Minilya1_rescaled = data.frame(Minilya1_interpolated$DEPT, Minilya1_scaled)

# Resampling the data using moving window statistics
Picard1_scaled = mwStats(Picard1_rescaled, cols = 2, win = 3, ends = T)
Picard1_standardized = data.frame(Picard1_scaled$Center_win, Picard1_scaled$Average)

Minilya1_scaled = mwStats(Minilya1_rescaled, cols = 2, win = 3, ends = T)
Minilya1_standardized = data.frame(Minilya1_scaled$Center_win, Minilya1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Picard1_standardized, type="l", xlim = c(150, 1300), ylim = c(-20, 20), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
plot(Minilya1_standardized, type="l", xlim = c(150, 1100), ylim = c(-20, 20), xlab = "Minilya1 Resampled Depth", ylab = "Normalized GR")

#### DTW with step pattern asymmetricP1 and no knowledge-based step pattern or window ####

# Perform dtw
system.time(al_m1_p1_p1 <- dtw(Minilya1_standardized$Minilya1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1, open.begin = F, open.end = T))
plot(al_m1_p1_p1, "threeway")

# Tuning the standardized data on reference depth scale
Minilya1_on_Picard1_depth = tune(Minilya1_standardized, cbind(Minilya1_standardized$Minilya1_scaled.Center_win[al_m1_p1_p1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_m1_p1_p1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (metres)", ylab = "Gamma Ray (gAPI)", cex.lab = 1.25)
lines(Minilya1_on_Picard1_depth, col = "red")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(-20,-10,0,10,20), cex.axis = 1.25, las = 1)
legend(x = max(Picard1_standardized$Picard1_scaled.Center_win)-320, y = max(Picard1_standardized$Picard1_scaled.Average)+10, legend = c("Picard-1", "Minilya-1"), col = c("black", "red"), lty = 1, lwd = 3, cex = 1.5, bty = "n")

# DTW Distance measure
al_m1_p1_p1$normalizedDistance
al_m1_p1_p1$distance


#### DTW with stratigraphy-optimized step pattern asymmetricP1.1 and knowledge-based window ####

# create matrix for the knowledge-based window

compare.window <- matrix(data=TRUE,nrow=nrow(Minilya1_standardized),ncol=nrow(Picard1_standardized))
image(x=Picard1_standardized[,1],y=Minilya1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- Closest(260, Picard1_standardized[,1],which=TRUE)
base_1_y <- Closest(300, Minilya1_standardized[,1],which=TRUE)

# Depth values for second datum
base_2_x <- Closest(370, Picard1_standardized[,1],which=TRUE)
base_2_y <- Closest(420, Minilya1_standardized[,1],which=TRUE)

# Depth values for third datum
base_3_x <- Closest(430, Picard1_standardized[,1],which=TRUE)
base_3_y <- Closest(510, Minilya1_standardized[,1],which=TRUE)

# Depth values for fourth datum
base_4_x <- Closest(545, Picard1_standardized[,1],which=TRUE)
base_4_y <- Closest(620, Minilya1_standardized[,1],which=TRUE)

# Depth values for fifth datum
base_5_x <- Closest(993, Picard1_standardized[,1],which=TRUE)
base_5_y <- Closest(790, Minilya1_standardized[,1],which=TRUE)

# Depth values for sixth datum
base_6_x <- Closest(1190, Picard1_standardized[,1],which=TRUE)
base_6_y <- Closest(940, Minilya1_standardized[,1],which=TRUE)

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Minilya1_standardized), ncol = nrow(Picard1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+100):nrow(Minilya1_standardized),1:(base_1_x-100)] <- 0
compare.window[1:(base_1_y-100),(base_1_x+100):ncol(compare.window)] <- 0

compare.window[(base_2_y+100):nrow(Minilya1_standardized),1:(base_2_x-100)] <- 0
compare.window[1:(base_2_y-100),(base_2_x+100):ncol(compare.window)] <- 0

compare.window[(base_3_y+100):nrow(Minilya1_standardized),1:(base_3_x-100)] <- 0
compare.window[1:(base_3_y-100),(base_3_x+100):ncol(compare.window)] <- 0

compare.window[(base_4_y+100):nrow(Minilya1_standardized),1:(base_4_x-100)] <- 0
compare.window[1:(base_4_y-100),(base_4_x+100):ncol(compare.window)] <- 0

compare.window[(base_5_y+200):nrow(Minilya1_standardized),1:(base_5_x-200)] <- 0
compare.window[1:(base_5_y-200),(base_5_x+200):ncol(compare.window)] <- 0

compare.window[(base_6_y+450):nrow(Minilya1_standardized),1:(base_6_x-450)] <- 0
compare.window[1:(base_6_y-200),(base_6_x+200):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Picard1_standardized[,1],y=Minilya1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Picard1_standardized[,1],y=Minilya1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a knowledge-based window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with knowledge-based window
system.time(al_m1_p1_ap1 <- dtw(Minilya1_standardized$Minilya1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_m1_p1_ap1, type = "threeway")

# DTW Distance measure
al_m1_p1_ap1$normalizedDistance
al_m1_p1_ap1$distance

# Dtw Density plot

plot(al_m1_p1_ap1, type = "density", xaxt = "n", yaxt = "n", xlab = "Minilya-1", ylab = "Picard-1", cex.lab = 1.25)
axis(1, at = c(113,1113,2113,3113,4113), labels = c(200,400,600,800,1000), cex.axis = 1.25)
axis(2, at = c(251,1251,2251,3251,4251,5251), labels = c(200,400,600,800,1000,1200), cex.axis = 1.25)
points(border_coords, cex = 0.3)

# Dtw knowledge-based window plot

image(y = Picard1_standardized[,1], x = Minilya1_standardized[,1], z = compare.window, useRaster = T, xlab = "Minilya-1", ylab = "Picard-1", cex.lab = 1.25, cex.axis = 1.25)
lines(Minilya1_standardized$Minilya1_scaled.Center_win[al_m1_p1_ap1$index1], Picard1_standardized$Picard1_scaled.Center_win[al_m1_p1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Minilya1_on_Picard1_depth_cw = tune(Minilya1_standardized, cbind(Minilya1_standardized$Minilya1_scaled.Center_win[al_m1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_m1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (metres)", ylab = "Gamma Ray (gAPI)", cex.lab = 1.25)
lines(Minilya1_on_Picard1_depth_cw, col = "red")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(-20,-10,0,10,20), cex.axis = 1.25, las = 1)
legend(x = max(Picard1_standardized$Picard1_scaled.Center_win)-320, y = max(Picard1_standardized$Picard1_scaled.Average)+10, legend = c("Picard-1", "Minilya-1"), col = c("black", "red"), lty = 1, lwd = 3, cex = 1.5, bty = "n")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Minilya1_originalGR_on_Picard1_depth = data.frame(Minilya1_on_Picard1_depth_cw$X1, Minilya1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (metres)", ylab = "Gamma Ray (gAPI)", cex.lab = 1.25)
lines(Minilya1_originalGR_on_Picard1_depth, col = "red")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.25, las = 1)
legend(x = max(Picard1_originalGR$Picard1_standardized.Picard1_scaled.Center_win)-340, y = max(Picard1_originalGR$Picard1_interpolated.GR)+6, legend = c("Picard-1", "Minilya-1"), col = c("black", "red"), lty = 1, lwd = 3, cex = 1.5, bty = "n")
