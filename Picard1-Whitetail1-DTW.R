# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Picard1 and Whitetail1 datasets

Picard1 <- read.csv("C:/Users/Rohit/OneDrive - Universität Münster/NOPIMS Data/NW_Australia_Digitized_Data/Carnarvon Basin/Picard 1/PICARD 1.csv", header=TRUE, stringsAsFactors=FALSE)
Picard1=Picard1[c(83:7750),] # Eocene-Miocene Unconformity
head(Picard1)
plot(Picard1, type="l", xlim = c(150, 1300), ylim = c(0, 50))

Whitetail1 <- read.csv("C:/Users/Rohit/OneDrive - Universität Münster/NOPIMS Data/NW_Australia_Digitized_Data/Carnarvon Basin/Whitetail 1/Whitetail1.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting attenuated signal
W1 = Gmean(Whitetail1[c(1:320),2])
W2 = Gmean(Whitetail1[c(300:500),2])
SD1 = Gsd(Whitetail1[c(1:320),2])
SD2 = Gsd(Whitetail1[c(300:500),2])
Whitetail1[c(1:320),2]=(Whitetail1[c(1:320),2]+(W2-W1))*(SD1/SD2)

Whitetail1=Whitetail1[c(1:4630),] # Oligocene-Miocene
head(Whitetail1)
plot(Whitetail1, type="l", xlim = c(1000, 1900), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Picard1_interpolated <- linterp(Picard1, dt = 0.2, genplot = F)
Whitetail1_interpolated <- linterp(Whitetail1, dt = 0.2, genplot = F)

# Scaling the data
Pmean = Gmean(Picard1_interpolated$GR)
Pstd = Gsd(Picard1_interpolated$GR)
Picard1_scaled = (Picard1_interpolated$GR - Pmean)/Pstd
Picard1_rescaled = data.frame(Picard1_interpolated$DEPT, Picard1_scaled)

Wmean = Gmean(Whitetail1_interpolated$GR)
Wstd = Gsd(Whitetail1_interpolated$GR)
Whitetail1_scaled = (Whitetail1_interpolated$GR - Wmean)/Wstd
Whitetail1_rescaled = data.frame(Whitetail1_interpolated$DEPT, Whitetail1_scaled)

# Resampling the data using moving window statistics
Picard1_scaled = mwStats(Picard1_rescaled, cols = 2, win=3, ends = T)
Picard1_standardized = data.frame(Picard1_scaled$Center_win, Picard1_scaled$Average)

Whitetail1_scaled = mwStats(Whitetail1_rescaled, cols = 2, win=3, ends = T)
Whitetail1_standardized = data.frame(Whitetail1_scaled$Center_win, Whitetail1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Picard1_standardized, type="l", xlim = c(150, 1300), ylim = c(-20, 20), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
plot(Whitetail1_standardized, type="l", xlim = c(1000, 1900), ylim = c(-20, 20), xlab = "Whitetail1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1 but no custom window ####

# Perform dtw
system.time(al_w1_p1_p1 <- dtw(Whitetail1_standardized$Whitetail1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1, open.begin = F, open.end = T))
plot(al_w1_p1_p1, "threeway")

# Tuning the standardized data on reference depth scale
Whitetail1_on_Picard1_depth = tune(Whitetail1_standardized, cbind(Whitetail1_standardized$Whitetail1_scaled.Center_win[al_w1_p1_p1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_w1_p1_p1$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (metres)", ylab = "Gamma Ray (gAPI)", cex.lab = 1.25)
lines(Whitetail1_on_Picard1_depth, col = "red")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(-20,-10,0,10,20), cex.axis = 1.25, las = 1)
legend(x = max(Picard1_standardized$Picard1_scaled.Center_win)-320, y = max(Picard1_standardized$Picard1_scaled.Average)+8, legend = c("Picard-1", "Whitetail-1"), col = c("black", "red"), lty = 1, lwd = 3, cex = 1.5, bty = "n")

# DTW Distance
al_w1_p1_p1$normalizedDistance
al_w1_p1_p1$distance


#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Whitetail1_standardized),ncol=nrow(Picard1_standardized))
image(x=Picard1_standardized[,1],y=Whitetail1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- Closest(260, Picard1_standardized[,1],which=TRUE)
base_1_y <- Closest(1100, Whitetail1_standardized[,1],which=TRUE)

# Depth values for second datum
base_2_x <- Closest(370, Picard1_standardized[,1],which=TRUE)
base_2_y <- Closest(1190, Whitetail1_standardized[,1],which=TRUE)

# Depth values for third datum
base_3_x <- Closest(430, Picard1_standardized[,1],which=TRUE)
base_3_y <- Closest(1230, Whitetail1_standardized[,1],which=TRUE)

# Depth values for fourth datum
base_4_x <- Closest(545, Picard1_standardized[,1],which=TRUE)
base_4_y <- Closest(1350, Whitetail1_standardized[,1],which=TRUE)

# Depth values for fifth datum
base_5_x <- Closest(1010, Picard1_standardized[,1],which=TRUE)
base_5_y <- Closest(1750, Whitetail1_standardized[,1],which=TRUE)

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Whitetail1_standardized), ncol = nrow(Picard1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+200):nrow(Whitetail1_standardized),1:(base_1_x-200)] <- 0
compare.window[1:(base_1_y-200),(base_1_x+200):ncol(compare.window)] <- 0

compare.window[(base_2_y+300):nrow(Whitetail1_standardized),1:(base_2_x-300)] <- 0
compare.window[1:(base_2_y-300),(base_2_x+300):ncol(compare.window)] <- 0

compare.window[(base_3_y+400):nrow(Whitetail1_standardized),1:(base_3_x-400)] <- 0
compare.window[1:(base_3_y-400),(base_3_x+400):ncol(compare.window)] <- 0

compare.window[(base_4_y+400):nrow(Whitetail1_standardized),1:(base_4_x-400)] <- 0
compare.window[1:(base_4_y-400),(base_4_x+400):ncol(compare.window)] <- 0

compare.window[(base_5_y+500):nrow(Whitetail1_standardized),1:(base_5_x-500)] <- 0
compare.window[1:(base_5_y-500),(base_5_x+500):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Picard1_standardized[,1],y=Whitetail1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Picard1_standardized[,1],y=Whitetail1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_w1_p1_ap1 <- dtw(Whitetail1_standardized$Whitetail1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_w1_p1_ap1, type = "threeway")

# DTW Distance measure
al_w1_p1_ap1$normalizedDistance
al_w1_p1_ap1$distance

image(y = Picard1_standardized[,1], x = Whitetail1_standardized[,1], z = compare.window, useRaster = T)
lines(Whitetail1_standardized$Whitetail1_scaled.Center_win[al_w1_p1_ap1$index1], Picard1_standardized$Picard1_scaled.Center_win[al_w1_p1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Whitetail1_on_Picard1_depth_cw = tune(Whitetail1_standardized, cbind(Whitetail1_standardized$Whitetail1_scaled.Center_win[al_w1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_w1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (metres)", ylab = "Gamma Ray (gAPI)", cex.lab = 1.25)
lines(Whitetail1_on_Picard1_depth_cw, col = "red")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(-20,-10,0,10,20), cex.axis = 1.25, las = 1)
legend(x = max(Picard1_standardized$Picard1_scaled.Center_win)-320, y = max(Picard1_standardized$Picard1_scaled.Average)+8, legend = c("Picard-1", "Whitetail-1"), col = c("black", "red"), lty = 1, lwd = 3, cex = 1.5, bty = "n")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Whitetail1_originalGR_on_Picard1_depth = data.frame(Whitetail1_on_Picard1_depth_cw$X1, Whitetail1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (metres)", ylab = "Gamma Ray (gAPI)", cex.lab = 1.25)
lines(Whitetail1_originalGR_on_Picard1_depth, col = "red")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.25, las = 1)
legend(x = max(Picard1_originalGR$Picard1_standardized.Picard1_scaled.Center_win)-340, y = max(Picard1_originalGR$Picard1_interpolated.GR)+6, legend = c("Picard-1", "Whitetail-1"), col = c("black", "red"), lty = 1, lwd = 3, cex = 1.5, bty = "n")
