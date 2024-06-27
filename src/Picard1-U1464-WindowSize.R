# Install packages if missing
install.packages(setdiff(c("dtw", "DescTools", "astrochron"), rownames(installed.packages())))

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

# Recorrecting attenuated signal
M1 = Gmean(U1464[c(1:551),2])
M2 = Gmean(U1464[c(535:900),2])
SD1 = Gsd(U1464[c(1:551),2])
SD2 = Gsd(U1464[c(535:900),2])
U1464[c(1:553),2]=(U1464[c(1:553),2]+(M2-M1))*(SD1/SD2)

head(U1464)
plot(U1464, type="l", xlim = c(0, 800), ylim = c(0, 60))

#### DTW with asymmetricP1 step pattern and different window size ####


#### Perform dtw with window size 500 ####
system.time(al_U1464_p1_ap5 <- dtw(U1464$HSGR, Picard1$GR, keep.internals = T, step.pattern = asymmetricP1, window.type = "sakoechiba", window.size = 500, open.begin = F, open.end = T))
plot(al_U1464_p1_ap5, "threeway")
plot(al_U1464_p1_ap5, "density")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth5 = tune(U1464, cbind(U1464$DEPTH_WMSF[al_U1464_p1_ap5$index1s], Picard1$DEPT[al_U1464_p1_ap5$index2s]), extrapolate = T)

dev.off()

#### Perform dtw with window size 1000 ####
system.time(al_U1464_p1_ap6 <- dtw(U1464$HSGR, Picard1$GR, keep.internals = T, step.pattern = asymmetricP1, window.type = "sakoechiba", window.size = 1000, open.begin = F, open.end = T))
plot(al_U1464_p1_ap6, "threeway")
plot(al_U1464_p1_ap6, "density")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth6 = tune(U1464, cbind(U1464$DEPTH_WMSF[al_U1464_p1_ap6$index1s], Picard1$DEPT[al_U1464_p1_ap6$index2s]), extrapolate = T)

dev.off()

#### Perform dtw with window size 2000 ####
system.time(al_U1464_p1_ap7 <- dtw(U1464$HSGR, Picard1$GR, keep.internals = T, step.pattern = asymmetricP1, window.type = "sakoechiba", window.size = 2000, open.begin = F, open.end = T))
plot(al_U1464_p1_ap7, "threeway")
plot(al_U1464_p1_ap7, "density")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth7 = tune(U1464, cbind(U1464$DEPTH_WMSF[al_U1464_p1_ap7$index1s], Picard1$DEPT[al_U1464_p1_ap7$index2s]), extrapolate = T)

dev.off()

pdf(file = "Figure 1.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 1.png", width = 6000, height = 4000, res = 600)

par(mar=c(2,5,1,1))
layout(matrix(c(1,2), 2, 1, byrow = TRUE), widths=c(1,1), heights=c(1,1))

plot(al_U1464_p1_ap5, "density", xaxt = "n", yaxt = "n", xlab = "U1464", ylab = "Picard-1", cex.lab = 1.5)
axis(1, at = c(1313,2625,3938), labels = c(200,400,600), cex.axis = 1.5)
axis(2, at = c(1723,4348,6973), labels = c(400,800,1200), cex.axis = 1.5)
text(x = 220, y = 7000,"(a)", cex = 1.75)

par(mar=c(5,5,0,1))

plot(Picard1, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1464_on_Picard1_depth5, col = "orange")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.5, las = 1)
legend(x = max(Picard1$DEPT)-200, y = max(Picard1$GR)+2, legend = c("Picard-1", "U1464"), col = c("black", "orange"), lty = 1, lwd = 3, cex = 1.25, bty = "n")
text(x = 160, y = 45,"(b)", cex = 1.75)

dev.off()

png(filename = "Figure 2.png", width = 6000, height = 4000, res = 600)

par(mar=c(2,5,1,1))
layout(matrix(c(1,2), 2, 1, byrow = TRUE), widths=c(1,1), heights=c(1,1))

plot(al_U1464_p1_ap6, "density", xaxt = "n", yaxt = "n", xlab = "U1464", ylab = "Picard-1", cex.lab = 1.5)
axis(1, at = c(1313,2625,3938), labels = c(200,400,600), cex.axis = 1.5)
axis(2, at = c(1723,4348,6973), labels = c(400,800,1200), cex.axis = 1.5)
text(x = 220, y = 7000,"(a)", cex = 1.75)

par(mar=c(5,5,0,1))

plot(Picard1, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1464_on_Picard1_depth6, col = "orange")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.5, las = 1)
legend(x = max(Picard1$DEPT)-200, y = max(Picard1$GR)+2, legend = c("Picard-1", "U1464"), col = c("black", "orange"), lty = 1, lwd = 3, cex = 1.25, bty = "n")
text(x = 160, y = 45,"(b)", cex = 1.75)

dev.off()

png(filename = "Figure 5.png", width = 6000, height = 4000, res = 600)

par(mar=c(2,5,1,1))
layout(matrix(c(1,2), 2, 1, byrow = TRUE), widths=c(1,1), heights=c(1,1))

plot(al_U1464_p1_ap7, "density", xaxt = "n", yaxt = "n", xlab = "U1464", ylab = "Picard-1", cex.lab = 1.5)
axis(1, at = c(1313,2625,3938), labels = c(200,400,600), cex.axis = 1.5)
axis(2, at = c(1723,4348,6973), labels = c(400,800,1200), cex.axis = 1.5)
text(x = 220, y = 7000,"(a)", cex = 1.75)

par(mar=c(5,5,0,1))

plot(Picard1, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1464_on_Picard1_depth7, col = "orange")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.5, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.5, las = 1)
legend(x = max(Picard1$DEPT)-200, y = max(Picard1$GR)+2, legend = c("Picard-1", "U1464"), col = c("black", "orange"), lty = 1, lwd = 3, cex = 1.25, bty = "n")
text(x = 160, y = 45,"(b)", cex = 1.75)

dev.off()
