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

#### DTW with different step pattern####

#### Perform dtw with asymmetric ####
system.time(al_U1464_p1_ap1 <- dtw(U1464$HSGR, Picard1$GR, keep.internals = T, step.pattern = asymmetric, open.begin = F, open.end = T))
plot(al_U1464_p1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth1 = tune(U1464, cbind(U1464$DEPTH_WMSF[al_U1464_p1_ap1$index1s], Picard1$DEPT[al_U1464_p1_ap1$index2s]), extrapolate = T)

dev.off()

#### Perform dtw with asymmetricP05 ####
system.time(al_U1464_p1_ap2 <- dtw(U1464$HSGR, Picard1$GR, keep.internals = T, step.pattern = asymmetricP05, open.begin = F, open.end = T))
plot(al_U1464_p1_ap2, "threeway")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth2 = tune(U1464, cbind(U1464$DEPTH_WMSF[al_U1464_p1_ap2$index1s], Picard1$DEPT[al_U1464_p1_ap2$index2s]), extrapolate = T)

dev.off()

#### Perform dtw with asymmetricP1 ####
system.time(al_U1464_p1_ap3 <- dtw(U1464$HSGR, Picard1$GR, keep.internals = T, step.pattern = asymmetricP1, open.begin = F, open.end = T))
plot(al_U1464_p1_ap3, "threeway")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth3 = tune(U1464, cbind(U1464$DEPTH_WMSF[al_U1464_p1_ap3$index1s], Picard1$DEPT[al_U1464_p1_ap3$index2s]), extrapolate = T)

dev.off()

#### Perform dtw with asymmetricP2 ####
system.time(al_U1464_p1_ap4 <- dtw(U1464$HSGR, Picard1$GR, keep.internals = T, step.pattern = asymmetricP2, open.begin = F, open.end = T))
plot(al_U1464_p1_ap4, "threeway")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth4 = tune(U1464, cbind(U1464$DEPTH_WMSF[al_U1464_p1_ap4$index1s], Picard1$DEPT[al_U1464_p1_ap4$index2s]), extrapolate = T)

dev.off()

png(filename = "Figure 4.png", width = 6000, height = 8000, res = 600)

par(mar=c(0,5,1,1))
layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE), widths=c(1,1), heights=c(1,1,1,1.25))

plot(Picard1, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1464_on_Picard1_depth1, col = "orange")
#axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.5, las = 1)
legend(x = max(Picard1$DEPT)-450, y = max(Picard1$GR)+2, legend = c("Picard-1", "U1464"), col = c("black", "orange"), lty = 1, lwd = 3, cex = 1.70, bty = "n")
text(220,45,"(a) asymmetric", cex = 1.6)

par(mar=c(0,5,0,1))

plot(Picard1, type = "l", ylim = c(0, 50), xlim =  c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1464_on_Picard1_depth2, col = "orange")
#axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.5, las = 1)
#legend(x = max(Picard1$DEPT)-250, y = max(Picard1$GR)+2, legend = c("Picard-1", "U1464"), col = c("black", "orange"), lty = 1, lwd = 3, cex = 1.5, bty = "n")
text(220,45,"(b) Asymmetric P05", cex = 1.6)

par(mar=c(0,5,0,1))

plot(Picard1, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1464_on_Picard1_depth3, col = "orange")
#axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.5, las = 1)
#legend(x = max(Picard1$DEPT)-250, y = max(Picard1$GR)+2, legend = c("Picard-1", "U1464"), col = c("black", "orange"), lty = 1, lwd = 3, cex = 1.5, bty = "n")
text(220,45,"(c) Asymmetric P1", cex = 1.6)

par(mar=c(5,5,0,1))

plot(Picard1, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(U1464_on_Picard1_depth4, col = "orange")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.5, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.5, las = 1)
#legend(x = max(Picard1$DEPT)-250, y = max(Picard1$GR)+2, legend = c("Picard-1", "U1464"), col = c("black", "orange"), lty = 1, lwd = 3, cex = 1.5, bty = "n")
text(220,45,"(d) Asymmetric P2", cex = 1.6)

dev.off()












