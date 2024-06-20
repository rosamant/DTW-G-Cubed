# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Picard1 and Goodwyn6 datasets

Picard1 <- read.csv("C:/Users/Rohit/OneDrive - Universität Münster/NOPIMS Data/NW_Australia_Digitized_Data/Carnarvon Basin/Picard 1/PICARD 1.csv", header=TRUE, stringsAsFactors=FALSE)
Picard1=Picard1[c(83:7750),] # Eocene-Miocene Unconformity
head(Picard1)
plot(Picard1, type="l", xlim = c(150, 1300), ylim = c(0, 50))

Goodwyn6 <- read.csv("C:/Users/Rohit/OneDrive - Universität Münster/NOPIMS Data/NW_Australia_Digitized_Data/Carnarvon Basin/Goodwyn 6/Goodwyn6.csv", header=TRUE, stringsAsFactors=FALSE)
Goodwyn6=Goodwyn6[c(1:10620),] # Oligocene-Miocene
head(Goodwyn6)
plot(Goodwyn6, type="l", xlim = c(150, 1750), ylim = c(0, 50))

#### DTW with custom step pattern asymmetricP1 and no window ####

# Perform dtw
system.time(al_g6_p1_ap1 <- dtw(Goodwyn6$GR, Picard1$GR, keep.internals = T, step.pattern = asymmetricP1, open.begin = T, open.end = F))
plot(al_g6_p1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Goodwyn6_on_Picard1_depth = tune(Goodwyn6, cbind(Goodwyn6$DEPT[al_g6_p1_ap1$index1s], Picard1$DEPT[al_g6_p1_ap1$index2s]), extrapolate = T)

dev.off()

# Plotting the data

setwd("C:/Users/Rohit/OneDrive - Universität Münster/Paper Drafts/Supplementary Figures/")

pdf(file = "Figure 6.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 6.png", width = 6000, height = 4000, res = 600)

plot(Picard1, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard1 Depth (meters)", ylab = "Natural Gamma Radiation (gAPI)", cex.lab = 1.25)
lines(Goodwyn6_on_Picard1_depth, col = "orange")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.25, las = 1)
legend(x = max(Picard1$DEPT)-340, y = max(Picard1$GR)+2, legend = c("Picard-1", "Goodwyn-6"), col = c("black", "orange"), lty = 1, lwd = 3, cex = 1.5, bty = "n")

dev.off()

# DTW Distance
al_g6_p1_ap1$normalizedDistance
al_g6_p1_ap1$distance
