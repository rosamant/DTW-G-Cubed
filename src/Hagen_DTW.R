install.packages(setdiff(c("align"), rownames(installed.packages())))

# Import Align library

library(align)

# Create dataframes to perform Hagen DTW algorithm

A_Picard1_standardized = data.frame(Picard1_scaled$Average, Picard1_scaled$Center_win)
write.csv(A_Picard1_standardized, file = "C:/Users/Rohit/Desktop/align/A_Picard1_standardized.csv", row.names = FALSE)

A_Angel2_standardized = data.frame(Angel2_scaled$Average, Angel2_scaled$Center_win)
write.csv(A_Angel2_standardized, file = "C:/Users/Rohit/Desktop/align/A_Angel2_standardized.csv", row.names = FALSE)

A_Minilya1_standardized = data.frame(Minilya1_scaled$Average, Minilya1_scaled$Center_win)
write.csv(A_Minilya1_standardized, file = "C:/Users/Rohit/Desktop/align/A_Minilya1_standardized.csv", row.names = FALSE)

A_Whitetail1_standardized = data.frame(Whitetail1_scaled$Average, Whitetail1_scaled$Center_win)
write.csv(A_Whitetail1_standardized, file = "C:/Users/Rohit/Desktop/align/A_Whitetail1_standardized.csv", row.names = FALSE)

A_U1464_standardized = data.frame(U1464_scaled$Average, U1464_scaled$Center_win)
write.csv(A_U1464_standardized, file = "C:/Users/Rohit/Desktop/align/A_U1464_standardized.csv", row.names = FALSE)

A_U1464_1_standardized = data.frame(U1464_1_scaled$Average, U1464_1_scaled$Center_win)
write.csv(A_U1464_1_standardized, file = "C:/Users/Rohit/Desktop/align/A_U1464_1_standardized.csv", row.names = FALSE)

A_U1482_standardized = data.frame(U1482_scaled$Average, U1482_scaled$Center_win)
write.csv(A_U1482_standardized, file = "C:/Users/Rohit/Desktop/align/A_U1482_standardized.csv", row.names = FALSE)

# Perform DTW between Picard-1 and U1464 for different g and edge paramters

system.time(A1_U1464_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_U1464_standardized, g = 0.99, edge = 0.15))
A1_U1464_Picard1 = data.frame(A1_U1464_Picard1_dtw[[2]], A1_U1464_Picard1_dtw[[1]])
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A1_U1464_Picard1, col = "red")

system.time(A2_U1464_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_U1464_standardized, g = 1.0, edge = 0.15))
A2_U1464_Picard1 = data.frame(A2_U1464_Picard1_dtw[[2]], A2_U1464_Picard1_dtw[[1]])
A2_U1464_Picard1 <- A2_U1464_Picard1[order(A2_U1464_Picard1$A2_U1464_Picard1_dtw..2..), ]
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A2_U1464_Picard1, col = "red")

system.time(A3_U1464_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_U1464_standardized, g = 1.05, edge = 0.25))
A3_U1464_Picard1 = data.frame(A3_U1464_Picard1_dtw[[2]], A3_U1464_Picard1_dtw[[1]])
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A3_U1464_Picard1, col = "red")

system.time(A4_U1464_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_U1464_standardized, g = 1.15, edge = 0.2))
A4_U1464_Picard1 = data.frame(A4_U1464_Picard1_dtw[[2]], A4_U1464_Picard1_dtw[[1]])
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A4_U1464_Picard1, col = "red")


############################
# Make Figure for Picard-1 U1464 correlation
###########################

setwd("Supplementary Figures/")

pdf(file = "Figure 8.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 8.png", width = 4000, height = 6000, res = 600)

par(mar=c(0,5,1,1))
layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE), widths=c(1,1,1,1), heights=c(1,1,1,1.5))

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A1_U1464_Picard1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
legend(x = max(A_Picard1_standardized$Picard1_scaled.Center_win)-450, y = max(A_Picard1_standardized$Picard1_scaled.Average)+6, legend = c("Picard-1", "U1464"), col = c("black", "red"), lty = 1, lwd = 2, cex = 1.40, bty = "n")
text(150,15,"(a)", cex = 1.6)
text(1200,15,"g=0.99", cex = 1.5)
text(1200,10,"edge=0.15", cex = 1.5)

par(mar=c(0,5,0,1))

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A2_U1464_Picard1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(150,15,"(b)", cex = 1.6)
text(1200,15,"g=1.0", cex = 1.5)
text(1200,10,"edge=0.15", cex = 1.5)

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A3_U1464_Picard1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(150,15,"(c)", cex = 1.6)
text(1200,15,"g=1.05", cex = 1.5)
text(1200,10,"edge=0.25", cex = 1.5)

par(mar=c(5,5,0,1))

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A4_U1464_Picard1, col = "red")
axis(1, at = c(150, seq(300,1300,100)), cex.axis = 1.5, las = 1)
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(150,15,"(d)", cex = 1.6)
text(1200,15,"g=1.15", cex = 1.5)
text(1200,10,"edge=0.2", cex = 1.5)

dev.off()


# Perform DTW between U1464 and U1482 for different g and edge paramters

system.time(A1_U1482_U1464_1_dtw <- dtw_r(A_U1464_1_standardized, A_U1482_standardized, g = 0.99, edge = 0.15))
A1_U1482_U1464_1 = data.frame(A1_U1482_U1464_1_dtw[[2]], A1_U1482_U1464_1_dtw[[1]])
plot(A_U1464_1_standardized$U1464_1_scaled.Center_win, A_U1464_1_standardized$U1464_1_scaled.Average, type = "l")
lines(A1_U1482_U1464_1, col = "red")

system.time(A2_U1482_U1464_1_dtw <- dtw_r(A_U1464_1_standardized, A_U1482_standardized, g = 1.0, edge = 0.15))
A2_U1482_U1464_1 = data.frame(A2_U1482_U1464_1_dtw[[2]], A2_U1482_U1464_1_dtw[[1]])
A2_U1482_U1464_1 <- A2_U1482_U1464_1[order(A2_U1482_U1464_1$A2_U1482_U1464_1_dtw..2..), ]
plot(A_U1464_1_standardized$U1464_1_scaled.Center_win, A_U1464_1_standardized$U1464_1_scaled.Average, type = "l")
lines(A2_U1482_U1464_1, col = "red")

system.time(A3_U1482_U1464_1_dtw <- dtw_r(A_U1464_1_standardized, A_U1482_standardized, g = 1.05, edge = 0.25))
A3_U1482_U1464_1 = data.frame(A3_U1482_U1464_1_dtw[[2]], A3_U1482_U1464_1_dtw[[1]])
plot(A_U1464_1_standardized$U1464_1_scaled.Center_win, A_U1464_1_standardized$U1464_1_scaled.Average, type = "l")
lines(A3_U1482_U1464_1, col = "red")

system.time(A4_U1482_U1464_1_dtw <- dtw_r(A_U1464_1_standardized, A_U1482_standardized, g = 1.15, edge = 0.4))
A4_U1482_U1464_1 = data.frame(A4_U1482_U1464_1_dtw[[2]], A4_U1482_U1464_1_dtw[[1]])
plot(A_U1464_1_standardized$U1464_1_scaled.Center_win, A_U1464_1_standardized$U1464_1_scaled.Average, type = "l")
lines(A4_U1482_U1464_1, col = "red")


############################
# Make Figure for U1464 U1482 correlation
###########################

setwd("Supplementary Figures/")

pdf(file = "Figure 10.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 10.png", width = 4000, height = 6000, res = 600)

par(mar=c(0,5,1,1))
layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE), widths=c(1,1,1,1), heights=c(1,1,1,1.5))

plot(A_U1464_1_standardized$U1464_1_scaled.Center_win, A_U1464_1_standardized$U1464_1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(0, 500), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A1_U1482_U1464_1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
legend(x = max(A_U1464_1_standardized$U1464_1_scaled.Center_win)-150, y = max(A_U1464_1_standardized$U1464_1_scaled.Average)+6, legend = c("U1464", "U1482"), col = c("black", "red"), lty = 1, lwd = 2, cex = 1.40, bty = "n")
text(20,15,"(a)", cex = 1.6)
text(450,15,"g=0.99", cex = 1.5)
text(450,10,"edge=0.15", cex = 1.5)

par(mar=c(0,5,0,1))

plot(A_U1464_1_standardized$U1464_1_scaled.Center_win, A_U1464_1_standardized$U1464_1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(0, 500), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A2_U1482_U1464_1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(20,15,"(b)", cex = 1.6)
text(450,15,"g=1.0", cex = 1.5)
text(450,10,"edge=0.15", cex = 1.5)

plot(A_U1464_1_standardized$U1464_1_scaled.Center_win, A_U1464_1_standardized$U1464_1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(0, 500), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A3_U1482_U1464_1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(20,15,"(c)", cex = 1.6)
text(450,15,"g=1.05", cex = 1.5)
text(450,10,"edge=0.25", cex = 1.5)

par(mar=c(5,5,0,1))

plot(A_U1464_1_standardized$U1464_1_scaled.Center_win, A_U1464_1_standardized$U1464_1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(0, 500), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A4_U1482_U1464_1, col = "red")
axis(1, at = c(500, seq(0,500,100)), cex.axis = 1.5, las = 1)
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(20,15,"(d)", cex = 1.6)
text(450,15,"g=1.15", cex = 1.5)
text(450,10,"edge=0.4", cex = 1.5)

dev.off()

# Perform DTW between Picard-1 and Angel-2 for different g and edge paramters

system.time(A1_Angel2_U1464_1_dtw <- dtw_r(A_Picard1_standardized, A_Angel2_standardized, g = 0.99, edge = 0.15))
A1_Angel2_Picard1 = data.frame(A1_Angel2_Picard1_dtw[[2]], A1_Angel2_Picard1_dtw[[1]])
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A1_Angel2_Picard1, col = "red")

system.time(A2_Angel2_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_Angel2_standardized, g = 1.0, edge = 0.15))
A2_Angel2_Picard1 = data.frame(A2_Angel2_Picard1_dtw[[2]], A2_Angel2_Picard1_dtw[[1]])
A2_Angel2_Picard1 <- A2_Angel2_Picard1[order(A2_Angel2_Picard1$A2_Angel2_Picard1_dtw..2..), ]
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A2_Angel2_Picard1, col = "red")

system.time(A3_Angel2_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_Angel2_standardized, g = 1.05, edge = 0.25))
A3_Angel2_Picard1 = data.frame(A3_Angel2_Picard1_dtw[[2]], A3_Angel2_Picard1_dtw[[1]])
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A3_Angel2_Picard1, col = "red")

system.time(A4_Angel2_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_Angel2_standardized, g = 1.15, edge = 0.2))
A4_Angel2_Picard1 = data.frame(A4_Angel2_Picard1_dtw[[2]], A4_Angel2_Picard1_dtw[[1]])
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A4_Angel2_Picard1, col = "red")


############################
# Make Figure for Picard-1 Angel-2 correlation
###########################

setwd("Supplementary Figures/")

pdf(file = "Figure 11.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 11.png", width = 4000, height = 6000, res = 600)

par(mar=c(0,5,1,1))
layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE), widths=c(1,1,1,1), heights=c(1,1,1,1.5))

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A1_Angel2_Picard1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
legend(x = max(A_Picard1_standardized$Picard1_scaled.Center_win)-470, y = max(A_Picard1_standardized$Picard1_scaled.Average)+6, legend = c("Picard-1", "Angel-2"), col = c("black", "red"), lty = 1, lwd = 2, cex = 1.40, bty = "n")
text(150,15,"(a)", cex = 1.6)
text(1200,15,"g=0.99", cex = 1.5)
text(1200,10,"edge=0.15", cex = 1.5)

par(mar=c(0,5,0,1))

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A2_Angel2_Picard1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(150,15,"(b)", cex = 1.6)
text(1200,15,"g=1.0", cex = 1.5)
text(1200,10,"edge=0.15", cex = 1.5)

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A3_Angel2_Picard1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(150,15,"(c)", cex = 1.6)
text(1200,15,"g=1.05", cex = 1.5)
text(1200,10,"edge=0.25", cex = 1.5)

par(mar=c(5,5,0,1))

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A4_Angel2_Picard1, col = "red")
axis(1, at = c(150, seq(300,1300,100)), cex.axis = 1.5, las = 1)
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(150,15,"(d)", cex = 1.6)
text(1200,15,"g=1.15", cex = 1.5)
text(1200,10,"edge=0.2", cex = 1.5)

dev.off()


# Perform DTW between Picard-1 and Minilya-1 for different g and edge paramters

system.time(A1_Minilya1_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_Minilya1_standardized, g = 0.99, edge = 0.15))
A1_Minilya1_Picard1 = data.frame(A1_Minilya1_Picard1_dtw[[2]], A1_Minilya1_Picard1_dtw[[1]])
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A1_Minilya1_Picard1, col = "red")

system.time(A2_Minilya1_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_Minilya1_standardized, g = 1.0, edge = 0.15))
A2_Minilya1_Picard1 = data.frame(A2_Minilya1_Picard1_dtw[[2]], A2_Minilya1_Picard1_dtw[[1]])
A2_Minilya1_Picard1 <- A2_Minilya1_Picard1[order(A2_Minilya1_Picard1$A2_Minilya1_Picard1_dtw..2..), ]
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A2_Minilya1_Picard1, col = "red")

system.time(A3_Minilya1_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_Minilya1_standardized, g = 1.05, edge = 0.25))
A3_Minilya1_Picard1 = data.frame(A3_Minilya1_Picard1_dtw[[2]], A3_Minilya1_Picard1_dtw[[1]])
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A3_Minilya1_Picard1, col = "red")

system.time(A4_Minilya1_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_Minilya1_standardized, g = 1.15, edge = 0.3))
A4_Minilya1_Picard1 = data.frame(A4_Minilya1_Picard1_dtw[[2]], A4_Minilya1_Picard1_dtw[[1]])
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A4_Minilya1_Picard1, col = "red")


############################
# Make Figure for Picard-1 Minilya-1 correlation
###########################

setwd("Supplementary Figures/")

pdf(file = "Figure 12.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 12.png", width = 4000, height = 6000, res = 600)

par(mar=c(0,5,1,1))
layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE), widths=c(1,1,1,1), heights=c(1,1,1,1.5))

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A1_Minilya1_Picard1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
legend(x = max(A_Picard1_standardized$Picard1_scaled.Center_win)-450, y = max(A_Picard1_standardized$Picard1_scaled.Average)+6, legend = c("Picard-1", "Minilya-1"), col = c("black", "red"), lty = 1, lwd = 2, cex = 1.40, bty = "n")
text(150,15,"(a)", cex = 1.6)
text(1200,15,"g=0.99", cex = 1.5)
text(1200,10,"edge=0.15", cex = 1.5)

par(mar=c(0,5,0,1))

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A2_Minilya1_Picard1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(150,15,"(b)", cex = 1.6)
text(1200,15,"g=1.0", cex = 1.5)
text(1200,10,"edge=0.15", cex = 1.5)

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A3_Minilya1_Picard1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(150,15,"(c)", cex = 1.6)
text(1200,15,"g=1.05", cex = 1.5)
text(1200,10,"edge=0.25", cex = 1.5)

par(mar=c(5,5,0,1))

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A4_Minilya1_Picard1, col = "red")
axis(1, at = c(150, seq(300,1300,100)), cex.axis = 1.5, las = 1)
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(150,15,"(d)", cex = 1.6)
text(1200,15,"g=1.15", cex = 1.5)
text(1200,10,"edge=0.3", cex = 1.5)

dev.off()


# Perform DTW between Picard-1 and Whitetail-1 for different g and edge paramters

system.time(A1_Whitetail1_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_Whitetail1_standardized, g = 0.99, edge = 0.15))
A1_Whitetail1_Picard1 = data.frame(A1_Whitetail1_Picard1_dtw[[2]], A1_Whitetail1_Picard1_dtw[[1]])
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A1_Whitetail1_Picard1, col = "red")

system.time(A2_Whitetail1_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_Whitetail1_standardized, g = 1.0, edge = 0.15))
A2_Whitetail1_Picard1 = data.frame(A2_Whitetail1_Picard1_dtw[[2]], A2_Whitetail1_Picard1_dtw[[1]])
A2_Whitetail1_Picard1 <- A2_Whitetail1_Picard1[order(A2_Whitetail1_Picard1$A2_Whitetail1_Picard1_dtw..2..), ]
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A2_Whitetail1_Picard1, col = "red")

system.time(A3_Whitetail1_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_Whitetail1_standardized, g = 1.05, edge = 0.25))
A3_Whitetail1_Picard1 = data.frame(A3_Whitetail1_Picard1_dtw[[2]], A3_Whitetail1_Picard1_dtw[[1]])
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A3_Whitetail1_Picard1, col = "red")

system.time(A4_Whitetail1_Picard1_dtw <- dtw_r(A_Picard1_standardized, A_Whitetail1_standardized, g = 1.15, edge = 0.3))
A4_Whitetail1_Picard1 = data.frame(A4_Whitetail1_Picard1_dtw[[2]], A4_Whitetail1_Picard1_dtw[[1]])
plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l")
lines(A4_Whitetail1_Picard1, col = "red")


############################
# Make Figure for Picard-1 Whitetail-1 correlation
###########################

setwd("Supplementary Figures/")

pdf(file = "Figure 13.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 13.png", width = 4000, height = 6000, res = 600)

par(mar=c(0,5,1,1))
layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE), widths=c(1,1,1,1), heights=c(1,1,1,1.5))

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A1_Whitetail1_Picard1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
legend(x = max(A_Picard1_standardized$Picard1_scaled.Center_win)-550, y = max(A_Picard1_standardized$Picard1_scaled.Average)+6, legend = c("Picard-1", "Whitetail-1"), col = c("black", "red"), lty = 1, lwd = 2, cex = 1.40, bty = "n")
text(150,15,"(a)", cex = 1.6)
text(1200,15,"g=0.99", cex = 1.5)
text(1200,10,"edge=0.15", cex = 1.5)

par(mar=c(0,5,0,1))

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "",ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A2_Whitetail1_Picard1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(150,15,"(b)", cex = 1.6)
text(1200,15,"g=1.0", cex = 1.5)
text(1200,10,"edge=0.15", cex = 1.5)

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A3_Whitetail1_Picard1, col = "red")
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(150,15,"(c)", cex = 1.6)
text(1200,15,"g=1.05", cex = 1.5)
text(1200,10,"edge=0.25", cex = 1.5)

par(mar=c(5,5,0,1))

plot(A_Picard1_standardized$Picard1_scaled.Center_win, A_Picard1_standardized$Picard1_scaled.Average, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)", ylab = "NGR (gAPI)", cex.lab = 1.5)
lines(A4_Whitetail1_Picard1, col = "red")
axis(1, at = c(150, seq(300,1300,100)), cex.axis = 1.5, las = 1)
axis(2, at = c(-15,0,15), cex.axis = 1.25, las = 1)
text(150,15,"(d)", cex = 1.6)
text(1200,15,"g=1.15", cex = 1.5)
text(1200,10,"edge=0.3", cex = 1.5)

dev.off()
