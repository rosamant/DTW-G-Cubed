
setwd("C:/Users/Rohit/OneDrive - Universität Münster/Paper Drafts/G-cubed/Figures/")

pdf(file = "Figure 11.pdf", width = 8.27, height = 11, paper = "a4")

png(filename = "Figure 11.png", width = 5000, height = 6000, res = 600)

# Result Plot of original GR values with custom window

par(mar=c(0,5,1,1))
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), widths=c(1,1,1), heights=c(1,1,1.25))

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)", ylab = "Natural Gamma Radiation (gAPI)", cex.lab = 1.5)
lines(Angel2_originalGR_on_Picard1_depth, col = "red")
#axis(1, at = c(150,200,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.5, las = 1)
legend(x = max(Picard1_originalGR$Picard1_standardized.Picard1_scaled.Center_win)-335, y = max(Picard1_originalGR$Picard1_interpolated.GR)+6, legend = c("Picard-1", "Angel-2"), col = c("black", "red"), lty = 1, lwd = 3, cex = 2, bty = "n")
text(200,50,"(a)", cex = 1.75)

par(mar=c(0,5,1,1))

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)", ylab = "Natural Gamma Radiation (gAPI)", cex.lab = 1.5)
lines(Minilya1_originalGR_on_Picard1_depth, col = "blue")
#axis(1, at = c(150,200,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.5, las = 1)
legend(x = max(Picard1_originalGR$Picard1_standardized.Picard1_scaled.Center_win)-335, y = max(Picard1_originalGR$Picard1_interpolated.GR)+6, legend = c("Picard-1", "Minilya-1"), col = c("black", "blue"), lty = 1, lwd = 3, cex = 2, bty = "n")
text(200,50,"(b)", cex = 1.75)

par(mar=c(5,5,0,1))

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard-1 Depth (meters)", ylab = "Natural Gamma Radiation (gAPI)", cex.lab = 1.5)
lines(Whitetail1_originalGR_on_Picard1_depth, col = "purple")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(0,10,20,30,40,50), cex.axis = 1.5, las = 1)
legend(x = max(Picard1_originalGR$Picard1_standardized.Picard1_scaled.Center_win)-335, y = max(Picard1_originalGR$Picard1_interpolated.GR)+6, legend = c("Picard-1", "Whitetail-1"), col = c("black", "purple"), lty = 1, lwd = 3, cex = 2, bty = "n")
text(200,50,"(c)", cex = 1.75)

dev.off()
