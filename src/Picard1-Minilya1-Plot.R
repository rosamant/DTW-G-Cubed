
setwd("C:/Users/Rohit/OneDrive - Universität Münster/Paper Drafts/G-cubed/Figures/")

pdf(file = "Figure 9.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 9.png", width = 5000, height = 6000, res = 600)

# Result Plot with no custom window and with custom window

par(mar=c(0,5,1,1))
layout(matrix(c(1,2), 2, 1, byrow = TRUE), widths=c(1,1), heights=c(1,1.25))

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard1 Depth (meters)", ylab = "Natural Gamma Radiation (gAPI)", cex.lab = 1.25)
lines(Minilya1_on_Picard1_depth, col = "blue")
title(main = "a) Correlation with no custom step pattern and window")
axis(2, at = c(-20,0,20), cex.axis = 1.25, las = 1)
legend(x = max(Picard1_standardized$Picard1_scaled.Center_win)-320, y = max(Picard1_standardized$Picard1_scaled.Average)+8, legend = c("Picard-1", "Minilya-1"), col = c("black", "blue"), lty = 1, lwd = 3, cex = 1.5, bty = "n")
arrows(x0=500, y0 =-12 ,x1 = 440, y1 = 0, length = 0.10, angle = 30, code = 2, lwd = 1.25)
arrows(x0=550, y0 =-12 ,x1 = 580, y1 = 0, length = 0.10, angle = 30, code = 2, lwd = 1.25)
arrows(x0=600, y0 =-12 ,x1 = 880, y1 = -2, length = 0.10, angle = 30, code = 2, lwd = 1.25)
text(x = 530, y = -14, "Over-compression", cex = 1.25)
arrows(x0=550, y0 =12 ,x1 = 650, y1 = 15, length = 0.10, angle = 30, code = 1, lwd = 1.25)
arrows(x0=850, y0 =15 ,x1 = 940, y1 = 10, length = 0.10, angle = 30, code = 2, lwd = 1.25)
text(x = 760, y = 17, "Inaccurate correlation", cex = 1.25)

par(mar=c(5,5,1,1))

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, yaxt = "n", 
     xlab = "Picard1 Depth (meters)", ylab = "Natural Gamma Radiation (gAPI)", cex.lab = 1.25)
lines(Minilya1_on_Picard1_depth_cw, col = "blue")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(-20,0,20), cex.axis = 1.25, las = 1)
title(main = "b) Correlation with with custom step pattern and window")

dev.off()
