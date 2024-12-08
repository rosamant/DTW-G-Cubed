
setwd("C:/Users/Rohit/OneDrive - Universität Münster/Paper Drafts/P&P/Figures/")

pdf(file = "Figure 9.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 9.png", width = 5000, height = 6000, res = 600)

# Result Plot with no knowledge-based window and with knowledge-based window

par(mar=c(0,5,1,1))
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), widths=c(1,1,1), heights=c(1,1.3,1.4))

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "Picard1 Depth (meters)", ylab = "Natural Gamma Radiation (gAPI)", cex.lab = 1.25)
lines(Minilya1_on_Picard1_depth, col = "blue")
title(main = "a) Correlation with asymmetricP1 step pattern and no window")
axis(2, at = c(-20,0,20), cex.axis = 1.25, las = 1)
legend(x = max(Picard1_standardized$Picard1_scaled.Center_win)-320, y = max(Picard1_standardized$Picard1_scaled.Average)+8, legend = c("Picard-1", "Minilya-1"), col = c("black", "blue"), lty = 1, lwd = 3, cex = 1.5, bty = "n")

par(mar=c(5,5,1,1))

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), axes = FALSE, yaxt = "n", 
     xlab = "Picard1 Depth (meters)", ylab = "Natural Gamma Radiation (gAPI)", cex.lab = 1.25)
lines(Minilya1_on_Picard1_depth_cw, col = "blue")
axis(1, at = c(150,400,600,800,1000,1200,1320), cex.axis = 1.25, las = 1)
axis(2, at = c(-20,0,20), cex.axis = 1.25, las = 1)
title(main = "b) Correlation with stratigraphy-optimized step pattern and knowledge-based window")

plot(al_m1_p1_ap1, type = "density", xaxt = "n", yaxt = "n", xlab = "Minilya-1", ylab = "Picard-1", cex.lab = 1.5)
axis(1, at = c(113,1113,2113,3113,4113), labels = c(200,400,600,800,1000), cex.axis = 1.5)
axis(2, at = c(251,1251,2251,3251,4251,5251), labels = c(200,400,600,800,1000,1200), cex.axis = 1.5)
#text(x = 1080, y = 5000,"(c) Cumulative cost within knowledge-based window", cex = 1.75)
title(main = "c) Cumulative cost within knowledge-based window")
points(border_coords, cex = 0.3)

dev.off()
