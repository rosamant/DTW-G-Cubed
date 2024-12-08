
# Dtw custom window plot

dev.off()

setwd("C:/Users/Rohit/OneDrive - Universität Münster/Paper Drafts/P&P/Figures/Figure 3/")

png(filename = "CustomDensityWindow.png", width = 6000, height = 8000, res = 600)

par(mar=c(5,5,1,1))
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), widths=c(1,1), heights=c(1,1,1))

dtwWindow.plot(sakoeChibaWindow, window.size = 30, cex.lab = 1.75, cex.axis = 1.5)
text(62,200,"(a) Sakoe-Chiba Window showing available space for warping", cex = 1.75)

par(mar=c(5,5,0,1))

image(y = Picard1_standardized[,1], x = Minilya1_standardized[,1], z = compare.window, useRaster = T, xlab = "Minilya-1", ylab = "Picard-1", cex.lab = 1.5, cex.axis = 1.5)
text(380,1200,"(b) Knowledge-based windows reflecting datums", cex = 1.75)
lines(Minilya1_standardized$Minilya1_scaled.Center_win[al_m1_p1_ap1$index1], Picard1_standardized$Picard1_scaled.Center_win[al_m1_p1_ap1$index2], col = "white", lwd = 2)

par(mar=c(5,5,0,1))

plot(al_m1_p1_ap1, type = "density", xaxt = "n", yaxt = "n", xlab = "Minilya-1", ylab = "Picard-1", cex.lab = 1.5)
axis(1, at = c(113,1113,2113,3113,4113), labels = c(200,400,600,800,1000), cex.axis = 1.5)
axis(2, at = c(251,1251,2251,3251,4251,5251), labels = c(200,400,600,800,1000,1200), cex.axis = 1.5)
text(x = 1080, y = 5000,"(c) Cumulative cost within knowledge-based window", cex = 1.75)
points(border_coords, cex = 0.3)


dev.off()
