# DTW-G-Cubed
The repository contains the R scripts for the DTW calculations performed and presented for the G-Cubed paper on the methodology of DTW custom technique. 

Figure 4 can be reproduced using the 'Picard1-U1464-StepPatterns.R' file. This file contains the DTW calculations using 4 different step patterns, followed by the codes to create the plot.

The 'Picard1-U1464-WindowSize.R' file contains codes to create Figure 5 and Supplementary Figures 1 & 2. Figure 5 is created using window size = 2000. Supplementary Figures 1 & 2 are created with window sizes - 500 & 1000. Standard 'Asymmetric P1' step pattern is used to perform DTW and the codes to plot individual figures are provided in the same file.

To reproduce Supplementary Figures 3 to 6, use files 'Picard1-Finucane1-Preliminary.R', 'Picard1-Angel2-Preliminary.R', 'Picard1-Goodwyn2-Preliminary.R', and 'Picard1-Goodwyn6-Preliminary.R'. DTW calculation and corresponding plot are available in the same R scripts.

To reproduce Figures 8 to 11, run the 'Custom_Step_Pattern.R' file first. To plot Figure 8, first use the 'Picard1-Angel2-DTW.R' file to compute DTW with the standard step pattern and then with the customized 'Stratigraphy-optimized' step pattern and 'Knowledge-based' windowing function. Individual plots can be created using the codes within the same file. To reproduce Figure 8, use 'Picard1-Angel2-Plot.R'.

Follow the same procedure for Angel-2 and Minilya-1 sites.

To create Figure 11, run the 3 scripts of 'Picard1-Angel2-DTW.R', 'Picard1-Minilya1-DTW.R', and 'Picard1-Whitetail1-DTW.R'. Then run the 'Picard1-Angel2-Minilya1-Whitetail1-Plot.R' file to create the plot.

To create Figure 3 of the standard and knowledge-based window, make sure to run the entire script of 'Picard1-Minilya1-DTW.R' file. This performs DTW with customized step pattern and windowing function. Also, run the 'CompareWindowBorders' to create the customized window borders. Then, run the 'CustomWindow.R' file.


