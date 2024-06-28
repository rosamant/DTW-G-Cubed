# DTW-G-Cubed

## Code supporting the manuscript "Enhancing the Accuracy of Dynamic Time Warping by Integrating Stratigraphic Constraints into the Automated Correlation of Sedimentary Sequences"

Natural Gamma Radiation data used for this paper is stored in .csv files and is freely available through the National Offshore Petroleum Information Management System (NOPIMS, Geoscience Australia) (https://www.ga.gov.au/nopims).

The repository contains the R scripts for the DTW calculations performed and presented for the G-Cubed paper on the methodology of DTW custom technique. 

Figure 4 can be reproduced using the `Picard1-U1464-StepPatterns.R` file. This file contains the DTW calculations using 4 different step patterns, followed by the codes to create the plot.

The `Picard1-U1464-WindowSize.R` file contains codes to create Figure 5 and Supplementary Figures 1 & 2. Figure 5 is created using window size = 2000. Supplementary Figures 1 & 2 are created with window sizes - 500 & 1000. Standard `Asymmetric P1` step pattern is used to perform DTW and the codes to plot individual figures are provided in the same file.

To reproduce Supplementary Figures 3 to 6, use files `Picard1-Finucane1-Preliminary.R`, `Picard1-Angel2-Preliminary.R`, `Picard1-Goodwyn2-Preliminary.R`, and `Picard1-Goodwyn6-Preliminary.R`. DTW calculation and corresponding plot are available in the same R scripts.

To reproduce Figures 8 to 11, run the `Custom_Step_Pattern.R` file first. To plot Figure 8, first use the `Picard1-Angel2-DTW.R` file to compute DTW with the standard step pattern and then with the customized 'Stratigraphy-optimized' step pattern and 'Knowledge-based' windowing function. Individual plots can be created using the codes within the same file. To reproduce Figure 8, use `Picard1-Angel2-Plot.R`.

Follow the same procedure for Angel-2 and Minilya-1 sites.

To create Figure 11, run the `Picard1-Angel2-Minilya1-Whitetail1-Plot.R` file to create the plot.

To create Figure 3 of the standard and knowledge-based window, make sure to run the entire script of `Picard1-Minilya1-DTW.R` file. This performs DTW with customized step pattern and windowing function. Also, run the `CompareWindowBorders` to create the customized window borders. Then, run the `CustomWindow.R` file.

## License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

## Dependencies

- Giorgino T (2009). “Computing and Visualizing Dynamic Time Warping Alignments in
  R: The dtw Package.” _Journal of Statistical Software_, *31*(7), 1-24.
  doi:10.18637/jss.v031.i07 <https://doi.org/10.18637/jss.v031.i07>.
  
- Meyers, S.R. (2014). Astrochron: An R Package for Astrochronology.
  https://cran.r-project.org/package=astrochron

- Signorell A (2024). _DescTools: Tools for Descriptive Statistics_. R package
  version 0.99.54, <https://CRAN.R-project.org/package=DescTools>.
