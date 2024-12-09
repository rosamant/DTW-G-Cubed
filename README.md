# DTW-G-Cubed

## Code supporting the manuscript "Enhancing the Accuracy of Dynamic Time Warping by Integrating Stratigraphic Constraints into the Automated Correlation of Sedimentary Sequences"

Natural Gamma Radiation data used for this paper is stored in .csv files and is freely available through the National Offshore Petroleum Information Management System (NOPIMS, Geoscience Australia) (https://www.ga.gov.au/nopims).

The repository contains the R scripts for the DTW calculations performed and presented for the Paleoceanography and Paleoclimatology paper on the methodology of DTW custom technique. 

Figure 6 can be reproduced using the `U1463_U1464_DTW.R` file. This file contains the DTW calculations using 5 different step patterns and 3 windowing functions, followed by the codes to create the plot. The script contains the calculation of RMSE and the codes to reproduce Supplementary Figure 1.

Figure 7 can be reproduced using the `U1464_U1482_DTW.R` file. This file contains the DTW calculations using 5 different step patterns and 3 windowing functions, followed by the codes to create the plot. The script contains the calculation of RMSE and the codes to reproduce Supplementary Figure 3.

To reproduce Figures 8 to 11, run the `Custom_Step_Pattern.R` file first. To plot Figure 8, first use the `Picard1-Angel2-DTW.R` file to compute DTW with the standard step pattern and then with the customized 'Stratigraphy-optimized' step pattern and 'Knowledge-based' windowing function. Individual plots can be created using the codes within the same file. Before running the density plot, run the `CompareWindowBorders` to create the customized window borders. To reproduce Figure 8, use `Picard1-Angel2-Plot.R`.

Follow the same procedure for Angel-2 and Minilya-1 sites.

To create Figure 11, run the `Picard1-Angel2-Minilya1-Whitetail1-Plot.R` file to create the plot.

To create Figure 4 of the standard and knowledge-based window, make sure to run the entire script of `Picard1-Minilya1-DTW.R` file. This performs DTW with customized step pattern and windowing function. Also, run the `CompareWindowBorders` to create the customized window borders. Then, run the `CustomWindow.R` file.

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

- Hagen, C., Creveling, J., & Huybers, P. (2024). Align: A User-Friendly App for Numerical Stratigraphic Correlation. GSA Today, 34(2). https://doi.org/10.1130/GSATG575A.1
