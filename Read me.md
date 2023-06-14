# Read me

This document describe the sripts used in Perdoux and al. for statistacal analysis and plots.

## Root growth

4 scripts are availlable for root growth analyses

- Dose_effect_1line.R

This script build a curve for root growth of one line (WT) with one inhibitor. Growth is expressed in percent of the growth on DMSO.

- Root_growth.R

This script take the growth data for many lines and one inhibitor. It performs the statistical analysis according to the normality status. On normally distributed data we performed analysis of variance with post hoc Tukey test and on nonparametric data we performed the Kruskal–Wallis test with post hoc Dunn test from the rstatix package ([Kassambara, 2021](https://elifesciences.org/articles/75041#bib34)). 

Two Graphs are produced using the package ggplot2 ([Wickham, 2009](https://elifesciences.org/articles/75041#bib74)) , one for the root growth over the time and one to compare root size at a defined day. Growth is expressed in percent of the growth on DMSO.

- Az_effect_barplot.R

This script take the growth data for many lines on DMSO and AZ. It performs the statistical analysis according to the normality status. On normally distributed data we performed analysis of variance with post hoc Tukey test and on nonparametric data we performed the Kruskal–Wallis test with post hoc Dunn test from the rstatix package ([Kassambara, 2021](https://elifesciences.org/articles/75041#bib34)). Growth is expressed in percent of the growth on DMSO.

A barplot is produced using the package ggplot2 ([Wickham, 2009](https://elifesciences.org/articles/75041#bib74))  to compare root size of all lines at a defined day.

- Az_effect_growthcurve.R

This script take the growth data for many lines on DMSO and AZ. It performs the statistical analysis according to the normality status. On normally distributed data we performed analysis of variance with post hoc Tukey test and on nonparametric data we performed the Kruskal–Wallis test with post hoc Dunn test from the rstatix package ([Kassambara, 2021](https://elifesciences.org/articles/75041#bib34)). 

A growth curve is produced using the package ggplot2 ([Wickham, 2009](https://elifesciences.org/articles/75041#bib74))  to compare root size of all lines with incresing concentration of AZ. Growth is expressed in percent of the growth on DMSO.

## FV/FM

2 scripts are availlable for Fv/FM analyses

- dose_effect_Fv-Fm.R

This script take the Fv/Fm data for many lines on DMSO and AZ. It performs the statistical analysis according to the normality status. On normally distributed data we performed analysis of variance with post hoc Tukey test and on nonparametric data we performed the Kruskal–Wallis test with post hoc Dunn test from the rstatix package ([Kassambara, 2021](https://elifesciences.org/articles/75041#bib34)). 

A  curve is produced using the package ggplot2 ([Wickham, 2009](https://elifesciences.org/articles/75041#bib74))  to compare root size of all lines with incresing concentration of AZ. Growth is expressed in percent of the growth on DMSO.

- Timecourse_Fv-Fm.R

This script take the Fv/Fm data for many lines on DMSO and AZ. It performs the statistical analysis according to the normality status. On normally distributed data we performed analysis of variance with post hoc Tukey test and on nonparametric data we performed the Kruskal–Wallis test with post hoc Dunn test from the rstatix package ([Kassambara, 2021](https://elifesciences.org/articles/75041#bib34)). 

A  curve is produced using the package ggplot2 ([Wickham, 2009](https://elifesciences.org/articles/75041#bib74))  to compare root size of all lines over the time. Growth is expressed in percent of the growth on DMSO.







