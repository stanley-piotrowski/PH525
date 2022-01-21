# PH525.1-- Statistics with R

This repository contains data, notes, and exercise solutions for HarvardX's [Statistics and R](https://www.edx.org/course/statistics-and-r), the first in the series of the [Data Analysis for the Life Sciences series](https://www.edx.org/professional-certificate/harvardx-data-analysis-for-life-sciences).  All data used in the course can be found in the `data` directory; notes and exercise solutions for each week, written primarily in the tidyverse style, can be found in the `scripts` directory.  A brief description of the materials covered each week and the content of the associated scripts can be found below.  

**Week 1:** introduction to R and exploratory data analysis (EDA)

  * Exercise solutions to wrangle data, compute summary statistics, and build a variety of EDA plots including QQ-plots, boxplots, and scatterplots. 
  
**Week 2:** random variables, probability distributions, and the central limit theorem

  * Simulating null hypotheses through repeated sampling using Monte Carlo simulations.
  * Exploring the central limit theorem using Monte Carlo simulations and QQ-plots.  
  * Verifying assumptions of t-tests on real and simulated data.
  * Using simulations and QQ-plots to explore how the normal approximation holds for data with varying sample sizes.
  
**Week 3:** statistical inference

  * Exploring the concept of confidence intervals with standard normal and t-distributions through Monte Carlo simulations using a data set with known population parameters.
  * Exploring the concept of statistical power by varying sample sizes and Type I error rates using Monte Carlo simulations.  
  * Verifying assumptions of normality with standard normal and t-distributions using QQ-plots.
  * Exploring the concept of permutation tests using custom functions when parametric approaches are inappropriate.
  * Performing association tests using 2x2 contingency tables and calculating exact p-values with Fisher's exact test and approximations using the chi-square test.  
  
**Week 4:** EDA and robust summaries

  * Exploring useful EDA plots covered previously (e.g., scatterplots, QQ-plots, barplots, and boxplots) and plots to avoid (e.g., donut plots, 3D plots, etc).
  * Comparing the influence of outliers on statistcal summaries like the mean and standard deviation to robust statistics like the median and median absolute deviation (MAD).  
  * Using Wilcoxon rank-sum tests when data contain outliers that could skew results of parametric approaches (e.g., t-tests are sensitive to outliers due to their influence on the mean and standard deviation). 
