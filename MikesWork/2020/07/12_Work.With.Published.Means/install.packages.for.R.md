# Install Packages I Need 
- Run this coded on a local machine and it will ensure I've got all the packages I need for these analyses
- Was originally in 'load.and.process.means.data.Rmd' prior to 12 Jul 2020

```{r}
if (!requireNamespace("readr", quietly = TRUE))
    install.packages("readr")
if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")
if (!requireNamespace("forcats", quietly = TRUE))
    install.packages("forcats")
## Distributions Hermite Polynomial Approximation: used to calculate moments for truncated distributions. using truncatedNormalMoment
## We actually have censored data, where values exist but their exact value is unknown beyond the fact they are less than some threshold, not truncated data, which ignores data outside of a range.
##if (!requireNamespace("hpa", quietly = TRUE))
##    install.packages("hpa")
if (!requireNamespace("EnvStats", quietly = TRUE))
    install.packages("EnvStats")
## eivreg
if (!requireNamespace("eivtools", quietly = TRUE))
    install.packages("eivtools")
if (!requireNamespace("Biostrings", quietly = TRUE))
    BiocManager::install("Biostrings", lib="/home/mikeg/R/x86_64-pc-linux-gnu-library/3.6")
if (!requireNamespace("rmarkdown", quietly = TRUE))
    install.packages("rmarkdown")
if (!requireNamespace("hash", quietly = TRUE))
    install.packages("hash")
if (!requireNamespace("purrr", quietly = TRUE))
    install.packages("purrr")
```
