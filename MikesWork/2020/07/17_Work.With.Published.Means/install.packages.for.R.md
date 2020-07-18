# Install Packages I Need 
- Run this coded on a local machine and it will ensure I've got all the packages I need for these analyses
- Was originally in 'load.and.process.means.data.Rmd' prior to 12 Jul 2020

```{r}

funcList<-
    c("readr",
      "ggplot2",
      "forcats", 
      "EnvStats",
      "eivtools",
      "Biostrings",
      "rmarkdown",
      "knitr",
      "hash",
      "purrr",
      "optimx",
      "ggpubr",
      "ggpmisc"
      )

tmp <- lapply(
    funcList,
    function(x) {if (!requireNamespace(x, quietly = TRUE)) install.packages(x, lib=Sys.getenv("R_LIBS_USER"))  }
    )
rm(tmp)
```
