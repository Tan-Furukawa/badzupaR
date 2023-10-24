# badzupaR (Beta Version)

[![DOI](https://zenodo.org/badge/497178337.svg)](https://zenodo.org/badge/latestdoi/497178337)

# installation
install badzupaR by your R console (https://github.com/Tan-Furukawa/badzupaR)

```
install.packages('remotes')
remotes::install_github('Tan-Furukawa/badzupaR')
```

# example

```
library(xbadzupaR)
# put your vector age data to dat
set.seed(1234); dat <- c(rnorm(200, 4), rnorm(50, 8, 2))
xbadzupa(dat)
```