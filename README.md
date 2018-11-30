phest - phenology estimation (in R) [![Build Status](https://travis-ci.org/willpearse/phest.svg?branch=master)](https://travis-ci.org/willpearse/phest)
===============================================================
Will Pearse (will.pearse@usu.edu)

A simple R package to calculate phenology metrics in R. It currently
has no non-base dependencies; to install, simply run:

```{r}
library(devtools)
install_github("willpearse/phest")
```

This contains an implementation of the new(ish) method outlined in
[Pearse, W. D., Davis, C. C., Inouye, D. W., Primack, R. B., & Davies,
T. J. (2017). A statistical estimator for determining the limits of
contemporary and historic phenology. Nature Ecology & Evolution,
1](https://www.nature.com/articles/s41559-017-0350-0). You can find
the details of this in the package's help under `?weibull`.
