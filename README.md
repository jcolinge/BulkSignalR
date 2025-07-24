
# BulkSignalR <img  width="120" height="139" src="man/figures/logo.png" align="right" />

<!-- badges: start -->
[![CRAN Version](https://www.r-pkg.org/badges/version/BulkSignalR)](https://cran.r-project.org/package=BulkSignalR)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/BulkSignalR)](https://cran.r-project.org/package=BulkSignalR)
<!-- badges: end -->

\

## Note about versions

**This repository contains our first version** of the library as published
in Nucleic Acids Research in 2023. Meanwhile, we added the library to
Bioconductor, which required slight modifications to fulfill Bioconductor
standards.

The current version of the library (the one in Bioconductor) can be found
[here](https://www.bioconductor.org/packages/release/bioc/html/BulkSignalR.html).
It is temporarily developed from another GitHub
[repository](https://github.com/ZheFrench/BulkSignalR) and it will ultimately
be integrated back into this one. Until then, this very repository should be
considered for the first version only.

\

## Overview

BulkSignalR is used to infer ligand-receptor (L-R) interactions from bulk
expression (transcriptomics/proteomics) data, or spatial
transcriptomics. Potential L-R interactions are taken from the
LR*db* database, which is  included in our other package SingleCellSignalR,
available from Bioconductor [here](https://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html).

Inferences rely on a statistical model linking potential
L-R interactions with biological pathways from Reactome or biological
processes from GO.

A number of visualization and data summary functions are proposed to
help navigating the predicted interactions.


<img   src="man/figures/workflow.png" align="center" width="85%" height="85%" />
  

## Installation

``` R

# BulkSignalR is not included in BioConductor yet.
# Installation goes via GitHub:
# install.packages("devtools")
devtools::install_github("jcolinge/BulkSignalR",build_vignettes = TRUE)

# To read the vignette
# browseVignettes("BulkSignalR")

```

## Notes

For a version history/change logs, see the [NEWS file](https://github.com/jcolinge/BulksignalR/blob/master/NEWS.md).


**BulkSignalR** has been successfully installed on Mac OS X, Linux, and Windows using R version 4.2.


The code in this repository is published with the [CeCILL](https://github.com/jcolinge/BulksignalR/blob/master/LICENSE.md) License.


<!-- badges: start -->
[![Generic badge](https://img.shields.io/badge/License-CeCILL-green.svg)](https://shields.io/)
<!-- badges: end -->



