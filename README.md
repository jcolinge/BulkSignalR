
# BulkSignalR <img  width="120" height="139" src="man/figures/logo.png" align="right" />

<!-- badges: start -->
[![CRAN Version](https://www.r-pkg.org/badges/version/BulkSignalR)](https://cran.r-project.org/package=BulkSignalR)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/BulkSignalR)](https://cran.r-project.org/package=BulkSignalR)
<!-- badges: end -->

## Overview

Inference of ligand-receptor (LR) interactions from `bulk`
(transcriptomic or proteomic) data or `spatial transcriptomics`.
 **BulkSignalR** bases its inferences
on the LRdb database included in our other package, **SingleCellSignalR**
available from Bioconductor [here](https://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html). It relies on a statistical model that
is specific to bulk data sets. Different visualization and data
summary functions are proposed to help navigating prediction results.
</br>
<img   src="man/figures/workflow.png" align="center" width="85%" height="85%" />
</br>

## Installation

``` R

# The easiest way to get BulkSignalR is to install :
# Not deployed under Bioconductor yet.
# install.packages("BiocManager")
# BiocManager::install("BulkSignalR")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("jcolinge/BulkSignalR")

```

## Notes

For a version history/changelog, please see the [NEWS file](https://github.com/zhefrench/BulksignalR/blob/master/NEWS.md).


**BulkSignalR** has been successfully installed on Mac OS X, Linux, and Windows using R 4.2 version.

<!-- badges: start -->
[![Linux](https://svgshare.com/i/Zhy.svg)](https://svgshare.com/i/Zhy.svg)
[![macOS](https://svgshare.com/i/ZjP.svg)](https://svgshare.com/i/ZjP.svg)
[![Windows](https://svgshare.com/i/ZhY.svg)](https://svgshare.com/i/ZhY.svg)
<!-- badges: end -->


The code in this repository is published with the [CeCILL](https://github.com/zhefrench/BulksignalR/blob/master/LICENSE.md) License.


<!-- badges: start -->
[![Generic badge](https://img.shields.io/badge/License-CeCILL-green.svg)](https://shields.io/)
<!-- badges: end -->



