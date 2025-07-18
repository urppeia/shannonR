# ShannonR

[![Website](https://github.com/urppeia/shannonR/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/urppeia/shannonR/actions/workflows/pkgdown.yaml)
[![Codecov test coverage](https://codecov.io/gh/urppeia/shannonR/branch/main/graph/badge.svg)](https://codecov.io/gh/urppeia/shannonR)
[![License](https://img.shields.io/github/license/urppeia/shannonR)](https://github.com/urppeia/shannonR/blob/main/LICENSE.md)
[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://urppeia.github.io/shannonR/)
[![Last Commit](https://img.shields.io/github/last-commit/urppeia/shannonR)](https://github.com/urppeia/shannonR/commits/main)
[![GitHub issues](https://img.shields.io/github/issues/urppeia/shannonR)](https://github.com/urppeia/shannonR/issues)
[![Dependencies](https://img.shields.io/badge/dependencies-analyzed-brightgreen.svg)](https://github.com/urppeia/shannonR/blob/main/DESCRIPTION)

## Overview
`ShannonR` is an `R` package designed for advanced analysis of 
methylation data, developed as part of the **URPP Evolution in Action**: 
_From Genomes to Ecosystems_ initiative. This package provides tools 
for processing, analyzing, and visualizing methylation data, 
with a focus on biological and evolutionary research applications.

## Features

**Data Processing:** Efficiently handle large methylation datasets with
robust preprocessing functions.
**Statistical Analysis:** Perform comprehensive statistical analyses 
tailored for methylation data.
**Visualization:** Generate high-quality plots and visualizations to 
interpret methylation patterns.
**Integration:** Seamlessly integrate with other bioinformatics tools 
and pipelines.

## Installation
To install `ShannonR`, you can use the following commands in R:

```
# Install devtools if not already installed
install.packages("devtools")

# Install ShannonR from GitHub
devtools::install_github("urppeia/shannonR")
```

## Usage
Here’s a basic example to get started with ShannonR:

```
# Load the library
library(shannonR)
```

## Using `renv`

1. Clone the repo locally
```
git clone git@github.com:urppeia/shannonR.git
```

2. Open the project by opening `shannonR.Rproj` file

3. Restore packages
```
install.packages("renv")
renv::restore()
```

## Contributing
Contributions are welcome! Please follow these steps:

1. `Fork` the repository.
2. Create a new `branch` (`git checkout -b feature-branch`).
3. `Commit` your changes (`git commit -m 'Add new feature'`).
4. `Push` to the branch (`git push origin feature-branch`).
5. Open a `Pull Request`.

Please ensure your code adheres to the Contributor Covenant Code of Conduct.

## License
This project is licensed under the MIT License. See the LICENSE file for details.
