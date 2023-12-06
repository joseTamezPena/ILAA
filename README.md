# ILAA and the ERT

This repository showcases the use of the **FRESA.CAD::ILAA(), FRESA.CAD::IDeA()** and the **FRESA.CAD::getLatentCoefficients()** for the discovery of latent variables from tabular data-sets.

## Table of Contents

-   [About](#about)
-   [Installation](#installation)
-   [Usage](#usage)
-   [Contributing](#contributing)
-   [License](#license)

## About {#about}

Examples of usage and applications of the **FRESA.CAD::ILAA()** and the associated exploratory residualization matrix (**ERT**).

This repository also holds the Source code of the Shiny App: [ERT Calculator](https://josetamezpena.shinyapps.io/ILAA/) the source code is at the ILAA folder.

The repository folder structure is:

+-------------+------------------------------------------------------------------------+
| Folder Name | Contents                                                               |
+=============+========================================================================+
| Data        | Data sets used in the examples of this repository                      |
|             |                                                                        |
|             | The data sources are described in: DataSpecsAndSource.xlsx             |
+-------------+------------------------------------------------------------------------+
| ILAA        | Shiny App code                                                         |
+-------------+------------------------------------------------------------------------+
| Main        | The ILAA Tutorial Code                                                 |
+-------------+------------------------------------------------------------------------+
| RMD         | The RMD scripts used for the validation and showcasing the ILAA method |
+-------------+------------------------------------------------------------------------+

: Repository Structure

## Installation {#installation}

ILAA is part of the FRESA.CAD package. To use ILAA first install FRESA.CAD You can install the official release of the package from CRAN using:

``` r
install.packages("FRESA.CAD")
```

To install the development version from GitHub, use:

``` r
# Install 'devtools' package if you haven't already
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("https://github.com/joseTamezPena/FRESA.CAD")
```

## Usage {#usage}

After installation you can test ILAA on the iris data set.

```         
library("FRESA.CAD") 
# The IRIS dataset
data('iris')  
##FCA Decorrelation at 0.25 threshold, pearson and fast estimation  
irisDecor <- IDeA(iris,thr=0.25)  
# Print the latent variables 
print(getLatentCoefficients(irisDecor))

# Lets model setosa using logistic regression
setosaData <- iris
setosaData$setosa <- 1.0*(as.character(iris$Species)=="setosa")
setosaData$Species <- NULL
setosaILAA <- ILAA(setosaData,thr=0.25)
modelSetosa <- glm(setosa~.,setosaILAA,family="binomial")
#The model cofficients in the ERT space
print(modelSetosa$coefficients)
#Get the observed Coefficients
observedCoef <- getObservedCoef(setosaILAA,modelSetosa)
#The model cofficients in the Observed space
print(observedCoef)
```

You can test ILAA using the ERT calculator at: <https://josetamezpena.shinyapps.io/ILAA/>

The app will compute the ERT transformation using the user-provided data set.

Also you can look at the output of the ILAA tutorial at: <https://rpubs.com/J_Tamez/ILAA_Tutorial>

## Contributing {#contributing}

Contributions are welcome! If you'd like to contribute to this project, please follow these guidelines:

-   Fork the repository.
-   Create a new branch: `git checkout -b feature/new-feature`.
-   Make your changes and commit them: `git commit -m 'Add new feature'`.
-   Push to the branch: `git push origin feature/new-feature`.
-   Submit a pull request.

## License {#license}

This project is licensed under the LGPL Licence 3.0 see the [LICENSE](LICENSE) file for details.

## Contact

Email: jose.tamezpena\@tec.mx

Twitter: [\@tamezpena](https://twitter.com/jtamezpena)
