# Regression Analysis of Group-Tested Current Status Data under Semiparametric Probit Model

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Description

This project includes the `Bayes.Probit.GT` function, which serves multiple purposes:

- Identifies significant risk factors and assesses their effects.
- Estimates the cumulative incidence functions for different subgroups under a semiparametric Probit model.
- Can be used for various group testing methodologies such as master pooling, Dorfman testing, and array testing.
- For Dorfman and array testing, it can also estimate the unknown array accuracies.

Basic simulation code is available in the `Simulate.R` file.

## Table of Contents

## Table of Contents

- [Installation](#installation)
  - [Download the Repository](#1-download-the-repository)
  - [Check Package Dependencies](#2-check-package-dependencies)
  - [Prepare Input Data](#3-prepare-input-data)
  - [Run Basic Examples](#4-run-basic-examples)
- [Usage](#usage)
  - [Source Necessary Files](#1-source-all-necessary-files)
  - [Input Formatting](#2-input-formatting)
  - [Arguments](#3-arguments)
  - [Output](#4-output)
- [Example](#example)
- [Contributing](#contributing)
- [License](#license)
- [Authors](#authors)
- [Acknowledgments](#acknowledgments)

## Installation

### 1. Download the Repository

- Click the green "Code" button at the top of the GitHub repository.
- Select "Download ZIP" to download the repository files to your local machine.
- Extract the contents of the ZIP file to a directory of your choice.

### 2. Check Package Dependencies

Before using the `Bayes.Probit.GT` function, ensure that the following R packages are installed:

```R
install.packages(c("mvtnorm", "truncnorm", "Rcpp"))
```
### 3. Prepare Input Data

The `Bayes.Probit.GT` function requires specific input data formats, which can be generated using functions provided in the Testing Functions.txt file.
Refer to the descriptions provided in the Testing Functions.txt file for details on each model and how to prepare the input data.

### 4. Run Basic Examples

Explore basic examples of how to use the `Bayes.Probit.GT` function by sourcing the necessary files and running simple simulations for each model.
Refer to the **Simulate.R** for a starting point and modify the code as needed for your specific use case.

By following these steps, you'll be ready to use the `Bayes.Probit.GT` package for Bayesian Probit regression for group testing data analysis.

## Usage

To use the Bayes.Probit.GT function, follow these steps:

### 1. Source all necessary files:

```R
Rcpp::sourceCpp("SampLatent.cpp")
source("Testing Functions.txt")
source("Bayes.Probit.GT.R")
```
### 2. Input formatting:

The input for the Bayes.Probit.GT function requires specific formatting, which can be generated using functions provided in the **Testing Functions.txt** file. Each model is explained in detail in the text file, so be sure to read it carefully. You can also read 'Simulate.R' file.

### 3. Arguments:

Define the arguments for the Bayes.Probit.GT function:
- Z: A matrix of testing responses. Each row represents a test, with columns indicating the individual's ID, the number of individuals in the test, the assay used, and the indices of the individuals assigned to the test pools. * can be produced by using functions in **Testing Functions.txt** file.
- X: Covariate matrix containing covariate information for each individual.
- Y: Matrix indicating the pools each individual was assigned to.  * can be produced by using functions in **Testing Functions.txt** file.
- c: Censoring or testing time for each individual.
- grid: Grid definition for calculating the baseline survival function. Default is NULL.
- n.grid: Length of the grid.
- init.theta: Initial value of theta (Default is 0 vector)
- eta: Initial value of hyper parameter eta (Default is 1)
- gam0: Initial value of gam0 for the splines (Default is -3)
- gam: Initial values of the spline coefficients (Defalut is rep(0.1, m+order))
- theta0, Sigma0, m0, v0, a0, b0, ae, be, ap, bp: Priors for the model. Default number is given.
- Se: Vector of sensitivity values, if known.
- Sp: Vector of specificity values, if known.
- order: Order for I splines (usually 3 or 4).
- knots: Interior knots for the spline functions. Default is NULL.
- m: Number of interior knots.
- quantile: If TRUE, knots are created based on quantiles. If FALSE, equally spaced knots are created.
- maxiter: Maximum number of iterations.
- burn.in: Burn-in period.
- na: Number of arrays.
- err.est: Set to TRUE if assay accuracies are unknown.

Be sure to use these arguments appropriately when calling the function.


### 4. Output:
- **theta.mat:** Matrix containing the theta chain.
- **g.mat:** Unknown function g(t) based on the grid.
- **grid:** Indicates the grid used for estimation.
- **St:** Estimated baseline survival function based on the provided grid.
- **summary.theta:** Table summarizing theta, including posterior mean, 2.5% quantile, and 97.5% quantile. If `err.est = TRUE`, it also provides information about array accuracy estimations.

### Examples

Here are some example codes from the `Simulate.R` file:

#### Simulate Dorfman Group Testing

```R
# Simulate Dorfman group testing
Dorf.Data = Dorfman.decode.diff.error(Y.true, Se.true, Sp.true, d)
Z2 = Dorf.Data$Z
Y2 = Dorf.Data$Y
```
The above code is for input formatting in the Dorfman testing scenario.

```R
df = Bayes.Probit.GT(Z2, X, Y2, c, Se = c(0.95, 0.98), Sp = c(0.98, 0.99), na = 2)
df$summary.theta
```
The above code is when using Bayes.Probit.GT, which is the main function. This is when assay accuracies are known.

```R
# Unknown assay accuracies
df.unknown = Bayes.Probit.GT(Z2, X, Y2, c, na = 2, err.est = TRUE)
df.unknown$summary.theta
```
This is when the assay accuracies are unknown.


## Contributing

Thank you for considering contributing to the project! Whether you're reporting a bug, proposing a new feature, or submitting a pull request, your contributions are welcome.

### How to Contribute

1. **Reporting Issues:** If you encounter any bugs or issues with the project, please [open a new issue](https://github.com/jihyunk1114/Bayes.Probit.GT/issues) on GitHub. Include detailed information about the problem and steps to reproduce it.

2. **Submitting Pull Requests:** If you'd like to contribute code changes or new features, feel free to submit a pull request. Before doing so, please make sure to fork the repository, create a new branch for your changes, and write clear commit messages.

3. **Feature Requests:** If you have ideas for new features or improvements, you can also [submit a feature request](https://github.com/jihyunk1114/Bayes.Probit.GT/issues) on GitHub. Provide as much detail as possible about the proposed feature and why it would be valuable.

Thank you for helping to improve the project!

## License

This project is licensed under the MIT License.

## Authors

Jihyun Kim: [jihyunk@email.sc.edu](mailto:jihyunk@email.sc.edu)

## Acknowledgments

This project utilizes several external packages and resources. We acknowledge the following:

- [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
- [mvtnorm](https://cran.r-project.org/web/packages/mvtnorm/index.html)
- [truncnorm](https://cran.r-project.org/web/packages/truncnorm/index.html)
- Dr. McMahan's C++ code containing the `SampLatent` function.

We are grateful for the contributions of these packages and resources to our project.
