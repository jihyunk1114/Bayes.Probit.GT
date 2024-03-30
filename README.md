# Regression Analysis of Group-Tested Current Status Data under Semiparametric Probit Model

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Description

This project includes the `Bayes.Probit.GT` function, which serves multiple purposes:

- Identifies significant risk factors and assesses their effects.
- Estimates the cumulative incidence functions for different subgroups under a semiparametric Probit model.
- Can be used for various group testing methodologies such as master pooling, Dorfman testing, and array testing.
- For Dorfman and array testing, it can also estimate the unknown array accuracies.

Basic simulation code is available in the `sample.R` file.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
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
## 3. Prepare Input Data

The `Bayes.Probit.GT` function requires specific input data formats, which can be generated using functions provided in the Testing Functions.txt file.
Refer to the descriptions provided in the Testing Functions.txt file for details on each model and how to prepare the input data.

## 4. Read Documentation

Familiarize yourself with the functionality and usage of the `Bayes.Probit.GT` function by reading the documentation provided in the Testing Functions.txt file.

## 5. Run Basic Examples

Explore basic examples of how to use the `Bayes.Probit.GT` function by sourcing the necessary files and running simple simulations for each model.
Refer to the Simulate.R file for a starting point and modify the code as needed for your specific use case.

By following these steps, you'll be ready to use the `Bayes.Probit.GT` package for Bayesian Probit regression for group testing data analysis.

## Usage

Instructions on how to use the project, including examples if applicable.

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
