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

Instructions for installing or setting up the project. Include any dependencies and how to install them.

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
