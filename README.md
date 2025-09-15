# WaterTAP-REFLO
## Renewable Energy and Flexible Load Optimization

Welcome to the code repository for **WaterTAP-REFLO**! 
WaterTAP-REFLO is an extension of the [WaterTAP package](https://watertap.readthedocs.io/en/stable/).

![GitHub issues](https://img.shields.io/github/issues/watertap-org/watertap-reflo)
![GitHub pull requests](https://img.shields.io/github/issues-pr/watertap-org/watertap-reflo)
<!-- ![CI status](https://img.shields.io/github/actions/workflow/status/watertap-org/watertap-reflo/test.yml?branch=main) -->

## Documentation

See the full [WaterTAP-REFLO documentation](<https://watertap-reflo.readthedocs.io>) for more information.

## Getting Started

**WaterTAP-REFLO** supports Python versions 3.9 through 3.11.

### Prerequisites

- The conda package and environment manager, for example by using the [Miniconda installer](https://docs.conda.io/en/latest/miniconda.html#miniconda) following the steps appropriate for your operating system

### Installation

To install **WaterTAP-REFLO**, run:

```sh
conda create --yes --name watertap-reflo python=3.10 && conda activate watertap-reflo
pip install watertap-reflo
```

To install **WaterTAP-REFLO** *for development*, run:

```sh
git clone https://github.com/watertap-org/watertap-reflo && cd watertap-reflo
conda create --yes --name watertap-reflo python=3.10 && conda activate watertap-reflo
pip install -r requirements-dev.txt
```

### Running tests

```sh
conda activate watertap-reflo
pytest --pyargs watertap_contrib.reflo
```

### Formatting code

Before committing, the Python code must be formatted with [Black](https://black.readthedocs.io).

Black is installed by default as part of the developer dependencies. To format the code, run the following command from the local repository root directory:

```sh
conda activate watertap-reflo
black .
```

