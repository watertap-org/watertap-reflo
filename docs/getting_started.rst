Getting Started
==============

These instructions are currently intended for contributors to WaterTAP-REFLO

**WaterTAP-REFLO** supports Python versions 3.8 through 3.10.

Prerequisites
-------------

- The conda package and environment manager, for example by using the `Miniconda installer <https://docs.conda.io/en/latest/miniconda.html#miniconda>`_ following the steps appropriate for your operating system

Installation
------------

To install **WaterTAP-REFLO**, run:

.. code-block:: shell

    git clone https://github.com/watertap-org/watertap-seto && cd watertap-seto
    conda create --yes --name watertap-seto-dev-env python=3.10 && conda activate watertap-seto-dev-env
    pip install -r requirements-dev.txt

Running tests
-------------

.. code-block:: shell
    conda activate watertap-seto-dev-env
    pytest --pyargs watertap_contrib.seto

Formatting code
---------------

Before committing, the Python code must be formatted with `Black <https://black.readthedocs.io>`_.

Black is installed by default as part of the developer dependencies. To format the code, run the following command from the local repository root directory:

.. code-block:: shell
    conda activate watertap-seto-dev-env
    black .

