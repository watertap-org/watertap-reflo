Getting Started
===============

The following sections provide information on how to setup your Conda environment and install WaterTAP-REFLO on your macOS, Linux, or Windows computer.

Using Conda environments
------------------------

Conda environments are a way to create and manage multiple sets of packages and/or Python versions on the same system without incurring conflicts. 
Each Conda environment is a dedicated directory, separate from other Conda environments and the operating system's own directories, containing 
its own collection of packages, executables, and Python installation, including the Python interpreter. Once a Conda environment is *activated*, 
using the ``conda activate`` command in a terminal/console, the environment's own version of Python will be used to run commands or interactive 
sessions for the remainder of the session. For these reasons, Conda environments are especially useful to install and manage multiple projects 
(and/or multiple *versions* of the same project) on the same computer with minimal effort, as they provide a way to seamlessly switch between 
different projects without conflicts.

Using Conda environments is not mandatory to be able to install and use WaterTAP; however, it is strongly recommended. 
To use Conda environments, the ``conda`` package manager is required. 
Refer to the `Conda installation guide <https://idaes-pse.readthedocs.io/en/stable/tutorials/getting_started/index.html#installation>`_ for detailed steps on how to install Conda for your operating system.

.. _install:


Install on Linux, Windows, and Mac (ARM)
----------------------------------------

This section is for Linux, Windows, and Apple computers with ARM processors (also known as Apple Silicon, e.g., M1, M2). 
To check if your Mac is using an ARM processor, click the apple icon in the upper left corner of your screen. 
Then select "About This Mac" to view the type of processor in your computer.

Create a Conda environment (in this example, named ``watertap-reflo``) where WaterTAP-REFLO and its runtime dependencies will be installed:

.. code-block:: shell

   conda create --name watertap-reflo --yes python=3.10

Activate the ``watertap-reflo`` environment using the command given below. 
If the environment was activated successfully, the environment's name will be displayed in the terminal prompt such as ``(watertap-reflo) project-directory $``.

.. code-block:: shell

   conda activate watertap-reflo

Install WaterTAP-REFLO in the Conda environment using ``pip``:

.. code-block:: shell

   pip install watertap-reflo

(Optional) See the :ref:`running-test-suite` section, if you want to verify that the installation was successful.

After installing WaterTAP-REFLO, the IDAES Extensions command can be used to automatically install the solvers distributed as part of the IDAES Extensions. 
Depending on your operating system, additional steps might be needed. 
For more information, refer to the `IDAES installation guide <https://idaes-pse.readthedocs.io/en/stable/tutorials/getting_started/index.html#installation>`_. 
From the same environment where WaterTAP-REFLO was installed, run:

.. code-block:: shell

   idaes get-extensions

.. important:: The ``conda activate`` command described above must be run each time a new terminal/console session is started.

.. note:: Typically, the ``idaes get-extensions`` command only needs to be run once for each system, as it will install the required files into a common, system-wide location.

Install on Mac (Intel)
----------------------

.. warning:: Intel Mac is currently not being tested, and should therefore be considered *unsupported*. The steps described here should not be expected to work without issues, and are here for historical and/or informational purposes only.

This section is for Apple's Mac computers with Intel processors. To check if your Mac is using an Intel processor, 
click the apple icon in the upper left corner of your screen. Then select About This Mac to view the type of processor in your computer.

Create a Conda environment (in this example, named ``watertap-reflo``) where WaterTAP-REFLO and its runtime dependencies will be installed:

.. code-block:: shell

   conda create --name watertap-reflo --yes python=3.10

Activate the ``watertap-reflo`` environment using the command given below. 
If the environment was activated successfully, the environment's name will be displayed in the terminal prompt such as ``(watertap-reflo) project-directory $``.

.. code-block:: shell

   conda activate watertap-reflo

Install WaterTAP-REFLO in the Conda environment using ``pip``:

.. code-block:: shell

   pip install watertap-reflo

(Optional) See the :ref:`running-test-suite` section, if you want to verify that the installation was successful.

After installing WaterTAP-REFLO, we need to ensure we have the Xcode toolkit, build the PyNumero Pyomo extensions, and obtain solvers from conda-forge. To install Xcode, run:

.. code-block:: shell

   xcode-select --install

To build PyNumero, from the same environment where WaterTAP-REFLO was installed, run the following commands:

.. code-block:: shell

   conda install --yes cmake
   pyomo build-extensions

The output of the second command should be something like:

.. code-block:: shell

   INFO: Finished building Pyomo extensions.
   INFO: The following extensions were built:
      [FAIL]  appsi
      [FAIL]  mcpp
      [ OK ]  pynumero

Next, we can obtain Ipopt and CBC from conda-forge:

.. code-block:: shell

   conda install --yes -c conda-forge ipopt coincbc

.. important:: The ``conda activate`` command described above must be run each time a new terminal/console session is started.

.. note:: The ``pyomo build-extensions`` command only needs to be run once for each system as it builds and installs the required libraries into a common, system-wide location. After building PyNumero, you should not need cmake. You can remove it by running ``conda uninstall cmake``.

.. _running-test-suite:

.. Running the test suite
.. ----------------------

.. To run the WaterTAP test suite, first install the ``pytest`` test framework:

.. .. code-block:: shell

..    pip install pytest

.. Then, run the following command to run the complete WaterTAP test suite:

.. .. code-block:: shell

..    pytest --pyargs watertap

.. (Optional) To see a list of available command-line options, run:

.. .. code-block:: shell

..    pytest --pyargs watertap --help

.. .. note:: Some tests will be skipped (denoted by an ``s`` symbol). This is to be expected, as some of the tests are only applicable within a developer environment.

.. _install-dev:


For WaterTAP-REFLO developers
-----------------------------

This section is for developers who plan to modify or contribute to REFLO's codebase. 
Contributing to REFLO will involve opening a Pull Request (PR) in REFLO's GitHub repository. For more information, refer to WaterTAP's documentation on 
`how to contribute to WaterTAP's development <https://watertap.readthedocs.io/en/stable/how_to_guides/how_to_contribute_development.html#developer-guide>`_. 
All of the same guidance presented there is applicable to REFLO.

Create a Conda environment (in this example, named ``watertap-reflo-dev``) where WaterTAP-REFLO and all dependencies needed for development will be installed, then activate it:

.. code-block:: shell

   conda create --name watertap-reflo-dev --yes python=3.11 && conda activate watertap-reflo-dev

Clone the WaterTAP-REFLO repository to your local development machine using ``git clone``, then enter the newly created ``watertap-reflo`` subdirectory:

.. code-block:: shell

   git clone https://github.com/watertap-org/watertap-reflo && cd watertap-reflo

Install WaterTAP-REFLO and the development dependencies using ``pip`` and the ``requirements-dev.txt`` file:

.. code-block:: shell

   pip install -r requirements-dev.txt

If needed, or if this is your first time installing IDAES, WaterTAP, or WaterTAP-REFLO on your machine, 
run the following line from the same environment where WaterTAP-REFLO was installed.

.. code-block:: shell

   idaes get-extensions

.. note:: Typically, the ``idaes get-extensions`` command only needs to be run once for each system, as it will install the required files into a common, system-wide location.  Depending on your operating system, you may need to follow additional steps described above to install solvers distributed through IDAES Extensions.
   
.. (Optional but recommended) `Pre-commit hooks <https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks>`_ are scripts that are automatically run by Git "client-side" (i.e. on a developer's local machine) whenever `git commit` is run. WaterTAP uses the `pre-commit <https://pre-commit.com/>`_ framework to manage a few hooks that are useful for WaterTAP developers. To install the WaterTAP pre-commit hooks, run:

.. .. code-block:: shell

..    pre-commit install

To verify that the installation was successful, try running the WaterTAP-REFLO test suite using ``pytest``:

.. code-block:: shell

   pytest

.. To view/change the generated documentation, see the :ref:`documentation-mini-guide` section.

.. Installation
.. ------------

.. To install **WaterTAP-REFLO**, run:

.. .. code-block:: shell

..     git clone https://github.com/watertap-org/watertap-reflo && cd watertap-reflo
..     conda create --yes --name watertap-reflo-dev python=3.10 && conda activate watertap-reflo-dev
..     pip install -r requirements-dev.txt

.. Running tests
.. -------------

.. .. code-block:: shell
    
..     conda activate watertap-reflo-dev
..     pytest --pyargs watertap_contrib.reflo

.. Formatting code
.. ---------------

.. Before committing, the Python code must be formatted with `Black <https://black.readthedocs.io>`_.

.. Black is installed by default as part of the developer dependencies. To format the code, run the following command from the local repository root directory:

.. .. code-block:: shell
    
..     conda activate watertap-reflo-dev
..     black .

