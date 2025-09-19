.. index::
   pair: watertap-reflo;Home

Welcome to WaterTAP-REFLO's documentation!
==========================================

WaterTAP-REFLO (Water treatment Technoeconomic Assessment Platform with Renewable Energy and Flexible Load Optimization) is an extension of 
`WaterTAP <https://watertap.readthedocs.io/en/stable/index.html>`_ that incorporates renewable energy models for RE-driven desalination systems.

WaterTAP is a National Alliance for Water Innovation (NAWI) funded initiative to create an open-source water treatment model library
that is compatible with the IDAES Platform, an advanced process systems engineering tool developed by the U.S. Department of Energy.
WaterTAP is a dependency of WaterTAP-REFLO, and users should refer to the WaterTAP documentation for information on installation,
modeling framework, and general unit models. All WaterTAP models and features are available in WaterTAP-REFLO.

Importantly, because REFLO is an extension of WaterTAP, REFLO imports are prefixed with ``watertap_contrib.reflo`` instead of e.g., ``reflo``.

.. code-block:: python

    import watertap_contrib.reflo as reflo


Collaborating Institutions
--------------------------

The WaterTAP-REFLO development team is composed of researchers from:

* National Renewable Energy Laboratory
* National Energy Technology Laboratory

Content
-------

.. toctree::
   :maxdepth: 2
   
   getting_started
   technical_reference/index
   how_to_guides/index
..
    license
