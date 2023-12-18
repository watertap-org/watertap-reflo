How to create a surrogate model
-------------------------------

Surrogate modeling toolboxes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* `PySMO <https://idaes-pse.readthedocs.io/en/1.5.1/surrogate/pysmo/index.html>`_
  
  Stuff and Text
* `ALAMO <https://idaes-pse.readthedocs.io/en/1.5.1/surrogate/alamopy.html>`_
  
  Stuff and Text


Example 1: Surrogate model for PV energy generation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
PySAM

Radial Basis Function

.. csv-table::
   :header: "Variables", "Variable name", "Symbol", "Valid range", "Unit"

      "Design Size", "feed_props.conc_mass_phase_comp['Liq', 'TDS']", ":math:`X_{f}`", "35 - 292", ":math:`\text{g/}\text{L}`"
      "Electricity", "feed_props.conc_mass_phase_comp['Liq', 'TDS']", ":math:`X_{f}`", "35 - 292", ":math:`\text{g/}\text{L}`"