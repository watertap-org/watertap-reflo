Low Temperature - Multi-effect Distillation (LT-MED)
====================
This Low Temperature Multi-effect Distillation (LT-MED) unit model
   * supports steady-state only
   * is a surrogate model
   * is verified against the operation data in Plataforma Solar de Almeria (PSA)

#TODO: Add index/reference to home page


Degrees of Freedom
------------------
The LT-MED model has at least 5 degrees of freedom that should be fixed for the unit to be fully specified.

Typically, the following variables are fixed, including the state variables at the inlet:
    * Feed salinity (30-60 ":math:`\text{g/}\text{L}`" = ":math:`\text{kg/}\text{m}^3`")
    * Feed temperature (15-35 ":math:`\text{C}`")
    * Heating steam temperature (60-85 ":math:`\text{C}`")
    * Recovery ratio (30%- 50%)
    * System capacity (> 2000 ":math:`\text{m}^3\text{/day}`")

The first four variables are independent input variables to the surrogate equations.


Model Structure
---------------
This LT-MED model consists of 4 StateBlocks (as 4 Ports in parenthesis below).

* Feed flow (feed)
* Distillate (distillate)
* Brine flow (brine)
* Heating steam (steam)


Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq', 'Vap']"
   "Components", ":math:`j`", "['H2O', 'TDS']"


Variables
---------
The system configuration variables should be fixed at the default values, with which the surrogate model was developed:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Value", "Units"

   "Temperature difference between the last and first effect", ":math:`\Delta\T_{last}`", "delta_T_last_effect", "10", ":math:`\text{K}`"
   "Temperature decrease in cooling reject water", ":math:`\Delta\T_{cooling}`", "delta_T_cooling_reject", "-3", ":math:`\text{K}`"


The following performance variables are derived from the surrogate equations:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Gain output ratio", ":math:`GOR`", "gain_output_ratio", "None", ":math:`\text{dimensionless}`"
   "Specific total area", ":math:`sA`", "specific_area_per_m3_day", "None", ":math:`\text{m}^2\text{/m}^3\text{/day}`"


Equations
---------
.. csv-table::
   :header: "Description", "Equation"

   "Surrogate equation foor calculating GOR", ":math:`GOR = a_{1}X_{f}+a_{2}RR+a_{3}RRX_{f}+a_{4}T_{N}+a_{5}T_{N}X_{f}+a_{6}T_{N}RR+a_{7}T_{s}X_{f}+a_{8}T_{s}X_{f}+a_{9}T_{s}RR+a_{10}T_{s}T_{N}+a_{11}+a_{12}{T_s}^2+a_{13}{T_N}^2+a_{14}{RR^2+a_{15}{X_f}^2`"


References
----------

[1] Palenzuela, P., Hassan, A. S., Zaragoza, G., & Alarcón-Padilla, D. C. (2014). Steady state model for
multi-effect distillation case study: Plataforma Solar de Almería MED pilot plant. Desalination, 337,
31-42.

[2] Ortega-Delgado, B., Garcia-Rodriguez, L., & Alarcón-Padilla, D. C. (2017). Opportunities of
improvement of the MED seawater desalination process by pretreatments allowing high-temperature
operation. Desalin Water Treat, 97, 94-108.