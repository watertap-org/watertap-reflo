Low Temperature - Multi-effect Distillation (LT-MED)
====================
This Low Temperature Multi-effect Distillation (LT-MED) unit model
   * supports steady-state only
   * is a surrogate model
   * is verified against the operation data in Plataforma Solar de Almeria (PSA)

#TODO: Add index/reference to home page


Ports
---------

The model provides four ports (Pyomo notation in parenthesis):

* Feed flow port (feed)
* Distillate port (distillate)
* Brine flow port (brine)
* Heating steam port (steam)


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





References
----------

[1] Palenzuela, P., Hassan, A. S., Zaragoza, G., & Alarcón-Padilla, D. C. (2014). Steady state model for
multi-effect distillation case study: Plataforma Solar de Almería MED pilot plant. Desalination, 337,
31-42.

[2] Ortega-Delgado, B., Garcia-Rodriguez, L., & Alarcón-Padilla, D. C. (2017). Opportunities of
improvement of the MED seawater desalination process by pretreatments allowing high-temperature
operation. Desalin Water Treat, 97, 94-108.