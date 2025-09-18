.. _solar_still_ref:

Solar Still
===========


.. code-block:: python

    from watertap_contrib.reflo.unit_models import SolarStill

Solar stills use solar energy to evaporate water, leaving contaminants behind.
They are often used in remote areas where access to electricity is limited but solar energy and space is abundant.


Model Structure
---------------

The solar still calculates the water yield based on local atmospheric weather conditions. For each hour in the day,
it calculates the water yield. Based on the user-defined initial depth of water in the solar still system, it determines how long it takes
for the water level to drop below a certain threshold, indicating the end of the solar still batch operation.
The duration of the batch operation is then used to calculate the total water yield for the system, and calculates the
total area required to meet a user-defined daily water production target.

The `MCAS property package <https://watertap.readthedocs.io/en/stable/technical_reference/property_models/mc_aq_sol.html>`_ is used for the solar still model.
The model consists of 3 StateBlocks (with 3 Ports in parenthesis below).

* Feed flow (``inlet``)
* Product water (``outlet``)
* Waste flow (``waste``)


.. Degrees of Freedom
.. ------------------

The only degrees of freedom that must be specified for the solar still model are the inlet state variables (i.e. temperature, pressure, component flowrates).
Users must also provide a weather data file that contains hourly weather data for the location of interest.
The weather data file can be downloaded from the `National Solar Radiation Database <https://nsrdb.nrel.gov/data-viewer>`_.
The weather data file must be a CSV file that contains the following columns:

.. csv-table::
   :header: "Column Name", "Description", "Units"

   "GHI", "Global horizontal irradiance", ":math:`\text{W/m}^2`"
   "Ambient Temperature", "Ambient air temperature", ":math:`\text{°C}`"
   "Wind Speed", "Wind speed at 10 m height", ":math:`\text{m/s}`"

The hourly weather data will be converted to arrays for every second of a year, with the hourly 
values repeated for each second of the hour.

The water yield calculation will calculate the following variables for each second of the year:

.. csv-table::
   :header: "Description", "Symbol", "Units"

   "Temperature of the sky", ":math:`T_{sky}`", ":math:`\text{K}`"
   "Temperature of the basin", ":math:`T_{basin}`", ":math:`\text{K}`"
   "Temperature of the glass", ":math:`T_{glass}`", ":math:`\text{K}`"
   "Temperature of saline water", ":math:`T_{sw}`", ":math:`\text{K}`"
   "Temperature difference inside the basin", ":math:`\Delta T_{inside}`", ":math:`\text{K}`"
   "Temperature difference outside the basin", ":math:`\Delta T_{outside}`", ":math:`\text{K}`"
   "Depth of the water in the basin", ":math:`A`", ":math:`\text{m}`"
   "Mass of water in the basin", ":math:`M_{sw}`", ":math:`\text{kg}`"
   "Mass of salt in the basin", ":math:`M_{salt}`", ":math:`\text{kg}`"
   "Density of saline water", ":math:`\rho`", ":math:`\text{kg/m}^3`"
   "Dynamic viscosity of saline water", ":math:`\mu`", ":math:`\text{kg/m/s}`"
   "Kinematic viscosity of saline water", ":math:`\nu`", ":math:`\text{m}^2/\text{s}`"
   "Specific heat of saline water", ":math:`c_{p}`", ":math:`\text{kJ/kg/K}`"
   "Thermal conductivity of saline water", ":math:`k`", ":math:`\text{W/m/K}`"
   "Latent heat of vaporization of water", ":math:`h_{fg}`", ":math:`\text{kJ/kg}`"
   "Activity of saltwater", ":math:`a_{sw}`", ":math:`\text{dimensionless}`"
   "Partial saturated vapor pressure of water at saline water temperature", ":math:`P_{sw}`", ":math:`\text{Pa}`"
   "Partial saturated vapor pressure of water at glass temperature", ":math:`P_{a}`", ":math:`\text{Pa}`"
   "Coefficient of volume expansion of water", ":math:`\beta`", ":math:`\text{K}^{-1}`"
   "Prandtl number", ":math:`\text{Pr}`", ":math:`\text{dimensionless}`"
   "Grashof number", ":math:`\text{Gr}`", ":math:`\text{dimensionless}`"
   "Heat transfer coefficient of water layer", ":math:`h_{sw}`", ":math:`\text{W/m}^2/\text{K}`"
   "Convective heat transfer coefficient between water and glass", ":math:`h_{c,water-glass}`", ":math:`\text{W/m}^2/\text{K}`"
   "Radiative heat transfer coefficient between water and glass", ":math:`h_{r,water-glass}`", ":math:`\text{W/m}^2/\text{K}`"
   "Evaporative heat transfer coefficient between water and glass", ":math:`h_{e,water-glass}`", ":math:`\text{W/m}^2/\text{K}`"
   "Total heat transfer coefficient between water and glass", ":math:`h_{water-glass}`", ":math:`\text{W/m}^2/\text{K}`"
   "Radiative heat transfer coefficient between glass and ambient", ":math:`h_{r,glass-amb}`", ":math:`\text{W/m}^2/\text{K}`"
   "Convective heat transfer coefficient between glass and ambient", ":math:`h_{c,glass-amb}`", ":math:`\text{W/m}^2/\text{K}`"
   "Convective heat transfer coefficient between basin and ambient", ":math:`h_{c,basin-amb}`", ":math:`\text{W/m}^2/\text{K}`"
   "Total heat transfer coefficient between glass and ambient", ":math:`h_{glass-amb}`", ":math:`\text{W/m}^2/\text{K}`"
   "Total heat transfer coefficient between basin and ambient", ":math:`h_{basin-amb}`", ":math:`\text{W/m}^2/\text{K}`"
   "Effective overall absorptivity of the basin", ":math:`\alpha_{eff}`", ":math:`\text{dimensionless}`"
   "Overall heat loss coefficient between glass and ambient", ":math:`U_{glass-amb}`", ":math:`\text{W/m}^2/\text{K}`"
   "Overall heat transfer coefficient between basin bottom and ambient", ":math:`U_{basin-amb}`", ":math:`\text{W/m}^2/\text{K}`"
   "Overall heat loss coefficient between basin bottom and ambient", ":math:`U_{overall}`", ":math:`\text{W/m}^2/\text{K}`"
   "Overall heat loss coefficient from basin sides", ":math:`U_{sides}`", ":math:`\text{W/m}^2/\text{K}`"
   "Overall heat transfer coefficient from basin to ambient", ":math:`U_{total}`", ":math:`\text{W/m}^2/\text{K}`"
   "Evaporated mass of freshwater per unit area", ":math:`m_{evap}`", ":math:`\text{kg/m}^2/\text{s}`"
   "Remaining mass of saline water in the basin", ":math:`m_{sw}`", ":math:`\text{kg}`"
   "Salt concentration in the basin", ":math:`C_{salt}`", ":math:`\text{kg/m}^3`"
   "Mass of salt precipitated", ":math:`m_{salt,precip}`", ":math:`\text{kg}`"

..    "Radiative heat transfer coefficient", ":math:`h_{r}`", ":math:`\text{W/m}^2/\text{K}`"
..    "Overall heat transfer coefficient", ":math:`U`", ":math:`\text{W/m}^2/\text{K}`"
..    "Mass transfer coefficient", ":math:`k_{m}`", ":math:`\text{kg/m}^2/\text{s}/\text{Pa}`"
..    "Water evaporation rate per unit area", ":math:`m_{evap}`", ":math:`\text{kg/m}^2/\text{s}`"
..    "Water yield per unit area", ":math:`Y_{evap}`", ":math:`\text{m}^3/\text{m}^2/\text{s}`"

The model uses the following constants in the water yield calculations:


.. csv-table::
    :header: "Description", "Symbol", "Value", "Units"

    "Gravitational constant", ":math:`g`", ":math:`9.81`", ":math:`\text{m/s}^2`"
    "Stefan-Boltzmann constant", ":math:`\sigma`", ":math:`5.6697 \times 10^{-8}`", ":math:`\text{W/m}^2\text{K}^4`"
    "Thickness of solar still insulation", ":math:`x_{insul}`", ":math:`0.005`", ":math:`\text{m}`"
    "Thermal conductivity of solar still insulation", ":math:`k_{insul}`", ":math:`0.033`", ":math:`\text{W/m/K}`"
    "Thickness of glass cover", ":math:`x_{glass}`", ":math:`0.004`", ":math:`\text{m}`"
    "Thermal conductivity of glass cover", ":math:`k_{glass}`", ":math:`1.03`", ":math:`\text{W/m/K}`"
    "Density of NaCl solid", ":math:`\rho_{salt}`", ":math:`2165`", ":math:`\text{g/L}`"
    "Maximum solubility of NaCl in water", ":math:`C_{salt,max}`", ":math:`365`", ":math:`\text{g/L}`"
    "Adsorptivity of glass cover", ":math:`\alpha_{glass}`", ":math:`0.047`", ":math:`\text{dimensionless}`"
    "Adsorptivity of water surface", ":math:`\alpha_{water}`", ":math:`0.20`", ":math:`\text{dimensionless}`"
    "Adsorptivity of basin", ":math:`\alpha_{basin}`", ":math:`0.65`", ":math:`\text{dimensionless}`"
    "Reflectivity of glass cover", ":math:`R_{glass}`", ":math:`0.047`", ":math:`\text{dimensionless}`"
    "Reflectivity of water surface", ":math:`R_{water}`", ":math:`0.08`", ":math:`\text{dimensionless}`"
    "Emissivity of glass cover", ":math:`\epsilon_{glass}`", ":math:`0.94`", ":math:`\text{dimensionless}`"
    "Emissivity of water surface", ":math:`\epsilon_{water}`", ":math:`0.95`", ":math:`\text{dimensionless}`"

These are used to calculate other constants used in the model:


.. csv-table::
    :header: "Description", "Symbol", "Units", "Equation"

    "Effective emissivity between glass and water surface", ":math:`\epsilon_{effective}`", ":math:`\text{dimensionless}`", ":math:`\cfrac{1}{\left( \cfrac{1}{\epsilon_{glass}} + \cfrac{1}{\epsilon_{water}} - 1 \right)}`"
    "Effective absorptivity of solar radiation absorbed by water", ":math:`\alpha_{water,eff}`", ":math:`\text{dimensionless}`", ":math:`\alpha_{water} \left(1 - \alpha_{glass}\right) \left(1 - R_{glass}\right) \left(1 - R_{water} \right)`"
    "Effective absorptivity of solar radiation absorbed by basin", ":math:`\alpha_{basin,eff}`", ":math:`\text{dimensionless}`", ":math:`\alpha_{basin} \left(1 - \alpha_{glass}\right) \left(1 - R_{glass} \right) \left(1 - \alpha_{water} \right) \left(1 - R_{water} \right)`"
    "Effective absorptivity of solar radiation absorbed by glass", ":math:`\alpha_{glass,eff}`", ":math:`\text{dimensionless}`", ":math:`\alpha_{glass} \left(1 - R_{glass} \right)`"

Finally, the default values for the dimension of the basin are:


.. csv-table::
    :header: "Description", "Symbol", "Units", "Equation"

    "Length of the basin", ":math:`L`", ":math:`\text{m}`", "0.6"
    "Width of the basin", ":math:`W`", ":math:`\text{m}`", "0.6"
    "Initial depth of water in the basin", ":math:`Z`", ":math:`\text{m}`", "0.1"


Equations
---------


.. csv-table::
    :header: "Description", "Equation"

    "Temperature difference inside the basin", ":math:`\Delta T_{inside,t} = T_{sw,t-1} - T_{glass,t-1}`"
    "Temperature difference outside the basin", ":math:`\Delta T_{outside,t} = T_{glass,t-1} - T_{ambient,t-1}`" 
    "Temperature of the sky", ":math:`T_{sky,t} = 0.0552 \left( T_{ambient}^{1.5} \right)`"
    "Area of the water on the basin side", ":math:`A_{side, t} = 2 \times (L + W) \times Z_{t-1}`"
    "Saltwater density", ":math:`\rho_{sw,t} = f\left(C_{salt,t-1}, T_{sw,t-1}\right)`" 
    "Saltwater dynamic viscosity", ":math:`\mu_{sw,t} = f\left(C_{salt,t-1}, T_{sw,t-1}\right)`" 
    "Saltwater specific heat", ":math:`c_{p,sw,t} = f\left(C_{salt,t-1}, T_{sw,t-1}\right)`"
    "Saltwater thermal conductivity", ":math:`k_{sw,t} = f\left(C_{salt,t-1}, T_{sw,t-1}\right)`"
    "Saltwater kinematic viscosity", ":math:`\nu_t = \cfrac{\mu_{sw,t}}{\rho_{sw,t}}`"
    "Prandtl number", ":math:`\text{Pr}_{t} = \cfrac{c_{p,sw,t} \mu_{sw,t}}{k_{sw,t}}`"
    "Latent heat of vaporization of pure water", ":math:`h_{fg} = \left( 2501.67 - 2.389 T_{sw,t-1}\right) \times 1000`"
    "Water activity", ":math:`a_{sw,t} = -0.000537 C_{salt,t-1} + 0.9985307`"
    "Partial saturated vapor pressure of water at saline water temperature", ":math:`P_{sw,t} = a_{sw,t} \times \text{exp}{\left(  25.317 - \cfrac{5144}{(T_{sw,t-1} + 273)} \right)}`"
    "Partial saturated vapor pressure of water at glass temperature", ":math:`P_{glass,t} = a_{sw,t} \times \text{exp}{\left(  25.317 - \cfrac{5144}{(T_{glass,t-1} + 273)} \right)}`"
    "Coefficient of volume expansion of water", ":math:`\beta_{t} = -0.000006 \times T_{sw,t-1}^4 + 0.001667 \times T_{sw,t-1}^3 - 0.197796 \times T_{sw,t-1}^2 + 16.862446 \times T_{sw,t-1} - 64.319951`"
    "Grashof number", ":math:`\text{Gr}_{t} = \cfrac{g \beta_{t} (T_{basin,t-1} - T_{sw,t-1}) Z_{t-1}^3}{\nu_t^2}`"
    "Heat transfer coefficient of water layer", ":math:`h_{sw,t} = \cfrac{k}{Z_{t-1}} \left(0.54 \text{Pr}_{t} \text{Gr}_{t}\right)^{0.25}`"
    "Convective heat transfer coefficient between water and glass", ":math:`h_{c,water-glass,t} = 0.884 \left( (T_{sw,t-1} - T_{glass,t-1}) + \left( \cfrac{(P_{sw,t} - P_{glass,t})(T_{sw,t-1} + 273.15)}{268900 - P_{sw,t}} \right)  \right)^{1/3}`"
    "Radiative heat transfer coefficient between water and glass", ":math:`h_{r,water-glass,t} = \epsilon_{water} \sigma \left( (T_{sw,t-1} + 273.15)^2 + (T_{glass,t-1} + 273.15)^2 \right) \left( T_{sw,t-1} + T_{glass,t-1} + 546 \right)`"
    "Evaporative heat transfer coefficient between water and glass", ":math:`h_{e,water-glass,t} = 16.273 \times h_{c,water-glass,t} \times \cfrac{P_{sw,t} - P_{glass,t}}{\Delta T_{inside, t}}`"
    "Total heat transfer coefficient between water and glass", ":math:`h_{water-glass,t} = h_{c,water-glass,t} + h_{r,water-glass,t} + h_{e,water-glass,t}`"
    "Radiative heat transfer coefficient between glass and ambient", ":math:`h_{r,glass-amb,t} = \epsilon_{glass} \sigma \left( \cfrac{(T_{glass,t-1} + 273.15)^4 - (T_{sky,t-1} + 273.15)^4}{\Delta T_{inside, t}} \right)`"
    "Convective heat transfer coefficient between glass and ambient", ":math:`h_{c,glass-amb,t} = 2.8 + \left( 3.0 \times V_{wind,t}\right)`"
    "Convective heat transfer coefficient between basin and ambient", ":math:`h_{c,basin-amb,t} = 2.8 + \left( 3.0 \times V_{wind,t}\right)`"
    "Total heat transfer coefficient between glass and ambient", ":math:`h_{glass-amb,t} = h_{r,glass-amb,t} + h_{c,glass-amb,t}`"
    "Total heat transfer coefficient between basin and ambient", ":math:`h_{basin-amb,t} = \cfrac{k_{insul}}{x_{insul}} + \cfrac{1}{h_{c,basin-amb,t}}`"
    "Effective overall absorptivity of the basin", ":math:`\alpha_{eff} = \alpha_{basin,eff} \cfrac{h_{sw,t}}{h_{sw,t} + h_{basin-amb,t} + h_{c,basin-amb,t}}+ \alpha_{water,eff} + \alpha_{glass,eff} \left( \cfrac{h_{water-glass,t}}{h_{water-glass,t} + h_{glass-amb,t}} \right)`"
    "Overall heat loss coefficient between glass and ambient", ":math:`U_{glass-amb,t} = \cfrac{ \cfrac{k_{glass}}{x_{glass}} h_{glass-amb,t} } { \cfrac{k_{glass}}{x_{glass}}+h_{glass-amb,t}}`"
    "Overall heat transfer coefficient between basin bottom and ambient", ":math:`U_{basin-amb,t} = \cfrac{ h_{glass-amb,t} U_{glass-amb,t}} { h_{glass-amb,t} + U_{glass-amb,t}}`"
    "Overall heat loss coefficient between basin bottom and ambient", ":math:`U_{overall,t} = \cfrac{ h_{sw,t} h_{basin-amb,t}} { h_{sw,t} + h_{basin-amb,t}}`"
    "Overall heat loss coefficient from basin sides", ":math:`U_{sides,t} = \cfrac{A_{side}}{A_{bottom}} U_{overall,t}`"
    "Overall heat transfer coefficient from basin to ambient", ":math:`U_{total,t} = U_{overall,t} + U_{sides,t}`"
    "Overall external heat transfer coefficient", ":math:`U_{external,t} = U_{basin-amb,t} + U_{total,t}`"
    "Grouping term for energy balance", ":math:`Q_{group,t} = \cfrac{U_{external,t}}{m_{sw,t-1} c_{p,sw,t}}`"
    "Time dependent term for energy balance", ":math:`X_t = \cfrac{\left( \alpha_{eff} * \text{GHI} \right) + \left(U_{external,t} T_{ambient,t} \right)}{m_{sw,t-1} c_{p,sw,t}}`"
    "Saline water temperature", ":math:`T_{sw,t} = \left( \cfrac{X_t}{Q_{group, t}} \right) \left( 1 - \text{exp} \left(-Q_{group,t} t \right) \right) + \left( T_{sw,t-1} \text{exp} \left(-Q_{group,t} t \right) \right)`"
    "Glass temperature", ":math:`T_{glass,t} = \cfrac{ \alpha_{glass,eff} \text{GHI} + h_{water-glass,t} T_{sw,t-1} +  U_{glass-amb,t} T_{ambient,t}}{h_{water-glass,t} + U_{glass-amb,t}}`"
    "Basin temperature", ":math:`T_{basin,t} = \cfrac{ \alpha_{basin,eff} \text{GHI} + h_{sw,t} T_{sw,t-1} + (h_{basin-amb,t} + h_{c,basin-amb,t}) T_{basin,t-1}}{h_{sw,t} + U_{basin-amb,t} + h_{c,basin-amb,t}}`"
    "Evaporated mass of freshwater per unit area", ":math:`m_{fw,evap,t} = \cfrac{A_{basin} h_{e,water-glass,t} \Delta T_{inside}}{h_{fg}}`"
    "Evaporated mass of saltwater", ":math:`m_{sw,evap,t} = \cfrac{m_{fw,evap,t}}{1 + \cfrac{C_{salt,t-1}}{\rho_{fw}}}`"
    "Remaining mass of freshwater in the basin", ":math:`m_{fw,t} = m_{fw,t-1} - m_{sw,evap,t}`"
    "Remaining mass of saltwater in the basin", ":math:`m_{sw,t} = m_{fw,t} + M_{salt}`"
    "Depth of water in the basin", ":math:`Z_t = \cfrac{m_{sw,t}}{\rho_{sw,t} A_{basin}}`"
    "Salt concentration in the basin :sup:`1`", ":math:`C_{salt,t} = \cfrac{M_{salt} \rho_{fw}}{m_{fw,t}}`"
    "Excess salinity that precipitates", ":math:`m_{salt,precip,t} = \text{max} \left( 0, (C_{salt,t} - C_{salt,max}) \times \cfrac{m_{fw,t}}{\rho_{fw}} \right)`"





:math:`k_{glass}` = conductivity of glass

.. note:: 
   :sup:`1` This is the break point in the calculation where the model checks if the water depth is \<= 0 or if all the water has evaporated.
   
   The salt concentration in the basin is capped at a maximum value (``max_salt_concentration``) to prevent numerical issues when all the water has evaporated.
      This means that the model assumes that any salt that would cause the concentration to exceed this value precipitates out of solution.
      The default value is 300,000 mg/L (i.e., 300 g/L), which is close to the solubility limit of NaCl in water at room temperature.
      Users can adjust this value as needed.
:sup:`1` BREAKPOINT 1

For each second in the year, the model will check if either the depth of water in the basin is \<= 0 of if all the water has been evaporated.
If either of these conditions are met, the model will stop the water yield calculation and the length of the batch is the number of seconds that have passed (i.e., :math:`t_{batch} = t`).
This is used to calculate the number of batch cycles per year:

.. math::

   n_{cycles} = \cfrac{31536000}{t_{batch}}

The total water yield per year is then calculated as:

.. math::

    Y = \cfrac{m_{fw,t0} n_{cycles}}{A_{basin}}

The model assumes that all of the water that evaporates is collected as product water:

.. math::

    m_{in,water} = m_{out,water}

And that all other components remain in the waste stream:

.. math::

    m_{in,component} = m_{waste,component}

This implicitly means that the salt concentration in the product water is 0 g/L and the mass of water in the waste stream is zero.
Then, the model calculates the required area of the solar still to meet a user-defined daily water production target:

.. math::

   A_{basin,tot} = \cfrac{m_{out,water}}{Y}



Costing
----------

References
----------
