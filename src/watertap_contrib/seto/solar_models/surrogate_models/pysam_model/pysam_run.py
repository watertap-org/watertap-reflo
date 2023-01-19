import json
from os.path import join, dirname
from math import floor, ceil, isnan
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import PySAM.TroughPhysicalProcessHeat as iph
import PySAM.IphToLcoefcr as iph_to_lcoefcr
import PySAM.Lcoefcr as lcoefcr

def load_config_files(file_names, modules):
    for file_name, module in zip(file_names, modules):
        with open(file_name, 'r') as file:
            data = json.load(file)
            missing_values = []             # for debugging
            for k, v in data.items():
                if k != "number_inputs":
                    try:
                        module.value(k, v)
                    except:
                        missing_values.append(k)
            pass

def tes_cost(tech_model):
    STORAGE_COST_SPECIFIC = 62                                  # [$/kWht] borrowed from physical power trough
    tes_thermal_capacity = tech_model.value('q_pb_design') * 1e3 \
                         * tech_model.value('tshours')          # [kWht]
    return tes_thermal_capacity * STORAGE_COST_SPECIFIC

def system_capacity(tech_model):
    return tech_model.value('q_pb_design') * tech_model.value('specified_solar_multiple') * 1e3  # [kW]

def setup_model(model_name, config_files, weather_file):
    tech_model = iph.new()
    post_model = iph_to_lcoefcr.from_existing(tech_model, model_name)
    cash_model = lcoefcr.from_existing(tech_model, model_name)
    modules = [tech_model, post_model, cash_model]
    load_config_files(config_files, modules)
    tech_model.Weather.file_name = weather_file

    # Determine storage cost component
    capital_cost_orig = cash_model.value('capital_cost')
    storage_cost_orig = tes_cost(tech_model)
    capital_cost_minus_storage = capital_cost_orig - storage_cost_orig
    capital_cost_minus_storage_per_kW = capital_cost_minus_storage / system_capacity(tech_model)

    fixed_operating_cost_per_kW = cash_model.value('fixed_operating_cost') / system_capacity(tech_model)

    return {
        'tech_model': tech_model,
        'post_model': post_model,
        'cash_model': cash_model,
        'capital_cost_minus_storage_per_kW': capital_cost_minus_storage_per_kW,
        'fixed_operating_cost_per_kW': fixed_operating_cost_per_kW
    }

def run_model(modules, heat_load=None, hours_storage=None):
    tech_model = modules['tech_model']
    post_model = modules['post_model']
    cash_model = modules['cash_model']

    if heat_load is not None:
        tech_model.value('q_pb_design', heat_load)
    if hours_storage is not None:
        tech_model.value('tshours', hours_storage)
    tech_model.execute()

    # NOTE: freeze_protection_field can sometimes be nan (when it should be 0) and this causes other nan's
    #  Thus, freeze_protection, annual_energy and capacity_factor must be calculated manually
    # annual_energy = tech_model.Outputs.annual_energy                            # [kWht] net, does not include that used for freeze protection
    # freeze_protection = tech_model.Outputs.annual_thermal_consumption           # [kWht]
    # capacity_factor = tech_model.Outputs.capacity_factor                        # [%]
    freeze_protection_field = tech_model.Outputs.annual_field_freeze_protection
    freeze_protection_field = 0 if isnan(freeze_protection_field) else freeze_protection_field      # occasionally seen to be nan
    freeze_protection_tes = tech_model.Outputs.annual_tes_freeze_protection
    freeze_protection_tes = 0 if isnan(freeze_protection_tes) else freeze_protection_tes
    freeze_protection = freeze_protection_field + freeze_protection_tes
    annual_energy = tech_model.Outputs.annual_gross_energy - freeze_protection  # [kWht] net, does not include that used for freeze protection
    capacity_factor = annual_energy / (tech_model.value('q_pb_design') * 1e3 * 8760) * 100 # [%] 
    electrical_load = tech_model.Outputs.annual_electricity_consumption         # [kWhe]
    
    post_model.value('fixed_operating_cost', modules['fixed_operating_cost_per_kW'] * system_capacity(tech_model))
    post_model.execute()

    cash_model.value('annual_energy', annual_energy)                            # override the linked value from the tech model, which could be nan
    cash_model.value('capital_cost', modules['capital_cost_minus_storage_per_kW'] * system_capacity(tech_model)
									 + tes_cost(tech_model))
    cash_model.execute()

    lcoh = cash_model.Outputs.lcoe_fcr                                          # [$/kWht]
    capital_cost = cash_model.value('capital_cost')                             # [$]
    fixed_operating_cost = cash_model.value('fixed_operating_cost')             # [$] more than shown in UI, includes electricity purchases
    variable_operating_cost = cash_model.value('variable_operating_cost')       # [$]

    return {
        'annual_energy': annual_energy,                                         # [kWh] annual net thermal energy in year 1
        'freeze_protection': freeze_protection,                                 # [kWht]
        'capacity_factor': capacity_factor,                                     # [%] capacity factor
        'electrical_load': electrical_load,                                     # [kWhe]
        'lcoh': lcoh,                                                           # [$/kWht] LCOH
        'capital_cost': capital_cost,                                           # [$]
        'fixed_operating_cost': fixed_operating_cost,                           # [$]
        'variable_operating_cost': variable_operating_cost,                     # [$]
        }

def plot_3d(df, x_index=0, y_index=1, z_index=2, grid=True, countour_lines=True):
    """
    index 0 = x axis
    index 1 = y axis
    index 2 = z axis
    """
    # 3D PLOT
    # fig = plt.figure(figsize=(8,6))
    # ax = fig.add_subplot(1, 1, 1, projection='3d')
    # surf = ax.plot_trisurf(df.iloc[:,0], df.iloc[:,1], df.iloc[:,2], cmap=plt.cm.viridis, linewidth=0.2)
    # modld_pts = ax.scatter(df.iloc[:,0], df.iloc[:,1], df.iloc[:,2], c='black', s=15)
    # ax.set_xlabel(df.columns[0])
    # ax.set_ylabel(df.columns[1])
    # ax.set_zlabel(df.columns[2])
    # plt.show()

    def _set_aspect(ax, aspect):
        x_left, x_right = ax.get_xlim()
        y_low, y_high = ax.get_ylim()
        ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*aspect)

    # CONTOUR PLOT
    levels = 25
    df2 = df.pivot(
        df.columns[y_index],
        df.columns[x_index],
        df.columns[z_index])
    y = df2.index.values
    x = df2.columns.values
    z = df2.values
    fig, ax = plt.subplots(1,1)
    cs = ax.contourf(x, y, z, levels=levels)
    if countour_lines:
        cl = ax.contour(x, y, z, colors='black', levels=levels)
        ax.clabel(cl, colors='black', fmt='%#.4g')
    if grid:
        ax.grid(color='black')
    _set_aspect(ax, 0.5)
    fig.colorbar(cs)
    ax.set_xlabel(df.columns[x_index])
    ax.set_ylabel(df.columns[y_index])
    ax.set_title(df.columns[z_index])
    plt.show()


#########################################################################################################
if __name__ == '__main__':
    model_name =        'PhysicalTroughIPHLCOHCalculator'
    config_files =      [join(dirname(__file__), 'untitled_trough_physical_process_heat.json'),
                        join(dirname(__file__), 'untitled_iph_to_lcoefcr.json'),
                        join(dirname(__file__), 'untitled_lcoefcr.json'),]
    weather_file =      join(dirname(__file__), 'tucson_az_32.116521_-110.933042_psmv3_60_tmy.csv')

    modules = setup_model(model_name, config_files, weather_file)

    # Run default model
    # result = run_model(modules, heat_load=None, hours_storage=None)

    # Run model at specific parameters
    # result = run_model(modules, heat_load=600, hours_storage=3)

    # Load and plot saved df
    # df = pd.read_pickle('pickle_filename2.pkl')
    # plot_3d(df, 0, 1, 2, grid=False, countour_lines=False)      # annual energy
    # plot_3d(df, 0, 1, 3, grid=False, countour_lines=False)      # capacity factor
    # plot_3d(df, 0, 1, 4, grid=False, countour_lines=False)      # LCOH

    # Run parametrics
    data = []
    # heat_loads =        np.arange(5, 115, 10)       # [MWt]
    heat_loads =        np.arange(100, 1100, 100)   # [MWt]
    hours_storages =    np.arange(0, 27, 3)         # [hr]
    comb = [(hl, hs) for hl in heat_loads for hs in hours_storages]
    for heat_load, hours_storage in comb:
        result = run_model(modules, heat_load, hours_storage)
        total_cost = result['capital_cost'] + result['fixed_operating_cost'] + result['variable_operating_cost']
        data.append([
            heat_load,
            hours_storage,
            result['annual_energy'],
            result['capacity_factor'],
            result['lcoh'],
            total_cost,
        ])
    df = pd.DataFrame(data, columns=[
        'heat_load',
        'hours_storage',
        'annual_energy',
        'capacity_factor',
        'lcoh',
        'total_cost'])

    df.to_pickle('pickle_filename4.pkl')
    plot_3d(df, 0, 1, 2, grid=False, countour_lines=False)    # annual energy
    plot_3d(df, 0, 1, 3, grid=False, countour_lines=False)    # capacity factor
    plot_3d(df, 0, 1, 4, grid=False, countour_lines=False)    # lcoh

    x = 1   # for breakpoint
    pass
