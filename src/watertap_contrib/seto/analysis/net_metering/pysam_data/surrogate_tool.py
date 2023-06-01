import json
import os
from os.path import join, dirname
import sys
from io import StringIO
import re
import datetime
from math import floor, ceil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing
import time
from watertap_contrib.seto.core import SETODatabase, PySAMWaterTAP
from idaes.core.surrogate.sampling.data_utils import split_training_validation, split_training_validation_testing
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate

absolute_path = os.path.dirname(__file__)
print(absolute_path)

tech_config_file = "/pvsamv1.json"
tech_config_file = absolute_path + tech_config_file
grid_config_file = "/grid.json"
grid_config_file = absolute_path + grid_config_file
rate_config_file = "/utilityrate5.json"
rate_config_file = absolute_path + rate_config_file
cash_config_file = "/singleowner.json"
cash_config_file = absolute_path + cash_config_file
weather_file = "/phoenix_az_33.450495_-111.983688_psmv3_60_tmy.csv"
weather_file = absolute_path + weather_file

dataset_filename = os.path.join(
        os.path.dirname(__file__), "dataset_2.pkl"
    )  # output dataset for surrogate training
daily_dataset_filename = os.path.join(
        os.path.dirname(__file__), "dataset_daily.pkl"
    )  # output dataset for surrogate training

pysam = PySAMWaterTAP(
        pysam_model="pv",
        tech_config_file=tech_config_file,
        grid_config_file=grid_config_file,
        rate_config_file=rate_config_file,
        cash_config_file=cash_config_file,
        weather_file=weather_file,
    )

output_labels = ['period_energy', 'land_req']

def run(desired_size):  
    pysam.run_pv_single_owner(desired_size=desired_size)
    return pd.DataFrame({'design_size':desired_size, 'period_energy':list(pysam.tech_model.Outputs.monthly_energy), "land_req": pysam.land_req, 'month':np.linspace(1,12,12)})

def get_hourly(desired_size):  
    pysam.run_pv_single_owner(desired_size=desired_size)
    return pd.DataFrame({'design_size':desired_size, 'hourly_gen':pysam.tech_model.Outputs.gen})

def get_rep_days(result_df):
    year_start = datetime.datetime(year=2020, month=1, day=1, hour=0, minute=0)
    sum_solstice = datetime.datetime(year=2020, month=6, day=20, hour=0, minute=0)
    win_solstice = datetime.datetime(year=2020, month=12, day=21, hour=0, minute=0)
    ver_eq = datetime.datetime(year=2020, month=3, day=20, hour=0, minute=0)
    aut_eq = datetime.datetime(year=2020, month=9, day=22, hour=0, minute=0)
    key_days = [sum_solstice, win_solstice, ver_eq, aut_eq]
    rep_days = [(x-year_start).days*24 for x in key_days]

    day_labels = ['Summer Solstice','Winter Solstice','Spring Eq', 'Fall Eq']
    key_days = dict.fromkeys(day_labels)
    frame_track = []
    for idx1, size in enumerate(result_df['design_size'].unique()):
        for idx, day in enumerate(rep_days):
            temp_df = result_df[result_df['design_size'] == size].copy()
            new_df = temp_df[day:day+24].copy()
            new_df.loc[:,'Hour'] = np.linspace(0,23,24)
            new_df['Key'] = day_labels[idx]
            frame_track.append(new_df)
    return pd.concat(frame_track)

def plot_2D(frame, x_var='design_size', y_var ='period_energy', x_label='Desired Size (kW)', y_label='Period Energy'):
    fig, ax = plt.subplots(figsize=(6,6))
    # ax2 = ax.twinx()

    ax.scatter(x=frame[x_var], y=frame[y_var])

    ax.set_xlabel(x_label, fontsize=18)
    ax.set_ylabel(y_label, fontsize=18)
    ax.tick_params(axis='x', labelsize = 18)
    ax.tick_params(axis='y', labelsize = 18)
    plt.locator_params(axis='x', nbins=5)
    plt.locator_params(axis='y', nbins=5)
    plt.show()

def get_training_validation(dataset_filename, sample_frac=1, training_fraction=0.8):
    print('Loading Training Data...\n')
    time_start = time.process_time()
    pkl_data = pd.read_pickle(dataset_filename)
    data = pkl_data.sample(n=int(sample_frac*len(pkl_data)))
    data_training, data_validation = split_training_validation(
        data, training_fraction, seed=len(data)
    )
    time_stop = time.process_time()
    print("Data Loading Time:", time_stop - time_start, "\n")
    return data_training, data_validation

def create_rbf_surrogate(
    training_dataframe, input_labels, output_labels, output_filename=None
):
    print('Creating Surrogate Model...\n')
    time_start = time.process_time()
    # Capture long output
    stream = StringIO()
    oldstdout = sys.stdout
    sys.stdout = stream

    # Create PySMO trainer object
    trainer = PysmoRBFTrainer(
        input_labels=input_labels,
        output_labels=output_labels,
        training_dataframe=training_dataframe,
    )

    # Set PySMO options
    trainer.config.basis_function = "gaussian"  # default = gaussian
    trainer.config.solution_method = "algebraic"  # default = algebraic
    trainer.config.regularization = True  # default = True

    # Train surrogate
    rbf_train = trainer.train_surrogate()

    # Remove autogenerated 'solution.pickle' file
    try:
        os.remove("solution.pickle")
    except FileNotFoundError:
        pass
    except Exception as e:
        raise e

    # Create callable surrogate object
    xmin, xmax = [0, 0], [1000000, 12]
    input_bounds = {
        input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
    }
    rbf_surr = PysmoSurrogate(rbf_train, input_labels, output_labels, input_bounds)

    # Save model to JSON
    if output_filename is not None:
        print(f'Writing surrogate model to {output_filename}')
        model = rbf_surr.save_to_file(output_filename, overwrite=True)

    # Revert back to standard output
    sys.stdout = oldstdout

    time_stop = time.process_time()
    print("Model Training Time:", time_stop - time_start, "\n")

    return rbf_surr

def load_surrogate(surrogate_filename=None):
    if surrogate_filename == None:
        surrogate_filename = join(dirname(__file__), "pv_multiperiod_surrogate.json")
    # Load surrogate model from file
    surrogate = PysmoSurrogate.load_from_file(surrogate_filename)

    return surrogate

def plot_parity(true_values, modeled_values, color='k', label=None, axes=None, title=None, R_val= None, fontsize = 16):
    data_scale = np.floor(np.log10(true_values.max()))
    data_bound = (10**data_scale)*np.ceil(true_values.max()/(10**data_scale))
    if axes == None:
        fig, ax = plt.subplots(figsize=(6,6))
    else:
        ax=axes
    ax.plot([0, data_bound], [0, data_bound], color='k', linewidth=2, linestyle='dashed', zorder=1)
    f = sns.scatterplot(x=true_values, y=modeled_values, ax=ax, s=50, alpha=0.7, facecolor=color, edgecolor='k', marker='o', zorder=2)
    if title!=None:
        ax.set_title(title, fontsize=fontsize)
    ax.set_xlim([0,data_bound])
    ax.set_ylim([0,data_bound])
    ax.set_xlabel(f'True Values', fontsize=fontsize)
    ax.set_ylabel(f'Model Values', fontsize=fontsize)
    ax.tick_params(axis='x', labelsize = fontsize)
    ax.tick_params(axis='y', labelsize = fontsize)
    ax.locator_params(axis='x', nbins=5)
    ax.locator_params(axis='y', nbins=5)
    ax.ticklabel_format(axis='both',style='scientific',useMathText=True,scilimits=(0,data_scale))

    ax.annotate(f'$R^2$ = {R_val}', (0.05*data_bound,0.9*data_bound), size=12)
    # ax.tight_layout()
    # plt.savefig(f'/Users/zbinger/watertap-seto/src/watertap_contrib/seto/solar_models/surrogate/plots/parity_plot_{label}.png', dpi=300)


def plot_training_validation(
    surrogate, data_training, data_validation, input_labels, output_labels, save=False
):
    fig, ax = plt.subplots(2,len(output_labels),figsize=(4*len(output_labels),8))
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    for idx, output_label in enumerate(output_labels):
        training_output = surrogate.evaluate_surrogate(data_training[input_labels])
        validation_output = surrogate.evaluate_surrogate(data_validation[input_labels])
        corr_matrix = np.corrcoef(data_validation[output_label], validation_output[output_label])
        corr = corr_matrix[0,1]
        r2_testing = round(corr**2,3)      
        r2_training=round(surrogate._trained._data[output_label].model.R2,3)
        rmse_training=round(surrogate._trained._data[output_label].model.rmse,3)

        label = re.sub("[^a-zA-Z0-9 \n\.]", " ", output_label.title())
        # Output fit metrics and create parity and residual plots
        print(f"{label}: \nR2: {r2_training} \nRMSE: {rmse_training}\n")
        if len(output_labels) > 1:
            train_ax = ax[0,idx]
            test_ax = ax[1,idx]
        else:
            train_ax = ax[0]
            test_ax = ax[1]
        plot_parity(
            true_values=np.array(data_training[output_label]),
            modeled_values=np.array(training_output[output_label]), color='#2A628F',
            label=output_label, axes=train_ax, title=f'{label}\n\nTraining', R_val=r2_training
        )
        plot_parity(
            true_values=np.array(data_validation[output_label]),
            modeled_values=np.array(data_validation[output_label]), color='#fc3c06',
            label=output_label, axes=test_ax, title='Testing', R_val=r2_testing
        )
    plt.tight_layout()
    if save == True:
        plt.savefig(f'/Users/zbinger/watertap-seto/src/watertap_contrib/seto/solar_models/surrogate/plots/parity_combined_2.png', dpi=900)
    plt.show()

def plot_period_parities(surrogate, data_training, data_validation, input_labels, selected_labels, output_label='period_energy'):
    fig, ax = plt.subplots(4,3,figsize=(8,10))
    fig.subplots_adjust(wspace=0.4, hspace=0.45)
    months = {0:'Jan',1:'Feb',2:'Mar',3:'Apr',4:'May',5:'Jun',6:'Jul',7:'Aug',8:'Sep',9:'Oct',10:'Nov',11:'Dec'}
    
    def sort(period):
        train = data_training.loc[data_training['month'] == period]
        val = data_validation.loc[data_training['month'] == period]
        return train, val
    
    def eval(train, val):
        train_eval = surrogate.evaluate_surrogate(train[input_labels])
        val_eval = surrogate.evaluate_surrogate(val[input_labels])
        corr_matrix = np.corrcoef(val[output_label], val_eval[output_label])
        corr = corr_matrix[0,1]
        r2_testing = round(corr**2,3)      
        r2_training=round(surrogate._trained._data[output_label].model.R2,3)
        rmse_training=round(surrogate._trained._data[output_label].model.rmse,3)
        return train_eval, val_eval, r2_testing, r2_training, rmse_training
    
    def plot(train, train_eval, idx, r2_training):
        plot_parity(
            true_values=np.array(train['period_energy']),
            modeled_values=np.array(train_eval['period_energy']), color='#2A628F',
            axes=ax[idx//3, idx%3], fontsize = 10,
            R_val=r2_training, title=months[idx]
        )
    
    for idx, period in enumerate(data_training['month'].unique()):
        train, val = sort(period)
        train_eval, val_eval, r2_testing, r2_training, rmse_training = eval(train, val)
        plot(train, train_eval, idx, r2_training)

    plt.tight_layout()
    plt.savefig(f'/Users/zbinger/watertap-seto/src/watertap_contrib/seto/solar_models/surrogate/plots/period_parity.png', dpi=900)
    plt.show()

def plot_overlap(true_values, modeled_values, color='#004E89', label=None, axes=None, title=None, R_val= None, fontsize = 16):
    data_scale = np.floor(np.log10(true_values['hourly_gen'].max()))
    data_bound = (10**data_scale)*np.ceil(true_values.max()/(10**data_scale))
    if axes == None:
        fig, ax = plt.subplots(figsize=(6,6))
    else:
        ax=axes
    # ax.plot([0, data_bound], [0, data_bound], color='k', linewidth=2, linestyle='dashed', zorder=1)
    f = sns.scatterplot(x=true_values['Hour'], y=true_values['hourly_gen'], ax=ax, s=50, alpha=0.7, facecolor='#4D4D4D', edgecolor='k', marker='o')
    f = sns.scatterplot(x=true_values['Hour'], y=modeled_values['hourly_gen'], ax=ax, s=50, alpha=0.7, facecolor=color, edgecolor='k', marker='o-')
    if title!=None:
        ax.set_title(title, fontsize=fontsize)
    ax.set_xlim([0,24])
    ax.set_ylim([0,100000])
    ax.set_xlabel(f'Hour', fontsize=fontsize)
    ax.set_ylabel(f'Energy', fontsize=fontsize)
    ax.tick_params(axis='x', labelsize = fontsize)
    ax.tick_params(axis='y', labelsize = fontsize)
    ax.locator_params(axis='x', nbins=5)
    ax.locator_params(axis='y', nbins=5)
    # ax.ticklabel_format(axis='both',style='scientific',useMathText=True,scilimits=(0,data_scale))

    # ax.annotate(f'$R^2$ = {R_val}', (0.05*data_bound,0.9*data_bound), size=12)

def plot_key_day_parities(idx, surrogate, data_training, data_validation, input_labels, selected_labels, axes=None, output_label='hourly_gen', title=''):
    ax=axes
    day_labels = ['Summer Solstice','Winter Solstice','Spring Eq', 'Fall Eq']
    
    def sort():
        train = data_training.copy()
        val = data_validation.copy()
        return train, val
    
    def eval(train, val):
        train_eval = surrogate.evaluate_surrogate(train[input_labels])
        val_eval = surrogate.evaluate_surrogate(val[input_labels])
        corr_matrix = np.corrcoef(val[output_label], val_eval[output_label])
        corr = corr_matrix[0,1]
        r2_testing = round(corr**2,3)      
        r2_training=round(surrogate._trained._data[output_label].model.R2,3)
        rmse_training=round(surrogate._trained._data[output_label].model.rmse,3)
        return train_eval, val_eval, r2_testing, r2_training, rmse_training
    
    train, val = sort()
    train_eval, val_eval, r2_testing, r2_training, rmse_training = eval(train, val)
    plot_parity(
            true_values=np.array(train['hourly_gen']),
            modeled_values=np.array(train_eval['hourly_gen']), color='#2A628F',
            axes=ax[0,idx], fontsize = 12,
            R_val=r2_training, title=title
        )
    
    plot_parity(
            true_values=np.array(data_validation[output_label]),
            modeled_values=np.array(val_eval[output_label]), color='#fc3c06',
            label=output_label, axes=ax[1,idx], fontsize = 12,
            R_val=r2_testing
        )

def create_yearly_surrogate():
        # Run parametrics via multiprocessing
    data = []
    pv_size = np.linspace(10,100000,5)
    df = pd.DataFrame()

    time_start = time.process_time()
    with multiprocessing.Pool(processes=8) as pool:
        # results = pool.map(run, pv_size)
        results = pool.map(run, pv_size)
    time_stop = time.process_time()
    print("Multiprocessing time:", time_stop - time_start, "\n")
    df = pd.concat(results)
    df.to_pickle(dataset_filename)
    print(f'Data written to {dataset_filename}')
    # plot_2D(df, y_var ='period_energy')
    # plot_2D(df, y_var ='land_req', y_label='Land Req (acre)')

    data_training, data_validation = get_training_validation(dataset_filename, sample_frac=1, training_fraction=0.9)

    surrogate_filename = join(dirname(__file__), "pv_multiperiod_surrogate.json")
    surrogate = create_rbf_surrogate(
        data_training, ['design_size', 'month'], output_labels, output_filename=surrogate_filename
    )

    # # Load surrogate model from file
    surrogate = load_surrogate()

    # Delete surrogate testing file
    # os.remove(surrogate_filename)

    selected_labels = ['period_energy', "land_req"]
    # Create parity and residual plots for training and validation
    plot_training_validation(
        surrogate, data_training, data_validation, ['design_size', 'month'], selected_labels
    )

    plot_period_parities(
    surrogate, data_training, data_validation, ['design_size', 'month'], selected_labels
    )

def create_daily_surrogate():
    # Run parametrics via multiprocessing
    data = []
    pv_size = np.linspace(10,100000,30) ################################
    pv_size = np.logspace(1,5,20) ################################
    df = pd.DataFrame()

    time_start = time.process_time()
    with multiprocessing.Pool(processes=8) as pool:
        # results = pool.map(run, pv_size)
        results = pool.map(get_hourly, pv_size)
    time_stop = time.process_time()
    print("Multiprocessing time:", time_stop - time_start, "\n")
    df = pd.concat(results)
    key_df = get_rep_days(df)
    # df.to_pickle(daily_dataset_filename)
    fig, ax = plt.subplots(2,4,figsize=(14,12))
    for idx, period in enumerate(key_df['Key'].unique()):
        period_df = key_df.loc[key_df['Key'] == period].copy()
        data = period_df.sample(n=int(1*len(period_df)))
        data_training, data_validation = split_training_validation(
            data, 0.95, seed=len(data) ######################################
        ) 

        surrogate_filename = join(dirname(__file__), "pv_"+period.replace(" ","_")+"_surrogate_w_land.json")
        surrogate = create_rbf_surrogate(
            data_training, ['design_size', 'Hour'], ['hourly_gen']
        )

        selected_labels = ['hourly_gen']
        # Create parity and residual plots for training and validation
        plot_key_day_parities(
            idx, surrogate, data_training, data_validation, ['design_size', 'Hour'], selected_labels, axes=ax, title=period
        )

    plt.tight_layout()
    plt.show()
        # plot_period_parities(
        # surrogate, data_training, data_validation, ['design_size', 'Hour'], selected_labels
        # )

if __name__ == "__main__":
    create_daily_surrogate()