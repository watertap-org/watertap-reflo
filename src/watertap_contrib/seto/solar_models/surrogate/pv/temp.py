import os
import time
import pandas as pd
from idaes.core.surrogate.sampling.data_utils import split_training_validation

dataset_filename = os.path.join(
    os.path.dirname(__file__), "data/pv_data.pkl"
)

sample_fraction = 1
training_fraction = 0.8

print('Loading Training Data...\n')
time_start = time.process_time()
pkl_data = pd.read_pickle(dataset_filename)
data = pkl_data.sample(n=int(sample_fraction*len(pkl_data)))
data_training, data_validation = split_training_validation(
    data, training_fraction, seed=len(data)
)
time_stop = time.process_time()
print("Data Loading Time:", time_stop - time_start, "\n")

print(data_training.head())