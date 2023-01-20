import os
import pandas as pd

def read_weather_data(filename):
    """Read a TMY3 weather datafile to a dict that can be directly set to the ssc table"""
    # Read header
    df = pd.read_csv(
        os.path.join(os.path.dirname(__file__), filename),
        sep=',',
        header=0,
        nrows=1,
        skipinitialspace=True,
        engine='python'
        )
    data = {}
    data['tz'] = int(df['Time Zone'][0])
    data['elev'] = float(df['Elevation'][0])
    data['lat'] = float(df['Latitude'][0])
    data['lon'] = float(df['Longitude'][0])

    # Read timeseries data
    df = pd.read_csv(
        os.path.join(os.path.dirname(__file__), filename),
        sep=',',
        skiprows=2,
        header=0,
        skipinitialspace=True
        )
    data['year'] = list(df['Year'])
    data['month'] = list(df['Month'])
    data['day'] = list(df['Day'])
    data['hour'] = list(df['Hour'])
    data['minute'] = list(df['Minute']) if 'Minute' in df.keys() else [0.0 for j in data['hour']]
    data['dn'] = list(df['DNI'])
    data['df'] = list(df['DHI'])
    data['gh'] = list(df['GHI'])
    data['wspd'] = list(df['Wind Speed'])
    data['tdry'] = list(df['Temperature'])
    data['pres'] = list(df['Pressure'])
    data['tdew'] = list(df['Dew Point'])
    
    return data
