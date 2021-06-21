import numpy as np
import matplotlib.pyplot as plt
import argparse
import time
from blimpy import Waterfall
from riptide import TimeSeries, ffa_search

start = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('--filename', type = str, help = "Location of data file.")
args = parser.parse_args()
filename = args.filename

obs = Waterfall(filename)
nchans = obs.header['nchans']
data = np.squeeze(obs.data)

best_periods = []
for i in range(int(0.1 * nchans), int(0.9 * nchans)):
    channeled = data[:, i]
    time_series = riptide.TimeSeries.from_numpy_array(channeled, tsamp = obs.header['tsamp'])
    ts, pgram = riptide.ffa_search(time_series, rmed_width=4.0, period_min=0.01, period_max=1.25, bins_min=2, bins_max=260)
    best_periods.append(pgram.periods[np.argmax(pgram.snrs.max(axis=1))])

ranked = pd.Series(np.round(np.array(best_periods), 4)).value_counts()
best_period = ranked.keys()[0]

end = time.time()
print('Best Period: ', best_period)
print('Time Taken: ', end - start)
