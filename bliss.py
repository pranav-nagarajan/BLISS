import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import argparse
import time

from blimpy import Waterfall
from riptide import TimeSeries, ffa_search

start = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('filename', type = str, help = "Location of data file.")
parser.add_argument('--cutoff', type = int, default = 40, help = 'SNR cutoff value.')
parser.add_argument('--alias', type = int, default = 1, help = 'Number of periods to check for harmonics.')
args = parser.parse_args()
filename = args.filename
cutoff = args.cutoff
num_periods = args.alias

obs = Waterfall(filename)
nchans = obs.header['nchans']
data = np.squeeze(obs.data)

periods = []
frequencies = []
snrs = []

best_periods = []
for i in range(int(0.1 * nchans), int(0.9 * nchans)):

    time_series = TimeSeries.from_numpy_array(data[:, i], tsamp = obs.header['tsamp'])
    ts, pgram = ffa_search(time_series, rmed_width=4.0, period_min=0.01, period_max=1.25, bins_min=2, bins_max=260)

    if i == int(0.5 * nchans):
        pgram.display()

    best_periods.append(pgram.periods[np.argmax(pgram.snrs.max(axis=1))])

    mask = pgram.snrs.T[0] >= cutoff
    periods.extend(pgram.periods[mask])
    frequencies.extend(np.ones(len(pgram.periods[mask])) * obs.freqs[i])
    snrs.extend(pgram.snrs.T[0][mask])

ranked = pd.Series(np.round(np.array(best_periods), 4)).value_counts()
harmonics = np.zeros(len(periods), dtype = bool)
for r in ranked.keys()[0:num_periods]:
    for i in range(len(periods)):
        check = (round(periods[i] / r) > 1)
        close = (abs((periods[i] / r) - round(periods[i] / r)) <= 1e-3)
        if not harmonics[i] and check and close:
            harmonics[i] = True

end = time.time()
print('Best Period: ', ranked.keys()[0])
print('Time Taken: ', end - start)
