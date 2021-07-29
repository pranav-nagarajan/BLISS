import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import argparse
import multiprocessing as mp
import time

from blimpy import Waterfall
from riptide import TimeSeries, ffa_search

start = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('signal', type = str, help = "Location of 'ON' data file.")
parser.add_argument('background', type = str, help = "Location of 'OFF' data file.")
args = parser.parse_args()
on_file = args.signal
off_file = args.background


def simulate_period(on_data, off_data, freqs, nchans, tsamp, start, stop, trial, channels):
    """Periodic analysis on a set of ON-OFF files with simulated periodic signal."""
    best_periods = []
    indicators = []

    for i in range(int(start * nchans), int(stop * nchans)):
        rts = TimeSeries.from_numpy_array(on_data[:, i], tsamp = tsamp).normalise()
        if freqs[i] in channels:
            fts = TimeSeries.generate(length=len(rts.data) * tsamp, tsamp=tsamp, period=trial, ducy=0.02, amplitude=20.0).normalise()
            rts = TimeSeries.from_numpy_array(rts.data + fts.data, tsamp = tsamp).normalise()
        ts, pgram = ffa_search(rts, rmed_width=4.0, period_min=2.5, period_max=10, bins_min=2, bins_max=260)
        best_period = pgram.periods[np.argmax(pgram.snrs.max(axis=1))]
        best_periods.append(best_period)

        ots = TimeSeries.from_numpy_array(off_data[:, i], tsamp = tsamp).normalise()
        ots, opgram = ffa_search(ots, rmed_width=4.0, period_min=2.5, period_max=10, bins_min=2, bins_max=260)
        off_period = opgram.periods[np.argmax(opgram.snrs.max(axis=1))]
        indicators.append(abs(best_period - off_period) > 1e-4)

    return np.array(best_periods), np.array(indicators)


def simulate_helper(trial):
    """Carries out simulation for a single trial period."""
    on_iterables = [(data, back_data, freqs, nchans, tsamp, 0.1 * i, 0.1 * (i + 1), trial, injection) for i in range(1, 9)]
    on_results = pool.starmap(simulate_period, on_iterables)

    best_periods = []
    indicators = []
    for package in on_results:
        best_periods.extend(package[0])
        indicators.extend(package[1])

    ranked = pd.Series(np.round(np.array(best_periods)[indicators], 4)).value_counts()
    best_period = ranked.keys()[0]
    return best_period


obs = Waterfall(on_file)
data = np.squeeze(obs.data)
freqs = np.array([obs.header['fch1'] + i * obs.header['foff'] for i in range(obs.header['nchans'])])
nchans, tsamp = obs.header['nchans'], obs.header['tsamp']
print("Progress: Read ON file.")

background = Waterfall(off_file)
back_data = np.squeeze(background.data)
injection = np.random.choice(freqs, int(len(data[:, 0]) * 0.55), replace = False)
print("Progress: Read OFF file.")

pool = mp.Pool(mp.cpu_count())

best_periods = []
for trial in np.linspace(1, 10, 10):
    best_periods.append(simulate_helper(trial))
best_periods = np.array(best_periods)
print("Progress: File processing complete.")

pool.close()
pool.join()

plt.figure()
plt.scatter(np.linspace(1, 10, 10), best_periods)
plt.xlabel('Injected Period')
plt.ylabel('Recovered Period')
plt.savefig('simulation.png')
plt.show()

end = time.time()
print('Time Taken: ', end - start)
