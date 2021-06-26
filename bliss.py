import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import pandas as pd

import argparse
import time

from blimpy import Waterfall
from riptide import TimeSeries, ffa_search

start = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('signal', type = str, help = "Location of 'ON' data file.")
parser.add_argument('background', type = str, help = "Location of 'OFF' data file.")
parser.add_argument('--cutoff', type = int, default = 40, help = 'SNR cutoff value.')
parser.add_argument('--alias', type = int, default = 1, help = 'Number of periods to check for harmonics.')
args = parser.parse_args()
on_file = args.signal
off_file = args.background
cutoff = args.cutoff
num_periods = args.alias

obs = Waterfall(on_file)
data = np.squeeze(obs.data)
freqs = np.array([obs.header['fch1'] + i * obs.header['foff'] for i in range(obs.header['nchans'])])

periods = []
frequencies = []
snrs = []

best_periods = []
for i in range(int(0.1 * obs.header['nchans']), int(0.9 * obs.header['nchans'])):

    time_series = TimeSeries.from_numpy_array(data[:, i], tsamp = obs.header['tsamp'])
    ts, pgram = ffa_search(time_series, rmed_width=4.0, period_min=0.01, period_max=1.25, bins_min=2, bins_max=260)

    best_periods.append(pgram.periods[np.argmax(pgram.snrs.max(axis=1))])

    mask = pgram.snrs.T[0] >= cutoff
    periods.extend(pgram.periods[mask])
    frequencies.extend(np.ones(len(pgram.periods[mask])) * freqs[i])
    snrs.extend(pgram.snrs.T[0][mask])

ranked = pd.Series(np.round(np.array(best_periods), 4)).value_counts()
harmonics = np.zeros(len(periods), dtype = bool)
for r in ranked.keys()[0:num_periods]:
    for i in range(len(periods)):
        check = (round(periods[i] / r) > 1)
        close = (abs((periods[i] / r) - round(periods[i] / r)) <= 1e-3)
        if not harmonics[i] and check and close:
            harmonics[i] = True

periods = np.array(periods)
frequencies = np.array(frequencies)
snrs = np.array(snrs)

background = Waterfall(off_file)
back_data = np.squeeze(background.data)
back_freqs = np.array([background.header['fch1'] + i * background.header['foff'] for i in range(background.header['nchans'])])

back_periods = []
back_frequencies = []

for i in range(int(0.1 * background.header['nchans']), int(0.9 * background.header['nchans'])):

    time_series = TimeSeries.from_numpy_array(back_data[:, i], tsamp = background.header['tsamp'])
    ts, pgram = ffa_search(time_series, rmed_width=4.0, period_min=0.01, period_max=1.25, bins_min=2, bins_max=26)

    mask = pgram.snrs.T[0] >= cutoff
    back_periods.extend(pgram.periods[mask])
    back_frequencies.extend(np.ones(len(pgram.periods[mask])) * back_freqs[i])

back_periods = np.array(back_periods)
back_frequencies = np.array(back_frequencies)

indicator = np.zeros(len(periods))
counter = 0
for (po, fo) in zip(periods, frequencies):
    for (pb, fb) in zip(back_periods, back_frequencies):
        if indicator[counter] == 0 and abs(po - pb) <= 1e-3 and abs(fo - fb) <= 1e-3:
            indicator[counter] = 1
    counter += 1

full_signal = list(zip(periods, frequencies, snrs))
on = []
alias = []
off = []
for i in range(len(harmonics)):
    if harmonics[i]:
        alias.append(full_signal[i])
    else:
        if indicator[i]:
            off.append(full_signal[i])
        else:
            on.append(full_signal[i])

cmap = plt.cm.viridis
norm = matplotlib.colors.Normalize(vmin = min(snrs), vmax = max(snrs))

plt.figure(figsize = (8, 6))
if len(on) > 0:
    plt.scatter(list(zip(*on))[0], list(zip(*on))[1], c = cmap(norm(list(zip(*on))[2])), marker = 'o')
if len(alias) > 0:
    plt.scatter(list(zip(*alias))[0], list(zip(*alias))[1], c = cmap(norm(list(zip(*alias))[2])), marker = '+')
if len(off) > 0:
    plt.scatter(list(zip(*off))[0], list(zip(*off))[1], c = cmap(norm(list(zip(*off))[2])), marker = '^')
cbar = plt.colorbar(plt.cm.ScalarMappable(cmap = cmap, norm = norm))
plt.xlabel('Periods')
plt.ylabel('Frequencies')
cbar.set_label('SNR', rotation = 270)
plt.savefig('output.png')
plt.show();

end = time.time()
print('Best Period: ', ranked.keys()[0])
print('Time Taken: ', end - start)
