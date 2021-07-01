import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import pandas as pd

import argparse
import time
import multiprocessing as mp

from blimpy import Waterfall
from riptide import TimeSeries, ffa_search

start = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('signal', type = str, help = "Location of 'ON' data file.")
parser.add_argument('--background', type = str, default = None, help = "Location of 'OFF' data file.")
parser.add_argument('--cutoff', type = int, default = 40, help = 'SNR cutoff value.')
parser.add_argument('--alias', type = int, default = 1, help = 'Number of periods to check for harmonics.')
args = parser.parse_args()
on_file = args.signal
off_file = args.background
cutoff = args.cutoff
num_periods = args.alias


def periodic_analysis(data, freqs, start, stop, cutoff, on = True):
    """Runs FFA on a blimpy Waterfall object."""
    periods = []
    frequencies = []
    snrs = []
    best_periods = []

    for i in range(int(start * obs.header['nchans']), int(stop * obs.header['nchans'])):

        time_series = TimeSeries.from_numpy_array(data[:, i], tsamp = obs.header['tsamp'])
        ts, pgram = ffa_search(time_series, rmed_width=4.0, period_min=0.01, period_max=100, bins_min=2, bins_max=260)
        if on:
            best_periods.append(pgram.periods[np.argmax(pgram.snrs.max(axis=1))])

        mask = pgram.snrs.T[0] >= cutoff
        periods.extend(pgram.periods[mask])
        frequencies.extend(np.ones(len(pgram.periods[mask])) * freqs[i])
        if on:
            snrs.extend(pgram.snrs.T[0][mask])

    if not on:
        return periods, frequencies
    return periods, frequencies, snrs, best_periods


def find_harmonics(best_periods, num_periods):
    """Finds and labels harmonics in periodograms."""
    harmonics = np.zeros(len(periods), dtype = bool)

    ranked = pd.Series(np.round(np.array(best_periods), 4)).value_counts()
    for r in ranked.keys()[0:num_periods]:
        for i in range(len(periods)):
            check = (round(periods[i] / r) > 1)
            close = (abs((periods[i] / r) - round(periods[i] / r)) <= 1e-3)
            if not harmonics[i] and check and close:
                harmonics[i] = True

    return ranked, harmonics


def compare_on_off(periods, frequencies, back_periods, back_frequencies):
    """Compares ON and OFF files."""
    indicator = np.zeros(len(periods))

    counter = 0
    for (po, fo) in zip(periods, frequencies):
        for (pb, fb) in zip(back_periods, back_frequencies):
            if indicator[counter] == 0 and abs(po - pb) <= 1e-3 and abs(fo - fb) <= 1e-3:
                indicator[counter] = 1
        counter += 1

    return indicator


def concat_helper(results):
    """Concatenates results from worker processes."""
    periods = []
    frequencies = []
    snrs = []
    best_periods = []

    for package in results:
        periods.extend(package[0])
        frequencies.extend(package[1])
        snrs.extend(package[2])
        best_periods.extend(package[3])

    return [np.array(periods), np.array(frequencies), np.array(snrs), np.array(best_periods)]


def plot_helper(periods, frequencies, snrs, harmonics, indicator):
    """Plots frequency channel vs. periodogram."""

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

    return full_signal


obs = Waterfall(on_file)
data = np.squeeze(obs.data)
freqs = np.array([obs.header['fch1'] + i * obs.header['foff'] for i in range(obs.header['nchans'])])
print("Progress: Read ON file.")

pool = mp.Pool(mp.cpu_count())
on_iterables = [(data, freqs, 0.1 * i, 0.1 * (i + 1), cutoff, True) for i in range(1, 9)]
on_results = pool.starmap(periodic_analysis, on_iterables)
on_results = concat_helper(on_results)
print("Progress: Processed ON file.")

ranked, harmonics = find_harmonics(on_results[3], num_periods)

if off_file is not None:

    background = Waterfall(off_file)
    background_data = np.squeeze(background.data)
    background_freqs = np.array([background.header['fch1'] + i * background.header['foff'] for i in range(background.header['nchans'])])
    print("Progress: Read OFF file.")

    off_iterables = [(background_data, background_freqs, 0.1 * i, 0.1 * (i + 1), cutoff, False) for i in range(1, 9)]
    off_results = pool.map(periodic_analysis, off_iterables)
    off_results = concat_helper(off_results)
    print("Progress: Processed OFF file.")

    indicator = compare_on_off(on_results[0], on_results[1], off_results[0], off_results[1])

else:

    indicator = np.zeros(len(on_results[0]))

plot_helper(on_results[0], on_results[1], on_results[2], harmonics, indicator)

end = time.time()
print('Best Period: ', ranked.keys()[0])
print('Time Taken: ', end - start)
