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
parser.add_argument('--cutoff', type = int, default = 20, help = 'SNR cutoff value.')
parser.add_argument('--alias', type = int, default = 1, help = 'Number of periods to check for harmonics.')
args = parser.parse_args()
on_file = args.signal
off_file = args.background
cutoff = args.cutoff
num_periods = args.alias


def periodic_analysis(on_data, off_data, freqs, nchans, tsamp, start, stop, cutoff):
    """Periodic analysis on a set of ON-OFF files."""
    periods = []
    frequencies = []
    snrs = []
    best_periods = []
    indicators = []

    for i in range(int(start * nchans), int(stop * nchans)):

        on_periods, on_freqs, sn_ratios, best_period = periodic_helper(on_data[:, i], freqs[i], tsamp, cutoff)
        if off_data is not None:
            off_periods = periodic_helper(off_data[:, i], freqs[i], tsamp, cutoff, False)
            indicator = compare_on_off(on_periods, off_periods)

        periods.extend(on_periods)
        frequencies.extend(on_freqs)
        snrs.extend(sn_ratios)
        best_periods.append(best_period)

        if off_data is not None:
            indicators.extend(indicator)

    if off_data is None:
        return periods, frequencies, snrs, best_periods
    return periods, frequencies, snrs, best_periods, indicators


def periodic_helper(data, frequency, tsamp, cutoff, on = True):
    """Runs FFA on data from a blimpy Waterfall object."""
    periods = []
    frequencies = []
    snrs = []
    best_periods = []

    time_series = TimeSeries.from_numpy_array(data, tsamp = tsamp)
    ts, pgram = ffa_search(time_series, rmed_width=4.0, period_min=0.01, period_max=10, bins_min=2, bins_max=260)
    mask = pgram.snrs.T[0] >= cutoff
    periods = pgram.periods[mask]

    if not on:
        return periods
    else:
        frequencies = np.ones(len(pgram.periods[mask])) * frequency
        snrs = pgram.snrs.T[0][mask]
        best_period = pgram.periods[np.argmax(pgram.snrs.max(axis=1))]
        return periods, frequencies, snrs, best_period


def compare_on_off(on_periods, off_periods):
    """"Compares ON and OFF files."""
    counter, tol = 0, 1e-4
    indicator = np.zeros(len(on_periods))
    for i in range(len(on_periods)):
        po = on_periods[i]
        for j in range(len(off_periods)):
            pf = off_periods[j]
            if indicator[counter] == 0 and abs(po - pf) <= tol:
                indicator[counter] = 1
        counter += 1
    return indicator


def find_harmonics(periods, best_periods, num_periods):
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


def plot_helper(periods, frequencies, snrs, harmonics, indicators):
    """Plots frequency channel vs. periodogram."""

    full_signal = list(zip(periods, frequencies, snrs))
    signal = []
    alias = []
    background = []
    for i in range(len(harmonics)):
        if harmonics[i]:
            alias.append(full_signal[i])
        else:
            if indicators[i]:
                background.append(full_signal[i])
            else:
                signal.append(full_signal[i])

    cmap = plt.cm.viridis
    norm = matplotlib.colors.Normalize(vmin = min(snrs), vmax = max(snrs))

    plt.figure(figsize = (8, 6))
    if len(signal) > 0:
        plt.scatter(list(zip(*signal))[0], list(zip(*signal))[1], c = cmap(norm(list(zip(*signal))[2])), marker = 'o')
    if len(alias) > 0:
        plt.scatter(list(zip(*alias))[0], list(zip(*alias))[1], c = cmap(norm(list(zip(*alias))[2])), marker = '+')
    if len(background) > 0:
        plt.scatter(list(zip(*background))[0], list(zip(*background))[1], c = cmap(norm(list(zip(*background))[2])), marker = '^')

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
nchans, tsamp = obs.header['nchans'], obs.header['tsamp']
print("Progress: Read ON file.")

#pool = mp.Pool(mp.cpu_count())

if off_file is not None:

    background = Waterfall(off_file)
    background_data = np.squeeze(background.data)
    on_iterables = [(data, background_data, freqs, nchans, tsamp, 0.1 * i, 0.1 * (i + 1), cutoff) for i in range(1, 9)]
    print("Progress: Read OFF file.")

else:

    on_iterables = [(data, None, freqs, nchans, tsamp, 0.1 * i, 0.1 * (i + 1), cutoff) for i in range(1, 9)]

#on_results = pool.starmap(periodic_analysis, on_iterables)
#on_results = concat_helper(on_results)
on_results = periodic_analysis(data, None, freqs, nchans, tsamp, 0.1, 0.9, cutoff)
ranked, harmonics = find_harmonics(on_results[0], on_results[3], num_periods)
print("Progress: File processing complete.")

if off_file is not None:
    plot_helper(on_results[0], on_results[1], on_results[2], harmonics, on_results[4])
else:
    plot_helper(on_results[0], on_results[1], on_results[2], harmonics, np.zeros(len(harmonics)))

#pool.close()
#pool.join()

end = time.time()
print('Best Period: ', ranked.keys()[0])
print('Time Taken: ', end - start)
