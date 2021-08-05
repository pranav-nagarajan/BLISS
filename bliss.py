import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import pandas as pd

import argparse
import time
from tqdm import tqdm
import multiprocessing as mp

from blimpy import Waterfall
from riptide import TimeSeries, ffa_search

start = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('signal', type = str, help = 'Location of ON data file.')
parser.add_argument('--background', action = "append", type = str, default = None, help = 'Location of OFF data file.')
parser.add_argument('--cutoff', type = int, default = 10, help = 'SNR cutoff value.')
parser.add_argument('--alias', type = int, default = 1, help = 'Number of periods to check for harmonics.')
parser.add_argument('--range', type = float, nargs = 2, default = [2.5, 10], help = 'Period range for FFA search.')
parser.add_argument('--multi', action = "store_true", help = 'Use multiprocessing.')
parser.add_argument('--simulate', action = "store_true", help = 'Turns on simulation of fake signal.')
parser.add_argument('--beam', action = "store_true", help = 'Creates a four-digit code summarizing ON-OFF comparison.')
parser.add_argument('--output', type = str, default = "signal.txt", help = 'Name of output file.')
args = parser.parse_args()
on_file = args.signal
off_files = args.background
cutoff = args.cutoff
num_periods = args.alias
period_range = args.range
multi = args.multi
simulate = args.simulate
beam = args.beam
output = args.output


def periodic_analysis(on_data, off_data, freqs, nchans, tsamp, start, stop, cutoff):
    """Periodic analysis on a set of ON-OFF files."""
    periods = []
    frequencies = []
    snrs = []
    best_periods = []
    indicators = []
    best_widths = []
    min_widths = []
    max_widths = []
    all_codes = []

    for i in range(int(start * nchans), int(stop * nchans)):

        if sum(on_data[:, i]) == 0:
            continue
        on_periods, on_freqs, sn_ratios, best_period, widths, mini, maxi = periodic_helper(on_data[:, i], freqs[i], tsamp, cutoff)
        if off_data is not None:
            indicator = np.zeros(len(on_periods))
            codes = np.zeros(len(on_periods), dtype = str)
            for j in range(len(off_data)):
                datum = off_data[j]
                off_periods = periodic_helper(datum[:, i], freqs[i], tsamp, cutoff, False)
                if beam:
                    indicator, codes = compare_on_off(on_periods, off_periods, indicator, codes)
                    prev = '' + '0' * j
                    if prev in codes:
                        codes = [s + '0' for s in codes]
                else:
                    indicator = compare_on_off(on_periods, off_periods, indicator)

        periods.extend(on_periods)
        frequencies.extend(on_freqs)
        snrs.extend(sn_ratios)
        best_periods.append(best_period)
        best_widths.extend(widths)
        min_widths.extend(mini)
        max_widths.extend(maxi)
        if off_data is not None:
            indicators.extend(indicator)
            if beam:
                all_codes.extend(codes)

    if off_data is None:
        return periods, frequencies, snrs, best_periods, best_widths, min_widths, max_widths
    else:
        if not beam:
            return periods, frequencies, snrs, best_periods, best_widths, min_widths, max_widths, indicators
        return periods, frequencies, snrs, best_periods, best_widths, min_widths, max_widths, indicators, all_codes


def periodic_helper(data, frequency, tsamp, cutoff, on = True):
    """Runs FFA on data from a blimpy Waterfall object."""
    periods = []
    frequencies = []
    snrs = []
    best_periods = []

    time_series = TimeSeries.from_numpy_array(data, tsamp = tsamp)
    if on and simulate and frequency in injection:
        time_series = time_series.normalise()
        fts = TimeSeries.generate(length=len(data)*tsamp, tsamp=tsamp, period=99.0, ducy=0.02, amplitude=100.0).normalise()
        time_series = TimeSeries.from_numpy_array(time_series.data + fts.data, tsamp = tsamp).normalise()
    ts, pgram = ffa_search(time_series, rmed_width=4.0, period_min=period_range[0], period_max=period_range[1], bins_min=2, bins_max=260)
    mask = pgram.snrs.max(axis = 1) >= cutoff
    periods = pgram.periods[mask]

    if not on:
        return periods
    else:
        frequencies = np.ones(len(pgram.periods[mask])) * frequency
        snrs = pgram.snrs.max(axis = 1)[mask]
        best_period = pgram.periods[np.argmax(pgram.snrs.max(axis=1))]
        best_widths = [pgram.widths[i] for i in np.argmax(pgram.snrs, axis = 1)[mask]]
        min_widths = [min(pgram.widths) for w in best_widths]
        max_widths = [max(pgram.widths) for w in best_widths]
        return periods, frequencies, snrs, best_period, best_widths, min_widths, max_widths


def compare_on_off(on_periods, off_periods, indicator, codes = None):
    """"Compares ON and OFF files."""
    counter, tol = 0, 1e-4
    for i in range(len(on_periods)):
        po = on_periods[i]
        for j in range(len(off_periods)):
            pf = off_periods[j]
            if abs(po - pf) <= tol:
                if indicator[counter] == 0:
                    indicator[counter] = 1
                if codes is not None:
                    codes[counter] += '1'
            else:
                if codes is not None:
                    codes[counter] += '0'
        counter += 1
    if codes is not None:
        return indicator, codes
    return indicator


def find_harmonics(periods, best_periods, num_periods):
    """Finds and labels harmonics in periodograms."""
    harmonics = np.zeros(len(periods), dtype = bool)

    ranked = pd.Series(np.round(np.array(best_periods), 4)).value_counts()
    inspect = ranked.keys()[0:num_periods]
    counts = np.zeros(len(inspect))
    for i in range(len(inspect)):
        for j in range(len(periods)):
            check = (round(periods[j] / inspect[i]) > 1)
            close = (abs((periods[j] / inspect[i]) - round(periods[i] / inspect[i])) <= 1e-3)
            if check and close:
                counts[i] += 1
                if not harmonics[j]:
                    harmonics[j] = True

    return ranked, counts, harmonics


def concat_helper(results):
    """Concatenates results from worker processes."""
    periods = []
    frequencies = []
    snrs = []
    best_periods = []
    widths = []
    min_widths = []
    max_widths = []
    indicators = []
    codes = []

    for package in results:
        periods.extend(package[0])
        frequencies.extend(package[1])
        snrs.extend(package[2])
        best_periods.extend(package[3])
        widths.extend(package[4])
        min_widths.extend(package[5])
        max_widths.extend(package[6])
        if off_files is not None:
            indicators.extend(package[7])
            if beam:
                codes.extend(package[8])

    final = [np.array(periods), np.array(frequencies), np.array(snrs), np.array(best_periods)]
    final.extend([np.array(widths), np.array(min_widths), np.array(max_widths)])
    if off_files is not None:
        final.append(np.array(indicators))
        if beam:
            final.append(np.array(codes))
    return final


def plot_helper(periods, frequencies, snrs, harmonics, indicators):
    """Plots frequency channel vs. periodogram."""

    full_signal = list(zip(periods, frequencies, snrs))
    filter = np.zeros(len(full_signal), dtype = bool)
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
                filter[i] = 1
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

    return filter, signal


obs = Waterfall(on_file)
data = np.squeeze(obs.data)

freqs = np.array([obs.header['fch1'] + i * obs.header['foff'] for i in range(obs.header['nchans'])])
nchans, tsamp = obs.header['nchans'], obs.header['tsamp']

if simulate:
    injection = np.random.choice(freqs, 10, replace = False)
print("Progress: Read ON file.")

background_data = None
if off_files is not None:
    background_data = []
    for off_file in off_files:
        background = Waterfall(off_file)
        back_data = np.squeeze(background.data)
        background_data.append(back_data)
    print("Progress: Read OFF files.")

if multi:
    pool = mp.Pool(mp.cpu_count())
    on_iterables = [(data, background_data, freqs, nchans, tsamp, 0.1 * i, 0.1 * (i + 1), cutoff) for i in range(1, 9)]
    on_results = pool.starmap(periodic_analysis, on_iterables)
    on_results = concat_helper(on_results)
else:
    on_results = periodic_analysis(data, background_data, freqs, nchans, tsamp, 0.1, 0.9, cutoff)
    on_results = concat_helper([on_results])

ranked, counts, harmonics = find_harmonics(on_results[0], on_results[3], num_periods)
print("Progress: File processing complete.")

if off_files is not None:
    filter, signal = plot_helper(on_results[0], on_results[1], on_results[2], harmonics, on_results[-2])
else:
    filter, signal = plot_helper(on_results[0], on_results[1], on_results[2], harmonics, np.zeros(len(harmonics)))

for i in tqdm(range(len(signal))):
    aliasing, round_period = 0, round(signal[i][0], 4)
    if round_period in ranked.keys()[0:num_periods]:
        aliasing = counts[np.where(ranked.keys() == round_period)]
    signal[i] += (np.where(freqs == signal[i][1])[0][0], abs(obs.header['foff']), aliasing, period_range[0], period_range[1])
    signal[i] += (on_results[4][filter][i], on_results[5][filter][i], on_results[6][filter][i])
    signal[i] += (obs.header['source_name'], obs.header['tstart'], on_file)
    if beam:
        signal[i] += (on_results[-1][i],)
np.savetxt(output, signal, fmt = "%s")

if multi:
    pool.close()
    pool.join()

end = time.time()
print('Best Period: ', ranked.keys()[0])
print('Time Taken: ', end - start)
