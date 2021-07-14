import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import time
from blimpy import Waterfall

start = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('signal', type = str, help = "Location of 'ON' data file.")
parser.add_argument('background', type = str, help = "Location of 'OFF' data file.")
parser.add_argument('--cutoff', type = int, default = 10, help = "SNR cutoff value.")
args = parser.parse_args()
on_file = args.signal
off_file = args.background
cutoff = args.cutoff


def prep_data(file, return_freqs = False):
    """Generates time-averaged spectrum."""
    obs = Waterfall(on_file)
    data = np.squeeze(obs.data)

    average = data[0]
    for i in range(len(data)):
        average += data[i]
    average = average / len(data)

    if return_freqs:
        freqs = np.array([obs.header['fch1'] + i * obs.header['foff'] for i in range(obs.header['nchans'])])
        return average, freqs
    return average


def calc_window(freqs, spectrum):
    """Finds means and standard deviations of 2.9 MHz windows of time-averaged spectrum."""
    current, start_freq, end_freq = min(freqs), min(freqs), max(freqs)
    store, means, sds = [], [], []
    index, counter, running = 0, 0, 0

    while current <= end_freq:

        running += spectrum[index]
        store.append(spectrum[index])
        counter += 1
        index += 1

        if current - start_freq > 2.9:
            means.append(running / counter)
            sds.append(np.std(store))
            counter, running = 0, 0
            on_store, off_store = [], []
            start_freq = current

        current += freqs[1] - freqs[0]

    return means, sds


def find_outlier(freqs, spectrum, means, sds, cutoff):
    """Identifies signals by comparing to noise background in 2.9 MHz window."""
    outliers = []
    for i in range(len(freqs)):
        bin = int((freqs[i] - min(freqs)) / 2.9)
        z_score = (spectrum[i] - means[bin]) / sds[bin]
        outliers.append(z_score > cutoff)
    return outliers


on_average, freqs = prep_data(on_file, True)
off_average = prep_data(off_file)

on_means, on_sds = calc_window(freqs, on_average)
off_means, off_sds = calc_window(freqs, off_average)

on_outliers = find_outlier(freqs, on_average, on_means, on_sds, cutoff)
off_outliers = find_outlier(freqs, off_average, off_means, off_sds, cutoff)

ignore = []
for i in range(len(on_outliers)):
    if on_outliers[i] and off_outliers[i]:
        ignore.append(freqs[i])

fig, ax = plt.subplots(2, 1, sharex = True)
ax[0].plot(freqs, on_average)
ax[1].plot(freqs, off_average)
ax[0].vlines(ignore, min(ax[0].get_ylim()), max(ax[0].get_ylim()), color = 'red')
ax[1].vlines(ignore, min(ax[1].get_ylim()), max(ax[1].get_ylim()), color = 'red')
plt.show()

np.save('ignore.txt', np.array(ignore))
end = time.time()
print('Time Taken: ', end - start)
