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


def prep_data(file):
    """Generates time-averaged spectrum."""
    obs = Waterfall(file)
    data = np.squeeze(obs.data)
    freqs = np.array([obs.header['fch1'] + i * obs.header['foff'] for i in range(obs.header['nchans'])])

    average = data[0].copy()
    for i in range(len(data)):
        average += data[i].copy()
    average = average / len(data)

    return freqs, average


def calc_window(freqs, spectrum):
    """Finds means and standard deviations of 2.9 MHz windows of time-averaged spectrum."""
    current, start_freq, end_freq = min(freqs), min(freqs), max(freqs)
    store, means, sds = [], [], []
    index, counter, running = 0, 0, 0
    diff, flag = 0, False

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
            if not flag:
                diff = current - start_freq
                flag = True
            start_freq = current

        current += freqs[0] - freqs[1]

    means.append(running / counter)
    sds.append(np.std(store))
    return means, sds, diff


def find_outlier(freqs, spectrum, means, sds, diff, cutoff):
    """Identifies signals by comparing to noise background in 2.9 MHz window."""
    outliers = []
    for i in range(len(freqs)):
        binned = int((freqs[i] - min(freqs)) / diff)
        z_score = (spectrum[i] - means[binned]) / sds[binned]
        outliers.append(z_score > cutoff)
    return outliers


on_freqs, on_average = prep_data(on_file)
print("Progress: Read ON file.")
off_freqs, off_average = prep_data(off_file)
print("Progress: Read OFF file.")

on_means, on_sds, on_diff = calc_window(on_freqs, on_average)
on_outliers = find_outlier(on_freqs, on_average, on_means, on_sds, on_diff, cutoff)
print("Progress: Processed ON file.")

off_means, off_sds, off_diff = calc_window(off_freqs, off_average)
off_outliers = find_outlier(on_freqs, off_average, off_means, off_sds, off_diff, cutoff)
print("Progress: Processed OFF file.")

ignore = []
for i in range(len(on_outliers)):
    if on_outliers[i] and off_outliers[i]:
        ignore.append(freqs[i])
print("Progress: ON-OFF comparison complete.")

fig, ax = plt.subplots(2, 1, sharex = True)
ax[0].plot(on_freqs, on_average)
ax[1].plot(off_freqs, off_average)
ax[0].vlines(ignore, min(ax[0].get_ylim()), max(ax[0].get_ylim()), color = 'red')
ax[1].vlines(ignore, min(ax[1].get_ylim()), max(ax[1].get_ylim()), color = 'red')
plt.title('ON-OFF Comparison')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Power')
plt.savefig('comparison.png')
plt.show()

np.savetxt('ignore.txt', np.array(ignore))
end = time.time()
print('Time Taken: ', end - start)
