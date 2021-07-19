import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import kurtosis

import argparse
import time
import multiprocessing as mp
from blimpy import Waterfall

start = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('signal', type = str, help = "Location of 'ON' data file.")
parser.add_argument('background', type = str, help = "Location of 'OFF' data file.")
parser.add_argument('--cutoff', type = int, default = 20, help = "SNR cutoff value.")
args = parser.parse_args()
on_file = args.signal
off_file = args.background
cutoff = args.cutoff


def prep_data(data, start, stop):
    """Generates time-averaged spectrum."""
    beg, end = int(start * len(data[0])), int(stop * len(data[0]))
    avg = data.mean(axis = 0)[beg:end]

    kurts = []
    for i in range(beg, end):
        kurt = kurtosis(data[:, i], nan_policy = 'omit')
        kurts.append(kurt)

    kurto = np.array(kurts)
    return avg, kurto


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

def concat_helper(results):
    """Concatenates results from worker processes."""
    averages = []
    kurtoses = []
    for package in results:
        averages.extend(package[0])
        kurtoses.extend(package[1])
    return np.array(averages), np.array(kurtoses)


pool = mp.Pool(mp.cpu_count())

signal = Waterfall(on_file)
on_data = np.squeeze(signal.data)
on_freqs = np.array([signal.header['fch1'] + i * signal.header['foff'] for i in range(signal.header['nchans'])])
on_iterables = [(on_data, 0.1 * k, 0.1 * (k + 1)) for k in range(10)]
on_results = pool.starmap(prep_data, on_iterables)
on_average, on_kurtosis = concat_helper(on_results)

plt.figure()
plt.hist(on_kurtosis)
plt.xlabel('Kurtosis')
plt.ylabel('Count')
plt.savefig('kurtosis.png')
print("Progress: Read ON file.")

background = Waterfall(off_file)
off_data = np.squeeze(background.data)
off_freqs = np.array([background.header['fch1'] + j * background.header['foff'] for j in range(background.header['nchans'])])
off_iterables = [(off_data, 0.1 * n, 0.1 * (n + 1)) for n in range(10)]
off_results = pool.starmap(prep_data, off_iterables)
off_average, off_kurtosis = concat_helper(off_results)
print("Progress: Read OFF file.")

# on_means, on_sds, on_diff = calc_window(on_freqs, on_average)
# on_outliers = find_outlier(on_freqs, on_average, on_means, on_sds, on_diff, cutoff)
on_outliers = (on_kurtosis >= cutoff)
print("Progress: Processed ON file.")

# off_means, off_sds, off_diff = calc_window(off_freqs, off_average)
# off_outliers = find_outlier(on_freqs, off_average, off_means, off_sds, off_diff, cutoff)
off_outliers = (off_kurtosis >= cutoff)
print("Progress: Processed OFF file.")

ignore = []
for i in range(len(on_outliers)):
    if on_outliers[i] and off_outliers[i]:
        ignore.append(on_freqs[i])
print("Progress: ON-OFF comparison complete.")

fig, ax = plt.subplots(2, 1, sharex = True)
ax[0].plot(on_freqs, on_average)
ax[1].plot(off_freqs, off_average)
ax[0].vlines(ignore, min(ax[0].get_ylim()), max(ax[0].get_ylim()), color = 'red')
ax[1].vlines(ignore, min(ax[1].get_ylim()), max(ax[1].get_ylim()), color = 'red')
fig.suptitle('ON-OFF Comparison')
plt.xlabel('Frequency (MHz)')
ax[0].set_ylabel('Power')
ax[1].set_ylabel('Power')
plt.savefig('comparison.png')
plt.show()

pool.close()
pool.join()

np.savetxt('ignore.txt', np.array(ignore))
end = time.time()
print('Time Taken: ', end - start)
