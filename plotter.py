import numpy as np
import matplotlib.pyplot as plt
import argparse
import time

from blimpy import Waterfall
from riptide import TimeSeries, ffa_search

start = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('--inputs', action = "append", type = str, help = "Location of text files.")
parser.add_argument('--signal', action = "append", type = str, help = "Location of 'ON' data.")
parser.add_argument('--background', action = "append", type = str, help = "Location of 'OFF' data.")
parser.add_argument('--range', type = float, nargs = 2, default = [2.5, 10], help = 'Period range for FFA search.')
parser.add_argument('--bins', type = int, default = 10, help = "Number of bins for folded profile.")
parser.add_argument('--beam', action = store_true,  help = "Creates a six-digit code summarizing ON-OFF comparison.")
args = parser.parse_args()
inputs = args.inputs
on_files = args.signal
off_files = args.background
period_range = args.range
num_bins = args.bins
beam = args.beam

full_input = []
for name in inputs:
    full_input.extend(np.loadtxt(name, dtype = str))
input = list(set([tuple(arr) for arr in full_input]))

if beam:
    codes = []
    for package in input:
        flags = ''
        for name in inputs:
            curr = np.loadtxt(name, dtype = str)
            if package in curr:
                flags += '1'
            else:
                flags += '0'
        semi = package[-1]
        codes.append("".join(map(lambda x, y: x + y, flags, semi)))


signal_data = []
for on_file in on_files:
    obs = Waterfall(on_file)
    data = np.squeeze(obs.data)
    signal_data.append(data)
print("Progress: Read ON files.")

background_data = []
for off_file in off_files:
    background = Waterfall(off_file)
    back_data = np.squeeze(background.data)
    background_data.append(back_data)
print("Progress: Read OFF files.")

freqs = np.array([obs.header['fch1'] + i * obs.header['foff'] for i in range(obs.header['nchans'])])
nchans, tsamp = obs.header['nchans'], obs.header['tsamp']

flag = False
for package in input:

    fig, axes = plt.subplots(len(signal_data) + len(background_data), 2, figsize = (20, 20))
    channel =  np.where(freqs == float(package[1]))[0][0]

    sig_index, back_index, counter = 0, 0, 0
    for i in range(len(signal_data) + len(background_data)):

        if counter % 2 == 0:
            plot_data = signal_data[sig_index][:, channel]
            plot_label = 'ON'
            sig_index += 1
        else:
            plot_data = background_data[back_index][:, channel]
            plot_label = 'OFF'
            back_index += 1

        rts = TimeSeries.from_numpy_array(plot_data, tsamp = tsamp)
        ts, pgram = ffa_search(rts, rmed_width=4.0, period_min=period_range[0], period_max=period_range[1], bins_min=2, bins_max=260)
        folded = ts.fold(float(package[0]), bins = num_bins, subints = 1)

        axes[counter][0].plot(pgram.periods, pgram.snrs.max(axis = 1))
        axes[counter][0].set_xlabel('Period (sec)')
        axes[counter][0].set_ylabel('SNR')

        axes[counter][1].plot(np.linspace(0, float(package[0]), num_bins), folded)
        axes[counter][1].set_xlabel('Time (sec)')
        axes[counter][1].set_ylabel('Power')

        axes[counter][0].set_xlim(float(package[0]) - 1.0, float(package[0]) + 1.0)
        axes[counter][0].set_title(plot_label + ' Periodogram')
        axes[counter][1].set_title(plot_label + ' Spectrum')

        plt.tight_layout()
        counter += 1

    if not flag:
        plt.savefig('plot.png')
        if beam:
            np.savetxt('codes.txt', codes)
        flag = True
        end = time.time()
        print('Time Taken: ', end - start)
        break

    plt.show()
