# Breakthrough Listen Periodic Spectral Signals Search

The Breakthrough Listen project aims to conduct a comprehensive search for radio technosignatures from potential extraterrestrial intelligence sources. Different classes of technosignatures require different search strategies. In the past, efforts have focused on searching for continuous narrowband signals that display Doppler drift due to relative motion between the source and the observer [1]. However, broadband signals with periodic spectral modulation patterns are also a promising indicator of ETI, considering that similar signals are common in communication systems on Earth.

Unfortunately, algorithms dedicated to search for such signals are computationally intensive. Furthermore, they typically rely on Fast Fourier Transform techniques, which are impacted by "red" (random walk) noise and perform poorly on modulation patterns with either a small duty cycle (<1%) or a long period (>500 milliseconds). This work aims to develop a pipeline for searching for periodic spectral signals based on the Fast Folding Algorithm, which displays improved performance in these regimes [2].

The initial stage of the project will be to write a python wrapper that uses the blimpy [3] and riptide-ffa [4] packages, along with necessary plotting tools, to search for and analyze periodic spectral signals across all frequency channels in a single set of 5-minute observations. The next stage of the project involves developing tools to compare a set of ON and OFF observations for a sky-localized signal. The final stage of the project will be to test the pipeline on data collected by the BL program, in one of the first comprehensive searches for this class of signals.

## Repository Organization

- `bliss.py`
  - The main script for this project, which carries out the FFA analysis using blimpy and riptide. It accepts as input an ON file and any number of optional OFF files, and outputs a text file containing all periodicities detected in each frequency channel above a certain SNR threshold.

- `comparison.py`
  - A script containing functionality for quick comparison of an ON and OFF file prior to analysis. The two techniques tested were based on windowing and kurtosis. This script is not currently used in the pipeline.

- `plotter.py`
  - The secondary script for this project, containing all necessary tools to generate a 6x2 plot for each detected periodic signal. The left column contains the periodograms for each file in the cadence, while the right column contains the corresponding folded profiles. 

- `simulation.py`
  - A script containing functionality for injecting periodic signals into an ON file, and testing to see if they are recovered by the FFA algorithm. Since this functionality has now been included as an option in the main script, this code is currently not used in the pipeline.

## FFA Algorithm

## Benchmarking

## Results

## Next Steps

## References

1. Enriquez, J. Emilio, et al. "The Breakthrough Listen search for intelligent life: 1.1–1.9 GHz observations of 692 nearby stars." https://arxiv.org/pdf/1709.03491.pdf
2. Morello, Vincent, et al. "Optimal periodicity searching: Revisiting the Fast Folding Algorithm for large-scale pulsar surveys." https://arxiv.org/pdf/2004.03701.pdf
3. https://blimpy.readthedocs.io/en/latest/
4. https://pypi.org/project/riptide-ffa/
