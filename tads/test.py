# test.py
# -------------------------------------------
# testing purposes

from scipy import signal
import numpy as np
import matplotlib.pyplot as plt

def negate(x):
    return -x

xs = np.arange(0, np.pi, 0.05)
data = np.sin(xs)
peakind = signal.find_peaks_cwt(data, np.arange(1,10))
peakind, xs[peakind], data[peakind]

data2 = map(negate, data)