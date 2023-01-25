#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 09:52:29 2021

@author: Katrine O. Bangsgaard, DTU
Modified version based on script from 
Evelina Ametova, The University of Manchester and Karlsruhe Institute of Technology (KIT), 
evelina.ametova@kit.edu

This file is distributed under the Apache 2.0 license available from
https://www.apache.org/licenses/LICENSE-2.0
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.io
#%% parse spectral info
# shutter values
tof_lim_left = np.array([15e-3, 27e-3, 44e-3, 53e-3])
tof_lim_right = np.array([26.68e-3, 43.68e-3, 52.68e-3, 72e-3])
tof_bin_width = np.array([10.24, 20.48, 20.48, 40.96])

# number of shutter intervals
n_intervals = 4

# binning in each shutter interval
binning = np.array([16, 8, 8, 4])

# calculate number of bins in each shutter interval
# TOF is in seconds, bins in microseconds
n_bins = np.int_(np.floor((tof_lim_right - tof_lim_left) / (tof_bin_width * binning * 1e-6)))
n_bins_total = np.sum(n_bins)

# prepare spectral axis for plots
tof_bins_left = np.zeros((n_bins_total), dtype = np.float32)
tof_bins_right = np.zeros((n_bins_total), dtype = np.float32)
counter = 0
for i in range(n_intervals):
    tof_bins_left[counter:(counter+n_bins[i])] = np.arange(tof_lim_left[i], tof_lim_right[i]-tof_bin_width[i]*binning[i]*1e-6, tof_bin_width[i]*binning[i]*1e-6, dtype = np.float32)
    tof_bins_right[counter:(counter+n_bins[i])] = tof_bins_left[counter:(counter+n_bins[i])] + (tof_bin_width[i]*binning[i]*1e-6)
    counter = counter+n_bins[i]

tof_bins_center = ((tof_bins_left + tof_bins_right) / 2)

l = 56.428

angstrom_bins_center = (tof_bins_center * 3957) / l


scipy.io.savemat("angstrom_bins.mat", {'angstrom':angstrom_bins_center})






