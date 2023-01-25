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
# load packages
import numpy as np
from astropy.io import fits
import scipy.io


# load angles
angles = range(0, 120)

## parse specral info
# shutter values
tof_lim_left = np.array([15e-3, 27e-3, 44e-3, 53e-3])
tof_lim_right = np.array([26.68e-3, 43.68e-3, 52.68e-3, 72e-3])
tof_bin_width = np.array([10.24, 20.48, 20.48, 40.96])

# number of shutter intervals
n_intervals = 4

# binning in each shutter interval
#binning = np.array([16, 8, 8, 4])
binning = np.array([16, 8, 8, 4])

# calculate number of bins in each shutter interval
# TOF is in seconds, bins in microseconds
n_bins = np.int_(np.floor((tof_lim_right - tof_lim_left) / (tof_bin_width * binning * 1e-6)))
n_bins_total = np.sum(n_bins)


# path to projections and flat images
pathname_projection = './data/tomo_nov_2019_for_meeting/'
projection_prefix = 'angle{}'
projection_channel_prefix = 'IMAT000{}_tomolg_000_'
projection_counter = 13470

pathname_flat_after = './data/tomo_nov_2019_for_meeting/'
flat_after_prefix = 'OpenBeam{}'
flat_after_channel_prefix = 'IMAT000{}_tomolg_000_'
flat_after_counter = 13591

pathname_flat_before = './data/tomo_nov_2019_for_meeting/'
flat_before_prefix = 'OpenBeam{}'
flat_before_channel_prefix = 'IMAT000{}_tomolg_000_'
flat_before_counter = 13466

image_format = 'fits'


# image geometry
n_angles   = 120
imsize     = 512
n_channels = n_bins_total
n_flat     = 4
cor_offset = 28
n_pix      = 460

#%% image geometry

# initialize structures
sino_raw = np.zeros((n_channels, n_angles, imsize+cor_offset), dtype = np.float32)
sino     = np.zeros((n_channels, n_angles, n_pix), dtype = np.float32)
flat_raw = np.zeros((n_channels, 2*n_flat, imsize+cor_offset), dtype = np.float32)
flat     = np.zeros((n_channels, 2*n_flat, n_pix), dtype = np.float32)

for i in range(n_channels):
    
    print(i)
    # flatfields before and after
    flat_before = np.zeros((imsize, imsize), dtype = np.float32)
    n_trigs_flat_before = 0
    
    flat_after = np.zeros((imsize, imsize), dtype = np.float32)
    n_trigs_flat_after = 0    
    for j in range(n_flat):
        # load flat field before 
        filename = (pathname_flat_before + 
                    flat_before_prefix + '/' + 
                    flat_before_channel_prefix+
                    '{:05d}_cor.fits').format(j + 1, flat_before_counter + j, i);
        flat_handler = fits.open(filename)
        tmp = np.array(flat_handler[0].data, dtype = np.float32)
                   
        # center of roation correction
        flat_raw[i, j, cor_offset:]  = tmp[:, 127] /flat_handler[0].header['N_TRIGS']
        flat[i,j,:] = flat_raw[i, j,40:500];
        
        # store mean flatfield for flux normalization
        flat_before += tmp
        n_trigs_flat_before += flat_handler[0].header['N_TRIGS']
    
    
        # load flat field after
        filename = (pathname_flat_after + 
                    flat_after_prefix + '/' + 
                    flat_after_channel_prefix+
                    '{:05d}_cor.fits').format(j + 1 + 4, flat_after_counter + j, i);
        flat_handler = fits.open(filename)
        tmp = np.array(flat_handler[0].data, dtype = np.float32)            
        
        # center of roation correction
        flat_raw[i, 4+j, cor_offset:]  = tmp[:, 127] /flat_handler[0].header['N_TRIGS']
        flat[i,4+j,:] = flat_raw[i, 4+j,40:500];
        
        
        # store mean flatfield for flux normalization
        flat_after += tmp
        n_trigs_flat_after += flat_handler[0].header['N_TRIGS']
    
    # normalize triggers
    flat_before /= n_trigs_flat_before        
    flat_after /= n_trigs_flat_after
    
    # projections 
    for j in range(n_angles):
        # load raw sinogram
        filename_projection = (pathname_projection + 
                               projection_prefix + '/' + 
                               projection_channel_prefix +
                               '{:05d}_cor.fits').format(j, projection_counter + j, i);    
        projection_handler = fits.open(filename_projection)
        projection         = np.array(projection_handler[0].data, dtype = np.float32)            

        # normalize triggers
        n_trigs_projection = flat_handler[0].header['N_TRIGS']
        projection        /= n_trigs_projection
        
        
        # flux normalization
        current_flat       = (flat_after + flat_before) / 2
        current_projection = projection[:, 127]
        coeff              = np.mean(np.mean(current_flat[6:10, :]))  / np.mean(np.mean(projection[6:10, :]))
        current_flat       = current_flat[:, 127]
        tmp                = np.zeros((imsize), dtype = np.float32)
             
        # flux normalization
        tmp[np.logical_and(current_flat > 0, current_projection > 0)] = coeff*current_projection[np.logical_and(current_flat > 0, current_projection > 0)] 
        

        # center of rotation correction
        sino_raw[i, j, cor_offset:] = tmp
        sino[i,j,:] = sino_raw[i, j,40:500];
        
                

scipy.io.savemat("neutrondata.mat", {'sino':sino,'flat':flat})













