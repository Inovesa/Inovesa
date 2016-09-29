#! /usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import h5py
import sys

fname = sys.argv[1]

hdf_f = h5py.File(fname, 'r')

inovesa_version = hdf_f['/Info/INOVESA_v']

current = (hdf_f['/Info/Parameters'].attrs['BunchCurrent'])
grid_size = (hdf_f['/Info/Parameters'].attrs['GridSize'])
time = hdf_f['/Info/AxisValues_t'][...]
bunchlength = hdf_f['/BunchLength/data'][...]
energyspread = hdf_f['/EnergySpread/data'][...]
csr_power = hdf_f['/CSR/Intensity/data'][...]
csr_factor = 1e-6
csr_label = r"CSR Intensity / a.u."

if inovesa_version >= 12:
    csr_factor = (hdf_f['/CSR/Intensity/data']).attrs['Factor4Watts']
    csr_label = r"CSR Intensity / W"

plt.figure()
ax1 = plt.subplot(111)
ax1.set_xlabel(r"Time / $T_s$")
ax1.set_ylabel(r"$\sigma_{p,q}$")
ax1.plot(time,bunchlength,label=r"$\sigma_q$")
ax1.plot(time,energyspread,label=r"$\sigma_p$")

ax2 = ax1.twinx()

ax2.plot(time,csr_power*csr_factor,"r-",label="CSR")
ax2.set_ylabel(csr_label)
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2,bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)


if inovesa_version >= 13:
    tracks = hdf_f['/Particles/data'][...]
    if tracks[1] > 0:
        plt.figure()
        plt.axes().set_aspect('equal')
        for i in range(tracks.shape[1]):
       	    plt.plot(tracks[:,i,0]/grid_size,tracks[:,i,1]/grid_size)


plt.show()
