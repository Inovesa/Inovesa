#! /usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import h5py
import sys

quantities = sys.argv[1]

fnames = sys.argv[2:]
fnames = fnames[::-1]
print(fnames)

inovesa_version = 0 

current = []
grid_size = []
time = np.array([])
bunchlength = np.array([])
energyspread = np.array([])
csr_power = np.array([])
csr_factor = 1e-6
csr_label = r"CSR Intensity (arb. units)"

for fname in fnames:
    hdf_f = h5py.File(fname, 'r')

    try:
        inovesa_version = hdf_f['/Info/Inovesa_v']
    except:
        inovesa_version = hdf_f['/Info/INOVESA_v']

    current.append(hdf_f['/Info/Parameters'].attrs['BunchCurrent'])
    grid_size = (hdf_f['/Info/Parameters'].attrs['GridSize'])
    if (len (time) > 0):
        time = np.append(time,hdf_f['/Info/AxisValues_t'][...] + time[-1])
    else:
        time = np.append(time,hdf_f['/Info/AxisValues_t'][...])
    bunchlength = np.append(bunchlength,hdf_f['/BunchLength/data'][...])
    energyspread = np.append(energyspread,hdf_f['/EnergySpread/data'][...])
    csr_power = np.append(csr_power,hdf_f['/CSR/Intensity/data'][...])
    if inovesa_version[1] >= 14:
        csr_factor = (hdf_f['/CSR/Intensity/data']).attrs['Watt']
        csr_label = r"CSR Intensity (W)"

plt.figure()
ax1 = plt.subplot(111)
ax1.set_xlabel(r"Time ($T_s$)")
ax1.set_ylabel(r"$\sigma_{p,q}$")
if ("sq" in quantities):
    ax1.plot(time,bunchlength,label=r"$\sigma_q$")
if ("sp" in quantities):
    ax1.plot(time,energyspread,label=r"$\sigma_p$")

h1, l1 = ax1.get_legend_handles_labels()


if ("sr" in quantities):
    ax2 = ax1.twinx()

    ax2.plot(time,csr_power*csr_factor,"r-",label="CSR")
    ax2.set_ylabel(csr_label)
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2,bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0.)
else:
    ax1.legend(h1, l1,bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0.)



if len(fnames) == 1 and inovesa_version[1] >= 13:
    tracks = hdf_f['/Particles/data'][...]
    if (not 0 in tracks.shape):
        plt.figure()
        plt.axes().set_aspect('equal')
        for i in range(tracks.shape[1]):
            if i < 16:
       	        plt.plot(tracks[:,i,0]/grid_size,tracks[:,i,1]/grid_size)


plt.show()
