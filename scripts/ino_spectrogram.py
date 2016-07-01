#! /usr/bin/env python

import numpy as np
import os
import h5py
import sys
import scipy.signal
import argparse

def minor_format(x, i=None):
    return x

parser = argparse.ArgumentParser(description='Generating CSR Spectrogramms from Inovesa result files')
parser.add_argument('directory', type=str, help='relativ path the Inovesa hdf5 files')
parser.add_argument('--filenameending', type=str, default='b.h5', help='only files with the specified filenameending are used')
parser.add_argument('--saveplot', type=str, nargs='?', default=False, const='replace by args.directory')
parser.add_argument('--saveformat', type=str, default='png')
parser.add_argument('--showplot', action='store_true')
parser.add_argument('--datalen', type=int, default=10000, help='spezifies how many simulations steps (starting from the end of the file) will be used')
parser.add_argument('--currentmin', type=float,default=None)
parser.add_argument('--currentmax', type=float,default=None)
parser.add_argument('--freqmin', type=float,default=None)
parser.add_argument('--freqmax', type=float,default=None)
parser.add_argument('--colmin', type=float,default=None)
parser.add_argument('--colmax', type=float,default=None)

args = parser.parse_args()
args.saveplot = args.directory if args.saveplot=='replace by args.directory' else args.saveplot

import matplotlib

if not args.showplot:
  matplotlib.use('Agg')

if args.saveformat == 'eps':
  matplotlib.use('PS')

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import colormaps as cmaps
plt.register_cmap(name='inferno', cmap=cmaps.inferno)


dataname='/CSR/Intensity/data'
timename='/Info/AxisValues_t'

fnames = []

for fname in os.listdir(args.directory):
  if fname[-4:] == args.filenameending:
    fnames.append(fname)

fnames.sort()

data = []
currents = []
deltat = 1

for fname in fnames:
  h5file1 = h5py.File(args.directory+'/'+fname, 'r')
  tmp=deltat
  deltat = h5file1['/Info/AxisValues_t'].attrs['Factor4Seconds']*h5file1['/Info/AxisValues_t'][1]  
  assert tmp==deltat or tmp==1, 'not same deltat in all files (at file %s)' %fname
  data.append(np.abs(np.fft.rfft(h5file1[dataname][-2*args.datalen:]))[1:])
  currents.append(h5file1['/Info/Parameters'].attrs['BunchCurrent'])
  h5file1.close()

scaling = 0.31621
freqs = np.linspace(0,1/(2*deltat),args.datalen)/1000.
freqs_scaled = np.linspace(0,scaling/(2*deltat),args.datalen)/1000.

data = np.array(data)
currents = np.array(currents)
currents = currents*1000
data = (data.T*np.power(currents,2)).T
data=data[currents.argsort(),:]
currents.sort()


#plt.rc('text', usetex=True)
#font = {'family' : 'normal',  'weight' : 'bold', 'size' : 16 }
font = {'size' : 19 }
matplotlib.rc('font', **font)

plt.figure(figsize=(10,5),tight_layout=True)
plt.xlabel("Frequency / kHz")
plt.ylabel("Bunch Current / mA")
spectrogram = plt.pcolormesh(freqs,currents,data,shading='flat',norm=matplotlib.colors.LogNorm(vmin=10,vmax=5e7))
plt.yscale('log', subsy=np.arange(1, 10, 1)[1:])
plt.gca().yaxis.set_minor_formatter(FuncFormatter(minor_format))
plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x, i=None: x))
plt.tick_params('x',which='major', direction='out', labelleft='on', width=2,length=5)
plt.tick_params('y',which='minor', direction='inout', labelleft='on', width=1,length=5)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom') 
plt.xlim((args.freqmin/scaling if args.freqmin else freqs[0], args.freqmax/scaling if args.freqmax else freqs[-1]))
plt.ylim((args.currentmin/scaling if args.currentmin else currents[0],args.currentmax/scaling if args.currentmax else currents[-1]))
plt.minorticks_on()
plt.colorbar(spectrogram).set_label("Spectral Intensity / a.u.")
spectrogram.set_cmap('inferno')

if args.saveplot:
  plt.savefig("%s/simulated-spectogram-isoring.png" %(args.directory),dpi=200)

currents = currents*scaling


plt.figure(figsize=(10,5),tight_layout=True)
plt.xlabel("Scaled Frequency / kHz")
plt.ylabel("Scaled Bunch Current / mA")
spectrogram = plt.pcolormesh(freqs_scaled,currents,data,shading='flat',norm=matplotlib.colors.LogNorm(vmin=10,vmax=5e7))
plt.yscale('log', subsy=np.arange(1, 10, 1)[1:])
plt.gca().yaxis.set_minor_formatter(FuncFormatter(minor_format))
plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x, i=None: x))
plt.tick_params('x',which='major', direction='out', labelleft='on', width=2,length=5)
plt.tick_params('y',which='minor', direction='inout', labelleft='on', width=1,length=5)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom') 
plt.xlim((args.freqmin if args.freqmin else freqs_scaled[0], args.freqmax if args.freqmax else freqs_scaled[-1]))
plt.ylim((args.currentmin if args.currentmin else currents[0],args.currentmax if args.currentmax else currents[-1]))
plt.minorticks_on()
plt.colorbar(spectrogram).set_label("Spectral Intensity / a.u.")
spectrogram.set_cmap('inferno')

if args.saveplot:
  plt.savefig("%s/simulated-spectogram.png" %args.directory,dpi=200)

plt.figure(figsize=(10,5),tight_layout=True)
plt.xlabel("Scaled Frequency / kHz")
plt.ylabel("Scaled Bunch Current / mA")
spectrogram = plt.pcolormesh(freqs_scaled,currents,data,shading='flat',norm=matplotlib.colors.LogNorm(vmin=10,vmax=5e8))
plt.xscale('log', subsy=np.arange(1, 10, 1)[1:])
plt.yscale('log', subsy=np.arange(1, 10, 1)[1:])
plt.gca().yaxis.set_minor_formatter(FuncFormatter(minor_format))
plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x, i=None: x))
plt.tick_params('x',which='major', direction='out', labelleft='on', width=2,length=5)
plt.tick_params('y',which='minor', direction='inout', labelleft='on', width=1,length=5)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom') 
plt.xlim((args.freqmin if args.freqmin else 0.07,args.freqmax if args.freqmax else freqs_scaled[-1]))
plt.ylim((args.currentmin if args.currentmin else currents[0],args.currentmax if args.currentmax else currents[-1]))
plt.minorticks_on()
plt.colorbar(spectrogram).set_label("Spectral Intensity / a.u.")
spectrogram.set_cmap('inferno')

if args.saveplot:
  plt.savefig("%s/simulated-spectogram_logscale.png" %args.directory,dpi=200)

if args.showplot:
  plt.show()
else:
  plt.close()

