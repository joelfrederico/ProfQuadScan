#!/usr/bin/env python

# System imports
import pdb
import argparse
import copy
import h5py as h5

# Math/science imports
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gs
import matplotlib as mpl
mpl.rcParams.update({'font.size':7})

# My module imports
import ButterflyEmittancePython as bt
import mytools.slactrac as sltr
import mytools as mt
import setQS

plt.close('all')

# ======================================
# Load and transfer matlab variables
# ======================================
infile     = 'tempfiles/forpython.mat'
f          = h5.File(infile);
data       = f['data']
imgs       = data['processed']['images']['CEGAIN']['dat']
stepvalues = data['raw']['scalars']['step_value']['dat']
stepvalues = mt.derefdataset(stepvalues,f)
stepvalues = np.unique(stepvalues)

# ======================================
# Set up image slices
# ======================================
res_y    = 10.3934e-6
res_x    = res_y / np.sqrt(2)
xstart = 550
xstop  = 700
# Get y center
ycent = np.linspace(510,415,6)
ystart = ycent - 25
ystop  = ycent + 25

# ======================================
# Check number of points
# and initialize arrays
# ======================================
# Choose step range to perform analysis on
stepstart  = 2
stepend    = 8
numsteps   = stepend-stepstart
pxmin = np.zeros(numsteps)
imgs       = imgs[stepstart:stepend]
stepvalues = stepvalues[stepstart:stepend]

# ======================================
# Set up PDF
# ======================================
# Create PDF
pp = PdfPages('output.pdf')
# Create figure
fig=mt.figure('Page 1',figsize=(8.5,11))
# Create larger gridspec
outergs= gs.GridSpec(4,2)
# pdb.set_trace()
# ======================================
# Find spot size for each step
# ======================================
# mt.figure('Shot')
for i,img in enumerate(imgs):
	# img=np.flipud(np.rot90(f[img[0]],3))
	img=np.flipud(np.rot90(f[img[0]]))
	ax=fig.add_subplot(outergs[i])
	ax.imshow(img,interpolation='none')
	fig.add_subplot(ax)
	# fig2=mt.figure('test')
	# ax2=fig2.add_subplot(111)
	# ax2.plot(np.linspace(1,10))
	# ax = fig.add_subplot(outergs[i+1])
	# ax.imshow(img)
	# plt.show()

	# pdb.set_trace()

	# Fit individual slices
	pxmin[i] = mt.findpinch(
			img,
			xbounds=(xstart,xstop),
			ybounds=(ystart[i],ystop[i]),
			step=2,
			verbose=False)
	
# ======================================
# Debugging code
# ======================================
# img=imgs[3];
# img=np.flipud(np.rot90(f[img[0]]))
# plt.imshow(img,interpolation='none')
# mt.figure('Std. Dev.')
# plt.plot(stepvalues,np.sqrt(variance),'.-')



outergs.tight_layout(fig)

ymin = pxmin*res_y*1e3

out=np.polyfit(stepvalues/20.35,ymin,1)

fig2=mt.plot_featured(stepvalues,ymin,'.-',stepvalues,np.polyval(out,stepvalues/20.35),'-',
		figlabel = 'Fit',
		toplabel='Dispersion Measurement: {}mm'.format(out[0]),
		xlabel = 'Pinch location (y) [mm]',
		ylabel = 'Energy offset [GeV]',
		legend= ('Data','Fit')
		)
# mt.addlabel('Dispersion Fit

pp.savefig(fig)
pp.savefig(fig2)
pp.close()

plt.show()
