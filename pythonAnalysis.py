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
mpl.rcParams.update({'font.size':9})

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
ystart = 455
ystop  = 457

# ======================================
# Check number of points
# and initialize arrays
# ======================================
# Choose step range to perform analysis on
stepstart  = 2
stepend    = 8
numsteps   = stepend-stepstart
variance   = np.zeros(numsteps)
stddev     = np.zeros(numsteps)
varerr     = np.zeros(numsteps)
LGPS_3261  = np.zeros(numsteps)
LGPS_3311  = np.zeros(numsteps)
imgs       = imgs[stepstart:stepend]
stepvalues = stepvalues[stepstart:stepend]
# print varerr.shape
chisq_red = np.zeros(numsteps)

# ======================================
# Set up PDF
# ======================================
# Create PDF
pp = PdfPages('output.pdf')
# Create figure
fig=mt.figure('Page 1',figsize=(8.5,11))
# Create larger gridspec
outergs= gs.GridSpec(3,2)
# pdb.set_trace()
# ======================================
# Find spot size for each step
# ======================================
# mt.figure('Shot')
for i,img in enumerate(imgs):
	# img=np.flipud(np.rot90(f[img[0]],3))
	img=np.flipud(np.rot90(f[img[0]]))
	ax=fig.add_subplot(outergs[i])
	ax.imshow(img[350:600,xstart:xstop],aspect='auto',interpolation='none')
	fig.add_subplot(ax)
	plt.figure(fig.number)
	outergs.tight_layout(fig,pad=5)
	mt.addlabel(toplabel='$\Delta E$={}GeV'.format(stepvalues[i]),xlabel='x [px]',ylabel='y [px]')
	# plt.show()

	# Fit individual slices
	popt,pcov,chisq_red[i] = mt.fitimageslice(
			img,
			res_x,
			res_y,
			(xstart,xstop),
			(ystart,ystop),
			plot=True
			)
	mt.addlabel(toplabel='Gaussian Fit to Spot Profile, $\Delta E$={}GeV'.format(stepvalues[i]),xlabel='x [m]',ylabel='Counts')
	
	variance[i] = popt[2]
	bact = setQS.set_QS_energy_ELANEX(stepvalues[i])
	LGPS_3261[i] = setQS.bdes2K1(bact[0],20.35)
	LGPS_3311[i] = setQS.bdes2K1(bact[1],20.35)

	print 'QS1 K1: {}\tQS2 K1: {}'.format(LGPS_3261[i],LGPS_3311[i])

pp.savefig(fig)
pp.close()

# ======================================
# Debugging code
# ======================================
# img=imgs[3];
# img=np.flipud(np.rot90(f[img[0]]))
# plt.imshow(img,interpolation='none')
# mt.figure('Std. Dev.')
# plt.plot(stepvalues,np.sqrt(variance),'.-')

# ======================================
# Set up initial conditions
# ======================================
B5D36_en = 20.35
gamma    = (B5D36_en/0.5109989)*1e3
emitx    = 0.000100
twiss    = sltr.Twiss(
		beta  = 0.5,
		alpha = 0
		)

# ======================================
# Create beamlines
# ======================================
beamline=bt.beamlines.IP_to_lanex(twiss_x=twiss,twiss_y=twiss,gamma=gamma)
beamline_array = np.array([])
for i,beam in enumerate(stepvalues):
	beamline.elements[1].K1 = LGPS_3261[i]
	beamline.elements[2].K1 = LGPS_3261[i]
	beamline.elements[4].K1 = LGPS_3311[i]
	beamline.elements[5].K1 = LGPS_3311[i]
	beamline_array = np.append(beamline_array,copy.deepcopy(beamline))

# ======================================
# Fudge error
# ======================================
chisq_factor = 1e-28
# used_error   = stddev*np.sqrt(chisq_factor)
used_error   = variance*np.sqrt(chisq_factor)

# ======================================
# Fit beamline scan
# ======================================
out = bt.fitBeamlineScan(beamline_array,
		variance,
		emitx,
		error=used_error,
		verbose=True)

# ======================================
# Plot results
# ======================================
mpl.rcParams.update({'font.size':12})
bt.plotfit(stepvalues/20.35,
		variance,
		out.beta,
		out.X_unweighted,
		top='Emittance/Twiss Fit to Quad Scan',
		figlabel='Quad Scan Fit',
		error=used_error)
mt.addlabel(xlabel='$\delta$')

# figchisquare = plt.figure()
# mt.plot_featured(stepvalues,chisq_red,'.-',
#                 toplabel='Chi-Squared for Each Gaussian Fit',
#                 xlabel='$E/E_0$',
#                 ylabel='$\chi^2$')

# plt.show()
