#!/usr/bin/env python
import argparse
import numpy as np
import scipy.io as sio
import ButterflyEmittancePython as bt
import mytools.slactrac as sltr
import matplotlib.pyplot as plt
import mytools as mt
import copy
import h5py as h5
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

res    = 10.3934
xstart = 455
xstop  = 457
ystart = 550
ystop  = 700
step   = 2
x_range = slice(xstart,xstop)
y_range = slice(ystart,ystop)

# ======================================
# Check number of points
# and initialize arrays
# ======================================
num_pts=(xstop-xstart)/step
sigs=np.zeros(num_pts)

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
# Find spot size for each step
# ======================================
mt.figure('Shot')
for i,img in enumerate(imgs):
	img=np.flipud(np.rot90(f[img[0]]))
	# plt.imshow(img,interpolation='none')

	# Fit individual slices
	popt,pcov,chisq_red[i] = mt.fitimageslice(img,res/1e6,x_range,y_range)
	
	variance[i] = popt[2]
	bact = setQS.set_QS_energy_ELANEX(stepvalues[i])
	LGPS_3261[i] = setQS.bdes2K1(bact[0],20.35)
	LGPS_3311[i] = setQS.bdes2K1(bact[1],20.35)

	print 'QS1 K1: {}\tQS2 K1: {}'.format(LGPS_3261[i],LGPS_3311[i])

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
bt.plotfit(stepvalues/20.35,
		variance,
		out.beta,
		out.X_unweighted,
		top='Data is real!',
		figlabel='Comparison',
		error=used_error)

# figchisquare = plt.figure()
# mt.plot_featured(stepvalues,chisq_red,'.-',
#                 toplabel='Chi-Squared for Each Gaussian Fit',
#                 xlabel='$E/E_0$',
#                 ylabel='$\chi^2$')

plt.show()
