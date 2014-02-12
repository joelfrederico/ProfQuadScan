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

# Load and transfer matlab variables
infile = 'tempfiles/forpython.mat'
f=h5.File(infile);
data=f['data']
imgs=data['processed']['images']['CEGAIN']['dat']
stepvalues = data['raw']['scalars']['step_value']['dat']
stepvalues=mt.derefdataset(stepvalues,f)
stepvalues = np.unique(stepvalues)

xstart=455
xstop=457
ystart=550
ystop=700
step=2
xrange = slice(xstart,xstop)
yrange = slice(ystart,ystop)

# These need to be fixed someday...
qs1_k_half = 3.077225846087095e-01;
qs2_k_half = -2.337527121004531e-01;

# Check number of points and initialize arrays
num_pts=(xstop-xstart)/step
print num_pts
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
chisq_red = np.zeros(num_pts)

# Find spot size for each step
mt.figure('Shot')
for i,img in enumerate(imgs):
	img=np.flipud(np.rot90(f[img[0]]))
	plt.imshow(img)
	# raw_input('Wait...')

	# Fit individual slices
	for val in mt.linspacestep(xstart,xstop-step,step):
		# Take a strip of the image
		# strip = img[slice(475,480),yrange]
		strip = img[slice(val,val+step),yrange]
		# fig=mt.figure('Strip')
		# plt.imshow(strip)
		
		# Sum over the strip to get an average of sorts
		histdata = sum(strip)
		xbins = len(histdata)
		x = np.linspace(1,xbins,xbins)*40.3
		# plt.plot(x,histdata)
		
		# Fit with a Gaussian to find spot size
		# fig=mt.figure('Gaussian Fit')
		# plotbool=True
		plotbool = False
		# varbool  = False
		varbool  = True
		popt,pcov,chisq_red = mt.gaussfit(x,histdata,sigma_y=np.ones(xbins),plot=plotbool,variance_bool=varbool,verbose=False,p0=[16000,3500,500,0])
		# raw_input('Press enter?')
		# plt.close(fig)
	
		variance[i] = popt[2]
		bact = setQS.set_QS_energy_ELANEX(stepvalues[i])
		LGPS_3261[i] = setQS.bdes2K1(bact[0],20.35)
		LGPS_3311[i] = setQS.bdes2K1(bact[1],20.35)

		print 'QS1 K1: {}\tQS2 K1: {}'.format(LGPS_3261[i],LGPS_3311[i])


mt.figure('Variance')
plt.plot(stepvalues,variance,'.-')

# Set up initial conditions

# RMS
# emitx = 0.00201/gamma
# betax = 11.2988573693
# alphax = 6.72697997971

# Gauss fit
B5D36_en  = 20.35
gamma     = (B5D36_en/0.5109989)*1e3

emitx  = 0.000100/gamma
betax  = .5
alphax = -1
gammax = (1+np.power(alphax,2))/betax
twiss  = np.array([betax,alphax,gammax])
T      = np.array([[betax,-alphax],[-alphax,gammax]])

# Create Beamline {{{
IP2QS1 = sltr.Drift(length = 5.4217)
QS1 = sltr.Quad(length= 5.000000000E-01,K1= qs1_k_half)
# QS1._change_E(gamma,new_gamma)
LQS12QS2 = sltr.Drift(length = 4.00E+00)
QS2 = sltr.Quad(length= 5.000000000E-01,K1=qs2_k_half)
# QS2._change_E(gamma,new_gamma)
LQS22BEND = sltr.Drift(length = 0.7428E+00)
B5D36 = sltr.Bend(
		length= 2*4.889500000E-01,          
               angle= 6.0E-03, 	     
	       order=1,
	       rotate=0
	       )
LBEND2ELANEX = sltr.Drift(length = 8.792573)
beamline = sltr.Beamline(
		element_list=[
			IP2QS1        ,
			QS1           ,
			QS1           ,
			LQS12QS2      ,
			QS2           ,
			QS2           ,
			LQS22BEND     ,
			B5D36         ,
			LBEND2ELANEX	
			],
		gamma= gamma
		)

beamline_array = np.array([])

for i,beam in enumerate(stepvalues):
	# Set qs1
	beamline.elements[1]._K1 = LGPS_3261[i]
	beamline.elements[2]._K1 = LGPS_3261[i]
	beamline.elements[4]._K1 = LGPS_3311[i]
	beamline.elements[5]._K1 = LGPS_3311[i]
	print beamline
	beamline_array = np.append(beamline_array,copy.deepcopy(beamline))

# Fit bowtie plot
chisq_factor = 1
# chisq_factor = 63.6632188
used_error   = stddev*np.sqrt(chisq_factor)

out          = bt.fitbowtie(beamline,davg,variance,T,twiss,emitx,error=used_error, verbose=True)
spotexpected = out.spotexpected
X            = out.X
beta         = out.beta
covar        = out.covar
# print covar

figcher=plt.figure()
top='Simulated Energy Emittance Measurement\nNOT PHYSICAL'
bt.plotfit(davg,variance,beta,out.X_unweighted,spotexpected,top,error=used_error)

figchisquare = plt.figure()
plt.plot(davg,chisq_red)
mt.plot_featured(davg,chisq_red,'.-',
		toplabel='Chi-Squared for Each Gaussian Fit',
		xlabel='$E/E_0$',
		ylabel='$\chi^2$')
# print davg
# print chisq_red


# if __name__ == '__main__':

#         parser=argparse.ArgumentParser(description=
#                         'Wrap python analysis to be called at the command line.')
#         parser.add_argument('-V',action='version',version='%(prog)s v0.2')
#         parser.add_argument('-v','--verbose',action='store_true',
#                         help='Verbose mode.')
#         parser.add_argument('-o','--output',action='store',
#                         help='Output filename. (Default: no file output.)')
#         parser.add_argument('inputfile',
#                         help='Input Matlab v7 file.')
#         parser.add_argument('-f','--fit', choices=['gauss', 'bigauss'], default='gauss', 
#                         help='Type of fit to spot size profile. (Default: %(default)s)')
#         arg=parser.parse_args()

#         out=wrap_analyze(arg.inputfile)
#         
#         if arg.verbose:
#                 plt.show()

