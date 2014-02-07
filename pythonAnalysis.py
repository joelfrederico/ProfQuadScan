#!/usr/bin/env python
import argparse
import numpy as np
import scipy.io as sio
import ButterflyEmittancePython as bt
import mytools.slactrac as sltr
import matplotlib.pyplot as plt
import mytools as mt
import copy
plt.close('all')

# Load and transfer matlab variables
infile = 'tempfiles/forpython.mat'
matvars        = sio.loadmat(infile)
transfer = matvars['transfer']
res = transfer['res']
res = res[0,0][0,0]
# print '{}'.format(res)
img = matvars['img']
img = img[0]
xstart=420
xstop=500
ystart=550
ystop=700
step=2
xrange = slice(xstart,xstop)
yrange = slice(ystart,ystop)
# fig=plt.figure()
fig=mt.figure('To process')
a=img[5]
plt.imshow(a[xrange,yrange])

# These need to be fixed someday...
qs1_k_half = 3.077225846087095e-01;
qs2_k_half = -2.337527121004531e-01;

# addpath('E200_Cam_Energy_Cal');

# Check number of points and initialize arrays
num_pts=(xstop-xstart)/step
sigs=np.zeros(num_pts)
# num_pts = len(sum_x)
variance  = np.zeros(num_pts)
stddev    = np.zeros(num_pts)
varerr    = np.zeros(num_pts)
# print varerr.shape
chisq_red = np.zeros(num_pts)

# Fit individual slices
for i,val in enumerate(mt.linspacestep(xstart,xstop-step,step)):
	print i
	# Take a strip of the image
	# strip = img[slice(475,480),yrange]
	strip = a[slice(val,val+step),yrange]
	# fig=mt.figure('Strip')
	# plt.imshow(strip)
	
	# Sum over the strip to get an average of sorts
	histdata = sum(strip)
	xbins = len(histdata)
	x = np.linspace(1,xbins,xbins)*res
	# plt.plot(x,histdata)
	
	# Fit with a Gaussian to find spot size
	# fig=mt.figure('Gaussian Fit')
	# plotbool=True
	plotbool = False
	varbool  = False
	popt,pcov,chisq_red = mt.gaussfit(x,histdata,sigma_y=np.ones(xbins),plot=plotbool,variance_bool=varbool,verbose=False,p0=[16000,3500,500,0])
	# raw_input('Press enter?')
	# plt.close(fig)

	variance[i] = popt[2]
	# varerr[i]   = pcov[2,2]
	# stddev[i]   = np.sqrt(pcov[2,2])

mt.figure('Varplot')
# plt.plot(variance)

xvar=np.shape(mt.linspacestep(xstart,xstop,step))[0]-1
xvar=mt.linspacestep(1,xvar)

out=np.polyfit(xvar,variance,2)

plt.plot(xvar,variance,'.-',xvar,np.polyval(out,xvar))

# Central energy is at 18.

print 'Minimum at {}'.format(-out[1]/(2*out[0]))

	# for i,el in enumerate(sum_x):
	#         y                   = sum_x[i,:]
	#         erry = np.sqrt(y)
	#         erry[erry==0] = 0.3
	#         # plt.plot(y)
	#         # plt.show()
	#         # popt,pcov,chisq_red[i] = mt.gaussfit(x_meter,y,sigma_y=erry,plot=True,variance_bool=True,verbose=False)
	#         popt,pcov,chisq_red[i] = mt.gaussfit(x_meter,y,sigma_y=np.ones(len(y)),plot=True,variance_bool=True,verbose=False)
	#         plt.show()
	#         variance[i]         = popt[2]
	#         varerr[i]           = pcov[2,2]
	#         stddev[i]           = np.sqrt(pcov[2,2])

# Set up initial conditions

# RMS
# emitx = 0.00201/gamma
# betax = 11.2988573693
# alphax = 6.72697997971

# Gauss fit
B5D36_en  = 20.35
gamma     = (B5D36_en/0.5109989)*1e3

emitx = 0.000100/gamma
betax = .5
alphax = -1 
gammax = (1+np.power(alphax,2))/betax
twiss=np.array([betax,alphax,gammax])
T = np.array([[betax,-alphax],[-alphax,gammax]])

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
beamline.calc_mat()
# beamline.change_energy(new_gamma);
# for el in beamline.elements:
#         print '---------------'
#         print 'Energy of {} = {}'.format(el._type,el._gamma)
#         try:
#                 print 'K1 = {}'.format(el._K1)
#         except:
#                 pass
#         try:
#                 print 'K2 = {}'.format(el._K2)
#         except:
#                 pass
#         try:
#                 print 'Length = {}'.format(el._length)
#         except:
#                 pass
# print '---------------'
#}}}

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

