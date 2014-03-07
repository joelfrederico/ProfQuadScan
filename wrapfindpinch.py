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
plt.close('all')


infile = 'tempfiles/forpython.mat'
f=h5.File(infile);
data=f['data']
imgs=data['processed']['images']['CEGAIN']['dat']
img6=f[imgs[5,0]]
img6=np.flipud(np.rot90(img6))
plt.imshow(img6)

ystart=420
ystop=500
xstart=550
xstop=700
step=2

mt.findpinch(img6,[xstart,xstop],[ystart,ystop],step=2)

plt.show()
