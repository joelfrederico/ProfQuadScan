import numpy as _np

def set_QS_energy_ELANEX(E):
	    VAL = (1+E/20.35)*_np.array((210.87, -164.95)); # Imaging condition for object = E200 IP or plasma exit and image = ELANEX in Fall of 2013
	    # control_magnetSet({'LGPS:LI20:3261', 'LGPS:LI20:3311'}, VAL,  'action', 'TRIM');
	    return VAL

def bdes2K1(bdes,E):
	Brho=E/0.029979
	K1 = bdes/(Brho)
	return K1
