#!/usr/bin/env python
import numpy as np
import pyfits
import sys

######################################
## Conver map.bin to map.fits
######################################

if (len(sys.argv) > 1): 
    dir_map = sys.argv[1]+'/'
else:
    dir_map = './'

f = open(dir_map+"mapmeta.dat","r") 
line0 = f.readline().split(' ')
line1 = f.readline().split(' ')
line2 = f.readline().split(' ')
line3 = f.readline().split(' ')

mu_av    = float(line0[0])
N_av     = float(line0[1])
res      =   int(line1[0])
map_size = float(line2[0])
ka       = float(line3[0])
ga       = float(line3[1])
#s_pa     = float(line3[2])


mu_th = 1./((1-ka)**2-ga**2)

f   = open(dir_map+"map.bin","rb")
mu = np.fromfile(f,'i',-1,"")
mu = mu.reshape(res,res)
print 'res:',res,
print '# of Nij=0:',np.count_nonzero(mu==0)
N_rays = N_av/mu_av
mu = mu*mu_av/N_av
print 'N_rays=',N_rays,'N_av=',N_av,'mu_av=',mu_av,'mu_th=',mu_th

pyfits.writeto(dir_map+'map.fits', mu, clobber=True)


