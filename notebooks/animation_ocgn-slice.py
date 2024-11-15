#!/usr/bin/env python
# coding: utf-8

import sys
import os
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(os.path.join(os.environ['HOME'], 'local', 'lesview'))
from lesview import *
from lesview.plots import plot_box_slice

casename = 'lsc_ymc22_bbl_N2p1_viz'
datapath = os.path.join(os.path.pardir, 'tests', '{:s}'.format(casename))

data_xys = OceananigansDataVolume(filepath=os.path.join(datapath, 'slices_xy.jld2'))
data_xyb = OceananigansDataVolume(filepath=os.path.join(datapath, 'slices_xy2.jld2'))
data_xz  = OceananigansDataVolume(filepath=os.path.join(datapath, 'slices_xz.jld2'))
data_yz  = OceananigansDataVolume(filepath=os.path.join(datapath, 'slices_yz2.jld2'))

var = 'w'
view = 'both'
levels = np.linspace(-0.04, 0.04, 41)
if var == 'u':
    U0 = 0.25
    if rf in casename:
        U0 *= -1.0
    levels += U0
mp4dir = 'animation_{:s}-{:s}'.format(casename, view)
os.makedirs(mp4dir, exist_ok=True)

da_xys = data_xys.dataset.data_vars[var][0,:,:,:]
da_xyb = data_xyb.dataset.data_vars[var][0,:,:,:]
da_xz  = data_xz.dataset.data_vars[var][:,0,:,:]
da_yz  = data_yz.dataset.data_vars[var][:,:,0,:]

for itime in np.arange(da_xz.time.size):
    fig = plot_box_slice(da_xz, da_yz, da_xys, -4, da_xy2=da_xyb, iz2=3, view=view, itime=itime, cmap='RdBu_r', levels=levels)
    figname = os.path.join(mp4dir, '{:s}_{:06d}'.format(var, itime))
    plt.tight_layout()
    fig.savefig(figname, dpi = 300, facecolor='w')
    fig.clear()
    plt.close(fig)

