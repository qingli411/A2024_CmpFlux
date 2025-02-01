#!/usr/bin/env python
# coding: utf-8

import sys
import os
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
import pandas as pd
sys.path.append(os.path.join(os.path.pardir, 'lesview'))
from lesview import *

casename = 'lsc_ymc22_sbl_bbl_v2'
datapath = os.path.join(os.path.pardir, 'oceananigans', '{:s}'.format(casename))
figpath  = 'overview_{:s}'.format(casename)
os.makedirs(figpath, exist_ok=True)

filepath = os.path.join(datapath, 'averages.jld2')
data_pfl = OceananigansDataProfile(filepath=filepath)

g = 9.81
H = 30
u10 = 8
N2 = 1.962e-4
bstar = N2 * H
cd = 1.25e-3
rhoa = 1.225
rhoo = 1026
tau = rhoa/rhoo*cd*u10*u10
ustar = np.sqrt(tau)
amplitude = 1.0
wavelength = 60
wavenumber = 2.*np.pi/wavelength
frequency = np.sqrt(g*wavenumber*np.tanh(wavenumber*H))
us0 = amplitude**2*wavenumber
la = np.sqrt(ustar/us0)
print(la)

tavg1 = dict(starttime='2000-01-04T00:00:00', endtime='2000-01-04T17:00:00', line_kw=dict(color='k', linestyle='--'))
tavg2 = dict(starttime='2000-01-10T00:00:00', endtime='2000-01-10T17:00:00', line_kw=dict(color='k', linestyle='-'))
tavgs = dict(T1=tavg1, T2=tavg2)

def plot_overview(das, levels, labels, tavgs):

    nv = len(das)
    fig, axarr = plt.subplots(nv, 2, gridspec_kw={'width_ratios': [1, 5]})
    fig.set_size_inches([8, 0.4+2*nv])
    rlcolor = {'RdBu_r': 'k', 'viridis': 'w'}
    date_form = DateFormatter("%d")
    for i, var in enumerate(das.keys()):
        ax = np.ravel(axarr)[i*2+1]
        cf = das[var].plot(ax=ax, levels=levels[var], cbar_kwargs={'label': labels[var]})
        cmap = cf.get_cmap().name
        for j, tag in enumerate(tavgs.keys()):
            ax.axvline(x=pd.Timestamp(tavgs[tag]['starttime']), linestyle=':', color=rlcolor[cmap])
            ax.axvline(x=pd.Timestamp(tavgs[tag]['endtime']), linestyle=':', color=rlcolor[cmap])
            ax.text(pd.Timestamp(tavgs[tag]['starttime']), 0, tag, va='bottom', ha='left')
        ax.set_title('')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.xaxis.set_major_formatter(date_form)
        for lb in ax.get_xticklabels(which='major'):
            lb.set(rotation=0, horizontalalignment='right')
        ax = np.ravel(axarr)[i*2+0]
        for j, tag in enumerate(tavgs.keys()):
            tslice = slice(tavgs[tag]['starttime'], tavgs[tag]['endtime'])
            das[var].sel(time=tslice).mean(dim='time').plot(ax=ax, y=das[var].dims[0], label=tag, **tavgs[tag]['line_kw'])
        ax.set_xlabel(labels[var])
        ax.set_ylabel('Depth [m]')
    axarr[0,0].legend()
    axarr[-1,-1].set_xlabel('Time [day]')

    plt.tight_layout()
    return fig


# plot buoyancy


das = dict(
    b  = ds.data_vars['b']/bstar,
    NN = ds.data_vars['b'].differentiate(coord='z')/N2,
    wb = ds.data_vars['wb']/ustar/bstar*1e3,
    bb = ds.data_vars['bb']/bstar**2*1e3,
)
labels = dict(
    b  = '$b/b_*$',
    NN = '$N^2/N_0^2$',
    wb = '$10^3${:s}$/u_*b_*$'.format(ds.data_vars['wb'].long_name),
    bb = '$10^3${:s}$/b_*^2$'.format(ds.data_vars['bb'].long_name),
)
levels = dict(
    b  = np.linspace(-1, 1, 41),
    NN = np.linspace(0, 4, 41),
    wb = np.linspace(-4, 4, 41),
    bb = np.linspace(0, 4),
)
fig = plot_overview(das, levels, labels, tavgs)
figname = os.path.join(figpath, 'buoyancy')
fig.savefig(figname, dpi = 300, facecolor='w')


# plot velocity


das = dict(
    u = ds.data_vars['u']/ustar,
    v = ds.data_vars['v']/ustar,
    dudz = ds.data_vars['u'].differentiate(coord='z')*H/ustar*0.01,
    dvdz = ds.data_vars['v'].differentiate(coord='z')*H/ustar*0.01,
)
labels = dict(
    u = '$u/u_*$',
    v = '$v/u_*$',
    dudz = '$0.01\partial_z u H/u_*$',
    dvdz = '$0.01\partial_z v H/u_*$',
)
levels = dict(
    u = np.linspace(-40, 40, 41),
    v = np.linspace(-15, 15, 31),
    dudz = np.linspace(-2, 2, 41),
    dvdz = np.linspace(-2, 2, 41),
)
fig = plot_overview(das, levels, labels, tavgs)
figname = os.path.join(figpath, 'velocity')
fig.savefig(figname, dpi = 300, facecolor='w')


# plot flux


das = dict(
    wb = ds.data_vars['wb']/ustar/bstar*1e3,
    wu = ds.data_vars['wu']/ustar**2,
    wv = ds.data_vars['wv']/ustar**2,
)
labels = dict(
    wb = '$10^3${:s}$/u_*b_*$'.format(ds.data_vars['wb'].long_name),
    wu = '{:s}$/u_*^2$'.format(ds.data_vars['wu'].long_name),
    wv = '{:s}$/u_*^2$'.format(ds.data_vars['wv'].long_name),
)
levels = dict(
    wb = np.linspace(-4, 4, 41),
    wu = np.linspace(-2, 2, 41),
    wv = np.linspace(-1.2, 1.2, 41),
)
fig = plot_overview(das, levels, labels, tavgs)
figname = os.path.join(figpath, 'flux')
fig.savefig(figname, dpi = 300, facecolor='w')


# plot TKE


das = dict(
    tke = 0.5*(ds.data_vars['uu']+ds.data_vars['vv']+ds.data_vars['ww'].interp(zi=ds.z))/ustar**2,
    uu  = ds.data_vars['uu']/ustar**2,
    vv  = ds.data_vars['vv']/ustar**2,
    ww  = ds.data_vars['ww']/ustar**2
)
labels = dict(
    tke = '$TKE/u_*^2$',
    uu  = '{:s}$/u_*^2$'.format(ds.data_vars['uu'].long_name),
    vv  = '{:s}$/u_*^2$'.format(ds.data_vars['vv'].long_name),
    ww  = '{:s}$/u_*^2$'.format(ds.data_vars['ww'].long_name),
)
levels = dict(
    tke = np.linspace(0, 8, 41),
    uu  = np.linspace(0, 8, 41),
    vv  = np.linspace(0, 8, 41),
    ww  = np.linspace(0, 8, 41),
)
fig = plot_overview(das, levels, labels, tavgs)
figname = os.path.join(figpath, 'tke')
fig.savefig(figname, dpi = 300, facecolor='w')


# plot skewness of vertical velocity


das = dict(
    ww = ds.data_vars['ww']/ustar**2,
    w3 = ds.data_vars['w3']/ustar**3,
    sk = ds.data_vars['w3']/ds.data_vars['ww']**1.5,
)
labels = dict(
    ww = '{:s}$/u_*^2$'.format(ds.data_vars['ww'].long_name),
    w3 = '{:s}$/u_*^3$'.format(ds.data_vars['w3'].long_name),
    sk = '{:s}$/(${:s}$)^{:s}$'.format(ds.data_vars['w3'].long_name, ds.data_vars['ww'].long_name,'{3/2}'),
)
levels = dict(
    ww = np.linspace(0, 4, 41),
    w3 = np.linspace(-4, 4, 41),
    sk = np.linspace(-1, 1, 41),
)
fig = plot_overview(das, levels, labels, tavgs)
figname = os.path.join(figpath, 'skewness')
fig.savefig(figname, dpi = 300, facecolor='w')


# plot N2 and S2


das = dict(
    NN = ds.data_vars['b'].differentiate(coord='z')/N2,
    SS = np.sqrt(ds.data_vars['u'].differentiate(coord='z')**2+ds.data_vars['v'].differentiate(coord='z')**2)*H**2/ustar**2*1e-5,
    Ri = ds.data_vars['b'].differentiate(coord='z')/np.sqrt(ds.data_vars['u'].differentiate(coord='z')**2+ds.data_vars['v'].differentiate(coord='z')**2)
)
labels = dict(
    NN = '$N^2/N_0^2$',
    SS = '$10^{-5}S^2 H^2/u_*^2$',
    Ri = '$N^2/S^2$',
)
levels = dict(
    NN = np.linspace(0, 4, 41),
    SS = np.linspace(0, 4, 41),
    Ri = np.linspace(0, 0.25, 26),
)
fig = plot_overview(das, levels, labels, tavgs)
figname = os.path.join(figpath, 'mixing')
fig.savefig(figname, dpi = 300, facecolor='w')


# plot tracer Cs


das = dict(
    Cs    = np.log10(ds.data_vars['Cs']),
    dCsdz = ds.data_vars['Cs'].differentiate(coord='z')*1e3,
    # wCs   = (ds.data_vars['wCs']+ds.data_vars['wCssb'])*1e7,
    wCs   = ds.data_vars['wCs']*1e7,

)
labels = dict(
    Cs    = '$\log_{10}(C_s)$',
    dCsdz = '$10^3\partial_z C_s$',
    wCs   = '$10^7\overline{w^\prime C_s^\prime}$',
)
levels = dict(
    Cs    = np.linspace(-2.5, -1.5, 41),
    dCsdz = np.linspace(-8, 8, 41),
    wCs   = np.linspace(-8, 8, 41),
)
fig = plot_overview(das, levels, labels, tavgs)
figname = os.path.join(figpath, 'tracer_Cs')
fig.savefig(figname, dpi = 300, facecolor='w')


# plot tracer Cb


das = dict(
    Cb    = np.log10(ds.data_vars['Cb']),
    dCbdz = ds.data_vars['Cb'].differentiate(coord='z')*1e3,
    # wCb   = (ds.data_vars['wCb']+ds.data_vars['wCbsb'])*1e7,
    wCb   = ds.data_vars['wCb']*1e7,

)
labels = dict(
    Cb    = '$\log_{10}(C_b)$',
    dCbdz = '$10^3\partial_z C_b$',
    wCb   = '$10^7\overline{w^\prime C_b^\prime}$',
)
levels = dict(
    Cb    = np.linspace(-2.5, -1.5, 41),
    dCbdz = np.linspace(-8, 8, 41),
    wCb   = np.linspace(-8, 8, 41),
)
fig = plot_overview(das, levels, labels, tavgs)
figname = os.path.join(figpath, 'tracer_Cb')
fig.savefig(figname, dpi = 300, facecolor='w')

