#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.dates import DateFormatter
import numpy as np
sys.path.append(os.path.join(os.environ['HOME'], 'local', 'lesview'))
from lesview import *
from lesview.plots import plot_box_field


# In[2]:


casename = 'lsc_ymc22_sbl_bbl_v2'
datapath = os.path.join(os.path.pardir, 'tests', '{:s}'.format(casename))


# In[3]:


filepath = os.path.join(datapath, 'fields.jld2')
data_fld = OceananigansDataVolume(filepath=filepath, fieldname=['b','w'])


# In[4]:


filepath = os.path.join(datapath, 'averages.jld2')
data_pfl = OceananigansDataProfile(filepath=filepath)


# In[5]:


ds = data_fld.dataset
ds


# In[6]:


dsp = data_pfl.dataset
dsp


# In[7]:


# g = 9.81
# N2 = 1.962e-4
# T0 = 20.0
# alphaT = 2.0e-4
# dTdz = N2/alphaT/g
# H = 30

# var = 'b'
# levels = np.linspace(0,6e-3,41)
# levels = np.linspace(17,20,31)
# da = ds.data_vars[var]/alphaT/g + (T0 - dTdz*H)
# da.attrs['long_name'] = '$T$'
# da.attrs['units'] = '$^\circ$C'
# # levels = None

# dap = dsp.data_vars['b']/alphaT/g + (T0 - dTdz*H)
# dap.attrs['long_name'] = '$T$'
# dap.attrs['units'] = '$^\circ$C'


# In[8]:


# date_form = DateFormatter("%d")

# itime = 10
# fig = plot_box_field(da, view='both', itime=itime, cmap='viridis', levels=levels)
# fig.set_size_inches(8,8)
# ax = plt.gca()
# gs = gridspec.GridSpec(5, 8)
# ax.set_subplotspec(gs[0:4, 0:-1])
# ax2 = fig.add_subplot(gs[4, 1:-1])
# dap.plot(ax=ax2, levels=levels)
# ax2.xaxis.set_major_formatter(date_form)
# for lb in ax2.get_xticklabels(which='major'):
#     lb.set(rotation=0, horizontalalignment='right')
# ax2.set_xlabel('Time [day]')
# ax2.axvline(x=da.time[itime].data, linestyle=':', color='w')
# plt.tight_layout()
# fig.savefig('test', dpi = 300, facecolor='w')


# In[9]:


# var = 'Cs'
# levels = np.linspace(0,0.14,281)
# # levels = None
# fig = plot_box_field(ds.data_vars[var], view='both', itime=-1, cmap='prism_r', levels=levels)


# In[10]:


# var = 'Cb'
# levels = np.linspace(0,0.14,281)
# # levels = None
# fig = plot_box_field(ds.data_vars[var], view='both', itime=-1, cmap='prism_r', levels=levels)


# In[11]:


# var = 'w'
# levels = np.linspace(-0.04, 0.04, 41)
# fig = plot_box_field(ds.data_vars[var], view='both', itime=-1, cmap='RdBu_r', levels=levels)


# In[12]:


# var = 'u'
# levels = np.linspace(-0.4, 0.4, 41)
# fig = plot_box_field(ds.data_vars[var], view='both', itime=-1, cmap='RdBu_r', levels=levels)


# In[13]:


# var = 'v'
# levels = np.linspace(-0.2, 0.2, 41)
# fig = plot_box_field(ds.data_vars[var], view='both', itime=-1, cmap='RdBu_r', levels=levels)


# In[14]:


# var = 'T'
# view = 'both'
# levels = np.linspace(17,20,61)
# mp4dir = './{:s}-{:s}'.format(casename, view)
# os.makedirs(mp4dir, exist_ok=True)

# g = 9.81
# N2 = 1.962e-4
# T0 = 20.0
# alphaT = 2.0e-4
# dTdz = N2/alphaT/g
# H = 30
# da = ds.data_vars['b']/alphaT/g + (T0 - dTdz*H)
# da.attrs['long_name'] = '$T$'
# da.attrs['units'] = '$^\circ$C'

# for itime in np.arange(da.time.size):
#     fig = plot_box_field(da, view=view, itime=itime, cmap='viridis', levels=levels)
#     figname = os.path.join(mp4dir, '{:s}_{:06d}'.format(var, itime))
#     plt.tight_layout()
#     fig.savefig(figname, dpi = 300, facecolor='w')
#     fig.clear()
#     plt.close(fig)


# In[15]:


var = 'T'
view = 'both'
date_form = DateFormatter("%d")
levels = np.linspace(17,20,61)
mp4dir = './{:s}-{:s}'.format(casename, view)
os.makedirs(mp4dir, exist_ok=True)

g = 9.81
N2 = 1.962e-4
T0 = 20.0
alphaT = 2.0e-4
dTdz = N2/alphaT/g
H = 30
da = ds.data_vars['b']/alphaT/g + (T0 - dTdz*H)
da.attrs['long_name'] = '$T$'
da.attrs['units'] = '$^\circ$C'

dap = dsp.data_vars['b']/alphaT/g + (T0 - dTdz*H)
dap.attrs['long_name'] = '$T$'
dap.attrs['units'] = '$^\circ$C'

for itime in np.arange(da.time.size):
    fig = plot_box_field(da, view=view, itime=itime, cmap='viridis', levels=levels)
    fig.set_size_inches(8,8)
    ax = plt.gca()
    gs = gridspec.GridSpec(5, 8)
    ax.set_subplotspec(gs[0:4, 0:-1])
    ax2 = fig.add_subplot(gs[4, 1:-1])
    dap.plot(ax=ax2, levels=levels)
    ax2.xaxis.set_major_formatter(date_form)
    for lb in ax2.get_xticklabels(which='major'):
        lb.set(rotation=0, horizontalalignment='right')
    ax2.set_xlabel('Time [day]')
    ax2.axvline(x=da.time[itime].data, linestyle=':', color='w')
    figname = os.path.join(mp4dir, '{:s}_v2_{:06d}'.format(var, itime))
    plt.tight_layout()
    fig.savefig(figname, dpi = 300, facecolor='w')
    fig.clear()
    plt.close(fig)


# In[17]:


# var = 'w'
# view = 'both'
# levels = np.linspace(-0.04, 0.04, 41)
# mp4dir = './{:s}-{:s}'.format(casename, view)
# os.makedirs(mp4dir, exist_ok=True)

# for itime in np.arange(ds.time.size):
#     fig = plot_box_field(ds.data_vars[var], view=view, itime=itime, cmap='RdBu_r', levels=levels)
#     figname = os.path.join(mp4dir, '{:s}_{:06d}'.format(var, itime))
#     plt.tight_layout()
#     fig.savefig(figname, dpi = 300, facecolor='w')
#     fig.clear()
#     plt.close(fig)

