import numpy as np
import xarray as xr
import pandas as pd

def get_edges(NN):
    z = NN.coords['z'].values
    nt = NN.time.shape[0]
    NNc = 0.2
    zup = np.zeros(nt)
    zdn = np.zeros(nt)
    for i in np.arange(nt):
        idxlst = np.where(NN[:,i]>NNc)
        if idxlst[0].size > 0:
            zup[i] = z[np.max(idxlst)]
            zdn[i] = z[np.min(idxlst)]
        else:
            zup[i] = np.nan
            zdn[i] = np.nan
    z0 = xr.DataArray(zup, dims=['time'], coords={'time': NN.time},
                      attrs={'long_name': 'upper edge', 'units': 'm'})
    z1 = xr.DataArray(zdn, dims=['time'], coords={'time': NN.time},
                      attrs={'long_name': 'lower edge', 'units': 'm'})
    return (z0, z1)


def get_flux(lam, num, gam):
    tmp = xr.zeros_like(num)
    nt = lam.shape[1]
    for i in np.arange(nt):
        tmp.data[1:-1,i] = (lam.data[:-1,i]-lam.data[1:,i])/(lam.z.data[:-1]-lam.z.data[1:])
    flux = - num * tmp + gam
    return flux

def inertial_period(lat=45):
    Omega = 2*np.pi/86400
    f = 2*Omega*np.sin(np.deg2rad(lat))
    return 2*np.pi/f

def nondim_da(da, H, Tf):
    zdim = da.dims[0]
    tdim = da.dims[1]
    z_norm = da.coords[zdim].data/H
    time_dtime = pd.to_datetime(da.coords[tdim].data)
    time_sec = (time_dtime-time_dtime[0]).total_seconds()
    time_norm = time_sec/Tf
    da_new = da.assign_coords({zdim: z_norm, tdim: time_norm})
    da_new.coords[zdim].attrs['long_name'] = '$z/H$'
    da_new.coords[zdim].attrs['units'] = ''
    da_new.coords[tdim].attrs['long_name'] = '$t/T_f$'
    da_new.coords[tdim].attrs['units'] = ''
    return da_new

def nondim_das(das, **kwargs):
    das_new = {}
    for var in das.keys():
        das_new[var] = nondim_da(das[var], **kwargs)
    return das_new

def get_merge(NN, Tf=inertial_period(lat=45), nTf0=1, nTf1=4):
    z0, z1 = get_edges(NN)
    time_merge = z0.dropna(dim='time').time[-1]
    if hasattr(time_merge, 'long_name') and time_merge.attrs['long_name'] == '$t/T_f$':
        tstart1, tend1 = time_merge-1-nTf0, time_merge-nTf0
        tstart2, tend2 = time_merge+nTf1, time_merge+nTf1+1
        print('Time of merge:')
        print('t/Tf = {:g}'.format(time_merge.data))
        print('T1:')
        print('t/Tf = {:g} -- {:g}'.format(tstart1.data, tend1.data))
        print('T2:')
        print('t/Tf = {:g} -- {:g}'.format(tstart2.data, tend2.data))
    else:
        tstart1, tend1 = time_merge - pd.Timedelta(Tf, 's') - pd.Timedelta(nTf0*Tf, 's'), time_merge - pd.Timedelta(nTf0*Tf, 's')
        tstart2, tend2 = time_merge + pd.Timedelta(nTf1*Tf, 's'), time_merge + pd.Timedelta(nTf1*Tf, 's') + pd.Timedelta(Tf, 's')
        print('Time of merge:')
        print(pd.to_datetime(time_merge.data))
        print('T1:')
        print(pd.to_datetime(tstart1.data))
        print(pd.to_datetime(tend1.data))
        print('T2:')
        print(pd.to_datetime(tstart2.data))
        print(pd.to_datetime(tend2.data))
    return z0, z1, tstart1.data, tend1.data, tstart2.data, tend2.data

def get_tslice(NN, Tf, **kwargs):
    _, _, tstart1, tend1, tstart2, tend2 = get_merge(NN, Tf=Tf, **kwargs)
    tslice1 = slice(tstart1, tend1)
    tslice2 = slice(tstart2, tend2)
    return tslice1, tslice2    