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

def get_tslice(N2, T_io, T_bm=1.5*86400, T_am=3.5*86400):
    z0, z1 = get_edges(N2)
    time_merge = z0.dropna(dim='time').time[-1]
    print('Time of merge:')
    print(pd.to_datetime(time_merge.data))
    tstart1, tend1 = time_merge - pd.Timedelta(T_io, 's') - pd.Timedelta(T_bm, 's'), time_merge - pd.Timedelta(T_bm, 's')
    tstart2, tend2 = time_merge + pd.Timedelta(T_am, 's'), time_merge + pd.Timedelta(T_am, 's') + pd.Timedelta(T_io, 's')
    print('T1:')
    print(pd.to_datetime(tstart1.data))
    print(pd.to_datetime(tend1.data))
    print('T2:')
    print(pd.to_datetime(tstart2.data))
    print(pd.to_datetime(tend2.data))
    tslice1 = slice(tstart1, tend1)
    tslice2 = slice(tstart2, tend2)
    return tslice1, tslice2