import numpy as np
import xarray as xr

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