import logging
logging.basicConfig(level=logging.INFO)
import CtdHexTrans as ctdtrans
from datetime import datetime, timedelta
import scipy.io as sio
import xarray as xr
import numpy as np
import seawater.eos80 as sw
import seawater
import glob as glob
import os.path
import subprocess
from jmkdata import bindata1d
from timeit import default_timer as timer
import getInletX
import sys

todo = sys.argv[1]

coffset = 1.7

innames = glob.glob(todo+'/*.hex')
cgrid = dict()

for name in innames:
    basename = os.path.splitext(os.path.basename(name))[0]
    dirname = os.path.dirname(name)
    logging.info(dirname)
    ncname = dirname+'/'+basename+'.nc'
    matname = dirname+'/'+basename+'New.mat'
    if not os.path.exists(ncname) or (os.path.getctime(name) >
            os.path.getctime(ncname)):
        logging.warning('Hi!')

        ctd = ctdtrans.ctd_hex_trans(name)
        ctd = ctdtrans.raw_to_science(ctd, coffset)

        alongx, acrossx = getInletX.getInletX(ctd['lon'], ctd['lat'])

        ctd['alongx'] = alongx
        ctd['acrossx'] = acrossx

        ctdtrans.tonetcdf(ncname, ctd, coffset)
        ctdtrans.tomatlab(matname, ctd)
        subprocess.call(['octave', 'saveasstruct.m', matname, 'ctd'])

        logging.info(ctd['lon'])
        logging.info(ctd['lat'])
        logging.info(ctd['time'])
        logging.info(ctd['matlabtime'])

        logging.debug(cgrid.keys())

# make cgrid

innames = glob.glob(todo+'/20*.nc')
zbins = np.arange(325)
cgrid = dict()
cgrid['depths'] = zbins[1:] - 0.5
cgrid['name'] = []
for name in innames:
    ds = xr.open_dataset(name)
    logging.debug(ds)
    logging.info('Starting %s', name)
    # find downcast...
    ind = np.where(ds.pres > 20)[0][0]
    start = ind
    stop = ind
    p = ds.pres.values
    while (p[start - 10] < p[start]) and (start > 10):
        start = start - 10

    while ((p[stop + 10] > p[stop])
            and (stop + 10 < len(p))):
        stop = stop + 10
    inds = range(start, stop)
    for tobin in ds.keys():
        logging.debug('Starting %s', tobin)
        start = timer()
        if ds[tobin].dims == ('scan',):
            dat, vvv, vv = bindata1d(zbins,
                    ds.pres.values[inds], ds[tobin].values[inds])
        else:
            dat = ds[tobin]
        end = timer()
        logging.debug('Done %s %1.2f', tobin, end-start)
        start = timer()
        if tobin in cgrid:
            logging.debug('append variable %s', tobin)
            if len(np.shape(dat)):
                ax = 1
                cgrid[tobin] = np.append(cgrid[tobin],
                        dat[:, np.newaxis], axis=ax)
            else:
                ax = 0
                dat = np.array([np.array(dat)])
                cgrid[tobin] = np.append(cgrid[tobin],
                                    dat, axis=ax)

        else:
            if len(np.shape(dat)):
                cgrid[tobin] = dat[:, np.newaxis]
            else:
                cgrid[tobin] = np.array([np.array(dat)])
        stop = timer()
        logging.debug('Concat time %1.2f', end-start)
    for towrite in ds.attrs.keys():
        logging.debug('Starting %s', towrite)
        if towrite in cgrid:
            cgrid[towrite] = np.append(cgrid[towrite], ds.attrs[towrite])
        else:
            cgrid[towrite] = ds.attrs[towrite]
    cgrid['name'] += [os.path.splitext(os.path.basename(name))[0]]
    logging.info('Done %s', name)

# sort by alongx:

ind = np.argsort(cgrid['alongx'])
for key in cgrid.keys():
    print(key)
    if isinstance(cgrid[key], np.ndarray):
        if len(cgrid[key].shape) == 2:
            cgrid[key] = cgrid[key][:, ind]
        elif cgrid[key].shape[0] == len(ind):
            cgrid[key] = cgrid[key][ind]
    else:
        cgrid[key] = [cgrid[key][i] for i in ind]

logging.debug(cgrid['alongx'])

# save to netcdf...
ctd = cgrid.copy()
ds = xr.Dataset(
    {
    'cond': (['depths', 'time'], ctd['cond'], {'units':'S/m', 'offset': '{} scans'.format(coffset)}),
    'cond0': (['depths', 'time'], ctd['cond0'], {'units':'S/m'}),
    'temp': (['depths', 'time'],ctd['temp'], {'units':'deg C'}),
    'pres': (['depths', 'time'],ctd['pres'], {'units':'dbar'}),
    'O2':  (['depths', 'time'],ctd['O2'], {'units':'mmol/kg'}),
    'O2sat':  (['depths', 'time'],ctd['O2sat'], {'units':'percent saturation'}),
    'Par':  (['depths', 'time'],ctd['Par'], {'units':'E/m^2'}),
    'Flu':  (['depths', 'time'],ctd['Flu'], {'units':'Flu'}),
    'sal':  (['depths', 'time'],ctd['sal'], {'units':'psu'}),
    'pden':  (['depths', 'time'],ctd['pden'], {'units':'kg/m^2'}),
    'lat': (['time'], ctd['lat'], {'units':'deg N'}),
    'lon': (['time'], ctd['lon'], {'units':'deg W'}),
    'alongx': (['time'], ctd['alongx'], {'units':'dist from S4 [km]'}),
    'acrossx': (['time'], ctd['acrossx'], {'units':'dist from S4 [km]'}),
    'id': (['time'], ctd['id'])
    },
    coords={'depths': (['depths'], ctd['depths']),
            'time': (['time'], ctd['time'])},
    attrs={}
)
ds.to_netcdf(dirname+'/CtdGrid.nc', 'w')
logging.debug(ds)
logging.info('Saved %s', dirname+'/CtdGrid.nc')

# save cgrid to matfile
cgridout = cgrid.copy()
cgridout['time'] = np.array([ctdtrans.datetime2matlab(
                    datetime.utcfromtimestamp(t.tolist()/1.e9))
        for t in cgrid['time']])
cgridout['den'] = cgridout['pden'].data
cgridout['t'] = cgridout['temp'].data
cgridout['c'] = cgridout['cond'].data
cgridout['c0'] = cgridout['cond0'].data
cgridout['flu'] = cgridout['Flu'].data
cgridout['p'] = cgridout['pres'].data
cgridout['O2'] = cgridout['O2'].data
cgridout['O2sat'] = cgridout['O2sat'].data
cgridout['par'] = cgridout['Par'].data
for topop in ('pden', 'temp', 'cond', 'cond0'):
    cgridout.pop(topop)
sio.savemat(dirname+'/CtdGridNew.mat', cgridout, format='5')
logging.info('Saved %s', dirname+'/CtdGridNew.mat')
ans = subprocess.call(['octave', 'saveasstruct.m', dirname+'/CtdGridNew.mat', 'cgrid'])
print(ans)
logging.info('Saved %s structured matlab', dirname+'/CtdGrid.mat')
