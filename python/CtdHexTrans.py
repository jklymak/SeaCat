import numpy as np
import logging
from datetime import datetime, timedelta
import seawater.eos80 as sw
import seawater
import xarray as xr
import scipy.io as sio



log = logging.getLogger(__name__)


def ctd_hex_trans(fname):
    """
    Translate a seacat hex file to a python dictionary

    Parameters
    ----------

    fname: string
        The name of the file to read

    Returns
    -------

    seacat: dict
        A dictionary with the raw data in the file as elements.
    """

    seacat = dict()
    seacat['headinfo']=''
    seacat['CTDCoeff'] = dict()
    seacat['t'] = np.zeros(0)
    seacat['c'] = np.zeros(0)
    seacat['p'] = np.zeros(0)
    seacat['v1'] = np.zeros(0)
    seacat['v2'] = np.zeros(0)
    seacat['v3'] = np.zeros(0)

    with open(fname,'rt') as fin:
        for l in fin:
            log.debug(len(l))
            if len(l) == 1:
                pass
            elif (l[0] == '*'):
                log.debug(l)
                seacat['headinfo'] += l+'\n'
                log.debug(l[0:5])
                if l[0:5] == '*    ':
                    i = l.find('=')
                    log.debug('Hi %d', i)
                    if i:
                        v=l[1:i]
                        val = float(l[i+1:])
                        seacat['CTDCoeff'][v.lower().strip()] = val
                        log.debug(val)
                        log.debug("Name %s", v.lower().strip())
            else:
                break

        for l in fin:
            # temperature:
            a0 = seacat['CTDCoeff']['ta0']
            a1 = seacat['CTDCoeff']['ta1']
            a2 = seacat['CTDCoeff']['ta2']
            a3 = seacat['CTDCoeff']['ta3']
            log.debug(l[0:6])
            rawT = int(l[0:6], 16)
            log.debug(rawT)
            mv = (rawT-524288)/1.6e7;
            r = (mv*2.9e9 + 1.024e8) / (2.048e4 - mv * 2e5);
            t = a0 + a1 * np.log(r) + a2 * np.log(r)**2 + a3 * np.log(r)**3;
            t = 1./t - 273.15;
            seacat['t'] = np.append(seacat['t'], t)

            # pressure

            pa0 = seacat['CTDCoeff']['pa0']
            pa1 = seacat['CTDCoeff']['pa1']
            pa2 = seacat['CTDCoeff']['pa2']
            ptempa0 = seacat['CTDCoeff']['ptempa0']
            ptempa1 = seacat['CTDCoeff']['ptempa1']
            ptempa2 = seacat['CTDCoeff']['ptempa2']
            ptca0 = seacat['CTDCoeff']['ptca0']
            ptca1 = seacat['CTDCoeff']['ptca1']
            ptca2 = seacat['CTDCoeff']['ptca2']
            ptcb0 = seacat['CTDCoeff']['ptcb0']
            ptcb1 = seacat['CTDCoeff']['ptcb1']
            ptcb2 = seacat['CTDCoeff']['ptcb2']

            rawP = int(l[12:18], 16)

            y = int(l[18:22], 16) / 13107.
            t = ptempa0 + ptempa1 * y + ptempa2 * y**2
            x = rawP - ptca0 - ptca1 * t - ptca2 * t**2
            n = x*ptcb0 / (ptcb0 + ptcb1 * t + ptcb2 * t**2)

            seacat['p'] = np.append(seacat['p'],
                        (pa0 + pa1 * n + pa2 * n**2 - 14.7) * 0.689476)

            # conductivity
            g = seacat['CTDCoeff']['g']
            h = seacat['CTDCoeff']['h']
            i = seacat['CTDCoeff']['i']
            j = seacat['CTDCoeff']['j']
            tcor = seacat['CTDCoeff']['ctcor']
            pcor = seacat['CTDCoeff']['cpcor']
            #
            f = int(l[6:12], 16) / 256. / 1000.

            seacat['c'] = np.append(seacat['c'],
                    (g+h*f**2 + i * f**3 + j *
                    f**4) / (1. + tcor * seacat['t'][-1]
                    + pcor * seacat['p'][-1]))

            # voltages:

            seacat['v1'] = np.append(seacat['v1'], int(l[22:26], 16)/13107.)
            seacat['v2'] = np.append(seacat['v2'], int(l[26:30], 16)/13107.)
            seacat['v3'] = np.append(seacat['v3'], int(l[30:34], 16)/13107.)

        log.debug(seacat['CTDCoeff'].keys())
        log.debug(seacat['t'][100:200])
        log.debug(seacat['p'][200:3000])
        log.debug(seacat['c'][200:3000])

        return seacat

def datetime2matlab(dt):
    mdn = dt + timedelta(days = 366)
    frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds

    return mdn.toordinal() + frac / (24.0 * 60.0 * 60.0)

def getFlu(ctd,vvar):

    scale = 15.5 # ug/L/V
    Voff = -0.053 # V
    ctd['Flu'] = scale * (ctd[vvar] - Voff)

    return ctd


def getPar(ctd, vvar):

    ConFac = 10111223458.03842
    offset = -0.15423974

    ctd['Par'] =  (ctd[vvar] - offset) *  ConFac

    return ctd

def getO2(ctd,vvar):

    if 'sal' not in ctd:
        ctd['sal'] = sw.salt(ctd['c'] * 10. / seawater.constants.c3515,
                ctd['t'], ctd['p'])

    ctd['O2sat'] = np.real( seawater.extras.satO2(ctd['sal'], ctd['t']))
    # from April 2008 calibration
    Voff = -0.4877
    Soc = 0.4401
    Boc = 0.0000
    Tcor = 0.0005
    Pcor = 1.35e-4

    #below numbers are true for CTD casts prior to April 2008
    #Voff = -0.4893;
    #Soc = 0.4227;
    #Boc = 0.0000;
    #Tcor = 0.001;
    #Pcor = 1.35e-4;

    ctd['O2'] = (Soc * (ctd[vvar] + Voff) * np.exp(Tcor * ctd['t']) *
            ctd['O2sat'] * np.exp(ctd['p'] * Pcor))

    if 'den' not in ctd:
        ctd['pden']=sw.dens(ctd['sal'], ctd['t'], ctd['p'])
    ctd['O2'] = 1.e6 * ctd['O2']/(ctd['pden'] * 22.3916);
    return ctd

def raw_to_science(ctd, coffset):
    """
    Parse the header and convert ctd data to scientifically useful data.

    Parameters
    ----------

    ctd: dict
        Dictionary created by `CtdHexTrans`

    coffset: float
        Amount to offset conductivity by in scans.  >0 means conductivity
        is moved backwards in time to match temperature (i.e. temperature
        is lagging)  <0 means conductivity is moved forward (i.e.
        conductivity is lagging)

    Returns
    -------
    ctd: dict
        with new info (sal, pden, voltages to data, lat, lon, time)

    """

    # parse some header stuff
    for l in ctd['headinfo'].splitlines():
        try:

            if (l.find('** Lon') == 0):
                lon, minute = (l.split(':')[1]).split(' ')
                ctd['lon'] = -float(lon) - float(minute)/60.
            if (l.find('** Lat') == 0):
                lat, minute = (l.split(':')[1]).split(' ')
                ctd['lat'] = float(lat) + float(minute)/60.
            if (l.find('** Sta') == 0):
                ctd['id'] = l.split(':')[1]
            if (l.find('* cast') == 0 ):
                timest = l[11:28]
                ctd['time'] = datetime.strptime(timest, '%d %b %Y %H:%M')
                ctd['matlabtime'] = datetime2matlab(ctd['time'])
        except:
            pass

    #
    # get salt and pden...
    ctd['c0'] = ctd['c']
    N = len(ctd['c'])
    ctd['coffset'] = coffset
    ctd['c'] = np.interp(np.arange(N)-coffset, np.arange(N), ctd['c0'])
    ctd['sal'] = sw.salt(ctd['c'] * 10. / seawater.constants.c3515,
            ctd['t'], ctd['p'])
    ctd['pden']=sw.dens(ctd['sal'], ctd['t'], ctd['p'])

    ctd = getO2(ctd,'v1');
    ctd = getFlu(ctd,'v2');
    ctd = getPar(ctd,'v3');

    return ctd

def tonetcdf(ncname, ctd, coffset):
    """
    Save ctd as a netcdf file.

    This routine is pretty hardwired to what came before it.
    """
    logging.getLogger(__name__)


    ds = xr.Dataset(
        {
        'cond': (['scan'], ctd['c'], {'units':'S/m', 'offset': '{} scans'.format(coffset)}),
        'cond0': (['scan'], ctd['c0'], {'units':'S/m'}),
        'temp': (['scan'],ctd['t'], {'units':'deg C'}),
        'pres': (['scan'],ctd['p'], {'units':'dbar'}),
        'O2':  (['scan'],ctd['O2'], {'units':'mmol/kg'}),
        'O2sat':  (['scan'],ctd['O2sat'], {'units':'percent saturation'}),
        'Par':  (['scan'],ctd['Par'], {'units':'E/m^2'}),
        'Flu':  (['scan'],ctd['Flu'], {'units':'Flu'}),
        'sal':  (['scan'],ctd['sal'], {'units':'psu'}),
        'pden':  (['scan'],ctd['pden'], {'units':'kg/m^2'}),
        'time': ctd['time'],
        },
        coords={'scan': (['scan'], np.arange(len(ctd['c'])))},
        attrs={
        'alongx':ctd['alongx'],
        'acrossx':ctd['acrossx'],
        'lat':ctd['lat'],
        'lon':ctd['lon'],
        'id':ctd['id'],
        'header':ctd['headinfo']}
    )
    logging.debug('ds:')
    logging.debug(ds)

    ds.to_netcdf(ncname)
    logging.info('saved nc')


def tomatlab(matname, ctd):
    ctdout = ctd.copy()
    ctdout.pop('time')
    ctdout.pop('matlabtime')
    ctdout['den'] = ctdout['pden']
    ctdout.pop('pden')
    ctdout['time'] = ctd['matlabtime']

    sio.savemat(matname, ctdout, format='5')
    # subprocess.call(['octave', 'saveasstruct.m', dirname+'/'+basename+'New.mat', 'ctd'])
