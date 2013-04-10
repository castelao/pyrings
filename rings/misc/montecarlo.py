#!/usr/bin/env python
# -*- coding: Latin-1 -*-
# vim: tabstop=4 shiftwidth=4 expandtab

""" Miscelanous functions for the Monte Carlo simulations with PyRings

"""

from datetime import datetime, timedelta
import multiprocessing

import numpy as np
from numpy import ma
from numpy.random import random, randn

from rings.ring import Ring
from rings.fitt import carton_uv, carton_scales
from rings.utils import xy2lonlat, lonlat2xy


def synthetic_CLring(x, y, t, cfg):
    """ Creates a synthetic ring like Carton n Legras

        The parameters are passed by a dictionary: cfg
    """
    x_c = cfg['u_c'] * t  # + cfg['xt0']
    y_c = cfg['v_c'] * t  # + cfg['yt0']
    u, v = carton_uv(x-x_c, y-y_c, cfg['omega0'], cfg['delta'], cfg['alpha'])
    return u+cfg['u_c'], v+cfg['v_c']

def drunken_walk(N, step=1, x0=0, y0=0):
    """ Simulates a random walker

    For N iterations, it turns in any direction, equally distributed. Each
    time moves one step
    """
    rad = random(N)*2*np.pi
    dx = ma.sin(rad)*step
    dy = ma.cos(rad)*step
    x = dx.cumsum() + x0
    y = dy.cumsum() + y0
    return x, y

def drunken_drive(N, step=1, x0=0, y0=0):
    """ Simulates a standard normal change in course

    For N iterations it changes the course with a normal probability around
      the previous course. Each time moves one step distance. On the drunken
      drive there is a memory on the course followed, different from the
      drunken walk, and so it is more similar to real sensors (ships, drifters
      etc).

      Variance between -90 to 90 deg.
    """
    rad = 0.25*np.pi*randn(N)
    rad[0] = random(1)*2*np.pi
    rad = rad.cumsum()
    dx = ma.sin(rad)*step
    dy = ma.cos(rad)*step
    x = dx.cumsum() + x0
    y = dy.cumsum() + y0
    return x, y

def random_sample_equal_area(N, Rlimit):
    """ Return x,y for N samples inside the Rlimit.

        The probability is equal to area.
    """
    x = ma.array(Rlimit*(random(N)-0.5)*2)
    y = ma.array(Rlimit*(random(N)-0.5)*2)
    r = (x**2+y**2)**0.5
    ind = np.nonzero(r>Rlimit)[0]
    while ind.any():
        x[ind] = Rlimit*(random(len(ind))-0.5)*2
        y[ind] = Rlimit*(random(len(ind))-0.5)*2
        r = (x**2+y**2)**0.5
        ind = np.nonzero(r>Rlimit)[0]
    return x, y

def regulargrid_sample(N, Rlimit):
    """
    """
    grid_side = int(np.ceil(N**0.5))
    x = np.linspace(-Rlimit,Rlimit,grid_side)
    x, y = np.meshgrid(x, x)
    r2 = (x**2 + y**2)
    #ind = np.argsort(r2)<N # Starts on 0, so < not <=
    ind = np.argsort(r2.flatten())[:N]
    return ma.array(x.flatten()[ind]), \
            ma.array(y.flatten()[ind])

def drifter_sample(N, cfg):
    dt = 6*3600 #cfg['montecarlo']['dt']
    # Initial position
    t = np.arange(N)*cfg['montecarlo']['dt']
    t = t - np.median(t)
    x, y, u, v = ma.masked_all(N), ma.masked_all(N),ma.masked_all(N),ma.masked_all(N),
    x[0], y[0] = (2*random(2)-1)*cfg['montecarlo']['Rlimit']
    u[0], v[0] = synthetic_CLring(x[0], y[0], t[0], cfg['ring'])
    for n in range(1,N): 
        x[n], y[n] = x[n-1]+u[n-1]*dt, y[n-1]+v[n-1]*dt
        u[n], v[n] = synthetic_CLring(x[n], y[n], t[n], cfg['ring'])


def error_estimate(cfg):
    N = cfg['montecarlo']['Nsamples'] 
    t = np.arange(N, dtype='i')*cfg['montecarlo']['dt']
    t = t - np.median(t)
    # Define the (x,y) sampling positions
    if cfg['montecarlo']['sampling_type'] == 'equal_area':
        x, y = random_sample_equal_area(N , cfg['montecarlo']['Rlimit'])
    elif cfg['montecarlo']['sampling_type'] == 'drunken_drive':
        x, y = drunken_drive(N, step=1) #, x0=0, y0=0)
    elif cfg['montecarlo']['sampling_type'] == 'regulargrid':
        x, y = regulargrid_sample(N, cfg['montecarlo']['Rlimit'])

    Rmedian = np.median((x**2+y**2)**0.5)
    # Estimate the measures
    u, v = synthetic_CLring(x, y, t, cfg['ring'])
    Vmedian = np.median((u**2+v**2)**0.5)
    # Add a noise
    u_noise = cfg['montecarlo']['Vnoise_sigma'] * randn(N)
    v_noise = cfg['montecarlo']['Vnoise_sigma'] * randn(N)
    noisesig_ratio = np.median(
            (u_noise**2 + v_noise**2)**0.5 / \
            (u**2 + v**2)**0.5
            )
    u = u + u_noise
    v = v + v_noise
    # Might have a better way to do the line below.
    d0 = datetime(2000,1,1)
    d = d0+ma.array([timedelta(seconds=dt) for dt in t])
    # Transform x,y into lon, lat
    lon, lat = xy2lonlat(x, y, cfg['ring']['lon_t0'], cfg['ring']['lat_t0'])
    # Creating input to Class Ring
    input = {'datetime': d, 'Lon': lon, 'Lat': lat, 'u':u, 'v':v}
    anel = Ring(input)
    # error estimate
    xc_err, yc_err = lonlat2xy(ma.array([anel['Lon_c']]), ma.array([anel['Lat_c']]), [cfg['ring']['lon_t0']], [cfg['ring']['lat_t0']])
    uc_err = anel['uc'] - cfg['ring']['u_c']
    vc_err = anel['vc'] - cfg['ring']['v_c']

    Vmax, Rmax = \
                carton_scales(cfg['ring']['omega0'], cfg['ring']['delta'],
                        cfg['ring']['alpha'])
    output = cfg.copy()
    output['output'] = {'xc_err': xc_err[0], 'yc_err': yc_err[0],
            'x_median': np.median(x), 'y_median': np.median(y),
            'uc_err': uc_err, 'vc_err': vc_err, 'Vmax': Vmax, 
            'Rmax': Rmax, 'noisesig_ratio': noisesig_ratio,
            'Rmedian': Rmedian, 'Vmedian': Vmedian,
            'opt_stat': anel.opt_stat}

    return output

def random_cfg(cfg):
    cfg2 = {}
    for k in cfg:
        if type(cfg[k]) == list:
            cfg2[k] = (cfg[k][1]-cfg[k][0])*random(1)[0]+cfg[k][0]
            if (type(cfg[k][0])==int) & (type(cfg[k][1])==int):
                    cfg2[k] = int(cfg2[k])
        elif type(cfg[k]) == dict:
            cfg2[k] = random_cfg(cfg[k])
        else:
            cfg2[k] = cfg[k]
    return cfg2

#import pandas as pd
def montecarlo(cfg_base, N):
    """
    """
    npes = 2 * multiprocessing.cpu_count()
    pool = multiprocessing.Pool(npes)
    results = []

    for n in range(N):
        cfg_tmp = random_cfg(cfg_base)
        # Temporary solution
        #output = error_estimate(cfg_tmp)
        results.append( pool.apply_async( error_estimate, (cfg_tmp,) ) )
    pool.close()

    data = []
    for i, r in enumerate(results):
        output = r.get()
        tmp = {}
        for k in output.keys():
            for kk in output[k].keys():
                tmp[kk] = output[k][kk]
        data.append(tmp)
    pool.terminate()

    #data = pd.DataFrame(data)
    return data
