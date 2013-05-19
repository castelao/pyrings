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

from rings.ring import RingCenter
from rings.fitt import carton_uv, carton_scales
from rings.utils import xy2lonlat, lonlat2xy


def synthetic_CLring(x, y, t, cfg):
    """ Creates a synthetic ring like Carton n Legras

        Input:
          x: x [m]
          y: y [m]
          t: t [s]

        The parameters are passed by a dictionary: cfg

        The center is at (0,0) at t=0
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

def drunken_drive(cfg):
    """ Simulates a standard normal change in course

    For N iterations it changes the course with a normal probability around
      the previous course. Each time moves one step distance. On the drunken
      drive there is a memory on the course followed, different from the
      drunken walk, and so it is more similar to real sensors (ships, drifters
      etc).

      Variance between -90 to 90 deg.
    """
    N = cfg['montecarlo']['Nsamples']
    Rlimit = cfg['montecarlo']['Rlimit']
    x0 = Rlimit*(random(1)[0]-0.5)*2
    y0 = Rlimit*(random(1)[0]-0.5)*2
    rad = 0.25*np.pi*randn(N)
    rad[0] = random(1)*2*np.pi
    rad = rad.cumsum()
    step = cfg['montecarlo']['dt'] * cfg['montecarlo']['VSampler']
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

def drifter_sample(cfg):
    N = cfg['montecarlo']['Nsamples']

    # Initial position
    t = np.arange(N)*cfg['montecarlo']['dt']
    x, y = ma.masked_all(N), ma.masked_all(N),
    u, v = ma.masked_all(N), ma.masked_all(N),

    x[0], y[0] = (2*random(2)-1)*cfg['montecarlo']['Rlimit']
    u[0], v[0] = synthetic_CLring(x[0], y[0], t[0], cfg['ring'])

    u_nonoise, v_nonoise = u.copy(), v.copy()

    for n in range(1,N): 
        dt = t[n] - t[n-1]
        x[n], y[n] = x[n-1]+u[n-1]*dt, y[n-1]+v[n-1]*dt
        u_nonoise[n], v_nonoise[n] = synthetic_CLring(x[n], y[n], t[n], cfg['ring'])
        u[n] = u_nonoise[n] + cfg['montecarlo']['Vnoise_sigma'] * randn()
        v[n] = v_nonoise[n] + cfg['montecarlo']['Vnoise_sigma'] * randn()
        #print x[n], y[n], u[n], v[n]
    #t = t - np.median(t)
    data= {'t': t, 'x':x, 'y':y, 'u':u, 'v':v}
    stats = {}
    return data, stats

def sampler(cfg):
    """ Simulates sampling

        Input:
            - t: Samples timing [s]
            - method:
                regulargrid
                equal_area
                drunken_drive
                drifter
            - Rlimit: Radial limit of the samples

    """
    # Number of samples
    N = cfg['montecarlo']['Nsamples'] 

    # Defines the timing of the samples
    if ('dt' not in cfg['montecarlo']) & \
    ('SamplingPeriod' in cfg['montecarlo']):
        cfg['montecarlo']['dt'] = cfg['montecarlo']['SamplingPeriod']/(N-1.)
    t = (np.arange(N, dtype='i')-N/2.)*cfg['montecarlo']['dt']

    # Define the (x,y) sampling positions
    if cfg['montecarlo']['sampling_type'] == 'regulargrid':
        x, y = regulargrid_sample(N, cfg['montecarlo']['Rlimit'])
    elif cfg['montecarlo']['sampling_type'] == 'equal_area':
        x, y = random_sample_equal_area(N , cfg['montecarlo']['Rlimit'])
    elif cfg['montecarlo']['sampling_type'] == 'drifter':
        return drifter_sample(cfg)
    elif cfg['montecarlo']['sampling_type'] == 'drunken_drive':
        x, y = drunken_drive(cfg) #, x0=0, y0=0)
    else:
        return

    #t = t - np.median(t)

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
    #d0 = datetime(2000,1,1)
    #d = d0+ma.array([timedelta(seconds=dt) for dt in t])
    # Transform x,y into lon, lat
    #lon, lat = xy2lonlat(x, y, cfg['ring']['lon_t0'], cfg['ring']['lat_t0'])
    Vmax, Rmax = \
                carton_scales(cfg['ring']['omega0'], cfg['ring']['delta'],
                        cfg['ring']['alpha'])
    # Creating input to Class Ring
    #input = {'datetime': d, 'Lon': lon, 'Lat': lat, 'u':u, 'v':v}
    #input = {'datetime': d, 'x': x, 'y': y, 'u':u, 'v':v}
    data = {'t': t, 'x': x, 'y': y, 'u':u, 'v':v}
    stats = {'Rmedian': Rmedian, 'Vmedian': Vmedian,
            'noisesig_ratio': noisesig_ratio,
            'x_median': np.median(x), 'y_median': np.median(y),
            'Vmax': Vmax, 'Rmax': Rmax,
            }

    return data, stats

def error_estimate(cfg):

    data, stats = sampler(cfg)
            #method=cfg['montecarlo']['sampling_type'],
            #Rlimit = cfg['montecarlo']['Rlimit'])

    anel = RingCenter(data)
    # error estimate
    xc_err = anel.center['x'] - cfg['ring']['u_c'] * anel.center['t']
    yc_err = anel.center['y'] - cfg['ring']['v_c'] * anel.center['t']
    #xc_err, yc_err = lonlat2xy(ma.array([anel['Lon_c']]), ma.array([anel['Lat_c']]), [cfg['ring']['lon_t0']], [cfg['ring']['lat_t0']])
    uc_err = anel.center['u'] - cfg['ring']['u_c']
    vc_err = anel.center['v'] - cfg['ring']['v_c']

    output = cfg.copy()
    output['output'] = {'xc_err': xc_err, 'yc_err': yc_err,
            'uc_err': uc_err, 'vc_err': vc_err,
            'opt_stat': anel.opt_stat,
            'dt': cfg['montecarlo']['dt']}
    for k in stats.keys():
        output['output'][k] = stats[k]

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
