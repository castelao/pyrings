#!/usr/bin/env python
# -*- coding: Latin-1 -*-
# vim: tabstop=4 shiftwidth=4 expandtab

""" Script to run Monte Carlo simulations with PyRings.

    The misc directory is probably a better place for this code.

"""

from datetime import datetime, timedelta

import numpy as np
from numpy import ma

from rings.misc.montecarlo import *
#from rings.misc.montecarlo import montecarlo

# ============================================================================
cfg_base = {'ring':
            {'omega0': [1.5e-5, 2.6e-5], #2.0e-5
                'delta': [120e3, 200e3], #170e3, 
                'alpha': [3.0, 5.0], #4.8
                'u_c': [-0.20, 0.20],
                'v_c': [-0.20, 0.20],
                'lat_t0': 8,
                'lon_t0': -50},
            'montecarlo':{
                'sampling_type': 'regulargrid',
                'dt': 0,
                'Vnoise_sigma': 0, #[0.0, 0.20],
                'Nsamples': [10, 500], 
                'Rlimit': [25e3, 300e3]
            }
        }


data = montecarlo(cfg_base, 10000)
import pandas as pd
data = pd.DataFrame(data)
data.describe()


store = pd.HDFStore('montecarlo.h5')
store.append('all_snapshot', data)

#data['Lc_err'] = (data['xc_err']**2 + data['yc_err']**2)**0.5
#data['Vc_err'] = (data['uc_err']**2 + data['vc_err']**2)**0.5
