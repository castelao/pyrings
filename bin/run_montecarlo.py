#!/usr/bin/env python
# -*- coding: Latin-1 -*-
# vim: tabstop=4 shiftwidth=4 expandtab

""" Script to run Monte Carlo simulations with PyRings.

"""

from datetime import datetime, timedelta

import numpy as np
from numpy import ma

from rings.misc.montecarlo import *
#from rings.misc.montecarlo import montecarlo

# ============================================================================
cfg_base = {'ring':
            {'omega0': [1.7e-5, 2.7e-5], 
                'delta': [120e3, 200e3], #170e3, 
                'alpha': 2.5,
                'u_c': [-0.25, 0.25],
                'v_c': [-0.25, 0.25],
                'lat_t0': 8,
                'lon_t0': -50},
            'montecarlo':{
                'dt': 0,
                'Vnoise_sigma': [0.0, 0.20],
                'Nsamples': [25, 1000], 
                'Rlimit': [25e3, 300e3]
            }
        }


data = montecarlo(cfg_base, 10000)
import pandas as pd
data = pd.DataFrame(data)
data.describe()


store = pd.HDFStore('montecarlo.h5')
store.append('all_snapshot', data)
