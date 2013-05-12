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
            {'omega0': [1.4e-5, 2.5e-5], #2.0e-5
                'delta': [110e3, 200e3], #170e3, 
                'alpha': [3.0, 5.0], #4.8
                'u_c': [-0.20, 0.20],
                'v_c': [-0.20, 0.20]},
            'montecarlo':{
                'sampling_type': 'equal_area',
                'SamplingPeriod': [3600, 3600*24*30],
                #'dt': 3600,
                'Vnoise_sigma': [0.0, 0.20],
                'Nsamples': [10, 400], 
                'Rlimit': [25e3, 300e3]
            }
        }


import pandas as pd
Niterations = 10
version = '0.8.0'

data = montecarlo(cfg_base, Niterations)
data = pd.DataFrame(data)
#data.describe()

store = pd.HDFStore('montecarlo_%s.h5' % version)
store.append('trans_noise_equalarea', data)
store.close()

#==========================
cfg_base['montecarlo']['Vnoise_sigma'] = 0.0
data = montecarlo(cfg_base, Niterations)
data = pd.DataFrame(data)

store = pd.HDFStore('montecarlo_%s.h5' % version)
store.append('trans_nonoise_equalarea', data)
store.close()

#==========================
cfg_base['ring']['u_c'] = 0.0
cfg_base['ring']['v_c'] = 0.0
data = montecarlo(cfg_base, Niterations)
data = pd.DataFrame(data)

store = pd.HDFStore('montecarlo_%s.h5' % version)
store.append('notrans_nonoise_equalarea', data)
store.close()

#==========================
cfg_base['montecarlo']['Vnoise_sigma'] = [0.0, 0.20]
data = montecarlo(cfg_base, Niterations)
data = pd.DataFrame(data)

store = pd.HDFStore('montecarlo_%s.h5' % version)
store.append('notrans_noise_equalarea', data)
store.close()

