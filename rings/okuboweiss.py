#!/usr/bin/env python
# -*- coding: Latin-1 -*-
# vim: tabstop=4 shiftwidth=4 expandtab

""" Estimate de Okubo-Weiss parameter and handle related variables/procedures.
"""

import logging
import multiprocessing
from UserDict import UserDict
#from UserDict import IterableUserDict

import numpy
from numpy import ma

import fluid

from rings.utils import basic_logger
from rings.utils import cfg2dict


def okuboweiss(input):
    """ Estimate W and zeta

        In theory I should be able to use the 2D simplification, otherwise
          if it's a compressible flow the whole concept wouldn't be appliable.
    """
    from fluid.common.common import FieldDiferentials
    output = {}

    differentials = FieldDiferentials(input, dims={'x': 1, 'y': 0})
    output['zeta'] = differentials['dvdx'] - differentials['dudy']

    output['W'] = differentials['dudx']**2 \
      - 2*differentials['dudx']*differentials['dvdy'] \
      + differentials['dvdy']**2 \
      + 4*differentials['dvdx']*differentials['dudy']

    # The simplification for the 2D incompressible flow
    #self.data['W'] = 4*(differentials['dudx']**2 + \
    #  differentials['dvdx']*differentials['dudy'])

    return output


# ----------------------------------------------------------------------------
class OkuboWeiss(UserDict):
    """ Class to estimate the Okubo-Weiss (OW) index from a regular
          velocity field

        Big change from the first version. OW was supposed to be applied
          in a field, therefore it don't deal anymore with 3rd dimension,
          i.e. 2D fields only. Whatever is using the OW, the higher level
          procedures that should deal with the extra dimensions (t,z...).
          Classes and functions should be as simple and objective as
          possible.

        The attribute parameters is ready to receive in the future
          extra variables, like if is to smooth, with each method and
          which cutoff frequency, or which method to estimate W0.
    """
    def __init__(self, input, metadata={'W0': 'chelton'}, logname=None, **keywords):
        """
        """
        self._set_logger(logname)

        self.input = input
        # If data is an input, use it, otherwise create and empty dicitonary
        if 'data' in keywords.keys():
            self.data = keywords['data']
        else:
            self.data = {}

        self.metadata = metadata

        self.go()

    def _set_logger(self, logname):
        self.logname = logname
        if logname is None:
            self.logger = basic_logger(logname="Generic")
        else:
            try:
                self.logger = logging.getLogger(logname)
            except:
                self.logger = basic_logger(logname)
        self.logger.info("Inside class OkuboWeiss")

    def go(self):
        """
        """
        for v in ['W', 'zeta']:
            if v not in self.data.keys():
                self.data[v] = ma.masked_all(self.input['u'].shape)

        if len(self.input['u'].shape) == 3 and len(self.input['v'].shape) == 3:
            self.logger.info("The fields u and v are 3D. \
                    I'll repeat along the first dimension.")
            nt, ni, nj = self.input['u'].shape

            npes = 2*multiprocessing.cpu_count()
            self.logger.info("I'll work with %s parallel processes" % npes)
            pool = multiprocessing.Pool(npes)
            results = []

            self.logger.debug("I'm about to fire the pool processes")
            for t in range(nt):
                snapshot = {'Lon': self.input['Lon'], 'Lat': self.input['Lat'],
                        'u': self.input['u'][t], 'v': self.input['v'][t]}
                results.append( pool.apply_async( okuboweiss, (snapshot,) ) )
            pool.close()

            for t, r in enumerate(results):
                self.logger.debug("Getting the OkuboWeiss parameter: %s/%s" % (t,nt))
                output = r.get()
                self.data['W'][t] = output['W']
                self.data['zeta'][t] = output['zeta']

        elif len(self.input['u'].shape)==2 and len(self.input['v'].shape)==2:
            self.logger.debug("Calculating the OkuboWeiss parameter")
            self.data['W'], self.data['zeta'] = okuboweiss(self.input)

        else:
            print "Error, u and v must be 2D or 3D"

        self.set_W0()

        self.logger.info("I'm done with class OkuboWeiss")

    def set_sigma_W(self):
        """ Pseudo-Standart deviation of W index

            ATENTION, THIS IS WRONG!!!
        """
        W = self.data['W'].compressed().data
        self.data['sigma_W'] = (W[round(len(W)*.75-1)]-W[round(len(W)*.25-1)])/1.349
        return

    def set_W0(self):
        """ W_0 is the reference for the background flow
        """
        if self.metadata['W0'] == 'chelton':
            self.data['W0'] = numpy.ones(N)*2e-12
        elif self.metadata['W0'] == 'fontanet':
            self.data['W0'] = 0.2*self.data['W'].std()
        elif type(self.metadata['W0']) == float:
            self.data['W0'] = self.metadata['W0']
        return
