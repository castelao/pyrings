#!/usr/bin/env python
# -*- coding: Latin-1 -*-
# vim: tabstop=4 shiftwidth=4 expandtab

"""Ocean Rings related stuff.

     Atention!!! Work on variables _r. Maybe I should use instead of _r
       for the ring variables, use _raw for the input, or a dictionary
       data_intput, and use x,r,u,v for the data corrected for the ring.


    Should strongly consider to use namedtuple, for efficiency and clarity
      purposes.
    http://docs.python.org/library/collections.html#collections.namedtuple

    Another idea, include uncertain on the data. That would allow to use
      different datasets and give different weights for each type of data.
"""

#import logging
from datetime import datetime, timedelta

import numpy as np
from numpy import ma
from scipy import optimize

from rings.utils import uv2nt
from rings import fitt
#from fluid.common.common import lonlat2xy, xy2lonlat
from rings.utils import lonlat2xy, xy2lonlat


class Ring(object):
    """ A class to deal with Rings

        Maybe should split in two classes, one for eta inputs and other for
          velocities input.

        Input:
          data:
              | Lat: Latitude
              | Lon: Longitude
              | or
              | X: x position
              | Y: y position

              | datetime: Time [datetime]
              | or
              | seconds
            u: u component of velocity
            v: v component of velocity
    """
    pass

class RingCenter(object):
    """ A class to find the ring center.

        Maybe should split in two classes, one for eta inputs and other for
          velocities input.

        Input:
          data:
            t: time [seconds]
            x: x position
            y: y position
            u: u component of velocity [m/s]
            v: v component of velocity [m/s [m/s]]
    """

    def __init__(self, input, metadata={}, auto=True, **keywords):
        """
        """
        self.input = input
        self.data = {}
        self.metadata = metadata

        self._set_default_values(keywords)

        if auto == True:
            #self.set_t()
            self.findcenter()
            self.set_ring_velocity()
            self.set_cilyndrical_components()

    def keys(self):
        """ .keys is composed from input + data elements
        """
        k = list(self.data.keys())
        k.extend(self.input.keys())
        return k

    def __getitem__(self, key):
        try:
            return self.data[key]
        except:
            return self.input[key]

    def _set_default_values(self, keywords):
        """
        """

        if 'center' in keywords:
            self.center = keywords['center']
        else:
            self.center = {}

    def set_cilyndrical_components(self):
        self.data['r'] = (self['x'] ** 2 + self['y'] ** 2) ** 0.5
        self.data['vrad'], self.data['vtan'] = uv2nt(self['x'], self['y'],
                self['u'], self['v'],
                x_c=self.center['x'], y_c=self.center['y'])

    def findcenter(self):
        """
optimize.tnc.RCSTRINGS
EINVAL       = -2 # Invalid parameters (n<1)
INFEASIBLE   = -1 # Infeasible (low > up)
LOCALMINIMUM =  0 # Local minima reach (|pg| ~= 0)
CONVERGED    =  1 # Converged (|f_n-f_(n-1)| ~= 0)
MAXFUN       =  2 # Max. number of function evaluations reach
LSFAIL       =  3 # Linear search failed
CONSTANT     =  4 # All lower bounds are equal to the upper bounds
NOPROGRESS   =  5 # Unable to progress
USERABORT    =  6 # User requested end of minimization
        """
        verbose = 0
        tref = np.median(self['t'])
        self.center['t'] = tref
        args = ((self['t'] - tref), self['x'], self['y'],
                self['u'], self['v'])
        f = fitt.v_circular()
        #f.lamb = 1e-2
        #f = fitt.v_circular_nontranslating()
        #f.set_p0(self.data['x'], self.data['y'],
        #        self.input['u'], self.input['v'])
        bounds = None  #[(None,None), (None,None), (None,None), (None,None)]
        op = optimize.fmin_tnc(f.cost, f.p0, fprime=None,  args=args,
                approx_grad=1,  bounds=bounds, epsilon=1e-08, scale=None,
                offset=None, messages=15, maxCGit=-1, maxfun=500, eta=-1,
                stepmx=0, accuracy=0, fmin=0, ftol=-1, xtol=-1, pgtol=-1,
                rescale=-1, disp=verbose)

        self.opt_stat = op[2]

        p = op[0]   # Output fitting parameters
        self.center['x'] = f.s[0] * p[0]
        self.center['y'] = f.s[1] * p[1]
        self.center['u'] = f.s[2] * p[2]
        self.center['v'] = f.s[3] * p[3]

    def set_ring_velocity(self):
        """
        """
        self.data['xr'] = self.input['x'] - self['t'] * self.center['u']
        self.data['yr'] = self.input['y'] - self['t'] * self.center['v']
        self.data['ur'] = self.input['u'] - self.center['u']
        self.data['vr'] = self.input['v'] - self.center['v']


class RingCenterFlex(RingCenter):
    """ A class to find the ring center, with flexible inputs.

        Input:
          data:
            | datetime: Time [datetime]
            or
            | t: time [seconds]

            | Lat: Latitude
            | Lon: Longitude
            or
            | x: x position
            | y: y position

            u: u component of velocity [m/s]
            v: v component of velocity [m/s]
    """
    def __init__(self, input, metadata={}, auto=True, **keywords):
        """
        """
        # I don't think this keywords will work properly
        super(RingCenterFlex, self).__init__(input, metadata, auto=False)
        self.name = 'RingCenterFlex'

        if auto == True:
            #self.set_xy()
            self.set_t()

            self.findcenter()
            #self.data['Lon_c'], self.data['Lat_c'] = xy2lonlat(self.data['xc'],
            #        self.data['yc'], self.lon_ref, self.lat_ref)


            #self.set_xy()   # Redefine the positions, discounting the uc|vc
            self.set_ring_velocity()
            self.set_cilyndrical_components()

    def _set_default_values(self, keywords):
        """
        """

        if 'center' in keywords:
            self.center = keywords['center']
        else:
            self.center = {}

        # Lat, Lon of the origin of the coordinate system, i.e.  (x,y)
        #if (not hasattr(self, 'lat_ref')) & (not hasattr(self, 'lon_ref')):
        #    try:
        #        self.lat_ref = self.center['lat']
        #        self.lon_ref = self.center['lon']
        #    #except AttributeError:
        #    except KeyError:
        #        self.lat_ref = ma.median(self.input['Lat'])
        #        self.lon_ref = ma.median(self.input['Lon'])

    def set_t(self):
        """ Create an array t as seconds from self.t_ref

            If t_ref is not defined, it uses the median of the
              input['datetime'], than create an array data['t'] with the
              seconds relative to t_ref.

            It's not the most efficient way, but is done to work well with any
              datetime dimension.

            Should I convert datetime to t using the very first date?
              I might avoid some errors doing on that way.
        """
        if ('t' not in self.keys()) & ('datetime' in self.input):
            t0 = min(self['datetime'])
            t = ma.array([dt.total_seconds() for dt in self['datetime'] - t0])
            t_median = ma.median(t)
            #self.t_ref = t0+timedelta(seconds = dt_median)
            self.data['t'] = t - t_median
        #self.data['t'] = self.input['t'] - np.median(self.input['t'])

        self.ref['datetime'] = t0 + timedelta(seconds=t_median)

    def set_xy(self):
        """ Set x, y coordinates from Lat, Lon.

            Most of the procedures are based on cartesian distances. This
              function define x and y distance from lat_ref, lon_ref.
        """
        if ('x' in self.input) & ('y' in self.input):
            self.data['x'] = self.input['x']
            self.data['y'] = self.input['y']
            print "It will be used the x and y given as input"

        else:
            self.data['x'], self.data['y'] = lonlat2xy(self.input['Lon'], 
                    self.input['Lat'], lat_0=self.lat_ref, lon_0=self.lon_ref)

        #if ('u' in self.center) & ('v' in self.center):
        #    self['x_r'] = self['x']-self.center['u']*self['dt']
        #    self['y_r'] = self['y']-self.center['v']*self['dt']


def plot_ring(self):
    import pylab
    pylab.quiver(self.data['x'], self.data['y'], self.input['u'],
            self.input['v'], color='r')
    pylab.quiver(self.data['xr'], self.data['yr'], self.data['ur'], self.data['vr'])
    pylab.plot(self.data['xc'], self.data['yc'], 'ro')
    pylab.figure()
    pylab.plot(self.data['r'], self.data['vtan'], 'r')
    pylab.plot(self.data['r'], self.data['vrad'], 'g')
    pylab.show()


# ============================================================================
def center_velocity_correction(dt, x, y, u, v, u_c, v_c):
    """Correct data due ring propagation velocity

    Input:
        -> u => Zonal velocity [m/s]
        -> v => Meridional velocity [m/s]
        -> u_c => Ring center zonal velocity [m/s]
        -> v_c => Ring center meridional velocity [m/s]

    Subtract the velocity field by the ring propagation and
      move the data position according to the ring propagation,
      so estimate what should be the sample if was done on one
      instant.
    """
    # Correct velocity field from center movement
    u_new = u-u_c
    v_new = v-v_c

    # Corrected positions due center movement
    x_new=x-dt*u_c
    y_new=y-dt*v_c

    return x_new, y_new, u_new, v_new

