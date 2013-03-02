#!/usr/bin/env python
# -*- coding: Latin-1 -*-
# vim: tabstop=4 shiftwidth=4 expandtab

""" A module to handle the fittings for the rings' class
"""

import numpy as np
from numpy import ma

from rings.utils import uv2nt, nt2uv

class v_circular(object):
    """

        Recomendation from the paper reviewer! Include the
          penalty function for the translation speed. Old eq. 7,
          eq. 11 on the new version.
          J = \frac{1}{2N} \sum_{n=1}^{N}
              \frac{{v_{r}^2}_n}{|V_n|}
            + \frac{\lambda}{2} (u_{cn} + v_{cn})^2.
    """
    def __init__(self):
        # Escalas
        self.s = [1e4, 1e4, 1e-1, 1e-1]
        # Improve it. Consider the data input to better suggest initial values
        #   like x0 as the median of x, or eta0 as marximum/minimum eta.
        # Consider to migrate for named tuples or dictionary.
        self.p0 = [-1, -1, .1, .1, 0]

    #def set_p0(self, p0=None):
    #    if p0 == None:
    #        self.p0 = [-1, -1, .1, .1]

    def cost(self, p, dt, x, y, u, v):
        """
        """
        # I'm not sure this is the best way to count the number of elements.
        n = len(u.compressed().flatten())

        s = self.s
        u_new = u - s[2]*p[2]
        v_new = v - s[3]*p[3]
        x_new = x - dt*s[2]*p[2]
        y_new = y - dt*s[2]*p[3]
        vr, vt = uv2nt(x_new, y_new, u_new, v_new, x_c=s[0]*p[0], y_c=s[1]*p[1])
        mag = (u**2+v**2)**0.5
        #vt = -u*ma.sin(ma.arctan2(y-p[1], x-p[0])) + \
        #      v*ma.cos(ma.arctan2(y-p[1], x-p[0]))
        #e = 1./(2*n)*ma.sum( (vr/mag)**2 )
        # Penalized version
        # Maybe
        #e = 1./(2*n)*ma.sum( vr**2/mag )
        #e = 1./n*ma.sum( vr**2 )
        e = 1./(2*n)*ma.sum( (vr/mag)**2 ) + \
            p[4]/2 * (s[2] + s[3])**0.5
        return e

class v_circular_nontranslating(object):
    def __init__(self):
        # Escalas
        self.s = [1e4, 1e4]

    def cost(self, p, x, y, u, v):
        """
        """
        # I'm not sure this is the best way to count the number of elements.
        n = len(u.compressed().flatten())

        s = self.s
        vr, vt = uv2nt(x, y, u, v, x_c=s[0]*p[0], y_c=s[1]*p[1])
        mag = (u**2+v**2)**0.5
        #vt = -u*ma.sin(ma.arctan2(y-p[1], x-p[0])) + \
        #      v*ma.cos(ma.arctan2(y-p[1], x-p[0]))
        e = 1./n*ma.sum( (vr/mag)**2 )
        #e = 1./n*ma.sum( vr**2 )
        return e

# =========================
# Moving carton_uv to here just to be able to run the Monte Carlo tests,
#   but in the future, create a new class for the Carton&Legras vortex.
def carton_uv(x, y,omega0, delta, alpha):
    """ Carton's model for azimuthal velocity
    """
    r = ((x)**2+(y)**2)**0.5
    mag = 0.5*omega0*(((x)**2+(y)**2)**0.5)*np.exp(-((((x)**2+(y)**2)**0.5)/delta)**alpha)
    u,v = nt2uv(x,y,0,mag)
    return u, v

