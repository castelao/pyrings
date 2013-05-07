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

        I'm not sure how would be the best way to include the
          penalty factor lambda. Might be a good idea to depend on
          the number of available observations. More observations,
          less succetible it is to the over fitting effect.

        I should probably do first a fit as a non-translating case
          and only then try to do a translating.
    """
    def __init__(self):
        # Escalas
        self.s = [1e4, 1e4, 1e-2, 1e-2]
        # xc|yc [10 km], uc|vc [1 cm], lambda

        # First guess for the parameters
        # Improve it. Consider the data input to better suggest initial values
        #   like x0 as the median of x, or eta0 as marximum/minimum eta.
        # Consider to migrate for named tuples or dictionary.
        self.p0 = [-1, -1, .01, .01]

        # I need to do a better job defining the lambda. It's working well,
        #   but I feel like it could be improved.
        self.lamb = 1e-2

    #def set_p0(self, p0=None):
    def set_p0(self, x, y, u, v):
        """

            Don't like to start with zeros. Improve this here.
        """
        #    if p0 == None:
        #        self.p0 = [-1, -1, .1, .1]
        self.p0[0] = np.median(x)/self.s[0] + 1e-2
        self.p0[1] = np.median(y)/self.s[1] + 1e-2
        #self.p0[2] = np.median(u)/self.s[2]
        #self.p0[3] = np.median(v)/self.s[3]

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
        # Maybe
        #e = 1./(2*n)*ma.sum( vr**2/mag )

        j = 1./(2*n)*ma.sum( (vr/mag)**2 ) + \
                self.lamb/2 * (s[2]*p[2]**2 + s[3]*p[3]**2)
        return j

class v_circular_nontranslating(object):
    def __init__(self):
        # Escalas
        self.s = [1e4, 1e4]

        self.p0 = [-1, -1]

    def set_p0(self, x, y, u, v):
        self.p0[0] = np.median(x)/self.s[0]
        self.p0[1] = np.median(y)/self.s[1]

    def cost(self, p, x, y, u, v):
        """
        """
        # I'm not sure this is the best way to count the number of elements.
        n = len(u.compressed().flatten())

        s = self.s
        vr, vt = uv2nt(x, y, u, v, x_c=s[0]*p[0], y_c=s[1]*p[1])
        mag = (u**2+v**2)**0.5
        j = 1./(2*n)*ma.sum( (vr/mag)**2 )
        return j

# =========================
# Moving carton_uv to here just to be able to run the Monte Carlo tests,
#   but in the future, create a new class for the Carton&Legras vortex.

# Moved from rings.py
#Saving this notes from a modified version of the paper.
#\note{Point the solid body inner core, hence CL would be a good model.}
#
#\citet{CartonLegras1994}
#
#\begin{displaymath}
#    v_{\theta CL} = \frac{1}{2} \omega_0 r \exp \left( - \left( \frac{r}{\delta} \right)^\alpha \right)
#\end{displaymath}
#\begin{displaymath}
#    %J = \frac{1}{2N} \sum_{n=1}^{N} \left( v_{\theta CL \ n} - v_{\theta n} \right)^2
#    J = \frac{1}{2N} \sum \left( v_{\theta CL} - v_{\theta} \right)^2
#    %J = \frac{1}{2N} \sum \left( v_{\theta CL}^2 + v_r^2 - |V|^2 \right)^2
#\end{displaymath}

def carton_uv(x, y,omega0, delta, alpha):
    """ Carton's model for azimuthal velocity
    """
    r = ((x)**2+(y)**2)**0.5
    mag = 0.5*omega0*(((x)**2+(y)**2)**0.5)*np.exp(-((((x)**2+(y)**2)**0.5)/delta)**alpha)
    u,v = nt2uv(x,y,0,mag)
    return u, v

def carton_scales(omega0, delta, alpha):
    r = np.arange(0, 400e3, 5e3)
    v = 0.5*omega0*(r)*np.exp(-(r/delta)**alpha)
    Vmax = v.max()
    Rmax = r[v.argmax()]
    return Vmax, Rmax
