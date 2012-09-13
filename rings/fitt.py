#!/usr/bin/env python
# -*- coding: Latin-1 -*-
# vim: tabstop=4 shiftwidth=4 expandtab

""" A module to handle the fittings for the rings' class
"""

import numpy import ma

from rings.utils import uv2nt

class v_circular(object):
    def __init__(self):
        # Escalas
        self.s = [1e4, 1e4, 1e-1, 1e-1]
        # Improve it. Consider the data input to better suggest initial values
        #  like x0 as the median of x, or eta0 as marximum/minimum eta.
        self.p0 = [-1, -1, .1, .1]

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
        e = 1./(2*n)*ma.sum( (vr/mag)**2 )
        # Maybe
        #e = 1./(2*n)*ma.sum( vr**2/mag )
        #e = 1./n*ma.sum( vr**2 )
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

