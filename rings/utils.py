#!/usr/bin/env python
# -*- coding: Latin-1 -*-
# vim: tabstop=4 shiftwidth=4 expandtab

""" Util functions for the module rings
"""

import ConfigParser
import re
from datetime import timedelta

import numpy
from numpy import ma

# ============================================================================
def lonlat2xy(lon,lat,lon_0=None,lat_0=None):
    """ Convert pairs of (Lat,Lon) into (x,y)

        Input:
                      Lon [deg]
          Lat [deg]
          Lon_0 [deg] => Lon of the origin of the cartesian system
          Lat_0 [deg] => Lat of the origin of the cartesian system
        Output:
                      x [m]
          y [m]

        The projection is deformed as get away from the center. Since the
          Latitudes don't deform, the y is estimated first, then for each
          point is estimated the distante to the meridian of reference
          (Lon_0) considering the Latitude of the measurement.
    """
    if (lat_0==None) or (lon_0==None):
        lat_0=numpy.median(lat)
        lon_0=numpy.median(lon)
    from fluid.common.distance import distance
    y=distance(lat,0,lat_0,0)
    y[lat<lat_0]=-1*y[lat<lat_0]
    x=distance(lat,lon,lat,lon_0)
    x[lon<lon_0]=-1*x[lon<lon_0]
    return x,y

# ============================================================================
def xy2lonlat(x,y,lon_0,lat_0):
    """ Think how to improve it.1
    """
    DEG2RAD = (2*numpy.pi/360)
    RAD2DEG = 1/DEG2RAD
    DEG2MIN = 60.
    DEG2NM  = 60.
    NM2M   = 1852.    # Defined in Pond & Pickard p303.
    lat = lat_0+y/(DEG2NM*NM2M)
    fac = numpy.cos((lat+lat_0)/2.*numpy.pi/180)
    lon = lon_0+x/(DEG2NM*NM2M)/fac
    return lon, lat

# ============================================================================
def uv2nt(x,y,u,v,x_c=0,y_c=0):
    """Convert orthogonal velocity components to normal and tangent ones

    Based on x_c and y_c, with default values (0,0) convert u and v to
    normal (n) and tangent(t) velocities.
    Input:
        - u =>
        - v =>
        - x =>
        - y =>
        - x_c =>
        - y_c =>
    Output:
        - n =>
        - t =>
    """
    dx = x-x_c
    dy = y-y_c
    #
    alpha = numpy.arctan2(dy,dx) # Angle between x-axis and a rotated x passing throught the velocity point
    #
    n = u*numpy.cos(alpha)+v*numpy.sin(alpha)
    t = v*numpy.cos(alpha)-u*numpy.sin(alpha)
    #
    return [n,t]

# ============================================================================
def nt2uv(x,y,n,t,x_c=0,y_c=0):
    """Convert normal and tangent velocities to Orthogonal ones

    Convert normal(n) and tangent(t) velocities to orthogonal ones,
    u and v, based on x,y of each (n,t) velocity and the circular center
    x_c and y_c
    Input:
        - n =>
        - t =>
        - x =>
        - y =>
        - x_c =>
        - y_c =>
    Output:
        - u =>
        - v =>
    """

    dx = x-x_c
    dy = y-y_c

    alpha = numpy.arctan2(dy,dx) # Angle between x-axis and a rotated x passing throught the velocity point

    u = n*numpy.cos(alpha) - t*numpy.sin(alpha);
    v = n*numpy.sin(alpha) + t*numpy.cos(alpha);

    return [u,v]


# ============================================================================
# ==== Logging System
# ============================================================================

# This basic logger should be out of here. At most, in a miscellaneous file.

def basic_logger(logname="Generic", logfile=None):
    """ Create a logger

        !!!ATENTION!!!
        In the future change it to use a config file, which would give 
          far away more flexibility.
    """
    import logging
    import logging.handlers

    # Creating another log level
    logging.VERBOSE = logging.DEBUG - 1
    logging.addLevelName(logging.VERBOSE, 'VERBOSE')

    #create logger
    logger = logging.getLogger(logname)
    #logger.setLevel(logging.DEBUG)
    logger.setLevel(logging.VERBOSE)

    #create console handler and set level to debug
    ch = logging.StreamHandler()
    #ch.setLevel(logging.WARN)
    ch.setLevel(logging.INFO)
    #create formatter
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s -  %(message)s")
    #add formatter to ch
    ch.setFormatter(formatter)
    #add ch to logger
    logger.addHandler(ch)

    if logfile is not None:
        #create a rotate file handler
        fh = logging.handlers.RotatingFileHandler(logfile,
               mode='a', maxBytes=1000000, backupCount=10)
        fh.setLevel(logging.DEBUG)
        #ch.setLevel(logging.WARN)
        #create formatter
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s -  %(message)s")
        #add formatter to ch
        fh.setFormatter(formatter)
        #add ch to logger
        logger.addHandler(fh)

    #"application" code
    #logger.debug("debug message")
    #logger.info("info message")
    #logger.warn("warn message")
    #logger.error("error message")
    #logger.critical("critical message")

    return logger


# ============================================================================
# ==== Auxiliary functions
# ============================================================================

def cfg2dict(cfg_file):
    """Extract all variables from a ConfigParser file and create a dict[section][option]

        !!ATENTION!!, migrate it to YAML.

    """
    cfg = ConfigParser.SafeConfigParser()
    try:
        cfg.readfp(open(cfg_file))
    except:
        print "Couldn't open the config file: %s" % cfg_file
        print "Improve this error handle"
        return
    parameters = {}
    for section in cfg.sections():
        parameters[section] = {}
        for option in cfg.options(section):
            par = cfg.get(section, option)
            if option == 'w0':
                option = 'W0'
            if re.match("^\ *\[(.*)\]\ *$",par):
                l = re.search("^\ *\[(.*)\]\ *$",par).group(1)
                tmp = []
                for i in re.split(',',l):
                    if re.match("'.*'",i):
                        tmp.append(re.search("'(.*)'",i).group(1))
                    elif re.match("^-?\d+$",i):
                        tmp.append(int(i))
                    elif re.match("^-?\d+\.\d*$",i):
                        tmp.append(float(i))
                parameters[section][option] = tmp
                del(l,tmp,i)
            #if re.search("\[('(\w+)')?((,\ ?'(\w+)')?)+\]",par):
            #    parameters[section][option] = re.findall("'(\w+)'",par)
            #elif re.search("^\[(\[(\w+)(,\ ?\w+)*\](,\ )?)+\]$",par):
            #    parameters[section][option] = [] 
            #    for x in re.findall("\[\w+(?:,\ ?\w+)*\]",par):
            #        parameters[section][option]
            elif re.search("^True$",par):
                parameters[section][option] = True
            elif re.search("^False$",par):
                parameters[section][option] = False
            elif re.search("^-?\d+$",par):
                parameters[section][option] = int(par)
            elif re.search("^-?\d+(\.\d+)?(e-?\d+)?$", par):
                parameters[section][option] = float(par)
            elif re.match('^(?P<value>\d+)(?P<unit>days|hours|minutes|seconds)$',par):
                tmp = re.match('^(?P<value>\d+)(?P<unit>days|hours|minutes|seconds)$',par)
                if tmp.group('unit')=='days':
                    parameters[section][option] = timedelta(days=float(tmp.group('value')))
            else:
                parameters[section][option] = cfg.get(section,option)
    ks = parameters.keys()
    for k in ks:
        if re.search('\w+\.(?:\w+)+',k):
            kk=re.split('\.',k)
            # Stupid way not generic
            if kk[0] not in parameters:
                parameters[kk[0]] = {}
            if kk[1] not in parameters[kk[0]]:
                parameters[kk[0]][kk[1]] = {}
            if len(kk)==2:
                for kkk in parameters[k].keys():
                    parameters[kk[0]][kk[1]][kkk]=parameters[k][kkk]
            elif len(kk)==3:
                parameters[kk[0]][kk[1]][kk[2]]=parameters[k]
            del(parameters[k])
    # ------- Default values
    # I need to think about an efficient way to load default values.
    # ---- Eddy
    if 'eddy' not in parameters:
        parameters['eddy'] = {}
    if 'force_reprocess' not in parameters['eddy']:
        parameters['eddy']['force_reprocess'] = False
    # ---- Eddy Tracking system
    if 'eddy_track' not in parameters:
        parameters['eddy_track'] = {}
    if 'force_reprocess' not in parameters['eddy_track']:
        parameters['eddy_track']['force_reprocess'] = False
    return parameters


# ============================================================================
def neighbor_cluster(ind, verbose=False):
    """ Define neighbor patches

        The input is a boolean array with True on
          the values of interest. The output is an
          integer array, with one number for each
          group. I.e. all gridpoints with 10 are
          part of the same patch, the patch #10.

       ATENTION:
        Think about it. Maybe I should consider as neighbors only the
          NSEW, and don't include the diagonals.
        Maybe I should consider a jumping cell, i.e. if a neighbor is W<0
          and W<|W0|, and the next one is OK, should take all them as one group. No,
          this case should be dealed with filters. Smooth the data before group it,
          so I have control on the scales that I'm manipulating.
    """
    output = ma.masked_all(ind.shape, dtype='i')
    I,J = numpy.nonzero(ind)
    g = -1
    while I.any():
        if verbose:
            print "Missing %s gridpoints to evaluate" % (I.shape[0])
        g += 1
        gi, gj = [I[0]], [J[0]]
        I, J = I[1:], J[1:]
        while len(gi)>0:
            i, j = gi[0], gj[0]
            ind_around = numpy.nonzero(((I==i-1) & (J==j)) |
              ((I==i+1) & (J==j)) | ((I==i) & (J==j-1)) |
              ((I==i) & (J==j+1)))[0]
            for n in ind_around:
                gi.append(I[n])
                gj.append(J[n])
            I = numpy.delete(I, ind_around)
            J = numpy.delete(J, ind_around)
            output[gi[0],gj[0]] = g
            gi, gj = gi[1:], gj[1:]
    return output

# ============================================================================
# ============================================================================


def prepare_masked_array(data, float_fillvalue=1e20, int_fillvalue=-999999):
    """ Prepare masked array to be saved.

        Force data os masked elements to be a certain value

            # Be sure to set all masked values to 1e+20, so the masked 
            # array can be recreated once read.
            # This happens because the masked_all just put a bogus data in .data, so
            #   who understands the masked array goes fine, otherwise.
    """
    for k in data.keys():
        if type(data[k]) == numpy.ma.core.MaskedArray:
            if data[k].dtype in ('>f4', 'float32', '>f8', 'float64'):
                if data[k].mask.any():
                    data[k].data[data[k].mask == True] = float_fillvalue
                data[k].set_fill_value(float_fillvalue)
            elif data[k].dtype in ('int32',):
                if data[k].mask.any():
                    data[k].data[data[k].mask == True] = int_fillvalue
                data[k].set_fill_value(int_fillvalue)
    return data
