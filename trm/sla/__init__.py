#!/usr/bin/env python

"""
an interface to slalib

This module defines a few routines to access Pat Wallace's sla through 
Python.

Functions
=========

amass       -- calculates observational parameters given position and time
cldj        -- compute MJD from a date
djcl        -- compute date from an MJD
dtt         -- gives TT-UTC
sun         -- computes Sun's position on the sky.
sun_at_elev -- works out when the Sun crosses a given elevation
utc2tdb     -- compute tdb, heliocentric corrections etc

"""
import sys
sys.path.append('.')
from _sla import *

import exceptions
import re

# Exception class
class SlaError(exceptions.Exception):
    """For throwing exceptions from the sla module"""
    def __init__(self, value):
        self.value = value
            
    def __str__(self):
        return repr(self.value)

def sun_at_elev(longitude, latitude, height, utc1, utc2, elev, wave=0.55, rh=0.2, acc=1.e-5, fast=False):
    """
    utc = sun_at_elev(utc1, utc2, elev, acc=1.e-5)

    Works out when the Sun is at elevation elev. The user supplies times before
    and after the time of interest and this routine narrows down the range. It is 
    up to the user to bracket the times closely enough that there is no ambiguity
    about which crossing of the critical elevation is returned. For sunrise and sunset,
    the standard definition is when the limb of the Sun touches the horizon so you 
    should set elev = -0.25

    longitude -- longitude of oberving site, East positive, degrees
    latitude  -- latitude of oberving site, degrees
    height    -- height of observing site, metres
    utc1      -- mjd before time of interest
    utc2      -- mjd after time of interest. Elevation of Sun must be on opposite
                 sides of elev at utc1 and utc2 or an SlaError is raised
    elev      -- elevation in degrees
    wave      -- wavelength of observation, microns
    rh        -- relative humidity, 0 to 1
    acc       -- accuracy of final answer in days.

    The utc is returned as a decimal MJD
    """

    (az1, el1, ref1, ra1, dec1) = sun(utc1,longitude,latitude,height,wave,rh,fast)
    (az2, el2, ref2, ra2, dec2) = sun(utc2,longitude,latitude,height,wave,rh,fast)

    if (el1 > elev and el2 > elev) or (el1 < elev and el2 < elev):
        raise SlaError('Initial times do not bracket the critical elevation, el1, el2: ' 
                       + str(el1) + ', ' + str(el2))

    
    while utc2 - utc1 > acc:
        utc = utc1 + (elev-el1)/(el2-el1)*(utc2-utc1)
        (az, el, ref, ra, dec) = sun(utc,longitude,latitude,height,wave,rh,fast)
        if (el > elev and el1 > elev) or (el < elev and el1 < elev):
            utc1 = utc
        else:
            utc2 = utc
    return utc

