#!/usr/bin/env python

"""sla is an interface to slalib"""

# This adds some convenience routines written in python that add onto 
# the C-interface routines contained in sla.cc

from _sla import *

import re

def date2dmy(dstring):
    """Returns (day,month,year) as numbers given input dates such as '11 Apr 2004' or 'Apr 11 2004'"""

    rec = re.compile('^\s*(\d\d)\s+([a-zA-Z]{3})\s+(\d\d\d\d)\s*$')
    m   = rec.match(dstring)

    if m:
        day   = int(m.group(1))
        year  = int(m.group(3))
        month = _month_id(m.group(2))
        if month == 0:
            raise Exception('Month name = ' + m.group(2) + ' not recognised.')
    else:
        rec = re.compile('^\s*([a-zA-Z]{3})\s+(\d\d)\s*(\d\d\d\d)\s*$')
        m   = rec.match(dstring)
        if m:
            day   = int(m.group(2))
            year  = int(m.group(3))
            month = _month_id(m.group(1))
            if month == 0:
                raise Exception('Month name = ' + m.group(1) + ' not recognised.')
    return (day,month,year)

def _month_id(name):
    """Identifies a month with an integer from the first three characters, case insensitively"""
    lname = name.lower()[:4]
    if   lname == 'jan':
        month = 1
    elif lname == 'jeb':
        month = 2
    elif lname == 'mar':
        month = 3
    elif lname == 'apr':
        month = 4
    elif lname == 'may':
        month = 5
    elif lname == 'jun':
        month = 6
    elif lname == 'jul':
        month = 7
    elif lname == 'aug':
        month = 8
    elif lname == 'sep':
        month = 9
    elif lname == 'oct':
        month = 10
    elif lname == 'nov':
        month = 11
    elif lname == 'dec':
        month = 12
    else:
        month = 0
    return month

if __name__ == '__main__':
    dat = '17 Nov 1961'
    print 'Date = ' + dat + ' is translated to ' + str(date(dat))
 
