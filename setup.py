from __future__ import absolute_import
from __future__ import print_function
from distutils.core import setup, Extension
import os, numpy

""" Setup script for the sla python extension"""

library_dirs = []
include_dirs = []

# need to direct to where includes and  libraries are
if 'TRM_SOFTWARE' in os.environ:
    library_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'lib'))
    include_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'include'))
else:
    print("Environment variable TRM_SOFTWARE pointing to location of shareable libraries and includes not defined!", file=sys.stderr)

include_dirs.append(numpy.get_include())

sla = Extension('trm.sla._sla',
                define_macros   = [('MAJOR_VERSION', '0'),
                                   ('MINOR_VERSION', '1')],
                undef_macros    = ['USE_NUMARRAY'],
                include_dirs    = include_dirs,
                library_dirs    = library_dirs,
                runtime_library_dirs = library_dirs,
                libraries       = ['csla'],
                sources         = [os.path.join('trm', 'sla', 'sla.cc')])

setup(name='trm.sla',
      version='0.1',
      packages = ['trm', 'trm.sla'],
      ext_modules=[sla],

      author='Tom Marsh',
      author_email='t.r.marsh@warwick.ac.uk',
      url='http://www.astro.warwick.ac.uk/',
      description='Python interface to slalib.',
      long_description="""
sla provides an interface to slalib, a set of routines for celestial coordinate and time 
manipulations written by Pat Wallace. The Python sla interface uses this to provide a few 
convenient routines for Python scripts. These provide light travel time corrections good to 
about 50 microsecs amongst other things. It is not a complete interface to all sla routines.
""",
      )

