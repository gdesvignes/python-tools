#!/usr/bin/env python

import struct, sys
import numpy as np
from ppgplot import *
from math import *
from optparse import OptionParser

full_usage = """
usage : pltasc.py [options] .asc files

  [-h, --help]        : Display this help
  [-r, --rot]         : Rotate the profile


"""
usage = "usage: %prog [options] files"

def main():
  parser = OptionParser(usage)
  parser.add_option("-x", "--xwin", action="store_true", dest="xwin",
                          default=False, help="Don't make a postscript plot, just use an X-window")

  parser.add_option("-m", "--m2", type="float", dest="m2", default=0.05,
                        help="Companion mass")
  parser.add_option("-s", "--sini", type="float", dest="sini", default=0.99,
                        help="Sin of the inclination of the orbit")
  parser.add_option("-o", "--omega", type="float", dest="omega", default=0.0,
                        help="Longitude of ascending node")
	
  (opts, args) = parser.parse_args()	

  if len(args)==0:
    print full_usage
    sys.exit(0)

  phase, resid, resid_e = np.loadtxt(args[0], usecols=(1,2,3), unpack=True)

  # Open PGPLOT
  if not opts.xwin:
    device='/xw'
  else:
    device='plot.ps/ps'
  pgopen(device)

  # 
  pgpap(0.0,0.7)


  T = 4.925e-6
  x = np.arange(0, 1., 0.001)
  phi = x * 2 * np.pi
  shapiro = - 2 * opts.m2 * T * np.log (1 - opts.sini * np.sin (phi + opts.omega*np.pi/180.))

  pgscf(2)
  pgsvp(0.1,0.9,0.1,0.9)
  pgswin(0.0, 1., -25 , 50)
  pgbox("BCNST", 0, 0, "BCNSTV", 0, 0)

  pgline(x, shapiro*1e6)
  pgpt(phase, resid*1e6, 17)
  pgerry(phase, resid*1e6 - resid_e, resid*1e6 + resid_e, 1.)

  pgsch(1.2)
  pgmtxt("B",3.0,0.5,0.5,"Rotational phase")
  pgmtxt("L",3.0,0.5,0.5,"Residual (\\gms)")
  #pgptxt(0.07,1.0,0.0,0.0,"(a)")

  pgend()

if __name__ == '__main__':
  main()
