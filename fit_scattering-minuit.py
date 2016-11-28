from ROOT import TH1D, TF1, TCanvas, TPad, TFile, TPaveLabel, TPaveText
from ROOT import gStyle,gROOT 
from ROOT import Math
from array import array

import psrchive
import sys
import time
import numpy as np
from optparse import OptionParser
import math

full_usage = """
"""
usage = "usage: %prog [options]"

def gaussian(x ,par):
    return  par[0] * math.exp (-1/2 *((x[0]-par[1])/par[2])**2)

def scat(x ,par):
    """
    A = par[0]
    sigma = par[1]
    tau = par[2]
    b = par[3]

    """
    t0=par[0]
    A=par[1]
    sigma=par[2]
    tau=par[3]
    b=par[4]
    dt=(x[0]-t0)

    t1 = tau * dt
    t3 = sigma * sigma
    t5 = tau * tau
    t9 = math.exp(-(0.4e1 * t1 - t3) / t5 / 0.4e1)
    t10 = 1.77245385090552;
    t19 = Math.erf((0.2e1 * t1 - t3) / sigma / tau / 0.2e1)
    return(A*t9 * t10 * sigma * (t19 + 0.1e1) / 0.2e1+b)


    #return  math.sqrt(math.pi) *  par[0] * Math.exp ( -x/par[2] )


parser = OptionParser(usage)

parser.add_option("-n", "--ncomp", type="int", dest="ncomp", default=1,
			help="Number of components")
parser.add_option("-t", "--template", type="string", dest="template",
			help="Template of gaussians")

(opts, args) = parser.parse_args()
if len(sys.argv) <= 1:
    print usage
    sys.exit()


fn = args[0]
ar = psrchive.Archive_load(fn)

#ar.rotate(0.5)

ar.dedisperse()
ar.tscrunch()
ar.fscrunch()
ar.pscrunch()
ar.remove_baseline()
try:
    ar.centre()
except:
    pass

# Get I
prof = ar.get_Profile(0,0,0)
x = prof.get_amps()


# Get stats
subint = ar.get_Integration(0)
(b_mean, b_var) = subint.baseline_stats()
rms = math.sqrt(b_var[0][0])

# compute
period = subint.get_folding_period()
max_phase = prof.find_max_phase()
max_phase = max_phase * period

maxval = prof.max() / rms
nbin = ar.get_nbin()

print b_var[0][0], rms, max_phase

x =  x / rms

gROOT.Reset()
gStyle.SetOptFit()



# Canvas
c1 = TCanvas( 'c1', 'Scattering fitting', 200, 10, 900, 500 )



npts = len(x)
h1 = TH1D( 'h1', 'Pulse profile', npts, 0, period)


for i in xrange(npts):
   h1.SetBinContent( i+1, x[i] )

## Function to fit
g1 = TF1( 'g1', scat,   max_phase-0.25, max_phase+0.4 , 5)

par = array( 'd', 5*[0.] )
par[0] = max_phase
par[1] = maxval
par[2] = .1
par[3] = .1
par[4] = 0.0001
print par

g1.SetParName(0, 't0')
g1.SetParName(1, 'A')
g1.SetParName(2, 'sigma')
g1.SetParName(3, 'tau')
g1.SetParName(4, 'b')

# Set initial parameters
g1.SetParameters( par )

# Do Fit
fit = h1.Fit( g1, 'WRE' )

h1.Draw()

raw_input("Press enter to continue...")


		     
