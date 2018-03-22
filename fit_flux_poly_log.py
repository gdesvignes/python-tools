from pylab import *
import scipy.optimize
import numpy as np
import sys


fn = sys.argv[1]
# File should be freq (GHz), flux (mJy)
xdata, ydata = np.loadtxt(fn, usecols=(0,1), unpack=True)

# Create fake errors for now
yerr = np.ones(len(xdata))/10.

# Fit 
out = scipy.optimize.curve_fit(lambda t,a,b,c,d: a+b*np.log10(t)+c*np.power(np.log10(t),2) +d*np.power(np.log10(t),3),  xdata,  np.log10(ydata))

pfinal = out[0]
covar = out[1]
print "Best fit parameters for flux calibration following"
print "[a0 a1 a2 a3]"
# Flux in Jy for a frequency in GHz is:
# log10(S) = a_0 + a_1*log10(f) + a_2*(log10(f))^2 + ...
print pfinal

a0 = pfinal[0]
a1 = pfinal[1]
a2 = pfinal[2]
a3 = pfinal[3]

# For reference, 3C48 values:
# see http://www.vla.nrao.edu/astro/calib/manual/baars.html
#
#a0 = 1.31752
#a1 = -0.74090
#a2 = -0.16708
#a3 = 0.01525

# Plotting data
clf()
subplot(1, 1, 1)
errorbar(xdata, ydata, yerr=yerr, fmt='k.')  # Data
plot(xdata, 10**(a0+a1*np.log10(xdata)+a2*np.log10(xdata)*np.log10(xdata) + a3*np.power(np.log10(xdata), 3)) , 'r')     # Fit
title("Best fit coefficients:\n a0=%f a1=%f a2=%f a3=%f"%(a0, a1, a2, a3))
xlabel('Frequency (GHz)')
ylabel('Flux (Jy)')
show()
