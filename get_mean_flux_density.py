import psrchive
import sys
import pylab as plt
import numpy as np
import uncertainties as u
from uncertainties import unumpy as up

ar = psrchive.Archive_load(sys.argv[1])

minb = 970
maxb = 1100
minb = 0
maxb = 2047

minb = 880
maxb = 1190

try:
    minb = int(sys.argv[2])
    maxb = int(sys.argv[3])
except:
    print "Using default bin range value"
    pass

ar.dedisperse()
ar.tscrunch()
ar.fscrunch()
#ar.bscrunch(2)
ar.remove_baseline()
ar.convert_state('Stokes')
ar.centre_max_bin()
#ar.rotate(.5)
data = ar.get_data()
nbin = ar.get_nbin()

mjd = ar.get_first_Integration().get_epoch().strtempo()


mean = np.average(np.append(data[0][0][0][0:minb],data[0][0][0][maxb:nbin-1])) 
std = np.std(np.append(data[0][0][0][0:minb],data[0][0][0][maxb:nbin-1])) 
mean = 0

I = data[0][0][0] - mean
I2 = data[0][0][0][minb:maxb] - mean
#print np.arange(minb,maxb), I2
#print len(I)
#print len(np.arange(minb,maxb)), len(I2)

plt.plot(I)
plt.plot(np.arange(minb,maxb), I2, color='red')
plt.ylabel("Flux (mJy)")
plt.xlabel("Bin")

# std *= 1.1 # Assume a 10% error on flux calibration itself


flux = up.uarray(I2, np.ones(I2.shape)*std*2)
sum_flux = np.sum(flux)

print "MJD:", mjd
print "Mean flux at %.1f MHz:"%(ar.get_centre_frequency()), sum_flux/ar.get_nbin(), "mJy"
plt.show()

