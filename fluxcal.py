#! /usr/bin/env python
import pyfits
import numpy as np
import psrchive as p


class Fluxcal(object):
    def __init__(self, filename):

        self.filename = filename
        self.f = pyfits.open(self.filename)

        self.get_fluxcal_params()

    def get_fluxcal_params(self):
        t = self.f['FLUX_CAL']
        nchan = t.header['NCHAN']
        nrcvr = t.header['NRCVR']
        self.mjd = float(t.header['EPOCH'])
        self.data_cal = t.data.field('S_CAL').reshape((nrcvr,nchan))
        self.err_cal = t.data.field('S_CALERR').reshape((nrcvr,nchan))
        self.data_sys = t.data.field('S_SYS').reshape((nrcvr,nchan))
        self.err_sys = t.data.field('S_SYSERR').reshape((nrcvr,nchan))
        self.wts = t.data.field('DAT_WTS').reshape((nchan,))
        self.freq = t.data.field('DAT_FREQ').reshape((nchan,))
        self.nchan = nchan
        self.nrcvr = nrcvr
        
        ar = p.Archive_load(self.filename)
        self.obsfreq = ar.get_centre_frequency()

    def get_scaling(self, calfilename):
        """Open a Stokes-I only CAL file to compute the scaling factors from counts to Jy"""
        self.calfile = calfilename
        
        #self.f = pyfits.open(self.calfile)
        cal = p.Archive_load(self.calfile)
        cal.tscrunch()

        #Check if CALfile matches the fluxcalibrator
        npol = cal.get_npol()
        nchan = cal.get_nchan()
        obsfreq = cal.get_centre_frequency()

        if npol!=1 or npol!=self.nrcvr:
            print "Error: Npol=%d in fluxcal does not match Npol=%d in CAL:"%(self.nrcvr, npol)
            sys.exit(-1)

        if self.nchan!=nchan:
            print "Error: Nchan=%d in fluxcal does not match Nchan=%d in CAL:"%(self.nchan, nchan)
            sys.exit(-1)

        # Check if same centre frequency
        if self.obsfreq != obsfreq:
            print "Error: obsfreq=%f in fluxcal does not match obsfreq=%f in CAL:"%(self.obsfreq, obsfreq)
            sys.exit(-1)

        i = cal.get_Integration(0)
        levels = i.cal_levels()
        self.cal_high = levels[0][0]
        self.cal_low = levels[1][0]
        self.cal_high_sig = levels[2][0]
        self.cal_low_sig = levels[3][0]
        self.dcal = self.cal_high - self.cal_low
        self.scale = self.data_cal[0] / self.dcal 

    def fcal_archive(self, ar_filename):
        """Flux calibrate the archive 'ar_filename'. The method 'get_scaling' needs to have been called beforehand"""


        ar = p.Archive_load(ar_filename)
        # Check Npol
        if ar.get_npol()!=1 or ar.get_npol()!=self.nrcvr:
            print "Error: Npol=%d in fluxcal does not match Npol=%d in archive:"%(self.nrcvr, ar.get_npol())
            sys.exit(-1)

        # Check Nchan
        if self.nchan!=ar.get_nchan():
            print "Error: Nchan=%d in fluxcal does not match Nchan=%d in archive:"%(self.nchan, ar.get_nchan())
            sys.exit(-1)

        # Check if same centre frequency
        if self.obsfreq != ar.get_centre_frequency():
            print "Error: obsfreq=%f in fluxcal does not match obsfreq=%f in archive:"%(self.obsfreq, ar.get_centre_frequency())
            sys.exit(-1)

        # Rescale the archive to Janksy
        for isub in range(ar.get_nsubint()):
            i = ar.get_Integration(isub)
            for ichan in range(ar.get_nchan()):
                prof = i.get_Profile(0, ichan)
                pdata = prof.get_amps()
                pdata *= self.scale[ichan]


        ar.set_scale("Jansky")
        ar.unload(ar_filename+ ".calib")
            
        
        

    def get_ssys(self, pol, reffreq):
        # Select the closest frequency to the one asked
	idx = (np.abs(self.freqs - reffreq)).argmin()
	freq = self.freqs[idx]

	if (np.abs(freq - reffreq)) > 50.:
	    return None

        return freq, self.ssys[pol]
        return freq, self.ssys.transpose()[pol][idx], self.ssys_err.transpose()[pol][idx]

    def get_scal(self, pol, reffreq):
        # Select the closest frequency to the one asked
	idx = (np.abs(self.freqs - reffreq)).argmin()
	freq = self.freqs[idx]

	if (np.abs(freq - reffreq)) > 50.:
	    return None

        return freq, self.scal.transpose()[pol][idx], self.scal_err.transpose()[pol][idx]



