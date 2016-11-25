#! /usr/bin/env python

import os.path
import sys
import pylab

import numpy as np
import scipy.constants
import numpy.ma as ma


from ROOT import TMinuit
from ROOT import Double
from ROOT import Long
from array import array as arr


from optparse import OptionParser

import psrchive

full_usage = """
usage : rmsyn.py [options] 
  [-h, --help]        : Display this help
  [-r, --rm]          : Estimate of the RM


"""

usage = "usage: %prog [options] "

# Has to be global to be accessible by ROOT
freqs_hz=0
pv=0
sigq=0
sigu=0

threshold = 4

def getPV(PV, RM):
    lbda = scipy.constants.c / freqs_hz
    return PV * np.exp(-2 * 1j * lbda * lbda * RM)

def get_chi2(par):
    Qmodel = par[0]
    Umodel = par[1]
    RM = par[2]
    a = par[3]
    b = par[4]

    datmodel = getPV((Qmodel + 1j* Umodel), RM) - (a +1j * b)

    sigq_ma = ma.masked_where(sigq <= 0, sigq)
    sigu_ma = ma.masked_where(sigu <= 0, sigu)

    chiq = np.real(pv-datmodel) / sigq_ma
    chiu = np.imag(pv-datmodel) / sigu_ma

    #print "RM = ", R
    #print np.sum(sigq_ma)
    #print np.sum(sigu_ma)
    #print sigq, sigu
    #print chiq, chiu
    #print np.sum(chiq*chiq + chiu*chiu) / ( 2*pv.size - 3 - 1)
    #return np.sum(chiq*chiq + chiu*chiu) / ( 2*pv.size - 3 - 1 -2) # -3 params - 1
    return np.sum(chiq*chiq + chiu*chiu)# -3 params - 1




def fcn(npar, gin, f, par, iflag):
    #    global ncount

    chisq = get_chi2(par)

    f[0] = chisq



def getL(rhos, pv):
    RM = rhos.reshape(rhos.size,1)
    return np.sum(getPV(pv, RM), axis=1)

def get_pulse_range(ar):
    prof = ar.get_Profile(0,0,0)
    subint = ar.get_Integration(0)
    (b_mean, b_var) = subint.baseline_stats()
    rms = np.sqrt(b_var[0][0])

    dat = prof.get_amps()
    lowbin = prof.find_max_bin()
    hibin = prof.find_max_bin()
    while dat[lowbin-1] > threshold * rms:
	lowbin -= 1
    while dat[hibin+1] > threshold * rms:
	hibin += 1
    return lowbin, hibin

def main():
    global freqs_hz, pv, sigq, sigu
    parser = OptionParser(usage)

    parser.add_option("-r", "--rm", type="float", dest="rm", default=0.,
                            help="Estimate of the RM")
    parser.add_option("-n", "--norm", action="store_true", dest="norm",
			    default=False, help="Verbose mode")
    parser.add_option("-b", "--baseline", action="store_true", dest="baseline",
			    default=False, help="Fit for a baseline")
    parser.add_option("-f", "--nofit", action="store_true", dest="nofit",
			    default=False, help="Don't fit for the RM with Minuit")
    parser.add_option("-c", "--clean", action="store_true", dest="clean",
			    default=False, help="Clean the baseline")

    (opts, args) = parser.parse_args()

    #if opts.help:
    #    print full_usage

    arch = psrchive.Archive_load(sys.argv[1])
    arch.tscrunch()
    arch.bscrunch(2)
    #arch.fscrunch(2)
    arch.dedisperse()
    arch.remove_baseline()
    arch.convert_state('Stokes')
    arch.centre_max_bin()
    #arch.rotate(0.5)
    data = arch.get_data()

    MJD = arch.get_first_Integration().get_epoch().strtempo()

    w = arch.get_weights().sum(0)
    I = data[:,0,:,:].mean(0)
    Q = data[:,1,:,:].mean(0)
    U = data[:,2,:,:].mean(0)
    #V = data[:,3,:,:].mean(0)

    I = ma.masked_array(I)
    Q = ma.masked_array(Q)
    U = ma.masked_array(U)
    #V = ma.masked_array(V)

    I[w==0] = np.ma.masked
    Q[w==0] = np.ma.masked
    U[w==0] = np.ma.masked


    # Get Freqs
    nchan = arch.get_nchan()
    freqs = np.zeros(nchan)
    for i in range(nchan):
	freqs[i] = arch.get_Profile(0,0,i).get_centre_frequency() 

    ## Determine phase range
    arch.fscrunch()
    arch.tscrunch()
    arch.pscrunch()
    prof = arch.get_Profile(0,0,0)
    #print prof.get_amps()

    lowbin, hibin = get_pulse_range(arch)

    #lowbin, hibin = 458,531
    #lowbin, hibin = 20,120



    if opts.clean:
      for n in range(nchan):
	if n%20: continue
        baseline_idx = np.append(np.arange(0, lowbin), np.arange(hibin, arch.get_nbin()))
	baseline_dat = I[n][baseline_idx]
	print n, baseline_idx, baseline_dat
	poly = np.polyfit(baseline_idx, baseline_dat, 1)
	p = np.poly1d(poly)
	pylab.plot(I[n])
	I[n] = I[n] - p(np.arange(0, arch.get_nbin()))
	pylab.plot(I[n])
	pylab.show()

    
    for ii, wx in enumerate(w):
	I[ii] *= wx/np.max(w)
	Q[ii] *= wx/np.max(w)
	U[ii] *= wx/np.max(w)
	#V[ii] *= wx

    if opts.norm:
        print "Will normalize by the rms of the noise"
	scale = np.std(I, axis=1) / np.max(np.std(I, axis=1))
	for ii, s in enumerate(scale):
	    I[ii,:] /= s
	    Q[ii,:] /= s
	    U[ii,:] /= s
    #Q = Q / scale	
    #U = U / scale	
	
    freq_lo = freqs[0]
    freq_hi = freqs[-1]


    lowphase = lowbin / float(arch.get_nbin())
    hiphase = hibin / float(arch.get_nbin())

    #lowphase= 0.48
    #hiphase = 0.515
    #lowphase = 0.515
    #hiphase = 0.54
    print "Pulse phase window: ", lowphase, hiphase

    #pylab.plot(np.sum(I), axis)
    #pylab.show()

    # 2D plot
    f = pylab.figure() 

    pylab.subplot(321)
    #print np.sum(I,axis=0)
    #pylab.plot(np.sum(I, axis=0))
    pylab.plot(prof.get_amps())
    pylab.xlim([0, arch.get_nbin()])

    f.text( .5, 0.95, r'%s'%os.path.split(sys.argv[1])[1], horizontalalignment='center') 
    pylab.subplot(323)
    pylab.axis([0,1,freq_lo,freq_hi])
    pylab.imshow(I,extent=(0,1,freq_lo,freq_hi), origin='lower', aspect='auto', vmax=I.max()/1.5, vmin=I.min()/1.5)
    pylab.xlabel('Pulse phase')
    pylab.ylabel('Frequency (MHz)')

    # Compute errors of Q and U
    sigq = np.std(Q[:,lowbin:hibin], axis=1) / np.sqrt(hibin-lowbin) 
    sigu = np.std(U[:,lowbin:hibin], axis=1) / np.sqrt(hibin-lowbin)
    #sigq = np.std(Q[:,:], axis=1) 
    #sigu = np.std(U[:,:], axis=1) 

    pylab.plot([lowphase,lowphase], [freq_lo, freq_hi], 'r--')
    pylab.plot([hiphase,hiphase], [freq_lo, freq_hi], 'r--')

    freqs_hz = freqs * 1e6
    # Select phase range
    lowbin = int(lowphase * arch.get_nbin())
    hibin = int(hiphase * arch.get_nbin())
    dat = Q + 1j * U
    pv = dat[:,lowbin:hibin].mean(1) # Select phase range and average

    #rhos = np.arange(-90000, 90000, 10)
    rhos = np.arange(-2000, 2000, 1)
    res = np.absolute(getL(rhos, pv))

    # Plot I f(RM)
    pylab.subplot(324)
    pylab.plot(rhos, res, 'b-')
    pylab.xlabel('RM')
    pylab.ylabel('I')


    # Plot Q
    pylab.subplot(325)
    pylab.errorbar(freqs, np.real(pv), yerr=sigq, ls='None')
    pylab.plot(freqs, np.real(pv), 'bo')
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Q')


    # Plot U
    pylab.subplot(326)
    pylab.errorbar(freqs, np.imag(pv), yerr=sigu, ls='None')
    pylab.plot(freqs, np.imag(pv), 'bo')
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('U')

    # Should get the initial RM
    if opts.rm:
	initial_RM = -opts.rm
    else:
        initial_RM = -rhos[np.argmax(res)]
	print "Will use initial RM", initial_RM
    initial_L = getL(np.array([initial_RM]), pv) / (hibin-lowbin)


    if not opts.nofit:
	# Minuit part
	gMinuit = TMinuit(5)
	#gMinuit.SetErrorDef(1) # 1-Sigma Error
	gMinuit.SetErrorDef(4) # 2-Sigma Error
	gMinuit.SetFCN(fcn)
	arglist = arr('d', 2*[0.01])
	ierflg = Long(0)

	#arglist[0] = 1
	#gMinuit.mnexcm("SET ERR", arglist ,1,ierflg)


	# Set initial parameter values for fit
	vstart = arr( 'd', (np.real(initial_L), np.imag(initial_L), initial_RM) )
	#vstart = arr( 'd', (2*np.real(initial_L), 2*np.imag(initial_L), initial_RM) )
	# Set step size for fit
	step =   arr( 'd', (0.001, 0.001, 0.001) )

	# Define the parameters for the fit
	gMinuit.mnparm(0, "Qf",  vstart[0], step[0], 0,0,ierflg)
	gMinuit.mnparm(1, "Uf",  vstart[1], step[1], 0,0,ierflg)
	gMinuit.mnparm(2, "RM",  vstart[2], step[2], 0,0,ierflg)

	gMinuit.mnparm(3, "a",  0.0, step[2], 0,0,ierflg)
	gMinuit.mnparm(4, "b",  0.0, step[2], 0,0,ierflg)
	
	#gMinuit.FixParameter(2)

	if not opts.baseline:
	    gMinuit.FixParameter(3)
	    gMinuit.FixParameter(4)

	arglist[0] = 6000 # Number of calls to FCN before giving up. 
	arglist[1] = 0.1  # Tolerance
	gMinuit.mnexcm("MIGRAD", arglist ,2,ierflg)

	amin, edm, errdef = Double(0.18), Double(0.19), Double(1.0)
	nvpar, nparx, icstat = Long(1), Long(2), Long(3)
	gMinuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);

	gMinuit.mnmnos();
	gMinuit.mnprin(3,amin);

	finalQ, finalQ_err = Double(0), Double(0)
	finalU, finalU_err = Double(0), Double(0)
	finalRM, finalRM_err = Double(0), Double(0)
	A, A_err = Double(0), Double(0)
	B, B_err = Double(0), Double(0)
	print "\n\n"

	gMinuit.GetParameter(0, finalQ, finalQ_err)
	gMinuit.GetParameter(1, finalU, finalU_err)
	gMinuit.GetParameter(2, finalRM, finalRM_err)
	gMinuit.GetParameter(3, A, A_err)
	gMinuit.GetParameter(4, B, B_err)

	finalRM *= -1.

	print finalQ, finalU
	line1 =  "Final RM: %.1f"%(finalRM) + "(+/-) %.1f "%(finalRM_err) + " (2-sig)    "
	line2 = "Chi**2r = %.2f"%(get_chi2((finalQ,finalU,-finalRM,A,B)) / ( 2*pv.size - 3 - 1 -2))
	print line1
	print line2
        f.text( .5, 0.92, line1 + line2, horizontalalignment='center') 

	print MJD, np.mean(freqs), finalRM, finalRM_err, get_chi2((finalQ,finalU,finalRM,A,B)), 2*pv.count()

	# Plot best fit model
	finalRM *= -1.
	pylab.subplot(325)
	pylab.plot(freqs, np.real( -A+getPV(finalQ+1j*finalU, finalRM) ))

	pylab.subplot(326)
	pylab.plot(freqs, np.imag( -B+getPV(finalQ+1j*finalU, finalRM) ))

    pylab.show()

if __name__ == '__main__':
    main()

