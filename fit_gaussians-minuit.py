from ROOT import TH1D, TF1, TCanvas, TPad, TFile, TPaveLabel, TPaveText
from ROOT import gROOT, gStyle
from array import array

import psrchive
import sys
import time
import numpy as np
from optparse import OptionParser


full_usage = """
"""
usage = "usage: %prog [options]"

gStyle.SetOptFit()



parser = OptionParser(usage)

parser.add_option("-n", "--ncomp", type="int", dest="ncomp", default=1,
			help="Number of components")
parser.add_option("-t", "--template", type="string", dest="template",
			help="Template of gaussians")
parser.add_option("-s", "--save", action="store_true", dest="save",
			 default=False, help="Save gaussians to file")
parser.add_option("-r", "--rotate", action="store_true", dest="rotate",
                         default=False, help="Rotate profile by 180deg for display and fit")
parser.add_option("-c", "--bscr", type="int", dest="bscr", default=1,
                        help="Bscrunch by this factor (power of 2)")

(opts, args) = parser.parse_args()
if len(sys.argv) <= 1:
    print usage
    sys.exit()

if opts.template:
    gaussian_id, phase, heigth, width = np.loadtxt(opts.template, unpack=True, dtype=float, skiprows=3)  
    print phase, heigth, width
			

fn = args[0]
ar = psrchive.Archive_load(fn)

#ar.dedisperse()
if opts.rotate:
    ar.rotate_phase(.5)
ar.tscrunch()
ar.fscrunch()
if opts.bscr>1:
    print "Will scrunch the profile by ", opts.bscr, " bins"
    ar.bscrunch(opts.bscr)
ar.pscrunch()
ar.remove_baseline()

# Get I
#prof = ar.get_Profile(0,0,0)
nbin = ar.get_nbin()
#x = prof.get_amps()

data = ar.get_data()[0][0][0]

# NO scrunching for flux calibrated data

#x = data / data.max()
x = data
#for a,b in enumerate(x):
#    print a,b
#max_bin = prof.find_max_bin()
#max_bin = nbin/2
#lowbin = max_bin-30
#hibin = max_bin+30

# For interpulse
#lowbin = max_bin-40
#hibin = max_bin+40

#low_phase = 0.47
#hi_phase = 0.6

low_phase = 0.42
hi_phase = 0.6

lowbin = int(low_phase * nbin)
hibin = int(hi_phase * nbin)

gROOT.Reset()


# Canvas
c1 = TCanvas( 'c1', 'Histogram Drawing Options', 200, 10, 700, 900 )

pad1 = TPad( 'pad1', 'The pad with the function',  0.03, 0.62, 0.50, 0.92, 21 )
pad2 = TPad( 'pad2', 'The pad with the histogram', 0.51, 0.62, 0.98, 0.92, 21 )
pad3 = TPad( 'pad3', 'The pad with the histogram', 0.03, 0.02, 0.97, 0.57, 21 )
pad1.Draw()
pad2.Draw()
pad3.Draw()



npts = len(x)
h1 = TH1D( 'h1', 'Full pulse profile', npts, 0, npts)
h2 = TH1D( 'h2', 'Diff with template', hibin-lowbin, lowbin, hibin)
h3 = TH1D( 'h3', 'Fitted pulse profile', hibin-lowbin, lowbin, hibin)

for i in xrange(npts):
   h1.SetBinContent( i+1, x[i] )
   #h2.SetBinContent( i+1, x[i] )

for i,xx in enumerate(x[lowbin:hibin]):
    h2.SetBinContent( i+1, xx )
    h3.SetBinContent( i+1, xx )

par = array( 'd', 3*[0.] )
total = TF1( 'total', 'gaus',  lowbin, hibin )
# Set initial conditions
par[0] = h2.GetMaximum()    
par[1] = h2.GetMaximumBin() + lowbin
par[2] = 5
print "Initial conditions:", par[0], par[1], par[2]
total.SetParameters( par )
#total.FixParameter(0, par[0])
#total.FixParameter(1, par[1])
#total.FixParameter(2, par[2])
#g2 = TF1( 'g2', 'gaus',  98, 108 )
#g3 = TF1( 'g3', 'gaus', 110, 121 )

#total = TF1( 'total', 'gaus(0)+gaus(3)+gaus(6)', 85, 125 )
#total.SetLineColor( 2 )

# Draw histogram hpx in first pad with the default option.
pad1.cd()
h1.Draw()


pad3.cd()
h3.Draw()
#fit = h3.Fit( total, 'RS')
fit = h3.Fit( total, 'RS' , '',lowbin, hibin)



par1 = total.GetParameters()
#print par1[0], par1[1], par1[2]


icomp = 1

# residual plot
for i in range(h3.GetNbinsX()):
  h2.SetBinContent( i+1, h3.GetBinContent(i) - total.Eval(h3.GetBinCenter(i)) )
  #print i+1, h3.GetBinContent(i) - total.Eval(h3.GetBinCenter(i))
pad2.cd()
h2.Draw()


if opts.template:
    ncomp = len(gaussian_id)
    print "Ncomponents: ", ncomp

    # init the model to a single gaussian  
    model = 'gaus(0)'
    for i in range(ncomp-1):
	model += '+gaus(%d)'%((i+1)*3)
    print model

    par = array( 'd', 3*ncomp*[0.] )
    for i in range(ncomp):
	par[3*(i)] = heigth[i]
	par[3*(i)+1] = phase[i] /360. *  nbin
	par[3*(i)+2] = width[i]
    print par    

    total = TF1( 'total', model,  lowbin, hibin )
    total.SetParameters( par )


    # fit and display
    pad3.Clear()
    h3.Fit(total, 'ER')
    par1 = total.GetParameters()

    # residual plot
    pad2.Clear()
    for i in range(h3.GetNbinsX()):
	h2.SetBinContent( i+1, h3.GetBinContent(i) - total.Eval(h3.GetBinCenter(i)) )
    pad2.cd()
    h2.Draw()

    icomp = ncomp

#if 1:
else:

    # init the model to a single gaussian  
    model = 'gaus(0)'
    ncomp = opts.ncomp

    while icomp < ncomp:
      par = array( 'd', 3*(icomp+1)*[0.] )

      model += '+gaus(%d)'%(icomp*3)
      #print model
      total = TF1( 'total', model,  lowbin, hibin )

      # Set params
      print "Max residual bin: ", h2.GetMaximumBin()   + lowbin, "Value: ", h2.GetMaximum()
      for i in range(icomp):
	  par[3*i],par[3*i+1],par[3*i+2] = par1[3*i],par1[3*i+1],par1[3*i+2]
      par[3*(i+1)] = h2.GetMaximum()    
      par[3*(i+1)+1] = h2.GetMaximumBin() + lowbin
      par[3*(i+1)+2] = 4.
      total.SetParameters( par )


      # fit and display
      pad3.Clear()
      h3.Fit(total, 'RS')
      par1 = total.GetParameters()

      # residual plot
      pad2.Clear()
      for i in range(h3.GetNbinsX()):
	  h2.SetBinContent( i+1, h3.GetBinContent(i) - total.Eval(h3.GetBinCenter(i)) )
	  #print i+1, h3.GetBinContent(i) - total.Eval(h3.GetBinCenter(i))
      pad2.cd()
      h2.Draw()

      icomp +=1     
      #time.sleep(1)

if opts.rotate:
    offset = 180.
else:
    offset = 0.

par = total.GetParameters()

if opts.save:
    pfo = open(fn + '.gaussian', 'w')    
    pfo.write('#### Gaussian Template #####\n')
    pfo.write('MJD: %s\n'%ar.start_time().strtempo())
    pfo.write('Components: %d\n'%ncomp)
    for icomp in range(ncomp):
	pfo.write('%d %f %f %f\n'%(icomp, par[3*icomp+1]*360./float(ar.get_nbin()) - offset, par[3*icomp], par[3*icomp+2]*360./float(ar.get_nbin())))
    pfo.close()    

print "\nSummary of fit"
for icomp in range(ncomp):
    print par[3*icomp+1]
    print 'Component #%d %f %f %fi (%f)\n'%(icomp, par[3*icomp+1]*360/float(ar.get_nbin()) - offset, par[3*icomp], par[3*icomp+2]*360/float(ar.get_nbin()), par[3*icomp+1]/float(ar.get_nbin()) - offset/360.)
	

raw_input("Press enter to continue...")
