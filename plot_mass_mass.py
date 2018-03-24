from ppgplot import *
import psr_utils as p
import numpy as np
import parfile as pf
import sys

pb = 0.16599304948208
x = 1.41995201704099
ecc = 0.0852970911069719753
omega = 2.15
omega_e = 0.08


par = pf.Parfile(sys.argv[1])

try:
    pb = par.PB
    x = par.A1
except:
    raise "Error: no binary parameters in parfile"

mp = np.arange(0.0001, 2.51, 0.01)

pgopen("?")
pgscf(2)
pgpap(0., 1.0)
pgsvp(0.1, 0.9, 0.1, 0.9)
pgswin(0., 2.5, 0., 2.5)

pgslw(2)

# Inclination angle
try:
    inc = np.arcsin(par.SINI)
    inc_e = par.INC_ERR
    mcs1 = np.array([])
    mcs2 = np.array([])
    for m in mp:
        mc1 = p.companion_mass(pb, x, inc=inc+inc_e, mpsr=m)
        mc2 = p.companion_mass(pb, x, inc=inc-inc_e, mpsr=m)
        mcs1 = np.append(mcs1, mc1)
        mcs2 = np.append(mcs2, mc2)    
    pgsls(2);pgsci(2)
    pgline(mp, mcs1)
    pgline(mp, mcs2)
    pgsls(1);pgsci(1)
except:
    print "No inclination constraint"

# Omdot
try:
    omdot = par.OMDOT
    omdot_e = par.OMDOT_ERR
    mc = np.arange(0.0001, 2.51, 0.01)
    mps1 = p.OMDOT_to_Mp(ecc, pb, omdot+omdot_e, mc)
    mps2 = p.OMDOT_to_Mp(ecc, pb, omdot-omdot_e, mc)
    pgsls(1);pgsci(1)
    pgline(mps1, mc)
    pgline(mps2, mc)
    pgsls(1);pgsci(1)
except:
    print "No omdot constraint"

# Pbdot
try:
    pbdot = par.PBDOT
    pbdot_e = par.PBDOT_ERR
    mcs1 = np.array([])
    mcs2 = np.array([])
    mp = np.arange(0.5, 2.51, 0.01)
    for m in mp:
        mc1 = p.PBDOT_to_Mc(pbdot-pbdot_e, pb, ecc, m)
        mc2 = p.PBDOT_to_Mc(pbdot+pbdot_e, pb, ecc, m)
        mcs1 = np.append(mcs1, mc1)
        mcs2 = np.append(mcs2, mc2)
    pgsls(5);pgsci(1)
    pgline(mp, mcs1)
    pgline(mp, mcs2)
    pgsls(1);pgsci(1)
except:
    print "No pbdot constraint"

# Gamma
try:
    gamma = par.GAMMA
    gamma_e = par.GAMMA_ERR
    mcs1 = np.array([])
    mcs2 = np.array([])
    mp = np.arange(0.0001, 2.51, 0.01)
    for m in mp:
        mc1 = p.GAMMA_to_Mc(gamma-gamma_e, pb, ecc, m)
        mc2 = p.GAMMA_to_Mc(gamma+gamma_e, pb, ecc, m)
        mcs1 = np.append(mcs1, mc1)
        mcs2 = np.append(mcs2, mc2)
    pgsls(3);pgsci(1)
    pgline(mp, mcs1)
    pgline(mp, mcs2)
    pgsls(1);pgsci(1)
except:
    print "No gamma constraint"

# Precession rate
try:
    omega = par.OMEGA
    omega_e = par.OMEGA_ERR
    mcs1 = np.array([])
    mcs2 = np.array([])
    mp = np.arange(0.0001, 2.51, 0.01)
    for m in mp:
        mc1 = p.OMEGA_to_Mc(omega-omega_e, pb, ecc, m)
        mc2 = p.OMEGA_to_Mc(omega+omega_e, pb, ecc, m)
        mcs1 = np.append(mcs1, mc1)
        mcs2 = np.append(mcs2, mc2)
    pgsls(4);pgsci(2)
    print mp, mcs1
    pgline(mp, mcs1)
    pgline(mp, mcs2)
    pgsls(1);pgsci(1)
except:
    print "No precession rate constraint"


# Mass function                                                                     
mcs1 = np.array([])
for m in mp:
    mc1 = p.companion_mass_limit(pb, x, mpsr=m)
    mcs1 = np.append(mcs1, mc1)
#pgline(mp, mcs1)
mp = np.append([0],mp)
mcs1 = np.append([0],mcs1)
mp = np.append(mp,[2.5])
mcs1 = np.append(mcs1, [0])
#pgsci(8)
pgsci(15)
pgpoly(mp, mcs1)
pgsci(1)


## Now plotting labels

# Pbdot
pgsch(0.8)
pgpt([0.94], [2.42], 17)
pgsch(1.5)
pgptxt(0.94, 2.3, 0., 0.5, "\\fiP\\d\\frb\\u")

# Omega-dot
pgsch(0.8)
pgpt([0.46], [2.35], 17)
pgsch(1.5)
pgptxt(0.46, 2.25, 0., 0.5, "\\gw")

# Gamma
pgptxt(0.19, 0.88, 0., 0.5, "\\gg")

# Omega
pgsci(2)
pgptxt(0.2, 1.3, 0., 0.5, "\\gW\\dp\\u")

# Inclination angle                                                  
pgptxt(1., 0.88, 0., 0.5, "\\fii")
pgsci(1)

pgsch(1)

pgbox("BCNST", 0, 0, "BCNST", 0, 0)

#pgmtxt("B", 2.7, 0.5, 0.5, "M\\dp\\u (M\\d\\(2281)\\u)")
#pgmtxt("L", 2.5, 0.5, 0.5, "M\\dc\\u (M\\d\\(2281)\\u)")
pgmtxt("B", 2.7, 0.5, 0.5, "Pulsar mass, M\\dp\\u (M\\d\\(2281)\\u)")
pgmtxt("L", 2.5, 0.5, 0.5, "Companion mass, M\\dc\\u (M\\d\\(2281)\\u)")
pgend()
