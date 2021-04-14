import numpy as np
import matplotlib.pyplot as P
import matplotlib.ticker as mticker
f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
fmt = mticker.FuncFormatter(g)
#dir='./'
p0dir='/home/angelos-89/runs/pinned/long_runs/pin_0/'
ddir = np.array(['tau_0.01/',
              'tau_0.1/',
              'tau_2.5/',
              'tau_10.0/',
              'tau_100.0/',
              'tau_1000.0/'])
tau = np.array([0.01,
              0.1,
              2.5,
              10,
              100,
              1000])
sym = np.array(['o',
              '*',
              'X', 
              'v',
              'd',
              '^'])
ndir = ddir.size
spec = np.empty(ndir,dtype=object)
for idir in range(ndir):
    print(p0dir+ddir[idir]+'hfield_spec_0.txt')
    spec[idir] = np.loadtxt(p0dir+ddir[idir]+'hfield_spec_0.txt')
dactive = '/home/dhruba/research/MC_membrane_runs/active_Apr20/r1/'
specA = np.loadtxt(dactive+'hfield_spec_0.txt')
N=80
aa = 1.
LL = aa*N
dk = 2*np.pi/LL
for idir in range(ndir):
#idir=0
    xx = spec[idir][1:,0]/dk
#P.plot( np.log10(xx),np.log10(spec[idir][1:,1]/spec[idir][1,1]), sym[idir] )
    P.loglog(xx,spec[idir][1:,1]/spec[idir][1,1], sym[idir], 
        label=r"$\tau$ = {}".format(fmt(tau[idir])) )
fit_tension = 1.*xx**(-1)
P.loglog(xx,fit_tension,'k')
fit_curv = 50*xx**(-3)
P.loglog(xx,fit_curv,'k:')
xx = specA[1:,0]/dk
P.loglog(xx,specA[1:,1]/specA[1,1], 's-',color='C10') 


#P.plot( np.log10(xx),np.log10(fit1), 'k' )
#P.plot( np.log10(xx),np.log10(fit2), 'k:' )
P.grid(True)
P.legend()
P.xlabel(r'$\log_{10}[\left(\frac{L}{2\pi}\right)k]$')
P.ylabel(r'$\log_{10}[S(k)]$')
P.xlim(1,80)
P.ylim(1e-4,10)
P.show()
