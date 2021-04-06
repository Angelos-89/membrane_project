import numpy as np
import matplotlib.pyplot as P
spec = np.loadtxt('hfield_spec_0.txt')
N=80
aa = 1.
LL = aa*N
dk = 2*np.pi/LL
xx = spec[1:,0]/dk
P.plot( np.log10(xx),np.log10(spec[1:,1]/spec[1,1]), 'o-' )
fit1 = 1.1e4*xx**(-3)
fit2 = 50*xx**(-1)
P.plot( np.log10(xx),np.log10(fit1), 'k' )
P.plot( np.log10(xx),np.log10(fit2), 'k:' )
P.grid(True)
P.xlabel(r'$\log_{10}[\left(\frac{L}{2\pi}\right)k]$')
P.ylabel(r'$\log_{10}[S(k)]$')
P.xlim(0,np.log10(N))
P.ylim(-2,1)
P.show()
