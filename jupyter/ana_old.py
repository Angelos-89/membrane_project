import numpy as np
import matplotlib.pyplot as P
#-----------------------------------------#
def pspec(dir='./',fname='hfield_spec_0.txt',N=80,aa=1):
#dir='/home/angelos-89/runs/pinned/long_runs/pin_0/tau_2.5/'
    spec = np.loadtxt(dir+fname)
    LL = aa*N
    dk = 2*np.pi/LL
    xx = spec[1:,0]/dk
    fig = P.figure()
    ax = fig.add_subplot(111)
    ax.plot( np.log10(xx),np.log10(spec[1:,1]/spec[1,1]), 'o-' )
    fit1 = 1.1e4*xx**(-3)
    fit2 = 50*xx**(-1)
    ax.plot( np.log10(xx),np.log10(fit1), 'k' )
    ax.plot( np.log10(xx),np.log10(fit2), 'k:' )
    ax.grid(True)
    ax.set_xlabel(r'$\log_{10}[\left(\frac{L}{2\pi}\right)k]$')
    ax.set_ylabel(r'$\log_{10}[S(k)]$')
    ax.set_xlim(0,np.log10(N))
    ax.set_ylim(-2,1)
    P.show()
    return(spec)
#-------------------------#
#************************************#
if __name__ == '__main__':
  main();
else:
  print('Importing Analysis module')

