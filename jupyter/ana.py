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
def load_spec( pdir='/home/angelos-89/runs/pinned/pin_0/',fname='hfield_spec_0.txt'):
    ddir = np.array(['tau_0.01/',
        'tau_0.1/',
        'tau_1.0/',
        'tau_10.0/',
        'tau_100.0/',
        'tau_1000.0/'])
    tau = np.array([0.01,
        0.1,
        1.0,
        10,
        100,
        1000])
    ndir = ddir.size
    spec = np.empty(ndir,dtype=object)
    for idir in range(ndir):
        print(pdir+ddir[idir]+'hfield_spec_0.txt')
        spec[idir] = np.loadtxt(pdir+ddir[idir]+fname)
    return(ndir,spec)
#----------------------------------------------#
def get_gamma(spec,ax,NN=80,aa=1,Nmax=32,sigma=1.45,kappa=2.5):
    LL = aa*NN
    dk = 2*np.pi/LL
    qq = spec[1:Nmax,0]
    Sq = spec[1:Nmax,1]
    fitk = 200*(1./(2*np.pi))*qq 
    fitk3 = 300*(kappa/(2*np.pi))*(qq**3) 
    ax.loglog( qq,1./(Sq*sigma), 'o-' )
    ax.plot(qq,fitk,'k')
    ax.plot(qq,fitk3,'k:')
    ax.set_ylim([1,1e4])
    ax.grid()
#-------------------------------
def speck(spec,ax,NN=80,aa=1,Nmax=32,sigma=1.45):
    LL = aa*NN
    dk = 2*np.pi/LL
    qq = spec[1:Nmax,0]
    Sq = spec[1:Nmax,1]
    fit = 30*(sigma/(2*np.pi))*qq 
    ax.loglog( qq,sigma*Sq, 'o-' )
    #ax.plot(qq,fit,'k')
    ax.grid()
# gamma is the fluctuating tension 
#--------------------------------------------------------#
sym = np.array(['o',
              '*',
              'X', 
              'v',
              'd',
              '^'])
sigma = np.array([0.5,
    0.57,
    1.45,
    10.28,
    100,
    1000])
#ndir,spec = load_spec() 
fig = P.figure()
ax = fig.add_subplot(111)
for idir in range(ndir):
#    get_gamma(spec[idir][:,:],ax,Nmax=40,sigma=sigma[idir])
     speck(spec[idir][:,:],ax,Nmax=40,sigma=sigma[idir])
ax.set_ylim([1e-4,1])
ax.grid(True)
#************************************#
