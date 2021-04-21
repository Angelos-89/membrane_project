import numpy as np
import matplotlib.pyplot as P
import re as re
header = np.array([ '%iter','total_moves','total_area',
                    'prj_area','alpha','curv_energy',
                    'entropic_corr', 'pinning_energy', 'tot_energy'])
#--------------------------------------------------------
def which_col(htext,string='tot_energy'):
    ncol = np.shape(htext)[0]
    for icol in range(ncol):
        match = re.match(htext[icol],string)
        if match:
            col = icol
            break
    return match,col
#------------------------------
def read_ts(ddir='./',fname='timeseries_0.txt',comment='%'):
    f = open(ddir+'/'+fname,'r')
    txt = f.readline();
    htxt = txt.split();
    ts = np.loadtxt(ddir+'/'+fname,comments=comment)
    return ts,htxt
#-----------------------------
def parea(ddir='./',fname='timeseries_0.txt',comment='%'):
    ts,htext = read_ts(ddir=ddir,fname=fname,comment=comment)
    match,col_area=which_col(htext,string='total_area')
    match,col_time=which_col(htext,string='total_moves')
    fig = P.figure()
    ax = fig.add_subplot(111)
    ax.plot(ts[:,col_time],ts[:,col_area])
    return(ax)
#-----------------------------
def rparam(ddir='./',nproc=16,nparam=10):
    param = np.zeros([nproc,nparam])
    for iproc in range(nproc):
        fname=ddir+'/'+'input_'+str(iproc)+'.txt'
        f = open(fname,'r')
        for iparam in range(nparam):
           param[iproc,iparam]=f.readline()
    return param
#-----------------------------
def rts_dir(ddir='./',nproc=16):
    rts = np.empty(nproc,dtype=object)
    rhtext = np.empty(nproc,dtype=object)
    for iproc in range(nproc):
        fname='timeseries_'+str(iproc)+'.txt'
        ts,htxt = read_ts(ddir=ddir,fname=fname)
        rts[iproc] = ts
        rhtext[iproc]=htxt
    return rts,rhtext
#-----------------------------
def rspec_dir(ddir='./',nproc=16,NN=80,aa=1):
    spec_dir = np.empty(nproc,dtype=object)
    LL = aa*NN
    dk = 2*np.pi/LL
    for iproc in range(nproc):
        spec = np.loadtxt(ddir+'/'+'hfield_spec_'+str(iproc)+'.txt')
        spec[:,0] = spec[:,0]/dk
        spec_dir[iproc] = spec
    return spec_dir
#-----------------------------
def rall(ddir,nproc,nparam):
    ndir = np.shape(ddir)[0]
    all_param = np.empty(ndir,dtype=object)
    all_ts = np.empty(ndir,dtype=object)
    all_htext = np.empty(ndir,dtype=object)
    all_spec=np.empty(ndir,dtype=object)
    for idir in range(ndir):
        param = rparam(ddir=ddir[0])
        all_param[idir] = param
        rts,rhtext = rts_dir(ddir=ddir[idir],nproc=nproc)
        all_ts[idir] = rts
        all_htext[idir] = rhtext
        spec_dir = rspec_dir(ddir=ddir[idir],nproc=nproc)
        all_spec[idir] = spec_dir 
    return all_param,all_ts,all_htext,all_spec
#----------------------------------------#
def get_ts_data(all_ts,all_htext,kdir=0,kproc=0,string='total_area',tstring='total_moves'):
    htext = all_htext[kdir][kproc]
    match,column=which_col(htext,string=string)
    match,col_time=which_col(htext,string=tstring)
    ts = all_ts[kdir][kproc]
    time=ts[:,col_time]
    data=ts[:,column]
    return time,data
#----------------------------------------#
def pArea(all_ts,all_htext,ax,kdir=0,kproc=0,window=0.8):
    time,area=get_ts_data(all_ts,all_htext,kdir=kdir,kproc=kproc,string='total_area')
    ntmax=np.shape(area)[0]
    ntmin=np.int(ntmax*0.8)
    print(ntmin,ntmax)
    mean_area=np.mean(area[ntmin:ntmax-1])
    ax.plot(time,area)
    return ax,mean_area
#----------------------------------------#
def pEnergy(all_ts,all_htext,ax,kdir=0,kproc=0):
    time,energy=get_ts_data(all_ts,all_htext,kdir=kdir,kproc=kproc,string='tot_energy')
    ntmax=np.shape(energy)[0]
    ntmin=np.int(ntmax*0.8)
    print(ntmin,ntmax)
    meanE=np.mean(energy[ntmin:ntmax-1])
    ax.plot(time,energy)
    return ax,meanE
#----------------------------------------#
def ptension(all_param,ax,nproc=16,kdir=0):
    sigma = np.empty(nproc)
    tau = np.empty(nproc)
    for iproc in range(0,nproc):
#    ax,area=pArea(all_ts,all_htext,ax,kdir=idir,kproc=iproc)
        sigma[iproc] = all_param[kdir][iproc][3]
        tau[iproc] = all_param[kdir][iproc][4]
    ax.loglog(tau,sigma,'o-')
    return ax,sigma,tau
#----------------------------------------#
def cArea(all_ts,all_htext,ax,nproc=16,kdir=0):
    for iproc in range(0,nproc):
        ax,meanA=pArea(all_ts,all_htext,ax,kdir=kdir,kproc=iproc)
        sigma = all_param[kdir][iproc][3]
        tau = all_param[kdir][iproc][4]
        print('tau,sigma,area=\t',tau,sigma,meanA)
    ax.set_xlabel('time')
    ax.set_ylabel('Area')
    return ax
#----------------------------------------#
def cEn(all_ts,all_htext,ax,nproc=16,kdir=0):
    for iproc in range(0,nproc):
        ax,meanE=pEnergy(all_ts,all_htext,ax,kdir=idir,kproc=iproc)
        sigma = all_param[idir][iproc][3]
        tau = all_param[idir][iproc][4]
        print('tau,sigma,Energy=\t',tau,sigma,meanE)
    ax.set_xlabel('time')
    ax.set_ylabel('Energy')
    return ax
#-------------------------------------------
def ctension(all_param,ndir,ax,nproc=16):
    sigma = np.zeros([ndir,nproc])
    tau = np.zeros([ndir,nproc])
    for idir in range(ndir):
        for iproc in range(0,nproc):
#    ax,area=pArea(all_ts,all_htext,ax,kdir=idir,kproc=iproc)
            sigma[idir][iproc] = all_param[idir][iproc][3]
            tau[idir][iproc] = all_param[idir][iproc][4]
    for idir in range(ndir):
        ax.loglog(tau[idir][:],sigma[idir][:],'o-')
    return ax
#-------------------------------------------
def pten():
    fig = P.figure()
    ax=fig.add_subplot(111)
    sigma0 = [0.5,  0.515, 0.525, 0.545, 0.58, 0.69, 0.85, 1.28, 1.48, 2.45, 4.1, 8.5, 10.4, 20.3, 40.05, 80.1]
    sigma2 = [0.28, 0.28, 0.28, 0.32, 0.34, 0.43, 0.6, 1.10, 1.1, 2.1, 4.1, 8.01, 10.01, 20.01, 40.01, 80.1]
    sigma4 = [0.1, 0.105, 0.15, 0.155, 0.18, 0.22, 0.25, 0.81, 1.1, 2.1, 4.1, 8.1, 10.1, 20.1, 40.1, 80.1]
    tau = [0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.4, 0.8, 1., 2., 4., 8., 10., 20., 40., 80.] 
    ax.loglog(tau,sigma0,'s-',label='equi')
    ax.loglog(tau,sigma2,'s-',label=r'$E_{\rm A} = 2$')
    ax.loglog(tau,sigma4,'s-',label=r'$E_{\rm A} = 4$')
    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel(r'$\sigma$')
    ax.legend()
    ax.grid(True)
    return ax
#-------------------------------------------
iread = 0
ddir=np.array(['active_2'])
#ddir=np.array(['active_4'])
nproc=16
nparam=10
if (iread == 0) :
    all_param,all_ts,all_htext,all_spec=rall(ddir,nproc,nparam)
    iread=1 
    print('read parameters, time series and spectra')
else:
    print('data already read')
ndir = np.shape(ddir)[0]
#--------------------------------------
fig = P.figure()
ax=fig.add_subplot(111)
ax = cArea(all_ts,all_htext,ax,nproc=nproc,kdir=0)
#ax.legend()
#ax=pten()
P.show()
#P.savefig('sigma_tau.pdf')
