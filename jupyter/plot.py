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
    ax.plot(time,energy)
    return ax
#----------------------------------------#
#iread = 0
ddir=np.array(['active_m2','active_zero','active_2'])
nproc=16
nparam=10
if (iread == 0) :
    all_param,all_ts,all_htext,all_spec=rall(ddir,nproc,nparam)
    iread=1 
    print('read parameters, time series and spectra')
else:
    print('data already read')
fig = P.figure()
ax=fig.add_subplot(111)
#ndir=np.shape(ddir)[0]
#iproc=4
#for idir in range(1,ndir):
#@    ax=pEnergy(all_ts,all_htext,ax,kdir=idir,kproc=iproc)
#    tau = all_param[idir][iproc][4]
#    print('tau=',tau)
#ax = pArea(all_ts,all_htext,ax,kdir=0)
#ax = pArea(all_ts,all_htext,ax,kdir=1)
#ax = pArea(all_ts,all_htext,ax,kdir=2)
idir=2
for iproc in range(0,nproc-1):
    ax,area=pArea(all_ts,all_htext,ax,kdir=idir,kproc=iproc)
    sigma = all_param[idir][iproc][3]
    tau = all_param[idir][iproc][4]
    print('tau,sigma,area=\t',tau,sigma,area)
P.grid()
P.show()
