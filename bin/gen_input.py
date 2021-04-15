# Filename : gen_input.py
"""
Generates the input file(s) for the membrane code.
"""
# Default values, change them as necessary.
# This file assumes that we are going to run 8 
# simulations in as many processors with all 
# parameters being the same except frame tension tau.
import numpy as np
iprefix = 'input_'
nproc = 8
# Every other parameters is an array of size nproc
# There are 10 parameters
nparam = 10
param = np.zeros([nparam,nproc])
# You can set some parameters from pre-defined arrays, e.g.
tau = np.zeros(nproc)
for iproc in range(0,nproc):
    tau[iproc] = 1e-1*(10**iproc)
sigma = tau
delta_hmax = np.array([0.3, 0.3, 0.3, 0.1, 0.1, 0.05, 0.05, 0.05])
delta_amin = np.array([0.98, 0.98, 0.98, 0.98, 0.95, 0.95, 0.94, 0.94])
delta_amax = np.array([1.02, 1.02, 1.02, 1.02, 1.12, 1.12, 1.12, 1.12])
# Now loop over the processors and set the parameters
for iproc in range(0,nproc):
    # first parameter nrun = 0 
    param[0,iproc] = 0
    # no. of iterations
    param[1,iproc] = 1e7
    # surface tension (sigma)
    param[2,iproc] = sigma[iproc] 
    # frame tansion (tau)
    param[3,iproc] = tau[iproc] 
    # maximum height perturbation
    param[4,iproc] = delta_hmax[iproc]
    # minimum lattice perturbation
    param[5,iproc] = delta_amin[iproc] 
    # maximum lattice perturbation
    param[6,iproc] = delta_amax[iproc] 
    # fraction of pinning
    param[7,iproc] = 0
    #how often we sample
    param[8,iproc] = 500 
    # Active energy
    param[9,iproc] = 0.
for iproc in range(nproc):
    fname = iprefix+str(iproc)+'.txt'
    ff = open(fname, 'w')
    for iparam in range(nparam):
        ff.write(str(param[iparam,iproc]))
        ff.write("\n")
    ff.close()
#if __name__ == "__main__":
#    import sys
#    fib(int(sys.argv[1]))
