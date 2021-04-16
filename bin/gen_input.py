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
nproc = 16 
# Every other parameters is an array of size nproc
# There are 10 parameters
nparam = 10
param = np.zeros([nparam,nproc])
# You can set some parameters from pre-defined arrays, e.g.
tau = np.zeros(nproc)
tau = [0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.4, 0.8, 1., 2., 4., 8., 10., 20., 40., 80.]
sigma = [0.5, 0.51, 0.52, 0.54, 0.57, 0.75, 0.8, 1.24, 1.44, 2.5, 4.8, 8.3, 10.3, 20.2, 40.2, 80.] 
#delta_hmax = np.array([0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3])
#delta_amin = np.array([0.98, 0.98, 0.98, 0.98, 0.95, 0.95, 0.94, 0.94])
#delta_amax = np.array([1.02, 1.02, 1.02, 1.02, 1.12, 1.12, 1.12, 1.12])
# Now loop over the processors and set the parameters
for iproc in range(0,nproc):
    # first parameter nrun = 0 
    param[0,iproc] = 0
    #how often we sample
    param[1,iproc] = 500 
    # no. of iterations
    param[2,iproc] = 1e7
    # surface tension (sigma)
    param[3,iproc] = sigma[iproc] 
    # frame tansion (tau)
    param[4,iproc] = tau[iproc] 
    # maximum height perturbation
    param[5,iproc] = 0.3 
    # minimum lattice perturbation
    param[6,iproc] = 0.98 
    # maximum lattice perturbation
    param[7,iproc] = 1.02 
    # fraction of pinning
    param[8,iproc] = 0
    # Active energy
    param[9,iproc] = 0.
for iproc in range(nproc):
    fname = iprefix+str(iproc)+'.txt'
    ff = open(fname, 'w')
    ff.write(str(np.int(param[0,iproc])))
    ff.write("\n")
    ff.write(str(np.int(param[1,iproc])))
    ff.write("\n")
    for iparam in range(2,nparam):
        ff.write(str(param[iparam,iproc]))
        ff.write("\n")
    ff.close()
#if __name__ == "__main__":
#    import sys
#    fib(int(sys.argv[1]))
