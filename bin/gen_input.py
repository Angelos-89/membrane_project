##Filename : gen_input.py

"""Generates the input file(s) for the membrane code."""

## Default values, change them as necessary.
## This file assumes that we are going to run
## 16 simulations with all parameters being
## the same, except frame tension tau and
## internal tension sigma.

import numpy as np
iprefix = 'input_'
nparam=11

def gen_param(nproc=1):

    ## Every other parameter is an array of size nproc. There are 10 parameters.
    param = np.zeros([nparam,nproc])

    ## You can set some parameters from pre-defined arrays, e.g.
    tau = np.zeros(nproc)
    tau = [0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.4, 0.8, 1., 2., 4., 8., 10., 20., 40., 80.]
    sigma = [0.5, 0.51, 0.53, 0.55, 0.6, 0.7, 0.9, 1.3, 1.5, 2.4, 4., 9., 10.5, 21.2, 41., 80.]
    #delta_hmax = np.array([0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3])
    #delta_amin = np.array([0.98, 0.98, 0.98, 0.98, 0.95, 0.95, 0.94, 0.94])
    #delta_amax = np.array([1.02, 1.02, 1.02, 1.02, 1.12, 1.12, 1.12, 1.12])

    # Now loop over the processors and set the parameters.
    for iproc in range(nproc):
        # New simulation (0) or start from earlier (1) 
        param[0,iproc] = 0
        # If we write snapshots of the membrane (1)
        param[2,iproc] = 0
        # How often we sample
        param[3,iproc] = 500 
        # No. of iterations
        param[4,iproc] = 1e7
        # Surface tension (sigma)
        param[5,iproc] = sigma[iproc] 
        # Frame tansion (tau)
        param[6,iproc] = tau[iproc] 
        # Maximum height perturbation
        param[7,iproc] = 0.2
        # Minimum lattice perturbation
        param[8,iproc] = 0.98 
        # Maximum lattice perturbation
        param[9,iproc] = 1.02 
        # Fraction of pinning
        param[10,iproc] = 0.
        # Activity energy
        param[11,iproc] = 0.
    return param

def pparam(nproc,param):
    
    for iproc in range(nproc):
        fname = iprefix+str(iproc)+'.txt'
        ff = open(fname, 'w')
        ff.write(str(np.int(param[0,iproc]))) # 1st entry is integer.
        ff.write("\n")
        ff.write(str(np.int(param[1,iproc]))) # This too.
        ff.write("\n")
        ff.write(str(np.int(param[2,iproc]))) # This too.
        ff.write("\n")
        for iparam in range(3,nparam):        # These are doubles.
            ff.write(str(param[iparam,iproc]))
            ff.write("\n")
        ff.close()

if __name__ == "__main__":
    nproc=16
    param=gen_param(nproc=nproc)
    pparam(nproc,param)
