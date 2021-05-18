##Filename : gen_input.py

"""Generates the input file(s) for the membrane code."""

## Default values, change them as necessary.
## This file assumes that we are going to run
## 64 simulations with all parameters being
## the same, except frame tension tau and
## internal tension sigma.

import numpy as np
iprefix = 'input_'
nparam = 11

def gen_param(nproc=1):

    ## Every other parameter is an array of size nproc. There are 10 parameters.
    param = np.zeros([nparam,nproc])

    ## You can set some parameters from pre-defined arrays, e.g.
    tau = np.zeros(nproc)

    #------------------------------
    #0% pinning (unpinned membrane)
    #------------------------------
    '''tau = [0.01, 0.01, 0.01, 0.01,   0.02, 0.02, 0.02, 0.02,   0.04, 0.04, 0.04, 0.04,     0.08, 0.08, 0.08, 0.08,
           0.1,  0.1,  0.1,  0.1,    0.2,  0.2,  0.2,  0.2,    0.4,  0.4,  0.4,  0.4,      0.8,  0.8,  0.8,  0.8,
           1.0,  1.0,  1.0,  1.0,    2.0,  2.0,  2.0,  2.0,    4.0,  4.0,  4.0,  4.0,      8.0,  8.0,  8.0,  8.0,
           10.0, 10.0, 10.0, 10.0,   40.0, 40.0, 40.0, 40.0,   100,  100,  100,  100,      1000, 1000, 1000, 1000]

    sigma = [0.506143, 0.506143, 0.506143, 0.506143,          0.515, 0.515, 0.515, 0.515,         0.5348571,  0.5348571,  0.5348571,  0.5348571,       0.5728571,  0.5728571, 0.5728571, 0.5728571,
             0.5938709, 0.5938709 , 0.5938709,  0.5938709,    0.69,  0.69,  0.69,   0.69,         0.89, 0.89, 0.89,  0.89,                             1.2805,  1.2805, 1.2805,  1.2805,
             1.479, 1.479, 1.479, 1.479,                      2.465,  2.465,  2.465,  2.465,      4.448, 4.448, 4.448,  4.448,                         8.41,  8.41, 8.41,  8.41,
             10.4, 10.4, 10.4,  10.4,                         40.28, 40.280, 40.280, 40.280,      100.18, 100.18, 100.18, 100.18,                      1000.03, 1000.03, 1000.03, 1000.03]
    
    delta_hmax = np.array( [0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,
                            0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,
                            0.16, 0.16, 0.16, 0.16,  0.16, 0.16, 0.16, 0.16,  0.16, 0.16, 0.16, 0.16,  0.15, 0.15, 0.15, 0.15,
                            0.15, 0.15, 0.15, 0.15,  0.13, 0.13, 0.13, 0.13,  0.11, 0.11, 0.11, 0.11,  0.04,0.04,0.04,0.04])
    
    delta_amin = np.array([0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,
                           0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,
                           0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,
                           0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.97, 0.97, 0.97, 0.97,  0.95, 0.95, 0.95, 0.95] )
    
    delta_amax = np.array([1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,
                           1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,
                           1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,
                           1.02, 1.02, 1.02, 1.02,  1.03, 1.03, 1.03, 1.03,  1.03, 1.03, 1.03, 1.03,  1.08, 1.08, 1.08, 1.08])
    pin = 0.0'''

    
    #----------
    #8% pinning
    #----------

    '''tau = [0.01, 0.01, 0.01, 0.01,   0.02, 0.02, 0.02, 0.02,   0.04, 0.04, 0.04, 0.04,     0.08, 0.08, 0.08, 0.08,
           0.1,  0.1,  0.1,  0.1,    0.2,  0.2,  0.2,  0.2,    0.4,  0.4,  0.4,  0.4,      0.8,  0.8,  0.8,  0.8,
           1.0,  1.0,  1.0,  1.0,    2.0,  2.0,  2.0,  2.0,    4.0,  4.0,  4.0,  4.0,      8.0,  8.0,  8.0,  8.0,
           10.0, 10.0, 10.0, 10.0,   40.0, 40.0, 40.0, 40.0,   100,  100,  100,  100,      1000, 1000, 1000, 1000]
    
    sigma = [0.46928,0.46928,0.46928,0.46928,            0.47666, 0.47666, 0.47666, 0.47666,     0.4972,0.4972,0.4972,0.4972,        0.5359677, 0.5359677, 0.5359677, 0.5359677,
             0.5594,0.5594,0.5594,0.5594,        	 0.653, 0.653, 0.653, 0.653,             0.85666,0.85666,0.85666,0.85666,    1.2496825, 1.2496825, 1.2496825, 1.2496825,
             1.451,1.451,1.451,1.451,            	 2.420, 2.420, 2.420, 2.420,             4.42,4.42,4.42,4.42,          	     8.39066, 8.39066, 8.39066, 8.39066,
             10.379677,10.379677,10.379677,10.379677,    40.264, 40.264, 40.264, 40.264,         100.17,100.17,100.17,100.17,  	     1000.03, 1000.03, 1000.03, 1000.03]
    
    delta_hmax = np.array( [0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,
                            0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,
                            0.17, 0.17, 0.17, 0.17,  0.17, 0.17, 0.17, 0.17,  0.16, 0.16, 0.16, 0.16,  0.16, 0.16, 0.16, 0.16,
                            0.15,0.15,0.15,0.15,     0.13,0.13,0.13,0.13,     0.10,0.10,0.10,0.10,     0.045, 0.045, 0.045, 0.045])
    
    delta_amin = np.array([0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,
                           0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,
                           0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,
                           0.98, 0.98, 0.98, 0.98,  0.97, 0.97, 0.97, 0.97,  0.97, 0.97, 0.97, 0.97,  0.96, 0.96, 0.96, 0.96])
    
    delta_amax = np.array([1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,
                           1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,
                           1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,
                           1.02, 1.02, 1.03, 1.03,  1.03, 1.03, 1.03, 1.03,  1.03, 1.03, 1.04, 1.04,  1.09, 1.09, 1.09, 1.09])
    pin = 0.08'''

    
    #-----------
    #16% pinning
    #-----------

    '''tau = [0.01, 0.01, 0.01, 0.01,   0.02, 0.02, 0.02, 0.02,   0.04, 0.04, 0.04, 0.04,     0.08, 0.08, 0.08, 0.08,
           0.1,  0.1,  0.1,  0.1,    0.2,  0.2,  0.2,  0.2,    0.4,  0.4,  0.4,  0.4,      0.8,  0.8,  0.8,  0.8,
           1.0,  1.0,  1.0,  1.0,    2.0,  2.0,  2.0,  2.0,    4.0,  4.0,  4.0,  4.0,      8.0,  8.0,  8.0,  8.0,
           10.0, 10.0, 10.0, 10.0,   40.0, 40.0, 40.0, 40.0,   100,  100,  100,  100,      1000, 1000, 1000, 1000]
    
    sigma = [0.427,0.427,0.427,0.427,		           0.438,0.438,0.438,0.438,        0.4602,0.4602,0.4602,0.4602,      0.496,0.496,0.496,0.496,
             0.519333,0.519333,0.519333,0.519333,          0.615,0.615,0.615,0.615,        0.816,0.816,0.816,0.816,          1.21333,1.21333, 1.21333,1.21333,
             1.41,1.41,1.41,1.41,         		   2.400,2.40,2.40,2.40,           4.388,4.3888,4.3888,4.3888,       8.365,8.365,8.365,8.365,
             10.3522580,10.3522580,10.3522580,10.3522580,  40.245,40.245,40.245,40.245,    100.16,100.16,100.16,100.16,      1000.0276984,1000.0276984,1000.0276984,1000.0276984]
    
    delta_hmax = np.array([0.15, 0.15, 0.14, 0.14,  0.15, 0.15, 0.14, 0.14,  0.15, 0.15, 0.14, 0.14,  0.15, 0.15, 0.14, 0.14,
                           0.15, 0.15, 0.14, 0.14,  0.15, 0.15, 0.14, 0.14,  0.15, 0.15, 0.14, 0.14,  0.15, 0.15, 0.14, 0.14,
                           0.15, 0.15, 0.14, 0.14,  0.15, 0.15, 0.14, 0.14,  0.15, 0.15, 0.15, 0.15,  0.15, 0.15, 0.15, 0.15,
                           0.13,0.13,0.12,0.12,     0.11, 0.11, 0.11, 0.11,  0.09, 0.09, 0.09, 0.09,  0.04, 0.04, 0.04, 0.04] )
    
    delta_amin = np.array([0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,
                           0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,
                           0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,
                           0.97, 0.97, 0.97, 0.97,  0.97, 0.97, 0.97, 0.97,  0.97, 0.97, 0.97, 0.97,  0.95, 0.95, 0.95, 0.95] )
    
    delta_amax = np.array([1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,
                           1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,
                           1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,  1.02, 1.02, 1.02, 1.02,
                           1.02, 1.02, 1.02, 1.03,  1.02, 1.02, 1.02, 1.02,  1.03, 1.03, 1.03, 1.03,  1.06, 1.06, 1.06, 1.06 ])

    pin = 0.16'''
    
    #-----------
    #32% pinning
    #-----------

    '''tau = [0.01, 0.01, 0.01, 0.01,   0.02, 0.02, 0.02, 0.02,   0.04, 0.04, 0.04, 0.04,     0.08, 0.08, 0.08, 0.08,
           0.1,  0.1,  0.1,  0.1,    0.2,  0.2,  0.2,  0.2,    0.4,  0.4,  0.4,  0.4,      0.8,  0.8,  0.8,  0.8,
           1.0,  1.0,  1.0,  1.0,    2.0,  2.0,  2.0,  2.0,    4.0,  4.0,  4.0,  4.0,      8.0,  8.0,  8.0,  8.0,
           10.0, 10.0, 10.0, 10.0,   40.0, 40.0, 40.0, 40.0,   100,  100,  100,  100,      1000, 1000, 1000, 1000]
    
    sigma = [0.350,0.350,0.350,0.350,                             0.36,0.360,0.36,0.36,                               0.38,0.38,0.38,0.38,                                  0.42,0.42,0.42,0.42,
             0.440,0.440,0.440,0.440,                             0.540,0.540,0.540,0.540,                            0.74,0.74,0.74,0.74,                                  1.139,1.139,1.139,1.139,
             1.335,1.335,1.335,1.335,                             2.329677, 2.329677, 2.329677, 2.329677,             4.3198412, 4.3198412, 4.3198412, 4.3198412,           8.3 ,8.3 , 8.30, 8.3,
             10.2913492, 10.2913492, 10.2913492, 10.2913492,      40.21079365,40.21079365,40.21079365,40.21079365,    100.13,100.13,100.13,100.13,                          1000.02419047,1000.02419047,1000.02419047, 1000.02419047]	
    
    delta_hmax = np.array( [0.10, 0.10, 0.10, 0.10,  0.10, 0.10, 0.10, 0.10,  0.10, 0.10, 0.10, 0.10,  0.10, 0.10, 0.10, 0.10,
                            0.10, 0.10, 0.10, 0.10,  0.10, 0.10, 0.10, 0.10,  0.10, 0.10, 0.10, 0.10,  0.10, 0.10, 0.10, 0.10,
                            0.10, 0.10, 0.10, 0.10,  0.10, 0.10, 0.10, 0.10,  0.10, 0.10, 0.10, 0.10,  0.09, 0.09, 0.09, 0.09,
                            0.09,0.09,0.09,0.09,     0.08,0.08,0.08,0.08,     0.08,0.08,0.08,0.08,     0.035, 0.035, 0.035, 0.035] )
    
    delta_amin = np.array([0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,
                           0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,
                           0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.98, 0.98, 0.98, 0.98,  0.97, 0.97, 0.97, 0.97,
                           0.97, 0.97, 0.97, 0.97,  0.97, 0.97, 0.97, 0.97,  0.97, 0.97, 0.97, 0.97,  0.94, 0.94, 0.94, 0.94] )
    
    delta_amax = np.array([1.02, 1.02, 1.02, 1.02,      1.02, 1.02, 1.02, 1.02,      1.02, 1.02, 1.02, 1.02,     1.02, 1.02, 1.02, 1.02,
                           1.02, 1.02, 1.02, 1.02,      1.02, 1.02, 1.02, 1.02,      1.02, 1.02, 1.02, 1.02,     1.02, 1.02, 1.02, 1.02,
                           1.02, 1.02, 1.02, 1.02,      1.02, 1.02, 1.02, 1.02,      1.02, 1.02, 1.02, 1.02,     1.03, 1.03, 1.03, 1.03,
                           1.035, 1.035, 1.035, 1.035,  1.045, 1.045, 1.045, 1.045,  1.05, 1.05, 1.05, 1.05,     1.12, 1.12, 1.12, 1.12])

    pin = 0.32'''
    
    # Now loop over the processors and set the parameters.
    for iproc in range(nproc):
        # New simulation (0) or start from earlier (1) 
        param[0,iproc] = 0
        # If we write snapshots of the membrane (1)
        param[1,iproc] = 0
        # How often we sample
        param[2,iproc] = 10000 
        # No. of iterations
        param[3,iproc] = 1e8
        # Surface tension (sigma)
        param[4,iproc] = sigma[iproc] 
        # Frame tansion (tau)
        param[5,iproc] = tau[iproc] 
        # Maximum height perturbation
        param[6,iproc] = delta_hmax[iproc]
        # Minimum lattice perturbation
        param[7,iproc] = delta_amin[iproc] 
        # Maximum lattice perturbation
        param[8,iproc] = delta_amax[iproc] 
        # Fraction of pinning
        param[9,iproc] = pin
        # Activity energy
        param[10,iproc] = 0.
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
    nproc=64
    param=gen_param(nproc=nproc)
    pparam(nproc,param)
