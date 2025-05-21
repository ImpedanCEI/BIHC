# copyright ################################# #
# This file is part of the BIHC Package.      #
# Copyright (c) CERN, 2024.                   #
# ########################################### #

'''
User defined filling schemes defined as
a sequence of true/falses. A True indicates
de presence of a bunch and False indicates
an empty slot
'''

# SPS user defined filling general
def fillingSchemeSPS(nbunches, ntrains, batchspacing=9):
    '''
    Returns the filling scheme for the SPS

    Parameters
    ----------
    ntrains: number of injections (batches)
    '''
    # Define filling scheme: parameters
    # ntrains = 4 # number of trains/batches
    nslots = 924 # Defining total number of slots for SPS
    # nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72.. 
    # batchspacing = 9 # Batch spacing in 25 ns slots 45/5

    # Defining the trains as lists of True/Falses
    bt = [True]*nbunches
    st = [False]*batchspacing
    sc = [False]*(nslots - (nbunches+batchspacing)*ntrains)
    an = (bt + st)*ntrains + sc

    return an

# LHC user defined filling general
def fillingSchemeLHC(ninj, nbunches, ntrains, batchspacing=7, injspacing=37):
    '''
    Returns the filling scheme for the LHC 
    using the standard pattern: 

    "nbunches x ntrains for ninj"

    Parameters
    ----------
    ninj: int
        number of injections of the "ntrains x nbunches" scheme
    nbunches: optional, default 72
        number of bunches per train
    ntrains: optional, default 4
        number of trains of bunches per single injection
    '''

    # Define filling scheme: parameters
    # ninj = 14 # Defining number of injections
    nslots = 3564 # Defining total number of slots for LHC
    # ntrain = 4 # Defining the number of trains
    # nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72.. 
    batchS = batchspacing # Batch spacing in 25 ns slots
    injspacing = injspacing # Injection spacing in 25 ns slots

    # Defining the trains as lists of True/Falses
    bt = [True]*nbunches
    st = [False]*batchS
    stt = [False]*injspacing
    sc = [False]*(nslots-(ntrains*nbunches*ninj+((ntrains-1)*(batchS)*ninj)+((1)*injspacing*(ninj))))
    an1 = (bt+st)*ntrains + stt
    an = an1 * ninj + sc # This is the final true false sequence that is the beam distribution

    return an

# -------------------------------------
# Case-specific filling schemes:

def fillingSchemeSPS_standard(ntrains):
    '''
    Returns standard the filling scheme for the SPS

    Parameters
    ----------
    ntrains: number of trains (batches)
    '''
    # Define filling scheme: parameters
    nslots = 920 # Defining total number of slots for SPS
    nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72.. 
    batchspacing = 9 # Batch spacing in 25 ns slots 45/5

    # Defining the trains as lists of True/Falses
    bt = [True]*nbunches
    st = [False]*batchspacing
    sc = [False]*(nslots - (nbunches+batchspacing)*ntrains)
    an = (bt + st)*ntrains + sc

    return an

def fillingSchemeSPS_BCMS(ninj):
    '''
    Returns the filling scheme for the SPS

    Parameters
    ----------
    ninj: number of injections (batches)
    nbunches: default, 48. number of bunches per train
    '''
    # Define filling scheme: parameters
    ntrain = 1 # SPS has 1 train per cycle
    nslots = 920 # Defining total number of slots for SPS
    nbunches = 48 # Defining a number of bunchs e.g. 18, 36, 72.. 
    batchspacing = 8 # Batch spacing in 25 ns slots (200 ns)

    # Defining the trains as lists of True/Falses
    bt = [True]*nbunches
    st = [False]*batchspacing
    sc = [False]*(nslots - (nbunches+batchspacing)*ninj)
    an = (bt + st)*ninj + sc

    return an

def fillingSchemeSPS_8b4e(ninj):
    '''
    Returns the filling scheme for the SPS 
    using the 8b4e pattern

    Parameters
    ----------
    ninj: number of injections (batches)
    '''
    # Define filling scheme: parameters
    nslots = 920 # Defining total number of slots for SPS
    nbunches = 8*7 # Defining number of bunches e.g. 18, 36, 72.. 
    nempty = 4*6  #Defining number of empty slots between bunches
    batchspacing = 8 # Batch spacing in 25 ns slots (200 ns)

    # Defining the trains as lists of True/Falses
    bt = ([True]*8+[False]*4)*6+[True]*8
    st = [False]*batchspacing
    sc = [False]*(nslots - (nbunches+nempty+batchspacing)*ninj)
    an = (bt + st)*ninj + sc

    return an

def fillingSchemeLHC_standard(ninj, nbunches=72, ntrains=4):
    '''
    Returns the filling scheme for the LHC 
    using the standard pattern: 

    "nbunches x ntrains for ninj"

    Parameters
    ----------
    ninj: int
        number of injections of the "ntrains x nbunches" scheme
    nbunches: optional, default 72
        number of bunches per train
    ntrains: optional, default 4
        number of trains of bunches per single injection
    '''

    # Define filling scheme: parameters
    ninj = ninj # Defining number of injections
    nslots = 3564 # Defining total number of slots for LHC
    ntrain = ntrains # Defining the number of trains
    nbunches = nbunches # Defining a number of bunchs e.g. 18, 36, 72.. 
    batchS = 7 # Batch spacing in 25 ns slots
    injspacing = 37 # Injection spacing in 25 ns slots

    # Defining the trains as lists of True/Falses
    bt = [True]*nbunches
    st = [False]*batchS
    stt = [False]*injspacing
    sc = [False]*(nslots-(ntrain*nbunches*ninj+((ntrain-1)*(batchS)*ninj)+((1)*injspacing*(ninj))))
    an1 = (bt+st)*ntrains + stt
    an = an1 * ninj + sc # This is the final true false sequence that is the beam distribution

    return an

# SPS user defined filling scheme: AWAKE single bunch
def fillingSchemeSPS_AWAKE(nbunches=1, ntrains=1):
    '''
    Returns the filling scheme for the SPS

    Parameters
    ----------
    ntrains: number of injections (batches)
    '''
    # Define filling scheme: parameters
    #ntrains = 4 # number of trains/batches
    nslots = 924 # Defining total number of slots for SPS
    #nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72.. 
    batchspacing = 9 # Batch spacing in 25 ns slots 45/5

    # Defining the trains as lists of True/Falses
    bt = [True]*nbunches
    st = [False]*batchspacing
    sc = [False]*(nslots - (nbunches+batchspacing)*ntrains)
    an = (bt + st)*ntrains + sc

    return an