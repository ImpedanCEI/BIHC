'''
User defined filling schemes defined as
a sequence of true/falses. A True indicates
de presence of a bunch and False indicates
an empty slot
'''

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