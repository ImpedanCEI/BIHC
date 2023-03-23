'''
Test for numeric spectrum with 
different profile shape

* date: 23/03/2023
* author: Elena de la Fuente
'''

import bihc
import os

def test_numeric_gaussian():

    prof = 'GAUSSIAN'
    file = '25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
    impedance_file = 'PillboxImpedance.txt'
    path = 'tests/'
    beam = bihc.Beam(LPCfile=path+file, bunchShape=prof, ppbk=1000, verbose=False)

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f) 
    Z.getImpedanceFromCST(path+impedance_file)

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == 2183.076195744422

def test_numeric_qgaussian():
    prof = 'q-GAUSSIAN'
    file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
    impedance_file = 'PillboxImpedance.txt'
    path = 'tests/'
    beam = bihc.Beam(LPCfile=path+file, bunchShape=prof, ppbk=1000, verbose=False)

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f) 
    Z.getImpedanceFromCST(path+impedance_file)

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == 2839.640870711469

def test_numeric_parabolic():
    prof = 'PARABOLIC'
    file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
    impedance_file = 'PillboxImpedance.txt'
    path = 'tests/'
    beam = bihc.Beam(LPCfile=path+file, bunchShape=prof, ppbk=1000, verbose=False)

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f) 
    Z.getImpedanceFromCST(path+impedance_file)

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == 4408.417172705872

def test_numeric_cos2():
    prof = 'COS2'
    file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
    impedance_file = 'PillboxImpedance.txt'
    path = 'tests/'
    beam = bihc.Beam(LPCfile=path+file, bunchShape=prof, ppbk=1000, verbose=False)

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f) 
    Z.getImpedanceFromCST(path+impedance_file)

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == 1652.8965230665121 #This value is probably off (Bunch length definition)

def test_numeric_binomial():
    prof = 'BINOMIAL'
    file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
    impedance_file = 'PillboxImpedance.txt'
    path = 'tests/'
    beam = bihc.Beam(LPCfile=path+file, bunchShape=prof, ppbk=1000, verbose=False)

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f) 
    Z.getImpedanceFromCST(path+impedance_file)

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == 445.484665898145 #This value is probably off (Bunch length definition)