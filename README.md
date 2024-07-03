# BIHC

[![Documentation Status](https://readthedocs.org/projects/bihc/badge/?version=latest)](https://bihc.readthedocs.io/en/latest/?badge=latest) [![PyPI version](https://badge.fury.io/py/bihc.svg)](https://badge.fury.io/py/bihc)

Beam Induced Heating Computation (BIHC) tool is a package that allows the estimation of the dissipated power due to the passage of a particle beam inside an accelerator component.

The dissipated power value depends on the characteristics of the particle beam (beam spectrum and intensity) and on the characteristics of the consdiered accelerator component (beam-coupling impedance).

Check :file_folder: `examples/` on how to use it!

Documentation is avaiable in [bihc.readthedocs.io](https://bihc.readthedocs.io/en/latest/). More practical information and code snippets in the `Users guide` section.

For specific needs, please contact the maintainers :woman_technologist: :man_technologist: :wave:
* elena.de.la.fuente.garcia@cern.ch
* leonardo.sito@cern.ch

:bookmark: Citing `bihc`
---
There is a [paper about `bihc`](10.18429/JACoW-HB2023-THBP52), presented at 8th ICFA Advanced Beam Dynamics Workshop on High-Intensity and High-Brightness Hadron Beam (0ct. 2023).
If you are using `bihc` in your scientific research, please help our scientific
visibility by citing our work:
> [1] E. de la Fuente, L. Sito, F. Giordano, G. Rumolo, B. Salvant, and C. Zannini, “A Python Package to Compute Beam-Induced Heating in Particle Accelerators and Applications,” JACoW, vol. HB2023, pp. 611–614, 2024, doi: [https://doi.10.18429/JACoW-HB2023-THBP52](https://jacow.org/hb2023/papers/thbp52.pdf). 

Bibtex:
```
@article{Sito:2024ywv,
    author = "Sito, Leonardo and de la Fuente, Elena and Giordano, Francesco and Rumolo, Giovanni and Salvant, Benoit and Zannini, Carlo ",
    title = "{A Python Package to Compute Beam-Induced Heating in Particle Accelerators and Applications}",
    doi = "10.18429/JACoW-HB2023-THBP52",
    journal = "JACoW",
    volume = "HB2023",
    pages = "611--614",
    year = "2024"
}
```
:mag_right: About `bihc` python package 
---

`bihc` is a computational package that integrates over a decade of experience in beam-induced heating calculations from the Impedance and Coherent Effects Section (see [^1], [^3], [^4], [^5], [^6], [^7]) into a comprehensive and flexible Python-based tool. 

The package has been presented at the 68th ICFA Advanced Beam Dynamics Workshop on High-Intensity and High-Brightness Hadron Beam (0ct. 2023)[^8], and is under continuous development to face the beam-induce heating challenges that become more relevant as the beam total intensity and bunch length is pushed.

`bihc` has been succesfully employed to assess the mitigation strategy for the CERN-SPS Beam Wire Scanners after the wire failure in 2023, and was extensively used to study the CERN-LHC Warm Vacuum modules limitations in intensity and bunch length for the 2024 run.

[^1]: B. Salvant et al., “Beam induced heating”, 2012, [Online]. Available: https://cds.cern.ch/record/1975499 

[^3]: C. Zannini, et al. "Power Loss Calculation in Separated and Common Beam Chambers of the LHC". Proceedings of the 5th Int. Particle Accelerator Conf., vol. IPAC2014, 2014, p. 3 pages, 1.928 MB. DOI.org (Datacite), https://doi.org/10.18429/JACOW-IPAC2014-TUPRI061. 

[^4]: C. Zannini, "Electromagnetic Simulation of CERN accelerator Components and Experimental Applications", 2013. [Online]. Available: https://cds.cern.ch/record/1561199 

[^5]: C. Zannini, “Multiphysics Simulations of Impedance Effects in Accelerators,” CERN Yellow Rep. Conf. Proc., vol. 1, pp. 141–144, 2018, doi: 10.23732/CYRCP-2018-001.141. 

[^6]: G. Rumolo, “Beam Instabilities”, 21 pages contribution to the CAS - CERN Accelerator School: Advanced Accelerator Physics Course, Trondheim, Norway, 2014, doi: 10.5170/CERN-2014-009.199. Available; https://cds.cern.ch/record/1982422

[^7]: F. Giordano, ‘Simulation Analysis and Machine Learning Based Detection of Beam-Induced Heating in Particle Accelerator at CERN’, University of Naples Federico II, 2020.

[^8]: L. Sito, E. de la Fuente, F. Giordano, G. Rumolo, B. Salvant, and C. Zannini, “A Python Package to Compute Beam-Induced Heating in Particle Accelerators and Applications,” in Proc. 68th Adv. Beam Dyn. Workshop High-Intensity High-Brightness Hadron Beams (HB’23), Geneva, Switzerland, Apr. 2024, no. 68, pp. 611–614. doi: 10.18429/JACoW-HB2023-THBP52. 

## :zap: Installation
This section explains how to set up the environment to start using BIHC package for power loss computations.
If you are using PyVista in your scientific research, please help our scientific
visibility by citing our work:


#### Developers: Download BIHC repository from Github
```
# SSH:
git clone git@github.com:ImpedanCEI/BIHC.git

# or HTTPS:
git clone https://github.com/ImpedanCEI/BIHC.git
```

#### Users: pip install 
```
pip install bihc
```
If already installed but want to have the newest version: `pip install bihc --upgrade`

### Installation with pytimber in CERN lxplus

Connect to CERN lxplus via ssh. Avoid connecting to lxplus8, the code will induce in Kerberos issues. Kerberos logging will expire 4h after each connection and needs to be renewed.
```
ssh -X user@lxplus.cern.ch
```
In your /user or /work directory, do:
```
# If miniconda is not installed
# Get, install and activate miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh 
source miniconda3/bin/activate

# Get standard packages 
# (to have all spark functionalities pandas needs to be installed before pytimber)
pip install numpy scipy matplotlib ipython pandas

# Change python package index to CERN index
pip install git+https://gitlab.cern.ch/acc-co/devops/python/acc-py-pip-config.git

# Install pytimber
pip install pytimber

# Change python package index back to default
pip uninstall acc-py-pip-config
```
Test the installation with 
```
$ ipython
import pytimber
ldb = pytimber.LoggingDB(source="nxcals") 
ldb.search('LHC%BEAM_ENERGY%')
ldb.get(ldb.search('LHC%BEAM_ENERGY%')[0], t1='2022-06-15 15:10:30.0000')
```


