# BIHC
Beam Induced Heating Computation (BIHC) tool is a package that allows the estimation of the dissipated power due to the passage of a particle beam inside an accelerator component.

The dissipated power value depends on the characteristics of the particle beam (beam spectrum and intensity) and on the characteristics of the consdiered accelerator component (beam-coupling impedance).

Check :file_folder: examples/ on how to use!

Documentation is avaiable in [bihc.readthedocs.io](https://bihc.readthedocs.io/en/latest/)

**Fist release coming soon** (January 2023)

## Installation
This section explains how to set up the environment to start using BIHC package for power loss computations.

### Installation in CERN lxplus

Connect to CERN lxplus via ssh. Avoid connecting to lxplus8, the code will induce in Kerberos issues. Kerberos logging will expire 4h after each connection and needs to be renewed.
```
ssh -XY user@lxplus.cern.ch
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
#### Setup Git in lxplus
```
git config --list
```
Look for user.name and user.email. If it is not yet set, run the following commands (which set the information for all your git repositories on lxplus)
```
git config --global user.name "Your Name"
git config --global user.email "your.name@cern.ch"
```
There are other settings recommended as well:
```
git config --global push.default simple
git config --global http.postBuffer 1048576000
git config --global http.emptyAuth true # Required on CC7
```
The push setting makes some operations more straightforward. The second addresses an issue with large pushes via plain http or krb5. The third addresses an issue with libcurl and krb5.

#### Developers: Download BIHC repository
```
git clone https://github.com/LeonardoSito/BIHC.git
```
#### Users: pip install 
(Available soon)
```
pip install bihc
```
