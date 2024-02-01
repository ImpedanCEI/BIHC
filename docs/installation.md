# Installation guide

This section explains how to set up the environment to start using BIHC package for power loss computations. For installing `pytimber` check the section *Installation with pytimber in CERN lxplus*. `pytimber` is also accesible from CERN SWAN python notebooks using the `102b NXCALS PRO` configuration.

<div class="warning">
<div class="admonition-title">
Warning
</div>
This documentation is currently under development
</div>

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



