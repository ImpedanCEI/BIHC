# copyright ################################# #
# This file is part of the BIHC Package.      #
# Copyright (c) CERN, 2024.                   #
# ########################################### #

from setuptools import setup, find_packages
from pathlib import Path

# read version
version_file = Path(__file__).parent / 'bihc/_version.py'
dd = {}
with open(version_file.absolute(), 'r') as fp:
    exec(fp.read(), dd)
__version__ = dd['__version__']

# read long_description
long_description = (Path(__file__).parent / "README.md").read_text(encoding="utf-8")

#########
# Setup #
#########

setup(
    name="bihc",
    version=__version__,
    description="BIHC: Beam-Induced Heating Computation package",
    author="Elena de la Fuente, Francesco Giordano, Leonardo Sito",
    author_email="elena.de.la.fuente.garcia@cern.ch", 
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    url="https://github.com/ImpedanCEI/BIHC",
    project_urls={"Bug Tracker": "https://github.com/ImpedanCEI/BIHC/issues"},
    install_requires = [
                    'numpy',
                    'matplotlib',
                    'scipy',
                       ],
    classifiers=[
        "Development Status :: 4 - Beta",
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
		"Topic :: Scientific/Engineering :: Physics",
    ],

)
