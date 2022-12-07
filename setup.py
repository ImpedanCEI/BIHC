from setuptools import setup

####################################
# Add long description from README #
####################################

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


#########
# Setup #
#########

setup(
    name="bihc",
    version="0.0.1",
    description="BIHC: Beam-Induced Heating Computation package",
    author="Francesco Giordano",
    author_email="francesco.giordano@cern.ch",
    packages=['bihc'],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LeonardoSito/BIHC",
    project_urls={"Bug Tracker": "https://github.com/LeonardoSito/BIHC/issues"},
    install_requires = [
					'numpy',
					'matplotlib',
					'scipy',
					'pandas',
					'pytimber',
						]
    ],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],

)
