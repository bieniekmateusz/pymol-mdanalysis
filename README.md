# PyMOL-MDAnalysis exploratory project

This is an exploratory project started as a PyMOL Fellowship in 2019 
which combines PyMOL and MDAnalysis together. 

If you need any help, feel free to create an issue. 


## Installation

This fork uses a dumbed-down version of the installation file setup.py 
which might be lacking in various ways. 

If you have a problem with missing dependencies, 
check the official instructions: please see [INSTALL](INSTALL) 
and the parent github page. 


### CONDA
It is easiest to use conda to install some necessary dependencies.
```
conda install -c conda-forge glew glm libpng freetype libxml2 netcdf4
```

Download the repository, e.g. with git:
`git clone git@github.com:bieniekmateusz/pymol-mdanalysis.git`
then go to the `pymol-mdanalysis` directory and enter:
```
pip install .
```
In this case `pip` should refer to conda's pip, which you can check
with the command `where pip`.

From now on, whenever your conda environment is active, 
type `pymol` to start the application. 


### Linux

Installation follows the same process as Windows. 
Similarly to Windows, you might need to install 
glew, glm or other implementing binaries.  

An example on ubuntu:
```
apt-get install build-essential python3-dev libglew-dev \
  libpng-dev libfreetype6-dev libxml2-dev \
  libmsgpack-dev python3-pyqt5.qtopengl libglm-dev libnetcdf-dev
```

Then the pip/conda instructions are applicable. 


## Errors / Debugging

A good test is to use `python setup.py build` and check the errors in the output. 

In a conda environment, gcc version 10 was tested and works fine. However, gcc 11.2 led to the following error:
```
ImportError: ... _cmd.cpython-39-x86_64-linux-gnu.so: undefined symbol: _ZSt28__throw_bad_array_new_lengthv
```
You can switch to an older compiler in that case by configuirng CC and CXX:
```
export CC=gcc-10
export CXX=g++-10
```
Then clear the files from build/* and then reinstall. 


## Citations

If you find this useful, please cite the SI of Mateusz Bieniek's thesis:
https://librarysearch.kcl.ac.uk/permalink/f/1fdleu2/TN_proquest2475221287

Direct Thesis Download via KCLPure:
https://kclpure.kcl.ac.uk/portal/en/theses/fibronectin-iii910-adsorption-to-selfassembled-monolayers-and-interdomain-orientation-in-the-context-of-materialdriven-fibrillogenesis-studied-with-molecular-dynamics-simulations(11b1d705-79f0-4d0a-8c53-b68bc2c3338e).html
 
## Official repo info 

### Contributing

See [DEVELOPERS](DEVELOPERS).

### License

Copyright (c) [Schrodinger, LLC](https://www.schrodinger.com/)

Published under a BSD-like license, see [LICENSE](LICENSE).
