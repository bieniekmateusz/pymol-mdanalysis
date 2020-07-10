# PyMOL-MDAnalysis exploratory project

This is an exploratory project started as a PyMOL Fellowship in 2019 
which combines PyMOL and MDAnalysis together. 

## Installation

This fork uses a dumbed-down version of the installation file setup.py 
which might be lacking in various ways. 

### Windows
It is easiest to use conda. Ensure you install some 
implementations of glew and glm. For example:
```
conda install -c menpo glew
conda install -c conda-forge glm
```

Use `git clone` or another way to download the repository,
go to the main directory and enter:
```
pip install .
```
In this case `pip` should refer to conda's pip. 
If you would like to use a startmenu shortcut PymolMda, 
use the following script:
```
python make_pymol_shortcut.py 
```
If you have a problem with missing dependencies, 
check the official instructions: please see [INSTALL](INSTALL) 
and the parent github page. 

### Linux

Use pip or conda. Conda follows the same format as Windows. 
Similarly to Window, you might need to install 
glew or glm implementing binaries.  

```
pip install .
```
If you would like a shortcut (added to the user space):
```
python make_pymol_shortcut.py
```
 

## Official repo info 

### Contributing

See [DEVELOPERS](DEVELOPERS).

### License

Copyright (c) [Schrodinger, LLC](https://www.schrodinger.com/)

Published under a BSD-like license, see [LICENSE](LICENSE).
