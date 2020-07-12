# PyMOL-MDAnalysis exploratory project

This is an exploratory project started as a PyMOL Fellowship in 2019 
which combines PyMOL and MDAnalysis together. 

## Installation

This fork uses a dumbed-down version of the installation file setup.py 
which might be lacking in various ways. 

### Windows
It is easiest to use conda to install some necessary dependencies.
```
conda install -c conda-forge glew glm libpng freetype libxml2
```

Download the repository, e.g. with git:
`git clone git@github.com:bieniekmateusz/pymol-mdanalysis.git`
then go to the `pymol-mdanalysis` directory and enter:
```
pip install .
```
In this case `pip` should refer to conda's pip, which you can check
with the command `where pip`.

If you would like to use a startmenu shortcut PymolMda, after
installation is complete, use the following script:
```
python make_pymol_shortcut.py 
```

If you have a problem with missing dependencies, 
check the official instructions: please see [INSTALL](INSTALL) 
and the parent github page. 

### Linux

Installation follows the same process as Windows. 
Similarly to Windows, you might need to install 
glew, glm or other implementing binaries.  

```
pip install .
```
If you would like a shortcut (added to the user space):
```
python make_pymol_shortcut.py
```

### MacOS

Installation on MacOS is not tested but should also be
possible via the same process described for Windows and
Linux.
 
## Official repo info 

### Contributing

See [DEVELOPERS](DEVELOPERS).

### License

Copyright (c) [Schrodinger, LLC](https://www.schrodinger.com/)

Published under a BSD-like license, see [LICENSE](LICENSE).
