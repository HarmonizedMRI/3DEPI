# Design CAIPI sampling pattern and write to .mat file

Python code for designing Cartesian ky-kz undersampling 
following the CAIPI pattern.
The sampling pattern is written to a .mat file.

See https://github.com/mrphysics-bonn/skipped-caipi.


## Usage example

Design CAIPI pattern suitable for SMS fMRI, with the following parameters:  
ny = 64  (y matrix size)  
nz = 6   (z matrix size)  
Ry = 1  
Rz = 6  (multiband factor)  
Dz = 2  (CAIPI shift along z)  
S = 1  (segmentation factor)

From linux command line:
```
$ python skippedcaipi_sampling.py 64 6 1 6 2 1
```
This creates the file `caipi.mat` that can be loaded into MATLAB.

```
caipi.mat:
indices   [64 2]  sampling indeces for one shot
mask      [64 6]  sampling pattern 
```


## Run test

Generate CAIPI pattern and plot.
To run from Linux command line:
```
$ python caipitest.py
```

From Python prompt (REPL):
```
>>> import caipitest
>>> # To run again in the same session:
>>> import importlib
>>> importlib.reload(caipitest)
```

