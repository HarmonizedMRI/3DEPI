# Design CAIPI sampling pattern and write to .mat file

## Usage example

From linux command line:
```
$ python skippedcaipi_sampling.py 1 6 2 1 64 30
```
This creates the file `caipi.mat` that can be loaded into MATLAB.


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

