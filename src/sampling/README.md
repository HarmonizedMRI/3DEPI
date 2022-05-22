
## Generate 2D CAIPI sampling patterns (python)

### Get submodule

The `skippedcaipi` repository (https://github.com/mrphysics-bonn/skipped-caipi)
is included here as a git submodule.
The contents of that submodule must be downloaded explicitly 
using the ``--recursive`` git directive:

If you've already cloned this 3DEPI repository, enter into the 3DEPI repository folder and do
```
git submodule update --init --recursive
```

If cloning this 3DEPI repository for the first time, do:
```
$ git clone --recursive 
```

### Run test

Generate CAIPI pattern and plot.

To run from Linux command line:
```
$ python test.py
```

From Python prompt (REPL):
```
>>> import test
```
To run again:
```
>>> import importlib
>>> importlib.reload(test)
```


