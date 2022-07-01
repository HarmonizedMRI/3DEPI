# Setup

Here we assume that you want to call the `recon` function provided
in this Julia package from some arbitrary folder on your computer.

To enable that we'll install the EPI3D package from a local folder
as well as any other packages you may need.
We'll also record the state of all packages and dependencies in a Julia environment
for later use; this information will be encapsulated
in two files residing in your working directory: 
'Project.toml' and 'Manifest.toml'.

The instructions below assumes the EPI3D Julia package is located in
~/github/HarmonizedMRI/3DEPI/recon/sense/, and that the 
working directory is 3DEPI_example.

From your working directory:
1. Create an empty Project.toml file
```
$ cd 3DEPI_example     # or wherever you want to work
$ touch Project.toml
```

3. Start Julia and press `]` to enter the Julia package manager

4. Activate the (empty) environment and add the EPI3D package
```
(@v1.7) pkg> activate .       # Set environment specified in Project.toml
(3DEPI_example) pkg> dev ~/github/HarmonizedMRI/3DEPI/recon/sense/  # adds local package
```
If you later wish to remove the EPI3D package, do:
```
(3DEPI_example) pkg> rm EPI3D 
```

4. Add any other packages you may need, e.g.,
``` 
(3DEPI_example) pkg> add MAT, JLD2, MIRTjim
``` 

4. Instantiate (create/update Manifest.toml)
```
(3DEPI_example) pkg> instantiate
```
This will create the file Manifest.toml that records the exact state of all
dependencies. This file (together with Project.toml) can later
be used to reproduce the exact same Julia environment.

4. Press `backspace` to exit the package manager and get back to the Julia prompt.

More info on packages and environments in Julia:
https://pkgdocs.julialang.org/v1/



## Julia tips

### Add Julia support in Vim
```
cd ~/.vim
mkdir -p pack/plugins/start && cd pack/plugins/start
git clone https://github.com/JuliaEditorSupport/julia-vim.git
```
