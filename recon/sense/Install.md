# Installation

We assume that you want to call the `recon` function provided
in this Julia package from some arbitrary folder on your computer.

To enable that we'll install the EPI3D package from a local folder
as well as any other packages you may need,
and record the state of all packages and dependencies in a Julia environment
(Project.toml and Manifest.toml files) for later use.

The instructions below assumes the EPI3D Julia package is located in
~/github/HarmonizedMRI/3DEPI/recon/sense/, and that the 
working directory is 3DEPI_example.


From your working directory:
1. Create an empty Project.toml file
```
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

4. Add any other packages you may need
``` 
(3DEPI_example) pkg> add MAT, JLD2, MIRTjim
``` 

4. Instantiate (create/update Manifest.toml)
```
(3DEPI_example) pkg> instantiate
```
This will create a file Manifest.toml that records the exact state of all
dependencies. This file (together with Project.toml) can later
be used to reproduce the exact same Julia environment.

4. Press `backspace` to exit the package manager and get back to the Julia prompt.

More info on packages and environments in Julia:
https://pkgdocs.julialang.org/v1/



## Workflow

1. Create `du.mat` containing undersampled data,
and do 3D SENSE recon in bart
```
>> main;
```

1. Reconstruct in Julia
```
pkg> activate .
(3DEPI_example) pkg> instantiate
julia> include("test2.jl")   # reconstruct real data
```
10 iterations takes about 17s on Jon's somewhat old HP EliteBook G3 
quad core laptop, Intel(R) Core(TM) i5-6300U CPU @ 2.40GHz.


1. View result (Matlab)
```
>> load result
>> se = norm(espiritreco(:));
>> sx = norm(x(:));
>> im(cat(1, espiritreco(:,:,30)/se, x(:,:,30)/sx))
```


## SMS/3D recon without B0 correction

The Julia function `recon` in the EPI3D package can be called as follows:
```
@time (x, costs) = recon(data, s, Ω; λ, niter)
```
where 
```
data        (nCoils,)  Acquired data (Cartesian)
data[1]     (Nk,)      Data for coil 1. Nk = number of acquired data samples,
                       e.g., in SMS EPI with no partial fourier, Nk = Nx*Ny.
s           (nCoils,)  Sensitivity maps
s[1]        (Nx, Ny, nSlices)   coil 1 sensitivity map
Ω           (Nx, Ny, nSlices)   Cartesian k-space sampling mask. sum(Ω) = Nk
```


## Test scripts

```
julia> include("test1.jl")   # numerical phantom
julia> include("test2.jl")   # real 3D data
```


## Add Julia support in Vim

```
cd ~/.vim
mkdir -p pack/plugins/start && cd pack/plugins/start
git clone https://github.com/JuliaEditorSupport/julia-vim.git
```
