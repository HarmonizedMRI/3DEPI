# SMS/3D EPI SENSE recon with/without B0 correction

Can also correct for through-voxel B0 gradients.  


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

## Installation
[Install.md](Install.md)


## Example workflow

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
