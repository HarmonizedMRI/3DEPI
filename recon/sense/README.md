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

Processing time:
For 6-fold undersampling, 6/8 partial Fourier, and matrix size [80 80 60],
10 iterations take about 17s on an HP EliteBook G3 
quad core laptop, Intel(R) Core(TM) i5-6300U CPU @ 2.40GHz.
The result is visually similar to SENSE reconstruction with BART, which
takes about 6s.



## Setup

To install this Julia package, see
[Setup.md](Setup.md).


