"""
    getWmn(tm, km, fn, gn, γn, Δ)

Compute the (m,n)th entry of a weighting matrix W
that incorporates off-resonance frequency,
the off-resonance frequency gradient,
and the gradient of phase imparted by a prephasing RF pulse.

# Arguments
- `tm`: Time (s)
- `km`: K-space location [kx, ky, kz] (cycles/cm)
- `fn`: Off-resonance frequency (cycles/s)
- `gn`: Off-resonance frequency gradient [gx, gy, gz] (cycles/s/cm)
- `γn`: Gradient of phase imparted by prephasing RF pulse [γx, γy, γz] (cycles/cm)
- `Δ`: Voxel size [Δx, Δy, Δz] (cm)

## Note
- `m` indexes k-space samples, `n` indexes spatial samples
- `km`, `gn`, `γn`, and `Δ` all must have the same number of dimensions
  (typically 2D or 3D)
- The units for time and distance
  of the arguments don't necessarily have to match those given above,
  but they do need to be consistent
  (e.g., don't have `km` be cycles/cm and `γn` be cycles/mm;
  both need to be cycles/cm or cycles/mm or cycles/m or whatever)
- The units for angle *do* need to be cycles (not rad or °)
"""
function getWmn(tm, km, fn, gn, γn, Δ)

    Omn = cispi(-2 * fn * tm)
    Gmn = prod(sinc((km[i] + gn[i] * tm + γn[i]) * Δ[i]) for i = 1:length(km))
    return Omn * Gmn

end

abstract type AbstractUVMethod end

struct ExactUV <: AbstractUVMethod
    L::Int
end

struct SketchU <: AbstractUVMethod
    L::Int
    Usamp::Int
end

struct SketchUV <: AbstractUVMethod
    L::Int
    Usamp::Int
    Vsamp::Int
end

"""
    getUV(method, t, k, f, g, γ, Δ; weights)

Compute `U` and `V` such that `W ≈ U * V`.
See [`getWmn`](@ref).

# Arguments
- `method`:
  - `::ExactUV`: Compute `U` and `V` by taking the first `method.L` singular
    vectors and values of the full matrix `W` (impractical for large problems)
  - `::SketchU`: Compute `U` approximately by taking the first `method.L` left
    singular vectors of a matrix consisting of `method.Usamp` random columns of
    `W`; then compute `V = U' * W`
  - `::SketchUV`: Compute `U` approximately by taking the first `method.L` left
    singular vectors and values of a matrix consisting of `method.Usamp` random
    columns of `W`; then compute `V` as the minimizer of ``||W - UV||_F`` but
    with only `method.Vsamp` randomly selected rows of `W` (and the
    corresponding rows of `U`)
- `t`: Vector of times indicating when the corresponding k-space data points
  were sampled
- `k`: Vector of [kx, ky, kz] k-space locations
- `f`: Vector of off-resonance frequencies at each spatial location
- `g`: Vector of [gx, gy, gz] off-resonance frequency gradients at each spatial
  location
- `γ`: Vector of [γx, γy, γz] prephasing RF pulse phase gradients at each
  spatial location
- `Δ`: Voxel size [Δx, Δy, Δz]
- `weights = trues(size(f))`: Weights (typically binary) to, e.g.,
  ignore background voxels
"""
function getUV(method::ExactUV, t, k, f, g, γ, Δ; weights = trues(size(f)))

    M = length(t)
    N = length(f)

    W = [getWmn(t[m], k[m], f[n], g[n], γ[n], Δ) for m = 1:M, n = 1:N]

    F = svd(W .* transpose(vec(weights)))

    U = F.U[:,1:method.L]
    V = U' * W

    return (U, V)

end

function getUV(method::SketchU, t, k, f, g, γ, Δ; weights = trues(size(f)))

    M = length(t)
    N = length(f)

    U = sketchU(method, t, k, f, g, γ, Δ, weights)
    V = [sum(conj(U[m,l]) * getWmn(t[m], k[m], f[n], g[n], γ[n], Δ) for m = 1:M) for l = 1:method.L, n = 1:N]

    return (U, V)

end

function getUV(method::SketchUV, t, k, f, g, γ, Δ; weights = trues(size(f)))

    M = length(t)
    N = length(f)

    U = sketchU(method, t, k, f, g, γ, Δ, weights)

    idx = randperm(M)[1:method.Vsamp]
    W = [getWmn(t[m], k[m], f[n], g[n], γ[n], Δ) for m in idx, n = 1:N]

    V = U[idx,:] \ W # Solves min ||W - UV||_F

    return (U, V)

end

function sketchU(method, t, k, f, g, γ, Δ, weights)

    M = length(t)
    N = length(f)

    idx = randperm(N)[1:method.Usamp]
    W = [getWmn(t[m], k[m], f[n], g[n], γ[n], Δ) * weights[n] for m = 1:M, n in idx]

    F = svd(W)

    U = F.U[:,1:method.L]

    return U

end
