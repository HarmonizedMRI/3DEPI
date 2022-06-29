"""
    signalmodel(x, s, ϕ, r, tm, km, f, g, γ, Δ)

Compute single-coil k-space data
for a given k-space location `km`
from image `x`
given the coil's sensitivity map `s`
and the complex RF response `ϕ`.

`r` is a collection of positions [x, y, z]
such that `r[n]` is the position of image sample `x[n]`.

See [`getWmn`](@ref) for a description of the other method arguments.
"""
function signalmodel(x, s, ϕ, r, tm, km, f, g, γ, Δ)

    dm = sum(s[n] * ϕ[n] * cispi(-2 * dot(km, r[n])) * x[n] *
             getWmn(tm, km, f[n], g[n], γ[n], Δ) for n in eachindex(x))
    return dm

end

"""
    signalmodel(x, s, r, km)

Compute k-space data ignoring B0 effects.
"""
function signalmodel(x, s, r, km)

    dm = sum(s[n] * cispi(-2 * dot(km, r[n])) * x[n] for n in eachindex(x))
    return dm

end

"""
    forwardmodel(x, s, ϕ, U, V, Ω)

Compute single coil k-space data from image `x`
given the coil's sensitivty map `s`,
the complex RF response `ϕ`,
and matrices `U` and `V` that incorporate information
about off-resonance frequency (gradients) (see [`getUV`](@ref)).

`Ω` is the k-space sampling mask.
`count(Ω) == size(U, 1)` should hold true.

Use [`adjointmodel`](@ref) to compute the adjoint of this operation.

## Note
This method assumes Cartesian k-space sampling;
make sure `U` and `V` were generated accordingly.
"""
function forwardmodel(x, s, ϕ, U, V, Ω)

    ϕsx = ϕ .* s .* x
    vϕsx = zeros(promote_type(eltype(ϕsx), eltype(V)), size(Ω))
    y = sum(1:size(U, 2)) do l
        @views vϕsx .= reshape(V[l,:], size(Ω)) .* ϕsx
        @views U[:,l] .* ffts(vϕsx)[Ω]
    end
    return y

end

"""
    adjointmodel(y, s, ϕ, U, V, Ω)

Compute single coil image data from k-space data `y`.

See [`forwardmodel`](@ref) for a description of the method arguments.
"""
function adjointmodel(y, s, ϕ, U, V, Ω)

    uy = zeros(promote_type(eltype(y), eltype(U)), size(Ω))
    sumv = sum(1:size(U, 2)) do l
        @views uy[Ω] = conj.(U[:,l]) .* y
        @views reshape(conj.(V[l,:]), size(Ω)) .* bffts(uy)
    end
    x = conj.(s) .* conj.(ϕ) .* sumv
    return x

end

"""
    forwardmodel(x, s, Ω)

Compute Cartesian k-space data ignoring B0 effects.
"""

function forwardmodel(x, s, Ω)

    y = ffts(s .* x)[Ω]
    return y

end

"""
    adjointmodel(y, s, Ω)

Compute image data ignoring B0 effects.
"""
function adjointmodel(y, s, Ω)

    y0 = zeros(eltype(y), size(Ω))
    y0[Ω] = y
    x = conj.(s) .* bffts(y0)
    return x

end

ffts(X) = fftshift(fft(ifftshift(X)))
bffts(Y) = fftshift(bfft(ifftshift(Y)))

encodingmatrix(s, ϕ, U, V, Ω) = LinearMapAA(
    x -> forwardmodel(x, s, ϕ, U, V, Ω),
    y -> adjointmodel(y, s, ϕ, U, V, Ω),
    (count(Ω), length(Ω));
    idim = size(Ω),
    T = ComplexF64
)

encodingmatrix(s, Ω) = LinearMapAA(
    x -> forwardmodel(x, s, Ω),
    y -> adjointmodel(y, s, Ω),
    (count(Ω), length(Ω));
    idim = size(Ω),
    T = ComplexF64
)
