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

struct SENSEWorkspace{T,N,P1,P2}
    x::Array{T,N}
    y::Vector{T}
    sx::Array{T,N}
    yfull::Array{T,N}
    fft_plan::P1
    bfft_plan::P2
    shift::Array{T,N}
end

function SENSEWorkspace(s, Ω)

    x = similar(s)
    y = similar(s, count(Ω))
    sx = similar(s)
    yfull = similar(s)
    fft_plan = plan_fft(s)
    bfft_plan = plan_bfft(s)
    shift = similar(s)
    return SENSEWorkspace(x, y, sx, yfull, fft_plan, bfft_plan, shift)

end

"""
    forwardmodel(x, s, Ω)

Compute Cartesian k-space data ignoring B0 effects.
"""

function forwardmodel(x, s, Ω)

    y = ffts(s .* x)[Ω]
    return y

end

function forwardmodel(x, s, Ω, workspace)

    workspace.sx .= s .* x

    # Compute `workspace.yfull = fftshift(fft(ifftshift(workspace.sx)))`
    ffts!(workspace.yfull, workspace.sx, workspace.fft_plan, workspace)

    # The following for loop computes
    # `workspace.y .= workspace.yfull[Ω]`
    # without allocating a temporary array for workspace.yfull[Ω]
    i = 0
    for j in eachindex(Ω, workspace.yfull)
        if Ω[j]
            i += 1
            workspace.y[i] = workspace.yfull[j]
        end
    end

    return workspace.y

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

function adjointmodel(y, s, Ω, workspace)

    # Zero-fill the k-space data
    i = 0
    for j in eachindex(Ω, workspace.yfull)
        if Ω[j]
            i += 1
            workspace.yfull[j] = y[i]
        else
            workspace.yfull[j] = 0
        end
    end

    # Compute `workspace.sx = fftshift(bfft(ifftshift(workspace.yfull)))`
    ffts!(workspace.sx, workspace.yfull, workspace.bfft_plan, workspace)

    workspace.x .= conj.(s) .* workspace.sx

    return workspace.x

end

ffts(X) = fftshift(fft(ifftshift(X)))
bffts(Y) = fftshift(bfft(ifftshift(Y)))

function ffts!(Y, X, plan, workspace)

    # Compute `Y = fftshift(fft(ifftshift(X)))`
    # (`fft` could be `bfft`, `ifft`, etc. depending on `plan`)
    ifftshift!(workspace.shift, X)
    mul!(Y, plan, workspace.shift)
    fftshift!(workspace.shift, Y)
    Y .= workspace.shift

end

encodingmatrix(s, ϕ, U, V, Ω) = LinearMapAA(
    x -> forwardmodel(x, s, ϕ, U, V, Ω),
    y -> adjointmodel(y, s, ϕ, U, V, Ω),
    (count(Ω), length(Ω));
    idim = size(Ω),
    T = ComplexF64
)

function encodingmatrix(s, Ω)

    workspace = SENSEWorkspace(s, Ω)
    return LinearMapAA(
        x -> forwardmodel(x, s, Ω, workspace),
        y -> adjointmodel(y, s, Ω, workspace),
        (count(Ω), length(Ω));
        idim = size(Ω),
        T = ComplexF64
    )

end
