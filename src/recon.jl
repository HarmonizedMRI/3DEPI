"""
    recon(y, s, ϕ, U, V, Ω; λ, kwargs...)

Reconstruct an image from multicoil Cartesian k-space samples `y`,
given sensitivity maps `s`,
the complex response `ϕ` of a (possibly prephasing) RF pulse,
matrices `U` and `V` that incorporate information
about off-resonance frequency (gradients) (see [`getUV`](@ref)),
and k-space sampling mask `Ω`.

`λ` is a regularization parameter that controls the tradeoff
between data fidelity and spatial smoothness.

`kwargs` include other keyword arguments accepted by
[`MIRT.ncg`](https://github.com/JeffFessler/MIRT.jl/blob/master/src/algorithm/general/ncg.jl).

## Note
This method assumes Cartesian k-space sampling;
make sure `U` and `V` were generated accordingly.
"""
function recon(y, s, ϕ, U, V, Ω; kwargs...)

    A = [encodingmatrix(s, ϕ, U, V, Ω) for s in s]
    _recon(y, s, A, Ω; kwargs...)

end

function _recon(y, s, A, Ω; λ = 0, kwargs...)

    ncoils = length(s)
    T = diffop(size(Ω))
    cost = x -> sum(0.5 * norm(A[c] * x - y[c])^2 for c = 1:ncoils) + λ * norm(T * x)^2
    fun = (x, iter) -> cost(x)
    gradf = [[v -> v - y for y in y]; v -> 2λ * v]
    curvf = [fill(v -> 1, ncoils); v -> 2λ]
    B = [A; T]
    x0 = conjugate_phase(y, s, A, Ω)
    (x, out) = ncg(B, gradf, curvf, x0; fun, kwargs...)
    costs = [out[] for out in out]

    return (x, costs)

end

function conjugate_phase(y, s, A, Ω)

    ncoils = length(s)
    sos = sum(abs2.(s) for s in s)
    x = div0.(sum(A[c]' * y[c] for c = 1:ncoils), sos) ./ count(Ω)
    return x

end

function div0(a, b)

    tmp = a / b
    return iszero(b) ? zero(tmp) : tmp

end
