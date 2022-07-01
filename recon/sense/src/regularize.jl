# Basically the same as MIRT.diff_map, but works with complex data and
# non-vector inputs
function diffop(N::Dims{D}; dims = 1:D) where {D}
    diff_check(N, dims)
    return LinearMapAA(
        x -> diff_forw(x; dims),
        d -> diff_adj(d, N; dims),
        (sum(diff_length(N, dim) for dim in dims), prod(N)),
        (; name = "diffop", dims);
        idim = N,
        T = ComplexF64
    )
end
