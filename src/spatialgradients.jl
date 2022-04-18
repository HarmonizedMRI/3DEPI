"""
    spatialgradients(x, r)

Compute the spatial (x, y, z) gradients of `x`.
(`x` is, e.g., a B0 field map.)
`r` is a map of spatial positions,
i.e., `r[n]` is the (x, y, z) position of `x[n]`.

The output is a map of spatial gradients,
i.e., `g[n]` contains the x, y, and z spatial gradients of `x`
at position `r[n]`.
"""
function spatialgradients(x, r)

    oneindex = zeros(Int, ndims(x))

    return map(CartesianIndices(x)) do i

        ntuple(ndims(x)) do d

            oneindex[d] = 1
            onecart = CartesianIndex(oneindex...)
            oneindex[d] = 0

            if i[d] == 1 # Forward difference
                i1 = i
                i2 = i + onecart
            elseif i[d] == size(x, d) # Backward difference
                i1 = i - onecart
                i2 = i
            else # Central difference
                i1 = i - onecart
                i2 = i + onecart
            end
            g = (x[i2] - x[i1]) / (r[i2][d] - r[i1][d])

        end

    end

end
