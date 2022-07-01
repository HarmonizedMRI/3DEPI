"""
    phasecorrect!(y, z_start, nz, Δ, Ω)

Correct the phase of multicoil k-space data `y`
for a set of simultaneously excited slices.
Such correction is needed because the center
of the set of simultaneously excited slices
does not necessarily align with isocenter
(i.e., the center of the full acquired volume).

# Arguments
- `y`: Multicoil k-space data
- `z_start`: z-index of first simultaneously excited slice
  (e.g., if slices 2, 5, and 8 are simultaneously excited,
  then `z_start = 2`)
- `nz`: Total number of slices in the imaging volume
  (*not* the number of simultaneously excited slices)
- `Δ`: Image resolution [Δx, Δy, Δz]
  (in particular, `Δ[3]` is *not* the interslice spacing
  of the simultaneously excited slices)
- `Ω`: SMS k-space sampling mask
"""
function phasecorrect!(y, z_start, nz, Δ, Ω)

    nexcited = size(Ω, 3) # Number of simultaneously excited slices
    nexcitations = nz ÷ nexcited # Number of sets of simultaneously excited slices
    nexcitations * nexcited == nz || error("number of simultaneously excited \
        slices (size(Ω, 3) = $nexcited) should evenly divide the total number \
        of slices ($nz)")
    1 <= z_start <= nexcitations || error("the starting z index ($z_start) \
        should be between 1 and the number of sets of simultaneously excited \
        slices (nz / nexcited = $nexcitations), inclusive")

    z_iso = nz ÷ 2 + 1 # Center of full volume
    # Center of set of simultaneously excited slices
    z_excited = z_start + nexcitations * (nexcited ÷ 2)
    z_offset = (z_excited - z_iso) * Δ[3]

    # k-space sample spacing
    Δk = (
        1 / (size(Ω, 1) * Δ[1]), # 1 / FOV
        1 / (size(Ω, 2) * Δ[2]), # 1 / FOV
        1 / (nz * Δ[3]) # 1 / FOV
    )

    # k-space sample locations
    kdim = d -> (-size(Ω, d)÷2:size(Ω, d)÷2-iseven(size(Ω, d))) * Δk[d]
    kx = kdim(1)
    ky = kdim(2)
    kz = kdim(3)
    kz = [kz for kx in kx, ky in ky, kz in kz][Ω]

    # Do phase correction for each coil
    for d in y
        d .*= cispi.(2 .* z_offset .* kz)
    end

    return y

end
