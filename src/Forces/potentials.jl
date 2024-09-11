module POTENTIALS
using LinearAlgebra

export F_lj_2d, elastic_f, harmonic_force

function elastic_f(r, k, d)
    r_norm = sqrt(sum(r .^ 2))
    force_direction = r / r_norm
    magnitude = -k*(r_norm-d)
    return force_direction * magnitude
end

function F_lj_2d(r, ϵ, σ)
    r_norm = sqrt(sum(r .^ 2))
    force_direction = r / r_norm
    magnitude = 4*ϵ * (12 * (σ / r_norm) ^ 12 - 6 * (σ / r_norm) ^ 6) / r_norm
    return force_direction * magnitude
end

f = open("./data.txt", "w")
function harmonic_force(r1, r2, angle, ka)
    r1_norm = sqrt(sum(r1 .^ 2))
    r2_norm = sqrt(sum(r2 .^ 2))
    r1_buffer = copy(r1)/r1_norm
    r2_buffer = copy(r2)/r2_norm
    push!(r1_buffer, 0.)
    push!(r2_buffer, 0.)
    force_direction = (cross(cross(r2_buffer, r1_buffer), r1_buffer))
    current_angle = acosd(dot(r1, r2) ./ (r1_norm * r2_norm))
    write(f, "force direction --> $force_direction, r1: $r1_buffer, r2: $r2_buffer, angle: $current_angle\n")
    magnitude = -ka * (current_angle - angle)/r1_norm
    pop!(force_direction)
    return force_direction .* magnitude
end

function dipole_magnetic_force()

end

function dipole_magnetic_force()

end

end
