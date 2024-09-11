module POTENTIALS
using LinearAlgebra

export F_lj_2d, elastic_f, harmonic_force, wca_f, dipole_magnetic_force

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

function wca_f(r, ϵ, σ)
    r_norm = sqrt(sum(r .^ 2))
    force_direction = r / r_norm
    magnitude = 0
    if r_norm < 2 ^ (1 / 6) * σ
        magnitude = 4 * ϵ * (12 * (σ / r_norm) ^ 12 - 6 * (σ / r_norm) ^ 6) / r_norm + ϵ
    end
    return magnitude * force_direction
end
function harmonic_force(r1, r2, angle, ka)
    r1_norm = sqrt(sum(r1 .^ 2))
    r2_norm = sqrt(sum(r2 .^ 2))
    r1_buffer = copy(r1)/r1_norm
    r2_buffer = copy(r2)/r2_norm
    push!(r1_buffer, 0.)
    push!(r2_buffer, 0.)
    force_direction = (cross(cross(r2_buffer, r1_buffer), r1_buffer))
    factor = dot(r1, r2) / (r1_norm * r2_norm)
    if factor > 1. && factor < 1.01
        factor = 1
    end
    if factor < -1. && factor > -1.1
        factor = -1
    end
    current_angle = acosd(factor)
    if current_angle == 0 || current_angle == 180
        force_direction = [0. , 1., 0.]
    end
    magnitude = -ka * (current_angle - angle)/r1_norm
    pop!(force_direction)
    return force_direction .* magnitude
end

function dipole_magnetic_force(r, m1, m2)
    μ = 4 * π * 10^(-7)
    r_norm = sqrt(sum(r .^ 2))
    force_direction = r / r_norm

    m1_dot_r = dot(m1, r)
    m2_dot_r = dot(m2, r)
    m1_dot_m2 = dot(m1, m2)

    prefactor = ( 3 * μ ) / ( 4 * π * r_norm ^ 3)
    magnitude = prefactor * ((m1_dot_r * m2) + (m2_dot_r * m1) + (m1_dot_m2 * r) - (5 * m1_dot_r * m2_dot_r)/(r_norm .^ 2) * r)
    return magnitude
end

end
