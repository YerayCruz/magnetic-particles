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

function harmonic_force(r1, r2, angle, ka)
  r1_norm = sqrt(sum(r1 .^ 2))
  r2_norm = sqrt(sum(r2 .^ 2))
  factor = dot(r1, r2)/(r1_norm*r2_norm)
  if factor > 1 && factor < 1.05
    factor = 1
  elseif factor < -1 && factor > -1.05
    factor = -1
  end
  current_angle = acosd(factor)
  absolute_angle = acosd(r1[1]/r1_norm)
  if r1[2] < 0
    absolute_angle = 360 - absolute_angle
  end
  magnitude = -ka*(current_angle - angle)/r1_norm
  force_direction = [-sind(absolute_angle), cosd(absolute_angle)]
  println(current_angle)
  return force_direction * magnitude
end

function dipole_magnetic_force()

end

end
