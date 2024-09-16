include("./main.jl")
using .MAIN_SIMULATION

payload = Dict("k" =>  0.5,
"ka"  => 0.1,
"ϵ" =>  1.0,
"σ"  => 1.0,
"d" =>  2.7,
"angles" =>  [180, 100],
"bonds" =>  [6, 3],
"m" => [[0.0, 1.0], [0.0, 1.0]],
"viscocity" =>  0.9,
"wall_dimension" =>  20.0)

brownian_motion(payload)
