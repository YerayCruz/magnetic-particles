include("./main.jl")
using .MAIN_SIMULATION

payload = Dict("k" =>  0.5, ##pN/nm
"ka"  => 100.0, #pN/°
"ϵ" =>  100.0, #pN*nm
"σ"  => 1000.0, #nm
"d" =>  2700.0, #nm
"angles" =>  [180, 100], #°
"bonds" =>  [6, 3],
"B" => [0.0, 1.0], #mT
"viscocity" =>  9.0, #pN*s/nm
"wall_dimension" =>  20000.0, #nm
"T" => 300.0 #K
)

brownian_motion(payload)
