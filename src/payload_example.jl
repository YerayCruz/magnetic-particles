include("./main.jl")
using .MAIN_SIMULATION

payload = Dict("k" =>  500.0, ##pN/nm
"ka"  => 1000.0, #pN/°
"ϵ" =>  1000.0, #pN*nm
"σ"  => 10.0, #μm
"d" =>  20.7, #μm
"angles" =>  [180, 100], #°
"bonds" =>  [6, 3],
"B" => [0.0, 1.0], #mT
"viscocity" =>  9.0, #pN*s/nm
"wall_dimension" =>  200.0, #μm
"T" => 300.0 #K
)

brownian_motion(payload)
