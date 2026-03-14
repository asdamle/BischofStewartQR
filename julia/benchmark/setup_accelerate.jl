if !Sys.isapple()
    error("Apple Accelerate setup is only relevant on macOS")
end

using AppleAccelerate
using LinearAlgebra

cfg = sprint(show, BLAS.get_config())
println("BLAS config after loading AppleAccelerate:")
println(cfg)

if !occursin("Accelerate", cfg)
    error("Apple Accelerate does not appear to be the active BLAS/LAPACK backend")
end

println("Apple Accelerate backend is active.")
