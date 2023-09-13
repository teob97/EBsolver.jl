using EBsolver
using Test, Statistics

@testset "BackgroundCosmology test" begin
  include("BackgroundCosmology.jl")
end

@testset "RecombinationHistory test" begin
  include("RecombinationHistory.jl")
end