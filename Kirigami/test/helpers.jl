# helpers.jl -- shared test fixtures.
using Kirigami, Test, LinearAlgebra
const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const CORPUS = joinpath(REPO, "data", "corpus")
