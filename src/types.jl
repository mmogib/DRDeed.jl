abstract type ModelData end
abstract type ModelSolution end
abstract type ModelResult end

struct SlimData <: ModelData
  customers::Int
  generators::Int
  periods::Int
  Caith::Matrix{Float64}
  Caitp::Matrix{Float64} # rows: customers, cols: periods
  Caits::Matrix{Float64}
  Caitf::Matrix{Float64}
  Cajtp::Matrix{Float64}
  Cajts::Matrix{Float64}
  Caijtx::Array{Float64,3}
  Caijty::Array{Float64,3}
  DemandDt::Matrix{Float64}
  CProductionLimitst::Matrix{Float64} # Lit
  CPowerUpperLimitst::Matrix{Float64} # Kit
  CStorageUpperLimitst::Matrix{Float64} # Cit
  PProductionLimitst::Matrix{Float64} # Ljt
  PPowerUpperLimitst::Matrix{Float64} # Kjt
  PStorageUpperLimitst::Matrix{Float64} # Cjt
end

struct SlimSolution <: ModelSolution
  U::Vector{Float64}
  V::Vector{Float64}
  x::Array{Float64,3}
  y::Array{Float64,3}
  p::Matrix{Float64}
  q::Matrix{Float64}
  f::Matrix{Float64}
  h::Matrix{Float64}
  si::Matrix{Float64}
  sj::Matrix{Float64}
end

# Deed
struct DeedData <: ModelData
  customers::Int
  generators::Int
  periods::Int
  Dt::Vector{Float64}
  B::Array{Float64,2}
  pjmin::Vector{Float64}
  pjmax::Vector{Float64}
  DR::Vector{Float64}
  UR::Vector{Float64}
  K1::Vector{Float64}
  K2::Vector{Float64}
  θ::Vector{Float64}
  CM::Vector{Float64}
  UB::Float64
  a::Vector{Float64}
  b::Vector{Float64}
  c::Vector{Float64}
  e::Vector{Float64}
  f::Vector{Float64}
  g::Vector{Float64}
  λ::Array{Float64,2}
end


struct DeedSolution <: ModelSolution
  P::Array{Float64,2}
  Cost::Float64
  Emission::Float64
  Loss::Vector{Float64}
  w::Vector{Float64}
end

Base.show(io::IO, ddata::DeedSolution) = println(
  io,
  """Deed Problem: solution
    Cost: $(ddata.Cost)
    Emission: $(ddata.Emission)
    Loss: $(ddata.Loss)
""",
)

struct DRDeedSolution <: ModelSolution
  χ::Array{Float64,2}
  ω::Array{Float64,2}
  q::Array{Float64,2}
  Cost::Float64
  Emission::Float64
  Utility::Float64
  Losst::Vector{Float64}
  w::Vector{Float64}
end

Base.show(io::IO, ddata::DRDeedSolution) = println(
  io,
  """Deed Problem: solution
    Cost: $(ddata.Cost)
    Emission: $(ddata.Emission)
    Utility: $(ddata.Utility)
    Loss: $(sum(ddata.Losst))
""",
)


# GameTheory Deed
struct GTDeedData <: ModelData
  deedData::DeedData
  ahdot_it::Matrix{Float64}
  adot_it::Matrix{Float64}
  bdot_it::Matrix{Float64}
  cdot_it::Matrix{Float64}
  ddot_it::Matrix{Float64}
  edot_it::Matrix{Float64}
  nudot_it::Matrix{Float64}
  afdot_it::Matrix{Float64}
  CDemandt::Matrix{Float64}
  pimin::Vector{Float64}
  pimax::Vector{Float64}
end
struct GTDRDeedSolution <: ModelSolution
  deedSolution::DRDeedSolution
  p::Matrix{Float64}
  s::Matrix{Float64}
  h::Matrix{Float64}
  f::Matrix{Float64}
  x::Matrix{Float64}
  y::Array{Float64,3}
  customerCost::Matrix{Float64}
end


## Results

struct SuccessResult <: ModelResult
  data::ModelData
  solution::ModelSolution
end


struct FailResult <: ModelResult
  data::ModelData
  messga::String
end