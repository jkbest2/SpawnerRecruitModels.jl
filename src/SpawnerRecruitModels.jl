module SpawnerRecruitModels

export
    AbstractSpawnerRecruitModel,
    BevertonHoltSpawnerRecruitModel,
    RickerSpawnerRecruitModel,
    LudwigWaltersSpawnerRecruitModel,
    CushingSpawnerRecruitModel,
    DerisoSchnuteSpawnerRecruitModel,
    ShepherdSpawnerRecruitModel,
    GammaSpawnerRecruitModel,
    recruit,
    max_recruits,
    max_spawnrecruits

abstract type AbstractSpawnerRecruitModel end

"""
    BevertonHoltSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModel

Parameters:

- α: density-independent productivity parameter
- β: density-dependence parameter

Follows equation 3.6 (p. 88) of Quinn and Deriso. Asymptotic recruitment with
increasing spawner biomass.
"""
struct BevertonHoltSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModel
    α::T
    β::T

    function BevertonHoltSpawnerRecruitModel(α::T, β::T) where T
        new{T}(α, β)
    end
end

function recruit(spawners, srmod::BevertonHoltSpawnerRecruitModel)
    srmod.α * spawners / (1 + srmod.β * spawners)
end

function max_recruits(srmod::BevertonHoltSpawnerRecruitModel)
    srmod.α / srmod.β
end

function max_spawnrecruits(srmod::BevertonHoltSpawnerRecruitModel)
    (Inf, max_recruits(srmod))
end

"""
    RickerSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModels

Parameters

- α: productivity parameter
- β: density-dependence parameter

Parameterized as in equation 3.8 of Quinn and Deriso. Dome-shaped recruitment
curve.
"""
struct RickerSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModels
    α::T
    β::T

    function RickerSpawnerRecruitModel(α::T, β::T) where T
        new{T}(α, β)
    end
end

function recruit(spawners, srmod::RickerSpawnerRecruitModel)
    srmod.α * spawners * exp(-srmod.β * spawners)
end

function max_recruits(srmod::RickerSpawnerRecruitModel)
    srmod.α / (srmod.β * ℯ)
end

function max_spawnrecruits(srmod::RickerSpawnerRecruitModel)
    (1 / srmod.β, max_recruits(srmod))
end

"""
    LudwigWaltersSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModels

Parameters

- α: productivity parameter
- β: density-dependence parameter
- γ: density-dependence power parameter

Generalization of Ricker spawner recruit model. Parameterized as in equation 3.9
of Quinn and Deriso. Dome-shaped recruitment curve.
"""
struct LudwigWaltersSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModel
    α::T
    β::T
    γ::T

    function LudwigWaltersSpawnerRecruitModel(α::T, β::T, γ::T) where T
        new{T}(α, β, γ)
    end
end

function recruit(spawners, srmod::LudwigWaltersSpawnerRecruitModel)
    srmod.α * spawners * exp(-srmod.β * spawners ^ srmod.γ)
end

function max_recruits(srmod::LudwigWaltersSpawnerRecruitModel)
    srmod.α * (srmod.β * srmod.γ) ^ -(1 / srmod.γ)
end

function max_spawnrecruits(srmod::LudwigWaltersSpawnerRecruitModel)
    (exp(-1 / srmod.γ), max_recruits(srmod))
end

"""
    CushingSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModels

Parameters

- α: productivity parameter
- γ: index of density dependence

Generalization of Ricker spawner recruit model. Parameterized as in equation 3.12
of Quinn and Deriso. Dome-shaped recruitment curve.
"""
struct CushingSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModel
    α::T
    γ::T

    function CushingSpawnerRecruitModel(α::T, γ::T) where T
        new{T}(α, γ)
    end
end

function recruit(spawners, srmod::CushingSpawnerRecruitModel)
    srmod.α * spawners ^ srmod.γ
end

max_recruits(srmod::CushingSpawnerRecruitModel) = Inf

max_spawnrecruits(srmod::CushingSpawnerRecruitModel) = (Inf, Inf)

"""
    DerisoSchnuteSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModels

Parameters

- α: productivity parameter
- β: optimality parameter
- γ: recruitment limitation or skewness parameter

Generalization of Beverton-Holt and Ricker spawner recruit models. Parameterized
as in equation 3.20 of Quinn and Deriso. May be asymptotic or dome-shaped
recruitment curve, depending on the value of γ.
"""
struct DerisoSchnuteSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModel
    α::T
    β::T
    γ::T

    function DerisoSchnuteSpawnerRecruitModel(α::T, β::T, γ::T) where T
        new{T}(α, β, γ)
    end
end

function recruit(spawners, srmod::DerisoSchnuteSpawnerRecruitModel)
    srmod.α * spawners * (1 - srmod.β * srmod.γ * spawners) ^ (1 / srmod.γ)
end

function max_recruits(srmod::DerisoSchnuteSpawnerRecruitModel)
    (srmod.α / srmod.β) * (1 + srmod.γ) ^ -((1 + srmod.γ) / srmod.γ)
end

function max_spawnrecruits(srmod::DerisoSchnuteSpawnerRecruitModels)
    (1 / (srmod.β * (1 + srmod.γ)), max_recruits(srmod))
end

"""
    ShepherdSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModels

Parameters

- α: productivity parameter
- γ: index of density dependence

Generalization of Beverton-Holt, Ricker, and Cushing spawner recruit models.
Parameterized as in equation 3.21 of Quinn and Deriso. Unlimited recruitment for
γ < 1 as in the Cushing model, asymptotic recruitment (Beverton-Holt)
recruitment at γ = 1, and dome-shaped recruitment at γ > 1.
"""
struct ShepherdSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModel
    α::T
    β::T
    γ::T

    function ShepherdSpawnerRecruitModel(α::T, β::T, γ::T) where T
        new{T}(α, β, γ)
    end
end

function recruit(spawners, srmod::ShepherdSpawnerRecruitModel)
    srmod.α * spawners / (1 + β * spawners ^ srmod.γ)
end

function max_recruits(srmod::ShepherdSpawnerRecruitModel)
    if srmod.γ > 1
        γm1 = srmod.γ - 1
        max_rec = srmod.α / (srmod.β * srmod.γ) * (srmod.β * γm1) ^ (γm1 / srmod.γ)
    elseif srmod.γ == 1
        max_rec = max_recruits(BevertonHoltSpawnerRecruitModel(srmod.α, srmod.β))
    else
        max_rec = Inf
    end
    max_rec
end

function max_spawnrecruits(srmod::ShepherdSpawnerRecruitModel)
    if srmod.γ > 1
        spawnmax = (srmod.β * (srmod.γ - 1)) ^ (-1 / srmod.γ)
    elseif srmod.γ ≤ 1
        spawnmax = Inf
    end
    (spawnmax, max_recruits(srmod))
end

"""
    GammaSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModel

Parameters:

- α: productivity parameter
- β:
- γ:

Takes the form of an unnomalized gamma function. Contains the Ricker (γ = 1) and
Cushing (β = 0) as special cases. Approaches Beverton-Holt as γ → 0
"""
struct GammaSpawnerRecruitModel{T} <: AbstractSpawnerRecruitModel
    α::T
    β::T
    γ::T

    function GammaSpawnerRecruitModel(α::T, β::T, γ::T) where T
        γ > 0 || DomainError("γ must be positive for a valid spawner-recruit model")
        new{T}(α, β, γ)
    end
end

function recruit(spawners, srmod::GammaSpawnerRecruitModel)
    srmod.α * spawners ^ srmod.γ * exp(-srmod.β * spawners)
end

function max_recruits(srmod::GammaSpawnerRecruitModel)
    srmod.α * (srmod.γ / srmod.β) ^ srmod.γ * exp(-srmod.γ)
end

function max_spawnrecruits(srmod::GammaSpawnerRecruitModel)
    (srmod.γ / srmod.β, max_recruits(srmod))
end

end
