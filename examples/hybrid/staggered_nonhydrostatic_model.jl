using LinearAlgebra: ×, norm, norm_sqr, dot

using ClimaCore: Operators, Fields

using Thermodynamics
using CLIMAParameters: AbstractEarthParameterSet, Planet, Atmos

const TD = Thermodynamics

include("schur_complement_W.jl")
include("hyperdiffusion.jl")

f_plane_coriolis_frequency(::AbstractEarthParameterSet) = 0

# Note: FT must be defined before `include("staggered_nonhydrostatic_model.jl")`

# Functions on which the model depends:
# Planet.R_d(params)         # dry specific gas constant
# Planet.kappa_d(params)     # dry adiabatic exponent
# Planet.T_triple(params)    # triple point temperature of water
# Planet.MSLP(params)        # reference pressure
# Planet.grav(params)        # gravitational acceleration
# Planet.Omega(params)       # rotation rate (only used if space is spherical)
# Planet.cv_d(params)        # dry isochoric specific heat capacity
# The value of cv_d is implied by the values of R_d and kappa_d

# The model also depends on f_plane_coriolis_frequency(params)
# This is a constant Coriolis frequency that is only used if space is flat

# To add additional terms to the explicit part of the tendency, define new
# methods for `additional_cache` and `additional_tendency!`.

const divₕ = Operators.Divergence()
const wdivₕ = Operators.WeakDivergence()
const gradₕ = Operators.Gradient()
const wgradₕ = Operators.WeakGradient()
const curlₕ = Operators.Curl()
const wcurlₕ = Operators.WeakCurl()

const ᶜinterp = Operators.InterpolateF2C()
const ᶠinterp = Operators.InterpolateC2F(
    bottom = Operators.Extrapolate(),
    top = Operators.Extrapolate(),
)
const ᶜdivᵥ = Operators.DivergenceF2C(
    top = Operators.SetValue(Geometry.Contravariant3Vector(FT(0))),
    bottom = Operators.SetValue(Geometry.Contravariant3Vector(FT(0))),
)
const ᶠgradᵥ = Operators.GradientC2F(
    bottom = Operators.SetGradient(Geometry.Covariant3Vector(FT(0))),
    top = Operators.SetGradient(Geometry.Covariant3Vector(FT(0))),
)
const ᶠcurlᵥ = Operators.CurlC2F(
    bottom = Operators.SetCurl(Geometry.Contravariant12Vector(FT(0), FT(0))),
    top = Operators.SetCurl(Geometry.Contravariant12Vector(FT(0), FT(0))),
)
const ᶜFC = Operators.FluxCorrectionC2C(
    bottom = Operators.Extrapolate(),
    top = Operators.Extrapolate(),
)

const ᶜinterp_stencil = Operators.Operator2Stencil(ᶜinterp)
const ᶠinterp_stencil = Operators.Operator2Stencil(ᶠinterp)
const ᶜdivᵥ_stencil = Operators.Operator2Stencil(ᶜdivᵥ)
const ᶠgradᵥ_stencil = Operators.Operator2Stencil(ᶠgradᵥ)

const C123 = Geometry.Covariant123Vector

partition(Yc) =
    TD.PhasePartition(Yc.ρq_tot / Yc.ρ, Yc.ρq_liq / Yc.ρ, Yc.ρq_ice / Yc.ρ)
function thermo_state_ρθ(ρθ, Yc, params) # Note: θ is liquid-ice potential temp
    if (
        :ρq_liq in propertynames(Yc) &&
        :ρq_ice in propertynames(Yc) &&
        :ρq_tot in propertynames(Yc)
    )
        return TD.PhaseNonEquil_ρθq(params, Yc.ρ, ρθ / Yc.ρ, partition(Yc))
    elseif :ρq_tot in propertynames(Yc)
        return TD.PhaseEquil_ρθq(params, Yc.ρ, ρθ / Yc.ρ, Yc.ρq_tot / Yc.ρ)
    else
        return TD.PhaseDry_ρθ(params, Yc.ρ, ρθ / Yc.ρ)
    end
end
function thermo_state_ρe_int(ρe_int, Yc, params)
    if (
        :ρq_liq in propertynames(Yc) &&
        :ρq_ice in propertynames(Yc) &&
        :ρq_tot in propertynames(Yc)
    )
        return TD.PhaseNonEquil(params, ρe_int / Yc.ρ, Yc.ρ, partition(Yc))
    elseif :ρq_tot in propertynames(Yc)
        return TD.PhaseEquil_ρeq(params, Yc.ρ, ρe_int / Yc.ρ, Yc.ρq_tot / Yc.ρ)
    else
        return TD.PhaseDry(params, ρe_int / Yc.ρ, Yc.ρ)
    end
end
thermo_state_ρe(ρe, Yc, K, Φ, params) =
    thermo_state_ρe_int(ρe - Yc.ρ * (K + Φ), Yc, params)

get_cache(Y, params, dt) =
    merge(default_cache(Y, params), additional_cache(Y, params, dt))

function default_cache(Y, params)
    ᶜcoord = Fields.local_geometry_field(Y.c).coordinates
    if eltype(ᶜcoord) <: Geometry.LatLongZPoint
        Ω = FT(Planet.Omega(params))
        ᶜf = @. 2 * Ω * sind(ᶜcoord.lat)
    else
        f = FT(f_plane_coriolis_frequency(params))
        ᶜf = map(_ -> f, ᶜcoord)
    end
    ᶜf = @. Geometry.Contravariant3Vector(Geometry.WVector(ᶜf))
    if (
        :ρq_liq in propertynames(Y.c) &&
        :ρq_ice in propertynames(Y.c) &&
        :ρq_tot in propertynames(Y.c)
    )
        ts_type = TD.PhaseNonEquil{FT}
    elseif :ρq_tot in propertynames(Y.c)
        ts_type = TD.PhaseEquil{FT}
    else
        ts_type = TD.PhaseDry{FT}
    end
    return (;
        ᶜuvw = similar(Y.c, Geometry.Covariant123Vector{FT}),
        ᶜK = similar(Y.c, FT),
        ᶜΦ = FT(Planet.grav(params)) .* ᶜcoord.z,
        ᶜts = similar(Y.c, ts_type),
        ᶜp = similar(Y.c, FT),
        ᶜω³ = similar(Y.c, Geometry.Contravariant3Vector{FT}),
        ᶠω¹² = similar(Y.f, Geometry.Contravariant12Vector{FT}),
        ᶠu¹² = similar(Y.f, Geometry.Contravariant12Vector{FT}),
        ᶠu³ = similar(Y.f, Geometry.Contravariant3Vector{FT}),
        ᶜf,
        ∂ᶜK∂ᶠw_data = similar(
            Y.c,
            Operators.StencilCoefs{-half, half, NTuple{2, FT}},
        ),
        params,
        ghost_buffer = (
            c = Spaces.create_ghost_buffer(Y.c),
            f = Spaces.create_ghost_buffer(Y.f),
            χ = Spaces.create_ghost_buffer(Y.c.ρ), # for hyperdiffusion
            χw = Spaces.create_ghost_buffer(Y.f.w.components.data.:1), # for hyperdiffusion
            χuₕ = Spaces.create_ghost_buffer(Y.c.uₕ), # for hyperdiffusion
        ),
    )
end

is_tracer(name) = !(name in (:ρ, :ρθ, :ρe, :ρe_int, :uₕ, :w))

function implicit_tendency!(Yₜ, Y, p, t)
    ᶜρ = Y.c.ρ
    ᶜuₕ = Y.c.uₕ
    ᶠw = Y.f.w
    (; ᶜK, ᶜΦ, ᶜts, ᶜp, params) = p

    ᶠupwind_product = Operators.UpwindBiasedProductC2F()

    # Used for automatically computing the Jacobian ∂Yₜ/∂Y. Currently requires
    # allocation because the cache is stored separately from Y, which means that
    # similar(Y, <:Dual) doesn't allocate an appropriate cache for computing Yₜ.
    if eltype(Y) <: Dual
        ᶜK = similar(ᶜρ)
        ᶜp = similar(ᶜρ)
    end

    @. ᶜK = norm_sqr(C123(ᶜuₕ) + C123(ᶜinterp(ᶠw))) / 2

    @. Yₜ.c.ρ = -(ᶜdivᵥ(ᶠinterp(ᶜρ) * ᶠw))

    if :ρθ in propertynames(Y.c)
        @. ᶜts = thermo_state_ρθ(Y.c.ρθ, Y.c, params)
        @. ᶜp = TD.air_pressure(params, ᶜts)
        @. Yₜ.c.ρθ =
            -(ᶜdivᵥ(ᶠinterp(Y.c.ρ) * ᶠupwind_product(ᶠw, Y.c.ρθ / Y.c.ρ)))
    elseif :ρe in propertynames(Y.c)
        @. ᶜts = thermo_state_ρe(Y.c.ρe, Y.c, ᶜK, ᶜΦ, params)
        @. ᶜp = TD.air_pressure(params, ᶜts)
        @. Yₜ.c.ρe =
            -(ᶜdivᵥ(
                ᶠinterp(Y.c.ρ) * ᶠupwind_product(ᶠw, (Y.c.ρe + ᶜp) / Y.c.ρ),
            ))
    elseif :ρe_int in propertynames(Y.c)
        @. ᶜts = thermo_state_ρe_int(Y.c.ρe_int, Y.c, params)
        @. ᶜp = TD.air_pressure(params, ᶜts)
        @. Yₜ.c.ρe_int =
            -(
                (ᶜdivᵥ(
                    ᶠinterp(Y.c.ρ) *
                    ᶠupwind_product(ᶠw, (Y.c.ρe_int + ᶜp) / Y.c.ρ),
                )) -
                ᶜinterp(dot(ᶠgradᵥ(ᶜp), Geometry.Contravariant3Vector(ᶠw)))
            )
        # or, equivalently,
        # @. Yₜ.c.ρe_int = -(ᶜdivᵥ(ᶠinterp(Y.c.ρe_int) * ᶠw) + ᶜp * ᶜdivᵥ(ᶠw))
    end

    # TODO: Add flux correction to the Jacobian
    # @. Yₜ.c.ρ += ᶜFC(ᶠw, ᶜρ)
    # if :ρθ in propertynames(Y.c)
    #     @. Yₜ.c.ρθ += ᶜFC(ᶠw, ᶜρθ)
    # elseif :ρe in propertynames(Y.c)
    #     @. Yₜ.c.ρe += ᶜFC(ᶠw, ᶜρe)
    # elseif :ρe_int in propertynames(Y.c)
    #     @. Yₜ.c.ρe_int += ᶜFC(ᶠw, ᶜρe_int)
    # end

    Yₜ.c.uₕ .= Ref(zero(eltype(Yₜ.c.uₕ)))

    @. Yₜ.f.w = -(ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ) + ᶠgradᵥ(ᶜK + ᶜΦ))

    # TODO: Add vertical advection of tracers to the Jacobian
    for tracer_name in filter(is_tracer, propertynames(Y.c))
        ᶜtracerₜ = getproperty(Yₜ.c, tracer_name)
        ᶜtracerₜ .= zero(eltype(ᶜtracerₜ))
    end

    return Yₜ
end

function remaining_tendency!(Yₜ, Y, p, t)
    Yₜ .= zero(eltype(Yₜ))
    default_remaining_tendency!(Yₜ, Y, p, t)
    additional_tendency!(Yₜ, Y, p, t)
    Spaces.weighted_dss!(Yₜ.c, p.ghost_buffer.c)
    Spaces.weighted_dss!(Yₜ.f, p.ghost_buffer.f)
    return Yₜ
end

function default_remaining_tendency!(Yₜ, Y, p, t)
    ᶜρ = Y.c.ρ
    ᶜuₕ = Y.c.uₕ
    ᶠw = Y.f.w
    (; ᶜuvw, ᶜK, ᶜΦ, ᶜts, ᶜp, ᶜω³, ᶠω¹², ᶠu¹², ᶠu³, ᶜf, params) = p
    point_type = eltype(Fields.local_geometry_field(axes(Y.c)).coordinates)

    @. ᶜuvw = C123(ᶜuₕ) + C123(ᶜinterp(ᶠw))
    @. ᶜK = norm_sqr(ᶜuvw) / 2

    # Mass conservation

    @. Yₜ.c.ρ -= divₕ(ᶜρ * ᶜuvw)
    @. Yₜ.c.ρ -= ᶜdivᵥ(ᶠinterp(ᶜρ * ᶜuₕ))

    # Energy conservation
    ᶠupwind_product = Operators.UpwindBiasedProductC2F()

    if :ρθ in propertynames(Y.c)
        @. ᶜts = thermo_state_ρθ(Y.c.ρθ, Y.c, params)
        @. ᶜp = TD.air_pressure(params, ᶜts)
        @. Yₜ.c.ρθ -= divₕ(Y.c.ρθ * ᶜuvw)
        @. Yₜ.c.ρθ -= ᶜdivᵥ(ᶠinterp(Y.c.ρθ * ᶜuₕ))
    elseif :ρe in propertynames(Y.c)
        @. ᶜts = thermo_state_ρe(Y.c.ρe, Y.c, ᶜK, ᶜΦ, params)
        @. ᶜp = TD.air_pressure(params, ᶜts)
        @. Yₜ.c.ρe -= divₕ((Y.c.ρe + ᶜp) * ᶜuvw)
        @. Yₜ.c.ρe -= ᶜdivᵥ(ᶠinterp((Y.c.ρe + ᶜp) * ᶜuₕ))
    elseif :ρe_int in propertynames(Y.c)
        @. ᶜts = thermo_state_ρe_int(Y.c.ρe_int, Y.c, params)
        @. ᶜp = TD.air_pressure(params, ᶜts)
        if point_type <: Geometry.Abstract3DPoint
            @. Yₜ.c.ρe_int -=
                divₕ((Y.c.ρe_int + ᶜp) * ᶜuvw) -
                dot(gradₕ(ᶜp), Geometry.Contravariant12Vector(ᶜuₕ))
        else
            @. Yₜ.c.ρe_int -=
                divₕ((Y.c.ρe_int + ᶜp) * ᶜuvw) -
                dot(gradₕ(ᶜp), Geometry.Contravariant1Vector(ᶜuₕ))
        end
        @. Yₜ.c.ρe_int -= ᶜdivᵥ(ᶠinterp((Y.c.ρe_int + ᶜp) * ᶜuₕ))
        # or, equivalently,
        # @. Yₜ.c.ρe_int -= divₕ(Y.c.ρe_int * ᶜuvw) + ᶜp * divₕ(ᶜuvw)
        # @. Yₜ.c.ρe_int -=
        #     ᶜdivᵥ(ᶠinterp(Y.c.ρe_int * ᶜuₕ)) + ᶜp * ᶜdivᵥ(ᶠinterp(ᶜuₕ))
    end

    # Momentum conservation

    if point_type <: Geometry.Abstract3DPoint
        @. ᶜω³ = curlₕ(ᶜuₕ)
        @. ᶠω¹² = curlₕ(ᶠw)
    elseif point_type <: Geometry.Abstract2DPoint
        ᶜω³ .= Ref(zero(eltype(ᶜω³)))
        @. ᶠω¹² = Geometry.Contravariant12Vector(curlₕ(ᶠw))
    end
    @. ᶠω¹² += ᶠcurlᵥ(ᶜuₕ)

    # TODO: Modify to account for topography
    @. ᶠu¹² = Geometry.Contravariant12Vector(ᶠinterp(ᶜuₕ))
    @. ᶠu³ = Geometry.Contravariant3Vector(ᶠw)

    @. Yₜ.c.uₕ -=
        ᶜinterp(ᶠω¹² × ᶠu³) + (ᶜf + ᶜω³) × Geometry.Contravariant12Vector(ᶜuₕ)
    if point_type <: Geometry.Abstract3DPoint
        @. Yₜ.c.uₕ -= gradₕ(ᶜp) / ᶜρ + gradₕ(ᶜK + ᶜΦ)
    elseif point_type <: Geometry.Abstract2DPoint
        @. Yₜ.c.uₕ -=
            Geometry.Covariant12Vector(gradₕ(ᶜp) / ᶜρ + gradₕ(ᶜK + ᶜΦ))
    end

    @. Yₜ.f.w -= ᶠω¹² × ᶠu¹²

    # Tracer conservation

    for tracer_name in filter(is_tracer, propertynames(Y.c))
        ᶜtracer = getproperty(Y.c, tracer_name)
        ᶜtracerₜ = getproperty(Yₜ.c, tracer_name)
        @. ᶜtracerₜ -= divₕ(ᶜtracer * ᶜuvw)
        @. ᶜtracerₜ -=
            ᶜdivᵥ(ᶠinterp(Y.c.ρ) * ᶠupwind_product(ᶠw, ᶜtracer / Y.c.ρ)) # TODO: put in implicit tend.
        @. ᶜtracerₜ -= ᶜdivᵥ(ᶠinterp(ᶜtracer * ᶜuₕ))
    end
end

additional_tendency!(Yₜ, Y, p, t) = nothing

# Allow one() to be called on vectors.
Base.one(::T) where {T <: Geometry.AxisTensor} = one(T)
Base.one(::Type{T}) where {T′, A, S, T <: Geometry.AxisTensor{T′, 1, A, S}} =
    T(axes(T), S(one(T′)))

function Wfact!(W, Y, p, dtγ, t)
    (; flags, dtγ_ref, ∂ᶜρₜ∂ᶠ𝕄, ∂ᶜ𝔼ₜ∂ᶠ𝕄, ∂ᶠ𝕄ₜ∂ᶜ𝔼, ∂ᶠ𝕄ₜ∂ᶜρ, ∂ᶠ𝕄ₜ∂ᶠ𝕄) = W
    ᶜρ = Y.c.ρ
    ᶜuₕ = Y.c.uₕ
    ᶠw = Y.f.w
    (; ᶜK, ᶜΦ, ᶜts, ᶜp, ∂ᶜK∂ᶠw_data, params) = p

    R_d = FT(Planet.R_d(params))
    κ_d = FT(Planet.kappa_d(params))
    cv_d = FT(Planet.cv_d(params))
    T_tri = FT(Planet.T_triple(params))
    MSLP = FT(Planet.MSLP(params))

    dtγ_ref[] = dtγ

    # If we let ᶠw_data = ᶠw.components.data.:1 and ᶠw_unit = one.(ᶠw), then
    # ᶠw == ᶠw_data .* ᶠw_unit. The Jacobian blocks involve ᶠw_data, not ᶠw.
    ᶠw_data = ᶠw.components.data.:1

    # ᶜinterp(ᶠw) =
    #     ᶜinterp(ᶠw)_data * ᶜinterp(ᶠw)_unit =
    #     ᶜinterp(ᶠw_data) * ᶜinterp(ᶠw)_unit
    # norm_sqr(ᶜinterp(ᶠw)) =
    #     norm_sqr(ᶜinterp(ᶠw_data) * ᶜinterp(ᶠw)_unit) =
    #     ᶜinterp(ᶠw_data)^2 * norm_sqr(ᶜinterp(ᶠw)_unit)
    # ᶜK =
    #     norm_sqr(C123(ᶜuₕ) + C123(ᶜinterp(ᶠw))) / 2 =
    #     norm_sqr(ᶜuₕ) / 2 + norm_sqr(ᶜinterp(ᶠw)) / 2 =
    #     norm_sqr(ᶜuₕ) / 2 + ᶜinterp(ᶠw_data)^2 * norm_sqr(ᶜinterp(ᶠw)_unit) / 2
    # ∂(ᶜK)/∂(ᶠw_data) =
    #     ∂(ᶜK)/∂(ᶜinterp(ᶠw_data)) * ∂(ᶜinterp(ᶠw_data))/∂(ᶠw_data) =
    #     ᶜinterp(ᶠw_data) * norm_sqr(ᶜinterp(ᶠw)_unit) * ᶜinterp_stencil(1)
    @. ∂ᶜK∂ᶠw_data =
        ᶜinterp(ᶠw_data) *
        norm_sqr(one(ᶜinterp(ᶠw))) *
        ᶜinterp_stencil(one(ᶠw_data))

    # ᶜρₜ = -ᶜdivᵥ(ᶠinterp(ᶜρ) * ᶠw)
    # ∂(ᶜρₜ)/∂(ᶠw_data) = -ᶜdivᵥ_stencil(ᶠinterp(ᶜρ) * ᶠw_unit)
    @. ∂ᶜρₜ∂ᶠ𝕄 = -(ᶜdivᵥ_stencil(ᶠinterp(ᶜρ) * one(ᶠw)))

    if :ρθ in propertynames(Y.c)
        ᶜρθ = Y.c.ρθ
        @. ᶜts = thermo_state_ρθ(Y.c.ρθ, Y.c, params)
        @. ᶜp = TD.air_pressure(params, ᶜts)

        if flags.∂ᶜ𝔼ₜ∂ᶠ𝕄_mode != :exact
            error("∂ᶜ𝔼ₜ∂ᶠ𝕄_mode must be :exact when using ρθ")
        end

        # ᶜρθₜ = -ᶜdivᵥ(ᶠinterp(ᶜρθ) * ᶠw)
        # ∂(ᶜρθₜ)/∂(ᶠw_data) = -ᶜdivᵥ_stencil(ᶠinterp(ᶜρθ) * ᶠw_unit)
        @. ∂ᶜ𝔼ₜ∂ᶠ𝕄 = -(ᶜdivᵥ_stencil(ᶠinterp(ᶜρθ) * one(ᶠw)))
    elseif :ρe in propertynames(Y.c)
        ᶜρe = Y.c.ρe
        @. ᶜK = norm_sqr(C123(ᶜuₕ) + C123(ᶜinterp(ᶠw))) / 2
        @. ᶜts = thermo_state_ρe(Y.c.ρe, Y.c, ᶜK, ᶜΦ, params)
        @. ᶜp = TD.air_pressure(params, ᶜts)

        if flags.∂ᶜ𝔼ₜ∂ᶠ𝕄_mode == :exact
            # ᶜρeₜ = -ᶜdivᵥ(ᶠinterp(ᶜρe + ᶜp) * ᶠw)
            # ∂(ᶜρeₜ)/∂(ᶠw_data) =
            #     -ᶜdivᵥ_stencil(ᶠinterp(ᶜρe + ᶜp) * ᶠw_unit) -
            #     ᶜdivᵥ_stencil(ᶠw) * ∂(ᶠinterp(ᶜρe + ᶜp))/∂(ᶠw_data)
            # ∂(ᶠinterp(ᶜρe + ᶜp))/∂(ᶠw_data) =
            #     ∂(ᶠinterp(ᶜρe + ᶜp))/∂(ᶜp) * ∂(ᶜp)/∂(ᶠw_data)
            # ∂(ᶠinterp(ᶜρe + ᶜp))/∂(ᶜp) = ᶠinterp_stencil(1)
            # ∂(ᶜp)/∂(ᶠw_data) = ∂(ᶜp)/∂(ᶜK) * ∂(ᶜK)/∂(ᶠw_data)
            # ∂(ᶜp)/∂(ᶜK) = -ᶜρ * R_d / cv_d
            @. ∂ᶜ𝔼ₜ∂ᶠ𝕄 =
                -(ᶜdivᵥ_stencil(ᶠinterp(ᶜρe + ᶜp) * one(ᶠw))) - compose(
                    ᶜdivᵥ_stencil(ᶠw),
                    compose(
                        ᶠinterp_stencil(one(ᶜp)),
                        -(ᶜρ * R_d / cv_d) * ∂ᶜK∂ᶠw_data,
                    ),
                )
        elseif flags.∂ᶜ𝔼ₜ∂ᶠ𝕄_mode == :no_∂ᶜp∂ᶜK
            # same as above, but we approximate ∂(ᶜp)/∂(ᶜK) = 0, so that ∂ᶜ𝔼ₜ∂ᶠ𝕄
            # has 3 diagonals instead of 5
            @. ∂ᶜ𝔼ₜ∂ᶠ𝕄 = -(ᶜdivᵥ_stencil(ᶠinterp(ᶜρe + ᶜp) * one(ᶠw)))
        else
            error("∂ᶜ𝔼ₜ∂ᶠ𝕄_mode must be :exact or :no_∂ᶜp∂ᶜK when using ρe")
        end
    elseif :ρe_int in propertynames(Y.c)
        ᶜρe_int = Y.c.ρe_int
        @. ᶜts = thermo_state_ρe_int(Y.c.ρe_int, Y.c, params)
        @. ᶜp = TD.air_pressure(params, ᶜts)

        if flags.∂ᶜ𝔼ₜ∂ᶠ𝕄_mode != :exact
            error("∂ᶜ𝔼ₜ∂ᶠ𝕄_mode must be :exact when using ρe_int")
        end

        # ᶜρe_intₜ =
        #     -ᶜdivᵥ(ᶠinterp(ᶜρe_int + ᶜp) * ᶠw) +
        #     ᶜinterp(dot(ᶠgradᵥ(ᶜp), Geometry.Contravariant3Vector(ᶠw)))
        # ∂(ᶜρe_intₜ)/∂(ᶠw_data) =
        #     -ᶜdivᵥ_stencil(ᶠinterp(ᶜρe_int + ᶜp) * ᶠw_unit) + ᶜinterp_stencil(
        #         dot(ᶠgradᵥ(ᶜp), Geometry.Contravariant3Vector(ᶠw_unit)),
        #     )
        @. ∂ᶜ𝔼ₜ∂ᶠ𝕄 =
            -(ᶜdivᵥ_stencil(ᶠinterp(ᶜρe_int + ᶜp) * one(ᶠw))) +
            ᶜinterp_stencil(dot(
                ᶠgradᵥ(ᶜp),
                Geometry.Contravariant3Vector(one(ᶠw)),
            ),)
    end

    # To convert ∂(ᶠwₜ)/∂(ᶜ𝔼) to ∂(ᶠw_data)ₜ/∂(ᶜ𝔼) and ∂(ᶠwₜ)/∂(ᶠw_data) to
    # ∂(ᶠw_data)ₜ/∂(ᶠw_data), we must extract the third component of each
    # vector-valued stencil coefficient.
    to_scalar_coefs(vector_coefs) =
        map(vector_coef -> vector_coef.u₃, vector_coefs)

    # TODO: If we end up using :gradΦ_shenanigans, optimize it to
    # `cached_stencil / ᶠinterp(ᶜρ)`.
    if flags.∂ᶠ𝕄ₜ∂ᶜρ_mode != :exact && flags.∂ᶠ𝕄ₜ∂ᶜρ_mode != :gradΦ_shenanigans
        error("∂ᶠ𝕄ₜ∂ᶜρ_mode must be :exact or :gradΦ_shenanigans")
    end
    if :ρθ in propertynames(Y.c)
        # ᶠwₜ = -ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ) - ᶠgradᵥ(ᶜK + ᶜΦ)
        # ∂(ᶠwₜ)/∂(ᶜρθ) = ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) * ∂(ᶠgradᵥ(ᶜp))/∂(ᶜρθ)
        # ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) = -1 / ᶠinterp(ᶜρ)
        # ∂(ᶠgradᵥ(ᶜp))/∂(ᶜρθ) =
        #     ᶠgradᵥ_stencil(
        #         R_d / (1 - κ_d) * (ᶜρθ * R_d / MSLP)^(κ_d / (1 - κ_d))
        #     )
        @. ∂ᶠ𝕄ₜ∂ᶜ𝔼 = to_scalar_coefs(
            -1 / ᶠinterp(ᶜρ) * ᶠgradᵥ_stencil(
                R_d / (1 - κ_d) * (ᶜρθ * R_d / MSLP)^(κ_d / (1 - κ_d)),
            ),
        )

        if flags.∂ᶠ𝕄ₜ∂ᶜρ_mode == :exact
            # ᶠwₜ = -ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ) - ᶠgradᵥ(ᶜK + ᶜΦ)
            # ∂(ᶠwₜ)/∂(ᶜρ) = ∂(ᶠwₜ)/∂(ᶠinterp(ᶜρ)) * ∂(ᶠinterp(ᶜρ))/∂(ᶜρ)
            # ∂(ᶠwₜ)/∂(ᶠinterp(ᶜρ)) = ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ)^2
            # ∂(ᶠinterp(ᶜρ))/∂(ᶜρ) = ᶠinterp_stencil(1)
            @. ∂ᶠ𝕄ₜ∂ᶜρ = to_scalar_coefs(
                ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ)^2 * ᶠinterp_stencil(one(ᶜρ)),
            )
        elseif flags.∂ᶠ𝕄ₜ∂ᶜρ_mode == :gradΦ_shenanigans
            # ᶠwₜ = (
            #     -ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ′) -
            #     ᶠgradᵥ(ᶜΦ) / ᶠinterp(ᶜρ′) * ᶠinterp(ᶜρ)
            # ), where ᶜρ′ = ᶜρ but we approximate ∂(ᶜρ′)/∂(ᶜρ) = 0
            @. ∂ᶠ𝕄ₜ∂ᶜρ = to_scalar_coefs(
                -(ᶠgradᵥ(ᶜΦ)) / ᶠinterp(ᶜρ) * ᶠinterp_stencil(one(ᶜρ)),
            )
        end
    elseif :ρe in propertynames(Y.c)
        # ᶠwₜ = -ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ) - ᶠgradᵥ(ᶜK + ᶜΦ)
        # ∂(ᶠwₜ)/∂(ᶜρe) = ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) * ∂(ᶠgradᵥ(ᶜp))/∂(ᶜρe)
        # ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) = -1 / ᶠinterp(ᶜρ)
        # ∂(ᶠgradᵥ(ᶜp))/∂(ᶜρe) = ᶠgradᵥ_stencil(R_d / cv_d)
        @. ∂ᶠ𝕄ₜ∂ᶜ𝔼 = to_scalar_coefs(
            -1 / ᶠinterp(ᶜρ) * ᶠgradᵥ_stencil(R_d / cv_d * one(ᶜρe)),
        )

        if flags.∂ᶠ𝕄ₜ∂ᶜρ_mode == :exact
            # ᶠwₜ = -ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ) - ᶠgradᵥ(ᶜK + ᶜΦ)
            # ∂(ᶠwₜ)/∂(ᶜρ) =
            #     ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) * ∂(ᶠgradᵥ(ᶜp))/∂(ᶜρ) +
            #     ∂(ᶠwₜ)/∂(ᶠinterp(ᶜρ)) * ∂(ᶠinterp(ᶜρ))/∂(ᶜρ)
            # ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) = -1 / ᶠinterp(ᶜρ)
            # ∂(ᶠgradᵥ(ᶜp))/∂(ᶜρ) =
            #     ᶠgradᵥ_stencil(R_d * (-(ᶜK + ᶜΦ) / cv_d + T_tri))
            # ∂(ᶠwₜ)/∂(ᶠinterp(ᶜρ)) = ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ)^2
            # ∂(ᶠinterp(ᶜρ))/∂(ᶜρ) = ᶠinterp_stencil(1)
            @. ∂ᶠ𝕄ₜ∂ᶜρ = to_scalar_coefs(
                -1 / ᶠinterp(ᶜρ) *
                ᶠgradᵥ_stencil(R_d * (-(ᶜK + ᶜΦ) / cv_d + T_tri)) +
                ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ)^2 * ᶠinterp_stencil(one(ᶜρ)),
            )
        elseif flags.∂ᶠ𝕄ₜ∂ᶜρ_mode == :gradΦ_shenanigans
            # ᶠwₜ = (
            #     -ᶠgradᵥ(ᶜp′) / ᶠinterp(ᶜρ′) -
            #     ᶠgradᵥ(ᶜΦ) / ᶠinterp(ᶜρ′) * ᶠinterp(ᶜρ)
            # ), where ᶜρ′ = ᶜρ but we approximate ∂ᶜρ′/∂ᶜρ = 0, and where
            # ᶜp′ = ᶜp but with ᶜK = 0
            @. ∂ᶠ𝕄ₜ∂ᶜρ = to_scalar_coefs(
                -1 / ᶠinterp(ᶜρ) *
                ᶠgradᵥ_stencil(R_d * (-(ᶜΦ) / cv_d + T_tri)) -
                ᶠgradᵥ(ᶜΦ) / ᶠinterp(ᶜρ) * ᶠinterp_stencil(one(ᶜρ)),
            )
        end
    elseif :ρe_int in propertynames(Y.c)
        # ᶠwₜ = -ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ) - ᶠgradᵥ(ᶜK + ᶜΦ)
        # ∂(ᶠwₜ)/∂(ᶜρe_int) = ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) * ∂(ᶠgradᵥ(ᶜp))/∂(ᶜρe_int)
        # ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) = -1 / ᶠinterp(ᶜρ)
        # ∂(ᶠgradᵥ(ᶜp))/∂(ᶜρe_int) = ᶠgradᵥ_stencil(R_d / cv_d)
        @. ∂ᶠ𝕄ₜ∂ᶜ𝔼 = to_scalar_coefs(
            -1 / ᶠinterp(ᶜρ) * ᶠgradᵥ_stencil(R_d / cv_d * one(ᶜρe_int)),
        )

        if flags.∂ᶠ𝕄ₜ∂ᶜρ_mode == :exact
            # ᶠwₜ = -ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ) - ᶠgradᵥ(ᶜK + ᶜΦ)
            # ∂(ᶠwₜ)/∂(ᶜρ) =
            #     ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) * ∂(ᶠgradᵥ(ᶜp))/∂(ᶜρ) +
            #     ∂(ᶠwₜ)/∂(ᶠinterp(ᶜρ)) * ∂(ᶠinterp(ᶜρ))/∂(ᶜρ)
            # ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) = -1 / ᶠinterp(ᶜρ)
            # ∂(ᶠgradᵥ(ᶜp))/∂(ᶜρ) = ᶠgradᵥ_stencil(R_d * T_tri)
            # ∂(ᶠwₜ)/∂(ᶠinterp(ᶜρ)) = ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ)^2
            # ∂(ᶠinterp(ᶜρ))/∂(ᶜρ) = ᶠinterp_stencil(1)
            @. ∂ᶠ𝕄ₜ∂ᶜρ = to_scalar_coefs(
                -1 / ᶠinterp(ᶜρ) * ᶠgradᵥ_stencil(R_d * T_tri * one(ᶜρe_int)) +
                ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ)^2 * ᶠinterp_stencil(one(ᶜρ)),
            )
        elseif flags.∂ᶠ𝕄ₜ∂ᶜρ_mode == :gradΦ_shenanigans
            # ᶠwₜ = (
            #     -ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ′) -
            #     ᶠgradᵥ(ᶜΦ) / ᶠinterp(ᶜρ′) * ᶠinterp(ᶜρ)
            # ), where ᶜp′ = ᶜp but we approximate ∂ᶜρ′/∂ᶜρ = 0
            @. ∂ᶠ𝕄ₜ∂ᶜρ = to_scalar_coefs(
                -1 / ᶠinterp(ᶜρ) * ᶠgradᵥ_stencil(R_d * T_tri * one(ᶜρe_int)) -
                ᶠgradᵥ(ᶜΦ) / ᶠinterp(ᶜρ) * ᶠinterp_stencil(one(ᶜρ)),
            )
        end
    end

    # ᶠwₜ = -ᶠgradᵥ(ᶜp) / ᶠinterp(ᶜρ) - ᶠgradᵥ(ᶜK + ᶜΦ)
    # ∂(ᶠwₜ)/∂(ᶠw_data) =
    #     ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) * ∂(ᶠgradᵥ(ᶜp))/∂(ᶠw_dataₜ) +
    #     ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜK + ᶜΦ)) * ∂(ᶠgradᵥ(ᶜK + ᶜΦ))/∂(ᶠw_dataₜ) =
    #     (
    #         ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) * ∂(ᶠgradᵥ(ᶜp))/∂(ᶜK) +
    #         ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜK + ᶜΦ)) * ∂(ᶠgradᵥ(ᶜK + ᶜΦ))/∂(ᶜK)
    #     ) * ∂(ᶜK)/∂(ᶠw_dataₜ)
    # ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜp)) = -1 / ᶠinterp(ᶜρ)
    # ∂(ᶠgradᵥ(ᶜp))/∂(ᶜK) =
    #     ᶜ𝔼_name == :ρe ? ᶠgradᵥ_stencil(-ᶜρ * R_d / cv_d) : 0
    # ∂(ᶠwₜ)/∂(ᶠgradᵥ(ᶜK + ᶜΦ)) = -1
    # ∂(ᶠgradᵥ(ᶜK + ᶜΦ))/∂(ᶜK) = ᶠgradᵥ_stencil(1)
    # ∂(ᶜK)/∂(ᶠw_data) =
    #     ᶜinterp(ᶠw_data) * norm_sqr(ᶜinterp(ᶠw)_unit) * ᶜinterp_stencil(1)
    if :ρθ in propertynames(Y.c) || :ρe_int in propertynames(Y.c)
        @. ∂ᶠ𝕄ₜ∂ᶠ𝕄 =
            to_scalar_coefs(compose(-1 * ᶠgradᵥ_stencil(one(ᶜK)), ∂ᶜK∂ᶠw_data))
    elseif :ρe in propertynames(Y.c)
        @. ∂ᶠ𝕄ₜ∂ᶠ𝕄 = to_scalar_coefs(compose(
            -1 / ᶠinterp(ᶜρ) * ᶠgradᵥ_stencil(-(ᶜρ * R_d / cv_d)) +
            -1 * ᶠgradᵥ_stencil(one(ᶜK)),
            ∂ᶜK∂ᶠw_data,
        ),)
    end

    if W.test
        # Checking every column takes too long, so just check one.
        i, j, h = 1, 1, 1
        if :ρθ in propertynames(Y.c)
            ᶜ𝔼_name = :ρθ
        elseif :ρe in propertynames(Y.c)
            ᶜ𝔼_name = :ρe
        elseif :ρe_int in propertynames(Y.c)
            ᶜ𝔼_name = :ρe_int
        end
        args = (implicit_tendency!, Y, p, t, i, j, h)
        @assert matrix_column(∂ᶜρₜ∂ᶠ𝕄, axes(Y.f), i, j, h) ==
                exact_column_jacobian_block(args..., (:c, :ρ), (:f, :w))
        @assert matrix_column(∂ᶠ𝕄ₜ∂ᶜ𝔼, axes(Y.c), i, j, h) ≈
                exact_column_jacobian_block(args..., (:f, :w), (:c, ᶜ𝔼_name))
        @assert matrix_column(∂ᶠ𝕄ₜ∂ᶠ𝕄, axes(Y.f), i, j, h) ≈
                exact_column_jacobian_block(args..., (:f, :w), (:f, :w))
        ∂ᶜ𝔼ₜ∂ᶠ𝕄_approx = matrix_column(∂ᶜ𝔼ₜ∂ᶠ𝕄, axes(Y.f), i, j, h)
        ∂ᶜ𝔼ₜ∂ᶠ𝕄_exact =
            exact_column_jacobian_block(args..., (:c, ᶜ𝔼_name), (:f, :w))
        if flags.∂ᶜ𝔼ₜ∂ᶠ𝕄_mode == :exact
            @assert ∂ᶜ𝔼ₜ∂ᶠ𝕄_approx ≈ ∂ᶜ𝔼ₜ∂ᶠ𝕄_exact
        else
            err = norm(∂ᶜ𝔼ₜ∂ᶠ𝕄_approx .- ∂ᶜ𝔼ₜ∂ᶠ𝕄_exact) / norm(∂ᶜ𝔼ₜ∂ᶠ𝕄_exact)
            @assert err < 1e-6
            # Note: the highest value seen so far is ~3e-7 (only applies to ρe)
        end
        ∂ᶠ𝕄ₜ∂ᶜρ_approx = matrix_column(∂ᶠ𝕄ₜ∂ᶜρ, axes(Y.c), i, j, h)
        ∂ᶠ𝕄ₜ∂ᶜρ_exact = exact_column_jacobian_block(args..., (:f, :w), (:c, :ρ))
        if flags.∂ᶠ𝕄ₜ∂ᶜρ_mode == :exact
            @assert ∂ᶠ𝕄ₜ∂ᶜρ_approx ≈ ∂ᶠ𝕄ₜ∂ᶜρ_exact
        else
            err = norm(∂ᶠ𝕄ₜ∂ᶜρ_approx .- ∂ᶠ𝕄ₜ∂ᶜρ_exact) / norm(∂ᶠ𝕄ₜ∂ᶜρ_exact)
            @assert err < 0.03
            # Note: the highest value seen so far for ρe is ~0.01, and the
            # highest value seen so far for ρθ is ~0.02
        end
    end
end
