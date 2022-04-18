using CLIMAParameters: Planet
using ClimaCore: Fields, DataLayouts
using RRTMGP
using NCDatasets

# Utility function for extracting a view of an array from a Field
function field2array(field::Fields.Field)
    if eltype(field) != eltype(parent(field))
        error("field2array only works when the element type of the Field is \
               also the base type of the underlying array")
    end
    return data2array(Fields.field_values(field))
end

const DataLayoutWithV =
    Union{DataLayouts.VF, DataLayouts.VIFH, DataLayouts.VIJFH}
const DataLayoutWithoutV =
    Union{DataLayouts.IF, DataLayouts.IFH, DataLayouts.IJF, DataLayouts.IJFH}
data2array(data::DataLayoutWithV) =
    reshape(parent(data), size(parent(data), 1), :)
data2array(data::DataLayoutWithoutV) = reshape(parent(data), :)
data2array(data::DataLayouts.AbstractData) =
    error("data2array not yet implemented for $(typeof(data).name.wrapper)")

function rrtmgp_artifact(subfolder, file_name)
    artifact_name = "RRTMGPReferenceData"
    artifacts_file =
        joinpath(dirname(dirname(pathof(RRTMGP))), "Artifacts.toml")
    data_folder = joinpath(
        Pkg.Artifacts.ensure_artifact_installed(artifact_name, artifacts_file),
        artifact_name,
    )
    return Dataset(joinpath(data_folder, subfolder, file_name), "r")
end

abstract type AbstractRadiationMode end
struct GrayRadiation <: AbstractRadiationMode end
struct ClearSkyRadiation <: AbstractRadiationMode end
struct FullRadiation <: AbstractRadiationMode end
struct FullRadiationWithClearSkyDiagnostics <: AbstractRadiationMode end

abstract type AbstractInterpolation end
struct NoInterpolation <: AbstractInterpolation end
struct ArithmeticMean <: AbstractInterpolation end
struct GeometricMean <: AbstractInterpolation end
struct UniformZ <: AbstractInterpolation end
struct UniformP <: AbstractInterpolation end
struct BestFit <: AbstractInterpolation end

abstract type AbstractBottomExtrapolation end
struct SameAsInterpolation <: AbstractBottomExtrapolation end
struct UseSurfaceTempAtBottom <: AbstractBottomExtrapolation end
struct HydrostaticBottom <: AbstractBottomExtrapolation end

# NOTE: GeometricMean, UniformZ, and UniformP are all special cases of BestFit
#       (they assume different values for z, rather than using the true values),
#       but ArithmeticMean, HydrostaticBottom, and UseSurfaceTempAtBottom are
#       not consistent with BestFit.

requires_z(::Any) = false
requires_z(::Union{BestFit, HydrostaticBottom}) = true

uniform_z_p(t, t₁, t₂, p₁, p₂) = t₁ == t₂ ?
    p₁ * (p₂ / p₁)^(1/2) : p₁ * (p₂ / p₁)^(log(t / t₁) / log(t₂ / t₁))
best_fit_p(t, t₁, t₂, p₁, p₂, z, z₁, z₂) = t₁ == t₂ ?
    p₁ * (p₂ / p₁)^((z - z₁) / (z₂ - z₁)) :
    p₁ * (p₂ / p₁)^(log(t / t₁) / log(t₂ / t₁))

function interp!(::ArithmeticMean, t, p, tꜜ, tꜛ, pꜜ, pꜛ)
    @. t = (tꜜ + tꜛ) / 2
    @. p = (pꜜ + pꜛ) / 2
end
function interp!(::GeometricMean, t, p, tꜜ, tꜛ, pꜜ, pꜛ)
    @. t = sqrt(tꜜ * tꜛ)
    @. p = sqrt(pꜜ * pꜛ)
end
function interp!(::UniformZ, t, p, tꜜ, tꜛ, pꜜ, pꜛ)
    @. t = (tꜜ + tꜛ) / 2
    @. p = uniform_z_p(t, tꜜ, tꜛ, pꜜ, pꜛ)
end
function interp!(::UniformP, t, p, tꜜ, tꜛ, pꜜ, pꜛ)
    @. p = (pꜜ + pꜛ) / 2
    @. t = tꜜ * (tꜛ / tꜜ)^(log(p / pꜜ) / log(pꜛ / pꜜ)) # assume that pꜜ != pꜛ
end
function interp!(::BestFit, t, p, tꜜ, tꜛ, pꜜ, pꜛ, z, zꜜ, zꜛ)
    @. t = tꜜ + (tꜛ - tꜜ) * (z - zꜜ) / (zꜛ - zꜜ)
    @. p = best_fit_p(t, tꜜ, tꜛ, pꜜ, pꜛ, z, zꜜ, zꜛ)
end

function extrap!(::ArithmeticMean, t, p, t⁺, t⁺⁺, p⁺, p⁺⁺, tₛ, params)
    @. t = (3 * t⁺ - t⁺⁺) / 2
    @. p = (3 * p⁺ - p⁺⁺) / 2
end
function extrap!(::GeometricMean, t, p, t⁺, t⁺⁺, p⁺, p⁺⁺, tₛ, params)
    @. t = sqrt(t⁺^3 / t⁺⁺)
    @. p = sqrt(p⁺^3 / p⁺⁺)
end
function extrap!(::UniformZ, t, p, t⁺, t⁺⁺, p⁺, p⁺⁺, tₛ, params)
    @. t = (3 * t⁺ - t⁺⁺) / 2
    @. p = uniform_z_p(t, t⁺, t⁺⁺, p⁺, p⁺⁺)
end
function extrap!(::UniformP, t, p, t⁺, t⁺⁺, p⁺, p⁺⁺, tₛ, params)
    @. p = (3 * p⁺ - p⁺⁺) / 2
    @. t = t⁺ * (t⁺⁺ / t⁺)^(log(p / p⁺) / log(p⁺⁺ / p⁺)) # assume that p⁺ != p⁺⁺
end
function extrap!(::BestFit, t, p, t⁺, t⁺⁺, p⁺, p⁺⁺, tₛ, params, z, z⁺, z⁺⁺)
    @. t = t⁺ + (t⁺⁺ - t⁺) * (z - z⁺) / (z⁺⁺ - z⁺)
    @. p = best_fit_p(t, t⁺, t⁺⁺, p⁺, p⁺⁺, z, z⁺, z⁺⁺)
end
function extrap!(::UseSurfaceTempAtBottom, p, t, p⁺, p⁺⁺, t⁺, t⁺⁺, tₛ, params)
    FT = eltype(p)
    cₚ = FT(Planet.cp_d(params))
    R = FT(Planet.R_d(params))
    @. t = tₛ
    @. p = p⁺ * (t / t⁺)^(cₚ / R)
end
function extrap!(
    ::HydrostaticBottom,
    t,
    p,
    t⁺,
    t⁺⁺,
    p⁺,
    p⁺⁺,
    tₛ,
    params,
    z,
    z⁺,
    z⁺⁺,
)
    FT = eltype(p)
    g = FT(Planet.grav(params))
    cₚ = FT(Planet.cp_d(params))
    R = FT(Planet.R_d(params))
    @. t = t⁺ + g / cₚ * (z⁺ - z)
    @. p = p⁺ * (t / t⁺)^(cₚ / R)
end

struct RRTMGPModel{R, I, B, L, P, S, V}
    radiation_mode::R
    interpolation::I
    bottom_extrapolation::B
    add_isothermal_boundary_layer::Bool
    disable_lw::Bool
    disable_sw::Bool
    lookups::L
    params::P
    max_threads::Int
    solver::S
    views::V # user-friendly views into the solver
end

function Base.getproperty(model::RRTMGPModel, s::Symbol)
    if s in fieldnames(typeof(model))
        return getfield(model, s)
    else
        return getproperty(getfield(model, :views), s)
    end
end

function Base.propertynames(model::RRTMGPModel, private::Bool = false)
    names = propertynames(getfield(model, :views))
    return private ? (names..., fieldnames(typeof(model))...) : names
end

# This sets array .= value, but it allows array to be to be a CuArray while
# value is an Array (in which case broadcasting throws an error).
# TODO: Should this be parallelized?
set_array!(array, value::Real, symbol) = fill!(array, value)
function set_array!(array, value::AbstractArray{<:Real}, symbol)
    if ndims(array) == 2
        if size(value) == size(array)
            copyto!(array, value)
        elseif size(value) == (size(array, 1),)
            for col in eachcol(array)
                copyto!(col, value)
            end
        elseif size(value) == (1, size(array, 2))
            for (icol, col) in enumerate(eachcol(array))
                fill!(col, value[1, icol])
            end
        else
            error("expected $symbol to be an array of size $(size(array)), \
                   ($(size(array, 1)),), or (1, $(size(array, 2))); received \
                   an array of size $(size(value))")
        end
    else
        if size(value) == size(array)
            copyto!(array, value)
        else
            error("expected $symbol to be an array of size $(size(array)); \
                   received an array of size $(size(value))")
        end
    end
end

function set_and_save!(
    array,
    name,
    views,
    domain_nlay,
    extension_nlay,
    dict = nothing,
)
    domain_symbol = Symbol(name)
    extension_symbol = Symbol("extension_$name")

    if isnothing(dict)
        domain_value = NaN
    else
        if !(domain_symbol in keys(dict))
            throw(UndefKeywordError(domain_symbol))
        end
        domain_value = pop!(dict, domain_symbol)
    end

    if (
        (startswith(name, "center_") || startswith(name, "face_")) &&
        extension_nlay > 0
    )
        if isnothing(dict)
            extension_value = NaN
        else
            if !(extension_symbol in keys(dict))
                if domain_value isa Real
                    extension_value = domain_value
                else
                    throw(UndefKeywordError(extension_symbol))
                end
            end
            extension_value = pop!(dict, extension_symbol)
        end

        if startswith(name, "center_")
            domain_range = 1:domain_nlay
            extension_range =
                (domain_nlay + 1):(domain_nlay + extension_nlay)
        else # startswith(name, "face_")
            domain_range = 1:(domain_nlay + 1)
            extension_range =
                (domain_nlay + 2):(domain_nlay + extension_nlay + 1)
        end
        domain_view = view(array, domain_range, :)
        extension_view = view(array, extension_range, :)

        set_array!(domain_view, domain_value, domain_symbol)
        push!(views, (domain_symbol, domain_view))
        set_array!(extension_view, extension_value, extension_symbol)
        push!(views, (extension_symbol, extension_view))
    else
        set_array!(array, domain_value, domain_symbol)
        push!(views, (domain_symbol, array))
    end
end

function RRTMGPModel(
    params::AbstractEarthParameterSet;
    ncol::Int,
    domain_nlay::Int,
    extension_nlay::Int = 0,
    FT::Type{<:AbstractFloat} = Float64,
    DA::Type{<:AbstractArray} = RRTMGP.Device.array_type(),
    radiation_mode::AbstractRadiationMode = ClearSkyRadiation(),
    interpolation::AbstractInterpolation = NoInterpolation(),
    bottom_extrapolation::AbstractBottomExtrapolation = SameAsInterpolation(),
    disable_longwave::Bool = false,
    disable_shortwave::Bool = false,
    use_one_scalar_mode::Bool = false,
    use_pade_cloud_optics_mode::Bool = false,
    use_global_means_for_trace_gases::Bool = false,
    add_isothermal_boundary_layer::Bool = false,
    max_threads::Int = 256,
    kwargs...,
)
    # turn kwargs into a Dict, so that values can be dynamically popped from it
    dict = Dict(kwargs)

    if disable_longwave && disable_shortwave
        error("either longwave or shortwave fluxes must be enabled")
    end
    if use_one_scalar_mode && !disable_shortwave
        @warn "upward shortwave fluxes are not computed when \
               use_one_scalar_mode is true"
    end
    if use_pade_cloud_optics_mode && (
        radiation_mode isa GrayRadiation ||
        radiation_mode isa ClearSkyRadiation
    )
        @warn "use_pade_cloud_optics_mode is ignored when using GrayRadiation \
               or ClearSkyRadiation"
    end
    if use_global_means_for_trace_gases && radiation_mode isa GrayRadiation
        @warn "use_global_means_for_trace_gases is ignored when using \
               GrayRadiation"
    end

    if (
        :center_pressure in keys(dict) &&
        :center_temperature in keys(dict) &&
        :face_pressure in keys(dict) &&
        :face_temperature in keys(dict)
    )
        if !(interpolation isa NoInterpolation)
            @warn "interpolation is ignored if both center and face pressures/\
                   temperatures are specified"
        end
        implied_values = :none
    elseif (
        :center_pressure in keys(dict) &&
        :center_temperature in keys(dict) &&
        !(:face_pressure in keys(dict)) &&
        !(:face_temperature in keys(dict))
    )
        if interpolation isa NoInterpolation
            error("interpolation cannot be NoInterpolation if only center \
                   pressures/temperatures are specified")
        end
        implied_values = :face
    elseif (
        !(:center_pressure in keys(dict)) &&
        !(:center_temperature in keys(dict)) &&
        :face_pressure in keys(dict) &&
        :face_temperature in keys(dict)
    )
        if interpolation isa NoInterpolation
            error("interpolation cannot be NoInterpolation if only face \
                   pressures/temperatures are specified")
        end
        implied_values = :center
    else
        error("please specify either center_pressure and center_temperature, \
               or face_pressure and face_temperature, or all four values")
    end

    if implied_values != :face
        if !(bottom_extrapolation isa SameAsInterpolation)
            @warn "bottom_extrapolation is ignored if face_pressure and \
                   face_temperature are specified"
        end
    end

    lookups = NamedTuple()
    views = []

    nlay = domain_nlay + extension_nlay + Int(add_isothermal_boundary_layer)
    op_symbol = use_one_scalar_mode ? :OneScalar : :TwoStream
    t = (views, domain_nlay, extension_nlay)

    if disable_longwave
        src_lw = flux_lw = fluxb_lw = bcs_lw = nothing
    else
        if radiation_mode isa GrayRadiation
            nbnd_lw = ngpt_lw = 1
        else
            ds_lw = rrtmgp_artifact("lookup_tables", "clearsky_lw.nc")
            lookup_lw, idx_gases =
                RRTMGP.LookUpTables.LookUpLW(ds_lw, Int, FT, DA)
            close(ds_lw)
            lookups = (; lookups..., lookup_lw, idx_gases)

            nbnd_lw = lookup_lw.n_bnd
            ngpt_lw = lookup_lw.n_gpt
            ngas = lookup_lw.n_gases

            if !(radiation_mode isa ClearSkyRadiation)
                ds_lw_cld = rrtmgp_artifact("lookup_tables", "cloudysky_lw.nc")
                lookup_lw_cld = RRTMGP.LookUpTables.LookUpCld(
                    ds_lw_cld,
                    Int,
                    FT,
                    DA,
                    !use_pade_cloud_optics_mode,
                )
                close(ds_lw_cld)
                lookups = (; lookups..., lookup_lw_cld)
            end
        end

        src_lw =
            RRTMGP.Sources.source_func_longwave(FT, ncol, nlay, op_symbol, DA)
        flux_lw = RRTMGP.Fluxes.FluxLW(ncol, nlay, FT, DA)
        fluxb_lw = radiation_mode isa GrayRadiation ? nothing :
            RRTMGP.Fluxes.FluxLW(ncol, nlay, FT, DA)
        set_and_save!(flux_lw2.flux_up, "face_lw_flux_up", t...)
        set_and_save!(flux_lw2.flux_dn, "face_lw_flux_dn", t...)
        set_and_save!(flux_lw2.flux_net, "face_lw_flux", t...)
        if radiation_mode isa FullRadiationWithClearSkyDiagnostics
            flux_lw2 = RRTMGP.Fluxes.FluxLW(ncol, nlay, FT, DA)
            set_and_save!(flux_lw2.flux_up, "face_clear_lw_flux_up", t...)
            set_and_save!(flux_lw2.flux_dn, "face_clear_lw_flux_dn", t...)
            set_and_save!(flux_lw2.flux_net, "face_clear_lw_flux", t...)
        end

        sfc_emis = DA{FT}(undef, nbnd_lw, ncol)
        set_and_save!(sfc_emis, "surface_emissivity", t..., dict)
        name = "top_of_atmosphere_lw_flux_dn"
        if Symbol(name) in keys(dict)
            inc_flux = DA{FT}(undef, ncol, ngpt_lw)
            set_and_save!(transpose(inc_flux), name, t..., dict)
        else
            inc_flux = nothing
        end
        bcs_lw = RRTMGP.BCs.LwBCs(sfc_emis, inc_flux)
    end

    if disable_shortwave
        src_sw = flux_sw = fluxb_sw = bcs_sw = nothing
    else
        if radiation_mode isa GrayRadiation
            nbnd_sw = ngpt_sw = 1
        else
            ds_sw = rrtmgp_artifact("lookup_tables", "clearsky_sw.nc")
            lookup_sw, idx_gases =
                RRTMGP.LookUpTables.LookUpSW(ds_sw, Int, FT, DA)
            close(ds_sw)
            lookups = (; lookups..., lookup_sw, idx_gases)

            nbnd_sw = lookup_sw.n_bnd
            ngpt_sw = lookup_sw.n_gpt
            ngas = lookup_sw.n_gases

            if !(radiation_mode isa ClearSkyRadiation)
                ds_sw_cld = rrtmgp_artifact("lookup_tables", "cloudysky_sw.nc")
                lookup_sw_cld = RRTMGP.LookUpTables.LookUpCld(
                    ds_sw_cld,
                    Int,
                    FT,
                    DA,
                    !use_pade_cloud_optics_method,
                )
                close(ds_sw_cld)
                lookups = (; lookups..., lookup_sw_cld)
            end
        end

        src_sw =
            RRTMGP.Sources.source_func_shortwave(FT, ncol, nlay, op_symbol, DA)
        flux_sw = RRTMGP.Fluxes.FluxSW(ncol, nlay, FT, DA)
        fluxb_sw = radiation_mode isa GrayRadiation ? nothing :
            RRTMGP.Fluxes.FluxSW(ncol, nlay, FT, DA)
        set_and_save!(flux_sw.flux_up, "face_sw_flux_up", t...)
        set_and_save!(flux_sw.flux_dn, "face_sw_flux_dn", t...)
        set_and_save!(flux_sw.flux_dn_dir, "face_sw_direct_flux_dn", t...)
        set_and_save!(flux_sw.flux_net, "face_sw_flux", t...)
        if radiation_mode isa FullRadiationWithClearSkyDiagnostics
            flux_sw2 = RRTMGP.Fluxes.FluxSW(ncol, nlay, FT, DA)
            set_and_save!(flux_sw2.flux_up, "face_clear_sw_flux_up", t...)
            set_and_save!(flux_sw2.flux_dn, "face_clear_sw_flux_dn", t...)
            set_and_save!(
                flux_sw2.flux_dn_dir,
                "face_clear_sw_direct_flux_dn",
                t...,
            )
            set_and_save!(flux_sw2.flux_net, "face_clear_sw_flux", t...)
        end

        zenith = DA{FT}(undef, ncol)
        set_and_save!(zenith, "solar_zenith_angle", t..., dict)
        toa_flux = DA{FT}(undef, ncol)
        set_and_save!(toa_flux, "weighted_irradiance", t..., dict)
        sfc_alb_direct = DA{FT}(undef, nbnd_sw, ncol)
        set_and_save!(sfc_alb_direct, "direct_sw_surface_albedo", t..., dict)
        sfc_alb_diffuse = DA{FT}(undef, nbnd_sw, ncol)
        set_and_save!(sfc_alb_diffuse, "diffuse_sw_surface_albedo", t..., dict)
        name = "top_of_atmosphere_diffuse_sw_flux_dn"
        if Symbol(name) in keys(init_dict)
            @warn "incoming diffuse shortwave fluxes are not yet implemented \
                   in RRTMGP.jl; the value of $name will be ignored"
            inc_flux_diffuse = DA{FT}(undef, ncol, ngpt_sw)
            set_and_save!(transpose(inc_flux_diffuse), name, t..., dict)
        else
            inc_flux_diffuse = nothing
        end
        bcs_sw = RRTMGP.BCs.SwBCs(
            zenith,
            toa_flux,
            sfc_alb_direct,
            inc_flux_diffuse,
            sfc_alb_diffuse,
        )
    end

    if disable_longwave
        set_and_save!(flux_sw.flux_net, "face_flux", t...)
        if radiation_mode isa FullRadiationWithClearSkyDiagnostics
            set_and_save!(flux_sw2.flux_net, "face_clear_flux", t...)
        end
    elseif disable_shortwave
        set_and_save!(flux_lw.flux_net, "face_flux", t...)
        if radiation_mode isa FullRadiationWithClearSkyDiagnostics
            set_and_save!(flux_lw2.flux_net, "face_clear_flux", t...)
        end
    else
        set_and_save!(similar(flux_lw.flux_net), "face_flux", t...)
        if radiation_mode isa FullRadiationWithClearSkyDiagnostics
            set_and_save!(similar(flux_lw2.flux_net), "face_clear_flux", t...)
        end
        if !(radiation_mode isa GrayRadiation)
            @assert lookup_lw.n_gases == lookup_sw.n_gases
            @assert lookup_lw.p_ref_min == lookup_sw.p_ref_min
        end
    end

    p_lay = DA{FT}(undef, nlay, ncol)
    p_lev = DA{FT}(undef, nlay + 1, ncol)
    t_lay = DA{FT}(undef, nlay, ncol)
    t_lev = DA{FT}(undef, nlay + 1, ncol)
    if implied_values != :center
        set_and_save!(p_lay, "center_pressure", t..., dict)
        set_and_save!(t_lay, "center_temperature", t..., dict)
    end
    if implied_values != :face
        set_and_save!(p_lev, "face_pressure", t..., dict)
        set_and_save!(t_lev, "face_temperature", t..., dict)
    end
    t_sfc = DA{FT}(undef, ncol)
    set_and_save!(t_sfc, "surface_temperature", t..., dict)

    if radiation_mode isa GrayRadiation
        z_lev = DA{FT}(undef, nlay + 1, ncol) # TODO: z_lev required but unused

        # lapse_rate is a constant, so don't use set_and_save! to get it
        if !(:lapse_rate in keys(dict))
            throw(UndefKeywordError(:lapse_rate))
        end
        α = pop!(dict, :lapse_rate)
        if !(α isa Real)
            error("lapse_rate must be a Real")
        end

        d0 = DA{FT}(undef, ncol)
        set_and_save!(d0, "optical_thickness_parameter", t..., dict)
        as = RRTMGP.AtmosphericStates.GrayAtmosphericState(
            p_lay,
            p_lev,
            t_lay,
            t_lev,
            z_lev,
            t_sfc,
            α,
            d0,
            nlay,
            ncol,
        )
    else
        if !(:latitude in keys(init_dict))
            lon = lat = nothing
        else
            lon = DA{FT}(undef, ncol) # TODO: lon required but unused
            lat = DA{FT}(undef, ncol)
            set_and_save!(lat, "latitude", t..., dict)
        end

        vmr_str = "volume_mixing_ratio_"
        gas_names = filter(
            gas_name -> !(gas_name in ("h2o", "h2o_frgn", "h2o_self", "o3")),
            keys(idx_gases),
        )
        # TODO: This gives the wrong types for CUDA 3.4 and above.
        # gm = use_global_means_for_trace_gases
        # vmr = RRTMGP.Vmrs.init_vmr(ngas, nlay, ncol, FT, DA; gm)
        if use_global_means_for_trace_gases
            vmr = RRTMGP.Vmrs.VmrGM(
                DA{FT}(undef, nlay, ncol),
                DA{FT}(undef, nlay, ncol),
                DA{FT}(undef, ngas),
            )
            vmr.vmr .= 0 # TODO: do we need this?
            set_and_save!(vmr.vmr_h2o, "center_$(vmr_str)h2o", t..., dict)
            set_and_save!(vmr.vmr_o3, "center_$(vmr_str)o3", t..., dict)
            for gas_name in gas_names
                gas_view = view(vmr.vmr, idx_gases[gas_name])
                set_and_save!(gas_view, "$vmr_str$gas_name", t..., dict)
            end
        else
            vmr = RRTMGP.Vmrs.Vmr(DA{FT}(undef, nlay, ncol, ngas))
            for gas_name in ["h2o", "o3", gas_names...]
                gas_view = view(vmr.vmr, :, :, idx_gases[gas_name])
                set_and_save!(gas_view, "center_$vmr_str$gas_name", t..., dict)
            end
        end

        if radiation_mode isa ClearSkyRadiation
            cld_r_eff_liq = cld_r_eff_ice = nothing
            cld_path_liq = cld_path_ice = cld_mask = nothing
            ice_rgh = 1
        else
            cld_r_eff_liq = DA{FT}(undef, nlay, ncol)
            name = "center_cloud_liquid_effective_radius"
            set_and_save!(cld_r_eff_liq, name, t..., dict)
            cld_r_eff_ice = DA{FT}(undef, nlay, ncol)
            name = "center_cloud_ice_effective_radius"
            set_and_save!(cld_r_eff_ice, name, t..., dict)
            cld_path_liq = DA{FT}(undef, nlay, ncol)
            name = "center_cloud_liquid_water_path"
            set_and_save!(cld_path_liq, name, t..., dict)
            cld_path_ice = DA{FT}(undef, nlay, ncol)
            name = "center_cloud_ice_water_path"
            set_and_save!(cld_path_ice, name, t..., dict)
            cld_mask = DA{Bool}(undef, nlay, ncol)
            set_and_save!(cld_mask, "center_cloud_boolean_mask", t..., dict)

            # ice_roughness is a constant, so don't use set_and_save! to get it
            if !(:ice_roughness in keys(dict))
                throw(UndefKeywordError(:ice_roughness))
            end
            ice_rgh = pop!(dict, :ice_roughness)
            if !(ice_rgh in (1, 2, 3))
                error("ice_roughness must be either 1, 2, or 3")
            end
        end

        as = RRTMGP.AtmosphericStates.AtmosphericState(
            lon,
            lat,
            p_lay,
            p_lev,
            t_lay,
            t_lev,
            t_sfc,
            DA{FT}(undef, nlay, ncol), # col_dry array is for internal use only
            vmr,
            cld_r_eff_liq,
            cld_r_eff_ice,
            cld_path_liq,
            cld_path_ice,
            cld_mask,
            ice_rgh,
            nlay,
            ncol,
            ngas,
        )
    end

    op_type = use_one_scalar ? RRTMGP.Optics.OneScalar : RRTMGP.Optics.TwoStream
    solver = RRTMGP.RTE.Solver(
        as,
        op_type(FT, ncol, nlay, DA),
        src_lw,
        src_sw,
        bcs_lw,
        bcs_sw,
        fluxb_lw,
        fluxb_sw,
        flux_lw,
        flux_sw,
    )

    if requires_z(interpolation) || requires_z(bottom_extrapolation)
        z_lay = DA{FT}(undef, nlay, ncol)
        set_and_save!(z_lay, "center_z", t..., dict)
        z_lev = DA{FT}(undef, nlay + 1, ncol)
        set_and_save!(z_lev, "face_z", t..., dict)
    end

    if length(dict) > 0
        @warn string(
            "unused keyword argument",
            length(dict) == 1 ? " " : "s ",
            join(keys(dict), ", ", length(dict) == 2 ? " and " : ", and "),
        )
    end

    return RRTMGPModel(
        radiation_mode,
        interpolation,
        bottom_extrapolation,
        add_isothermal_boundary_layer,
        disable_lw,
        disable_sw,
        lookups,
        params,
        max_threads,
        solver,
        NamedTuple(views),
    )
end

get_p_min(model) = get_p_min(model.solver.as, model.lookups)
get_p_min(as::RRTMGP.AtmosphericStates.GrayAtmosphericState, lookups) =
    zero(eltype(as.p_lay))
get_p_min(as::RRTMGP.AtmosphericStates.AtmosphericState, lookups) =
    lookups[1].p_ref_min

function update_implied_values!(model)
    (; p_lay, p_lev, t_lay, t_lev) = model.solver.as
    if requires_z(model.interpolation) || requires_z(model.bottom_extrapolation)
        z_lay = parent(model.center_z)
        z_lev = parent(model.face_z)
    end
    nlay = size(p_lay, 1) - Int(model.add_isothermal_boundary_layer)
    if !(:center_pressure in propertynames(model))
        if requires_z(model.interpolation)
            z_args = (
                view(z_lay, 1:nlay, :),
                view(z_lev, 1:nlay, :),
                view(z_lev, 2:(nlay + 1), :),
            )
        else
            z_args = ()
        end
        interp!(
            model.interpolation,
            view(t_lay, 1:nlay, :),
            view(p_lay, 1:nlay, :),
            view(t_lev, 1:nlay, :),
            view(t_lev, 2:(nlay + 1), :),
            view(p_lev, 1:nlay, :),
            view(p_lev, 2:(nlay + 1), :),
            z_args...,
        )
    elseif !(:face_pressure in propertynames(model))
        if requires_z(model.interpolation)
            z_args = (
                view(z_lev, 2:nlay, :),
                view(z_lay, 1:(nlay - 1), :),
                view(z_lay, 2:nlay, :),
            )
        else
            z_args = ()
        end
        interp!(
            model.interpolation,
            view(t_lev, 2:nlay, :),
            view(p_lev, 2:nlay, :),
            view(t_lay, 1:(nlay - 1), :),
            view(t_lay, 2:nlay, :),
            view(p_lay, 1:(nlay - 1), :),
            view(p_lay, 2:nlay, :),
            z_args...,
        )
        if requires_z(model.interpolation)
            z_args = (
                view(z_lev, nlay + 1, :),
                view(z_lay, nlay, :),
                view(z_lay, nlay - 1, :),
            )
        else
            z_args = ()
        end
        extrap!(
            model.interpolation,
            view(t_lev, nlay + 1, :),
            view(p_lev, nlay + 1, :),
            view(t_lay, nlay, :),
            view(t_lay, nlay - 1, :),
            view(p_lay, nlay, :),
            view(p_lay, nlay - 1, :),
            model.solver.as.t_sfc,
            model.params,
            z_args...,
        )
        bottom_extrapolation =
            model.bottom_extrapolation isa SameAsInterpolation ?
            model.interpolation : model.bottom_extrapolation
        if requires_z(bottom_extrapolation)
            z_args = (view(z_lev, 1, :), view(z_lay, 1, :), view(z_lay, 2, :))
        else
            z_args = ()
        end
        extrap!(
            bottom_extrapolation,
            view(t_lev, 1, :),
            view(p_lev, 1, :),
            view(t_lay, 1, :),
            view(t_lay, 2, :),
            view(p_lay, 1, :),
            view(p_lay, 2, :),
            model.solver.as.t_sfc,
            model.params,
            z_args...,
        )
    end
end

function update_boundary_layer!(model)
    as = model.solver.as
    p_min = get_p_min(model)
    @views as.p_lay[end, :] .= (as.p_lev[end - 1, :] .+ p_min) ./ 2
    @views as.p_lev[end, :] .= p_min
    @views as.t_lay[end, :] .= as.t_lev[end - 1, :]
    @views as.t_lev[end, :] .= as.t_lev[end - 1, :]
    update_boundary_layer_vmr!(model.radiation_mode, as)
end
update_boundary_layer_vmr!(::GrayRadiation, as) = nothing
update_boundary_layer_vmr!(radiation_mode, as) =
    update_boundary_layer_vmr!(as.vmr)
function update_boundary_layer_vmr!(vmr::RRTMGP.Vmrs.VmrGM)
    @views vmr.vmr_h2o[end, :] .= vmr.vmr_h2o[end - 1, :]
    @views vmr.vmr_o3[end, :] .= vmr.vmr_o3[end - 1, :]
end
update_boundary_layer_vmr!(vmr::RRTMGP.Vmrs.Vmr) =
    @views vmr.vmr[end, :, :] .= vmr.vmr[end - 1, :, :]

update_concentrations!(::GrayRadiation, model) = nothing
update_concentrations!(radiation_mode, model) =
    RRTMGP.Optics.compute_col_dry!(
        model.solver.as.p_lev,
        model.solver.as.col_dry,
        model.params,
        get_vmr_h2o(model.solver.as.vmr, model.lookups.idx_gases),
        model.solver.as.lat,
        model.max_threads,
    )
get_vmr_h2o(vmr::RRTMGP.Vmrs.VmrGM, idx_gases) = vmr.vmr_h2o
get_vmr_h2o(vmr::RRTMGP.Vmrs.Vmr, idx_gases) =
    view(vmr.vmr, :, :, idx_gases["h2o"])

update_lw_fluxes!(::GrayRadiation, model) =
    RRTMGP.RTESolver.solve_lw!(model.solver, model.max_threads)
update_lw_fluxes!(::ClearSkyRadiation, model) =
    RRTMGP.RTESolver.solve_lw!(
        model.solver,
        model.max_threads,
        model.lookups.lookup_lw,
    )
update_lw_fluxes!(::FullRadiation, model) =
    RRTMGP.RTESolver.solve_lw!(
        model.solver,
        model.max_threads,
        model.lookups.lookup_lw,
        model.lookups.lookup_lw_cld,
    )
function update_lw_fluxes!(::FullRadiationWithClearSkyDiagnostics, model)
    RRTMGP.RTESolver.solve_lw!(
        model.solver,
        model.max_threads,
        model.lookups.lookup_lw,
    )
    parent(model.face_clear_lw_flux_up) .= parent(model.face_lw_flux_up)
    parent(model.face_clear_lw_flux_dn) .= parent(model.face_lw_flux_dn)
    parent(model.face_clear_lw_flux) .= parent(model.face_lw_flux)
    RRTMGP.RTESolver.solve_lw!(
        model.solver,
        model.max_threads,
        model.lookups.lookup_lw,
        model.lookups.lookup_lw_cld,
    )
end

update_sw_fluxes!(::GrayRadiation, model) =
    RRTMGP.RTESolver.solve_sw!(model.solver, model.max_threads)
update_sw_fluxes!(::ClearSkyRadiation, model) =
    RRTMGP.RTESolver.solve_sw!(
        model.solver,
        model.max_threads,
        model.lookups.lookup_sw,
    )
update_sw_fluxes!(::FullRadiation, model) =
    RRTMGP.RTESolver.solve_sw!(
        model.solver,
        model.max_threads,
        model.lookups.lookup_sw,
        model.lookups.lookup_sw_cld,
    )
function update_sw_fluxes!(::FullRadiationWithClearSkyDiagnostics, model)
    RRTMGP.RTESolver.solve_sw!(
        model.solver,
        model.max_threads,
        model.lookups.lookup_sw,
    )
    parent(model.face_clear_sw_flux_up) .= parent(model.face_sw_flux_up)
    parent(model.face_clear_sw_flux_dn) .= parent(model.face_sw_flux_dn)
    parent(model.face_clear_sw_direct_flux_dn) .=
        parent(model.face_sw_direct_flux_dn)
    parent(model.face_clear_sw_flux) .= parent(model.face_sw_flux)
    RRTMGP.RTESolver.solve_sw!(
        model.solver,
        model.max_threads,
        model.lookups.lookup_sw,
        model.lookups.lookup_sw_cld,
    )
end

update_net_flux!(radiation_mode, model) =
    parent(model.face_flux) .=
        parent(model.face_lw_flux) .+ parent(model.face_sw_flux)
function update_net_flux!(::FullRadiationWithClearSkyDiagnostics, model)
    parent(model.face_clear_flux) .=
        parent(model.face_clear_lw_flux) .+ parent(model.face_clear_sw_flux)
    parent(model.face_flux) .=
        parent(model.face_lw_flux) .+ parent(model.face_sw_flux)
end

function update_fluxes!(model)
    model.implied_values != :none && update_implied_values!(model)
    model.add_isothermal_boundary_layer && update_boundary_layer!(model)

    p_min = get_p_min(model)
    @. model.solver.as.p_lay = max(model.solver.as.p_lay, p_min)
    @. model.solver.as.p_lev = max(model.solver.as.p_lev, p_min)

    update_concentrations!(model.radiation_mode, model)
    !model.disable_lw && update_lw_fluxes!(model.radiation_mode, model)
    !model.disable_sw && update_sw_fluxes!(model.radiation_mode, model)
    !(model.disable_lw || model.disable_sw) &&
        update_net_fluxes!(model.radiation_mode, model)
    return model.face_flux
end

# Overriding the definition of optical depth for GrayRadiation.
function τ_lw_gray(p, pꜜ, pꜛ, p₀, τ₀)
    FT = eltype(p)
    f = FT(0.2)
    α = 1
    return α * τ₀ * (f / p₀ + 4 * (1 - f) / p₀ * (p / p₀)^3) * (pꜜ - pꜛ)
end
τ_sw_gray(p, pꜜ, pꜛ, p₀, τ₀) = 2 * τ₀ * (p / p₀) / p₀ * (pꜜ - pꜛ)
# Note: the original value for both of these functions was
#     abs(α * τ₀ * (p / p₀)^α / p * (p⁺ - pꜜ)),
# where α is the lapse rate and τ₀ is the optical thickness parameter.

import RRTMGP.Optics: compute_optical_props_kernel!
function compute_optical_props_kernel!(
    op::RRTMGP.Optics.AbstractOpticalProps{FT},
    as::RRTMGP.AtmosphericStates.GrayAtmosphericState{FT},
    glaycol,
    source::RRTMGP.Sources.AbstractSourceLW{FT},
) where {FT<:AbstractFloat}
    compute_optical_props_kernel_lw!(op, as, glaycol)
    RRTMGP.Optics.compute_sources_gray_kernel!(source, as, glaycol)
end
function compute_optical_props_kernel_lw!(
    op::RRTMGP.Optics.AbstractOpticalProps{FT},
    as::RRTMGP.AtmosphericStates.GrayAtmosphericState{FT},
    glaycol,
) where {FT<:AbstractFloat}
    glay, gcol = glaycol
    (; p_lay, p_lev, d0) = as
    @inbounds op.τ[glay, gcol] = τ_lw_gray(
        p_lay[glay, gcol],
        p_lev[glay, gcol],
        p_lev[glay + 1, gcol],
        p_lev[1, gcol],
        d0[gcol],
    )
    if op isa RRTMGP.Optics.TwoStream
        op.ssa[glaycol...] = FT(0)
        op.g[glaycol...] = FT(0)
    end
end
function compute_optical_props_kernel!(
    op::RRTMGP.Optics.AbstractOpticalProps{FT},
    as::RRTMGP.AtmosphericStates.GrayAtmosphericState{FT},
    glaycol,
) where {FT<:AbstractFloat}
    glay, gcol = glaycol
    (; p_lay, p_lev) = as
    @inbounds op.τ[glay, gcol] = τ_sw_gray(
        p_lay[glay, gcol],
        p_lev[glay, gcol],
        p_lev[glay + 1, gcol],
        p_lev[1, gcol],
        FT(0.22),
    )
    if op isa RRTMGP.Optics.TwoStream
        op.ssa[glaycol...] = FT(0)
        op.g[glaycol...] = FT(0)
    end
end
