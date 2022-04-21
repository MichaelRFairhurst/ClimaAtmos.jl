using Dates: Second, DateTime
using Insolation: instantaneous_zenith_angle
using CLIMAParameters: AbstractEarthParameterSet, Planet, astro_unit

include("../rrtmgp_model.jl")

function radiation_cache(Y, params)
    input_data = rrtmgp_artifact("atmos_state", "clearsky_as.nc")
    rrtmgp_model = RRTMGPModel(
        params;
        FT = Float64,
        ncol = length(Spaces.all_nodes(axes(Spaces.level(Y.c, 1)))),
        domain_nlay = Spaces.nlevels(axes(Y.c)),
        radiation_mode = ClearSkyRadiation(),
        interpolation = ArithmeticMean(),
        bottom_extrapolation = SameAsInterpolation(),
        use_global_means_for_well_mixed_gases = true,
        add_isothermal_boundary_layer = true,
        center_pressure = NaN, # initialize in tendency
        center_temperature = NaN, # initialize in tendency
        surface_temperature = 280,
        surface_emissivity = mean(input_data["surface_emissivity"]),
        direct_sw_surface_albedo = mean(input_data["surface_albedo"]),
        diffuse_sw_surface_albedo = mean(input_data["surface_albedo"]),
        solar_zenith_angle = NaN, # initialize in tendency
        weighted_irradiance = NaN, # initialize in tendency
        center_volume_mixing_ratio_h2o = NaN, # initialize in tendency
        center_volume_mixing_ratio_o3 = mean(input_data["ozone"]),
        volume_mixing_ratio_co2 = mean(input_data["carbon_dioxide_GM"]),
        volume_mixing_ratio_n2o = mean(input_data["nitrous_oxide_GM"]),
        volume_mixing_ratio_co = mean(input_data["carbon_monoxide_GM"]),
        volume_mixing_ratio_ch4 = mean(input_data["methane_GM"]),
        volume_mixing_ratio_o2 = mean(input_data["oxygen_GM"]),
        volume_mixing_ratio_n2 = mean(input_data["nitrogen_GM"]),
        volume_mixing_ratio_ccl4 = mean(input_data["carbon_tetrachloride_GM"]),
        volume_mixing_ratio_cfc11 = mean(input_data["cfc11_GM"]),
        volume_mixing_ratio_cfc12 = mean(input_data["cfc12_GM"]),
        volume_mixing_ratio_cfc22 = mean(input_data["hcfc22_GM"]),
        volume_mixing_ratio_hfc143a = mean(input_data["hfc143a_GM"]),
        volume_mixing_ratio_hfc125 = mean(input_data["hfc125_GM"]),
        volume_mixing_ratio_hfc23 = mean(input_data["hfc23_GM"]),
        volume_mixing_ratio_hfc32 = mean(input_data["hfc32_GM"]),
        volume_mixing_ratio_hfc134a = mean(input_data["hfc134a_GM"]),
        volume_mixing_ratio_cf4 = mean(input_data["cf4_GM"]),
        volume_mixing_ratio_no2 = 0, # not available in input_data
        latitude =
            field2array(Fields.coordinate_field(Spaces.level(Y.c, 1)).lat),
        # center_z = field2array(Fields.coordinate_field(Y.c).z),
        # face_z = field2array(Fields.coordinate_field(Y.f).z),
    )
    close(input_data)
    return (;
        ᶜT = similar(Y.c, FT),
        ᶜvmr_h2o = similar(Y.c, FT),
        insolation_tuple = similar(Spaces.level(Y.c, 1), Tuple{FT, FT, FT}),
        zenith_angle = similar(Spaces.level(Y.c, 1), FT),
        weighted_irradiance = similar(Spaces.level(Y.c, 1), FT),
        ᶠradiation_flux = similar(Y.f, FT),
        rrtmgp_model,
    )
end
function rrtmgp_model_tendency!(Yₜ, Y, p, t)
    (; ᶜts, ᶜp, params) = p # assume ᶜts and ᶜp have been updated
    (; ᶜT, ᶜvmr_h2o, insolation_tuple, zenith_angle, weighted_irradiance) = p
    (; ᶠradiation_flux, rrtmgp_model) = p

    date_time = DateTime(2022) + Second(round(Int, t)) # t secs into 2022
    max_zenith_angle = FT(π) / 2 - eps(FT)
    irradiance = FT(Planet.tot_solar_irrad(params))
    au = FT(astro_unit())

    @. ᶜT = TD.air_temperature(ᶜts)
    @. ᶜvmr_h2o = TD.vol_vapor_mixing_ratio(params, TD.PhasePartition(ᶜts))
    bottom_coords = Fields.coordinate_field(Spaces.level(Y.c, 1))
    @. insolation_tuple = instantaneous_zenith_angle(
        date_time,
        bottom_coords.long,
        bottom_coords.lat,
        params,
    ) # each tuple contains (zenith angle, azimuthal angle, earth-sun distance)
    @. zenith_angle = min(first(insolation_tuple), max_zenith_angle)
    @. weighted_irradiance = irradiance * (au / last(insolation_tuple))^2

    rrtmgp_model.center_pressure .= field2array(ᶜp)
    rrtmgp_model.center_temperature .= field2array(ᶜT)
    rrtmgp_model.center_volume_mixing_ratio_h2o .= field2array(ᶜvmr_h2o)
    rrtmgp_model.solar_zenith_angle .= field2array(zenith_angle)
    rrtmgp_model.weighted_irradiance .= field2array(weighted_irradiance)
    update_fluxes!(rrtmgp_model)

    field2array(ᶠradiation_flux) .= rrtmgp_model.face_flux
    if :ρe in propertynames(Y.c)
        @. Yₜ.c.ρe -= ᶜdivᵥ(Geometry.WVector(ᶠradiation_flux))
    elseif :ρe_int in propertynames(Y.c)
        @. Yₜ.c.ρe_int -= ᶜdivᵥ(Geometry.WVector(ᶠradiation_flux))
    elseif :ρθ in propertynames(Y.c)
        error("rrtmgp_model_tendency! not implemented for ρθ")
    end
end
