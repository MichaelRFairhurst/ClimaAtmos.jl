using Statistics: mean
using Dierckx: Spline1D
using Dates: Second, DateTime
using Insolation: instantaneous_zenith_angle
using CLIMAParameters: AbstractEarthParameterSet, Planet, astro_unit

include("../rrtmgp_model.jl")

function rrtmgp_model_cache(
    Y,
    params;
    radiation_mode = ClearSkyRadiation(),
    interpolation = BestFit(),
    bottom_extrapolation = SameAsInterpolation(),
    idealized_insolation = true,
)
    bottom_coords = Fields.coordinate_field(Spaces.level(Y.c, 1))
    if eltype(bottom_coords) <: Geometry.LatLongZPoint
        latitude = field2array(bottom_coords.lat)
    else
        latitude = field2array(zero(bottom_coords.z)) # flat space is on Equator
    end
    input_data = rrtmgp_artifact("atmos_state", "clearsky_as.nc")
    if radiation_mode isa GrayRadiation
        kwargs = (;
            lapse_rate = 3.5,
            optical_thickness_parameter = (@. (
                (300 + 60 * (FT(1 / 3) - sind(latitude)^2)) / 200
            )^4 - 1),
        )
    else
        # the pressure and ozone concentrations are provided for each of 100
        # sites, which we average across
        n = input_data.dim["layer"]
        input_center_pressure =
            vec(mean(reshape(input_data["pres_layer"][:, :], n, :); dims = 2))
        # the first values along the third dimension of the ozone concentration
        # data are the present-day values
        input_center_volume_mixing_ratio_o3 =
            vec(mean(reshape(input_data["ozone"][:, :, 1], n, :); dims = 2))

        # interpolate the ozone concentrations to our initial pressures (set the
        # kinetic energy to 0 when computing the pressure using total energy)
        pressure2ozone =
            Spline1D(input_center_pressure, input_center_volume_mixing_ratio_o3)
        if :ρθ in propertynames(Y.c)
            ᶜts = @. thermo_state_ρθ(Y.c.ρθ, Y.c, params)
        elseif :ρe in propertynames(Y.c)
            ᶜΦ = FT(Planet.grav(params)) .* Fields.coordinate_field(Y.c).z
            ᶜts = @. thermo_state_ρe(Y.c.ρe, Y.c, 0, ᶜΦ, params)
        elseif :ρe_int in propertynames(Y.c)
            ᶜts = @. thermo_state_ρe_int(Y.c.ρe_int, Y.c, params)
        end
        ᶜp = @. TD.air_pressure(params, ᶜts)
        center_volume_mixing_ratio_o3 = field2array(@. FT(pressure2ozone(ᶜp)))

        # the first value for each global mean volume mixing ratio is the
        # present-day value
        input_vmr(name) =
            input_data[name][1] * parse(FT, input_data[name].attrib["units"])
        kwargs = (;
            use_global_means_for_well_mixed_gases = true,
            center_volume_mixing_ratio_h2o = NaN, # initialize in tendency
            center_volume_mixing_ratio_o3,
            volume_mixing_ratio_co2 = input_vmr("carbon_dioxide_GM"),
            volume_mixing_ratio_n2o = input_vmr("nitrous_oxide_GM"),
            volume_mixing_ratio_co = input_vmr("carbon_monoxide_GM"),
            volume_mixing_ratio_ch4 = input_vmr("methane_GM"),
            volume_mixing_ratio_o2 = input_vmr("oxygen_GM"),
            volume_mixing_ratio_n2 = input_vmr("nitrogen_GM"),
            volume_mixing_ratio_ccl4 = input_vmr("carbon_tetrachloride_GM"),
            volume_mixing_ratio_cfc11 = input_vmr("cfc11_GM"),
            volume_mixing_ratio_cfc12 = input_vmr("cfc12_GM"),
            volume_mixing_ratio_cfc22 = input_vmr("hcfc22_GM"),
            volume_mixing_ratio_hfc143a = input_vmr("hfc143a_GM"),
            volume_mixing_ratio_hfc125 = input_vmr("hfc125_GM"),
            volume_mixing_ratio_hfc23 = input_vmr("hfc23_GM"),
            volume_mixing_ratio_hfc32 = input_vmr("hfc32_GM"),
            volume_mixing_ratio_hfc134a = input_vmr("hfc134a_GM"),
            volume_mixing_ratio_cf4 = input_vmr("cf4_GM"),
            volume_mixing_ratio_no2 = 1e-8, # not available in input_data
            latitude,
        )
        if !(radiation_mode isa ClearSkyRadiation)
            error("rrtmgp_model_cache not yet implemented for $radiation_mode")
        end
    end

    if requires_z(interpolation) || requires_z(bottom_extrapolation)
        kwargs = (;
            kwargs...,
            center_z = field2array(Fields.coordinate_field(Y.c).z),
            face_z = field2array(Fields.coordinate_field(Y.f).z),
        )
    end

    if idealized_insolation
        # perpetual equinox with no diurnal cycle
        solar_zenith_angle = FT(π) / 3
        weighted_irradiance =
            @. 1360 * (1 + FT(1.2) / 4 * (1 - 3 * sind(latitude)^2)) /
               (4 * cos(solar_zenith_angle))
    else
        solar_zenith_angle = weighted_irradiance = NaN # initialized in tendency
    end

    # surface_emissivity and surface_albedo are provided for each of 100 sites,
    # which we average across
    rrtmgp_model = RRTMGPModel(
        params;
        FT = Float64,
        ncol = length(Spaces.all_nodes(axes(Spaces.level(Y.c, 1)))),
        domain_nlay = Spaces.nlevels(axes(Y.c)),
        radiation_mode,
        interpolation,
        bottom_extrapolation,
        add_isothermal_boundary_layer = true,
        center_pressure = NaN, # initialized in tendency
        center_temperature = NaN, # initialized in tendency
        surface_temperature = (@. 29 * exp(-(latitude / 26)^2 / 2) + 271),
        surface_emissivity = mean(input_data["surface_emissivity"]),
        direct_sw_surface_albedo = mean(input_data["surface_albedo"]),
        diffuse_sw_surface_albedo = mean(input_data["surface_albedo"]),
        solar_zenith_angle,
        weighted_irradiance,
        kwargs...,
    )
    close(input_data)
    return (;
        ᶜT = similar(Y.c, FT),
        ᶜvmr_h2o = similar(Y.c, FT),
        insolation_tuple = similar(Spaces.level(Y.c, 1), Tuple{FT, FT, FT}),
        zenith_angle = similar(Spaces.level(Y.c, 1), FT),
        weighted_irradiance = similar(Spaces.level(Y.c, 1), FT),
        ᶠradiation_flux = similar(Y.f, Geometry.WVector{FT}),
        idealized_insolation,
        rrtmgp_model,
    )
end
function rrtmgp_model_tendency!(Yₜ, Y, p, t)
    (; ᶠradiation_flux) = p
    ᶜdivᵥ = Operators.DivergenceF2C()
    if :ρθ in propertynames(Y.c)
        error("rrtmgp_model_tendency! not implemented for ρθ")
    elseif :ρe in propertynames(Y.c)
        @. Yₜ.c.ρe -= ᶜdivᵥ(ᶠradiation_flux)
    elseif :ρe_int in propertynames(Y.c)
        @. Yₜ.c.ρe_int -= ᶜdivᵥ(ᶠradiation_flux)
    end
end
function rrtmgp_model_callback!(integrator)
    Y = integrator.u
    p = integrator.p
    t = integrator.t

    (; ᶜK, ᶜΦ, ᶜts, ᶜp, params) = p
    (; ᶜT, ᶜvmr_h2o, insolation_tuple, zenith_angle, weighted_irradiance) = p
    (; ᶠradiation_flux, idealized_insolation, rrtmgp_model) = p

    if :ρθ in propertynames(Y.c)
        @. ᶜts = thermo_state_ρθ(Y.c.ρθ, Y.c, params)
    elseif :ρe in propertynames(Y.c)
        @. ᶜK = norm_sqr(C123(Y.c.uₕ) + C123(ᶜinterp(Y.f.w))) / 2
        @. ᶜts = thermo_state_ρe(Y.c.ρe, Y.c, ᶜK, ᶜΦ, params)
    elseif :ρe_int in propertynames(Y.c)
        @. ᶜts = thermo_state_ρe_int(Y.c.ρe_int, Y.c, params)
    end

    @. ᶜp = TD.air_pressure(params, ᶜts)
    @. ᶜT = TD.air_temperature(params, ᶜts)
    rrtmgp_model.center_pressure .= field2array(ᶜp)
    rrtmgp_model.center_temperature .= field2array(ᶜT)

    if !(rrtmgp_model.radiation_mode isa GrayRadiation)
        @. ᶜvmr_h2o =
            TD.vol_vapor_mixing_ratio(params, TD.PhasePartition(params, ᶜts))
        rrtmgp_model.center_volume_mixing_ratio_h2o .= field2array(ᶜvmr_h2o)
    end

    if !idealized_insolation
        date_time = DateTime(2022) + Second(round(Int, t)) # t secs into 2022
        max_zenith_angle = FT(π) / 2 - eps(FT)
        irradiance = FT(Planet.tot_solar_irrad(params))
        au = FT(astro_unit())

        bottom_coords = Fields.coordinate_field(Spaces.level(Y.c, 1))
        @. insolation_tuple = instantaneous_zenith_angle(
            date_time,
            bottom_coords.long,
            bottom_coords.lat,
            params,
        ) # each tuple is (zenith angle, azimuthal angle, earth-sun distance)
        @. zenith_angle = min(first(insolation_tuple), max_zenith_angle)
        @. weighted_irradiance = irradiance * (au / last(insolation_tuple))^2
        rrtmgp_model.solar_zenith_angle .= field2array(zenith_angle)
        rrtmgp_model.weighted_irradiance .= field2array(weighted_irradiance)
    end

    update_fluxes!(rrtmgp_model)
    field2array(ᶠradiation_flux) .= rrtmgp_model.face_flux
end
