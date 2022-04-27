@inline function calculate_pressure(Y, Ya, _...)
    error("not implemented for this model configuration.")
end

@inline function calculate_pressure(
    Y,
    Ya,
    ::AbstractBaseModelStyle,
    ::PotentialTemperature,
    ::Dry,
    params,
    FT,
)
    ρ = Y.base.ρ
    ρθ = Y.thermodynamics.ρθ

    p = @. TD.air_pressure(params, TD.PhaseDry_ρθ(params, ρ, ρθ / ρ))

    return p
end

@inline function calculate_pressure(
    Y,
    Ya,
    ::AdvectiveForm,
    ::InternalEnergy,
    ::Dry,
    params,
    FT,
)
    ρ = Y.base.ρ
    ρe_int = Y.thermodynamics.ρe_int

    e_int = @. ρe_int / ρ
    p = TD.air_pressure.(Ref(params), TD.PhaseDry.(Ref(params), e_int, ρ))

    return p
end

@inline function calculate_pressure(
    Y,
    Ya,
    ::AdvectiveForm,
    ::TotalEnergy,
    ::Dry,
    params,
    FT,
)
    ρ = Y.base.ρ
    uh = Y.base.uh
    w = Y.base.w
    ρe_tot = Y.thermodynamics.ρe_tot

    interp_f2c = Operators.InterpolateF2C()

    z = Fields.coordinate_field(axes(ρ)).z
    uvw = @. Geometry.Covariant123Vector(uh) +
       Geometry.Covariant123Vector(interp_f2c(w))
    Φ = calculate_gravitational_potential(Y, Ya, params, FT)

    e_int = @. ρe_tot / ρ - Φ - norm(uvw)^2 / 2
    p = TD.air_pressure.(Ref(params), TD.PhaseDry.(Ref(params), e_int, ρ))

    return p
end

@inline function calculate_pressure(
    Y,
    Ya,
    ::AdvectiveForm,
    ::TotalEnergy,
    ::EquilibriumMoisture,
    params,
    FT,
)
    ρ = Y.base.ρ
    uh = Y.base.uh
    w = Y.base.w
    ρe_tot = Y.thermodynamics.ρe_tot
    ρq_tot = Y.moisture.ρq_tot

    interp_f2c = Operators.InterpolateF2C()

    z = Fields.coordinate_field(axes(ρ)).z
    uvw = @. Geometry.Covariant123Vector(uh) +
       Geometry.Covariant123Vector(interp_f2c(w))
    Φ = calculate_gravitational_potential(Y, Ya, params, FT)
    e_int = @. ρe_tot / ρ - Φ - norm(uvw)^2 / 2
    q_tot = @. ρq_tot / ρ

    # saturation adjustment
    p =
        TD.air_pressure.(
            Ref(params),
            TD.PhaseEquil_ρeq.(Ref(params), ρ, e_int, q_tot),
        )

    return p
end
