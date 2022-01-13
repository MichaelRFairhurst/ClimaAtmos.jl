"""
    Nonhydrostatic3DModel <: AbstractModel

A three-dimensional non-hydrostatic model, which is typically used for simulating
the Euler equations. Required fields are `domain`, `boundary_conditions`, and
`parameters`.
"""
Base.@kwdef struct Nonhydrostatic3DModel{D <: AbstractHybridDomain, BC, P} <:
                   AbstractModel
    domain::D
    boundary_conditions::BC
    parameters::P
    name::Symbol = :nhm
    varnames::Tuple = (:ρ, :ρuh, :ρw, :ρθ) # ρuh is the horizontal momentum
end

function Models.default_initial_conditions(model::Nonhydrostatic3DModel)
    space_c, space_f = make_function_space(model.domain)
    local_geometry_c = Fields.local_geometry_field(space_c)
    local_geometry_f = Fields.local_geometry_field(space_f)

    # functions that make zeros for this model
    zero_val = zero(Spaces.undertype(space_c))
    zero_scalar(lg) = zero_val # .
    zero_12vector(lg) = Geometry.Covariant12Vector(zero_val, zero_val) # ---->
    zero_3vector(lg) = Geometry.Covariant3Vector(zero_val)

    ρ = zero_scalar.(local_geometry_c)
    uh = zero_12vector.(local_geometry_c)
    w = zero_3vector.(local_geometry_f) # faces
    ρe_tot = zero_scalar.(local_geometry_c)

    return Fields.FieldVector(
        nhm = Fields.FieldVector(ρ = ρ, uh = uh, w = w, ρe_tot = ρe_tot),
    )
end

function Models.make_ode_function(model::Nonhydrostatic3DModel)
    FT = eltype(model.domain)

    # relevant parameters
    Ω::FT = CLIMAParameters.Planet.Omega(model.parameters)
    g::FT = CLIMAParameters.Planet.grav(model.parameters)
    κ₄::FT = model.hyperdiffusivity

    # gravitational potential
    Φ(z::FT) where {FT} = g * z

    # operators
    # spectral horizontal operators
    hdiv = Operators.Divergence()
    hwdiv = Operators.Divergence()
    hgrad = Operators.Gradient()
    hwgrad = Operators.Gradient()
    hcurl = Operators.Curl()
    hwcurl = Operators.Curl()

    # vertical FD operators with BC's
    # interpolators
    interp_c2f = Operators.InterpolateC2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    interp_f2c = Operators.InterpolateF2C()

    # gradients
    scalar_vgrad_c2f = Operators.GradientC2F(
        bottom = Operators.SetGradient(Geometry.Covariant3Vector(FT(0))),
        top = Operators.SetGradient(Geometry.Covariant3Vector(FT(0))),
    )

    # divergences
    vector_vdiv_f2c = Operators.DivergenceF2C(
        top = Operators.SetValue(Geometry.Contravariant3Vector(FT(0))),
        bottom = Operators.SetValue(Geometry.Contravariant3Vector(FT(0))),
    )

    # curls
    vector_vcurl_c2f = Operators.CurlC2F(
        bottom = Operators.SetCurl(Geometry.Contravariant12Vector(
            FT(0),
            FT(0),
        )),
        top = Operators.SetCurl(Geometry.Contravariant12Vector(FT(0), FT(0))),
    )

    # flux correction aka upwinding
    flux_correction_center = Operators.FluxCorrectionC2C(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )

    function rhs!(dY, Y, Ya, t)
        dYm = dY.nhm
        dρ = dYm.ρ # scalar on centers
        duh = dYm.uh # Covariant12Vector on centers
        dw = dYm.w # Covariant3Vector on faces
        dρe_tot = dYm.ρe_tot # scalar on centers
        Ym = Y.nhm
        ρ = Ym.ρ
        uh = Ym.uh
        w = Ym.w
        ρe_tot = Ym.ρe_tot

        # TODO!: Initialize all tendencies to zero for good practice!

        # calculate relevant thermodynamic quantities
        z = Fields.coordinate_field(axes(ρ)).z
        uvw = @. Geometry.Covariant123Vector(uh) +
           Geometry.Covariant123Vector(interp_f2c(w))
        e_int = @. ρe_tot / ρ - Φ(z) - norm(uvw)^2 / 2
        ts = Thermodynamics.PhaseDry.(model.parameters, e_int, ρ)
        p = Thermodynamics.air_pressure.(ts)

        # hyperdiffusion
        χe = @. dρe_tot = hwdiv(hgrad(ρe_tot / ρ))
        χuh = @. duh =
            hwgrad(hdiv(uh)) -
            Geometry.Covariant12Vector(hwcurl(Geometry.Covariant3Vector(hcurl(
                uh,
            ))),)
        Spaces.weighted_dss!(dρe_tot)
        Spaces.weighted_dss!(duh)
        @. dρe_tot = -κ₄ * hwdiv(ρ * hgrad(χe))
        @. duh =
            -κ₄ * (
                hwgrad(hdiv(χuh)) -
                Geometry.Covariant12Vector(hwcurl(Geometry.Covariant3Vector(hcurl(
                    χuh,
                ))),)
            )

        # density
        # the vector is split into horizontal and vertical components so that they can be
        # applied individually
        uvw = @. Geometry.Covariant123Vector(uh) +
           Geometry.Covariant123Vector(interp_f2c(w))
        @. dρ = -hdiv(ρ * uvw) # horizontal divergence
        @. dρ -= vector_vdiv_f2c(interp_c2f(ρ * uh)) # explicit vertical part
        @. dρ -= vector_vdiv_f2c(interp_c2f(ρ) * w) # TODO: implicit vertical part

        # horizontal momentum
        ω³ = @. hcurl(uh) # Contravariant3Vector
        ω¹² = @. hcurl(w) # Contravariant12Vector
        @. ω¹² += vector_vcurl_c2f(uh) # Contravariant12Vector
        u¹² = # TODO!: Will need to be changed with topography
            @. Geometry.Contravariant12Vector(Geometry.Covariant123Vector(interp_c2f(
                uh,
            )),) # Contravariant12Vector in 3D
        u³ = @. Geometry.Contravariant3Vector(Geometry.Covariant123Vector(w))  # TODO!: Will need to be changed with topography
        # coriolis
        if Ω != 0
            lat = Fields.coordinate_field(axes(ρ)).lat
            f = @. Geometry.Contravariant3Vector(Geometry.WVector(
                2 * Ω * sind(lat),
            ))  # TODO!: Will need to be changed with topography
        else
            y = Fields.coordinate_field(axes(ρ)).y
            f = @. Geometry.Contravariant3Vector(Geometry.WVector(Ω * y))
        end
        E = @. (norm(uvw)^2) / 2 + Φ(z)

        @. duh -= interp_f2c(ω¹² × u³)
        @. duh -=
            (f + ω³) ×
            Geometry.Contravariant12Vector(Geometry.Covariant123Vector(uh))
        @. duh -= hgrad(p) / ρ
        @. duh -= hgrad(E)

        # vertical momentum
        @. dw = ω¹² × u¹² # Covariant3Vector on faces
        @. dw -= scalar_vgrad_c2f(p) / interp_c2f(ρ)
        @. dw -= scalar_vgrad_c2f(E)

        # thermodynamic variable
        @. dρe_tot -= hdiv(uvw * (ρe_tot + p))
        @. dρe_tot -= vector_vdiv_f2c(w * interp_c2f(ρe_tot + p))
        @. dρe_tot -= vector_vdiv_f2c(interp_c2f(uh * (ρe_tot + p)))

        if model.flux_corr
            @. dρ += flux_correction_center(w, ρ)
            @. dρe_tot += flux_correction_center(w, ρe_tot)
        end

        # discrete stiffness summation for spectral operations
        Spaces.weighted_dss!(dρ)
        Spaces.weighted_dss!(duh)
        Spaces.weighted_dss!(dw)
        Spaces.weighted_dss!(dρe_tot)

        return dY
    end
end
