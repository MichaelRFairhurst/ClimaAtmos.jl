function set_thermo_state!(ts, Yc, params, ::Val{:ρθ}, ::Val{:dry})
    (; ρ, ρθ) = Yc
    @. ts = TD.PhaseDry_ρθ(params, ρ, ρθ / ρ)
end
function set_thermo_state!(ts, Yc, params, ::Val{:ρθ}, ::Val{:equil})
    (; ρ, ρθ, ρq_tot) = Yc
    @. ts = TD.PhaseEquil_ρθq(params, ρ, ρθ / ρ, ρq_tot / ρ)
end
function set_thermo_state!(ts, Yc, params, ::Val{:ρθ}, ::Val{:nonequil})
    (; ρ, ρθ, ρq_tot, ρq_liq, ρq_ice) = Yc
    @. ts = TD.PhaseNonEquil_ρθq(
        params,
        ρ,
        ρθ / ρ,
        TD.PhasePartition(ρq_tot / ρ, ρq_liq / ρ, ρq_ice / ρ),
    )
end

set_thermo_state!(ts, Yc, params, args...) = set_thermo_state!(
    ts,
    Yc,
    params,
    Val(energy_name()),
    Val(moisture_mode()),
    args...,
)

function set_thermo_state!(
    ts,
    Yc,
    params,
    ::Val{:ρe},
    moisture_mode,
    K,
    Φ,
    ρe_int,
)
    (; ρ, ρe) = Yc
    @. ρe_int = ρe - ρ * (K + Φ)
    set_thermo_state!(ts, Yc, params, Val(:ρe_int), moisture_mode, ρe_int)
end

# Dispatcher:
set_thermo_state!(ts, Yc, params, ::Val{:ρe_int}, moisture_mode) =
    set_thermo_state!(ts, Yc, params, Val(:ρe_int), moisture_mode, Yc.ρe_int)

function set_thermo_state!(ts, Yc, params, ::Val{:ρe_int}, ::Val{:dry}, ρe_int)
    (; ρ) = Yc
    @. ts = TD.PhaseDry(params, ρe_int / ρ, ρ)
end
function set_thermo_state!(
    ts,
    Yc,
    params,
    ::Val{:ρe_int},
    ::Val{:equil},
    ρe_int,
)
    (; ρ, ρq_tot) = Yc
    @. ts = TD.PhaseEquil_ρeq(params, ρ, ρe_int / ρ, ρq_tot / ρ)
end
function set_thermo_state!(
    ts,
    Yc,
    params,
    ::Val{:ρe_int},
    ::Val{:nonequil},
    ρe_int,
)
    (; ρ, ρθ, ρq_tot, ρq_liq, ρq_ice) = Yc
    @. ts = TD.PhaseNonEquil(
        params,
        ρe_int / ρ,
        ρ,
        TD.PhasePartition(ρq_tot / ρ, ρq_liq / ρ, ρq_ice / ρ),
    )
end
