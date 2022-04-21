include("cli_options.jl")

const FT = parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
TEST_NAME = parsed_args["TEST_NAME"]


# Test-specific definitions (may be overwritten in each test case file)
# TODO: Allow some of these to be environment variables or command line arguments
params = nothing
horizontal_mesh = nothing # must be object of type AbstractMesh
npoly = 0
z_max = 0
z_elem = 0
t_end = if isnothing(parsed_args["t_end"])
    FT(60 * 60 * 24 * 10)
else
    parsed_args["t_end"]
end
dt = if isnothing(parsed_args["dt"])
    FT(400)
else
    parsed_args["dt"]
end
dt_save_to_sol = parsed_args["dt_save_to_sol"]
dt_save_to_disk = 0 # 0 means don't save to disk
ode_algorithm = nothing # must be object of type OrdinaryDiffEqAlgorithm
jacobian_flags = () # only required by implicit ODE algorithms
max_newton_iters = 10 # only required by ODE algorithms that use Newton's method
show_progress_bar = true
additional_callbacks = () # e.g., printing diagnostic information
additional_solver_kwargs = () # e.g., abstol and reltol
test_implicit_solver = false # makes solver extremely slow when set to `true`
additional_cache(Y, params, dt) = NamedTuple()
additional_tendency!(Yₜ, Y, p, t) = nothing
center_initial_condition(local_geometry, params) = NamedTuple()
face_initial_condition(local_geometry, params) = NamedTuple()
postprocessing(sol, output_dir) = nothing

################################################################################
is_distributed = haskey(ENV, "CLIMACORE_DISTRIBUTED")

using Logging
if is_distributed
    using ClimaComms
    if ENV["CLIMACORE_DISTRIBUTED"] == "MPI"
        using ClimaCommsMPI
        const comms_ctx = ClimaCommsMPI.MPICommsContext()
    else
        error("ENV[\"CLIMACORE_DISTRIBUTED\"] only supports the \"MPI\" option")
    end
    const pid, nprocs = ClimaComms.init(comms_ctx)
    logger_stream = ClimaComms.iamroot(comms_ctx) ? stderr : devnull
    prev_logger = global_logger(ConsoleLogger(logger_stream, Logging.Info))
    @info "Setting up distributed run on $nprocs \
        processor$(nprocs == 1 ? "" : "s")"
else
    using TerminalLoggers: TerminalLogger
    prev_logger = global_logger(TerminalLogger())
end
atexit() do
    global_logger(prev_logger)
end

using OrdinaryDiffEq
using DiffEqCallbacks
using JLD2

include("../implicit_solver_debugging_tools.jl")
include("../ordinary_diff_eq_bug_fixes.jl")
include("../common_spaces.jl")

test_dir, test_file_name = split(TEST_NAME, '/')
include(joinpath(test_dir, "$test_file_name.jl"))

import ClimaCore: enable_threading
enable_threading() = parsed_args["enable_threading"]

# TODO: When is_distributed is true, automatically compute the maximum number of
# bytes required to store an element from Y.c or Y.f (or, really, from any Field
# on which gather() or weighted_dss!() will get called). One option is to make a
# non-distributed space, extract the local_geometry type, and find the sizes of
# the output types of center_initial_condition() and face_initial_condition()
# for that local_geometry type. This is rather inefficient, though, so for now
# we will just hardcode the value of 4.
max_field_element_size = 4 # ρ = 1 byte, 𝔼 = 1 byte, uₕ = 2 bytes

if haskey(ENV, "RESTART_FILE")
    restart_file_name = ENV["RESTART_FILE"]
    if is_distributed
        restart_file_name =
            split(restart_file_name, ".jld2")[1] * "_pid$pid.jld2"
    end
    restart_data = jldopen(restart_file_name)
    t_start = restart_data["t"]
    Y = restart_data["Y"]
    close(restart_data)
    ᶜlocal_geometry = Fields.local_geometry_field(Y.c)
    ᶠlocal_geometry = Fields.local_geometry_field(Y.f)
else
    t_start = FT(0)
    if is_distributed
        h_space, comms_ctx =
            make_distributed_horizontal_space(horizontal_mesh, npoly, comms_ctx)
    else
        h_space = make_horizontal_space(horizontal_mesh, npoly)
        comms_ctx = nothing
    end
    center_space, face_space = make_hybrid_spaces(h_space, z_max, z_elem)
    ᶜlocal_geometry = Fields.local_geometry_field(center_space)
    ᶠlocal_geometry = Fields.local_geometry_field(face_space)
    Y = Fields.FieldVector(
        c = center_initial_condition.(ᶜlocal_geometry, params),
        f = face_initial_condition.(ᶠlocal_geometry, params),
    )
end
p = get_cache(Y, params, dt)

if ode_algorithm <: Union{
    OrdinaryDiffEq.OrdinaryDiffEqImplicitAlgorithm,
    OrdinaryDiffEq.OrdinaryDiffEqAdaptiveImplicitAlgorithm,
}
    use_transform = !(ode_algorithm in (Rosenbrock23, Rosenbrock32))
    W = SchurComplementW(Y, use_transform, jacobian_flags, test_implicit_solver)
    jac_kwargs = use_transform ? (; jac_prototype = W, Wfact_t = Wfact!) :
        (; jac_prototype = W, Wfact = Wfact!)

    alg_kwargs = (; linsolve = linsolve!)
    if ode_algorithm <: Union{
        OrdinaryDiffEq.OrdinaryDiffEqNewtonAlgorithm,
        OrdinaryDiffEq.OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    }
        alg_kwargs =
            (; alg_kwargs..., nlsolve = NLNewton(; max_iter = max_newton_iters))
    end
else
    jac_kwargs = alg_kwargs = ()
end

job_id = if isnothing(parsed_args["job_id"])
    job_id_from_parsed_args(s, parsed_args)
else
    parsed_args["job_id"]
end
output_dir = if isnothing(parsed_args["output_dir"])
    job_id
else
    parsed_args["output_dir"]
end
@info "Output directory: $output_dir"
mkpath(output_dir)

function make_save_to_disk_func(output_dir, test_file_name, is_distributed)
    function save_to_disk_func(integrator)
        day = floor(Int, integrator.t / (60 * 60 * 24))
        @info "Saving prognostic variables to JLD2 file on day $day"
        suffix = is_distributed ? "_pid$pid.jld2" : ".jld2"
        output_file = joinpath(output_dir, "$(test_file_name)_day$day$suffix")
        jldsave(output_file; t = integrator.t, Y = integrator.u)
        return nothing
    end
    return save_to_disk_func
end

save_to_disk_func =
    make_save_to_disk_func(output_dir, test_file_name, is_distributed)

dss_callback = FunctionCallingCallback(func_start = true) do Y, t, integrator
    p = integrator.p
    Spaces.weighted_dss!(Y.c, p.ghost_buffer.c)
    Spaces.weighted_dss!(Y.f, p.ghost_buffer.f)
end
if dt_save_to_disk == 0
    save_to_disk_callback = nothing
else
    save_to_disk_callback = PeriodicCallback(
        save_to_disk_func,
        dt_save_to_disk;
        initial_affect = true,
    )
end
callback =
    CallbackSet(dss_callback, save_to_disk_callback, additional_callbacks...)

problem = SplitODEProblem(
    ODEFunction(
        implicit_tendency!;
        jac_kwargs...,
        tgrad = (∂Y∂t, Y, p, t) -> (∂Y∂t .= FT(0)),
    ),
    remaining_tendency!,
    Y,
    (t_start, t_end),
    p,
)
integrator = OrdinaryDiffEq.init(
    problem,
    ode_algorithm(; alg_kwargs...);
    saveat = dt_save_to_sol == 0 ? [] : dt_save_to_sol,
    callback = callback,
    dt = dt,
    adaptive = false,
    progress = show_progress_bar,
    progress_steps = 1000,
    additional_solver_kwargs...,
)

if haskey(ENV, "CI_PERF_SKIP_RUN") # for performance analysis
    throw(:exit_profile)
end

@info "Running `$test_dir/$test_file_name` test case"
sol = @timev OrdinaryDiffEq.solve!(integrator)

if is_distributed # replace sol.u on the root processor with the global sol.u
    if ClimaComms.iamroot(comms_ctx)
        global_h_space = make_horizontal_space(horizontal_mesh, npoly)
        global_center_space, global_face_space =
            make_hybrid_spaces(global_h_space, z_max, z_elem)
        global_Y_c_type = Fields.Field{
            typeof(Fields.field_values(Y.c)),
            typeof(global_center_space),
        }
        global_Y_f_type = Fields.Field{
            typeof(Fields.field_values(Y.f)),
            typeof(global_face_space),
        }
        global_Y_type = Fields.FieldVector{
            FT,
            NamedTuple{(:c, :f), Tuple{global_Y_c_type, global_Y_f_type}},
        }
        global_sol_u = similar(sol.u, global_Y_type)
    end
    for i in 1:length(sol.u)
        global_Y_c =
            DataLayouts.gather(comms_ctx, Fields.field_values(sol.u[i].c))
        global_Y_f =
            DataLayouts.gather(comms_ctx, Fields.field_values(sol.u[i].f))
        if ClimaComms.iamroot(comms_ctx)
            global_sol_u[i] = Fields.FieldVector(
                c = Fields.Field(global_Y_c, global_center_space),
                f = Fields.Field(global_Y_f, global_face_space),
            )
        end
    end
    if ClimaComms.iamroot(comms_ctx)
        sol = DiffEqBase.sensitivity_solution(sol, global_sol_u, sol.t)
    end
end

import JSON
using Test
import OrderedCollections
if !is_distributed
    ENV["GKSwstype"] = "nul" # avoid displaying plots
    if TEST_NAME == "sphere/baroclinic_wave_rhoe" ||
       TEST_NAME == "sphere/baroclinic_wave_rhotheta"
        paperplots(sol, output_dir, p, FT(90), FT(180))
    else
        postprocessing(sol, output_dir)
    end
end

if !is_distributed || (is_distributed && ClimaComms.iamroot(Context))
    include(joinpath(@__DIR__, "..", "..", "post_processing", "mse_tables.jl"))

    if parsed_args["regression_test"]

        Y_last = sol.u[end]
        # This is helpful for starting up new tables
        @info "Job-specific MSE table format:"
        println("all_best_mse[$job_id] = OrderedCollections.OrderedDict()")
        for prop_chain in Fields.property_chains(Y_last)
            println("all_best_mse[$job_id][$prop_chain] = 0.0")
        end
        @info "Solution variables:"
        for prop_chain in Fields.property_chains(Y_last)
            println(prop_chain)
        end

        # Extract best mse for this job:
        best_mse = all_best_mse[job_id]

        include(joinpath(
            @__DIR__,
            "..",
            "..",
            "post_processing",
            "compute_mse.jl",
        ))

        ds_filename_computed = joinpath(output_dir, "prog_state.nc")

        function process_name(s::AbstractString)
            # "c_ρ", "c_ρe", "c_uₕ_1", "c_uₕ_2", "f_w_1"
            s = replace(s, "components_data_" => "")
            s = replace(s, "ₕ" => "_h")
            s = replace(s, "ρ" => "rho")
            return s
        end
        varname(pc::Tuple) = process_name(join(pc, "_"))

        export_nc(Y_last; nc_filename = ds_filename_computed, varname)
        computed_mse =
            regression_test(; job_id, best_mse, ds_filename_computed, varname)

        computed_mse_filename = joinpath(job_id, "computed_mse.json")

        open(computed_mse_filename, "w") do io
            JSON.print(io, computed_mse)
        end
        NCRegressionTests.test_mse(computed_mse, best_mse)
    end

end
