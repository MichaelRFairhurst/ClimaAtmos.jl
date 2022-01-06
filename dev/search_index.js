var documenterSearchIndex = {"docs":
[{"location":"running_instructions/#Running-instructions","page":"Running instructions","title":"Running instructions","text":"","category":"section"},{"location":"running_instructions/","page":"Running instructions","title":"Running instructions","text":"We store conventional benchmark simulation examples in the test folder, where they have access to standardized initial conditions. These benchmarks are used and reused in several tests ranging from unit, over regression, to complex validation tests. This guarantees that the same code piece can be efficiently reused.","category":"page"},{"location":"running_instructions/","page":"Running instructions","title":"Running instructions","text":"Run all the test cases with:","category":"page"},{"location":"running_instructions/","page":"Running instructions","title":"Running instructions","text":"$ julia --project test/runtests.jl","category":"page"},{"location":"contributor_guide/#Contributors-Guide","page":"Contributor Guide","title":"Contributors Guide","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Thank you for considering contributions to ClimaAtmos! We hope this guide helps you make a contribution.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Feel free to ask us questions and chat with us at any time about any topic at all by:","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Opening a GitHub issue","category":"page"},{"location":"contributor_guide/#Creating-issues","page":"Contributor Guide","title":"Creating issues","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"The simplest way to contribute to ClimaAtmos is to create or comment on issues.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"The most useful bug reports:","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Provide an explicit code snippet –- not just a link –- that reproduces the bug in the latest tagged version of ClimaAtmos. This is sometimes called the \"minimal working example\". Reducing bug-producing code to a minimal example can dramatically decrease the time it takes to resolve an issue.\nPaste the entire error received when running the code snippet, even if it's unbelievably long.\nUse triple backticks (e.g., ```some_code; and_some_more_code;```) to enclose code snippets, and other markdown formatting syntax to make your issue easy and quick to read.\nReport the ClimaAtmos version, Julia version, machine (especially if using a GPU) and any other possibly useful details of the computational environment in which the bug was created.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Discussions are recommended for asking questions about (for example) the user interface, implementation details, science, and life in general.","category":"page"},{"location":"contributor_guide/#But-I-want-to-*code*!","page":"Contributor Guide","title":"But I want to code!","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"New users help write ClimaAtmos code and documentation by forking the ClimaAtmos repository, using git to edit code and docs, and then creating a pull request. Pull requests are reviewed by ClimaAtmos collaborators.\nA pull request can be merged once it is reviewed and approved by collaborators. If the pull request author has write access, they have the reponsibility of merging their pull request. Otherwise, ClimaAtmos.jl collabators will execute the merge with permission from the pull request author.\nNote: for small or minor changes (such as fixing a typo in documentation), the GitHub editor is super useful for forking and opening a pull request with a single click.\nWrite your code with love and care. In particular, conform to existing ClimaAtmos style and formatting conventions. For example, we love verbose and explicit variable names, use TitleCase for types, snake_case for objects, and always.put.spaces.after.commas. For formatting decisions we loosely follow the YASGuide. It's worth few extra minutes of our time to leave future generations with well-written, readable code.","category":"page"},{"location":"contributor_guide/#General-coding-guidelines","page":"Contributor Guide","title":"General coding guidelines","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Keep the number of members of Julia structs small if possible (less than 8 members).\nCode should reflect \"human intuition\" if possible. This mean abstraction should reflect how humans reason about the problem under consideration.\nCode with small blast radius. If your code needs to be modified or extendended, the resulting required changes should be as small and as localized as possible.\nWhen you write code, write it with testing and debugging in mind.\nIdeally, the lowest level structs have no defaults for their member fields. Nobody can remember all the defaults, so it is better to introduce them at the high-level API only.\nMake sure that module imports are specific so that it is easy to trace back where functions that are used inside a module are coming from.\nConsider naming abstract Julia types \"AbstractMyType\" in order to avoid confusion for the reader of your code.\nComments in your code should explain why the code exists and clarify if necessary, not just restate the line of code in words.\nBe mindful of namespace issues when writing functional code, especially when writing function code that represents mathematical or physical concepts.\nCondider using keywords in your structs to allow readers to more effectively reason about your code.","category":"page"},{"location":"contributor_guide/#What-is-a-\"collaborator\"-and-how-can-I-become-one?","page":"Contributor Guide","title":"What is a \"collaborator\" and how can I become one?","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Collaborators have permissions to review pull requests and status allows a contributor to review pull requests in addition to opening them. Collaborators can also create branches in the main ClimaAtmos repository.\nWe ask that new contributors try their hand at forking ClimaAtmos, and opening and merging a pull request before requesting collaborator status.","category":"page"},{"location":"contributor_guide/#What's-a-good-way-to-start-developing-ClimaAtmos?","page":"Contributor Guide","title":"What's a good way to start developing ClimaAtmos?","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Tackle an existing issue. We keep a list of good first issues that are self-contained and suitable for a newcomer to try and work on.\nTry to run ClimaAtmos and play around with it to simulate your favorite fluids and atmosphere physics. If you run into any problems or find it difficult to use or understand, please open an issue!\nWrite up an example or tutorial on how to do something useful with ClimaAtmos, like how to set up a new physical configuration.\nImprove documentation or comments if you found something hard to use.\nImplement a new feature if you need it to use ClimaAtmos.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"If you're interested in working on something, let us know by commenting on existing issues or  by opening a new issue. This is to make sure no one else is working on the same issue and so  we can help and guide you in case there is anything you need to know beforehand.","category":"page"},{"location":"contributor_guide/#Ground-Rules","page":"Contributor Guide","title":"Ground Rules","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Each pull request should consist of a logical collection of changes. You can include multiple bug fixes in a single pull request, but they should be related. For unrelated changes, please submit multiple pull requests.\nDo not commit changes to files that are irrelevant to your feature or bugfix (eg: .gitignore).\nBe willing to accept criticism and work on improving your code; we don't want to break other users' code, so care must be taken not to introduce bugs. We discuss pull requests and keep working on them until we believe we've done a good job.\nBe aware that the pull request review process is not immediate, and is generally proportional to the size of the pull request.","category":"page"},{"location":"contributor_guide/#Reporting-a-bug","page":"Contributor Guide","title":"Reporting a bug","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"The easiest way to get involved is to report issues you encounter when using ClimaAtmos or by requesting something you think is missing.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Head over to the issues page.\nSearch to see if your issue already exists or has even been solved previously.\nIf you indeed have a new issue or request, click the \"New Issue\" button.\nPlease be as specific as possible. Include the version of the code you were using, as well as what operating system you are running. The output of Julia's versioninfo() and ] status is helpful to include. Try your best to include a complete, \"minimal working example\" that reproduces the issue.","category":"page"},{"location":"contributor_guide/#Setting-up-your-development-environment","page":"Contributor Guide","title":"Setting up your development environment","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Install Julia on your system.\nInstall git on your system if it is not already there (install XCode command line tools on a Mac or git bash on Windows).\nLogin to your GitHub account and make a fork of the ClimaAtmos repository by clicking the \"Fork\" button.\nClone your fork of the ClimaAtmos repository (in terminal on Mac/Linux or git shell/ GUI on Windows) in the location you'd like to keep it.\ngit clone https://github.com/your-user-name/ClimaAtmos.jl.git\nNavigate to that folder in the terminal or in Anaconda Prompt if you're on Windows.\nConnect your repository to the upstream (main project).\ngit remote add ClimaAtmos https://github.com/CLiMA/ClimaAtmos.jl.git\nCreate the development environment by opening Julia via julia --project then typing in ] instantiate. This will install all the dependencies in the Project.toml file.\nYou can test to make sure ClimaAtmos works by typing in ] test. Doing so will run all the tests (and this can take a while).","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Your development environment is now ready!","category":"page"},{"location":"contributor_guide/#Pull-Requests","page":"Contributor Guide","title":"Pull Requests","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"We follow the ColPrac guide for collaborative practices. We ask that new contributors read that guide before submitting a pull request.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Changes and contributions should be made via GitHub pull requests against the main branch.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"When you're done making changes, commit the changes you made. Chris Beams has written a  guide on how to write good commit messages.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"When you think your changes are ready to be merged into the main repository, push to your fork and submit a pull request.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Working on your first Pull Request? You can learn how from this free video series How to Contribute to an Open Source Project on GitHub, Aaron Meurer's tutorial on the git workflow, or the guide “How to Contribute to Open Source\".","category":"page"},{"location":"contributor_guide/#Documentation","page":"Contributor Guide","title":"Documentation","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Generally, we follow the Julia conventions for documentation https://docs.julialang.org/en/v1/manual/documentation/.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Now that you've made your awesome contribution, it's time to tell the world how to use it. Writing documentation strings is really important to make sure others use your functionality properly. Didn't write new functions? That's fine, but be sure that the documentation for the code you touched is still in great shape. It is not uncommon to find some strange wording or clarification that you can take care of while you are here.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Here is an example of a docstring:","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"\"\"\"\n    Column([FT=Float64]; zlim, nelements)\n    \nCreates a column domain of type `FT`,\nwith extents zlim[1] < zlim[2] and `nelements` cells. \n\nExample:\nGenerate a Column{Float64} with extents (0,1) and 10 elements.\n```julia-repl\njulia> using ClimaAtmos.Domains\njulia> z_domain = Column(Float64, \n                            zlim = (0,1), \n                            nelements = 10)\n```\n\"\"\"\nfunction Column(FT::DataType = Float64; zlim, nelements)\n    @assert zlim[1] < zlim[2]\n    return Column{FT}(zlim, nelements)\nend","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"You can preview how the Documentation will look like after merging by building the documentation  locally. From the main directory of your local repository call","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"julia --project -e 'using Pkg; Pkg.instantiate()'\njulia --project=docs/ -e 'using Pkg; Pkg.instantiate()'\nJULIA_DEBUG=Documenter julia --project=docs/ docs/make.jl","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"and then open docs/build/index.html in your favorite browser. Providing the environment variable  JULIA_DEBUG=Documenter will provide with more information in the documentation build process and thus help figuring out a potential bug.","category":"page"},{"location":"contributor_guide/#Credits","page":"Contributor Guide","title":"Credits","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"This contributor's guide is heavily based on the excellent Oceananigans.jl contributor's guide which is heavily based on the excellent MetPy contributor's guide.","category":"page"},{"location":"function_index/#API","page":"Function Index","title":"API","text":"","category":"section"},{"location":"function_index/#Domains","page":"Function Index","title":"Domains","text":"","category":"section"},{"location":"function_index/","page":"Function Index","title":"Function Index","text":"Modules = [ClimaAtmos.Domains]","category":"page"},{"location":"function_index/#ClimaAtmos.Domains.AbstractDomain","page":"Function Index","title":"ClimaAtmos.Domains.AbstractDomain","text":"Supertype for all domains.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.Domains.Column-Union{Tuple{}, Tuple{Type{FT}}, Tuple{FT}} where FT","page":"Function Index","title":"ClimaAtmos.Domains.Column","text":"Column([FT = Float64]; zlim, nelements)\n\nConstruct a domain of type FT that represents a column along the z-axis with limits zlim (where zlim[1] < zlim[2]) and nelements elements. This domain is not periodic.\n\nExample\n\njulia> Column(zlim = (0, 1), nelements = 10)\nDomain set-up:\n\tz-column:\t[0.0, 1.0]\n\t# of elements:\t10\n\n\n\n\n\n","category":"method"},{"location":"function_index/#ClimaAtmos.Domains.HybridBox-Union{Tuple{}, Tuple{Type{FT}}, Tuple{FT}} where FT","page":"Function Index","title":"ClimaAtmos.Domains.HybridBox","text":"HybridBox([FT = Float64]; xlim, ylim, zlim, nelements, npolynomial, xperiodic = true, yperiodic = true)\n\nConstruct a domain of type FT that represents an xz-plane with limits xlim ylim and zlim (where xlim[1] < xlim[2],ylim[1] < ylim[2], and zlim[1] < zlim[2]), nelements elements of polynomial order npolynomial, x-axis periodicity xperiodic, and y-axis periodicity yperiodic. nelements must be a tuple with two values, with the first value corresponding to the x-axis, the second corresponding to the y-axis, and the third corresponding to the z-axis.  This domain is not periodic along the z-axis.\n\nExample\n\njulia> HybridBox(\n            xlim = (0, π),\n            ylim = (0, π),\n            zlim = (0, 1),\n            nelements = (5, 5, 10),\n            npolynomial = 5,\n            xperiodic = true,\n            yperiodic = true,\n            )\nDomain set-up:\n\txyz-box:\t[0.0, 3.1) × [0.0, 3.1) × [0.0, 1.0]\n\t# of elements:\t(5, 5, 10)\n\tpoly order:\t5\n\n\n\n\n\n","category":"method"},{"location":"function_index/#ClimaAtmos.Domains.HybridPlane-Union{Tuple{}, Tuple{Type{FT}}, Tuple{FT}} where FT","page":"Function Index","title":"ClimaAtmos.Domains.HybridPlane","text":"HybridPlane([FT = Float64]; xlim, zlim, nelements, npolynomial, xperiodic = true)\n\nConstruct a domain of type FT that represents an xz-plane with limits xlim and zlim (where xlim[1] < xlim[2] and zlim[1] < zlim[2]), nelements elements of polynomial order npolynomial, and x-axis periodicity xperiodic. nelements must be a tuple with two values, with the first value corresponding to the x-axis and the second corresponding to the z-axis. This domain is not periodic along the z-axis.\n\nExample\n\njulia> HybridPlane(\n            xlim = (0, π),\n            zlim = (0, 1),\n            nelements = (5, 10),\n            npolynomial = 5,\n            xperiodic = true,\n        )\nDomain set-up:\n\txz-plane:\t[0.0, 3.1) × [0.0, 1.0]\n\t# of elements:\t(5, 10)\n\tpoly order:\t5\n\n\n\n\n\n","category":"method"},{"location":"function_index/#ClimaAtmos.Domains.SphericalShell-Union{Tuple{}, Tuple{Type{FT}}, Tuple{FT}} where FT","page":"Function Index","title":"ClimaAtmos.Domains.SphericalShell","text":"SphericalShell([FT = Float64]; radius, height, nelements, npolynomial)\n\nConstruct a domain of type FT that represents a spherical shell with radius radius, height height, and nelements elements of polynomial order npolynomial.\n\nExample\n\njulia> SphericalShell(radius = 1, height = 1, nelements = (6, 10), npolynomial = 5)\nDomain set-up:\n\tsphere radius:\t1.0\n\tsphere height:\t1.0\n\t# of elements:\t(6, 10)\n\tpoly order:\t5\n\n\n\n\n\n","category":"method"},{"location":"function_index/#ClimaAtmos.Domains.make_function_space","page":"Function Index","title":"ClimaAtmos.Domains.make_function_space","text":"make_function_space(domain)\n\nConvert an AbstractDomain into a ClimaCore.Spaces.AbstactSpace.\n\n\n\n\n\n","category":"function"},{"location":"function_index/#Boundary-Conditions","page":"Function Index","title":"Boundary Conditions","text":"","category":"section"},{"location":"function_index/","page":"Function Index","title":"Function Index","text":"Modules = [ClimaAtmos.BoundaryConditions]","category":"page"},{"location":"function_index/#ClimaAtmos.BoundaryConditions.AbstractBoundaryCondition","page":"Function Index","title":"ClimaAtmos.BoundaryConditions.AbstractBoundaryCondition","text":"Supertype for all boundary conditions.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.BoundaryConditions.BulkFormulaCondition","page":"Function Index","title":"ClimaAtmos.BoundaryConditions.BulkFormulaCondition","text":"BulkFormulaCondition{C, T} <: AbstractBoundaryCondition\n\nComputes the boundary flux using the bulk formula and constant or Mohnin-Obukhov-based heat transfer coefficient Ch. Specific to potential temperature density.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.BoundaryConditions.CustomFluxCondition","page":"Function Index","title":"ClimaAtmos.BoundaryConditions.CustomFluxCondition","text":"CustomFluxCondition{F} <: AbstractBoundaryCondition\n\nComputes a user-defined boundary flux. The user is charged with making the custom flux function consistent wth the numerics of the model that invokes this boundary condition.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.BoundaryConditions.DragLawCondition","page":"Function Index","title":"ClimaAtmos.BoundaryConditions.DragLawCondition","text":"DragLawCondition{C} <: AbstractBoundaryCondition\n\nComputes the boundary flux using the bulk formula and constant or Mohnin-Obukhov-based drag coefficient Cd. Specific to momentum density.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.BoundaryConditions.NoFluxCondition","page":"Function Index","title":"ClimaAtmos.BoundaryConditions.NoFluxCondition","text":"NoFluxCondition <: AbstractBoundaryCondition\n\nComputes a fixed boundary flux of 0.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.BoundaryConditions.get_boundary_flux","page":"Function Index","title":"ClimaAtmos.BoundaryConditions.get_boundary_flux","text":"get_boundary_flux(model, bc, var, Ym, Ya)\n\nGet the flux of variable var across the boundary of model using boundary condition bc.\n\n\n\n\n\n","category":"function"},{"location":"function_index/#Models-Interface","page":"Function Index","title":"Models Interface","text":"","category":"section"},{"location":"function_index/","page":"Function Index","title":"Function Index","text":"Modules = [ClimaAtmos.Models]","category":"page"},{"location":"function_index/#ClimaAtmos.Models.AbstractModel","page":"Function Index","title":"ClimaAtmos.Models.AbstractModel","text":"Supertype for all models.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.Models.default_initial_conditions-Tuple{ClimaAtmos.Models.AbstractModel}","page":"Function Index","title":"ClimaAtmos.Models.default_initial_conditions","text":"default_initial_conditions(model)\n\nConstruct the initial conditions for model.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#ClimaAtmos.Models.make_ode_function-Tuple{ClimaAtmos.Models.AbstractModel}","page":"Function Index","title":"ClimaAtmos.Models.make_ode_function","text":"make_ode_function(model)\n\nConstruct the ordinary differential equations for model.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#Models","page":"Function Index","title":"Models","text":"","category":"section"},{"location":"function_index/","page":"Function Index","title":"Function Index","text":"Modules = (m = ClimaAtmos.Models; filter(x -> x isa Module && x != m, map(name -> getproperty(m, name), names(m; all = true)))) # all submodules of ClimaAtmos.Models","category":"page"},{"location":"function_index/#ClimaAtmos.Models.Nonhydrostatic2DModels.Nonhydrostatic2DModel","page":"Function Index","title":"ClimaAtmos.Models.Nonhydrostatic2DModels.Nonhydrostatic2DModel","text":"Nonhydrostatic2DModel <: AbstractModel\n\nA two-dimensional non-hydrostatic model, which is typically used for simulating the Euler equations. Required fields are domain, boundary_conditions, and parameters.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.Models.Nonhydrostatic3DModels.Nonhydrostatic3DModel","page":"Function Index","title":"ClimaAtmos.Models.Nonhydrostatic3DModels.Nonhydrostatic3DModel","text":"Nonhydrostatic3DModel <: AbstractModel\n\nA three-dimensional non-hydrostatic model, which is typically used for simulating the Euler equations. Required fields are domain, boundary_conditions, and parameters.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.Models.SingleColumnModels.SingleColumnModel","page":"Function Index","title":"ClimaAtmos.Models.SingleColumnModels.SingleColumnModel","text":"SingleColumnModel <: AbstractModel\n\nA single column model. Required fields are domain, boundary_conditions, and parameters.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#Callbacks","page":"Function Index","title":"Callbacks","text":"","category":"section"},{"location":"function_index/","page":"Function Index","title":"Function Index","text":"Modules = [ClimaAtmos.Callbacks]","category":"page"},{"location":"function_index/#ClimaAtmos.Callbacks.AbstractCallback","page":"Function Index","title":"ClimaAtmos.Callbacks.AbstractCallback","text":"Supertype for all callbacks.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.Callbacks.JLD2Output","page":"Function Index","title":"ClimaAtmos.Callbacks.JLD2Output","text":"JLD2Output{M, I} <: AbstractCallback\n\nSpecifies that a DiffEqCallbacks.PeriodicCallback should be constructed that extracts the model state from the integrator and stores it in a .jld2 file. \n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.Callbacks.generate_callback","page":"Function Index","title":"ClimaAtmos.Callbacks.generate_callback","text":"generate_callback(callback; kwargs...)\n\nConvert an AbstractCallback to a SciMLBase.DECallback.\n\n\n\n\n\n","category":"function"},{"location":"function_index/#Simulations","page":"Function Index","title":"Simulations","text":"","category":"section"},{"location":"function_index/","page":"Function Index","title":"Function Index","text":"Modules = [ClimaAtmos.Simulations]","category":"page"},{"location":"function_index/#ClimaAtmos.Simulations.AbstractRestart","page":"Function Index","title":"ClimaAtmos.Simulations.AbstractRestart","text":"Supertype for all restart modes.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.Simulations.NoRestart","page":"Function Index","title":"ClimaAtmos.Simulations.NoRestart","text":"NoRestart <: AbstractRestart\n\nSpecifies that a simulation should use its original initial conditions and end time.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.Simulations.Restart","page":"Function Index","title":"ClimaAtmos.Simulations.Restart","text":"Restart <: AbstractRestart\n\nSpecifies that a simulation should begin from the state recorded in a restart file and end at a specific time. The restart file must be a .jld2 file containing the simulation.integrator and simulation.model objects. Users must set! the simulation prior to restart.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.Simulations.Simulation-Tuple{ClimaAtmos.Models.AbstractModel, Any}","page":"Function Index","title":"ClimaAtmos.Simulations.Simulation","text":"Simulation(\n    model,\n    method;\n    Y_init = nothing,\n    dt,\n    tspan,\n    callbacks = nothing,\n    restart = NoRestart(),\n)\n\nConstruct a Simulation for a model with a time stepping method, initial conditions Y_init, and time step dt for a time interval of tspan. If Y_init is not provided, the model's default initial conditions are used. ```\n\n\n\n\n\n","category":"method"},{"location":"function_index/#ClimaAtmos.Simulations.run!-Tuple{ClimaAtmos.Simulations.Simulation, Vararg{Any, N} where N}","page":"Function Index","title":"ClimaAtmos.Simulations.run!","text":"run!(simulation, args...; kwargs...)\n\nRun the simulation to completion.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#ClimaAtmos.Simulations.set!","page":"Function Index","title":"ClimaAtmos.Simulations.set!","text":"set!(simulation, submodel_name = nothing; kwargs...)\n\nSet the simulation state to a new state, either through an array or a function.\n\n\n\n\n\n","category":"function"},{"location":"function_index/#ClimaAtmos.Simulations.step!-Tuple{ClimaAtmos.Simulations.Simulation, Vararg{Any, N} where N}","page":"Function Index","title":"ClimaAtmos.Simulations.step!","text":"step!(simulation, args...; kwargs...)\n\nAdvance the simulation by one time step.\n\n\n\n\n\n","category":"method"},{"location":"#ClimaAtmos.jl","page":"Home","title":"ClimaAtmos.jl","text":"","category":"section"},{"location":"installation_instructions/#Installation-instructions","page":"Installation instructions","title":"Installation instructions","text":"","category":"section"},{"location":"installation_instructions/","page":"Installation instructions","title":"Installation instructions","text":"Download the ClimaAtmos source with:","category":"page"},{"location":"installation_instructions/","page":"Installation instructions","title":"Installation instructions","text":"$ git clone https://github.com/CliMA/ClimaAtmos.jl.git","category":"page"},{"location":"installation_instructions/","page":"Installation instructions","title":"Installation instructions","text":"Now change into the ClimaAtmos.jl directory with ","category":"page"},{"location":"installation_instructions/","page":"Installation instructions","title":"Installation instructions","text":"$ cd ClimaAtmos.jl","category":"page"},{"location":"installation_instructions/","page":"Installation instructions","title":"Installation instructions","text":"To use ClimaAtmos, you need to instantiate all dependencies with:","category":"page"},{"location":"installation_instructions/","page":"Installation instructions","title":"Installation instructions","text":"$ julia --project\njulia>]\n(v1.6) pkg> instantiate","category":"page"}]
}
