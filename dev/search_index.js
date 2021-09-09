var documenterSearchIndex = {"docs":
[{"location":"running_instructions/#Running-instructions","page":"Running instructions","title":"Running instructions","text":"","category":"section"},{"location":"running_instructions/","page":"Running instructions","title":"Running instructions","text":"We store conventional benchmark simulation examples in the test folder, where they have access to standardized initial conditions. These benchmarks are used and reused in several tests ranging from unit, over regression, to complex validation tests. This guarantees that the same code piece can be efficiently reused.","category":"page"},{"location":"running_instructions/","page":"Running instructions","title":"Running instructions","text":"Run all the test cases with:","category":"page"},{"location":"running_instructions/","page":"Running instructions","title":"Running instructions","text":"$ julia --project test/runtests.jl","category":"page"},{"location":"contributor_guide/#Contributors-Guide","page":"Contributor Guide","title":"Contributors Guide","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Thank you for considering contributions to ClimaAtmos! We hope this guide helps you make a contribution.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Feel free to ask us questions and chat with us at any time about any topic at all by:","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Opening a GitHub issue","category":"page"},{"location":"contributor_guide/#Creating-issues","page":"Contributor Guide","title":"Creating issues","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"The simplest way to contribute to ClimaAtmos is to create or comment on issues.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"The most useful bug reports:","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Provide an explicit code snippet –- not just a link –- that reproduces the bug in the latest tagged version of ClimaAtmos. This is sometimes called the \"minimal working example\". Reducing bug-producing code to a minimal example can dramatically decrease the time it takes to resolve an issue.\nPaste the entire error received when running the code snippet, even if it's unbelievably long.\nUse triple backticks (e.g., ```some_code; and_some_more_code;```) to enclose code snippets, and other markdown formatting syntax to make your issue easy and quick to read.\nReport the ClimaAtmos version, Julia version, machine (especially if using a GPU) and any other possibly useful details of the computational environment in which the bug was created.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Discussions are recommended for asking questions about (for example) the user interface, implementation details, science, and life in general.","category":"page"},{"location":"contributor_guide/#But-I-want-to-*code*!","page":"Contributor Guide","title":"But I want to code!","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"New users help write ClimaAtmos code and documentation by forking the ClimaAtmos repository, using git to edit code and docs, and then creating a pull request. Pull requests are reviewed by ClimaAtmos collaborators.\nA pull request can be merged once it is reviewed and approved by collaborators. If the pull request author has write access, they have the reponsibility of merging their pull request. Otherwise, ClimaAtmos.jl collabators will execute the merge with permission from the pull request author.\nNote: for small or minor changes (such as fixing a typo in documentation), the GitHub editor is super useful for forking and opening a pull request with a single click.\nWrite your code with love and care. In particular, conform to existing ClimaAtmos style and formatting conventions. For example, we love verbose and explicit variable names, use TitleCase for types, snake_case for objects, and always.put.spaces.after.commas. For formatting decisions we loosely follow the YASGuide. It's worth few extra minutes of our time to leave future generations with well-written, readable code.","category":"page"},{"location":"contributor_guide/#General-coding-guidelines","page":"Contributor Guide","title":"General coding guidelines","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Keep the number of members of Julia structs small if possible (less than 8 members).\nCode should reflect \"human intuition\" if possible. This mean abstraction should reflect how humans reason about the problem under consideration.\nCode with small blast radius. If your code needs to be modified or extendended, the resulting required changes should be as small and as localized as possible.\nWhen you write code, write it with testing and debugging in mind.\nIdeally, the lowest level structs have no defaults for their member fields. Nobody can remember all the defaults, so it is better to introduce them at the high-level API only.\nMake sure that module imports are specific so that it is easy to trace back where functions that are used inside a module are coming from.\nConsider naming abstract Julia types \"AbstractMyType\" in order to avoid confusion for the reader of your code.\nComments in your code should explain why the code exists and clarify if necessary, not just restate the line of code in words.\nBe mindful of namespace issues when writing functional code, especially when writing function code that represents mathematical or physical concepts.\nCondider using keywords in your structs to allow readers to more effectively reason about your code.","category":"page"},{"location":"contributor_guide/#What-is-a-\"collaborator\"-and-how-can-I-become-one?","page":"Contributor Guide","title":"What is a \"collaborator\" and how can I become one?","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Collaborators have permissions to review pull requests and status allows a contributor to review pull requests in addition to opening them. Collaborators can also create branches in the main ClimaAtmos repository.\nWe ask that new contributors try their hand at forking ClimaAtmos, and opening and merging a pull request before requesting collaborator status.","category":"page"},{"location":"contributor_guide/#What's-a-good-way-to-start-developing-ClimaAtmos?","page":"Contributor Guide","title":"What's a good way to start developing ClimaAtmos?","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Tackle an existing issue. We keep a list of good first issues that are self-contained and suitable for a newcomer to try and work on.\nTry to run ClimaAtmos and play around with it to simulate your favorite fluids and atmosphere physics. If you run into any problems or find it difficult to use or understand, please open an issue!\nWrite up an example or tutorial on how to do something useful with ClimaAtmos, like how to set up a new physical configuration.\nImprove documentation or comments if you found something hard to use.\nImplement a new feature if you need it to use ClimaAtmos.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"If you're interested in working on something, let us know by commenting on existing issues or  by opening a new issue. This is to make sure no one else is working on the same issue and so  we can help and guide you in case there is anything you need to know beforehand.","category":"page"},{"location":"contributor_guide/#Ground-Rules","page":"Contributor Guide","title":"Ground Rules","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Each pull request should consist of a logical collection of changes. You can include multiple bug fixes in a single pull request, but they should be related. For unrelated changes, please submit multiple pull requests.\nDo not commit changes to files that are irrelevant to your feature or bugfix (eg: .gitignore).\nBe willing to accept criticism and work on improving your code; we don't want to break other users' code, so care must be taken not to introduce bugs. We discuss pull requests and keep working on them until we believe we've done a good job.\nBe aware that the pull request review process is not immediate, and is generally proportional to the size of the pull request.","category":"page"},{"location":"contributor_guide/#Reporting-a-bug","page":"Contributor Guide","title":"Reporting a bug","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"The easiest way to get involved is to report issues you encounter when using ClimaAtmos or by requesting something you think is missing.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Head over to the issues page.\nSearch to see if your issue already exists or has even been solved previously.\nIf you indeed have a new issue or request, click the \"New Issue\" button.\nPlease be as specific as possible. Include the version of the code you were using, as well as what operating system you are running. The output of Julia's versioninfo() and ] status is helpful to include. Try your best to include a complete, \"minimal working example\" that reproduces the issue.","category":"page"},{"location":"contributor_guide/#Setting-up-your-development-environment","page":"Contributor Guide","title":"Setting up your development environment","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Install Julia on your system.\nInstall git on your system if it is not already there (install XCode command line tools on a Mac or git bash on Windows).\nLogin to your GitHub account and make a fork of the ClimaAtmos repository by clicking the \"Fork\" button.\nClone your fork of the ClimaAtmos repository (in terminal on Mac/Linux or git shell/ GUI on Windows) in the location you'd like to keep it.\ngit clone https://github.com/your-user-name/ClimaAtmos.jl.git\nNavigate to that folder in the terminal or in Anaconda Prompt if you're on Windows.\nConnect your repository to the upstream (main project).\ngit remote add ClimaAtmos https://github.com/CLiMA/ClimaAtmos.jl.git\nCreate the development environment by opening Julia via julia --project then typing in ] instantiate. This will install all the dependencies in the Project.toml file.\nYou can test to make sure ClimaAtmos works by typing in ] test. Doing so will run all the tests (and this can take a while).","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Your development environment is now ready!","category":"page"},{"location":"contributor_guide/#Pull-Requests","page":"Contributor Guide","title":"Pull Requests","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"We follow the ColPrac guide for collaborative practices. We ask that new contributors read that guide before submitting a pull request.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Changes and contributions should be made via GitHub pull requests against the main branch.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"When you're done making changes, commit the changes you made. Chris Beams has written a  guide on how to write good commit messages.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"When you think your changes are ready to be merged into the main repository, push to your fork and submit a pull request.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Working on your first Pull Request? You can learn how from this free video series How to Contribute to an Open Source Project on GitHub, Aaron Meurer's tutorial on the git workflow, or the guide “How to Contribute to Open Source\".","category":"page"},{"location":"contributor_guide/#Documentation","page":"Contributor Guide","title":"Documentation","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Generally, we follow the Julia conventions for documentation https://docs.julialang.org/en/v1/manual/documentation/.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Now that you've made your awesome contribution, it's time to tell the world how to use it. Writing documentation strings is really important to make sure others use your functionality properly. Didn't write new functions? That's fine, but be sure that the documentation for the code you touched is still in great shape. It is not uncommon to find some strange wording or clarification that you can take care of while you are here.","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"Here is an example of a docstring:","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"\"\"\"\n    Column([FT=Float64]; zlim, nelements)\n    \nCreates a column domain of type `FT`,\nwith extents zlim[1] < zlim[2] and `nelements` cells. \n\nExample:\nGenerate a Column{Float64} with extents (0,1) and 10 elements.\n```julia-repl\njulia> using ClimaAtmos.Domains\njulia> z_domain = Column(Float64, \n                            zlim = (0,1), \n                            nelements = 10)\n```\n\"\"\"\nfunction Column(FT::DataType = Float64; zlim, nelements)\n    @assert zlim[1] < zlim[2]\n    return Column{FT}(zlim, nelements)\nend","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"You can preview how the Documentation will look like after merging by building the documentation  locally. From the main directory of your local repository call","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"julia --project -e 'using Pkg; Pkg.instantiate()'\njulia --project=docs/ -e 'using Pkg; Pkg.instantiate()'\nJULIA_DEBUG=Documenter julia --project=docs/ docs/make.jl","category":"page"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"and then open docs/build/index.html in your favorite browser. Providing the environment variable  JULIA_DEBUG=Documenter will provide with more information in the documentation build process and thus help figuring out a potential bug.","category":"page"},{"location":"contributor_guide/#Credits","page":"Contributor Guide","title":"Credits","text":"","category":"section"},{"location":"contributor_guide/","page":"Contributor Guide","title":"Contributor Guide","text":"This contributor's guide is heavily based on the excellent Oceananigans.jl contributor's guide which is heavily based on the excellent MetPy contributor's guide.","category":"page"},{"location":"function_index/#API","page":"Function Index","title":"API","text":"","category":"section"},{"location":"function_index/#Simulations","page":"Function Index","title":"Simulations","text":"","category":"section"},{"location":"function_index/","page":"Function Index","title":"Function Index","text":"Modules = [ClimaAtmos.Simulations]","category":"page"},{"location":"function_index/#ClimaAtmos.Simulations.AbstractSimulation","page":"Function Index","title":"ClimaAtmos.Simulations.AbstractSimulation","text":"AbstractSimulation\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.Simulations.Simulation","page":"Function Index","title":"ClimaAtmos.Simulations.Simulation","text":"struct Simulation <: AbstractSimulation\n\nA simulation wraps an abstract ClimaAtmos model containing  equation specifications and an instance of an integrator used for time integration of the discretized model PDE.\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.Simulations.Simulation-Tuple{ClimaAtmos.Models.AbstractModel, Any}","page":"Function Index","title":"ClimaAtmos.Simulations.Simulation","text":"Simulation(model::AbstractModel, method; Y_init = nothing, dt, tspan)\n\nConstruct a Simulation for a model with a time stepping method, initial conditions Y_init, time step Δt for tspan time interval.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#ClimaAtmos.Simulations.run!-Tuple{ClimaAtmos.Simulations.AbstractSimulation, Vararg{Any, N} where N}","page":"Function Index","title":"ClimaAtmos.Simulations.run!","text":"run!(simulation::AbstractSimulation, args...; kwargs...)\n\nRun a simulation to the end.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#ClimaAtmos.Simulations.set!","page":"Function Index","title":"ClimaAtmos.Simulations.set!","text":"set!(\n    simulation::AbstractSimulation,\n    submodel_name = nothing;\n    kwargs...,\n)\n\nSet the simulation state to a new state, either through  an array or a function.\n\n\n\n\n\n","category":"function"},{"location":"function_index/#SciMLBase.step!-Tuple{ClimaAtmos.Simulations.AbstractSimulation, Vararg{Any, N} where N}","page":"Function Index","title":"SciMLBase.step!","text":"step!(simulation::AbstractSimulation, args...; kwargs...)\n\nStep forward a simulation one time step.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#Models","page":"Function Index","title":"Models","text":"","category":"section"},{"location":"function_index/","page":"Function Index","title":"Function Index","text":"Modules = [ClimaAtmos.Models]\n","category":"page"},{"location":"function_index/#ClimaAtmos.Models.AbstractModel","page":"Function Index","title":"ClimaAtmos.Models.AbstractModel","text":"AbstractModel\n\n\n\n\n\n","category":"type"},{"location":"function_index/#ClimaAtmos.Models.default_initial_conditions-Tuple{ClimaAtmos.Models.AbstractModel}","page":"Function Index","title":"ClimaAtmos.Models.default_initial_conditions","text":"default_initial_conditions(model)\n\nConstruct the initial conditions for model.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#ClimaAtmos.Models.get_boundary_flux","page":"Function Index","title":"ClimaAtmos.Models.get_boundary_flux","text":"get_boundary_flux\n\nConstruct the boundary fluxes.\n\n\n\n\n\n","category":"function"},{"location":"function_index/#ClimaAtmos.Models.make_ode_function-Tuple{ClimaAtmos.Models.AbstractModel}","page":"Function Index","title":"ClimaAtmos.Models.make_ode_function","text":"make_ode_function(model)\n\nConstruct the ordinary differential equations for model.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#ShallowWaterModels","page":"Function Index","title":"ShallowWaterModels","text":"","category":"section"},{"location":"function_index/","page":"Function Index","title":"Function Index","text":"Modules = [ClimaAtmos.ShallowWaterModels]","category":"page"},{"location":"#ClimaAtmos.jl","page":"Home","title":"ClimaAtmos.jl","text":"","category":"section"},{"location":"installation_instructions/#Installation-instructions","page":"Installation instructions","title":"Installation instructions","text":"","category":"section"},{"location":"installation_instructions/","page":"Installation instructions","title":"Installation instructions","text":"Download the ClimaAtmos source with:","category":"page"},{"location":"installation_instructions/","page":"Installation instructions","title":"Installation instructions","text":"$ git clone https://github.com/CliMA/ClimaAtmos.jl.git","category":"page"},{"location":"installation_instructions/","page":"Installation instructions","title":"Installation instructions","text":"Now change into the ClimaAtmos.jl directory with ","category":"page"},{"location":"installation_instructions/","page":"Installation instructions","title":"Installation instructions","text":"$ cd ClimaAtmos.jl","category":"page"},{"location":"installation_instructions/","page":"Installation instructions","title":"Installation instructions","text":"To use ClimaAtmos, you need to add the ClimaCore package and instantiate all dependencies with:","category":"page"},{"location":"installation_instructions/","page":"Installation instructions","title":"Installation instructions","text":"$ julia --project\njulia>]\n(v1.6) pkg> add https://github.com/CliMA/ClimaCore.jl.git\n(v1.6) pkg> instantiate","category":"page"}]
}
