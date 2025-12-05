[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# On the approximation of separable non-convex optimization programs to an arbitrary numerical tolerance

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported in the paper [On the approximation of separable non-convex optimization programs to an arbitrary numerical tolerance](https://doi.org/10.1287/ijoc.2022.0221) by Claudio Contardo and Sandra U. Ngueveu.

The upstream repository (including a potentially more up-to-date version of this code) can be found at https://redmine.laas.fr/laas/users/ngueveu/r/ipwlb.git. This snapshot has been built from [this commit](https://redmine.laas.fr/laas/users/ngueveu/r/ipwlb.git/commit/f483ce379b24d49e842db33cc32df46cddca3ae5).

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2022.0221

https://doi.org/10.1287/ijoc.2022.0221.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{ContardoNgueveuIPWLB2025,
  author =        {Claudio Contardo and Sandra U. Ngueveu},
  publisher =     {INFORMS Journal on Computing},
  title =         {On the approximation of separable non-convex optimization programs to an arbitrary numerical tolerance},
  year =          {2025},
  doi =           {10.1287/ijoc.2022.0221.cd},
  url = {https://github.com/INFORMSJoC/2022.0221},
  note =          {Available for download at https://github.com/INFORMSJoC/2022.0221},
}
```

## Description

This repository provides the source code, data files and detailed results as reported in the article.

The main folders are 'data', 'results', and 'src'.
- '[instances](instances)': This folder includes instances of the multiple MINLPs benchmark problems.
- '[data](data)': This folder contains toy instances used by the method for precompilation.
- '[src](src)': The source code.
- '[results](results)': Detailed results, please see the README.md included in that folder for more information.


## Installation

To install this Julia package, simply execute from a Julia Pkg REPL (by pressing `]` within a regular Julia REPL) the following:
```julia
(@v1.10) pkg> add https://github.com/INFORMSJoC/2022.0221
```

This package depends on `LinA.jl`, a greedy method for constructing piecewise linear approximations to univariate functions of minimum size. We make use of the code available at https://github.com/claud10cv/LinA.jl. The official repository and a more up-to-date implementation of this method can be found at https://github.com/LICO-labs/PiecewiseLinApprox.jl. LinA.jl and PiecewiseLinApprox.jl implement the method described in [this paper](https://link.springer.com/article/10.1007/s12532-024-00274-8).

# Parameters available

The following struct is passed to the optimizer and can be used to control its behavior

```julia
mutable struct Parameters
    with_solver::DataType # Optimization solver, default = GLPK.Optimizer,
    with_solver_parameters::Any # Parameters passed to the solver, default = empty_parameters,
    with_lina_solver::DataType # Type of LinA method, either LinA.HeuristicLin or LinA.ExactLin
    with_cost_type::Symbol # type of cost function
    ptype::Symbol # problem type
    eps0::Float64 # INitial tolerance, default = 1e-1,
    eps::Float64 # final desired relative tolerance
    tilim::Float64 # Time limit, default = 3600,
    printon::Bool # Verbosity level, default = false,
    stepdiv::Float64 # Divisor of the tolrance when upadting a linear section, default =2.0,
    dynstep::Bool # Not relevant, not used, default = false,
    activeUpdateLargePwlf::Bool # Not relevant, not used, default =false,
    warmstartMILP::Bool # Solution of previous MILP is passed as warmstart the next call, default =false,
    miptype::Symbol # Not used, default = :grb,
    with_fixed_charge::Bool # Whether a fixed charge term is necessary at the origin, default = false
end
```
# Basic usage

To solve a nonlinear transportation problem, the following script reads the instance name and the type of cost function from the command line and solves the associated problem using our solver
```julia
using IterativePWLB
using LinA
using Gurobi

INSTANCE = ARGS[1]
CTYPE = Symbol(ARGS[2])

with_fixed_charge = (CTYPE in [:exp, :divexp, :pardivexp, :cubic, :sincos])

params = IterativePWLB.Parameters(Gurobi.Optimizer,
                                IterativePWLB.gurobi_parameters,
                                LinA.ExactLin,
                                CTYPE,
                                :transp,
                                1e+1,
                                1e-4,
                                3600.0,
                                false,
                                2.0,
                                false,
                                false,
                                true,
                                :grb,
                                with_fixed_charge)

res = solve_transp_from_file(INSTANCE,
                             params)

println("WRITING LOG...")
println.(res)
```
