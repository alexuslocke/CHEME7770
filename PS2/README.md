### Simulation model for dual emission green fluorescent protein (deGFP) expression in a TX/TL 2.0 reaction.

This model describes the mRNA and protein concentration as a function of time as a system of
Ordinary Differential Equations (ODEs). The model equations are implemented in the [Julia](https://www.julialang.org) programming language. Required packaged are described in the ``Include.jl`` file.

#### Model parameters
Infrastructure parameters (elongation rates, RNAP/Ribosome concentrations, etc) are contained in the ``CellFree.json`` file. Kinetic parameter sets (N=20) are contained in the ``Ensemble-T4.dat`` file (more on how we got these later).

#### Solution of model equations
The model equations are solved using the [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) package, and the solutions are visualized using the
[PyPlot](https://github.com/JuliaPy/PyPlot.jl) package. The solution is generated using the ``sample_parameter_ensemble.jl`` script, while the solutions are plotted versus data (contained in the [data](https://github.com/varnerlab/CHEME-5440-7770-S20/tree/master/problem_sets/PS2/data) directory) using the ``visualize*.jl`` scripts.
