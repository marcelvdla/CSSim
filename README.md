# CSSim

Reproducing [Networks of forest ecosystems: Mathematical modeling of their biotic pump mechanism and resilience to certain patch deforestation](https://www.sciencedirect.com/science/article/pii/S1476945X20300386) as part of UvA's Complex System Simulation 2023 course.

Members:

* Eva Lampret
* Marcel van de Lagemaat
* Jared Frazier
* Dominique Weltevreden

# Python Setup

If you have conda installed, do the following:

```shell
conda create --name CSSim pip
conda activate CSSim 
```

Then do the following regardless of whether you have conda:

```shell
pip install -r requirements.txt
pip install -e .
```


# Julia Setup

## Installation

Go to [julialang.org](https://julialang.org/downloads/).

*On Linux (this is the most stable setup, so please use this if possible):*

```
cd ~
wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.1-linux-x86_64.tar.gz
tar zxvf julia-1.9.1-linux-x86_64.tar.gz
```

Then edit your `~/.bashrc` file and add `export PATH="$PATH:~/julia-1.9.1/bin"`
to the end of the file. Note that you could replace the `~` in the path you
export to the `PATH` variable with `/home/<INSERT YOUR USERNAME>`. Finally,
`source ~/.bashrc` and then check to see that the Julia REPL work
by typing `julia`.

## Usage

If Julia is already installed, change to the directory of this project
and do the following in a command prompt:

```
julia> using Pkg
julia> Pkg.install("DrWatson") # only if not already installed
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```

You may notice that Julia scripts start with the commands:

```julia
using DrWatson
@quickactivate :CSSim
```

which auto-activate the project and enable local path handling from DrWatson.

The above syntax is based on [making your DrWatson project a usable module](https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1),
which was done for the purposes of using `Revise.jl`. See the [discussion](https://discourse.julialang.org/t/best-debug-workflow-for-dr-watson/97234/5)
here for more information on this.

*For developers only:*

To setup `Revise.jl`, follow the instructions described [here](https://timholy.github.io/Revise.jl/stable/) and make sure
to modify your `~/.julia/config/startup.jl` to include `using Revise`.

# Description of Codebase

Here we briefly describe relevant files to the project.

## `src`

`n_forest.jl & n_forest.py`: N-forest models (systems of ODEs).

`common.jl & common.py`: Helper functions for models.

`perturbation.py`: Functions for modeling perturbation of N-forest system.

## `scripts`

`deforested_position_sim.jl`: Sanity check code for deforestation simulation (corresponding in-part to table 4 of Cantin2020)

`plots.jl`: Used to compute `n` system phase portraits.

`sensitivity.py`: Conducts sensitivity analysis of parameters described in Cantin2020.

## `notebooks`

`jl_fig*.ipynb`: Reproduces figures the corresponding figures in Cantin2020.

`two_forest_model.ipynb`: Reproduces figure 8 from Cantin2020; base for scripts n_forest.py

`perturbed_networks.ipynb`: Reproduces figure 9 from Cantin2020.

`simulation_perturbed.ipynb`: Simulates multiple runs of the perturbed networks to gather data.

# References

[Cantin2020](https://www.sciencedirect.com/science/article/pii/S1476945X20300386)

[Antonovsky1990](https://www.sciencedirect.com/science/article/abs/pii/004058099090043U?via%3Dihub)

[Tutorial: Scipy odeint/solve_vip](https://danielmuellerkomorowska.com/2021/02/16/differential-equations-with-scipy-odeint-or-solve_ivp/)

[Julia Debugger.jl Tutorial](https://www.educative.io/answers/how-to-debug-script-in-julia)
