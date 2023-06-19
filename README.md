# CSSim

Reproducing [Networks of forest ecosystems: Mathematical modeling of their biotic pump mechanism and resilience to certain patch deforestation](https://www.sciencedirect.com/science/article/pii/S1476945X20300386) as part of UvA's Complex System Simulation 
2023 course.

Members:

* Eva Lampret
* Marcel van de Lagemaat
* Jared Frazier
* Dominique Weltevreden

# Python Setup (if applicable)

...

# Julia Setup 

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> CSSim

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "CSSim"
```
which auto-activate the project and enable local path handling from DrWatson.
