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

## Installation
Go to [julialang.org](https://juliailang.org/downloads/). 

*On Linux:*

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

*On Windows:*

```
TODO
```

*On MacOS:*

```
TODO
```

*For developers only:*

To setup `Revise.jl`, follow the instructions
described [here](https://timholy.github.io/Revise.jl/stable/) and make sure
to modify your `~/.julia/config/startup.jl` to include `using Revise`.

## Usage

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
@quickactivate :CSSim
```
which auto-activate the project and enable local path handling from DrWatson.

The above syntax is based on [making your DrWatson project a usable module](https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1),
which was done for the purposes of using `Revise.jl`. See the [discussion](https://discourse.julialang.org/t/best-debug-workflow-for-dr-watson/97234/5)
here for more information on this. 
