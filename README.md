# CSSim

Reproducing [Networks of forest ecosystems: Mathematical modeling of their biotic pump mechanism and resilience to certain patch deforestation](https://www.sciencedirect.com/science/article/pii/S1476945X20300386) as part of UvA's Complex System Simulation
2023 course.

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

Functionally, what this means is that in the `scripts/` and `notebooks/`
folder, you can import modules directly from the `src/` folder. For example,

```python
# In `notebooks/example_import.ipynb`
from src.penalty import alpha
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

*On Windows (unstable):*

```
TODO
```

*On MacOS (semi-stable):*

On macOS, a julia-1.9.1-mac64.dmg file is provided, which contains Julia-1.9.app. Installation is the same as any other Mac software: drag the Julia-1.9.app to Applications Folder's Shortcut. The Julia download runs on macOS 10.9 Mavericks and later releases. You can build from source for macOS 10.6 Snow Leopard (possibly earlier versions as well) and 32-bit but neither are fully supported.

You can launch Julia by opening the Julia app like any other application.
Optional: Add Julia to PATH

If you want to launch Julia from the command line, first open a new terminal window, then run the following snippet from your shell (e.g., using the Terminal app, not inside the Julia prompt).

```
sudo mkdir -p /usr/local/bin
sudo rm -f /usr/local/bin/julia
sudo ln -s /Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia
```

This code creates a symlink to a Julia version (here 1.9) of your choosing. To launch Julia, simply type julia inside your shell and press return.

## Usage

If Julia is already installed, change to the directory of this project
and do the following in a command prompt:

```
julia> using Pkg
julia> Pkg.install("DrWatson") # only if not already installed
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```

You may notice that most scripts start with the commands:

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

# References

[Cantin2020](https://www.sciencedirect.com/science/article/pii/S1476945X20300386)

[Antonovsky1990](https://www.sciencedirect.com/science/article/abs/pii/004058099090043U?via%3Dihub)

[Tutorial: Scipy odeint/solve_vip](https://danielmuellerkomorowska.com/2021/02/16/differential-equations-with-scipy-odeint-or-solve_ivp/)
