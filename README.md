

# cbp: branch-and-price algorithms for conic submodular binpacking

Binpacking is an important combinatorial problem in operations research. Conic submodular binpacking generalizes the classical binpacking problem to nonlinear setting, and it can formulate the chance-constrained binpacking, the distributionally robust binpacking etc. cbp implements several branch-and-price algorithms to solve the conic submodular binpacking problems.

cbp is written in C++ based on [SCIP](https://www.scipopt.org/) and CPLEX.


## Installation
cbp has been developped and tested in *Linux* System. 

To use cbp, it requires to install the following libraries:
- SCIP,
- CPLEX,
- CMake and C++ Compilors.

### Build
First, you need to modify the `CMakeLists.txt` file such that it matches with SCIP and CPLEX's include path and library path in your system.

You may use the command line to configure and build cbp by creating a `build` directory at the project directory and then building the configuration:
```
mkdir build
```
```
cd  build
```
If your SCIP 's build type is *Release*,  the following command builds and sets the compiler's optimization level to **O3** : 
```
camke .. -DCMAKE_BUILD_TYPE=Release
```
Otherwise your SCIP's build type is *Debug*,  the following command builds but does not enable compilor's optimization : 
```
camke .. -DCMAKE_BUILD_TYPE=Debug
```
The complier generates a binary exectuable `cbp`.



## Benchmarks
The `data` folder contains three benchmarks `CloudSmall`, `CloudMedium`, `CloudBig`. Each benchmark contains instances in `.cbp` file format.

The capacity constraint of conic submodular binpacking problems is the following
```
\sum_{i=1}^{n} a_i x_i + \sigma \sqrt{\sum_{i=1}^{n} b_i x_i} \le c
```

The `.cbp` file organizes the instance data in the following way:

In the first line, the statistics of the instance is given:
```
  type c n alpha
```
  * `type`: one of three data generation types in `d` (distributinally robust case), `g` (Gaussian case) and `h` (Hoeffding inequality case).
  * `c`:  a positive real of capacity.
  * `n`:  a positive integer of item size.
  * `alpha`: a real in [0,1] of risk level. 

The second line gives the `\sigma`.

Then, the third line gives the `a` in a list, and the fourth line gives the `b` in a list.


## Usage

You can run `cbp`  via the following command:
```
./cbp 
```
use the `read` command to read ab `.cbp` instance,

use the `set` command to set the parameters of the cbp algorithms:

  * `limits/time `: a positve real of the time limit in CPU seconds.
  * `cbp/is_misocp`:  a Boolean (TRUE/FALSE) value indicating wether to use the LP-BC algorithm / PWL-BC algorithm to solve the pricing problem (default: FALSE).
  * `cbp/is_parallelscplex`: a Boolean (TRUE/FALSE) value indicating wether to use CPLEX in parallelism mode (default: FALSE).
  * `cbp/is_bd_tight`: a Boolean (TRUE/FALSE) value indicating wether to use the bound tightening procedure for the PWL-BC algorithm (default: TRUE).
  * `cbp/is_adapt_points`: a Boolean (TRUE/FALSE) value indicating wether to use the two-stage adaptive breakpoints the PWL-BC algorithm (default: FALSE).
  * `cbp/is_heur`: a Boolean (FALSE/TRUE) value indicating wether to use the heuritic pricing algorithm  and the hybrid pricing strategy for the pricing problem (default: TRUE).
  * `heuristics/rmp/freq`: an integer value of the frequency of the column selection heuristic (-1: disable, default: 1).
  * `cbp/is_stablize`: the stablization mode of column generation, current in test and disabled.

To reproduce the computational results in the accompanied paper, run the script `testdir.py` to clear the `results` folder, and execute the following command in the terminal (in *Linux*)
```
/bin/bash run.sh
```
Note that you should change the [GNU parallel](https://www.gnu.org/software/parallel/) test option `gnuparalleltest` according to its availability in your system.


### Algorithm Option
We provide `.set` files in the `bp_setting` folder. These files can automatically set the options of `cbp` to the algorithms in the accompanied paper.

The algorithm option is summarized in the following table.


|                    | name in the paper| cbp/is_misocp |  cbp/is_bd_tight     |  cbp/is_adapt_points| cbp/is_heur   |heuristics/rmp/freq| cbp/is_parallelscplex| 
|--------------------|------------------|---------------|----------------------|---------------------|---------------|-------------------|-------|
| `bpsocp.set`       |      DW-BC    |  TRUE         |  FALSE               | FALSE               |  FALSE         | -1            | FALSE  |
| `bppl.set`         |      DW-PWL     |  FALSE        |  TRUE                | FALSE               |  FALSE          | -1            | FALSE  |
| `bphybrid.set`     |    DW-Hybrid    |  FALSE        |  TRUE                | FALSE              |  TRUE         | -1             | FALSE  |
| `bpheur.set`   |    DW-Hybrid*   |  FALSE        |  TRUE                | FALSE               |  TRUE            | 1             | FALSE  |
| `bphybridphadd.set`|    DW-Hybrid**  |  FALSE        |  TRUE                | TRUE                | TRUE             | 1          | FALSE  |





## References

If you find cbp useful in your work, we kindly request that you cite the following paper draft ([arXiv preprint](https://arxiv.org/abs/2204.00320)), which is recommended reading for advanced users:


    @misc{https://doi.org/10.48550/arxiv.2204.00320,
      doi = {10.48550/ARXIV.2204.00320},
      
      url = {https://arxiv.org/abs/2204.00320},
      
      author = {Xu, Liding and D'Ambrosio, Claudia and Vanier, Sonia Haddad and Traversi, Emiliano},
      
      keywords = {Optimization and Control (math.OC), FOS: Mathematics, FOS: Mathematics},
      
      title = {Branch and Price for Submodular Bin Packing},
      
      publisher = {arXiv},
      
      year = {2022},
      
      copyright = {Creative Commons Attribution Non Commercial No Derivatives 4.0 International}
    }



