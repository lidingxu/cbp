

# cbp: a solver for second order conic sub-modular bin pakcing problems

cbp is a specialized **branch-and-price** (BP) based solver for **second order conic sub-modular bin pakcing problems** (SOC SMBPs). SOC SMBPs are used to formulate chance constrained/robust bin packing problems. cbp is written in C++ based on [SCIP](https://www.scipopt.org/) and CPLEX.



## Benchmarks
We create a data format `cbp` to describe the network and interference. `ntrt` format has 5 groups of descriptors, and these 5 groups are listed one by one in the data file.

The first line (group) lists the numbers of vertices, edges,  interference cliques and demands.
```
#nv #ne #nc #nf
```
In the second group, each line describes the cost of the  i-th vertex, from the 1st to  the nv-th vertex.  Currently, the vertex cost is not used.

In the third group,  each line describes the tail vertex, end vertex, cost, capacity of the i-th edge, from the 1st to the ne-th vertex. 

In the fourth group, every two lines desribe the size of interference clique and the vertices in the i-th clique, from the 1st to the nc-th clique.

In the fifth group, each line describes the source vertex, the target vertex and the flow of the i-th demand, from the 1st to the nf-th demand. 

Instances are located in `data/benchmark1`,  `data/benchmark2` and  `data/benchmark1_rlp`.
Files in `data/benchmark1` and `data/benchmark1`are in `ntrt` format. `data/benchmark1_rlp` consists of `rlp` format files converted from `data/benchmark1`, `data/benchmark1` and `data/benchmark1_rlp` are used to test wUMCFC with CPLEX; `data/benchmark2` is used to test the effect of network coding.


## Installation
cbp has been developped and tested in *Linux* System. 

*Building* from sources requires:
- a SCIP installation.
- CMake.
- C++ compiler.
- CPLEX (version >= 20).

### Build
You may use the command line to configure and build wUMCFC by creating a `build` directory at the project directory and then building the configuration:
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
Otherwise your SCIP 's build type is *Debug*,  the following command builds but does not enable compilor's optimization : 
```
camke .. -DCMAKE_BUILD_TYPE=Debug
```
The complier generates a binary exectuable `wUMCFC`.

## Usage

You can run `cbp`  in the `build` directory via the following command:
```
./cbp
```
Set the #timelimit:
```
set  limits/time #timelimit
```

### cbp options
To modify the behavior of `cbp` (see [cppmain.cpp](https://github.com/mlubin/Pajarito.jl/blob/master/src/solver.jl) for default values), set parameters of `cbp` via command:
```
set  cbp #option
```
where '#option' can be set to the following option:

  * `is_misocp` Verbosity flag: 0 for quiet, 1 for basic solve info, 2 for iteration info, 3 for detailed timing and cuts and solution feasibility info
  * `is_bd_tight` Time limit for algorithm (in seconds)
  * `is_heur` Relative optimality gap termination condition
  * `is_stablize` Let MIP solver manage convergence ("branch and cut")
  * `is_parallelscplex`




Read the test #instance in 'ntrt' format:
```
read #instance
```
Optimize:
```
optimize
```



## References

If you find cbp useful in your work, we kindly request that you cite the following paper draft, which is recommended reading for advanced users:


    @unpublished{haddadvanier:hal-03242888,
      TITLE = {{A Branch-and-Price Algorithm for Energy Optimization in Multi-Hop Wireless Networks}},
      AUTHOR = {Haddad-Vanier, Sonia and Xu, Liding},
      URL = {https://hal.archives-ouvertes.fr/hal-03242888},
      NOTE = {working paper or preprint},
      TYPE = {Technical Report},
      INSTITUTION = {{OptimiX Team at LIX, Laboratoire d'Informatique de l'Ecole Polytechnique}},
      YEAR = {2021},
      MONTH = May,
      HAL_ID = {hal-03242888},
      HAL_VERSION = {v1},
    }
