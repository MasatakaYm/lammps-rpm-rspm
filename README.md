# lammps-rpm-rspm

## Description

LAMMPS implementation of the replica-permutation method (RPM) and replica sub-permutation method (RSPM). These methods are improved alternatives to the replica-exchange method.

The below is a list of codes that are put in `src` directory:
1. `temper_permutation.h`, `temper_permutation.cpp`: RPM in the *NVT* ensemble.
2. `temper_subpermutation.h`, `temper_subpermutation.cpp`: RSPM in the  *NVT* ensemble.
3. `temper_permutation_npt.h`, `temper_permutation_npt.cpp`: RPM in the *NPT* ensemble.
4. `temper_subpermutation_npt.h`, `temper_subpermutation_npt.cpp`: RSPM in the *NPT* ensemble.


## Build

First of all, clone a repository of LAMMPS and change the branch as follows:
```sh
# Clone repository
$ git clone https://github.com/lammps/lammps.git
$ cd lammps
$ git checkout stable_3Mar2020
```

After changing the branch, put the codes of the RPM and RSPM (i.e., ` temper_*permutation*.h` and ` temper_*permutation*.cpp`) in the `lammps/src/USER-MISC/`. Then, install required packages and compile codes.

```sh
$ cd lammps/src

# Installing packages
$ make yes-kspace
$ make yes-molecule
$ make yes-rigid
$ make yes-replica
$ make yes-user-molfile
$ make yes-user-misc
$ make yes-opt
$ make yes-user-omp
$ make yes-user-intel

# Compile
$ make serial -j 10
$ make mpi -j 10
```

After compile, you will find two executable binary files; `lmp_serial` and `lmp_mpi`.

## LAMMPS version

Our source codes were tested with LAMMPS stable version **3Mar20**. We should remark that there were significant changes in the release of stable version **29Oct20** (In details, please see https://github.com/lammps/lammps/releases/tag/stable_29Oct2020). Because we have not been tested our codes with this version yet, we recommend that you use the stable 3Mar20 version of LAMMPS in order to perform RPM and RSPM.

## Using codes

### RPM and RSPM in the *NVT* ensemble

The RPM and RSPM in the *NVT* ensemble are provided as commands `temper/permutation` and `temper/subpermutation`, respectively. The syntax is as follows:

```txt
temper/permutation N M temp fix-ID seed1 seed2 index
```

```txt
temper/subpermutation N M temp fix-ID seed1 seed2 index
```

- `N`: the total number of simulation time steps.
- `M`: a time step interval of the temperature permutation (sub-permutation).
- `fix-ID`: a fix command ID of a thermostat.
- `seed1`: random number seed to decide on which set of subsets are used in the permutation attempts.
- `seed2`: random number seed for the Suwa-Todo algorithm.
- `index` (optional): specifies the temperature label of the replica.


Note that these commands have the same syntax as the `temper` command (https://docs.lammps.org/temper.html).


### RPM and RSPM in the *NPT* ensemble

The RPM and RSPM in the *NPT* ensemble are provided as commands `temper/permutation/npt` and `temper/subpermutation/npt`, respectively. The syntax is as follows:

```txt
temper/permutation/npt N M temp fix-ID seed1 seed2 pressure index
```

```txt
temper/subpermutation/npt N M temp fix-ID seed1 seed2 pressure index
```

- `pressure`: an pressure value of the replica.

Note that these commands have the same syntax as the `temper/npt` command (https://docs.lammps.org/temper_npt.html).

### Running simulations
Sample input files are put in `inputs` directory.

```sh
# An example simulation using temper/permutation command
$ cd inputs/poly-ethlene-chain/rpmd
mpirun -np 8 ${LMP} -partition 8x1 -in input-seed1.in > run-seed1.log

# An exmaple simulation using temper/subpermutation command
$ cd inputs/poly-ethlene-chain/rspmd
mpirun -np 8 ${LMP} -partition 8x1 -in input-seed1.in > run-seed1.log

# An example simulation using temper/permutation/npt
$ cd inpust/ala-dipeptide/rpmd
mpirun -np 80 ${LMP} -partition 10x8 -in md-rpmd-seed1.in > md-rpmd-seed1.log

# An example simulation using temper/subpermutation/npt
$ cd inputs/ala-dipeptide/rspmd
mpirun -np 80 ${LMP} -partition 10x8 -in md-rspmd-seed1.in > md-rspmd-seed1.lo
```

## References
1. [S. G. Itoh and H. Okumura, JCTC (2013)](https://pubs.acs.org/doi/10.1021/ct3007919)
2. [M. Yamauchi and H. Okumura, JCP (2017)](https://aip.scitation.org/doi/full/10.1063/1.4996431)
3. [M. Yamauchi and H. Okumura, JCC (2019)](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.26030)
4. [M. Yamauchi, G. La Penna, S. G. Itoh, and H. Okumura, Comput. Phys. Commun. (2022)](https://www.sciencedirect.com/science/article/pii/S0010465522000819)