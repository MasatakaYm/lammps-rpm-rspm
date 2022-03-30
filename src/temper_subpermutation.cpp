/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Masataka Yamauchi (SOKENDAI, IMS)
------------------------------------------------------------------------- */

#include "temper_subpermutation.h"
#include <cmath>
#include <cstring>
#include "universe.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "integrate.h"
#include "modify.h"
#include "compute.h"
#include "force.h"
#include "fix.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

// #define TEMPER_SUBPERMUTATION_DEBUG 1

/* --------------------------q-------------------------------------------- */

TemperSubPermutation::TemperSubPermutation(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

TemperSubPermutation::~TemperSubPermutation()
{
  MPI_Comm_free(&roots);
  if (ranswap) delete ranswap;
  delete ranboltz;
  delete [] set_temp;
  delete [] temp2world;
  delete [] world2temp;
  delete [] world2root;

  delete [] nworlds_subsets;
  delete [] npermutations;

  for (int i = 0; i < max_npermutations; i++) {
    delete[] permutations[i];
  }
  delete[] permutations;
}

/* ----------------------------------------------------------------------
   perform tempering with inter-world swaps
------------------------------------------------------------------------- */

void TemperSubPermutation::command(int narg, char **arg)
{

  // Load arguments and check it.

  if (universe->nworlds == 1)
    error->all(FLERR,"Must have more than one processor partition to temper");
  if (domain->box_exist == 0)
    error->all(FLERR,"temper/subpermutation command before simulation box is defined");
  if (narg != 6 && narg != 7)
    error->universe_all(FLERR,"Illegal temper/subpermutation command");

  int nsteps = force->inumeric(FLERR,arg[0]);
  nevery = force->inumeric(FLERR,arg[1]);
  double temp = force->numeric(FLERR,arg[2]);

  // ignore temper command, if walltime limit was already reached

  if (timer->is_timeout()) return;

  for (whichfix = 0; whichfix < modify->nfix; whichfix++)
    if (strcmp(arg[3],modify->fix[whichfix]->id) == 0) break;
  if (whichfix == modify->nfix)
    error->universe_all(FLERR,"Tempering fix ID is not defined");

  seed_swap = force->inumeric(FLERR,arg[4]);
  seed_boltz = force->inumeric(FLERR,arg[5]);

  my_set_temp = universe->iworld;
  if (narg == 7) my_set_temp = force->inumeric(FLERR,arg[6]);
  if ((my_set_temp < 0) || (my_set_temp >= universe->nworlds))
    error->universe_one(FLERR,"Illegal temperature index");

  // swap frequency must evenly divide total # of timesteps

  if (nevery <= 0)
    error->universe_all(FLERR,"Invalid frequency in temper/subpermutation command");
  nswaps = nsteps/nevery;
  if (nswaps*nevery != nsteps)
    error->universe_all(FLERR,"Non integer # of swaps in temper/subpermutation command");

  // fix style must be appropriate for temperature control, i.e. it needs
  // to provide a working Fix::reset_target() and must not change the volume.

  if ((!utils::strmatch(modify->fix[whichfix]->style,"^nvt")) &&
      (!utils::strmatch(modify->fix[whichfix]->style,"^langevin")) &&
      (!utils::strmatch(modify->fix[whichfix]->style,"^gl[de]$")) &&
      (!utils::strmatch(modify->fix[whichfix]->style,"^rigid/nvt")) &&
      (!utils::strmatch(modify->fix[whichfix]->style,"^temp/")))
    error->universe_all(FLERR,"Tempering temperature fix is not supported");

  // setup for long tempering run

  update->whichflag = 1;
  timer->init_timeout();

  update->nsteps = nsteps;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + nsteps;
  if (update->laststep < 0)
    error->all(FLERR,"Too many timesteps");

  lmp->init();

  // local storage

  me_universe = universe->me;    // My place in universe
  MPI_Comm_rank(world,&me);      // Get my lank in world
  nworlds = universe->nworlds;   // # of worlds in universe
  iworld = universe->iworld;     // Which world I am in
  boltz = force->boltz;          // Boltzmann constant (eng/degree-K)

  // pe_compute = ptr to thermo_pe compute
  // notify compute it will be called at first swap

  int id = modify->find_compute("thermo_pe");
  if (id < 0) error->all(FLERR,"Tempering could not find thermo_pe compute");
  Compute *pe_compute = modify->compute[id];
  pe_compute->addstep(update->ntimestep + nevery);

  // create MPI communicator for root proc from each world

  int color;
  if (me == 0) color = 0;
  else color = 1;
  MPI_Comm_split(universe->uworld,color,0,&roots);

  // RNGs for swaps and Boltzmann test
  // warm up Boltzmann RNG
  // all procs have the same random seed so as to realize a transition to
  // the next state without MPI communications when parameter-permutation
  // trials are performed

  if (seed_swap) ranswap = new RanPark(lmp,seed_swap);
  else ranswap = NULL;
  ranboltz = new RanPark(lmp,seed_boltz);
  for (int i = 0; i < 100; i++) ranboltz->uniform();

  // world2root[i] = global proc that is root proc of world i

  world2root = new int[nworlds];
  if (me == 0)
    MPI_Allgather(&me_universe,1,MPI_INT,world2root,1,MPI_INT,roots);
  MPI_Bcast(world2root,nworlds,MPI_INT,0,world);

  // create static list of set temperatures
  // allgather tempering arg "temp" across root procs
  // bcast from each root to other procs in world

  set_temp = new double[nworlds];
  if (me == 0) MPI_Allgather(&temp,1,MPI_DOUBLE,set_temp,1,MPI_DOUBLE,roots);
  MPI_Bcast(set_temp,nworlds,MPI_DOUBLE,0,world);

  // create world2temp only on root procs from my_set_temp
  // create temp2world on root procs from world2temp,
  //   then bcast to all procs within world

  world2temp = new int[nworlds];
  temp2world = new int[nworlds];
  if (me == 0) {
    MPI_Allgather(&my_set_temp,1,MPI_INT,world2temp,1,MPI_INT,roots);
    for (int i = 0; i < nworlds; i++) temp2world[world2temp[i]] = i;
  }
  MPI_Bcast(temp2world,nworlds,MPI_INT,0,world);

  // if restarting tempering, reset temp target of Fix to current my_set_temp

  if (narg == 7) {
    double new_temp = set_temp[my_set_temp];
    modify->fix[whichfix]->reset_target(new_temp);
  }

  // setup tempering/subpermutation runs

  int which,swap,partner_set_temp;
  double pe,new_temp;

  const int limit_nworlds_subsets = 28;
  const int nsubsets = ((nworlds-1)/limit_nworlds_subsets) + 1;
  const int nworlds_average = nworlds/nsubsets;
  int nworlds_remainder = nworlds%nsubsets;

  nworlds_subsets = new int[nsubsets];
  npermutations = new int[nsubsets];

  for (int isub = 0; isub < nsubsets; ++isub)
      nworlds_subsets[isub] = nworlds_average;

  for (int isub = 0; isub < nsubsets; ++isub) {
    if (nworlds_remainder == 0) break;
    nworlds_subsets[isub]++;
    nworlds_remainder--;
  }
  max_nworlds_subsets = nworlds_subsets[0];

  max_npermutations = calc_nsubpermutations(max_nworlds_subsets);
  for (int isub = 0; isub < nsubsets; isub++) {
    npermutations[isub] = calc_nsubpermutations(nworlds_subsets[isub]);
  }

  permutations = new int *[max_npermutations];
  for (int iperm = 0; iperm < max_npermutations; iperm++) {
    permutations[iperm] = new int[max_nworlds_subsets];
  }

  init_subpermutations(max_nworlds_subsets, permutations);

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting up tempering with sub-permutation ...\n");

  if (me_universe == 0 && universe->uscreen) {
    // Subsets for which = 0
    fprintf(universe->uscreen, "Subsets:");
    for (int isub = 0; isub < nsubsets; isub++) {
      int iparam = 0;
      for (int jsub = 0; jsub < isub; jsub++) {
        iparam += nworlds_subsets[jsub];
      }

      int iworld = 0;
      int jparam = iparam;
      fprintf(universe->uscreen, "X%d=[", 2*isub);
      for (iworld = 0; iworld < nworlds_subsets[isub] - 1; iworld++) {
        jparam = iparam + iworld;
        jparam = jparam < nworlds ? jparam : jparam - nworlds;
        fprintf(universe->uscreen, "T%d,", jparam);
      }
      iworld = nworlds_subsets[isub] - 1;
      jparam = iparam + iworld;
      jparam = jparam < nworlds ? jparam : jparam - nworlds;
      fprintf(universe->uscreen, "T%d] ", jparam);
     }
     fprintf(universe->uscreen, "\n");
  }
  if (me_universe == 0 && universe->ulogfile) {
    // Subsets for which = 0
    fprintf(universe->ulogfile, "Subsets:");
    for (int isub = 0; isub < nsubsets; isub++) {
      int iparam = 0;
      for (int jsub = 0; jsub < isub; jsub++) {
        iparam += nworlds_subsets[jsub];
      }

      int iworld = 0;
      int jparam = iparam;
      fprintf(universe->ulogfile, "X%d=[", 2*isub);
      for (iworld = 0; iworld < nworlds_subsets[isub] - 1; iworld++) {
        jparam = iparam + iworld;
        jparam = jparam < nworlds ? jparam : jparam - nworlds;
        fprintf(universe->ulogfile, "T%d,", jparam);
      }
      iworld = nworlds_subsets[isub] - 1;
      jparam = iparam + iworld;
      jparam = jparam < nworlds ? jparam : jparam - nworlds;
      fprintf(universe->ulogfile, "T%d] ", jparam);
     }
     fprintf(universe->ulogfile, "\n");
  }

  if (nsubsets > 1) {
    if (me_universe == 0 && universe->uscreen) {
      // Subsets for which = 1
      fprintf(universe->uscreen, "Subsets:");
      for (int isub = 0; isub < nsubsets; isub++) {
        int iparam = nworlds_subsets[0] / 2;
        for (int jsub = 0; jsub < isub; ++jsub) {
          iparam += nworlds_subsets[jsub];
        }

        int iworld = 0;
        int jparam = iparam;
        fprintf(universe->uscreen, "X%d=[", 2*isub + 1);
        for (int iworld = 0; iworld < nworlds_subsets[isub]-1; iworld++) {
          jparam = iparam + iworld;
          jparam = jparam < nworlds ? jparam : jparam - nworlds;
          fprintf(universe->uscreen, "T%d,", jparam);
        }
        iworld = nworlds_subsets[isub] - 1;
        jparam = iparam + iworld;
        jparam = jparam < nworlds ? jparam : jparam - nworlds;
        fprintf(universe->uscreen, "T%d] ", jparam);
      }
      fprintf(universe->uscreen, "\n");
    }
  }
  if (nsubsets > 1) {
    if (me_universe == 0 && universe->ulogfile) {
      // Subsets for which = 1
      fprintf(universe->ulogfile, "Subsets:");
      for (int isub = 0; isub < nsubsets; isub++) {
        int iparam = nworlds_subsets[0] / 2;
        for (int jsub = 0; jsub < isub; ++jsub) {
          iparam += nworlds_subsets[jsub];
        }

        int iworld = 0;
        int jparam = iparam;
        fprintf(universe->ulogfile, "X%d=[", 2*isub + 1);
        for (int iworld = 0; iworld < nworlds_subsets[isub]-1; iworld++) {
          jparam = iparam + iworld;
          jparam = jparam < nworlds ? jparam : jparam - nworlds;
          fprintf(universe->ulogfile, "T%d,", jparam);
        }
        iworld = nworlds_subsets[isub] - 1;
        jparam = iparam + iworld;
        jparam = jparam < nworlds ? jparam : jparam - nworlds;
        fprintf(universe->ulogfile, "T%d] ", jparam);
      }
      fprintf(universe->ulogfile, "\n");
    }
  }

  update->integrate->setup(1);

  if (me_universe == 0) {
    if (universe->uscreen) {
      fprintf(universe->uscreen,"Step");
      for (int i = 0; i < nworlds; i++)
        fprintf(universe->uscreen," T%d",i);
      fprintf(universe->uscreen,"\n");
    }
    if (universe->ulogfile) {
      fprintf(universe->ulogfile,"Step");
      for (int i = 0; i < nworlds; i++)
        fprintf(universe->ulogfile," T%d",i);
      fprintf(universe->ulogfile,"\n");
    }
    print_status();
  }

  timer->init();
  timer->barrier_start();

  for (int iswap = 0; iswap < nswaps; iswap++) {

    // run for nevery timesteps

    timer->init_timeout();
    update->integrate->run(nevery);
    if (timer->is_timeout()) break;

    // compute PE
    // notify compute it will be called at next swap

    pe = pe_compute->compute_scalar();
    pe_compute->addstep(update->ntimestep + nevery);


    // which = which of 2 kinds of subsets
    // which = 0: even numberd subsets are used in permutation trials
    // which = 1: odd numberd subsets are used in permutation trials

    if (nsubsets == 1) which = 0;
    else if (!ranswap) which = iswap % 2;
    else if (ranswap->uniform() < 0.5) which = 0;
    else which = 1;

    // all procs (worlds) use the same random seed
    // results of permutation trials (i.e., istate_next value) must be the same
    // between all worlds; and therefore, it is not necessary to communicate
    // the istate_next value among worlds after determining it

    if (me == 0) {

      double pe_worlds[nworlds];

      // gather potential energy from the whole worlds.

      MPI_Allgather(&pe,1,MPI_DOUBLE,pe_worlds,1,MPI_DOUBLE,roots);

      for (int isub = 0; isub < nsubsets; isub++) {

        const int nparam_subset = nworlds_subsets[isub];
        int subtemp2temp[nparam_subset];

        // setup mapping function: subtemp2temp

        int iparam_sub_begin;
        if (which%2 == 0) iparam_sub_begin = 0;
        else iparam_sub_begin = nworlds_subsets[0] / 2;

        for (int jsub = 0; jsub < isub; jsub++) {
          iparam_sub_begin += nworlds_subsets[jsub];
        }

        for (int iparam_sub = 0; iparam_sub < nparam_subset; iparam_sub++) {
          int iparam = iparam_sub_begin + iparam_sub;
          if (iparam >= nworlds) iparam -= nworlds;
          subtemp2temp[iparam_sub] = iparam;
        }

        // calculate weight factors for all possible candidates:
        //   w(X) = prod_{m=1}^{M} exp[-beta_{m} U(x_{m}^{i(m)})]

        double w[npermutations[isub]];
        for (int i = 0; i < npermutations[isub]; i++) w[i] = 0.0e0;

        for (int icand = 0; icand < npermutations[isub]; icand++) {
          for (int jparam_sub = 0; jparam_sub < nparam_subset; jparam_sub++) {
            const int kparam_sub = permutations[icand][jparam_sub];
            const int jparam = subtemp2temp[jparam_sub];
            const int kparam = subtemp2temp[kparam_sub];
            const int kworld = temp2world[kparam];
            w[icand] -= pe_worlds[kworld] / (boltz*set_temp[jparam]);
          }
        }

        // search the maximum weight value

        int iw_max = 0;
        double w_max = w[0];
        for (int icand = 1; icand < npermutations[isub]; icand++)
        {
          if (w_max < w[icand]) {
            iw_max = icand;
            w_max  = w[icand];
          }
        }

        for (int icand = 0; icand < npermutations[isub]; icand++)
        {
          w[icand] = exp(w[icand] - w_max);
        }

        // calculate transition probability by Suwa-Todo algorithm

        const int icurrent_state = 0;
        double probs[npermutations[isub]];
        suwa_todo_algorithm(icurrent_state, npermutations[isub], w, probs);

        // select a next state following obtained probabilities

        int istate_next = 0;
        double sum_probs = 0.0;
        double rand_mcmc = ranboltz->uniform();

        for (int icand = 0; icand < npermutations[isub]; ++icand)
        {
          sum_probs += probs[icand];
          if (rand_mcmc <= sum_probs)
          {
            istate_next = icand;
            break;
          }
        }

        // update world2temp

        for (int iparam_sub = 0; iparam_sub < nparam_subset; iparam_sub++) {
          const int jparam_sub = permutations[istate_next][iparam_sub];
          const int jparam = subtemp2temp[jparam_sub];
          const int jworld = temp2world[jparam];
          world2temp[jworld] = subtemp2temp[iparam_sub];
        }


#ifdef TEMPER_SUBPERMUTATION_DEBUG
        if (me == 0 && me_universe == 0) {
          printf("Debug) Calculate weight factor.\n");
          printf("  Cand: ");
          for (int iparam_sub = 0; iparam_sub < nparam_subset; iparam_sub++) {
            const int iparam = subtemp2temp[iparam_sub];
            printf("[Rep%02d|T%02d] ", temp2world[iparam], iparam);
          }
          printf("\n");

          printf("--------"); // "Cand: "
          for (int iparam_sub = 0; iparam_sub < nparam_subset; iparam_sub++) {
            printf("-----------");
          }
          printf("--------------------------------------------------\n");

          double p_sum = 0.0;
          for (int icand = 0; icand < npermutations[isub]; icand++) {
            printf("%6d: ", icand);
            for (int jparam_sub = 0; jparam_sub < nparam_subset; jparam_sub++) {
              const int kparam_sub = permutations[icand][jparam_sub];
              const int jparam = subtemp2temp[jparam_sub];
              const int kparam = subtemp2temp[kparam_sub];
              const int kworld = temp2world[kparam];
              const int jworld = temp2world[jparam];
              printf("[Rep%02d|T%02d] ", jworld, kparam);
            }
            p_sum += probs[icand];
            printf("w = %6.4lf, p = %6.4lf, p_sum = %6.4lf\n"
                   , w[icand], probs[icand], p_sum);
          }
          printf("Select:\n");
          printf("%6d: ", istate_next);
          for (int jparam_sub = 0; jparam_sub < nparam_subset; jparam_sub++) {
              const int kparam_sub = permutations[istate_next][jparam_sub];
              const int jparam = subtemp2temp[jparam_sub];
              const int kparam = subtemp2temp[kparam_sub];
              const int kworld = temp2world[kparam];
              const int jworld = temp2world[jparam];
              printf("[Rep%02d|T%02d] ", jworld, kparam);
          }
          printf("w = %6.4lf, p = %6.4lf, rand  = %6.4lf\n"
                , w[istate_next], probs[istate_next], rand_mcmc);
        }
#endif
      }

      // partner_set_temp = which set temp my temp transits to

      partner_set_temp = world2temp[iworld];

      swap = 1;
      if (my_set_temp == partner_set_temp) swap = 0;
    }

    // bcast swap and partner_set_temp to other procs in my world

    MPI_Bcast(&partner_set_temp,1,MPI_INT,0,world);
    MPI_Bcast(&swap,1,MPI_INT,0,world);

    // rescale kinetic energy via velocities
    // if the temp transits to the other temperatures

    if (swap) scale_velocities(partner_set_temp,my_set_temp);

    // if my world swapped, all procs in world reset temp target of Fix

    if (swap) {
      new_temp = set_temp[partner_set_temp];
      modify->fix[whichfix]->reset_target(new_temp);
    }

    // update my_set_temp and temp2world on every proc
    // root procs update their value if swap took place
    // allgather across root procs
    // bcast within my world

    if (swap) my_set_temp = partner_set_temp;

    if (me == 0) {
      for (int iworld = 0; iworld < nworlds; iworld++) {
        temp2world[world2temp[iworld]] = iworld;
      }
    }
    MPI_Bcast(temp2world,nworlds,MPI_INT,0,world);

    // print out current swap status

    if (me_universe == 0) print_status();

  }

  timer->barrier_stop();

  update->integrate->cleanup();

  Finish finish(lmp);
  finish.end(1);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}

/* ----------------------------------------------------------------------
   scale kinetic energy via velocities a la Sugita
------------------------------------------------------------------------- */

void TemperSubPermutation::scale_velocities(int t_partner, int t_me)
{
  double sfactor = sqrt(set_temp[t_partner]/set_temp[t_me]);

  double **v = atom->v;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    v[i][0] = v[i][0]*sfactor;
    v[i][1] = v[i][1]*sfactor;
    v[i][2] = v[i][2]*sfactor;
  }
}

/* ----------------------------------------------------------------------
   proc 0 prints current tempering status
------------------------------------------------------------------------- */

void TemperSubPermutation::print_status()
{
  if (universe->uscreen) {
    fprintf(universe->uscreen,BIGINT_FORMAT,update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->uscreen," %d",world2temp[i]);
    fprintf(universe->uscreen,"\n");
  }
  if (universe->ulogfile) {
    fprintf(universe->ulogfile,BIGINT_FORMAT,update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->ulogfile," %d",world2temp[i]);
    fprintf(universe->ulogfile,"\n");
    fflush(universe->ulogfile);
  }
}

/* ----------------------------------------------------------------------
   construction of permutations
------------------------------------------------------------------------- */

int TemperSubPermutation::calc_nsubpermutations(const int nelems)
{
    int fibo_0 = 0;
    int fibo_1 = 1;
    int fibo_n = 1;

    // the Fibonacci sequence: n = 1,2,3,5,8,13,21,34...

    for (int ielem = 1; ielem <= nelems; ielem++)
    {
        fibo_n = fibo_0 + fibo_1;
        fibo_0 = fibo_1;
        fibo_1 = fibo_n;
    }
    return fibo_n;
}

bool TemperSubPermutation::check_range
(const int range, const int nelems, const int *array)
{
    for (int i = 0; i < nelems; i++)
      if (abs(i - array[i]) > range)
        return false;

    return true;
}

void TemperSubPermutation::swap_elements(const int ielem, const int jelem,
const int nelems, const int* array, int* ret_array)
{
  for (int i = 0; i < nelems; i++) ret_array[i] = array[i];

  ret_array[ielem] = array[jelem];
  ret_array[jelem] = array[ielem];
}

void TemperSubPermutation::init_subpermutations(const int nelems, int **permutations)
{
  // initialization

  const int range = 1;

  for (int iperm = 0; iperm < calc_nsubpermutations(nelems); iperm++) {
    for (int ielem = 0; ielem < nelems; ielem++) {
      permutations[iperm][ielem] = 0;
    }
  }

  // construct permutation function

  int nperm = 1;
  for (int ielem = 1; ielem < nelems; ielem++) {

    // # of permutation elements currently being constructed

    const int jelem = ielem + 1;

    // add "ielem" at the end of permutation function

    for (int jperm = 0; jperm < nperm; jperm++) {
      permutations[jperm][ielem] = ielem;
    }

    // exchange process

    int begin = 0;
    int ncandidates = nperm;

    for (int jswap = 0; jswap < ielem; jswap++) {
      int npermutation_add = 0;

      for (int kperm = begin; kperm < begin + ncandidates; kperm++) {
        int new_permutation[nelems];
        const int i = ielem - jswap;
        const int j = ielem - jswap - 1;
        swap_elements(i, j, nelems, permutations[kperm], new_permutation);

        // add new permutation if the candidate satisfies range condition
        if (check_range(range, jelem, new_permutation)) {
          for (int lelem = 0; lelem < nelems; lelem++) {
            permutations[nperm][lelem] = new_permutation[lelem];
          }
          nperm++;
          npermutation_add++;
        }
      }
      begin += ncandidates;
      ncandidates = npermutation_add;
    }
  }
}

int TemperSubPermutation::next_state(const int istate, const int nstate)
{
  return (istate + 1) < nstate ? istate + 1 : 0;
}

void TemperSubPermutation::suwa_todo_algorithm
(const int icurrent, const int nstates, const double* w, double* probs)
{
  double v [nstates];  // stochastic flows from the current state
  double ww[nstates];  // copy of weight factors

  for (int istate = 0; istate < nstates; istate++) {
    v [istate] = 0.0;
    ww[istate] = w[istate];
  }

  // search a state that has the maximum weight value

  int i_max  = 0;
  double w_max = w[0];
  for (int istate = 1; istate < nstates; istate++) {
    if (w_max < w[istate]) {
      i_max = istate;
      w_max = w[istate];
    }
  }

  // weight allocation process

  int istate, jstate, icount;

  icount = 0;
  istate = i_max;
  jstate = next_state(istate, nstates);

  while (icount < nstates) {

    double wr = w[istate];

    while (wr > 0.0e0) {
      if (wr >= ww[jstate]) {
        if (istate == icurrent) v[jstate] = ww[jstate];
        wr -= ww[jstate];
        jstate = next_state(jstate, nstates);
      }
      else {
        if (istate == icurrent) v[jstate] = wr;
        ww[jstate] = ww[jstate] - wr;
        wr = 0.0e0;
      }
    }
    istate = next_state(istate, nstates);
    icount++;
  }

  for (jstate = 0; jstate < nstates; ++jstate) {
    probs[jstate] = v[jstate] / w[icurrent];
  }
}