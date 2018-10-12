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
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
                        Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_spring_rg_cylindrical_capped.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSpringRGCylindricalCapped::FixSpringRGCylindricalCapped(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix spring/rg/cylindricalCapped command");
  scalar_flag = 1;
  global_freq = 1;

  k = force->numeric(FLERR,arg[3]);
  rg0_flag = 0;
  rg0 = force->numeric(FLERR,arg[4]);
  k_bin = force->numeric(FLERR,arg[5]);
  dynamic_group_allow = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
}

/* ---------------------------------------------------------------------- */

int FixSpringRGCylindricalCapped::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpringRGCylindricalCapped::init()
{
  masstotal = group->mass(igroup);
  force_flag = 0;
  // if rg0 was specified as NULL, compute current Rg
  // only occurs on 1st run

  if (rg0_flag) {
    double xcm[3];
    group->xcm(igroup,masstotal,xcm);
    rg0 = group->gyration(igroup,masstotal,xcm);
    rg0_flag = 0;
  }

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringRGCylindricalCapped::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringRGCylindricalCapped::post_force(int vflag)
{
  // compute current Rg and center-of-mass
  force_flag = 0;
  double xcm[3];
  if (group->dynamic[igroup])
    masstotal = group->mass(igroup);
  group->xcm(igroup,masstotal,xcm);
  double rg = group->cylindricalgyration(igroup,masstotal,xcm);

  // apply restoring force to atoms in group
  // f = -k*(r-r0)*mass/masstotal

  double dx,dy,term1;

  double **f = atom->f;
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double massfrac;
  double unwrap[3];
  double theta, theta_prime;
  theta = 0.5* (1 + tanh(k_bin*(rg - rg0)));
  theta_prime = (1 - pow(tanh(k_bin*(rg - rg0)), 2))*k_bin;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xcm[0];
      dy = unwrap[1] - xcm[1];
      term1 = k * (2*(1.0 - theta) + (rg0 * rg0 / rg - rg)*theta_prime);
      if (masstotal > 0.0) {
        if (rmass) massfrac = rmass[i]/masstotal;
        else  massfrac = mass[type[i]]/masstotal;

        f[i][0] -= term1*dx*massfrac;
        f[i][1] -= term1*dy*massfrac;
      }
    }
  total_energy = k * rg * rg * (1 - theta) + k * rg0 * rg0 * theta;
  // fprintf(stderr,"Rg: %f, total energy: %f\n", rg, total_energy);
}

/* ---------------------------------------------------------------------- */

void FixSpringRGCylindricalCapped::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}
double FixSpringRGCylindricalCapped::compute_scalar()
{
  return total_energy;
}
