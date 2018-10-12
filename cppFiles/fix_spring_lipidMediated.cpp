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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_spring_lipidMediated.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "force.h"
#include "group.h"
#include "error.h"
#include <fstream>

using namespace LAMMPS_NS;
using namespace FixConst;
using std::ifstream;
#define SMALL 1.0e-10

enum{TETHER,COUPLE};

/* ---------------------------------------------------------------------- */

FixSpringLipidMediated::FixSpringLipidMediated(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  group2(NULL)
{
  if (narg < 10) error->all(FLERR,"Illegal fix spring command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  dynamic_group_allow = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  if (strcmp(arg[3],"couple") == 0) {
    if (narg != 12) error->all(FLERR,"Illegal fix spring command");
    styleflag = COUPLE;

    int n = strlen(arg[4]) + 1;
    group2 = new char[n];
    strcpy(group2,arg[4]);
    igroup2 = group->find(arg[4]);
    if (igroup2 == -1)
      error->all(FLERR,"Fix spring couple group ID does not exist");
    if (igroup2 == igroup)
      error->all(FLERR,"Two groups cannot be the same in fix spring couple");
    group2bit = group->bitmask[igroup2];

    k_spring = force->numeric(FLERR,arg[5]);
    xflag = yflag = zflag = 1;
    if (strcmp(arg[6],"NULL") == 0) xflag = 0;
    else xc = force->numeric(FLERR,arg[6]);
    if (strcmp(arg[7],"NULL") == 0) yflag = 0;
    else yc = force->numeric(FLERR,arg[7]);
    if (strcmp(arg[8],"NULL") == 0) zflag = 0;
    else zc = force->numeric(FLERR,arg[8]);
    r0 = force->numeric(FLERR,arg[9]);
    k_bin = force->numeric(FLERR,arg[10]);

    n = strlen(arg[11]) + 1;
    lipidName = new char[n];
    strcpy(lipidName,arg[11]);
    // if (r0 < 0) error->all(FLERR,"R0 < 0 for fix spring command");
    ifstream in_lipid(lipidName);
    if (!in_lipid) error->all(FLERR,"Lipid file doesn't exist");
    in_lipid >> left_bound >> right_bound;
    in_lipid >> a5 >> a4 >> a3 >> a2 >> a1 >> a0;
    in_lipid.close();

  } else error->all(FLERR,"Illegal fix spring command");

  ftotal[0] = ftotal[1] = ftotal[2] = ftotal[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixSpringLipidMediated::~FixSpringLipidMediated()
{
  delete [] group2;
}

/* ---------------------------------------------------------------------- */

int FixSpringLipidMediated::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpringLipidMediated::init()
{
  // recheck that group 2 has not been deleted

  if (group2) {
    igroup2 = group->find(group2);
    if (igroup2 == -1)
      error->all(FLERR,"Fix spring couple group ID does not exist");
    group2bit = group->bitmask[igroup2];
  }

  masstotal = group->mass(igroup);
  if (styleflag == COUPLE) masstotal2 = group->mass(igroup2);

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

}

/* ---------------------------------------------------------------------- */

void FixSpringLipidMediated::setup(int vflag)
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

void FixSpringLipidMediated::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringLipidMediated::post_force(int vflag)
{
  spring_couple();
}


/* ---------------------------------------------------------------------- */

void FixSpringLipidMediated::spring_couple()
{
  double xcm[3],xcm2[3];

  if (group->dynamic[igroup])
    masstotal = group->mass(igroup);

  if (group->dynamic[igroup2])
    masstotal2 = group->mass(igroup2);

  group->xcm(igroup,masstotal,xcm);
  group->xcm(igroup2,masstotal2,xcm2);

  // fx,fy,fz = components of k * (r-r0) / masstotal
  // fx2,fy2,fz2 = components of k * (r-r0) / masstotal2

  // double dx,dy,dz,fx,fy,fz,fx2,fy2,fz2,r,dr;
  //
  // dx = xcm2[0] - xcm[0] - xc;
  // dy = xcm2[1] - xcm[1] - yc;
  // dz = xcm2[2] - xcm[2] - zc;
  // if (!xflag) dx = 0.0;
  // if (!yflag) dy = 0.0;
  // if (!zflag) dz = 0.0;
  // r = sqrt(dx*dx + dy*dy + dz*dz);
  // r = MAX(r,SMALL);
  // dr = r - r0;
  double dx,dy,dz,fx,fy,fz,fx2,fy2,fz2,r,dr;

  double fprime;
  dx = xcm2[0] - xcm[0];
  dy = xcm2[1] - xcm[1];
  dz = xcm2[2] - xcm[2];
  if (!xflag) dx = 0.0;
  if (!yflag) dy = 0.0;
  if (!zflag) dz = 0.0;
  r = sqrt(dx*dx + dy*dy + dz*dz);
  r = r + r0; // r0 is the offset
  r = MAX(r,SMALL);
  // fprintf(stderr,"G1: %d, G2: %d, Center of masses distance: %f\n", igroup-5, igroup2-5, r);
  // dis 15
  // a0 = -2.1417e1;
  // a1 = 3.06674;
  // a2 = -1.46955e-01;
  // a3 = 2.88432e-03;
  // a4 = -1.99082e-5;
  // fprime = a0 + a1*r + a2*pow(r,2) + a3*pow(r,3) + a4*pow(r,4);
  // if(r > 43.03356301){
  //   fprime = 0.0;
  //   r = 43.03356301;
  // }
  // if(r < 14.32902487){
  //   fprime = 0.0;
  //   r = 14.32902487;
  // }
  // fprime = a0 + a1*r + a2*pow(r,2) + a3*pow(r,3) + a4*pow(r,4);
  // if(r > 43.03356301){
  //   fprime = 0.0;
  //   r = 43.03356301;
  // }
  // if(r < 14.32902487){
  //   fprime = 0.0;
  //   r = 14.32902487;
  // }
  // dis 10

  // a0 = -2.1417e1;
  // a1 = 3.06674;
  // a2 = -1.46955e-01;
  // a3 = 2.88432e-03;
  // a4 = -1.99082e-5;
  // fprime = a0 + a1*r + a2*pow(r,2) + a3*pow(r,3) + a4*pow(r,4);
  fprime = a1 + 2*a2*r + 3*a3*pow(r,2) + 4*a4*pow(r,3) + 5*a5*pow(r,4);
  // if(r > 38.03356){
  //   fprime = 0.0;
  //   r = 38.03356;
  // }
  // if(r < 9.32874381){
  //   fprime = 0.0;
  //   r = 9.32874381;
  // }
  if(r > right_bound){
    fprime = 0.0;
    r = right_bound;
  }
  if(r < left_bound){
    fprime = 0.0;
    r = left_bound;
  }
  // kt to kcal/mol 0.593
  double theta, theta_prime;
  double memb_b, z;
  double theta_1, theta_2, theta_prime_1, theta_prime_2;
  memb_b = zc/2.0;
  // k_bin = 0.2;
  z = xcm2[2];
  theta_1 = 0.5*(tanh(k_bin*(z+memb_b))-tanh(k_bin*(z-memb_b)));
  theta_prime_1 = 0.5*k_bin*(-pow(tanh(k_bin*(z+memb_b)),2) + pow(tanh(k_bin*(z-memb_b)),2));
  z = xcm[2];
  theta_2 = 0.5*(tanh(k_bin*(z+memb_b))-tanh(k_bin*(z-memb_b)));
  theta_prime_2 = 0.5*k_bin*(-pow(tanh(k_bin*(z+memb_b)),2) + pow(tanh(k_bin*(z-memb_b)),2));
  // // dis 15
  // a0 = 1.004588e2;
  // a1 = -2.1417e1;
  // a2 = 1.53337;
  // a3 = -4.8985e-2;
  // a4 = 7.2108e-4;
  // a5 = -3.98164e-6;
  // dis 10
  // a0 = 2.602268e1;
  // a1 = -9.40916215;
  // a2 = 9.017854e-1;
  // a3 = -3.5558475e-2;
  // a4 = 6.21541e-4;
  // a5 = -3.98164e-6;
  double ef;
  ef = 0.593*k_spring*(a0 + a1*r + a2*pow(r,2) + a3*pow(r,3) + a4*pow(r,4) + a5*pow(r,5));
  espring = ef*theta_1*theta_2;

  fx = 0.593*k_spring*dx*fprime/r*theta_1*theta_2;
  fy = 0.593*k_spring*dy*fprime/r*theta_1*theta_2;
  fz = 0.593*k_spring*dz*fprime/r*theta_1*theta_2;
  ftotal[0] = fx;
  ftotal[1] = fy;
  ftotal[2] = fz;
  ftotal[3] = sqrt(fx*fx + fy*fy + fz*fz);
  // if (dr < 0.0) ftotal[3] = -ftotal[3];



  if (masstotal2 > 0.0) {
    fx2 = fx/masstotal2;
    fy2 = fy/masstotal2;
    fz2 = (fz + ef*theta_prime_1*theta_2)/masstotal2;
  } else fx2 = fy2 = fz2 = 0.0;

  if (masstotal > 0.0) {
    fx /= masstotal;
    fy /= masstotal;
    fz = (fz - ef*theta_prime_2*theta_1)/ masstotal;
  } else fx = fy = fz = 0.0;

  // apply restoring force to atoms in group
  // f = -k*(r-r0)*mass/masstotal

  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double massone;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        massone = rmass[i];
        f[i][0] += fx*massone;
        f[i][1] += fy*massone;
        f[i][2] += fz*massone;
      }
      if (mask[i] & group2bit) {
        massone = rmass[i];
        f[i][0] -= fx2*massone;
        f[i][1] -= fy2*massone;
        f[i][2] -= fz2*massone;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        massone = mass[type[i]];
        f[i][0] += fx*massone;
        f[i][1] += fy*massone;
        f[i][2] += fz*massone;
      }
      if (mask[i] & group2bit) {
        massone = mass[type[i]];
        f[i][0] -= fx2*massone;
        f[i][1] -= fy2*massone;
        f[i][2] -= fz2*massone;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringLipidMediated::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringLipidMediated::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of stretched spring
------------------------------------------------------------------------- */

double FixSpringLipidMediated::compute_scalar()
{
  return espring;
}

/* ----------------------------------------------------------------------
   return components of total spring force on fix group
------------------------------------------------------------------------- */

double FixSpringLipidMediated::compute_vector(int n)
{
  return ftotal[n];
}
