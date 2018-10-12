/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(spring/lipidMediated,FixSpringLipidMediated)

#else

#ifndef LMP_FIX_SPRING_LIPIDMEDIATED_H
#define LMP_FIX_SPRING_LIPIDMEDIATED_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSpringLipidMediated : public Fix {
 public:
  FixSpringLipidMediated(class LAMMPS *, int, char **);
  ~FixSpringLipidMediated();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 private:
  double xc,yc,zc,r0;
  double k_spring;
  double k_bin;
  int xflag,yflag,zflag;
  int styleflag;
  char *group2;
  char *lipidName;
  int igroup2,group2bit;
  double masstotal,masstotal2;
  double left_bound, right_bound;
  double a0,a1,a2,a3,a4,a5;
  int ilevel_respa;
  double espring,ftotal[4];
  void spring_couple();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: R0 < 0 for fix spring command

Equilibrium spring length is invalid.

E: Fix spring couple group ID does not exist

Self-explanatory.

E: Two groups cannot be the same in fix spring couple

Self-explanatory.

*/
