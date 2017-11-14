/* parallel remap functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level directory of the distribution.
*/

#include "stdio.h"
#include "mpi.h"

#include "pack_3d.h"
#include "remap_3d.h"

#ifndef BLUEGENE
   #define remap_3d_interface remap_3d_interface_
   #define remap_3d_create_plan_interface remap_3d_create_plan_interface_
   #define remap_3d_destroy_plan_interface remap_3d_destroy_plan_interface_
#endif

/* ------------------------------------------------------------------- */
/* F77 wrapper on remap */

void remap_3d_interface(double *in, double *out, double *buf,
	       struct remap_plan_3d **plan)

{
  remap_3d(in,out,buf,*plan);
}

/* ------------------------------------------------------------------- */
/* F77 wrapper on remap_create_plan */

void remap_3d_create_plan_interface(
       MPI_Comm *comm,
       int *in_ilo, int *in_ihi, int *in_jlo, int *in_jhi,
       int *in_klo, int *in_khi,
       int *out_ilo, int *out_ihi, int *out_jlo, int *out_jhi,
       int *out_klo, int *out_khi,
       int *nqty, int *permute, int *memory, int *precision,
       struct remap_plan_3d **plan)

{
  int me;

/* convert F77 indices to C */

  *plan = remap_3d_create_plan(*comm,
			       *in_ilo-1,*in_ihi-1,*in_jlo-1,*in_jhi-1,
			       *in_klo-1,*in_khi-1,
			       *out_ilo-1,*out_ihi-1,*out_jlo-1,*out_jhi-1,
			       *out_klo-1,*out_khi-1,
			       *nqty, *permute, *memory, *precision);
  if (*plan == NULL) {
    MPI_Comm_rank(*comm,&me);
    printf("ERROR: REMAP 3d plan is NULL on proc %d\n",me);
  }
}

/* ------------------------------------------------------------------- */
/* F77 wrapper on remap_destroy_plan */

void remap_3d_destroy_plan_interface(struct remap_plan_3d **plan)

{
  remap_3d_destroy_plan(*plan);
}
