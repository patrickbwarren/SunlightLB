/*
 * This file is part of SunlightLB - a 3D Lattice Boltzmann code
 * Copyright (C) 2005 Unilever UK Central Resources Ltd.
 *  
 * SunlightLB is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2 
 * of the License, or (at your option) any later version. 
 *
 * SunlightLB is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
 */

#define MPMISC_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "sunlightlb.h"

/* Local variables and function prototypes */

static time_t tstart;

static LBReal *ralloc(int);
static int *ialloc(int);

/* Deliver a new MPSys structure with all space initialised */

MPSys newmp(LBSys lb, char *name) {
  int inode;
  MPSys mp;
  LBReal *p, *newp, *rho, rho0;
  mp = (MPSys) malloc(sizeof(MPSys_type));
  if (mp == NULL) error("No room for mpsys structure\n");
  mp->lb = lb;
  mp->p = ralloc(lb->nnodes);
  mp->newp = ralloc(lb->nnodes);
  if (mp->lb->rho != NULL) free(mp->lb->rho);
  mp->lb->rho = ralloc(lb->nnodes);
  p = mp->p; newp = mp->newp; rho = mp->lb->rho; rho0 = lb->rho0;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    p[inode] = newp[inode] = 0.0; rho[inode] = rho0;
  }
  mp->objects_moving = lb->objects_moving;
  mp->to_be_initialised = 1;
  if (myscmp("null", name)) mp->name = NULL; else mp->name = mycopy(name);
  printf("MP '%s' initialised for %i * %i * %i lattice\n",
	 mp->name, mp->lb->esi, mp->lb->esj, mp->lb->esk);
  return mp;
}


/* (Re)initialise the density field required for moment propagation */

void mprhoinit(MPSys mp) {
  int i, inode;
  LBReal *a, *rho, rho0;
  LBObj *obj;
  rho0 = mp->lb->rho0; rho = mp->lb->rho;
  a = mp->lb->a; obj = mp->lb->obj;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    rho[inode] = rho0;
    if (obj[inode] == NULL) {
      for (i=0; i<15; i++) rho[inode] += a[i+15*inode];
    }
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mprhoinit done\n", mp->name, mp->step);
  }
}

/* Initialise the propagating scalar field, and other stuff. */

void mpinit(MPSys mp) {
  int inode;
  LBReal *p;
/* First calculate density */
  mprhoinit(mp);
/* Clear out the arrays */
  p = mp->p;
  for (inode=0; inode<(mp->lb->nnodes); inode++) p[inode] = 0.0;
/* Reset the step count and set the status */
  mp->step = 0; mp->to_be_initialised = 0;
/* Restart the clock, so monitoring can get a reasonable idea of execution time */
  tstart = time(NULL); 
  printf("MP '%s' initialised and ready to go.\n", mp->name);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpinit done\n", mp->name, mp->step);
  }
}


/* Write out parameters for info */

void mpwrpars(MPSys mp) {
  LBReal d;
  d = (1.0 - mp->alpha) / 6.0;
  printf("\nMP '%s' parameter list\n\n", mp->name);
  printf(" | LB size  = %i * %i * %i = %i\n",
	 mp->lb->esi, mp->lb->esj, mp->lb->esk, mp->lb->nnodes);
  printf(" | alpha, D = %0.4f %0.6f\n", mp->alpha, d);
  if (mp->objects_moving) printf(" | Objects allowed to move\n");
  else printf(" | All objects are stationary\n");
  printf("\n");
}


/* Check the total propagating scalar field, which should be conserved,
 * at least if objects are stationary.
 */

void mpcheck(MPSys mp) {
  int inode;
  LBReal *p;
  LBReal pfl, psl;
  LBObj *obj;
  obj = mp->lb->obj; p = mp->p;
  pfl = psl = 0.0;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    if (obj[inode] == NULL) pfl += p[inode]; else psl += p[inode];
  }
  printf("MP '%s' %i : fluid + solid scalar = %10.6f + %10.6f = %10.6f\n",
         mp->name, mp->step, pfl, psl, pfl+psl);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpcheck done\n", mp->name, mp->step);
  }
}


/* Return the mean value of the propagating scalar field.
 */

LBReal mpmean(MPSys mp) {
  int inode;
  LBReal ptot, pmean;
  LBReal *p;
  LBObj *obj;
  obj = mp->lb->obj; p = mp->p;
  ptot = 0.0;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    if (obj[inode] == NULL) ptot += p[inode];
  }
  pmean = ptot / (LBReal)(countfluid(mp->lb));
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpnorm done, ptot = %f, pmean = %f\n",
	   mp->name, mp->step, ptot, pmean);
  }
  return pmean;
}


/* Return the maximum value of the propagating scalar field.
 */

LBReal mpmax(MPSys mp) {
  int inode;
  int found = 0;
  LBReal pval, pmax = 0.0;
  LBReal *p;
  LBObj *obj;
  obj = mp->lb->obj; p = mp->p;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    if (obj[inode] == NULL) {
      pval = p[inode];
      if (found) { if (pmax < pval) pmax = pval; }
      else { pmax = pval; found = 1; }
    }
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpnorm done, pmax = %f\n",
	   mp->name, mp->step, pmax);
  }
  return pmax;
}


/* Return the minimum value of the propagating scalar field.
 */

LBReal mpmin(MPSys mp) {
  int inode;
  int found = 0;
  LBReal pval, pmin = 0.0;
  LBReal *p;
  LBObj *obj;
  obj = mp->lb->obj; p = mp->p;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    if (obj[inode] == NULL) {
      pval = p[inode];
      if (found) { if (pmin > pval) pmin = pval; }
      else { pmin = pval; found = 1; }
    }
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpnorm done, pmin = %f\n",
	   mp->name, mp->step, pmin);
  }
  return pmin;
}


/* Rescale the propagating scalar.
 */

void mprescale(MPSys mp, LBReal fac) {
  int inode;
  LBReal *p;
  LBObj *obj;
  obj = mp->lb->obj; p = mp->p;
  for (inode=0; inode<(mp->lb->nnodes); inode++) if (obj[inode] == NULL) {
    p[inode] *= fac;
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mprescale done\n", mp->name, mp->step);
  }
}


/* Copy the propagating scalar from one MP system to another
 * Assumes that both belong to the same LB systems, but this
 * is not checked.
 */

void mpcopy(MPSys mp, MPSys mp_dest) {
  int inode;
  LBReal *p;
  LBReal *p_dest;
  LBObj *obj;
  obj = mp->lb->obj; p = mp->p; p_dest = mp_dest->p;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    if (obj[inode] == NULL) p_dest[inode] += p[inode];
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpcopy done, into MP '%s'\n", 
	   mp->name, mp->step, mp_dest->name);
  }
}


/* Save out the propagating scalar field, (very similar to 
 * mpconwr).  For use with companion routine mpload.
 */

void mpsave(MPSys mp, char *file, char *type) {
  int inode;
  int saveasbinary = 0;
  LBReal *p;
  LBStoreReal val;
  FILE *fp;
  if (type[0] == 'b' || type[0] == 'B') saveasbinary = 1;
  if ((fp = fopen(file,"w")) == NULL) {
    printf("mpsave: %s could not be opened\n", file);
    return;
  }
  printf("MP '%s' %6i: Saving mp data to %s\n", mp->name, mp->step, file);
  p = mp->p;
/* Save on a per-node basis.  Equivalent to traversing */
/* p[i][j][k] array with k varying fastest */
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    val = p[inode];
    if (saveasbinary) fwrite(&val, sizeof(LBStoreReal),1, fp);
    else fprintf(fp,"%f\n", val);
  }
  fclose(fp);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpsave done\n", mp->name, mp->step);
  }
}


/* Load the propagating scalar field.
 * For use with companion routine mpsave.
 */

void mpload(MPSys mp, char *file, char *type) {
  int inode;
  int loadfrombinary = 0;
  LBReal *p;
  LBStoreReal val;
  FILE *fp;
  if (type[0] == 'b' || type[0] == 'B') loadfrombinary = 1;
  if ((fp = fopen(file,"r")) == NULL) {
    printf("mpload: %s could not be opened\n", file);
    exit(1);
  }
  printf("MP '%s' %6i: Loading mp data from %s\n", mp->name, mp->step, file);
/* load on a per-node basis.  Equivalent to filling */
/* p[i][j][k] array with k varying fastest */
  p = mp->p;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    if (loadfrombinary) fread(&val, sizeof(LBStoreReal),1, fp);
    else fscanf(fp,"%f",&val);
    p[inode] = val;
  }
  fclose(fp);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpload done\n", mp->name, mp->step);
  }
}


/* Write out the propagating scalar field, 'method' codes as follows: 
 *  method[0] = 'u' --> unique file name constructed.
 *  method[0] = 'n' --> basic file name used as supplied.
 * In the former case, the file name must be a printf format string containing %i.
 * In the latter case (basic file name used) then:
 *  method[1] = 'w' --> overwrite the file (default).
 *  method[1] = 'a' --> append to the file.
 * The final letter in method is used to determine the file type.
 *  method[1|2] = 'a' --> ascii storage (default).
 *  method[1|2] = 'b' --> binary storage.
 *  method[1|2] = 'c' --> character storage.
 * In the last case, the concentration (a floating point variable)
 * is converted to a byte in the range 0x00 to 0xff before writing. 
 * Compared to binary, this reduces the file size by a factor four.  The byte
 * value is constructed by linearly interpolating the value of p between p2bytemin
 * and p2bytemax, clamping to this range if necessary.
 */

void mpconwr(MPSys mp, char *mpconfile, char *method) {
  int inode, ival;
  char cval;
  char file[50];
  char meth[8] = "w";
  char filetype = 'a';
  LBReal pmin, pmax, *p;
  LBStoreReal val;
  FILE *fp;
  LBObj *obj;
  switch (method[0]) {
    case 'u':
      sprintf(file, mpconfile, mp->step);
      filetype = method[1];
      break;
    case 'n':
      strcpy(file, mpconfile);
      if (method[1] != '\0') {
	if (method[1] == 'a') strcpy(meth,"a");
	filetype = method[2];
      } else filetype = 'a';
      break;
    default:
      error("Unrecognised method code (1st letter) in mpconwr");
  }
  if (filetype != 'b' && filetype != 'c') filetype = 'a';
/* Make sure p2bytemin and p2bytemax are sensible */
  pmin = mp->p2bytemin; pmax = mp->p2bytemax;
  if (pmin >= pmax) {
    printf("Switching to p2bytemin = 0.0 and p2bytemax = 1.0\n");
    pmin = 0.0; pmax = 1.0; 
  }
  if ((fp = fopen(file, meth)) == NULL) {
    printf("mpconfile: %s could not be opened\n", file);
    return;
  }
  printf("MP '%s' %6i: ", mp->name, mp->step);
  if (meth[0] == 'w') printf("Writing");
  if (meth[0] == 'a') printf("Appending");
  switch (filetype) {
  case 'a' : printf(" ascii"); break;
  case 'b' : printf(" binary"); break;
  case 'c' : printf(" byte-converted"); break;
  }
  printf(" scalar density to %s\n", file);
/* Save on a per-node basis.  Equivalent to traversing */
/* p[i][j][k] array with k varying fastest */
  obj = mp->lb->obj; p = mp->p;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    if (obj[inode] == NULL) val = p[inode]; else val = 0.0;
    switch (filetype) {
    case 'c' :
      ival = (int)(256.0*(val-pmin)/(pmax-pmin));
      if (ival > 255) ival = 255; if (ival < 0) ival = 0;
      cval = (char)ival;
      fwrite(&cval, sizeof(char),1, fp);
      break;
    case 'b' :
      fwrite(&val, sizeof(LBStoreReal),1, fp);
      break;
    case 'a' :
      fprintf(fp,"%f\n", val);
      break;
    }
  }
  fclose(fp);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpconwr done\n", mp->name, mp->step);
  }
}


/* Set aside storage for an isosurface field */

MPIso newmpiso(MPSys mp, LBReal val) {
  int inode;
  MPIso iso;
  int *count;
  iso = (MPIso) malloc(sizeof(MPIso_type));
  if (iso == NULL) {
    fprintf(stderr,"No room at the inn for another isosurface\n");
    return NULL;
  }
  count = ialloc(mp->lb->nnodes);
  for (inode=0; inode<(mp->lb->nnodes); inode++) count[inode] = 0;
  iso->mp = mp; 
  iso->val = val;
  iso->count = count;
  return iso;
}

/* Store an isoslice through the concentration field */

void mpiso(MPIso iso) {
  int inode;
  int *count;
  LBReal *p, val;
  p = iso->mp->p; val = iso->val; count = iso->count;
  for (inode=0; inode<(iso->mp->lb->nnodes); inode++) if (p[inode] > val) count[inode]++;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpiso done\n", iso->mp->name, iso->mp->step);
  }
}


/* Write out the isosurface array */

void mpisowr(MPIso iso, char *isofile, char *type) {
  int inode, val,*count;
  int saveasbinary = 0;
  FILE *fp;
  if (type[0] == 'b' || type[0] == 'B') saveasbinary = 1;
  if ((fp = fopen(isofile,"w")) == NULL) {
    printf("mpisowr: %s could not be opened\n", isofile);
    return;
  }
  printf("MP '%s' %6i: Saving cumulative isosurface data to %s\n",
         iso->mp->name, iso->mp->step, isofile);
  count = iso->count;
/* Save on a per-node basis.  Equivalent to traversing */
/* count[i][j][k] array with k varying fastest */
  for (inode=0; inode<(iso->mp->lb->nnodes); inode++) {
    val = count[inode];
    if (saveasbinary) fwrite(&val, sizeof(int),1, fp);
    else fprintf(fp,"%i\n", val);
  }
  fclose(fp);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpisowr done\n", iso->mp->name, iso->mp->step);
  }
}


/* Measure the spread (ie variance) of the concentration profile in
 * a desired direction, specified by dirn.
 */

LBReal mpspread(MPSys mp, char *dirn) {
  int inode, i, j, k, esi, esj, esk, fx, fy, fz, pos;
  LBReal tot, mean, var, ptemp;
  LBReal *p;
  LBObj *obj;
  fx = fy = fz = 0;
  switch (dirn[0]) {
    case 'x': case 'i': fx = 1; break;
    case 'y': case 'j': fy = 1; break;
    case 'z': case 'k': fz = 1; break;
    default: printf("Unrecognised direction for measuring spread\n"); return -1.0;
  }
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; p = mp->p;
  tot = mean = var = 0.0;
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    if (obj[inode] == NULL) {
      ptemp = p[inode];
      pos = fx*i + fy*j + fz*k;
      tot  += ptemp;
      mean += ptemp * pos;
      var  += ptemp * pos * pos;
    }
  }
  if (tot > 0.0) {
    mean = mean / tot;
    var = var / tot - mean*mean;
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpspread done\n", mp->name, mp->step);
  }
  return var;
}

/* Routines to allocate space */

LBReal *ralloc(int n) {
  LBReal *rarr;
  rarr = (LBReal *) malloc(n*sizeof(LBReal));
  if (rarr == NULL) error("No room at the inn");
  return rarr;
}

int *ialloc(int n) {
  int *iarr;
  iarr = (int *) malloc(n*sizeof(int));
  if (iarr == NULL) error("No room at the inn");
  return iarr;
}

/* End of mpmisc.c */
