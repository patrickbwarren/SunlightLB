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

#define MPOBJS_C

#include <stdio.h>
#include <stdlib.h>
#include "sunlightlb.h"

/* Insert diffusing material on all fluid nodes */

void mpall(MPSys mp, LBReal pval) {
  int i, j, k, esi, esj, esk, inode, ns;
  LBObj *obj;
  LBReal *p;
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; p = mp->p; ns = 0;
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    if (obj[inode] == NULL) { p[inode] = pval; ns++; }
  }
  if (slb_verbosity & MP_OBJECT_DETAILS) {
    printf("Diffusing material slab inserted everywhere, ");
    printf("%i sites set\n", ns);
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : mpall done\n", mp->step);
}




/* Insert a slab of diffusing material from x1, y1, z1 to x2, y2, z2 */

void mpslab(MPSys mp, LBReal x1, LBReal x2, LBReal y1, LBReal y2, LBReal z1, LBReal z2, LBReal pval) {
  int i, j, k, esi, esj, esk, inode, ns;
  LBObj *obj;
  LBReal *p;
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; p = mp->p; ns = 0;
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    if (i > x1 && i < x2 && j > y1 && j < y2 && k > z1 && k < z2) { 
      if (obj[inode] == NULL) { p[inode] = pval; ns++; }
    }
  }
  if (slb_verbosity & MP_OBJECT_DETAILS) {
    printf("Diffusing material slab from %5.2f %5.2f %5.2f ", x1, y1, z1);
    printf("to %5.2f %5.2f %5.2f, ", x2, y2, z2);
    printf("%i sites set\n", ns);
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : mpslab done\n", mp->step);
}


/* Insert a wall of diffusing material,
 * direction determined by first letter of type.
 */

void mpwall(MPSys mp, char *type, int wallpos, LBReal pval) {
  int i, j, k, esi, esj, esk, inode, ns;
  LBObj *obj;
  LBReal *p;
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; p = mp->p; ns = 0;
  switch (type[0]) {
    case 'x': case 'i':
      if ((i = wallpos) < 0 || i > esi-1) {
        printf("wallpos out of range\n"); return;
      }
      for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
	inode = esk*(esj*i+j)+k;
	if (obj[inode] == NULL) { p[inode] = pval; ns++; }
      }
      break;
    case 'y': case 'j':
      if ((j = wallpos) < 0 || j > esj-1) {
        printf("wallpos out of range\n"); return;
      }
      for (i=0; i<esi; i++) for (k=0; k<esk; k++) {
	inode = esk*(esj*i+j)+k;
	if (obj[inode] == NULL) { p[inode] = pval; ns++; }
      }
      break;
    case 'z': case 'k':
      if ((k = wallpos) < 0 || k > esk-1) {
        printf("wallpos out of range\n"); return;
      }
      for (i=0; i<esi; i++) for (j=0; j<esj; j++) {
	inode = esk*(esj*i+j)+k;
	if (obj[inode] == NULL) { p[inode] = pval; ns++; }
      }
      break;
    default: 
      printf("Unrecognised direction for walls\n"); return;
  }
  if (slb_verbosity & MP_OBJECT_DETAILS) {
    printf("Diffusing material wall at %c = %i, ", type[0], wallpos);
    printf("%i sites set\n", ns);
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : mpwall done\n", mp->step);
}


/* Insert a sphere of diffusing material radius rad, centre ox, y, z */

void mpsph(MPSys mp, LBReal rad, LBReal ox, LBReal oy, LBReal oz, LBReal pval) {
  int i, j, k, esi, esj, esk, inode, ns;
  LBObj *obj;
  LBReal *p;
  LBReal rsq, radsq;
  radsq = rad*rad;
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; p = mp->p; ns = 0;
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    rsq = (i-ox)*(i-ox) + (j-oy)*(j-oy) + (k-oz)*(k-oz);
    if (rsq < radsq) { 
      if (obj[inode] == NULL) { p[inode] = pval; ns++; }
    }
  }
  if (slb_verbosity & MP_OBJECT_DETAILS) {
    printf("Diffusing material sphere, radius %5.2f at %5.2f %5.2f %5.2f, ",
           rad, ox, oy, oz);
    printf("%i sites set\n", ns);
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : mpsph done\n", mp->step);
}


/* Insert a cylinder of diffusing material radius rad,
 * axis (a, b) = (ox, oy) or (oy, oz) or (ox, oz),
 * direction determined by first letter of type,
 * whether solid or a cylindrical hole by the second letter.
 */

void mpcyl(MPSys mp, char *type, LBReal rad, LBReal a, LBReal b, LBReal pval) {
  int i, j, k, esi, esj, esk, inode, ns, within;
  LBObj *obj;
  LBReal *p;
  LBReal rsq, radsq;
  radsq = rad*rad; rsq = 0.0;
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; p = mp->p; ns = 0;
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    switch (type[0]) {
      case 'x': case 'i': rsq = (j-a)*(j-a) + (k-b)*(k-b); break;
      case 'y': case 'j': rsq = (i-a)*(i-a) + (k-b)*(k-b); break;
      case 'z': case 'k': rsq = (i-a)*(i-a) + (j-b)*(j-b); break;
      default: error("Unrecognised cylinder direction");
    }
    within = (rsq < radsq); if (type[1] != '\0') within = !within;
    if (within) { 
      if (obj[inode] == NULL) { p[inode] = pval; ns++; }
    }
  }
  if (slb_verbosity & MP_OBJECT_DETAILS) {
    if (type[1] == '\0') printf("Diffusing material cylinder, \n");
    else printf("Diffusing material cylindrical hole, \n");
    printf("radius %5.2f, axis %5.2f %5.2f along %c, ", rad, a, b, type[0]);
    printf("%i sites set\n", ns);
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : mpcyl done\n", mp->step);
}


/* Insert the substrate field specified in subsfile. The format is
 * assumed to be a byte image, 0x00 -> 0.0, 0xff -> pval, and
 * other values are by interpolation.
 */

void mpsubs(MPSys mp, char *subsfile, int ilo, int ihi, int jlo, int jhi, int klo, int khi, LBReal pval) {
  unsigned char byte;
  int i, j, k, esj, esk, inode, ns;
  LBObj *obj;
  LBReal *p, ptot, val;
  FILE *fp;
  if ((fp = fopen(subsfile,"r")) == NULL) {
    printf("mpsubs: %s could not be opened\n", subsfile); return;
  }
  if (slb_verbosity & MP_OBJECT_DETAILS) {
    printf("\nInserting mp substrate data from %s\n", subsfile);
    printf("Substrate region is (%i-%i,%i-%i,%i-%i) inclusive\n",
           ilo, ihi, jlo, jhi, klo, khi);
  }
/* Read in data and insert in map, */
/* note array size must match that of data file */
/* also note k varies fastest */
  esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; p = mp->p; ns = 0; ptot = 0.0;
  for (i=ilo; i<=ihi; i++)
    for (j=jlo; j<=jhi; j++)
      for (k=klo; k<=khi; k++) {
	inode = esk*(esj*i+j)+k;
	if (fscanf(fp,"%c",&byte) != 1) {
	  printf("mpsubs: ran out of data in %s\n", subsfile);
	  break;
	}
	if (obj[inode] == NULL) {
	  val = pval * (double)byte / 255.0;
	  p[inode] = val; ns++; ptot += val;
	}
      }
  fclose(fp);
  if (slb_verbosity & MP_OBJECT_DETAILS) {
    printf("Inserted mp substrate field, %i sites set\n", ns);
    printf("Total amount of added material is %f\n", ptot);
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : mpsubs done\n", mp->step);
}


/* Insert the substrate field specified in subsfile. The format is
 * assumed to be a list of points to be initialised to pval.
 */

void mplist(MPSys mp, char *subsfile, LBReal pval) {
  int i, j, k, esj, esk, inode, ns;
  LBObj *obj;
  LBReal *p, ptot;
  FILE *fp;
  if ((fp = fopen(subsfile,"r")) == NULL) {
    printf("mpsubs: %s could not be opened\n", subsfile); return;
  }
  if (slb_verbosity & MP_OBJECT_DETAILS) {
    printf("\nInserting mp list data from %s\n", subsfile);
  }
/* Read in data and insert in map, */
/* also note k varies fastest */
  esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; p = mp->p; ns = 0; ptot = 0.0;
  for (;;) {
    if (fscanf(fp,"%i%i%i",&i,&j,&k) != 3) break;
    inode = esk*(esj*i+j)+k;
    if (obj[inode] == NULL) {
      p[inode] = pval; ns++; ptot += pval;
    }
  }
  fclose(fp);
  if (slb_verbosity & MP_OBJECT_DETAILS) {
    printf("Inserted mp substrate field, %i sites set\n", ns);
    printf("Total amount of added material is %f\n", ptot);
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : mpsubs done\n", mp->step);
}

/* End of mpobjs.c */
