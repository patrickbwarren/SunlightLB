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

#define MPADS_C

#include <stdio.h>
#include <stdlib.h>
#include "sunlightlb.h"

/* SUGGEST USING LB SYSTEM RATHER THAN MP SYSTEM IN THESE CALLS */

/* Adsorption sites are indicated by bit 0x4000 in mask.
 * Mask is only non-zero for solid nodes with boundary links,
 * this is preserved by the following suite of routines for
 * specifying adsorption sites.
 */


/* Make all surface nodes into adsorption sites */

void adall(MPSys mp) {
  int inode, ns;
  int *mask;
  mask = mp->lb->mask; ns = 0;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    if (mask[inode] != 0) { mask[inode] |= 0x4000; ns++; }
  }
  printf ("All nodes are adsorption sites, %i sites set\n", ns);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : adall done\n", mp->step);
}


/* Make all surface nodes into non-adsorption sites */

void adnone(MPSys mp) {
  int inode, ns;
  int *mask;
  mask = mp->lb->mask; ns = 0;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    if (mask[inode] & 0x4000) { mask[inode] &= ~0x4000; ns++; }
  }
  printf ("Reset %i adsorption sites to non-adsorbing sites\n", ns);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP %i : adnone done\n", mp->step);
  }
}


/* Insert a slab of adsorption sites from x1, y1, z1 to x2, y2, z2 */

void adslab(MPSys mp, LBReal x1, LBReal x2, LBReal y1, LBReal y2, LBReal z1, LBReal z2) {
  int i, j, k, esi, esj, esk, inode, ns;
  int *mask;
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  mask = mp->lb->mask; ns=0;
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    if (i > x1 && i < x2 && j > y1 && j < y2 && k > z1 && k < z2) { 
      if (mask[inode] != 0) { mask[inode] |= 0x4000; ns++; }
    }
  }
  printf("Adsorption site slab from %5.2f %5.2f %5.2f ", x1, y1, z1);
  printf("to %5.2f %5.2f %5.2f, ", x2, y2, z2);
  printf("%i sites set\n", ns);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : adslab done\n", mp->step);
}


/* Insert a sphere of adsorption sites radius rad, centre ox, y, z */

void adsph(MPSys mp, LBReal rad, LBReal ox, LBReal oy, LBReal oz) {
  int i, j, k, esi, esj, esk, inode, ns;
  int *mask;
  LBReal rsq, radsq;
  radsq = rad*rad;
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  mask = mp->lb->mask; ns=0;
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    rsq = (i-ox)*(i-ox) + (j-oy)*(j-oy) + (k-oz)*(k-oz);
    if (rsq < radsq) { 
      if (mask[inode] != 0) { mask[inode] |= 0x4000; ns++; }
    }
  }
  printf("Adsorption site sphere, radius %5.2f at %5.2f %5.2f %5.2f, ",
         rad, ox, oy, oz);
  printf("%i sites set\n", ns);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : adsph done\n", mp->step);
}


/* Insert a cylinder of adsorption sites radius rad,
 * axis (a, b) = (ox, oy) or (oy, oz) or (ox, oz),
 * direction determined by first letter of type,
 * whether solid or a cylindrical hole by the second letter.
 */

void adcyl(MPSys mp, char *type, LBReal rad, LBReal a, LBReal b) {
  int i, j, k, esi, esj, esk, inode, within, ns;
  int *mask;
  LBReal rsq, radsq;
  radsq = rad*rad; rsq = 0.0;
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  mask = mp->lb->mask; ns=0;
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
      if (mask[inode] != 0) { mask[inode] |= 0x4000; ns++; }
    }
  }
  if (type[1] == '\0') printf("Adsorbtion site cylinder, \n");
  else printf("Adsorbtion site cylindrical hole, \n");
  printf("radius %5.2f, axis %5.2f %5.2f along %c, ", rad, a, b, type[0]);
  printf("%i sites set\n", ns);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : adcyl done\n", mp->step);
}


/* Compute total adsorbed amount, write out, and return value 
 * If called with "null" for the file name, then do not attempt to write.
 * Use private static int nadt to count how many times this routine is
 * called, and write out multiplied by a given factor.  This allows for 
 * multiple steps to be taken between calls.
 */

LBReal adtotal(MPSys mp, char *adtfile, int imult) {
  int inode;
  int *mask;
  static int nadt = 0;
  char method[2] = " ";
  LBReal *p, ptot;
  FILE *fp;
  mask = mp->lb->mask; p = mp->p; ptot = 0.0;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    if (mask[inode] & 0x4000) ptot += p[inode];
  }
  if (nadt == 0) method[0] = 'w'; else method[0] = 'a';
  if (!myscmp("null", adtfile)) {
    if ((fp = fopen(adtfile, method)) == NULL) {
      printf("adtotal: %s could not be opened to write/append data\n", adtfile);
    } else {
      fprintf(fp,"%8i %f\n", nadt*imult, ptot);
      fclose(fp);
    }
  }
  nadt++;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : adtotal done\n", mp->step);
  return ptot;
}


/* Reset the adsorption site totals to zero */

void adzero(MPSys mp) {
  int inode;
  int *mask;
  LBReal *p;
  mask = mp->lb->mask; p = mp->p;
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    if (mask[inode] & 0x4000) p[inode] = 0.0;
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("MP %i : adzero done\n", mp->step);
}

/* End of mpads.c */
