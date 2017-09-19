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

#define LBCORE_C

#include <stdio.h>
#include <stdlib.h>
#include "sunlightlb.h"

/* Local function prototypes */

static void bblstat(LBSys);
static void bblstatforce(LBSys);
static void bblmove(LBSys);
static void bblmoveforce(LBSys);

/* Two generic calls to bounce-back routines above */

void bbl(LBSys lb) {
  if (lb->changed_objects) syncobjs(lb);
  if (lb->objects_moving) bblmove(lb);
  else bblstat(lb); 
}

void bblforce(LBSys lb) {
  if (lb->changed_objects) syncobjs(lb);
  if (lb->objects_moving) bblmoveforce(lb);
  else bblstatforce(lb); 
}

/* Take an LB update step. */

void lbstep(LBSys lb) {
  bbl(lb); propagate(lb); collide(lb); 
}

/* Implement bounce back on the links connecting interior and
 * exterior parts of the fluid for stationary objects.
 * This version is for stationary objects
 */

static void bblstat(LBSys lb) {
  int i, j, k, inode1, inode2, ip, jp, kp, im, jm, km, m;
  int esi, esj, esk;
  int *mask;
  LBReal *aa;
  LBReal a, b;
/* Check for correct place in cycle */
  if (lb->status != POST_COLLIDE) {
    printf("Bounceback apparently out of position\n"); wrstatus(lb, stdout);
  }
  esi = lb->esi; esj = lb->esj; esk = lb->esk;
  aa = lb->a; mask = lb->mask;
/* Run through all boundary links */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode1 = esk*(esj*i+j)+k;
    m = mask[inode1];
    if (m > 0) {
      ip = i + 1; if (ip >= esi) ip = 0;
      jp = j + 1; if (jp >= esj) jp = 0;
      kp = k + 1; if (kp >= esk) kp = 0;
      im = i - 1; if (im < 0) im = esi-1;
      jm = j - 1; if (jm < 0) jm = esj-1;
      km = k - 1; if (km < 0) km = esk-1;
/* Check for boundary link c_01 = ( 1, 0, 0) */
      if (m & 0x0001) {
	inode2 = esk*(esj*ip+j)+k;
        a = aa[ 1+15*inode1]; b = aa[ 2+15*inode2];
        aa[ 1+15*inode1] = b; aa[ 2+15*inode2] = a;
      }
/* Check for boundary link c_02 = (-1, 0, 0) */
      if (m & 0x0002) {
	inode2 = esk*(esj*im+j)+k;
        a = aa[ 2+15*inode1]; b = aa[ 1+15*inode2];
        aa[ 2+15*inode1] = b; aa[ 1+15*inode2] = a;
      }
/* Check for boundary link c_03 = ( 0, 1, 0) */
      if (m & 0x0004) {
	inode2 = esk*(esj*i+jp)+k;
        a = aa[ 3+15*inode1]; b = aa[ 4+15*inode2];
        aa[ 3+15*inode1] = b; aa[ 4+15*inode2] = a;
      }
/* Check for boundary link c_04 = ( 0,-1, 0) */
      if (m & 0x0008) {
	inode2 = esk*(esj*i+jm)+k;
        a = aa[ 4+15*inode1]; b = aa[ 3+15*inode2];
        aa[ 4+15*inode1] = b; aa[ 3+15*inode2] = a;
      }
/* Check for boundary link c_05 = ( 0, 0, 1) */
      if (m & 0x0010) {
	inode2 = esk*(esj*i+j)+kp;
        a = aa[ 5+15*inode1]; b = aa[ 6+15*inode2];
        aa[ 5+15*inode1] = b; aa[ 6+15*inode2] = a;
      }
/* Check for boundary link c_06 = ( 0, 0,-1) */
      if (m & 0x0020) {
	inode2 = esk*(esj*i+j)+km;
        a = aa[ 6+15*inode1]; b = aa[ 5+15*inode2];
        aa[ 6+15*inode1] = b; aa[ 5+15*inode2] = a;
      }
/* Check for boundary link c_07 = ( 1, 1, 1) */
      if (m & 0x0040) {
	inode2 = esk*(esj*ip+jp)+kp;
        a = aa[ 7+15*inode1]; b = aa[14+15*inode2];
        aa[ 7+15*inode1] = b; aa[14+15*inode2] = a;
      }
/* Check for boundary link c_08 = (-1, 1, 1) */
      if (m & 0x0080) {
	inode2 = esk*(esj*im+jp)+kp;
        a = aa[ 8+15*inode1]; b = aa[13+15*inode2];
        aa[ 8+15*inode1] = b; aa[13+15*inode2] = a;
      }
/* Check for boundary link c_09 = ( 1,-1, 1) */
      if (m & 0x0100) {
	inode2 = esk*(esj*ip+jm)+kp;
        a = aa[ 9+15*inode1]; b = aa[12+15*inode2];
        aa[ 9+15*inode1] = b; aa[12+15*inode2] = a;
      }
/* Check for boundary link c_10 = (-1,-1, 1) */
      if (m & 0x0200) {
	inode2 = esk*(esj*im+jm)+kp;
        a = aa[10+15*inode1]; b = aa[11+15*inode2];
        aa[10+15*inode1] = b; aa[11+15*inode2] = a;
      }
/* Check for boundary link c_11 = ( 1, 1,-1) */
      if (m & 0x0400) {
	inode2 = esk*(esj*ip+jp)+km;
        a = aa[11+15*inode1]; b = aa[10+15*inode2];
        aa[11+15*inode1] = b; aa[10+15*inode2] = a;
      }
/* Check for boundary link c_12 = (-1, 1,-1) */
      if (m & 0x0800) {
	inode2 = esk*(esj*im+jp)+km;
        a = aa[12+15*inode1]; b = aa[ 9+15*inode2];
        aa[12+15*inode1] = b; aa[ 9+15*inode2] = a;
      }
/* Check for boundary link c_13 = ( 1,-1,-1) */
      if (m & 0x1000) {
	inode2 = esk*(esj*ip+jm)+km;
        a = aa[13+15*inode1]; b = aa[ 8+15*inode2];
        aa[13+15*inode1] = b; aa[ 8+15*inode2] = a;
      }
/* Check for boundary link c_14 = (-1,-1,-1) */
      if (m & 0x2000) {
	inode2 = esk*(esj*im+jm)+km;
        a = aa[14+15*inode1]; b = aa[ 7+15*inode2];
        aa[14+15*inode1] = b; aa[ 7+15*inode2] = a;
      }
    }
  }
/* The status needs updating */
  lb->status = POST_BBL;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : bbl (STAT) done\n", lb->step);
}


/* Implement bounce back on the links connecting interior and
 * exterior parts of the fluid for stationary objects, 
 * collecting forces and torques.
 * This version is for stationary objects
 */

static void bblstatforce(LBSys lb) {
  int i, j, k, inode1, inode2, ip, jp, kp, im, jm, km, m;
  int esi, esj, esk;
  int *mask;
  LBReal *aa;
  LBReal *fx, *fy, *fz;
  LBReal fac, a, b;
/* Check for correct place in cycle */
  if (lb->status != POST_COLLIDE) {
    printf("Bounceback apparently out of position\n"); wrstatus(lb, stdout);
  }
  esi = lb->esi; esj = lb->esj; esk = lb->esk;
  fx = lb->fx; fy = lb->fy; fz = lb->fz;
  aa = lb->a; mask = lb->mask;
/* Run through all boundary links */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode1 = esk*(esj*i+j)+k;
    m = mask[inode1];
    if (m > 0) {
      ip = i + 1; if (ip >= esi) ip = 0;
      jp = j + 1; if (jp >= esj) jp = 0;
      kp = k + 1; if (kp >= esk) kp = 0;
      im = i - 1; if (im < 0) im = esi-1;
      jm = j - 1; if (jm < 0) jm = esj-1;
      km = k - 1; if (km < 0) km = esk-1;
/* Check for boundary link c_01 = ( 1, 0, 0) */
      if (m & 0x0001) {
	inode2 = esk*(esj*ip+j)+k;
        a = aa[ 1+15*inode1]; b = aa[ 2+15*inode2];
        aa[ 1+15*inode1] = b; aa[ 2+15*inode2] = a;
        fac = 2.0*(a - b);
        fx[inode1] += fac;
      }
/* Check for boundary link c_02 = (-1, 0, 0) */
      if (m & 0x0002) {
	inode2 = esk*(esj*im+j)+k;
        a = aa[ 2+15*inode1]; b = aa[ 1+15*inode2];
        aa[ 2+15*inode1] = b; aa[ 1+15*inode2] = a;
        fac = 2.0*(a - b);
        fx[inode1] -= fac;
      }
/* Check for boundary link c_03 = ( 0, 1, 0) */
      if (m & 0x0004) {
	inode2 = esk*(esj*i+jp)+k;
        a = aa[ 3+15*inode1]; b = aa[ 4+15*inode2];
        aa[ 3+15*inode1] = b; aa[ 4+15*inode2] = a;
        fac = 2.0*(a - b);
        fy[inode1] += fac;
      }
/* Check for boundary link c_04 = ( 0,-1, 0) */
      if (m & 0x0008) {
	inode2 = esk*(esj*i+jm)+k;
        a = aa[ 4+15*inode1]; b = aa[ 3+15*inode2];
        aa[ 4+15*inode1] = b; aa[ 3+15*inode2] = a;
        fac = 2.0*(a - b);
        fy[inode1] -= fac;
      }
/* Check for boundary link c_05 = ( 0, 0, 1) */
      if (m & 0x0010) {
	inode2 = esk*(esj*i+j)+kp;
        a = aa[ 5+15*inode1]; b = aa[ 6+15*inode2];
        aa[ 5+15*inode1] = b; aa[ 6+15*inode2] = a;
        fac = 2.0*(a - b);
        fz[inode1] += fac;
      }
/* Check for boundary link c_06 = ( 0, 0,-1) */
      if (m & 0x0020) {
	inode2 = esk*(esj*i+j)+km;
        a = aa[ 6+15*inode1]; b = aa[ 5+15*inode2];
        aa[ 6+15*inode1] = b; aa[ 5+15*inode2] = a;
        fac = 2.0*(a - b);
        fz[inode1] -= fac;
      }
/* Check for boundary link c_07 = ( 1, 1, 1) */
      if (m & 0x0040) {
	inode2 = esk*(esj*ip+jp)+kp;
        a = aa[ 7+15*inode1]; b = aa[14+15*inode2];
        aa[ 7+15*inode1] = b; aa[14+15*inode2] = a;
        fac = 2.0*(a - b);
        fx[inode1] += fac; fy[inode1] += fac; fz[inode1] += fac;
      }
/* Check for boundary link c_08 = (-1, 1, 1) */
      if (m & 0x0080) {
	inode2 = esk*(esj*im+jp)+kp;
        a = aa[ 8+15*inode1]; b = aa[13+15*inode2];
        aa[ 8+15*inode1] = b; aa[13+15*inode2] = a;
        fac = 2.0*(a - b);
        fx[inode1] -= fac; fy[inode1] += fac; fz[inode1] += fac;
      }
/* Check for boundary link c_09 = ( 1,-1, 1) */
      if (m & 0x0100) {
	inode2 = esk*(esj*ip+jm)+kp;
        a = aa[ 9+15*inode1]; b = aa[12+15*inode2];
        aa[ 9+15*inode1] = b; aa[12+15*inode2] = a;
        fac = 2.0*(a - b);
        fx[inode1] += fac; fy[inode1] -= fac; fz[inode1] += fac;
      }
/* Check for boundary link c_10 = (-1,-1, 1) */
      if (m & 0x0200) {
	inode2 = esk*(esj*im+jm)+kp;
        a = aa[10+15*inode1]; b = aa[11+15*inode2];
        aa[10+15*inode1] = b; aa[11+15*inode2] = a;
        fac = 2.0*(a - b);
        fx[inode1] -= fac; fy[inode1] -= fac; fz[inode1] += fac;
      }
/* Check for boundary link c_11 = ( 1, 1,-1) */
      if (m & 0x0400) {
	inode2 = esk*(esj*ip+jp)+km;
        a = aa[11+15*inode1]; b = aa[10+15*inode2];
        aa[11+15*inode1] = b; aa[10+15*inode2] = a;
        fac = 2.0*(a - b);
        fx[inode1] += fac; fy[inode1] += fac; fz[inode1] -= fac;
      }
/* Check for boundary link c_12 = (-1, 1,-1) */
      if (m & 0x0800) {
	inode2 = esk*(esj*im+jp)+km;
        a = aa[12+15*inode1]; b = aa[ 9+15*inode2];
        aa[12+15*inode1] = b; aa[ 9+15*inode2] = a;
        fac = 2.0*(a - b);
        fx[inode1] -= fac; fy[inode1] += fac; fz[inode1] -= fac;
      }
/* Check for boundary link c_13 = ( 1,-1,-1) */
      if (m & 0x1000) {
	inode2 = esk*(esj*ip+jm)+km;
        a = aa[13+15*inode1]; b = aa[ 8+15*inode2];
        aa[13+15*inode1] = b; aa[ 8+15*inode2] = a;
        fac = 2.0*(a - b);
        fx[inode1] += fac; fy[inode1] -= fac; fz[inode1] -= fac;
      }
/* Check for boundary link c_14 = (-1,-1,-1) */
      if (m & 0x2000) {
	inode2 = esk*(esj*im+jm)+km;
        a = aa[14+15*inode1]; b = aa[ 7+15*inode2];
        aa[14+15*inode1] = b; aa[ 7+15*inode2] = a;
        fac = 2.0*(a - b);
        fx[inode1] -= fac; fy[inode1] -= fac; fz[inode1] -= fac;
      }
    }
  }
/* The status needs updating */
  lb->status = POST_BBL;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : bblforce (STAT) done\n", lb->step);
}


/* Implement bounce back on the links connecting interior and
 * exterior parts of the fluid for moving objects.
 * This version is for moving objects
 * In this version, the mean density rho0 is used.
 */

static void bblmove(LBSys lb) {
  int i, j, k, inode1, inode2, ip, jp, kp, im, jm, km, m;
  int esi, esj, esk;
  int *mask;
  LBReal *aa;
  LBReal *ux, *uy, *uz;
  LBReal a, b, fac;
/* Check for correct place in cycle */
  if (lb->status != POST_COLLIDE) {
    printf("Bounceback apparently out of position\n"); wrstatus(lb, stdout);
  }
  esi = lb->esi; esj = lb->esj; esk = lb->esk;
  ux = lb->ux; uy = lb->uy; uz = lb->uz;
  aa = lb->a; mask = lb->mask;
/* Run through all boundary links */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode1 = esk*(esj*i+j)+k;
    m = mask[inode1];
    if (m > 0) {
      ip = i + 1; if (ip >= esi) ip = 0;
      jp = j + 1; if (jp >= esj) jp = 0;
      kp = k + 1; if (kp >= esk) kp = 0;
      im = i - 1; if (im < 0) im = esi-1;
      jm = j - 1; if (jm < 0) jm = esj-1;
      km = k - 1; if (km < 0) km = esk-1;
/* The velocity components are precalculated but when added to the */
/* a_i, also need to include the weight w_i */
/* Check for boundary link c_01 = ( 1, 0, 0) */
      if (m & 0x0001) {
	inode2 = esk*(esj*ip+j)+k;
        a = aa[ 1+15*inode1]; b = aa[ 2+15*inode2]; fac = 8.0 * ux[inode1];
        aa[ 2+15*inode2] = a - fac; aa[ 1+15*inode1] = b + fac;
      }
/* Check for boundary link c_02 = (-1, 0, 0) */
      if (m & 0x0002) {
	inode2 = esk*(esj*im+j)+k;
        a = aa[ 2+15*inode1]; b = aa[ 1+15*inode2]; fac = - 8.0 * ux[inode1];
        aa[ 1+15*inode2] = a - fac; aa[ 2+15*inode1] = b + fac;
      }
/* Check for boundary link c_03 = ( 0, 1, 0) */
      if (m & 0x0004) {
	inode2 = esk*(esj*i+jp)+k;
        a = aa[ 3+15*inode1]; b = aa[ 4+15*inode2]; fac = 8.0 * uy[inode1];
        aa[ 4+15*inode2] = a - fac; aa[ 3+15*inode1] = b + fac;
      }
/* Check for boundary link c_04 = ( 0,-1, 0) */
      if (m & 0x0008) {
	inode2 = esk*(esj*i+jm)+k;
        a = aa[ 4+15*inode1]; b = aa[ 3+15*inode2]; fac = - 8.0 * uy[inode1];
        aa[ 3+15*inode2] = a - fac; aa[ 4+15*inode1] = b + fac;
      }
/* Check for boundary link c_05 = ( 0, 0, 1) */
      if (m & 0x0010) {
	inode2 = esk*(esj*i+j)+kp;
        a = aa[ 5+15*inode1]; b = aa[ 6+15*inode2]; fac = 8.0 * uz[inode1];
        aa[ 6+15*inode2] = a - fac; aa[ 5+15*inode1] = b + fac;
      }
/* Check for boundary link c_06 = ( 0, 0,-1) */
      if (m & 0x0020) {
	inode2 = esk*(esj*i+j)+km;
        a = aa[ 6+15*inode1]; b = aa[ 5+15*inode2]; fac = - 8.0 * uz[inode1];
        aa[ 5+15*inode2] = a - fac; aa[ 6+15*inode1] = b + fac;
      }
/* Check for boundary link c_07 = ( 1, 1, 1) */
      if (m & 0x0040) {
	inode2 = esk*(esj*ip+jp)+kp;
        a = aa[ 7+15*inode1]; b = aa[14+15*inode2];
	fac = ux[inode1] + uy[inode1] + uz[inode1];
        aa[14+15*inode2] = a - fac; aa[ 7+15*inode1] = b + fac;
      }
/* Check for boundary link c_08 = (-1, 1, 1) */
      if (m & 0x0080) {
	inode2 = esk*(esj*im+jp)+kp;
        a = aa[ 8+15*inode1]; b = aa[13+15*inode2];
	fac = - ux[inode1] + uy[inode1] + uz[inode1];
        aa[13+15*inode2] = a - fac; aa[ 8+15*inode1] = b + fac;
      }
/* Check for boundary link c_09 = ( 1,-1, 1) */
      if (m & 0x0100) {
	inode2 = esk*(esj*ip+jm)+kp;
        a = aa[ 9+15*inode1]; b = aa[12+15*inode2];
	fac = ux[inode1] - uy[inode1] + uz[inode1];
        aa[12+15*inode2] = a - fac; aa[ 9+15*inode1] = b + fac;
      }
/* Check for boundary link c_10 = (-1,-1, 1) */
      if (m & 0x0200) {
	inode2 = esk*(esj*im+jm)+kp;
        a = aa[10+15*inode1]; b = aa[11+15*inode2];
	fac = - ux[inode1] - uy[inode1] + uz[inode1];
        aa[11+15*inode2] = a - fac; aa[10+15*inode1] = b + fac;
      }
/* Check for boundary link c_11 = ( 1, 1,-1) */
      if (m & 0x0400) {
	inode2 = esk*(esj*ip+jp)+km;
        a = aa[11+15*inode1]; b = aa[10+15*inode2];
	fac = ux[inode1] + uy[inode1] - uz[inode1];
        aa[10+15*inode2] = a - fac; aa[11+15*inode1] = b + fac;
      }
/* Check for boundary link c_12 = (-1, 1,-1) */
      if (m & 0x0800) {
	inode2 = esk*(esj*im+jp)+km;
        a = aa[12+15*inode1]; b = aa[ 9+15*inode2];
	fac = - ux[inode1] + uy[inode1] - uz[inode1];
        aa[ 9+15*inode2] = a - fac; aa[12+15*inode1] = b + fac;
      }
/* Check for boundary link c_13 = ( 1,-1,-1) */
      if (m & 0x1000) {
	inode2 = esk*(esj*ip+jm)+km;
        a = aa[13+15*inode1]; b = aa[ 8+15*inode2];
	fac = ux[inode1] - uy[inode1] - uz[inode1];
        aa[ 8+15*inode2] = a - fac; aa[13+15*inode1] = b + fac;
      }
/* Check for boundary link c_14 = (-1,-1,-1) */
      if (m & 0x2000) {
	inode2 = esk*(esj*im+jm)+km;
        a = aa[14+15*inode1]; b = aa[ 7+15*inode2];
	fac = - ux[inode1] - uy[inode1] - uz[inode1];
        aa[ 7+15*inode2] = a - fac; aa[14+15*inode1] = b + fac;
      }
    }
  }
/* The status needs updating */
  lb->status = POST_BBL;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : bbl (MOV) done\n", lb->step);
}


/* Implement bounce back on the links connecting interior and
 * exterior parts of the fluid for moving objects.
 * Measure the forces and torques thus generated.
 * This version is for moving objects
 * In this version, the mean density is used.
 */

static void bblmoveforce(LBSys lb) {
  int i, j, k, inode1, inode2, ip, jp, kp, im, jm, km, m;
  int esi, esj, esk;
  int *mask;
  LBReal *aa;
  LBReal *fx, *fy, *fz;
  LBReal *ux, *uy, *uz;
  LBReal a, b, fac;
/* Check for correct place in cycle */
  if (lb->status != POST_COLLIDE) {
    printf("Bounceback apparently out of position\n"); wrstatus(lb, stdout);
  }
  esi = lb->esi; esj = lb->esj; esk = lb->esk;
  fx = lb->fx; fy = lb->fy; fz = lb->fz;
  ux = lb->ux; uy = lb->uy; uz = lb->uz;
  aa = lb->a; mask = lb->mask;
/* Run through all boundary links */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode1 = esk*(esj*i+j)+k;
    m = mask[inode1];
    if (m > 0) {
      ip = i + 1; if (ip >= esi) ip = 0;
      jp = j + 1; if (jp >= esj) jp = 0;
      kp = k + 1; if (kp >= esk) kp = 0;
      im = i - 1; if (im < 0) im = esi-1;
      jm = j - 1; if (jm < 0) jm = esj-1;
      km = k - 1; if (km < 0) km = esk-1;
/* The velocity components are precalculated but when added to the */
/* a_i, also need to include the weight w_i */
/* Check for boundary link c_01 = ( 1, 0, 0) */
      if (m & 0x0001) {
	inode2 = esk*(esj*ip+j)+k;
        a = aa[ 1+15*inode1]; b = aa[ 2+15*inode2]; fac = 8.0 * ux[inode1];
        aa[ 2+15*inode2] = a - fac; aa[ 1+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fx[inode1] += fac;
      }
/* Check for boundary link c_02 = (-1, 0, 0) */
      if (m & 0x0002) {
	inode2 = esk*(esj*im+j)+k;
        a = aa[ 2+15*inode1]; b = aa[ 1+15*inode2]; fac = - 8.0 * ux[inode1];
        aa[ 1+15*inode2] = a - fac; aa[ 2+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fx[inode1] -= fac;
      }
/* Check for boundary link c_03 = ( 0, 1, 0) */
      if (m & 0x0004) {
	inode2 = esk*(esj*i+jp)+k;
        a = aa[ 3+15*inode1]; b = aa[ 4+15*inode2]; fac = 8.0 * uy[inode1];
        aa[ 4+15*inode2] = a - fac; aa[ 3+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fy[inode1] += fac;
      }
/* Check for boundary link c_04 = ( 0,-1, 0) */
      if (m & 0x0008) {
	inode2 = esk*(esj*i+jm)+k;
        a = aa[ 4+15*inode1]; b = aa[ 3+15*inode2]; fac = - 8.0 * uy[inode1];
        aa[ 3+15*inode2] = a - fac; aa[ 4+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fy[inode1] -= fac;
      }
/* Check for boundary link c_05 = ( 0, 0, 1) */
      if (m & 0x0010) {
	inode2 = esk*(esj*i+j)+kp;
        a = aa[ 5+15*inode1]; b = aa[ 6+15*inode2]; fac = 8.0 * uz[inode1];
        aa[ 6+15*inode2] = a - fac; aa[ 5+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fz[inode1] += fac;
      }
/* Check for boundary link c_06 = ( 0, 0,-1) */
      if (m & 0x0020) {
	inode2 = esk*(esj*i+j)+km;
        a = aa[ 6+15*inode1]; b = aa[ 5+15*inode2]; fac = - 8.0 * uz[inode1];
        aa[ 5+15*inode2] = a - fac; aa[ 6+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fz[inode1] -= fac;
      }
/* Check for boundary link c_07 = ( 1, 1, 1) */
      if (m & 0x0040) {
	inode2 = esk*(esj*ip+jp)+kp;
        a = aa[ 7+15*inode1]; b = aa[14+15*inode2];
	fac = ux[inode1] + uy[inode1] + uz[inode1];
        aa[14+15*inode2] = a - fac;  aa[ 7+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fx[inode1] += fac; fy[inode1] += fac; fz[inode1] += fac;
      }
/* Check for boundary link c_08 = (-1, 1, 1) */
      if (m & 0x0080) {
	inode2 = esk*(esj*im+jp)+kp;
        a = aa[ 8+15*inode1]; b = aa[13+15*inode2];
	fac = - ux[inode1] + uy[inode1] + uz[inode1];
        aa[13+15*inode2] = a - fac; aa[ 8+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fx[inode1] -= fac; fy[inode1] += fac; fz[inode1] += fac;
      }
/* Check for boundary link c_09 = ( 1,-1, 1) */
      if (m & 0x0100) {
	inode2 = esk*(esj*ip+jm)+kp;
        a = aa[ 9+15*inode1]; b = aa[12+15*inode2];
	fac = ux[inode1] - uy[inode1] + uz[inode1];
        aa[12+15*inode2] = a - fac; aa[ 9+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fx[inode1] += fac; fy[inode1] -= fac; fz[inode1] += fac;
      }
/* Check for boundary link c_10 = (-1,-1, 1) */
      if (m & 0x0200) {
	inode2 = esk*(esj*im+jm)+kp;
        a = aa[10+15*inode1]; b = aa[11+15*inode2];
	fac = - ux[inode1] - uy[inode1] + uz[inode1];
        aa[11+15*inode2] = a - fac; aa[10+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fx[inode1] -= fac; fy[inode1] -= fac; fz[inode1] += fac;
      }
/* Check for boundary link c_11 = ( 1, 1,-1) */
      if (m & 0x0400) {
	inode2 = esk*(esj*ip+jp)+km;
        a = aa[11+15*inode1]; b = aa[10+15*inode2];
	fac = ux[inode1] + uy[inode1] - uz[inode1];
        aa[10+15*inode2] = a - fac; aa[11+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fx[inode1] += fac; fy[inode1] += fac; fz[inode1] -= fac;
      }
/* Check for boundary link c_12 = (-1, 1,-1) */
      if (m & 0x0800) {
	inode2 = esk*(esj*im+jp)+km;
        a = aa[12+15*inode1]; b = aa[ 9+15*inode2];
	fac = - ux[inode1] + uy[inode1] - uz[inode1];
        aa[ 9+15*inode2] = a - fac; aa[12+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fx[inode1] -= fac; fy[inode1] += fac; fz[inode1] -= fac;
      }
/* Check for boundary link c_13 = ( 1,-1,-1) */
      if (m & 0x1000) {
	inode2 = esk*(esj*ip+jm)+km;
        a = aa[13+15*inode1]; b = aa[ 8+15*inode2];
	fac = ux[inode1] - uy[inode1] - uz[inode1];
        aa[ 8+15*inode2] = a - fac; aa[13+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fx[inode1] += fac; fy[inode1] -= fac; fz[inode1] -= fac;
      }
/* Check for boundary link c_14 = (-1,-1,-1) */
      if (m & 0x2000) {
	inode2 = esk*(esj*im+jm)+km;
        a = aa[14+15*inode1]; b = aa[ 7+15*inode2];
	fac = - ux[inode1] - uy[inode1] - uz[inode1];
        aa[ 7+15*inode2] = a - fac; aa[14+15*inode1] = b + fac;
        fac = 2.0*(a - b - fac);
        fx[inode1] -= fac; fy[inode1] -= fac; fz[inode1] -= fac;
      }
    }
  }
/* The status needs updating */
  lb->status = POST_BBL;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : bblforce (MOV) done\n", lb->step);
}


/* Propogate occupation densities */

void propagate(LBSys lb) { 
  int i, j, k, inode1, inode2;
  int esi, esj, esk;
  LBReal *aa;
  LBReal t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14;
  /* Check for correct place in cycle */
  if (lb->status != POST_BBL) {
    printf("Propagation apparently out of position\n"); wrstatus(lb, stdout);
  }
  esi = lb->esi; esj = lb->esj; esk = lb->esk; aa = lb->a;
/* Stream a01, a07, a09, a11, a13 in positive i direction
      and a02, a08, a10, a12, a14 in negative i direction */
 for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode1 = esk*(esj*(esi-1)+j)+k; inode2 = esk*j+k;
    t01 = aa[ 1+15*inode1];
    t07 = aa[ 7+15*inode1];
    t09 = aa[ 9+15*inode1];
    t11 = aa[11+15*inode1];
    t13 = aa[13+15*inode1];
    t02 = aa[ 2+15*inode2];
    t08 = aa[ 8+15*inode2];
    t10 = aa[10+15*inode2];
    t12 = aa[12+15*inode2];
    t14 = aa[14+15*inode2];
    for (i=esi-1; i>0; i--) {
      inode1 = esk*(esj*i+j)+k; inode2 = esk*(esj*(i-1)+j)+k;
      aa[ 1+15*inode1] = aa[ 1+15*inode2];
      aa[ 7+15*inode1] = aa[ 7+15*inode2];
      aa[ 9+15*inode1] = aa[ 9+15*inode2];
      aa[11+15*inode1] = aa[11+15*inode2];
      aa[13+15*inode1] = aa[13+15*inode2];
    }
    for (i=0; i<esi-1; i++) {
      inode1 = esk*(esj*i+j)+k; inode2 = esk*(esj*(i+1)+j)+k;
      aa[ 2+15*inode1] = aa[ 2+15*inode2];
      aa[ 8+15*inode1] = aa[ 8+15*inode2];
      aa[10+15*inode1] = aa[10+15*inode2];
      aa[12+15*inode1] = aa[12+15*inode2];
      aa[14+15*inode1] = aa[14+15*inode2];
    }
    inode1 = esk*j+k; inode2 = esk*(esj*(esi-1)+j)+k;
    aa[ 1+15*inode1] = t01;
    aa[ 7+15*inode1] = t07;
    aa[ 9+15*inode1] = t09;
    aa[11+15*inode1] = t11;
    aa[13+15*inode1] = t13;
    aa[ 2+15*inode2] = t02;
    aa[ 8+15*inode2] = t08;
    aa[10+15*inode2] = t10;
    aa[12+15*inode2] = t12;
    aa[14+15*inode2] = t14;
  }
  /* Stream a03, a07, a08, a11, a12 in positive j direction
      and a04, a09, a10, a13, a14 in negative j direction */
  for (i=0; i<esi; i++) for (k=0; k<esk; k++) {
    inode1 = esk*(esj*i+esj-1)+k; inode2 = esk*esj*i+k;
    t03 = aa[ 3+15*inode1];
    t07 = aa[ 7+15*inode1];
    t08 = aa[ 8+15*inode1];
    t11 = aa[11+15*inode1];
    t12 = aa[12+15*inode1];
    t04 = aa[ 4+15*inode2];
    t09 = aa[ 9+15*inode2];
    t10 = aa[10+15*inode2];
    t13 = aa[13+15*inode2];
    t14 = aa[14+15*inode2];
    for (j=esj-1; j>0; j--) {
      inode1 = esk*(esj*i+j)+k; inode2 = esk*(esj*i+j-1)+k;
      aa[ 3+15*inode1] = aa[ 3+15*inode2];
      aa[ 7+15*inode1] = aa[ 7+15*inode2];
      aa[ 8+15*inode1] = aa[ 8+15*inode2];
      aa[11+15*inode1] = aa[11+15*inode2];
      aa[12+15*inode1] = aa[12+15*inode2];
    }
    for (j=0; j<esj-1; j++) {
      inode1 = esk*(esj*i+j)+k; inode2 = esk*(esj*i+j+1)+k;
      aa[ 4+15*inode1] = aa[ 4+15*inode2];
      aa[ 9+15*inode1] = aa[ 9+15*inode2];
      aa[10+15*inode1] = aa[10+15*inode2];
      aa[13+15*inode1] = aa[13+15*inode2];
      aa[14+15*inode1] = aa[14+15*inode2];
    }
    inode1 = esk*esj*i+k; inode2 = esk*(esj*i+esj-1)+k;
    aa[ 3+15*inode1] = t03;
    aa[ 7+15*inode1] = t07;
    aa[ 8+15*inode1] = t08;
    aa[11+15*inode1] = t11;
    aa[12+15*inode1] = t12;
    aa[ 4+15*inode2] = t04;
    aa[ 9+15*inode2] = t09;
    aa[10+15*inode2] = t10;
    aa[13+15*inode2] = t13;
    aa[14+15*inode2] = t14;
  }
  /* Stream a05, a07, a08, a09, a10 in positive k direction
      and a06, a11, a12, a13, a14 in negative k direction */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) {
    inode1 = esk*(esj*i+j)+esk-1; inode2 = esk*(esj*i+j);
    t05 = aa[ 5+15*inode1];
    t07 = aa[ 7+15*inode1];
    t08 = aa[ 8+15*inode1];
    t09 = aa[ 9+15*inode1];
    t10 = aa[10+15*inode1];
    t06 = aa[ 6+15*inode2];
    t11 = aa[11+15*inode2];
    t12 = aa[12+15*inode2];
    t13 = aa[13+15*inode2];
    t14 = aa[14+15*inode2];
    for (k=esk-1; k>0; k--) {
      inode1 = esk*(esj*i+j)+k; inode2 = esk*(esj*i+j)+k-1;
      aa[ 5+15*inode1] = aa[ 5+15*inode2];
      aa[ 7+15*inode1] = aa[ 7+15*inode2];
      aa[ 8+15*inode1] = aa[ 8+15*inode2];
      aa[ 9+15*inode1] = aa[ 9+15*inode2];
      aa[10+15*inode1] = aa[10+15*inode2];
    }
    for (k=0; k<esk-1; k++) {
      inode1 = esk*(esj*i+j)+k; inode2 = esk*(esj*i+j)+k+1;
      aa[ 6+15*inode1] = aa[ 6+15*inode2];
      aa[11+15*inode1] = aa[11+15*inode2];
      aa[12+15*inode1] = aa[12+15*inode2];
      aa[13+15*inode1] = aa[13+15*inode2];
      aa[14+15*inode1] = aa[14+15*inode2];
    }
    inode1 = esk*(esj*i+j); inode2 = esk*(esj*i+j)+esk-1;
    aa[ 5+15*inode1] = t05;
    aa[ 7+15*inode1] = t07;
    aa[ 8+15*inode1] = t08;
    aa[ 9+15*inode1] = t09;
    aa[10+15*inode1] = t10;
    aa[ 6+15*inode2] = t06;
    aa[11+15*inode2] = t11;
    aa[12+15*inode2] = t12;
    aa[13+15*inode2] = t13;
    aa[14+15*inode2] = t14;
  }
  /* The status needs updating */
  lb->status = POST_PROPAGATE;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : propagate done\n", lb->step);
}


/* Do node collisions everywhere unless skip_solid is defined.
 * Also use stokes to determine appropriate equilibrium stress.
 */

void collide(LBSys lb) {
  int inode;
  int skip_solid, stokes_flow;
  LBReal *aa;
  LBReal *bx, *by, *bz;
  LBReal drho, rux, ruy, ruz, ruxyz, sxx, sxy, sxz, syy, syz, szz, sxyz, trsby3;
  LBReal t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14;
  LBReal t1p2, t3p4, t5p6, t0p6;
  LBReal tlppp, tlmpm, tlpmm, tlmmp, thppp, thmpm, thpmm, thmmp, tlhp;
  LBReal twoby3, dbybc2, edbybc2, dd2bybc4, fdd2bybc4;
  LBReal rhobyb, erhobyb, b;
  LBReal comega, omega;
  LBReal obyr, rho0;
  LBObj *obj;
/* Check for correct place in cycle */
  if (lb->status != POST_PROPAGATE) {
    printf("Collision apparently out of position\n"); wrstatus(lb, stdout);
  }
  aa = lb->a; bx = lb->bx; by = lb->by; bz = lb->bz;
  skip_solid = lb->skip_solid;
  stokes_flow = lb->stokes_flow;
  obj = lb->obj;
/* Generate some constant parameters */
  omega = 1.0/(lb->tau); comega = 1.0 - omega;
  twoby3 = 2.0/3.0; b = 72.0; rho0 = lb->rho0;
  dbybc2 = 1.0/24.0; dd2bybc4 = 1.0/8.0;
  fdd2bybc4 = 4.0*dd2bybc4; edbybc2 = 8.0*dbybc2;
/* Run through as many nodes as are required */
  for (inode = 0; inode<(lb->nnodes); inode++) {
    if (skip_solid && (obj[inode] != NULL)) continue;
    t00 = aa[ 0+15*inode]; t01 = aa[ 1+15*inode]; t02 = aa[ 2+15*inode];
    t03 = aa[ 3+15*inode]; t04 = aa[ 4+15*inode]; t05 = aa[ 5+15*inode];
    t06 = aa[ 6+15*inode]; t07 = aa[ 7+15*inode]; t08 = aa[ 8+15*inode];
    t09 = aa[ 9+15*inode]; t10 = aa[10+15*inode]; t11 = aa[11+15*inode];
    t12 = aa[12+15*inode]; t13 = aa[13+15*inode]; t14 = aa[14+15*inode];
/* Calculate some partial sums */
    t1p2 = t01 + t02; t3p4 = t03 + t04; t5p6 = t05 + t06;
    tlppp = t07 + t08 + t09 + t10; thppp = t11 + t12 + t13 + t14;
    tlmpm = t07 - t08 + t09 - t10; thmpm = t11 - t12 + t13 - t14;
    tlpmm = t07 + t08 - t09 - t10; thpmm = t11 + t12 - t13 - t14;
    tlmmp = t07 - t08 - t09 + t10; thmmp = t11 - t12 - t13 + t14;
    t0p6 = t00 + t1p2 + t3p4 + t5p6; tlhp = tlppp + thppp;
/* Calculate density and momentum */
    drho = t0p6 + tlhp;
    rux = t01 - t02 + tlmpm + thmpm;
    ruy = t03 - t04 + tlpmm + thpmm;
    ruz = t05 - t06 + tlppp - thppp;
/* Add body force to drive flow, only to external fluid nodes */
    if (obj[inode] == NULL) {
      rux += bx[inode]; ruy += by[inode]; ruz += bz[inode];
    }
/* Project out stress */
    sxyz = twoby3*drho - t0p6;
    sxx = sxyz + t1p2;    syy = sxyz + t3p4;    szz = sxyz + t5p6;
    sxy = tlmmp + thmmp;  sxz = tlmpm - thmpm;  syz = tlpmm - thpmm;
/* relax stress distribution */
    if (stokes_flow) {
      sxx = comega*sxx; sxy = comega*sxy; sxz = comega*sxz;
      syy = comega*syy; syz = comega*syz; szz = comega*szz;
    } else {
      obyr = omega/(rho0+drho);
      sxx = comega*sxx + obyr*rux*rux;
      sxy = comega*sxy + obyr*rux*ruy;
      sxz = comega*sxz + obyr*rux*ruz;
      syy = comega*syy + obyr*ruy*ruy;
      syz = comega*syz + obyr*ruy*ruz;
      szz = comega*szz + obyr*ruz*ruz;
    }
/* Rebuild occupation probabilities: a_i = w_i (f_i - rho_0 / b) */
    rhobyb = drho/b;  erhobyb = 8.0*rhobyb;
    trsby3 = (sxx+syy+szz)/3.0;
    aa[ 0+15*inode] = 2.0*erhobyb - 2.0*fdd2bybc4 * trsby3;
    ruxyz = edbybc2*rux;  sxyz = fdd2bybc4*(sxx - trsby3);
    aa[ 1+15*inode] = erhobyb + ruxyz + sxyz;
    aa[ 2+15*inode] = erhobyb - ruxyz + sxyz;
    ruxyz = edbybc2*ruy;  sxyz = fdd2bybc4*(syy - trsby3);
    aa[ 3+15*inode] = erhobyb + ruxyz + sxyz;
    aa[ 4+15*inode] = erhobyb - ruxyz + sxyz;
    ruxyz = edbybc2*ruz;  sxyz = fdd2bybc4*(szz - trsby3);
    aa[ 5+15*inode] = erhobyb + ruxyz + sxyz;
    aa[ 6+15*inode] = erhobyb - ruxyz + sxyz;
    ruxyz = dbybc2*(rux + ruy + ruz);
    sxyz = dd2bybc4*(trsby3 + sxy + sxz + syz);
    aa[ 7+15*inode] = rhobyb + ruxyz + sxyz;
    aa[14+15*inode] = rhobyb - ruxyz + sxyz;
    ruxyz = dbybc2*(rux - ruy - ruz);
    sxyz = dd2bybc4*(trsby3 - sxy - sxz + syz);
    aa[ 8+15*inode] = rhobyb - ruxyz + sxyz;
    aa[13+15*inode] = rhobyb + ruxyz + sxyz;
    ruxyz = dbybc2*(rux - ruy + ruz);
    sxyz = dd2bybc4*(trsby3 - sxy + sxz - syz);
    aa[ 9+15*inode] = rhobyb + ruxyz + sxyz;
    aa[12+15*inode] = rhobyb - ruxyz + sxyz;
    ruxyz = dbybc2*(rux + ruy - ruz);
    sxyz = dd2bybc4*(trsby3 + sxy - sxz - syz);
    aa[10+15*inode] = rhobyb - ruxyz + sxyz;
    aa[11+15*inode] = rhobyb + ruxyz + sxyz;
  }
/* Increment LB step counter, since collide is regarded as the */
/* end of a lattice Boltzmann update step */
  lb->step += 1;
/* The status needs updating */
  lb->status = POST_COLLIDE;
/* Lastly, handle a signal if one arrives */
#ifdef SIGNAL
  if (slb_signal_flag) sighandle(lb);
#endif /* SIGNAL */
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : collide done\n", lb->step);
}

/* End of lbcore.c */
