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

#define MPCORE_C

#include <stdio.h>
#include <stdlib.h>
#include "sunlightlb.h"

/* Internal function prototypes */

static void mpstat(MPSys);
static void mpmove(MPSys);

/* One-liner to take a moment propagation step, handle signals */

void mpstep(MPSys mp) {
  if (mp->to_be_initialised) mpinit(mp);
  if (mp->objects_moving) mpmove(mp); else mpstat(mp); 
#ifdef SIGNAL
  if (slb_signal_flag) sighandle(mp->lb);
#endif /* SIGNAL */
  mp->step += 1;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpdo done\n", mp->name, mp->step);
  }
}

/* Perform a moment propagation step for stationary objects. */

void mpstat(MPSys mp) {
  int i, j, k, ip, jp, kp, im, jm, km, esi, esj, esk, m, inode, inode1, inode2;
  int *mask;
  LBReal alpha, rho0, rhot, ptemp, fac, piece, piece8;
  LBReal piece16, onebyb, eightbyb;
  LBReal *a, *p, *newp, *rho;
  LBObj *obj;
/* Calculate some constants appropriate to D3Q15(r=16) */
  onebyb = 1.0/72.0; eightbyb = 8.0/72.0;
/* Check for correct place in LB update cycle, should be post-collision */
  if (mp->lb->status != POST_COLLIDE) {
    printf("Moment propagation apparently out of position\n");
    wrstatus(mp->lb, stdout);
  }
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; mask = mp->lb->mask; a = mp->lb->a;
  newp = mp->newp; p = mp->p; rho = mp->lb->rho;
  rho0 = mp->lb->rho0; alpha = mp->alpha;
/* Clear array that stores new scalar density */
  for (inode = 0; inode<(mp->lb->nnodes); inode++) newp[inode] = 0.0;
/* Run through all sites on the lattice */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode1 = esk*(esj*i+j)+k;
    ip = i + 1; if (ip >= esi) ip = 0;
    jp = j + 1; if (jp >= esj) jp = 0;
    kp = k + 1; if (kp >= esk) kp = 0;
    im = i - 1; if (im < 0) im = esi-1;
    jm = j - 1; if (jm < 0) jm = esj-1;
    km = k - 1; if (km < 0) km = esk-1;
/* Check for a fluid node, move stuff in appropriate directions */
/* Stuff that propagates along boundary links is taken care of below */
    if (obj[inode1] == NULL) {
      rhot = rho[inode1]; ptemp = p[inode1];
      fac = ptemp/rhot; piece = ptemp*onebyb*(rho0/rhot - alpha);
      piece8 = 8.0*piece; piece16 = 16.0*piece;
      newp[esk*(esj*i+j)+k] += fac*a[15*inode1] + piece16 + alpha*ptemp;
      newp[esk*(esj*ip+j)+k] += fac*a[1+15*inode1] + piece8;
      newp[esk*(esj*im+j)+k] += fac*a[2+15*inode1] + piece8;
      newp[esk*(esj*i+jp)+k] += fac*a[3+15*inode1] + piece8;
      newp[esk*(esj*i+jm)+k] += fac*a[4+15*inode1] + piece8;
      newp[esk*(esj*i+j)+kp] += fac*a[5+15*inode1] + piece8;
      newp[esk*(esj*i+j)+km] += fac*a[6+15*inode1] + piece8;
      newp[esk*(esj*ip+jp)+kp] += fac*a[7+15*inode1] + piece;
      newp[esk*(esj*im+jp)+kp] += fac*a[8+15*inode1] + piece;
      newp[esk*(esj*ip+jm)+kp] += fac*a[9+15*inode1] + piece;
      newp[esk*(esj*im+jm)+kp] += fac*a[10+15*inode1] + piece;
      newp[esk*(esj*ip+jp)+km] += fac*a[11+15*inode1] + piece;
      newp[esk*(esj*im+jp)+km] += fac*a[12+15*inode1] + piece;
      newp[esk*(esj*ip+jm)+km] += fac*a[13+15*inode1] + piece;
      newp[esk*(esj*im+jm)+km] += fac*a[14+15*inode1] + piece;
    } else if ((m = mask[inode1]) && !(m & 0x4000)) {
/* Interior node with boundary links found */
/* For stationary objects, account for these by bounce back into the */
/* appropriate fluid node. */
/* Implement adsorption here by simple adsorbing the total amount */
      if (m & 0x0001) { /* Check for boundary link c_01 = ( 1, 0, 0) */
	inode2 = esk*(esj*ip+j)+k;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[2+15*inode2] + piece;
      }
      if (m & 0x0002) { /* Check for boundary link c_02 = (-1, 0, 0) */
	inode2 = esk*(esj*im+j)+k;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[1+15*inode2] + piece;
      }
      if (m & 0x0004) { /* Check for boundary link c_03 = ( 0, 1, 0) */
	inode2 = esk*(esj*i+jp)+k;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[4+15*inode2] + piece;
      }
      if (m & 0x0008) { /* Check for boundary link c_04 = ( 0,-1, 0) */
	inode2 = esk*(esj*i+jm)+k;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[3+15*inode2] + piece;
      }
      if (m & 0x0010) { /* Check for boundary link c_05 = ( 0, 0, 1) */
	inode2 = esk*(esj*i+j)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[6+15*inode2] + piece;
      }
      if (m & 0x0020) { /* Check for boundary link c_06 = ( 0, 0,-1) */
	inode2 = esk*(esj*i+j)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[5+15*inode2] + piece;
      }
      if (m & 0x0040) { /* Check for boundary link c_07 = ( 1, 1, 1) */
	inode2 = esk*(esj*ip+jp)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[14+15*inode2] + piece;
      }
      if (m & 0x0080) { /* Check for boundary link c_08 = (-1, 1, 1) */
	inode2 = esk*(esj*im+jp)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[13+15*inode2] + piece;
      }
      if (m & 0x0100) { /* Check for boundary link c_09 = ( 1,-1, 1) */
	inode2 = esk*(esj*ip+jm)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[12+15*inode2] + piece;
      }
      if (m & 0x0200) { /* Check for boundary link c_10 = (-1,-1, 1) */
	inode2 = esk*(esj*im+jm)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[11+15*inode2] + piece;
      }
      if (m & 0x0400) { /* Check for boundary link c_11 = ( 1, 1,-1) */
	inode2 = esk*(esj*ip+jp)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[10+15*inode2] + piece;
      }
      if (m & 0x0800) { /* Check for boundary link c_12 = (-1, 1,-1) */
	inode2 = esk*(esj*im+jp)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[9+15*inode2] + piece;
      }
      if (m & 0x1000) { /* Check for boundary link c_13 = ( 1,-1,-1) */
	inode2 = esk*(esj*ip+jm)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[8+15*inode2] + piece;
      }
      if (m & 0x2000) { /* Check for boundary link c_14 = (-1,-1,-1) */
	inode2 = esk*(esj*im+jm)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[7+15*inode2] + piece;
      }
    }
  }
/* Update the scalar density in fluid sites, use p to store */
/* the adsorbed scalar density on adsorption sites */
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    ptemp = newp[inode];
    if (obj[inode] == NULL) p[inode] = ptemp;
    else {
      if (mask[inode] & 0x4000) p[inode] += ptemp;
      else p[inode] = 0.0;
    }
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpstat done\n", mp->name, mp->step);
  }
}


/* Perform a moment propagation step, for moving objects. */

void mpmove(MPSys mp) {
  int i, j, k, ip, jp, kp, im, jm, km, m, esi, esj, esk, inode, inode1, inode2;
  int *mask;
  LBReal alpha, rho0, rhot, ptemp, fac, piece, piece8;
  LBReal piece16, onebyb, eightbyb;
  LBReal rx, ry, rz, ux, uy, uz;
  LBReal *a, *p, *newp, *rho;
  LBObj *obj, lp;
/* Calculate some constants appropriate to D3Q15(r=16) */
  onebyb = 1.0/72.0; eightbyb = 8.0/72.0;
/* Check for correct place in LB update cycle, should be post-collision */
  if (mp->lb->status != POST_COLLIDE) {
    printf("Moment propagation apparently out of position\n");
    wrstatus(mp->lb, stdout);
  }
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; mask = mp->lb->mask; a = mp->lb->a;
  newp = mp->newp; p = mp->p; rho = mp->lb->rho;
  rho0 = mp->lb->rho0; alpha = mp->alpha;
/* Clear array that stores new scalar density */
  for (inode = 0; inode<(mp->lb->nnodes); inode++) newp[inode] = 0.0;
/* Run through all sites on the lattice */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode1 = esk*(esj*i+j)+k;
    ip = i + 1; if (ip >= esi) ip = 0;
    jp = j + 1; if (jp >= esj) jp = 0;
    kp = k + 1; if (kp >= esk) kp = 0;
    im = i - 1; if (im < 0) im = esi-1;
    jm = j - 1; if (jm < 0) jm = esj-1;
    km = k - 1; if (km < 0) km = esk-1;
/* Check for a fluid node, move stuff in appropriate directions */
/* Stuff that propagates along boundary links is taken care of below */
    if (obj[inode1] == NULL) {
      rhot = rho[inode1]; ptemp = p[inode1];
      fac = ptemp/rhot; piece = ptemp*onebyb*(rho0/rhot - alpha);
      piece8 = 8.0*piece; piece16 = 16.0*piece;
      newp[esk*(esj*i+j)+k] += fac*a[15*inode1] + piece16 + alpha*ptemp;
      newp[esk*(esj*ip+j)+k] += fac*a[1+15*inode1] + piece8;
      newp[esk*(esj*im+j)+k] += fac*a[2+15*inode1] + piece8;
      newp[esk*(esj*i+jp)+k] += fac*a[3+15*inode1] + piece8;
      newp[esk*(esj*i+jm)+k] += fac*a[4+15*inode1] + piece8;
      newp[esk*(esj*i+j)+kp] += fac*a[5+15*inode1] + piece8;
      newp[esk*(esj*i+j)+km] += fac*a[6+15*inode1] + piece8;
      newp[esk*(esj*ip+jp)+kp] += fac*a[7+15*inode1] + piece;
      newp[esk*(esj*im+jp)+kp] += fac*a[8+15*inode1] + piece;
      newp[esk*(esj*ip+jm)+kp] += fac*a[9+15*inode1] + piece;
      newp[esk*(esj*im+jm)+kp] += fac*a[10+15*inode1] + piece;
      newp[esk*(esj*ip+jp)+km] += fac*a[11+15*inode1] + piece;
      newp[esk*(esj*im+jp)+km] += fac*a[12+15*inode1] + piece;
      newp[esk*(esj*ip+jm)+km] += fac*a[13+15*inode1] + piece;
      newp[esk*(esj*im+jm)+km] += fac*a[14+15*inode1] + piece;
    } else if ((m = mask[inode1]) && !(m & 0x4000)) {
/* Interior node with boundary links found */
/* Calculate velocity components and factor in 2D rho0 / b c^2. */
/* Note that the offset 1/2 c_i can be omitted from the angular part. */
/* When added to the a_i, also need to include the weight w_i */
/* A factor rho0 / 12.0 = 2D rho0 / bc^2 is included in the velocities */
      lp = obj[inode1];
      rx = i - (lp->posx); ry = j - (lp->posy); rz = k - (lp->posz);
      ux = rho0 * ((lp->velx) + (lp->angy)*rz - (lp->angz)*ry) / 12.0;
      uy = rho0 * ((lp->vely) + (lp->angz)*rx - (lp->angx)*rz) / 12.0;
      uz = rho0 * ((lp->velz) + (lp->angx)*ry - (lp->angy)*rx) / 12.0;
      if (m & 0x0001) { /* Check for boundary link c_01 = ( 1, 0, 0) */
	inode2 = esk*(esj*ip+j)+k;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[2+15*inode2] + piece + 8.0*fac*ux;
      }
      if (m & 0x0002) { /* Check for boundary link c_02 = (-1, 0, 0) */
	inode2 = esk*(esj*im+j)+k;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[1+15*inode2] + piece - 8.0*fac*ux;
      }
      if (m & 0x0004) { /* Check for boundary link c_03 = ( 0, 1, 0) */
	inode2 = esk*(esj*i+jp)+k;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[4+15*inode2] + piece + 8.0*fac*uy;
      }
      if (m & 0x0008) { /* Check for boundary link c_04 = ( 0,-1, 0) */
	inode2 = esk*(esj*i+jm)+k;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[3+15*inode2] + piece - 8.0*fac*uy;
      }
      if (m & 0x0010) { /* Check for boundary link c_05 = ( 0, 0, 1) */
	inode2 = esk*(esj*i+j)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[6+15*inode2] + piece + 8.0*fac*uz;
      }
      if (m & 0x0020) { /* Check for boundary link c_06 = ( 0, 0,-1) */
	inode2 = esk*(esj*i+j)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[5+15*inode2] + piece - 8.0*fac*uz;
      }
      if (m & 0x0040) { /* Check for boundary link c_07 = ( 1, 1, 1) */
	inode2 = esk*(esj*ip+jp)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[14+15*inode2] + piece + fac*(ux+uy+uz);
      }
      if (m & 0x0080) { /* Check for boundary link c_08 = (-1, 1, 1) */
	inode2 = esk*(esj*im+jp)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[13+15*inode2] + piece - fac*(ux-uy-uz);
      }
      if (m & 0x0100) { /* Check for boundary link c_09 = ( 1,-1, 1) */
	inode2 = esk*(esj*ip+jm)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[12+15*inode2] + piece + fac*(ux-uy+uz);
      }
      if (m & 0x0200) { /* Check for boundary link c_10 = (-1,-1, 1) */
	inode2 = esk*(esj*im+jm)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[11+15*inode2] + piece - fac*(ux+uy-uz);
      }
      if (m & 0x0400) { /* Check for boundary link c_11 = ( 1, 1,-1) */
	inode2 = esk*(esj*ip+jp)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[10+15*inode2] + piece + fac*(ux+uy-uz);
      }
      if (m & 0x0800) { /* Check for boundary link c_12 = (-1, 1,-1) */
	inode2 = esk*(esj*im+jp)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[9+15*inode2] + piece - fac*(ux-uy+uz);
      }
      if (m & 0x1000) { /* Check for boundary link c_13 = ( 1,-1,-1) */
	inode2 = esk*(esj*ip+jm)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[8+15*inode2] + piece + fac*(ux-uy-uz);
      }
      if (m & 0x2000) { /* Check for boundary link c_14 = (-1,-1,-1) */
	inode2 = esk*(esj*im+jm)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[7+15*inode2] + piece - fac*(ux+uy+uz);
      }
    }
  }
/* Update the scalar density in fluid sites, use p to store */
/* the adsorbed scalar density on adsorption sites */
/* The latter should be used with care if objects are physically moved */
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    ptemp = newp[inode];
    if (obj[inode] == NULL) p[inode] = ptemp;
    else {
      if (mask[inode] & 0x4000) p[inode] += ptemp;
      else p[inode] = 0.0;
    }
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpmove done\n", mp->name, mp->step);
  }
}

/* Perform an adjoint moment propagation step for stationary objects. 
 * The only difference here is that the propagation probabilities are
 * used to gather scalar density rather than scatter scalar density.
 * Adjoint in the sense sum_r psi(r) P(r, r') phi(r') is preserved.
 */

void mpadjoint(MPSys mp) {
  int i, j, k, ip, jp, kp, im, jm, km, esi, esj, esk, m, inode, inode1, inode2;
  int *mask;
  LBReal alpha, rho0, rhot, ptemp, fac, piece, piece8;
  LBReal piece16, onebyb, eightbyb;
  LBReal *a, *p, *newp, *rho;
  LBObj *obj;
/* Calculate some constants appropriate to D3Q15(r=16) */
  onebyb = 1.0/72.0; eightbyb = 8.0/72.0;
/* Check for correct place in LB update cycle, should be post-collision */
  if (mp->lb->status != POST_COLLIDE) {
    printf("Moment propagation apparently out of position\n");
    wrstatus(mp->lb, stdout);
  }
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; mask = mp->lb->mask; a = mp->lb->a;
  newp = mp->newp; p = mp->p; rho = mp->lb->rho;
  rho0 = mp->lb->rho0; alpha = mp->alpha;
/* Clear array that stores new scalar density */
  for (inode = 0; inode<(mp->lb->nnodes); inode++) newp[inode] = 0.0;
/* Run through all sites on the lattice */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode1 = esk*(esj*i+j)+k;
    ip = i + 1; if (ip >= esi) ip = 0;
    jp = j + 1; if (jp >= esj) jp = 0;
    kp = k + 1; if (kp >= esk) kp = 0;
    im = i - 1; if (im < 0) im = esi-1;
    jm = j - 1; if (jm < 0) jm = esj-1;
    km = k - 1; if (km < 0) km = esk-1;
/* Check for a fluid node, move stuff in appropriate directions */
/* Stuff that propagates along boundary links is taken care of below */
    if (obj[inode1] == NULL) {

      rhot = rho[inode1];
      fac = 1.0/rhot; piece = onebyb*(rho0/rhot - alpha);
      piece8 = 8.0*piece; piece16 = 16.0*piece;

      newp[inode1] += p[esk*(esj*i+j)+k]*(fac*a[15*inode1] + piece16 + alpha);
      newp[inode1] += p[esk*(esj*ip+j)+k]*(fac*a[1+15*inode1] + piece8);
      newp[inode1] += p[esk*(esj*im+j)+k]*(fac*a[2+15*inode1] + piece8);
      newp[inode1] += p[esk*(esj*i+jp)+k]*(fac*a[3+15*inode1] + piece8);
      newp[inode1] += p[esk*(esj*i+jm)+k]*(fac*a[4+15*inode1] + piece8);
      newp[inode1] += p[esk*(esj*i+j)+kp]*(fac*a[5+15*inode1] + piece8);
      newp[inode1] += p[esk*(esj*i+j)+km]*(fac*a[6+15*inode1] + piece8);
      newp[inode1] += p[esk*(esj*ip+jp)+kp]*(fac*a[7+15*inode1] + piece);
      newp[inode1] += p[esk*(esj*im+jp)+kp]*(fac*a[8+15*inode1] + piece);
      newp[inode1] += p[esk*(esj*ip+jm)+kp]*(fac*a[9+15*inode1] + piece);
      newp[inode1] += p[esk*(esj*im+jm)+kp]*(fac*a[10+15*inode1] + piece);
      newp[inode1] += p[esk*(esj*ip+jp)+km]*(fac*a[11+15*inode1] + piece);
      newp[inode1] += p[esk*(esj*im+jp)+km]*(fac*a[12+15*inode1] + piece);
      newp[inode1] += p[esk*(esj*ip+jm)+km]*(fac*a[13+15*inode1] + piece);
      newp[inode1] += p[esk*(esj*im+jm)+km]*(fac*a[14+15*inode1] + piece);

    } else if ((m = mask[inode1]) && !(m & 0x4000)) {
/* Interior node with boundary links found */
/* For stationary objects, account for these by bounce back into the */
/* appropriate fluid node. */
/* Implement adsorption here by simple adsorbing the total amount */

/*  THIS NEEDS SORTING OUT */

      if (m & 0x0001) { /* Check for boundary link c_01 = ( 1, 0, 0) */
	inode2 = esk*(esj*ip+j)+k;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[2+15*inode2] + piece;
      }
      if (m & 0x0002) { /* Check for boundary link c_02 = (-1, 0, 0) */
	inode2 = esk*(esj*im+j)+k;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[1+15*inode2] + piece;
      }
      if (m & 0x0004) { /* Check for boundary link c_03 = ( 0, 1, 0) */
	inode2 = esk*(esj*i+jp)+k;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[4+15*inode2] + piece;
      }
      if (m & 0x0008) { /* Check for boundary link c_04 = ( 0,-1, 0) */
	inode2 = esk*(esj*i+jm)+k;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[3+15*inode2] + piece;
      }
      if (m & 0x0010) { /* Check for boundary link c_05 = ( 0, 0, 1) */
	inode2 = esk*(esj*i+j)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[6+15*inode2] + piece;
      }
      if (m & 0x0020) { /* Check for boundary link c_06 = ( 0, 0,-1) */
	inode2 = esk*(esj*i+j)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = eightbyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[5+15*inode2] + piece;
      }
      if (m & 0x0040) { /* Check for boundary link c_07 = ( 1, 1, 1) */
	inode2 = esk*(esj*ip+jp)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[14+15*inode2] + piece;
      }
      if (m & 0x0080) { /* Check for boundary link c_08 = (-1, 1, 1) */
	inode2 = esk*(esj*im+jp)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[13+15*inode2] + piece;
      }
      if (m & 0x0100) { /* Check for boundary link c_09 = ( 1,-1, 1) */
	inode2 = esk*(esj*ip+jm)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[12+15*inode2] + piece;
      }
      if (m & 0x0200) { /* Check for boundary link c_10 = (-1,-1, 1) */
	inode2 = esk*(esj*im+jm)+kp;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[11+15*inode2] + piece;
      }
      if (m & 0x0400) { /* Check for boundary link c_11 = ( 1, 1,-1) */
	inode2 = esk*(esj*ip+jp)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[10+15*inode2] + piece;
      }
      if (m & 0x0800) { /* Check for boundary link c_12 = (-1, 1,-1) */
	inode2 = esk*(esj*im+jp)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[9+15*inode2] + piece;
      }
      if (m & 0x1000) { /* Check for boundary link c_13 = ( 1,-1,-1) */
	inode2 = esk*(esj*ip+jm)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[8+15*inode2] + piece;
      }
      if (m & 0x2000) { /* Check for boundary link c_14 = (-1,-1,-1) */
	inode2 = esk*(esj*im+jm)+km;
        rhot = rho[inode2]; ptemp = p[inode2];
        fac = ptemp/rhot; piece = onebyb*ptemp*(rho0/rhot - alpha);
        newp[inode2] += fac*a[7+15*inode2] + piece;
      }
    }
  }
/* Update the scalar density in fluid sites, use p to store */
/* the adsorbed scalar density on adsorption sites */
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    ptemp = newp[inode];
    if (obj[inode] == NULL) p[inode] = ptemp;
    else {
      if (mask[inode] & 0x4000) p[inode] += ptemp;
      else p[inode] = 0.0;
    }
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpstat done\n", mp->name, mp->step);
  }
}

/* End of mpcore.c */
