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

#define LBMEASURE_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sunlightlb.h"


/* Measure some quantities. */

void measure(LBSys lb) {
/* Should be called after a collision step. */
  if (lb->status != POST_COLLIDE) {
    printf("Measurements apparently out of position\n"); wrstatus(lb, stdout);
  }
/* Start accumulating measurements from zero */
  zeroforce(lb);
  lb->vmeanx = lb->vmeany = lb->vmeanz = 0.0;
/* Accumulate data for two LB time steps */
  accvmean(lb, 0.25);
  bblforce(lb); propagate(lb); collide(lb);
  accvmean(lb, 0.5);
  bblforce(lb); propagate(lb); collide(lb);
  accvmean(lb, 0.25);
/* Halve node force to account for the two steps */
  normforce(lb, 0.5);
/* Finally distribute the forces to the particles */
  distribf(lb);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : measure done\n", lb->step);
}


/* Measure instantaneous bulk flow rate, corrected for half
 * added pressure gradient.  This is added to the stored data 
 * with a weight fac.  A scalar probability amplitude can be used
 * to bias the average.  It does not have to be normalised.
 */

void accvmean(LBSys lb, LBReal prefac) {
  int inode, nfluid;
  LBReal t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14;
  LBReal fac, rho, ux, uy, uz, uxtot, uytot, uztot, wfac, wnorm;
  LBReal *aa, *psi;
  LBReal *bx, *by, *bz;
  LBReal rho0;
  LBObj *obj;
  rho0 = lb->rho0; obj = lb->obj; psi = lb->psi;
  aa = lb->a; bx = lb->bx; by = lb->by; bz = lb->bz;
  fac = uxtot = uytot = uztot = wnorm = 0.0;
/* This should be called pre- or post-collision (not post-bounceback) */
  fac = 0.0;
  if (lb->vel_corrected) switch (lb->status) {
    case POST_PROPAGATE: fac = 0.5; break;
    case POST_COLLIDE: fac = -0.5; break;
    default:
      printf("Accvmean apparently of position\n"); wrstatus(lb, stdout);
      break;
  }
/* Note how this just sums over fluid nodes */
  wfac = 1.0; nfluid = 0;
  for (inode=0; inode<(lb->nnodes); inode++) {
    if (obj[inode] != NULL) continue;
    t00 = aa[ 0+15*inode]; t01 = aa[ 1+15*inode]; t02 = aa[ 2+15*inode];
    t03 = aa[ 3+15*inode]; t04 = aa[ 4+15*inode]; t05 = aa[ 5+15*inode];
    t06 = aa[ 6+15*inode]; t07 = aa[ 7+15*inode]; t08 = aa[ 8+15*inode];
    t09 = aa[ 9+15*inode]; t10 = aa[10+15*inode]; t11 = aa[11+15*inode];
    t12 = aa[12+15*inode]; t13 = aa[13+15*inode]; t14 = aa[14+15*inode];
/* Calculate density and velocity, correct for body force */
    rho = rho0 + t00+t01+t02+t03+t04+t05+t06
                     +t07+t08+t09+t10+t11+t12+t13+t14;
    ux = (t01-t02+t07-t08+t09-t10+t11-t12+t13-t14 + fac*bx[inode])/rho;
    uy = (t03-t04+t07+t08-t09-t10+t11+t12-t13-t14 + fac*by[inode])/rho;
    uz = (t05-t06+t07+t08+t09+t10-t11-t12-t13-t14 + fac*bz[inode])/rho;
    if (psi != NULL) { wfac = psi[inode]; wfac = wfac * wfac; }
/* Sum up results to get mean */
    uxtot += wfac*ux; uytot += wfac*uy; uztot += wfac*uz; wnorm += wfac;
    nfluid++;
  }
/* Take out the normalisation */
  wfac = nfluid / wnorm;
  uxtot *= wfac; uytot *= wfac; uztot *= wfac;
/* Add to stored mean flow rate */
  lb->vmeanx += prefac * uxtot / (lb->nnodes);
  lb->vmeany += prefac * uytot / (lb->nnodes);
  lb->vmeanz += prefac * uztot / (lb->nnodes);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : accvmean done\n", lb->step);
}

/* Zero the stored velocity field at all the fluid nodes */

void vzero(LBSys lb) {
  int inode;
  LBReal *ux, *uy, *uz;
  LBObj *obj;
  ux = lb->ux; uy = lb->uy; uz = lb->uz; obj = lb->obj;
  for (inode=0; inode<(lb->nnodes); inode++) {
    if (obj[inode] == NULL) {
      ux[inode] = uy[inode] = uz[inode] = 0.0;
    }
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : vzero done\n", lb->step);
}


/* Capture instantaneous flow field on the fluid nodes, corrected for
 * half added pressure gradient.  This is added to the stored data
 * with a weight fac.
 */

void vgrab(LBSys lb, LBReal prefac) {
  int inode;
  LBReal t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14;
  LBReal fac, rho;
  LBReal *ux, *uy, *uz;
  LBReal *aa;
  LBReal *bx, *by, *bz;
  LBReal rho0;
  LBObj *obj;
  ux = lb->ux; uy = lb->uy; uz = lb->uz;
  rho0 = lb->rho0; obj = lb->obj;
  aa = lb->a; bx = lb->bx; by = lb->by; bz = lb->bz;
/* This should be called pre- or post-collision (not post-bounceback) */
  fac = 0.0;
  if (lb->vel_corrected) switch (lb->status) {
    case POST_PROPAGATE: fac = 0.5; break;
    case POST_COLLIDE: fac = -0.5; break;
    default:
      printf("vgrab apparently of position\n"); wrstatus(lb, stdout);
      break;
  }
/* Note how this just includes fluid nodes */
  for (inode=0; inode<(lb->nnodes); inode++) {
    if (obj[inode] != NULL) continue;
    t00 = aa[ 0+15*inode]; t01 = aa[ 1+15*inode]; t02 = aa[ 2+15*inode];
    t03 = aa[ 3+15*inode]; t04 = aa[ 4+15*inode]; t05 = aa[ 5+15*inode];
    t06 = aa[ 6+15*inode]; t07 = aa[ 7+15*inode]; t08 = aa[ 8+15*inode];
    t09 = aa[ 9+15*inode]; t10 = aa[10+15*inode]; t11 = aa[11+15*inode];
    t12 = aa[12+15*inode]; t13 = aa[13+15*inode]; t14 = aa[14+15*inode];
/* Calculate density and velocity, correct for body force */
    rho = rho0 + t00+t01+t02+t03+t04+t05+t06
                     +t07+t08+t09+t10+t11+t12+t13+t14;
    ux[inode] += prefac * (t01-t02+t07-t08+t09-t10+t11-t12+t13-t14 + fac*bx[inode])/rho;
    uy[inode] += prefac * (t03-t04+t07+t08-t09-t10+t11+t12-t13-t14 + fac*by[inode])/rho;
    uz[inode] += prefac * (t05-t06+t07+t08+t09+t10-t11-t12-t13-t14 + fac*bz[inode])/rho;
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : vgrab done\n", lb->step);
}


/* Write the measured velocity field at all nodes including solid
 nodes.  Some thought is still required to treat the skip_solid = 1 case
 properly */

void vsave(LBSys lb, char *savefile) {
  int i, j, k, inode;  
  int esi, esj, esk;
  int skip_solid;
  LBReal ux, uy, uz, rho0, fac;
  LBObj *obj;
  FILE *fp;
  esi = lb->esi; esj = lb->esj; esk = lb->esk;
  rho0 = lb->rho0; fac = 12.0 / rho0;
  if ((fp = fopen(savefile,"w")) == NULL) {
    printf("savefile: %s could not be opened\n", savefile); return;
  }
  obj = lb->obj; skip_solid = lb->skip_solid;
/* Compute and save velocity on a per-node basis */
/* This is equivalent to traversing v[i][j][k] arrays with k varying fastest */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    if (obj[inode] != NULL) {
      if (skip_solid) {
	ux = uy = uz = 0.0;
      } else {
	ux = fac*lb->ux[inode]; uy = fac*lb->uy[inode]; uz = fac*lb->uz[inode];
      }
    } else {
      ux = lb->ux[inode]; uy = lb->uy[inode]; uz = lb->uz[inode];
    }
    fprintf(fp,"%g %g %g %i %i %i\n", ux, uy, uz, i, j, k);
  }
  fclose(fp);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : lbsavevel done\n", lb->step);
}


/* This sets the individual node force array to zero */

void zeroforce(LBSys lb) {
  int inode;
  LBReal *fx, *fy, *fz;
  fx = lb->fx; fy = lb->fy; fz = lb->fz;
  for (inode=0; inode<(lb->nnodes); inode++) {
    fx[inode] = fy[inode] = fz[inode] = 0.0;
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : zeroforce done\n", lb->step);
}


/* This multiplies the individual force array by a given factor */

void normforce(LBSys lb, LBReal fac) {
  int inode;
  LBReal *fx, *fy, *fz;
  fx = lb->fx; fy = lb->fy; fz = lb->fz;
  for (inode=0; inode<(lb->nnodes); inode++) {
    fx[inode] *= fac; fy[inode] *= fac; fz[inode] *= fac;
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("LB %i : normforce(%f) done\n", lb->step, fac);
  }
}


/* This distributes the node forces across the objects */

void distribf(LBSys lb) {
  int i, j, k, inode, esi, esj, esk;
  LBReal *fx, *fy, *fz;
  LBReal rx, ry, rz, forcex, forcey, forcez;
  LBObj p,*obj;
  fx = lb->fx; fy = lb->fy; fz = lb->fz;
  obj = lb->obj; esi = lb->esi; esj = lb->esj; esk = lb->esk;
/* First zero the object force and torque arrays */
  for (p=lb->obj_head; p!=NULL; p=p->next_obj) {
    p->frcx = p->frcy = p->frcz = 0.0;
    p->trqx = p->trqy = p->trqz = 0.0;
  }
/* Then distribute node forces across objects */
/* Need i, j, k resolution in here to work out the node position */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    if ((p = obj[inode]) == NULL) continue;
    rx = i - (p->posx); ry = j - (p->posy); rz = k - (p->posz);
    forcex = fx[inode]; forcey = fy[inode]; forcez = fz[inode];
    p->frcx += forcex;
    p->frcy += forcey;
    p->frcz += forcez;
    p->trqx += ry * forcez - rz * forcey;
    p->trqy += rz * forcex - rx * forcez;
    p->trqz += rx * forcey - ry * forcex;
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : distribf done\n", lb->step);
}


/* Print out basic results. */

void wrmeas(LBSys lb) {
  int inode;
  double discrx, discry, discrz, discr, norm;
  double fxtot, fytot, fztot;
  double fxsum, fysum, fzsum;
  LBObj p;
  printf("\n\nProperties at lb steps %i to %i\n",((lb->step)-2),(lb->step));
  printf(" mean flow rate = %0.6f %0.6f %0.6f\n", lb->vmeanx, lb->vmeany, lb->vmeanz);
  fxtot = fytot = fztot = 0.0;
  for (p=lb->obj_head; p!=NULL; p=p->next_obj) {
    if (slb_verbosity & FORCE_DETAILS) {
      printf(" object id %3i (base id %3i), force = %8.3f %8.3f %8.3f\n",
            p->id, p->base_obj->id, p->frcx, p->frcy, p->frcz);
      printf("                             torque = %8.3f %8.3f %8.3f\n",
              p->trqx, p->trqy, p->trqz);
    }
    fxtot += p->frcx; fytot += p->frcy; fztot += p->frcz; 
  }
  printf("            total force = %8.3f %8.3f %8.3f\n",
            fxtot, fytot, fztot);
  fxsum = fysum = fzsum = 0.0;
  for (inode=0; inode<(lb->nnodes); inode++) {
    if (lb->obj[inode] == NULL) {
      fxsum += lb->bx[inode]; fysum += lb->by[inode]; fzsum += lb->bz[inode];
    }
  }
  printf("  total force should be = %8.3f %8.3f %8.3f\n", fxsum, fysum, fzsum);
  discrx = fxtot - fxsum; discry = fytot - fysum; discrz = fztot - fzsum;
  discr = sqrt(discrx*discrx + discry*discry + discrz*discrz);
  norm = sqrt(fxsum*fxsum + fysum*fysum + fzsum*fzsum);
  if (norm != 0.0) discr = discr/norm;
  printf(" force discrepancy = %0.4f (norm = %0.4f)\n", discr, norm);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : wrforce done\n", lb->step);
}


/* Return the fractional discrepancy between the total applied force
 * and the force on the objects.
 */

LBReal discrf(LBSys lb) {
  int inode;
  double discrx, discry, discrz, norm;
  double fxtot, fytot, fztot;
  double fxsum, fysum, fzsum;
  LBReal discr;
  LBObj p;
  fxtot = fytot = fztot = 0.0;
  for (p=lb->obj_head; p!=NULL; p=p->next_obj) {
    fxtot += p->frcx; fytot += p->frcy; fztot += p->frcz; 
  }
  fxsum = fysum = fzsum = 0.0;
  for (inode=0; inode<(lb->nnodes); inode++) {
    if (lb->obj[inode] == NULL) {
      fxsum += lb->bx[inode]; fysum += lb->by[inode]; fzsum += lb->bz[inode];
    }
  }
  discrx = fxtot - fxsum; discry = fytot - fysum; discrz = fztot - fzsum;
  discr = sqrt(discrx*discrx + discry*discry + discrz*discrz);
  norm = sqrt(fxsum*fxsum + fysum*fysum + fzsum*fzsum) ;
  if (norm != 0.0) discr = discr/norm;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : wrforce done, discr = %0.4f\n", 
					    lb->step, discr);
  return discr;
}

/* End of lbmeasure.c */
