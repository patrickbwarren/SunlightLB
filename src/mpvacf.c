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

#define MPVACF_C

#include <stdio.h>
#include <stdlib.h>
#include "sunlightlb.h"

/* Local variables */

static int vacfl1, vacfl2, nvacf, maxnvacf;
static LBReal *vacfarr;

/* Set up to calculate tagged particle velocity autocorrelation function.
 * The first letter of vacfdir indicates the start direction, the second
 * letter indicates the end direction, thus "xx" and so on.  Storage space
 * for reqnvacf vacf steps is also allocated.  Two tracer fields are used,
 * the first is tagged with the velocity, and the second keeps track
 * of the number of tracer particles.
 */

void mpvacfinit(MPSys mp, MPSys mp2, char *vacfdir, int reqnvacf) {
  int i, j, k, esi, esj, esk, ip, jp, kp, im, jm, km, m, inode, inode1, inode2;
  int *mask, *c;
  LBReal t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14;
  LBReal sxx, sxy, sxz, syy, syz, szz, sxyz, t1p2, t3p4, t5p6, t0p6;
  LBReal tlppp, tlmpm, tlpmm, tlmmp, thppp, thmpm, thpmm, thmmp, tlhp;
  LBReal piece, piece8, onebyb, eightbyb, css, vacftot, vacfval, wfac;
  LBReal rho0, alpha, ntracers;
  LBReal *a, *p, *p2, *rho, *psi;
  LBObj *obj;
/* Allocate storage for (reqnvacf+1) reals */
  if (vacfarr != NULL) free(vacfarr);
  vacfarr = (LBReal *) malloc((reqnvacf+1)*sizeof(LBReal));
  if (vacfarr == NULL) {
    fprintf(stderr,"Not enought room at the inn for %i reals\n",
	    reqnvacf+1);
    return;
  } 
  maxnvacf = reqnvacf;
/* Calculate some constants for D3Q15(r=16) system */
  onebyb = 1.0/72.0; eightbyb = 8.0/72.0; css = 1.0/3.0;
/* Map initial and final directions onto vacfl1 and vacfl2 */
  switch (vacfdir[0]) {
    case 'x': case 'i': vacfl1 = 0; break;
    case 'y': case 'j': vacfl1 = 1; break;
    case 'z': case 'k': vacfl1 = 2; break;
    case '\0': error("mpvacfinit: please supply initial direction");
    default: error("mpvacfinit: unrecognised initial direction");
  }
  switch (vacfdir[1]) {
    case 'x': case 'i': vacfl2 = 0; break;
    case 'y': case 'j': vacfl2 = 1; break;
    case 'z': case 'k': vacfl2 = 2; break;
    case '\0': error("mpvacfinit: please supply final direction");
    default: error("mpvacfinit: unrecognised final direction");
  }
/* Calculate P(r,1) = sum_i w_i c_ia p_i(r-c_i) rho(r-c_i) psi(r-c_i) */
/* and P2(r,1) = sum_i w_i p_i(r-c_i) rho(r-c_i) psi(r-c_i) */
/* where the tracer propagation probabilities are given by */
/* p_i = (f_i^+ / rho - alpha / b) + alpha delta_{i0} / w_0 */
/* Note that for this to work, routine should be called post collision */
  if (mp->lb->status != POST_COLLIDE) {
    fprintf(stderr,"vacf initialisation apparently out of position\n");
    wrstatus(mp->lb, stdout);
  }
  esi = mp->lb->esi; esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; mask = mp->lb->mask; a = mp->lb->a; c = mp->lb->c;
  rho0 = mp->lb->rho0; alpha = mp->alpha;
  p = mp->p; p2 = mp2->p; wfac = 1.0;
  rho = mp->lb->rho; psi = mp->lb->psi;

  ntracers = 0.0;

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
      piece = onebyb*(rho0 - alpha*rho[inode1]); piece8 = 8.0*piece;
      if (psi != NULL) wfac = psi[inode1];
      /* Omit i=0 here since c_0 = 0 */
      p[esk*(esj*ip+j)+k]   += wfac*c[3+vacfl1]*(a[1+15*inode1] + piece8);
      p[esk*(esj*im+j)+k]   += wfac*c[6+vacfl1]*(a[2+15*inode1] + piece8);
      p[esk*(esj*i+jp)+k]   += wfac*c[9+vacfl1]*(a[3+15*inode1] + piece8);
      p[esk*(esj*i+jm)+k]   += wfac*c[12+vacfl1]*(a[4+15*inode1] + piece8);
      p[esk*(esj*i+j)+kp]   += wfac*c[15+vacfl1]*(a[5+15*inode1] + piece8);
      p[esk*(esj*i+j)+km]   += wfac*c[18+vacfl1]*(a[6+15*inode1] + piece8);
      p[esk*(esj*ip+jp)+kp] += wfac*c[21+vacfl1]*(a[7+15*inode1] + piece);
      p[esk*(esj*im+jp)+kp] += wfac*c[24+vacfl1]*(a[8+15*inode1] + piece);
      p[esk*(esj*ip+jm)+kp] += wfac*c[27+vacfl1]*(a[9+15*inode1] + piece);
      p[esk*(esj*im+jm)+kp] += wfac*c[30+vacfl1]*(a[10+15*inode1] + piece);
      p[esk*(esj*ip+jp)+km] += wfac*c[33+vacfl1]*(a[11+15*inode1] + piece);
      p[esk*(esj*im+jp)+km] += wfac*c[36+vacfl1]*(a[12+15*inode1] + piece);
      p[esk*(esj*ip+jm)+km] += wfac*c[39+vacfl1]*(a[13+15*inode1] + piece);
      p[esk*(esj*im+jm)+km] += wfac*c[42+vacfl1]*(a[14+15*inode1] + piece);
      p2[esk*(esj*ip+j)+k]   += wfac*(a[15*inode1] + 16.0*piece
				                 + alpha*rho[inode1]);
      p2[esk*(esj*ip+j)+k]   += wfac*(a[1+15*inode1] + piece8);
      p2[esk*(esj*im+j)+k]   += wfac*(a[2+15*inode1] + piece8);
      p2[esk*(esj*i+jp)+k]   += wfac*(a[3+15*inode1] + piece8);
      p2[esk*(esj*i+jm)+k]   += wfac*(a[4+15*inode1] + piece8);
      p2[esk*(esj*i+j)+kp]   += wfac*(a[5+15*inode1] + piece8);
      p2[esk*(esj*i+j)+km]   += wfac*(a[6+15*inode1] + piece8);
      p2[esk*(esj*ip+jp)+kp] += wfac*(a[7+15*inode1] + piece);
      p2[esk*(esj*im+jp)+kp] += wfac*(a[8+15*inode1] + piece);
      p2[esk*(esj*ip+jm)+kp] += wfac*(a[9+15*inode1] + piece);
      p2[esk*(esj*im+jm)+kp] += wfac*(a[10+15*inode1] + piece);
      p2[esk*(esj*ip+jp)+km] += wfac*(a[11+15*inode1] + piece);
      p2[esk*(esj*im+jp)+km] += wfac*(a[12+15*inode1] + piece);
      p2[esk*(esj*ip+jm)+km] += wfac*(a[13+15*inode1] + piece);
      p2[esk*(esj*im+jm)+km] += wfac*(a[14+15*inode1] + piece);
/* Check for an interior node with boundary links */
/* If this is found, account for it by bounce back into the */
/* appropriate fluid node. Adsorbing sites are discounted here. */
    } else if ((m = mask[inode1]) && !(m & 0x4000)) {
      if (m & 0x0001) { /* Check for boundary link c_01 = ( 1, 0, 0) */
	inode2 = esk*(esj*ip+j)+k;
        piece8 = eightbyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[6+vacfl1]*(a[2+15*inode2] + piece8);
        p2[inode2] += wfac*(a[2+15*inode2] + piece8);
      }
      if (m & 0x0002) { /* Check for boundary link c_02 = (-1, 0, 0) */
	inode2 = esk*(esj*im+j)+k;
        piece8 = eightbyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[3+vacfl1]*(a[1+15*inode2] + piece8);
        p2[inode2] += wfac*(a[1+15*inode2] + piece8);
      }
      if (m & 0x0004) { /* Check for boundary link c_03 = ( 0, 1, 0) */
	inode2 = esk*(esj*i+jp)+k;
        piece8 = eightbyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[12+vacfl1]*(a[4+15*inode2] + piece8);
        p2[inode2] += wfac*(a[4+15*inode2] + piece8);
      }
      if (m & 0x0008) { /* Check for boundary link c_04 = ( 0,-1, 0) */
	inode2 = esk*(esj*i+jm)+k;
        piece8 = eightbyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[9+vacfl1]*(a[3+15*inode2] + piece8);
        p2[inode2] += wfac*(a[3+15*inode2] + piece8);
      }
      if (m & 0x0010) { /* Check for boundary link c_05 = ( 0, 0, 1) */
	inode2 = esk*(esj*i+j)+kp;
        piece8 = eightbyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[18+vacfl1]*(a[6+15*inode2] + piece8);
        p2[inode2] += wfac*(a[6+15*inode2] + piece8);
      }
      if (m & 0x0020) { /* Check for boundary link c_06 = ( 0, 0,-1) */
	inode2 = esk*(esj*i+j)+km;
        piece8 = eightbyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[15+vacfl1]*(a[5+15*inode2] + piece8);
        p2[inode2] += wfac*(a[5+15*inode2] + piece8);
      }
      if (m & 0x0040) { /* Check for boundary link c_07 = ( 1, 1, 1) */
	inode2 = esk*(esj*ip+jp)+kp;
        piece = onebyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[42+vacfl1]*(a[14+15*inode2] + piece);
        p2[inode2] += wfac*(a[14+15*inode2] + piece);
      }
      if (m & 0x0080) { /* Check for boundary link c_08 = (-1, 1, 1) */
	inode2 = esk*(esj*im+jp)+kp;
        piece = onebyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[39+vacfl1]*(a[13+15*inode2] + piece);
        p2[inode2] += wfac*(a[13+15*inode2] + piece);
      }
      if (m & 0x0100) { /* Check for boundary link c_09 = ( 1,-1, 1) */
	inode2 = esk*(esj*ip+jm)+kp;
        piece = onebyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[36+vacfl1]*(a[12+15*inode2] + piece);
        p2[inode2] += wfac*(a[12+15*inode2] + piece);
      }
      if (m & 0x0200) { /* Check for boundary link c_10 = (-1,-1, 1) */
	inode2 = esk*(esj*im+jm)+kp;
        piece = onebyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*(a[11+15*inode2] + piece);
        p2[inode2] += wfac*(a[11+15*inode2] + piece);
      }
      if (m & 0x0400) { /* Check for boundary link c_11 = ( 1, 1,-1) */
	inode2 = esk*(esj*ip+jp)+km;
        piece = onebyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[30+vacfl1]*(a[10+15*inode2] + piece);
        p2[inode2] += wfac*(a[10+15*inode2] + piece);
      }
      if (m & 0x0800) { /* Check for boundary link c_12 = (-1, 1,-1) */
	inode2 = esk*(esj*im+jp)+km;
        piece = onebyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[27+vacfl1]*(a[9+15*inode2] + piece);
        p2[inode2] += wfac*(a[9+15*inode2] + piece);
      }
      if (m & 0x1000) { /* Check for boundary link c_13 = ( 1,-1,-1) */
	inode2 = esk*(esj*ip+jm)+km;
        piece = onebyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[24+vacfl1]*(a[8+15*inode2] + piece);
        p2[inode2] += wfac*(a[8+15*inode2] + piece);
      }
      if (m & 0x2000) { /* Check for boundary link c_14 = (-1,-1,-1) */
	inode2 = esk*(esj*im+jm)+km;
        piece = onebyb*(rho0 - alpha*rho[inode2]);
	if (psi != NULL) wfac = psi[inode2];
        p[inode2] += wfac*c[21+vacfl1]*(a[7+15*inode2] + piece);
        p2[inode2] += wfac*(a[7+15*inode2] + piece);
      }
    }
  }
/* Remove the stuff that ended up in solid object nodes */
  for (inode=0; inode<(mp->lb->nnodes); inode++) {
    if (obj[inode] != NULL) { p[inode] = 0.0; p2[inode] = 0.0; }
  }
/* Compute and save the zero time value of the vacf, namely */
/* 1/Ntracers * sum_r [ rho css (1-alpha) delta_ab + S_ab] psi^2 */
/* where Ntracers = sum_r rho psi^2 */
  vacftot = ntracers = 0.0; wfac = 1.0;
  for (inode=0; inode<(mp->lb->nnodes); inode++) if (obj[inode] == NULL) {
    t00 = a[15*inode]; t01 = a[1+15*inode]; t02 = a[2+15*inode];
    t03 = a[3+15*inode]; t04 = a[4+15*inode]; t05 = a[5+15*inode];
    t06 = a[6+15*inode]; t07 = a[7+15*inode]; t08 = a[8+15*inode];
    t09 = a[9+15*inode]; t10 = a[10+15*inode]; t11 = a[11+15*inode];
    t12 = a[12+15*inode]; t13 = a[13+15*inode]; t14 = a[14+15*inode];
/* Calculate some partial sums */
    t1p2 = t01 + t02; t3p4 = t03 + t04; t5p6 = t05 + t06;
    tlppp = t07 + t08 + t09 + t10; thppp = t11 + t12 + t13 + t14;
    tlmpm = t07 - t08 + t09 - t10; thmpm = t11 - t12 + t13 - t14;
    tlpmm = t07 + t08 - t09 - t10; thpmm = t11 + t12 - t13 - t14;
    tlmmp = t07 - t08 - t09 + t10; thmmp = t11 - t12 - t13 + t14;
    t0p6 = t00 + t1p2 + t3p4 + t5p6; tlhp = tlppp + thppp;
/* Project out stress */
    sxyz = 2.0*(t0p6 + tlhp)/3.0 - t0p6;
    sxx = sxyz + t1p2;    syy = sxyz + t3p4;    szz = sxyz + t5p6;
    sxy = tlmmp + thmmp;  sxz = tlmpm - thmpm;  syz = tlpmm - thpmm;
/* Add appropriate quantities to ntracers and vacftot */
    if (psi != NULL) wfac = psi[inode];
    ntracers += wfac * p2[inode];
    wfac = wfac * wfac;
    ntracers += wfac * rho[inode];
    piece = rho[inode] * css * (1.0 - alpha);
    switch (10*vacfl1 + vacfl2) {
      case 00: vacftot += wfac*(piece + sxx); break; 
      case 11: vacftot += wfac*(piece + syy); break; 
      case 22: vacftot += wfac*(piece + szz); break; 
      case 01: case 10: vacftot += wfac*sxy; break;
      case 02: case 20: vacftot += wfac*sxz; break;
      case 12: case 21: vacftot += wfac*syz; break;
    }
  }
  vacfval = vacftot / ntracers;
  vacfarr[0] = vacfval; nvacf = 1;
  printf("MP '%s' %i : Initialised for vacf, init drn %c, final drn %c\n",
         mp->name, mp->step, vacfdir[0], vacfdir[1]);
  printf("MP '%s' %i : Initialised storage for %i vacf values\n",
	 mp->name, mp->step, maxnvacf);
  printf("MP '%s' %i : ntracers = %f\n",
         mp->name, mp->step, ntracers);
  printf("MP '%s' %i : vacf(0) = %f\n",
         mp->name, mp->step, vacfval);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpvacfinit done\n", mp->name, mp->step);
  }
}


/* Alternative initialisation method */

void mpvacfinit2(MPSys mp, MPSys mp2) {
  int inode;
  LBReal t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14;
  LBReal rhot, u, wfac;
  LBReal *a, *p, *p2, *rho, *psi;
  LBObj *obj;
/* Initialise P(r,0) = u_a(r) psi(r), P_2(r,0) = psi(r) */
  obj = mp->lb->obj; a = mp->lb->a; p = mp->p; p2 = mp2->p;
  psi = mp->lb->psi; rho = mp->lb->rho;
  u = 0.0; wfac = 1.0;
  for (inode=0; inode<(mp->lb->nnodes); inode++) if (obj[inode] == NULL) {
    rhot = rho[inode];
    t01 = a[1+15*inode]; t02 = a[2+15*inode];
    t03 = a[3+15*inode]; t04 = a[4+15*inode]; t05 = a[5+15*inode];
    t06 = a[6+15*inode]; t07 = a[7+15*inode]; t08 = a[8+15*inode];
    t09 = a[9+15*inode]; t10 = a[10+15*inode]; t11 = a[11+15*inode];
    t12 = a[12+15*inode]; t13 = a[13+15*inode]; t14 = a[14+15*inode];
    switch (vacfl1) {
      case 0: u = (t01-t02+t07-t08+t09-t10+t11-t12+t13-t14)/rhot; break;
      case 1: u = (t03-t04+t07+t08-t09-t10+t11+t12-t13-t14)/rhot; break;
      case 2: u = (t05-t06+t07+t08+t09+t10-t11-t12-t13-t14)/rhot; break;
    }
    if (psi != NULL) wfac = psi[inode];
    p[inode] = wfac * u;
    p2[inode] = wfac;
  }
/* Don't need a special calculation for the vacf at t=0 */
  nvacf = 0;
  printf("MP '%s' %i : Re-initialised vacf calculation\n", mp->name, mp->step);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpvacfinit2 done\n", mp->name, mp->step);
  }
}


/* Calculate current vacf.  Note the presence of the adjoint
 * probability amplitude in all spatial averages.
 */

void mpvacf(MPSys mp, MPSys mp2) {
  int inode;
  LBReal t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14;
  LBReal rhot, u, wfac, vacftot, ntracers, vacfval, ptot, pmean, psitot;
  LBReal *a, *p, *p2, *rho, *psi;
  LBObj *obj;
/* Check there is enough storage space otherwise wasting effort */
  if (nvacf > maxnvacf) {
    fprintf(stderr,"Ran out of storage after %i vacf values\n", maxnvacf);
    return;
  }
  obj = mp->lb->obj; a = mp->lb->a; p = mp->p; p2 = mp2->p;
  psi = mp->lb->psi; rho = mp->lb->rho;
/* Make sure P(r) is orthogonal to psi(r), ie sum_r P(r) psi(r) = 0 */
  ptot = psitot = 0.0; wfac = 1.0;
  for (inode=0; inode<(mp->lb->nnodes); inode++) if (obj[inode] == NULL) {
    if (psi != NULL) wfac = psi[inode];
    ptot += wfac * p[inode]; psitot += wfac*wfac;
  }
  pmean = ptot / psitot;
  for (inode=0; inode<(mp->lb->nnodes); inode++) if (obj[inode] == NULL) {
    if (psi != NULL) wfac = psi[inode];
    p[inode] -= wfac*pmean;
  }
/* Calculate sum_r P(r) u_b(r) and sum_r P2(r) */
  ptot = vacftot = ntracers = 0.0; u = 0.0; wfac = 1.0;
  for (inode=0; inode<(mp->lb->nnodes); inode++) if (obj[inode] == NULL) {
    rhot = rho[inode];
    t01 = a[1+15*inode]; t02 = a[2+15*inode];
    t03 = a[3+15*inode]; t04 = a[4+15*inode]; t05 = a[5+15*inode];
    t06 = a[6+15*inode]; t07 = a[7+15*inode]; t08 = a[8+15*inode];
    t09 = a[9+15*inode]; t10 = a[10+15*inode]; t11 = a[11+15*inode];
    t12 = a[12+15*inode]; t13 = a[13+15*inode]; t14 = a[14+15*inode];
    switch (vacfl2) {
      case 0: u = (t01-t02+t07-t08+t09-t10+t11-t12+t13-t14)/rhot; break;
      case 1: u = (t03-t04+t07+t08-t09-t10+t11+t12-t13-t14)/rhot; break;
      case 2: u = (t05-t06+t07+t08+t09+t10-t11-t12-t13-t14)/rhot; break;
    }
    if (psi != NULL) wfac = psi[inode];
    ptot += wfac * p[inode];
    vacftot += wfac * p[inode] * u;
    ntracers += wfac * p2[inode];
  }
/* Compute and save the value of the vacf */
  vacfval = vacftot / ntracers;
  vacfarr[nvacf] = vacfval; nvacf++;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpvacf done\n", mp->name, mp->step);
  }
}


/* Write out vacf and partial sums, return final diffusion coefficient */

LBReal mpvacfwr(MPSys mp, MPSys mp2, char *vacffile, int nvals) {
  int i, ntrck;
  LBReal diffcon, vacfval;
  FILE *fp;
  if (!myscmp("null", vacffile)) {
    if ((fp = fopen(vacffile,"w")) == NULL) {
      fprintf(stderr,"%s could not be opened\n", vacffile);
      return 0.0;
    }
    printf("MP '%s' %i : Writing vacf data to %s\n",
	   mp->name, mp->step, vacffile);
  } else {
    fp = NULL;
    printf("MP '%s' %i : Analyzing vacf data\n", mp->name, mp->step);
  }
/* Now compute c_vv and other info */
  ntrck = nvacf / nvals;
  diffcon = 0.5 * vacfarr[0];
/* For t>0 simply correct vacf and sum to find diffusion constant */
  for (i=1; i<nvacf; i++) {
    vacfval = vacfarr[i];
    diffcon += vacfval;
    if (fp != NULL && i % ntrck == 0) {
      fprintf(fp,"%6i %g %g\n", i, vacfval, diffcon);
    }
    if (i > nvacf-4) {
      printf("MP '%s' %i : step, d = %8i %g %g\n",
	     mp->name, mp->step, i, vacfval, diffcon);
    }
  }
  fclose(fp);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP '%s' %i : mpvacfwr done\n", mp->name, mp->step);
  }
  return diffcon;
}

/* End of mpvacf.c */
