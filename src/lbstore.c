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

#define LBSTORE_C

#include <stdio.h>
#include <stdlib.h>
#include "sunlightlb.h"


/* Save density, momentum, stress field to a file in binary.
 * This is done just before the collision step as the collapse
 * to (density, momentum, stress) is also an integral part of
 * the collision step, thus no essential information is lost.
 * The routine therefore takes a bounceback and propagation step
 * before saving, and a collision step after saving.
 * Intended for use with the companion routine lbload.
 */

void lbsave(LBSys lb, char *savefile) {
  int inode;
  LBReal t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14;
  LBReal tlppp, tlmpm, tlpmm, tlmmp, thppp, thmpm, thpmm, thmmp, tlhp;
  LBReal t1p2, t3p4, t5p6, t0p6, sxyz, twoby3;
  LBReal *aa;
  LBStoreReal drho, rux, ruy, ruz, sxx, syy, szz, sxy, sxz, syz;
  FILE *fp;
  twoby3 = 2.0/3.0;
  if ((fp = fopen(savefile,"wb")) == NULL) {
    printf("savefile: %s could not be opened\n", savefile); return;
  }
  printf("\nSaving density, momentum, stress to %s\n", savefile);
/* Take a bounceback and a propagation step to get into correct position */
  bbl(lb); propagate(lb);
  aa = lb->a;
/* Save the density, momentum and stress on a per-node basis */
/* This is equivalent to traversing axx[i][j][k] arrays with k varying fastest */
  for (inode=0; inode<(lb->nnodes); inode++) {
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
/* Calculate density (actually drho = rho - rho0), momentum and stress */
    drho = t0p6 + tlhp;
    rux = t01 - t02 + tlmpm + thmpm;
    ruy = t03 - t04 + tlpmm + thpmm;
    ruz = t05 - t06 + tlppp - thppp;
    sxyz = twoby3*drho - t0p6;
    sxx = sxyz + t1p2;    syy = sxyz + t3p4;    szz = sxyz + t5p6;
    sxy = tlmmp + thmmp;  sxz = tlmpm - thmpm;  syz = tlpmm - thpmm;
/* Write out data in binary */
    fwrite(&drho, sizeof(LBStoreReal),1, fp);
    fwrite(&rux, sizeof(LBStoreReal),1, fp);
    fwrite(&ruy, sizeof(LBStoreReal),1, fp);
    fwrite(&ruz, sizeof(LBStoreReal),1, fp);
    fwrite(&sxx, sizeof(LBStoreReal),1, fp);
    fwrite(&syy, sizeof(LBStoreReal),1, fp);
    fwrite(&szz, sizeof(LBStoreReal),1, fp);
    fwrite(&sxy, sizeof(LBStoreReal),1, fp);
    fwrite(&sxz, sizeof(LBStoreReal),1, fp);
    fwrite(&syz, sizeof(LBStoreReal),1, fp);
  }
  fclose(fp);
/* Take a collision step to get back into correct position in cycle */
  collide(lb);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : lbsave done\n", lb->step);
}


/* Load density, momentum, stress field from a file in binary,
 * and construct occupation densities appropriately.
 * The load is pre-collision step as the collapse to (density,
 * momentum, stress) is an integral part of the collision step.
 * Therefore a collision step is taken at the end.
 * Intended for use with the companion routine lbsave.
 * The parameter 'lb->rescale' is used to multiply the velocity field.
 * It is *reset* to 1.0 at the routine end.
 */

void lbload(LBSys lb, char *loadfile) {
  int i, inode;
  LBReal ruxyz, sxyz, trsby3, b, rhobyb, erhobyb;
  LBReal dbybc2, edbybc2, dd2bybc4, fdd2bybc4;
  LBReal *aa;
  LBStoreReal drho, rux, ruy, ruz, sxx, syy, szz, sxy, sxz, syz;
  FILE *fp;
  LBObj *obj;
/* Synchronise objects */
  if (lb->changed_objects) syncobjs(lb);
/* Generate some constant parameters */
  b = 72.0; dbybc2 = 1.0/24.0; dd2bybc4 = 1.0/8.0;
  fdd2bybc4 = 4.0*dd2bybc4; edbybc2 = 8.0*dbybc2;
/* Initialise LB either from file or to uniform state */
  if ((fp = fopen(loadfile,"r")) == NULL) {
    printf("lbload: %s could not be opened\n", loadfile);
    exit(EXIT_FAILURE);
  }
  if (lb->load_average) {
    printf("\nAveraging density, momentum, stress with %s\n", loadfile);
  } else {
    printf("\nLoading density, momentum, stress from %s\n", loadfile);
  }
  if (lb->rescale != 1.0) {
    printf("Multiplying velocity field (only) by %f\n", lb->rescale);
    printf("*** this factor will be reset to unity at end of loading ***\n");
  }
  aa = lb->a; obj = lb->obj;
/* Load density, momentum and stress on a per-node basis */
  for (inode=0; inode<(lb->nnodes); inode++) {
    fread(&drho, sizeof(LBStoreReal),1, fp);
    fread(&rux, sizeof(LBStoreReal),1, fp);
    fread(&ruy, sizeof(LBStoreReal),1, fp);
    fread(&ruz, sizeof(LBStoreReal),1, fp);
    fread(&sxx, sizeof(LBStoreReal),1, fp);
    fread(&syy, sizeof(LBStoreReal),1, fp);
    fread(&szz, sizeof(LBStoreReal),1, fp);
    fread(&sxy, sizeof(LBStoreReal),1, fp);
    fread(&sxz, sizeof(LBStoreReal),1, fp);
    fread(&syz, sizeof(LBStoreReal),1, fp);
/* Rescale velocity components */
    rux = (lb->rescale)*rux; ruy = (lb->rescale)*ruy; ruz = (lb->rescale)*ruz;
/* Rebuild occupation probabilities: a_i = w_i (f_i - rho_0 / b) */
    rhobyb = drho/b;  erhobyb = 8.0*rhobyb;
    trsby3 = (sxx+syy+szz)/3.0;
    if (lb->load_average) {
      aa[ 0+15*inode] += 0.5 * (2.0*erhobyb - 2.0*fdd2bybc4 * trsby3);
      ruxyz = edbybc2*rux;  sxyz = fdd2bybc4*(sxx - trsby3);
      aa[ 1+15*inode] += 0.5 * (erhobyb + ruxyz + sxyz);
      aa[ 2+15*inode] += 0.5 * (erhobyb - ruxyz + sxyz);
      ruxyz = edbybc2*ruy;  sxyz = fdd2bybc4*(syy - trsby3);
      aa[ 3+15*inode] += 0.5 * (erhobyb + ruxyz + sxyz);
      aa[ 4+15*inode] += 0.5 * (erhobyb - ruxyz + sxyz);
      ruxyz = edbybc2*ruz;  sxyz = fdd2bybc4*(szz - trsby3);
      aa[ 5+15*inode] += 0.5 * (erhobyb + ruxyz + sxyz);
      aa[ 6+15*inode] += 0.5 * (erhobyb - ruxyz + sxyz);
      ruxyz = dbybc2*(rux + ruy + ruz);
      sxyz = dd2bybc4*(trsby3 + sxy + sxz + syz);
      aa[ 7+15*inode] += 0.5 * (rhobyb + ruxyz + sxyz);
      aa[14+15*inode] += 0.5 * (rhobyb - ruxyz + sxyz);
      ruxyz = dbybc2*(rux - ruy - ruz);
      sxyz = dd2bybc4*(trsby3 - sxy - sxz + syz);
      aa[ 8+15*inode] += 0.5 * (rhobyb - ruxyz + sxyz);
      aa[13+15*inode] += 0.5 * (rhobyb + ruxyz + sxyz);
      ruxyz = dbybc2*(rux - ruy + ruz);
      sxyz = dd2bybc4*(trsby3 - sxy + sxz - syz);
      aa[ 9+15*inode] += 0.5 * (rhobyb + ruxyz + sxyz);
      aa[12+15*inode] += 0.5 * (rhobyb - ruxyz + sxyz);
      ruxyz = dbybc2*(rux + ruy - ruz);
      sxyz = dd2bybc4*(trsby3 + sxy - sxz - syz);
      aa[10+15*inode] += 0.5 * (rhobyb - ruxyz + sxyz);
      aa[11+15*inode] += 0.5 * (rhobyb + ruxyz + sxyz);
    } else {
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
  }
  fclose(fp);
/* Explicitly reset the occupancies in non-fluid regions */
/* if relaxation on the solid nodes is not done. */
  if (lb->skip_solid) {
    for (inode=0; inode<(lb->nnodes); inode++) {
      if (obj[inode] == NULL) continue;
      for (i=0; i<15; i++) aa[i+15*inode] = 0.0;
    }
  }
/* Finally, take a collision step */
  lb->status = POST_PROPAGATE; collide(lb);
  lb->rescale = 1.0;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : lbload done\n", lb->step);
}


/* Write velocity field in ascii or binary to a file, selected by
 * type.  A simulation can be approximately restarted with this data
 * too, by the companion routine lbloadvel below.  This routine uses
 * double precision for improved accuracy.
 */

void lbsavevel(LBSys lb, char *savefile, char *type) {
  int i, j, k, inode;  
  int esi, esj, esk;
  int saveasbinary = 0;
  int skip_solid;
  LBReal t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14;
  LBReal rho, fac = 0.0;
  LBReal rho0;
  LBReal *aa;
  LBReal *fx, *fy, *fz;
  LBReal ux, uy, uz;
  FILE *fp;
  LBObj *obj;
  esi = lb->esi; esj = lb->esj; esk = lb->esk; 
  if (type[0] == 'b' || type[0] == 'B') saveasbinary = 1;
  if ((fp = fopen(savefile,"w")) == NULL) {
    printf("savefile: %s could not be opened\n", savefile); return;
  }
  rho0 = lb->rho0; obj = lb->obj; skip_solid = lb->skip_solid;
  fac = 0.0;
  if (lb->vel_corrected) {
    printf("\nSaving corrected velocity field in %s to %s",
	   saveasbinary ? "binary" : "ascii", savefile);
    switch (lb->status) {
      case POST_BBL: printf(" [post-bounceback - not a good idea]\n"); break;
      case POST_PROPAGATE: fac = 0.5; printf("\n"); break;
      case POST_COLLIDE: fac = -0.5; printf("\n"); break;
    }
  } else {
    printf("\nSaving ");
    switch (lb->status) {
      case POST_BBL: printf("post-bounceback [not a good idea]"); break;
      case POST_PROPAGATE: printf("post-propagation"); break;
      case POST_COLLIDE: printf("post-collision"); break;
    }
    printf(" uncorrected velocity field in %s to %s\n",
	   saveasbinary ? "binary" : "ascii", savefile);
  }
/* Compute and save velocity on a per-node basis */
/* This is equivalent to traversing v[i][j][k] arrays with k varying fastest */
  aa = lb->a; fx = lb->fx; fy = lb->fy; fz = lb->fz;
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
/* Explicitly set velocity on the solid nodes if relaxation is not done. */
    if ((obj[inode] != NULL) && skip_solid) {
      ux = uy = uz = 0.0;
    } else {
/* Calculate density and velocity, corrected for added force */            
      t00 = aa[ 0+15*inode]; t01 = aa[ 1+15*inode]; t02 = aa[ 2+15*inode];
      t03 = aa[ 3+15*inode]; t04 = aa[ 4+15*inode]; t05 = aa[ 5+15*inode];
      t06 = aa[ 6+15*inode]; t07 = aa[ 7+15*inode]; t08 = aa[ 8+15*inode];
      t09 = aa[ 9+15*inode]; t10 = aa[10+15*inode]; t11 = aa[11+15*inode];
      t12 = aa[12+15*inode]; t13 = aa[13+15*inode]; t14 = aa[14+15*inode];
      rho = rho0 + t00+t01+t02+t03+t04+t05+t06+t07+t08+t09+t10+t11+t12+t13+t14;
      ux = (t01-t02+t07-t08+t09-t10+t11-t12+t13-t14 + fac*fx[inode])/rho;
      uy = (t03-t04+t07+t08-t09-t10+t11+t12-t13-t14 + fac*fy[inode])/rho;
      uz = (t05-t06+t07+t08+t09+t10-t11-t12-t13-t14 + fac*fz[inode])/rho;
    }
    if (saveasbinary) {
      fwrite(&rho, sizeof(LBReal),1, fp);
      fwrite(&ux, sizeof(LBReal),1, fp);
      fwrite(&uy, sizeof(LBReal),1, fp);
      fwrite(&uz, sizeof(LBReal),1, fp);
    } else {
      fprintf(fp,"%g %g %g %g %i %i %i\n", rho, ux, uy, uz, i, j, k);
    }
  }
  fclose(fp);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : lbsavevel done\n", lb->step);
}


/* Load velocity field from a velocity file, and construct
 * occupation densities approximately.  Intended for use with the
 * companion routine lbsavevel.  The LBReal parameter 'lb->rescale'
 * multiplies the velocity field.  It is *reset* to 1.0 at the routine end.
 */

void lbloadvel(LBSys lb, char *loadfile, char *type) {
  int i, j, k, inode;
  int stokes_flow;
  int loadfrombinary = 0;
  float tr, tx, ty, tz;
  LBReal ruxyz, sxyz, trsby3, b, rhobyb, erhobyb;
  LBReal dbybc2, edbybc2, dd2bybc4, fdd2bybc4;
  LBReal rho0, rho, drho, rux, ruy, ruz, sxx, syy, szz, sxy, sxz, syz;
  LBReal *aa;
  LBReal ux, uy, uz;
  FILE *fp;
  LBObj *obj;
  if (type[0] == 'b' || type[0] == 'B') loadfrombinary = 1;
/* Synchronise objects */
  if (lb->changed_objects) syncobjs(lb);
/* Generate some constant parameters */
  b = 72.0; dbybc2 = 1.0/24.0; dd2bybc4 = 1.0/8.0;
  fdd2bybc4 = 4.0*dd2bybc4; edbybc2 = 8.0*dbybc2;
/* Initialise LB either from file or to uniform state */
  if (myscmp("null", loadfile)) {
    printf("\nNo load file\n"); fp = NULL; return;
  }
  printf("\nLoading velocity field from %s file %s\n",
	 loadfrombinary ? "binary" : "ascii", loadfile);
  printf(" (approximate LB initialisation)\n");
  if ((fp = fopen(loadfile,"r")) == NULL) {
    printf("lbloadvel: %s could not be opened\n", loadfile);
    exit(EXIT_FAILURE);
  }
  if (lb->rescale != 1.0) {
    printf("Multiplying velocity field by %f\n", lb->rescale);
    printf("*** this factor will be reset to unity at end of loading ***\n");
  }
  aa = lb->a; obj = lb->obj; stokes_flow = lb->stokes_flow; rho0 = lb->rho0;
/* Load velocity on a per-node basis and determine occupation probabilities */
  for (inode=0; inode<(lb->nnodes); inode++) {
    if (loadfrombinary) {
      fread(&rho, sizeof(LBReal),1, fp);
      fread(&ux, sizeof(LBReal),1, fp);
      fread(&uy, sizeof(LBReal),1, fp);
      fread(&uz, sizeof(LBReal),1, fp);
    } else {
      fscanf(fp,"%f%f%f%f%i%i%i",&tr,&tx,&ty,&tz,&i,&j,&k);
      rho = tr; ux = tx; uy = ty; uz = tz;
    }
/* Rescale velocity fields */
    ux = (lb->rescale)*ux; uy = (lb->rescale)*uy; uz = (lb->rescale)*uz;
    /* Rconstruct density, momentum and stress */
    drho = rho - rho0;
    rux = rho*ux; ruy = rho*uy; ruz = rho*uz;
    if (stokes_flow) {
      sxx = sxy = sxz = syy = syz = szz = 0.0;
    } else {
      sxx = ux*rux; sxy = ux*ruy; sxz = ux*ruz;
      syy = uy*ruy; syz = uy*ruz; szz = uz*ruz;
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
  fclose(fp);
/* Explicitly reset the occupancies in non-fluid regions */
/* if relaxation on the solid nodes is not done. */
  if (lb->skip_solid) {
    for (inode=0; inode<(lb->nnodes); inode++) {
      if (obj[inode] == NULL) continue;
      for (i=0; i<15; i++) aa[i+15*inode] = 0.0;
    }
  }
/* Finally, take a collision step */
  lb->status = POST_PROPAGATE; collide(lb);
  lb->rescale = 1.0;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : lbloadvel done\n", lb->step);
}


/* Save stress field to a file in ascii or binary.  The average of the
 * post-propagation and post-collision stress is saved, thus this
 * routine takes one lb step.
 */

void lbsavestr(LBSys lb, char *savefile, char *type) {
  int inode;
  int saveasbinary = 0;
  int skip_solid, void_space;
  LBReal t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14;
  LBReal tlppp, tlmpm, tlpmm, tlmmp, thppp, thmpm, thpmm, thmmp, tlhp;
  LBReal t1p2, t3p4, t5p6, t0p6, sxyz, twoby3, drho;
  LBReal *aa;
  LBReal *sxx, *syy, *szz, *sxy, *sxz, *syz;
  LBReal txx, tyy, tzz, txy, txz, tyz;
  LBObj *obj;
  FILE *fp;
  twoby3 = 2.0/3.0;
  aa = lb->a;
  obj = lb->obj; skip_solid = lb->skip_solid;
/* Obtain temporary stress field storage */
  sxx = (LBReal *) malloc((lb->nnodes)*sizeof(LBReal));
  syy = (LBReal *) malloc((lb->nnodes)*sizeof(LBReal));
  szz = (LBReal *) malloc((lb->nnodes)*sizeof(LBReal));
  sxy = (LBReal *) malloc((lb->nnodes)*sizeof(LBReal));
  sxz = (LBReal *) malloc((lb->nnodes)*sizeof(LBReal));
  syz = (LBReal *) malloc((lb->nnodes)*sizeof(LBReal));
  if (sxx == NULL || syy == NULL || szz == NULL
      || sxy == NULL || sxz == NULL || szz == NULL) {
    printf("lbsavestr: ran out of memory for temporary storage\n"); return;
  }    
/* Take a bounceback and a propagation step */
  bbl(lb); propagate(lb);
/* Sample stress field */
  for (inode=0; inode<(lb->nnodes); inode++) {
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
/* Calculate density (actually drho = rho - rho0) and stress */
    drho = t0p6 + tlhp;
    sxyz = twoby3*drho - t0p6;
    sxx[inode] = sxyz + t1p2;    syy[inode] = sxyz + t3p4;    szz[inode] = sxyz + t5p6;
    sxy[inode] = tlmmp + thmmp;  sxz[inode] = tlmpm - thmpm;  syz[inode] = tlpmm - thpmm;
  }
/* Take a collision step */
  collide(lb);
/* Sample stress field again */
  for (inode=0; inode<(lb->nnodes); inode++) {
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
/* Calculate density (actually drho = rho - rho0) and stress */
    drho = t0p6 + tlhp;
    sxyz = twoby3*drho - t0p6;
    sxx[inode] += sxyz + t1p2;    syy[inode] += sxyz + t3p4;    szz[inode] += sxyz + t5p6;
    sxy[inode] += tlmmp + thmmp;  sxz[inode] += tlmpm - thmpm;  syz[inode] += tlpmm - thpmm;
  }
/* Save stress field */
  if (type[0] == 'b' || type[0] == 'B') saveasbinary = 1;
  if ((fp = fopen(savefile,"wb")) == NULL) {
    printf("savefile: %s could not be opened\n", savefile); return;
  }
  printf("\nSaving stress field in %s to %s\n",
	 saveasbinary ? "binary" : "ascii", savefile);
  for (inode=0; inode<(lb->nnodes); inode++) {
/* Explicitly set velocity on the solid nodes if relaxation is not done. */
    if ((obj[inode] != NULL) && skip_solid) {
      txx = tyy = tzz = txy = txz = tyz = 0.0;
    } else {
      txx = 0.5*sxx[inode]; tyy = 0.5*syy[inode]; tzz = 0.5*szz[inode];
      txy = 0.5*sxy[inode]; txz = 0.5*sxz[inode]; tyz = 0.5*syz[inode];
    }
    if (obj[inode] == NULL) void_space = 1; else void_space = 0;
    if (saveasbinary) {
      fwrite(&txx, sizeof(LBStoreReal), 1, fp);
      fwrite(&tyy, sizeof(LBStoreReal), 1, fp);
      fwrite(&tzz, sizeof(LBStoreReal), 1, fp);
      fwrite(&txy, sizeof(LBStoreReal), 1, fp);
      fwrite(&txz, sizeof(LBStoreReal), 1, fp);
      fwrite(&tyz, sizeof(LBStoreReal), 1, fp);
      fwrite(&void_space, sizeof(int), 1, fp);
    } else {
      fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\t%i\n", 
	      txx, tyy, tzz, txy, txz, tyz, void_space);
    }
  }
  fclose(fp);
/* Release temporary stress field storage */
  free(sxx); free(syy); free(szz); free(sxy); free(sxz); free(syz);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : lbsavestr done\n", lb->step);
}


/* Write object prositions, velocities and forces in ascii to a file.
 * These can then be analysed offline, or restored with the companion
 * routine lbloadobj below.
 */

void lbsaveobj(LBSys lb, char *savefile) {
  int moreinfo = 0;
  FILE *fp;
  LBObj p;
  if (myscmp("stdout", savefile)) { fp = stdout; moreinfo = 1; }
  else if (myscmp("stderr", savefile)) { fp = stderr; moreinfo = 1; }
  else if ((fp = fopen(savefile,"w")) == NULL) {
    printf("lbsaveobj: %s could not be opened\n", savefile); return;
  }
  printf("\nSaving object data to %s\n", savefile);
  if (moreinfo) fprintf(fp,"START of particle data: TYPE, id, base id, x, y, z\n\n");
  for (p=lb->obj_head; p!=NULL; p=p->next_obj) {
    fprintf(fp,"POS %3i %3i %f %f %f\n", p->id, p->base_obj->id, p->posx, p->posy, p->posz);
    fprintf(fp,"VEL %3i %3i %f %f %f\n", p->id, p->base_obj->id, p->velx, p->vely, p->velz);
    fprintf(fp,"ANG %3i %3i %f %f %f\n", p->id, p->base_obj->id, p->angx, p->angy, p->angz);
    fprintf(fp,"FRC %3i %3i %f %f %f\n", p->id, p->base_obj->id, p->frcx, p->frcy, p->frcz);
    fprintf(fp,"TRQ %3i %3i %f %f %f\n", p->id, p->base_obj->id, p->trqx, p->trqy, p->trqz);
    if (moreinfo) fprintf(fp,"\n");
  }
  if (moreinfo) fprintf(fp,"END of particle data\n");
  fclose(fp);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : lbsaveobj done\n", lb->step);
}


/* Read object positions and velocities in ascii from a file.
 * The factor 'lb->rescale' is used to multiply velocities and angular
 * velocities.  It is *reset* to unity at the end of the routine.
 * Assumes that object map and object data structures are all initialised.
 */

void lbloadobj(LBSys lb, char *loadfile) {
  int irc = 0;
  int i1, i2;
  float t1, t2, t3;
  FILE *fp;
  LBObj p1, p2;
  if ((fp = fopen(loadfile,"r")) == NULL) {
    printf("lbloadobj: %s could not be opened\n", loadfile); return;
  } else printf("Loading object data from %s\n", loadfile);
  if (lb->rescale != 1.0) {
    printf("Multiplying velocities and angular velocities by %f\n", lb->rescale);
    printf("*** this factor will be reset to unity at end of loading ***\n");
  }
  while (fscanf(fp,"%*s%i%i%f%f%f",&i1,&i2,&t1,&t2,&t3) == 5) {
    if ((p1 = lookupobj(lb, i1)) == NULL) {irc = 1; break; }
    if (i2 != i1) {
      if ((p2 = lookupobj(lb, i2)) == NULL) { irc = 2; break; }
      p1->base_obj->id = i2;
    }
    p1->posx = t1; p1->posy = t2; p1->posz = t3;
    if (fscanf(fp,"%*s%i%i%f%f%f",&i1,&i2,&t1,&t2,&t3) != 5) { irc = 3; break; }
    p1->velx = (lb->rescale)*t1; p1->vely = (lb->rescale)*t2; p1->velz = (lb->rescale)*t3;
    if (fscanf(fp,"%*s%i%i%f%f%f",&i1,&i2,&t1,&t2,&t3) != 5) { irc = 4; break; }
    p1->angx = (lb->rescale)*t1; p1->angy = (lb->rescale)*t2; p1->angz = (lb->rescale)*t3;
    if (fscanf(fp,"%*s%i%i%f%f%f",&i1,&i2,&t1,&t2,&t3) != 5) { irc = 5; break; }
    if (fscanf(fp,"%*s%i%i%f%f%f",&i1,&i2,&t1,&t2,&t3) != 5) { irc = 6; break; }
  }
  if (irc == 1 || irc == 2) printf("UFO found (unidentified floating object)\n");
  if (irc > 3) printf("Ran out of data apparently\n");
  fclose(fp);
  lb->rescale = 1.0;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : lbloadobj done\n", lb->step);
}


/*
 * This writes out the instantaneous node force data, and
 * should be called immediately after measuring the forces.
 * Type can be "binary", "ascii" or "list".
 */

void lbsavefrc(LBSys lb, char *frcfile, char *type) {
  int i, j, k, esi, esj, esk, inode;
  int savebinary = 0;
  int saveaslist = 0;
  int *mask;
  LBReal *fx, *fy, *fz;
  LBStoreReal forcex, forcey, forcez;
  LBObj p;
  FILE *fp;
  if (type[0] == 'b' || type[0] == 'B') savebinary = 1;
  if (type[0] == 'l' || type[0] == 'L') saveaslist = 1;
  if ((fp = fopen(frcfile,"w")) == NULL) {
    printf("lbsavefrc: %s could not be opened\n", frcfile);
    exit(EXIT_FAILURE);
  }
  fx = lb->fx; fy = lb->fy; fz = lb->fz;
  mask = lb->mask; esi = lb->esi; esj = lb->esj; esk = lb->esk; 
  if (saveaslist) {
    printf("\nSaving force data as a list to %s\n", frcfile);
    for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
      inode = esk*(esj*i+j)+k;
      if (mask[inode] > 0) {
	p = lb->obj[inode];
	fprintf(fp,"%i %i %i %f %f %f %i %f %f %f\n", 
		i, j, k, fx[inode], fy[inode], fz[inode], 
		p->id, p->posx, p->posy, p->posz);
      }
    }
  } else {
    printf("\nSaving %s force data to %s\n", savebinary?"binary":"ascii", frcfile);
/* Write out forces acting on individual nodes on a per-node basis */
/* This is equivalent to traversing f[i][j][k] arrays with k varying fastest */
    for (inode=0; inode<(lb->nnodes); inode++) {
      forcex = fx[inode]; forcey = fy[inode]; forcez = fz[inode];
      if (savebinary) {
	fwrite(&forcex, sizeof(LBStoreReal),1, fp);
	fwrite(&forcey, sizeof(LBStoreReal),1, fp);
	fwrite(&forcez, sizeof(LBStoreReal),1, fp);
      } else {
	fprintf(fp,"%f %f %f\n", forcex, forcey, forcez);
      }
    }
  }
  fclose(fp);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : lbsavefrc done\n", lb->step);
}


/* This writes the solid regions to a mapfile, either as a list or
 * as a character array (substrate file).
 */

void lbsavemap(LBSys lb, char *mapfile, char *type) {
  int i, j, k, esi, esj, esk, inode;
  int saveaslist = 0;
  FILE *fp;
  LBObj *obj;
  if ((fp = fopen(mapfile,"w")) == NULL) {
    printf("lbsavemap: %s could not be opened\n", mapfile);
    exit(EXIT_FAILURE);
  }
  if (type[0] == 'l' || type[0] == 'L') saveaslist = 1;
  obj = lb->obj; esi = lb->esi; esj = lb->esj; esk = lb->esk; 
  if (saveaslist) {
    printf("\nSaving list of solid nodes to %s\n", mapfile);
    for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
      inode = esk*(esj*i+j)+k;
      if (obj[inode] != NULL) fprintf(fp,"%i %i %i\n", i, j, k);
    }
  } else {
    printf("\nSaving object map as character array to %s\n", mapfile);
/* Write out map acting on individual nodes on a per-node basis */
/* This is equivalent to traversing obj[i][j][k] with k varying fastest */
    for (inode=0; inode<(lb->nnodes); inode++) {
      if (obj[inode] == NULL) {
	fprintf(fp,"%c",(char)0);
      } else {
	fprintf(fp,"%c",(char)255);
      }
    }
  }
  fclose(fp);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : lbsavemap done\n", lb->step);
}


/* End of lbstore.c */
