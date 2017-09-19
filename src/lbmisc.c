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

#define LBMISC_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "sunlightlb.h"

/* Local variables and routines */

static time_t tstart;

static LBReal *ralloc(int);
static int *ialloc(int);
static void wrhms(FILE *, char *, int);

/* Routine used when there is a non-recoverable error */

void error(char *s1) {
  fprintf(stderr,"Opps!  %s\n", s1);
  exit(EXIT_FAILURE); 
}


/* Start the clock */

void startclock() {
  if (tstart == 0) printf("Clock started\n");
  else printf("Clock restarted\n");
  tstart = time(NULL); 
}

/* Set the verbosity level */

void set_verbosity(int what) { 
  slb_verbosity |= what; 
  report_verbosity();
  return;
}

void unset_verbosity(int what) { 
  if (slb_verbosity & what) slb_verbosity -= what; 
  report_verbosity();
  return;
}

void report_verbosity() {
  printf("Verbosity =");
  if (slb_verbosity == 0) { printf(" NONE\n"); return; }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf(" TRACE_ROUTINE_EXITS");
  if (slb_verbosity & OBJECT_DETAILS) printf(" OBJECT_DETAILS");
  if (slb_verbosity & FORCE_DETAILS) printf(" FORCE_DETAILS");
  printf("\n");
  return;
}

/* Deliver up a newly initialised lattice Boltzmann simulation system,
 * this is a hook for stuff that needs to be done right at the start.
 */

LBSys newlb(int reqesi, int reqesj, int reqesk, char *name) {
  int i, inode;
  LBSys lb;
/* Make the new LB system and allocate variables */
  lb = (LBSys) malloc(sizeof(LBSys_type));
  if (lb == NULL) error("No room for lbsys structure\n");
/* These are apropriate for D3Q15(r=16) model */
/* Note that i'th component of velocity for link l is c[3*l+i] */
  lb->c = ialloc(45); lb->w = ialloc(15);
  lb->c[ 0] =  0; lb->c[ 1] =  0; lb->c[ 2] =  0; lb->w[ 0] = 16;
  lb->c[ 3] =  1; lb->c[ 4] =  0; lb->c[ 5] =  0; lb->w[ 1] =  8;
  lb->c[ 6] = -1; lb->c[ 7] =  0; lb->c[ 8] =  0; lb->w[ 2] =  8;
  lb->c[ 9] =  0; lb->c[10] =  1; lb->c[11] =  0; lb->w[ 3] =  8;
  lb->c[12] =  0; lb->c[13] = -1; lb->c[14] =  0; lb->w[ 4] =  8;
  lb->c[15] =  0; lb->c[16] =  0; lb->c[17] =  1; lb->w[ 5] =  8;
  lb->c[18] =  0; lb->c[19] =  0; lb->c[20] = -1; lb->w[ 6] =  8;
  lb->c[21] =  1; lb->c[22] =  1; lb->c[23] =  1; lb->w[ 7] =  1;
  lb->c[24] = -1; lb->c[25] =  1; lb->c[26] =  1; lb->w[ 8] =  1;
  lb->c[27] =  1; lb->c[28] = -1; lb->c[29] =  1; lb->w[ 9] =  1;
  lb->c[30] = -1; lb->c[31] = -1; lb->c[32] =  1; lb->w[10] =  1;
  lb->c[33] =  1; lb->c[34] =  1; lb->c[35] = -1; lb->w[11] =  1;
  lb->c[36] = -1; lb->c[37] =  1; lb->c[38] = -1; lb->w[12] =  1;
  lb->c[39] =  1; lb->c[40] = -1; lb->c[41] = -1; lb->w[13] =  1;
  lb->c[42] = -1; lb->c[43] = -1; lb->c[44] = -1; lb->w[14] =  1;
  lb->esi = reqesi; lb->esj = reqesj; lb->esk = reqesk;
  lb->nnodes = (lb->esi)*(lb->esj)*(lb->esk);
  lb->a = ralloc(15*(lb->nnodes));
  lb->mask = ialloc(lb->nnodes);
  lb->obj = (LBObj *) malloc((lb->nnodes)*sizeof(LBObj));
  if (lb->obj == NULL) error("No room for LBObj array\n");
  for (i=0; i<(lb->nnodes); i++) { lb->obj[i] = NULL; lb->mask[i] = 0; }
  for (i=0; i<15*(lb->nnodes); i++) lb->a[i] = 0.0;
  lb->ux = ralloc(lb->nnodes); lb->uy = ralloc(lb->nnodes); lb->uz = ralloc(lb->nnodes);
  lb->fx = ralloc(lb->nnodes); lb->fy = ralloc(lb->nnodes); lb->fz = ralloc(lb->nnodes);
  lb->bx = ralloc(lb->nnodes); lb->by = ralloc(lb->nnodes); lb->bz = ralloc(lb->nnodes);
  for (inode=0; inode<(lb->nnodes); inode++) lb->ux[inode] = lb->uy[inode] = lb->uz[inode] = 0.0;
  for (inode=0; inode<(lb->nnodes); inode++) lb->fx[inode] = lb->fy[inode] = lb->fz[inode] = 0.0;
  for (inode=0; inode<(lb->nnodes); inode++) lb->bx[inode] = lb->by[inode] = lb->bz[inode] = 0.0;
  lb->tau = 0.75; 
  lb->rho0 = 0.0; for (i=0; i<15; i++) lb->rho0 += lb->w[i];
  lb->objects_moving = 0;
  lb->skip_solid = 0;
  lb->stokes_flow = 0;
  lb->vel_corrected = 1;
  lb->step = 0;
  lb->next_id = 0; lb->obj_head = NULL; 
  lb->status = POST_COLLIDE;
  lb->changed_objects = 0;
/* Initialise rescale and load_average used in storage routines */
  lb->rescale = 1.0;
  lb->load_average = 0;
/* Mean velocity weight is absent at first */
  lb->psi = NULL;
/* Density field is also only calculated if required */
  lb->rho = NULL;
/* Golbal variable - minimal verbosity at first */
  slb_verbosity = 0;
/* Global variable - clock starting time */
  tstart = 0;
/* Finish off with a few housekeeping items */
  if (myscmp("null", name)) lb->name = NULL; else lb->name = mycopy(name);
  lb->monitoring = (lb->name == NULL) ? 0 : 1;
  printf("LB system '%s' initialised for %i * %i * %i lattice\n", 
	 lb->name, lb->esi, lb->esj, lb->esk);
#ifdef SIGNAL
  siginit();
  printf("Signalling implemented, process id = %i\n",(int)getpid());
#else
  printf("Signalling not implemented\n");
#endif /* SIGNAL */
  return lb;
}

/* Writes out LB info */

void about() {
  printf("\nSunlightLB - a 3D Lattice Boltzmann code - version 1.0\n");
  printf("Copyright (C) 2005 Unilever UK Central Resources Ltd.\n\n");
  printf("SunlightLB is free software; you can redistribute it and/or\n");
  printf("modify it under the terms of the GNU General Public License\n");
  printf("as published by the Free Software Foundation; either version 2\n");
  printf("of the License, or (at your option) any later version.\n\n");
  printf("SunlightLB is distributed in the hope that it will be useful,\n");
  printf("but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
  printf("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
  printf("GNU General Public License for more details.\n\n");
}

/* Write out information on link weights and velocities, as a check. */

void wrlinks(LBSys lb) {
  printf(" w00 = %2i  c00 = (%2i,%2i,%2i)\n", lb->w[ 0], lb->c[ 0], lb->c[ 1], lb->c[ 2]);
  printf(" w01 = %2i  c01 = (%2i,%2i,%2i)\n", lb->w[ 1], lb->c[ 3], lb->c[ 4], lb->c[ 5]);
  printf(" w02 = %2i  c02 = (%2i,%2i,%2i)\n", lb->w[ 2], lb->c[ 6], lb->c[ 7], lb->c[ 8]);
  printf(" w03 = %2i  c03 = (%2i,%2i,%2i)\n", lb->w[ 3], lb->c[ 9], lb->c[10], lb->c[11]);
  printf(" w04 = %2i  c04 = (%2i,%2i,%2i)\n", lb->w[ 4], lb->c[12], lb->c[13], lb->c[14]);
  printf(" w05 = %2i  c05 = (%2i,%2i,%2i)\n", lb->w[ 5], lb->c[15], lb->c[16], lb->c[17]);
  printf(" w06 = %2i  c06 = (%2i,%2i,%2i)\n", lb->w[ 6], lb->c[18], lb->c[19], lb->c[20]);
  printf(" w07 = %2i  c07 = (%2i,%2i,%2i)\n", lb->w[ 7], lb->c[21], lb->c[22], lb->c[23]);
  printf(" w08 = %2i  c08 = (%2i,%2i,%2i)\n", lb->w[ 8], lb->c[24], lb->c[25], lb->c[26]);
  printf(" w09 = %2i  c09 = (%2i,%2i,%2i)\n", lb->w[ 9], lb->c[27], lb->c[28], lb->c[29]);
  printf(" w10 = %2i  c10 = (%2i,%2i,%2i)\n", lb->w[10], lb->c[30], lb->c[31], lb->c[32]);
  printf(" w11 = %2i  c11 = (%2i,%2i,%2i)\n", lb->w[11], lb->c[33], lb->c[34], lb->c[35]);
  printf(" w12 = %2i  c12 = (%2i,%2i,%2i)\n", lb->w[12], lb->c[36], lb->c[37], lb->c[38]);
  printf(" w13 = %2i  c13 = (%2i,%2i,%2i)\n", lb->w[13], lb->c[39], lb->c[40], lb->c[41]);
  printf(" w14 = %2i  c14 = (%2i,%2i,%2i)\n", lb->w[14], lb->c[42], lb->c[43], lb->c[44]);
  printf("\n");
}


/* Print out parameters for information.  */

void wrpars(LBSys lb) {
  int inode, nobjs, nfluid, msize;
  double fxsum, fysum, fzsum;
  LBReal nu, eta;
  LBObj p;
/* Count number of objects */
  nobjs = 0;
  for (p=lb->obj_head; p!=NULL; p=p->next_obj) nobjs++;
/* Count fluid nodes and the expected force */
  fxsum = fysum = fzsum = 0.0; nfluid = 0;
  for (inode=0; inode<(lb->nnodes); inode++) {
    if (lb->obj[inode] == NULL) {
      fxsum += lb->bx[inode]; fysum += lb->by[inode]; fzsum += lb->bz[inode]; nfluid++;
    }
  }
/* Calculate viscosity */
  nu = ((lb->tau)-0.5)/3.0; eta = (lb->rho0)*nu;
/* Write out parameters */
  printf("---------------------------------------------------\n");
  printf(" Operating parameters\n");
  printf("---------------------------------------------------\n");
  printf(" LB size = %i * %i * %i = %i\n", lb->esi, lb->esj, lb->esk, lb->nnodes);
  if (sizeof(LBReal) == sizeof(double)) printf(" Using double precision\n");
  else printf(" Using single precision\n");
  msize = (24*sizeof(LBReal) + sizeof(LBObj) + sizeof(int)) * (lb->nnodes);
  printf(" Memory requirement approx %0.2f MByte\n",(double)(msize)/(1024*1024));
  printf(" rho0 = %0.4f\n", lb->rho0);
  printf(" tau, nu, eta = %0.4f %0.4f %0.4f\n", lb->tau, nu, eta);
  printf(" Number of objects = %i\n", nobjs);
  printf(" Number of fluid + solid nodes = %i + %i = %i\n",
         nfluid,(lb->nnodes)-nfluid,(lb->nnodes));
  printf(" Expected total Fx, y, z = %0.2f %0.2f %0.2f\n", fxsum, fysum, fzsum);
  if (lb->stokes_flow) printf(" Stokes equations simulated\n");
  else printf(" Navier-Stokes equations simulated\n");
  if (lb->vel_corrected) printf(" Corrected velocity field used\n");
  else printf(" Uncorrected velocity field used\n");
  if (lb->skip_solid) printf(" Relaxation only on fluid nodes\n");
  else printf(" Relaxation on fluid and solid nodes\n");
  if (lb->objects_moving) printf(" Objects allowed to move\n");
  else printf(" All objects are stationary\n");
#ifdef SIGNAL
  printf(" Signalling implemented, process id = %i\n",(int)getpid());
#else
  printf(" Signalling not implemented\n");
#endif /* SIGNAL */
  if (lb->monitoring != 0) printf(" Monitoring data will be written to %s.mon\n", lb->name);
  else printf(" No monitoring\n");
  printf("---------------------------------------------------\n");
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : wrpars done\n", lb->step);
}


/* Short routine to write out position in (bbl, propagate, collide) cycle */

void wrstatus(LBSys lb, FILE *fp) {
  fprintf(fp,"Currently at ");
  switch (lb->status) {
    case POST_BBL: fprintf(fp,"post-bounceback"); break;
    case POST_PROPAGATE: fprintf(fp,"post-propagation"); break;
    case POST_COLLIDE: fprintf(fp,"post-collision"); break;
    default: fprintf(fp,"(help, I'm stymied!)"); break;
  }
  fprintf(fp," stage in the LB update cycle\n");
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : wrstatus done\n", lb->step);
}


/* Writes out current stats to a file handle.
 * This is also used to update the monitor file.
 */

void fwrstats(FILE *fp, int lbstep, int lbsize) {
  double rate = 0.0;
  time_t telapse;
  telapse = time(NULL) - tstart;
  if (telapse > 0) rate = (double)lbstep / (double)(telapse);
  fprintf(fp,"LB step %6i: ", lbstep);
  wrhms(fp,"long", telapse); 
  fprintf(fp," elapsed\n");
  fprintf(fp," %6.2f updates/sec\n", rate);
  fprintf(fp," %6.1f ksite-updates/sec\n", lbsize*rate/1000.0);
}

/* Use this to write to standard output */

void wrstats(LBSys lb) { fwrstats(stdout, lb->step, lb->nnodes); }


/* This routine writes out stats to the monitor file.
 * The second argument should be one of 'started' or 'initiated', 'running',
 * and 'finished' or 'ended'.  Only the first letter counts.
 */

void wrmon(LBSys lb, char *what) {
#ifdef SIGNAL
  int lbepid = 0;
#endif /* SIGNAL */
  double rate = 0.0;
  time_t telapse = 0;
  char *monfile;
  char w0;
  FILE *fp;
  if (lb->monitoring == 0) return;
  w0 = what[0];
  monfile = mycat(lb->name, ".mon");
  if ((fp = fopen(monfile, "w")) == NULL) {
    fprintf(stderr,"wrmon: I couldn't open %s.\n", monfile); 
    return;
  }
  switch (w0) {
  case 'i': case 's': 
    startclock();
    break;
  case 'r': case 'f': case 'e':
    telapse = time(NULL) - tstart;
    if (telapse > 0) rate = (double)(lb->step) / (double)(telapse);
    break;
  }
/* Write out a header line which can be picked up by sunlightlb_monitor */
#ifdef SIGNAL
  lbepid = (int)getpid();
  fprintf(fp,"%i %s ", lbepid, lb->name);
  switch (w0) {
  case 'i': case 's': 
    fprintf(fp,"started 0:00:00 ? ?\n");
    break;
  case 'r': case 'f': case 'e':
    if (w0 == 'r') fprintf(fp, "running ");
    else fprintf(fp, "finished ");
    wrhms(fp,"short", telapse);
    if (telapse > 0) {
      fprintf(fp," %0.2f %0.2f\n", rate, (lb->nnodes)*rate/1000.0);
    } else fprintf(fp," ? ?\n");
    break;
  }
#endif /* SIGNAL */
/* Write out rest of file in humanly digestible form */
  switch (w0) {
  case 'i': case 's': 
    break;
  case 'r': case 'f': case 'e':
#ifndef SIGNAL
    fprintf(fp, "%s ", lb->name);
    if (w0 == 'r') fprintf(fp, "running\n");
    else fprintf(fp, "finished\n");
#endif /* !SIGNAL */
    fwrstats(fp, lb->step, lb->nnodes);
    break;
  }
#ifdef SIGNAL
  if (w0 == 'i' || w0 == 's' || w0 == 'r') {
    fprintf(fp,"To update this file type : kill -USR1 %i\n", lbepid);
  }
#endif /* SIGNAL */
  if (fclose(fp) != 0) {
    fprintf(stderr,"wrmon: I couldn't close %s.\n", monfile); 
    exit(EXIT_FAILURE);
  }
  free(monfile);
  printf("Monitor file %s.mon %s at %i\n", lb->name, what, lb->step);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : wrmon done\n", lb->step);
}

/* Create space for a new string, copy 's' and 't' into it, return a
 * pointer to the new string.  This is equivalent to strdup
 * but provided for compatibility.  Returns pointer to new string.
 */

char *mycat(char *s, char *t) {
  char *v;
  v = (char *) malloc(strlen(s) + strlen(t) + 1);
  strcpy(v, s); strcat(v, t);
  return v;
}

/* Copy a string into a new space and return pointer */

char *mycopy(char *s) { return mycat(s,""); }

/* Compare two strings and return true
 * if 's' matches the first part of 't'.
 */

int myscmp(char *s, char *t) {
  return (strncmp(t, s, strlen(s)) == 0); 
}

/* Convert a time in seconds to hours, mins, secs and print to fp 
 * The conversion is set by fmt : "long" or "short"
 */

void wrhms(FILE *fp, char *fmt, int secs) {
  int hours, mins;
  if (secs < 0) return;
  hours = secs / 3600; secs -= 3600 * hours;
  mins = secs / 60; secs -= 60 * mins;
  if (fmt[0] == 'l') {
    if (hours == 0 && mins == 0 && secs == 0) {
      fprintf(fp,"no time at all");
    } else {
      switch (hours) {
      case 0: break;
      case 1: fprintf(fp,"1 hour "); break;
      default: fprintf(fp,"%i hours ", hours);
      }
      switch (mins) {
      case 0: break;
      case 1: fprintf(fp,"1 min "); break;
      default: fprintf(fp,"%i mins ", mins);
      }
      switch (secs) {
      case 0: fprintf(fp,"exactly"); break;
      case 1: fprintf(fp,"1 sec"); break;
      default: fprintf(fp,"%i secs", secs);
      }
    }
  } else {
    fprintf(fp,"%i:%02i:%02i", hours, mins, secs);
  }
}

/* Routines to allocate space */

LBReal *ralloc(int n) {
  LBReal *rarr;
  rarr = (LBReal *) malloc(n*sizeof(LBReal));
  if (rarr == NULL) error("No room for LBReal array");
  return rarr;
}

int *ialloc(int n) {
  int *iarr;
  iarr = (int *) malloc(n*sizeof(int));
  if (iarr == NULL) error("No room for int array");
  return iarr;
}

/* End of lbmisc.c */
