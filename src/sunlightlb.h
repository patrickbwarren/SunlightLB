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

#ifndef SUNLIGHTLB_H_INC
#define SUNLIGHTLB_H_INC

/* Use these to switch between single and double precision */

#ifdef _SINGLE
typedef float LBReal;
#else
typedef double LBReal;
#endif

#ifdef _STORE_DOUBLE
typedef double LBStoreReal;
#else
typedef float LBStoreReal;
#endif

/* Verbosity level setting */

#define TRACE_ROUTINE_EXITS   0x02
#define TRACING_ROUTINE_EXITS 0x02
#define OBJECT_DETAILS        0x04
#define MP_OBJECT_DETAILS     0x04
#define FORCE_DETAILS         0x08

/* For specifying the status...arbitrary numbering */

#define POST_BBL       42
#define POST_PROPAGATE 47
#define POST_COLLIDE   49

/* A structure for lattice objects */

typedef struct LBObj_struct {
  int id;
  LBReal posx, posy, posz;
  LBReal velx, vely, velz;
  LBReal angx, angy, angz;
  LBReal frcx, frcy, frcz;
  LBReal trqx, trqy, trqz;
  struct LBObj_struct *next_obj;
  struct LBObj_struct *base_obj;
} LBObj_type;

typedef LBObj_type *LBObj;

/* A structure for lattice Boltzmann things */

typedef struct LBSys_struct {
/* Weights and velocities associated with D3Q15(r=16) model */
  int *w, *c;
/* Storage for occupation probabilities D3Q15 (r=16) model */
/* These are defined by a_i = w_i (f_i - rho_0 / b ) */
  LBReal *a;
/* Storage for object node velocity data, and measured fluid velocity */
  LBReal *ux, *uy, *uz;
/* Storage for object node force data */
  LBReal *fx, *fy, *fz;
/* Storage for local body forces */
  LBReal *bx, *by, *bz;
/* General operating parameters */
  int esi, esj, esk, nnodes, step, status, changed_objects;
  int objects_moving, skip_solid, stokes_flow, vel_corrected, monitoring;
  LBReal tau, rho0;
/* Storage for measured mean velocities */
  LBReal vmeanx, vmeany, vmeanz;
/* Storage for a weight function psi(r) for averaging, NULL if absent */
  LBReal *psi;
/* Storage for the density, may be NULL */
  LBReal *rho;
/* Used for name of monitoring file, etc */
  char *name;
/* These are kludges in storage routines */
  LBReal rescale;
  int load_average;
/* This is for objects embedded in the lattice */
  LBObj *obj;
  LBObj obj_head;
  int *mask;
  int next_id;
} LBSys_type;

typedef LBSys_type *LBSys;

/* A structure for moment propagation things */

typedef struct MPSys_struct {
/* Pointer to the associated LB system */
  LBSys lb;
/* Storage for the tracer particle density, etc */
  LBReal *p, *newp;
/* General operating parameters */
  int step;
  int objects_moving;
  int to_be_initialised;
  LBReal alpha, isoslice, poffset;
  LBReal p2bytemin, p2bytemax;
/* Name used for monitoring */
  char *name;
} MPSys_type;

typedef MPSys_type *MPSys;

/* A structure for isosurfaces */

typedef struct MPIso_struct {
  MPSys mp;
  LBReal val;
  int *count;
} MPIso_type;

typedef MPIso_type *MPIso;

/* Structures for holding colours, colour maps, and images */

typedef struct LBColour_struct {
  int red, green, blue;
} LBColour_type;

typedef LBColour_type *LBColour;

typedef struct LBColourMap_struct {
/* Colour range */
  int red_lo, green_lo, blue_lo;
  int red_hi, green_hi, blue_hi;
/* Range of values that this maps into */
  LBReal v_lo, v_hi;
} LBColourMap_type;

typedef LBColourMap_type *LBColourMap;

typedef struct LBImage_struct {
/* Image size */
  int width, height;
/* Image colours */
  int *red, *green, *blue;
/* Background and object cols */
  int back_red, back_green, back_blue;
  int object_red, object_green, object_blue;
} LBImage_type;

typedef LBImage_type *LBImage;

/* There is one external variable that is defined in lbmisc.c */

#ifdef LBMISC_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int slb_verbosity;


/* Function prototypes associated with basic LB */

LBSys newlb(int, int, int, char *);
void set_verbosity(int);
void unset_verbosity(int);
void report_verbosity();
void lbstep(LBSys);
void bbl(LBSys);
void bblforce(LBSys);
void collide(LBSys);
void propagate(LBSys);
LBObj lookupobj(LBSys, int);
LBObj addobj(LBSys);
void zeroobj(LBObj);
LBObj insertsubs(LBSys, char *, int, int, int, int, int, int);
LBObj insertlist(LBSys, char *);
LBObj insertslab(LBSys, LBReal, LBReal, LBReal, LBReal, LBReal, LBReal);
LBObj insertsph(LBSys, LBReal, LBReal, LBReal, LBReal);
LBObj insertcyl(LBSys, char *, LBReal, LBReal, LBReal);
LBObj insertwall(LBSys, char *, int);
LBObj flattenobj(LBSys);
LBObj invertobj(LBSys);
void removeobjs(LBSys);
void wrobjs(LBSys);
void wrmap(LBSys, char *);
void makemask(LBSys);
void vspread(LBSys);
void zerobf(LBSys);
void addbf(LBSys, double, double, double);
void vset(LBSys, int, int, int, LBReal, LBReal, LBReal);
int countmask(LBSys, char *);
int countfluid(LBSys);
void syncobjs(LBSys);
void lbsave(LBSys, char *);
void lbload(LBSys, char *);
void lbsavevel(LBSys, char *, char *);
void lbloadvel(LBSys, char *, char *);
void lbsavestr(LBSys, char *, char *);
void lbsaveobj(LBSys, char *);
void lbloadobj(LBSys, char *);
void lbsavefrc(LBSys, char *, char *);
void lbsavemap(LBSys, char *, char *);
void error(char *);
void about();
void wrpars(LBSys);
void wrstatus(LBSys, FILE *);
void wrlinks(LBSys);
void wrstats(LBSys);
void fwrstats(FILE *, int, int);
void startclock();
char *mycat(char *, char *);
char *mycopy(char *);
int myscmp(char *, char *);
void accvmean(LBSys, LBReal);
void vzero(LBSys);
void vgrab(LBSys, LBReal);
void vsave(LBSys, char *);
void measure(LBSys);
void zeroforce(LBSys);
void normforce(LBSys, LBReal);
void distribf(LBSys);
void wrmeas(LBSys);
LBReal discrf(LBSys);
void wrmon(LBSys, char *);

/* Function prototypes associated with moment propagation */

void mpstep(MPSys);
MPSys newmp(LBSys, char *);
MPIso newmpiso(MPSys, LBReal);
void mpiso(MPIso);
void mpisowr(MPIso, char *, char *);
void mpinit(MPSys);
void mpwrpars(MPSys);
void mprhoinit(MPSys);
void mpcheck(MPSys);
LBReal mpmean(MPSys);
LBReal mpmax(MPSys);
LBReal mpmin(MPSys);
void mprescale(MPSys, LBReal);
void mpcopy(MPSys, MPSys);
void mpsave(MPSys, char *, char *);
void mpload(MPSys, char *, char *);
void mpconwr(MPSys, char *, char *);
LBReal mpspread(MPSys, char *);
void mpall(MPSys, LBReal);
void mpslab(MPSys, LBReal, LBReal, LBReal, LBReal, LBReal, LBReal, LBReal);
void mpwall(MPSys, char *, int, LBReal);
void mpsph(MPSys, LBReal, LBReal, LBReal, LBReal, LBReal);
void mpcyl(MPSys, char *, LBReal, LBReal, LBReal, LBReal);
void mpsubs(MPSys, char *, int, int, int, int, int, int, LBReal);
void mplist(MPSys, char *, LBReal);
void mpvacfinit(MPSys, MPSys, char *, int);
void mpvacfinit2(MPSys, MPSys);
void mpvacf(MPSys, MPSys);
LBReal mpvacfwr(MPSys, MPSys, char *, int);
void adall(MPSys);
void adnone(MPSys);
void adslab(MPSys, LBReal, LBReal, LBReal, LBReal, LBReal, LBReal);
void adsph(MPSys, LBReal, LBReal, LBReal, LBReal);
void adcyl(MPSys, char *, LBReal, LBReal, LBReal);
LBReal adtotal(MPSys, char *, int);
void adzero(MPSys);

/* Function prototypes associated with making images for LB systems */

LBColour newlbcolour(int, int, int);
void lbcolourset(LBColour, int, int, int);
void lbcolourcopy(LBColour, LBColour);
LBColourMap newlbcolourmap(LBReal, LBReal, LBColour, LBColour);
void lbcolourmapset(LBColourMap, LBReal, LBReal, LBColour, LBColour);
LBImage newlbimage(int, int, LBColour, LBColour);
void lbimageclear(LBImage);
void lbimagemap(LBImage, LBSys, char *, int);
void lbimagevel(LBImage, LBColourMap, LBSys, char *, char *, int);
void lbcolourget(LBImage, LBColour, int, int);
int lbgetred(LBImage, int, int);
int lbgetgreen(LBImage, int, int);
int lbgetblue(LBImage, int, int);
void lbcolourclamp(int *, int *, int *);
void mpimageset(LBImage, MPSys, LBColourMap, char *, char *, int);

/* Function prototypes dealing with signal handling, compiled in only if required */

#ifndef _NO_SIGNALLING

#define SIGNAL

#ifdef LBSIGNAL_C
#define EXTERN_LBSIGNAL
#else
#define EXTERN_LBSIGNAL extern
#endif

EXTERN_LBSIGNAL int slb_signal_flag;

void sighandle(LBSys);
void sighook(int);
void siginit();

#endif /* _NO_SIGNALLING */

#endif /* SUNLIGHTLB_H_INC */
