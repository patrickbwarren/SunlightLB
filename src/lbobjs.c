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

#define LBOBJS_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sunlightlb.h"

/* Local function prototypes */

static LBObj newobj();

/* Return a pointer to a new object structure */

LBObj newobj() {
  LBObj p;
  p = (LBObj) malloc(sizeof(LBObj_type));
  if (p == NULL) error("Ran out of memory for object creation\n");
  return p;
}


/* Add a new object to the top of the chain for object list in lb.
 * Increment the id count.
 */

LBObj addobj(LBSys lb) {
  LBObj p;
  p = newobj();
  p->base_obj = p;
  p->id = lb->next_id; lb->next_id++;
  p->next_obj = lb->obj_head; lb->obj_head = p;
  return p;
}


/* Zero object properties */

void zeroobj(LBObj p) {
  p->posx = p->posy = p->posz = 0.0;
  p->velx = p->vely = p->velz = 0.0;
  p->angx = p->angy = p->angz = 0.0;
  p->frcx = p->frcy = p->frcz = 0.0;
  p->trqx = p->trqy = p->trqz = 0.0;
}


/* Return a pointer to the object indexed by tid, or NULL if not found */

LBObj lookupobj(LBSys lb, int tid) {
  LBObj p;
  for (p=lb->obj_head; p!=NULL; p=p->next_obj) if (tid == p->id) return p;
  return NULL;
}


/* Insert substrate data specified in subsfile as a character array */

LBObj insertsubs(LBSys lb, char *subsfile, int ilo, int ihi, int jlo, int jhi, int klo, int khi) {
  int i, j, k, ns, inode, tot = 0, overlap = 0;
  int esj, esk;
  unsigned char byte;
  FILE *fp;
  LBObj p,*obj;
  if ((fp = fopen(subsfile,"r")) == NULL) {
    printf("insertsubs: %s could not be opened\n", subsfile); return NULL;
  }
  if (slb_verbosity & OBJECT_DETAILS) {
    printf("\nInserting substrate data from %s\n", subsfile);
    printf("Substrate region is (%i-%i, %i-%i, %i-%i) inclusive\n", 
	   ilo, ihi, jlo, jhi, klo, khi);
  }
/* Read in data and insert in map, */
/* note array size must match that of data file */
/* also note k varies fastest */
  obj = lb->obj; esj = lb->esj; esk = lb->esk; 
  p = addobj(lb); zeroobj(p); ns = 0;
  for (i=ilo; i<=ihi; i++)
    for (j=jlo; j<=jhi; j++)
      for (k=klo; k<=khi; k++) {
	inode = esk*(esj*i+j)+k;
	if (fscanf(fp,"%c",&byte) != 1) {
	  if (slb_verbosity & OBJECT_DETAILS) printf("Ran out of data in %s\n", subsfile);
	  break;
	}
	if (byte & 1) {
	  if (obj[inode] != NULL) overlap = 1;
	  obj[inode] = p; ns++;
	}
	if (byte & 2) { 
	  tot++;
	}
      }
  if (slb_verbosity & OBJECT_DETAILS) {
    printf("Inserted substrate data, id no = %i, # nodes = %i\n", p->id, ns);
    if (overlap) printf("Overlaps occurred\n");
  }
  if (ns == 0) {
    if (slb_verbosity & OBJECT_DETAILS) printf("No sites set in object\n");
    lb->obj_head = p->next_obj; free(p); return NULL;
  } else lb->changed_objects = 1;
  fclose(fp);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : insertsubs done\n", lb->step);
  return p;
}


/* Insert substrate data specified in subsfile as a list of nodes */

LBObj insertlist(LBSys lb, char *listfile) {
  int i, j, k, ns, inode, overlap = 0;
  int esi, esj, esk;
  FILE *fp;
  LBObj p,*obj;
  if ((fp = fopen(listfile,"r")) == NULL) {
    printf("insertlist: %s could not be opened\n", listfile); return NULL;
  }
  if (slb_verbosity & OBJECT_DETAILS) {
    printf("\nInserting substrate node list data from %s\n", listfile);
  }
/* Read in data and insert in map, */
  obj = lb->obj; esi = lb->esi; esj = lb->esj; esk = lb->esk; 
  p = addobj(lb); zeroobj(p); ns = 0;
  for (;;) {
    if (fscanf(fp,"%i%i%i",&i,&j,&k) != 3) break;
    if (i >= 0 && i < esi && j >= 0 && j < esj && k >= 0 && k < esk) {
      inode = esk*(esj*i+j)+k;
      if (obj[inode] != NULL) overlap = 1;
      obj[inode] = p; ns++; 
    } else {
      if (slb_verbosity & OBJECT_DETAILS) {
        printf("insertlist: ignoring node (%i,%i,%i), it is out of range\n", i, j, k);
      }
    }
  }
  if (slb_verbosity & OBJECT_DETAILS) {
    printf("Inserted node list data, id no = %i, # nodes = %i\n", p->id, ns);
    if (overlap) printf("Overlaps occurred\n");
  }
  if (ns == 0) {
    if (slb_verbosity & OBJECT_DETAILS) printf("No sites set in object\n");
    lb->obj_head = p->next_obj; free(p); return NULL;
  } else lb->changed_objects = 1;
  fclose(fp);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : insertlist done\n", lb->step);
  return p;
}


/* Insert a slab from x1, y1, z1 to x2, y2, z2 */

LBObj insertslab(LBSys lb, LBReal x1, LBReal x2, LBReal y1, LBReal y2, LBReal z1, LBReal z2) {
  int i, j, k, ns, inode;
  int esi, esj, esk;
  int overlap = 0;
  LBObj p,*obj;
  obj = lb->obj; esi = lb->esi; esj = lb->esj; esk = lb->esk; 
  p = addobj(lb); zeroobj(p); ns = 0;
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    if (i > x1 && i < x2 && j > y1 && j < y2 && k > z1 && k < z2) { 
      if (obj[inode] != NULL) overlap = 1;
      obj[inode] = p; ns++;
    }
  }
  if (slb_verbosity & OBJECT_DETAILS) {
    printf("Inserted slab ");
    printf("%5.2f < x < %5.2f, ", x1, x2);
    printf("%5.2f < y < %5.2f, ", y1, y2);
    printf("%5.2f < z < %5.2f, ", z1, z2);
    printf("   id no = %i, # nodes = %i\n", p->id, ns);
    if (overlap) printf("Overlaps occurred\n");   
  }
  if (ns == 0) {
    if (slb_verbosity & OBJECT_DETAILS) printf("No sites set in object\n");
    lb->obj_head = p->next_obj; free(p); return NULL;
  } else lb->changed_objects = 1;
  p->posx = 0.5*(x1+x2-1); p->posy = 0.5*(y1+y2-1); p->posz = 0.5*(z1+z2-1);
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : insertslab done\n", lb->step);
  return p;
}


/* Insert a sphere radius rad, centre ox, y, z */

LBObj insertsph(LBSys lb, LBReal rad, LBReal ox, LBReal oy, LBReal oz) {
  int i, j, k, ns, inode;
  int esi, esj, esk;
  int overlap = 0;
  LBReal rsq, radsq;
  LBObj p,*obj;
  obj = lb->obj; esi = lb->esi; esj = lb->esj; esk = lb->esk; 
  p = addobj(lb); zeroobj(p); ns = 0; radsq = rad*rad;
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    rsq = (i-ox)*(i-ox) + (j-oy)*(j-oy) + (k-oz)*(k-oz);
    if (rsq < radsq) { 
      if (obj[inode] != NULL) overlap = 1;
      obj[inode] = p; ns++;
    }
  }
  if (slb_verbosity & OBJECT_DETAILS) {
    printf("Sphere radius %5.2f at %6.2f %6.2f %6.2f ", rad, ox, oy, oz);
    printf("id = %5i, # nodes = %i", p->id, ns);
    if (overlap) printf(" (overlaps)\n"); else printf("\n");
  }
  if (ns == 0) {
    if (slb_verbosity & OBJECT_DETAILS) printf("No sites set in object\n");
    lb->obj_head = p->next_obj; free(p); lb->next_id--; return NULL; 
  } else lb->changed_objects = 1;
  p->posx = ox; p->posy = oy; p->posz = oz;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : insertsph done\n", lb->step);
  return p;
}


/* Insert a cylinder radius rad,
 * axis (a, b) = (ox, oy) or (oy, oz) or (ox, oz),
 * direction determined by first letter of type,
 * whether solid or a cylindrical hole by the second letter.
 */

LBObj insertcyl(LBSys lb, char *type, LBReal rad, LBReal a, LBReal b) {
  int i, j, k, ns, inode;
  int esi, esj, esk;
  int within, overlap = 0;
  LBReal rsq, radsq;
  LBObj p,*obj;
  obj = lb->obj; esi = lb->esi; esj = lb->esj; esk = lb->esk; 
  p = addobj(lb); zeroobj(p); ns = 0; radsq = rad*rad;
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    switch (type[0]) {
      case 'x': case 'i': rsq = (j-a)*(j-a) + (k-b)*(k-b); break;
      case 'y': case 'j': rsq = (i-a)*(i-a) + (k-b)*(k-b); break;
      case 'z': case 'k': rsq = (i-a)*(i-a) + (j-b)*(j-b); break;
      default: error("Unrecognised cylinder direction"); return NULL;
    }
    within = (rsq < radsq); if (type[1] != '\0') within = !within;
    if (within) { 
      if (obj[inode] != NULL) overlap = 1;
      obj[inode] = p; ns++;
    }
  }
  if (slb_verbosity & OBJECT_DETAILS) {
    if (type[1] == '\0') printf("Cylinder ");
    else printf("Cylindrical hole ");
    printf("radius %5.2f, axis %5.2f %5.2f along %c ", rad, a, b, type[0]);
    printf("id = %i, # nodes = %i", p->id, ns);
    if (ns == 0) printf(" (no sites set)");
    if (overlap) printf(" (overlaps)\n"); else printf("\n");
  }
  if (ns == 0) {
    lb->obj_head = p->next_obj; free(p); lb->next_id--; return NULL; 
  } else lb->changed_objects = 1;
  switch (type[0]) {
    case 'x': case 'i':
      p->posx = 0.5*(esi-1); p->posy = a; p->posz = b;
      break;
    case 'y': case 'j':
      p->posx = a; p->posy = 0.5*(esj-1); p->posz = b;
      break;
    case 'z': case 'k':
      p->posx = a; p->posy = b; p->posz = 0.5*(esk-1);
      break;
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : insertcyl done\n", lb->step);
  return p;
}


/* Insert a wall, position wallpos, 
 * direction determined by first letter of type,
 */

LBObj insertwall(LBSys lb, char *type, int wallpos) {
  int i, j, k, ns, inode;
  int esi, esj, esk;
  int overlap = 0;
  LBObj p,*obj;
  obj = lb->obj; esi = lb->esi; esj = lb->esj; esk = lb->esk; 
  p = addobj(lb); zeroobj(p); ns = 0;
  switch (type[0]) {
    case 'x': case 'i':
      if ((i = wallpos) < 0 || i > esi-1) {
        printf("insertwall: wallpos out of range\n"); return NULL;
      }
      for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
	inode = esk*(esj*i+j)+k;
        if (obj[inode] != NULL) overlap = 1;
        obj[inode] = p; ns++;
      }
      break;
    case 'y': case 'j':
      if ((j = wallpos) < 0 || j > esj-1) {
        printf("insertwall: wallpos out of range\n"); return NULL;
      }
      for (i=0; i<esi; i++) for (k=0; k<esk; k++) {
	inode = esk*(esj*i+j)+k;
        if (obj[inode] != NULL) overlap = 1;
        obj[inode] = p; ns++;
      }
      break;
    case 'z': case 'k':
      if ((k = wallpos) < 0 || k > esk-1) {
        printf("insertwall: wallpos out of range\n"); return NULL;
      }
      for (i=0; i<esi; i++) for (j=0; j<esj; j++) {
	inode = esk*(esj*i+j)+k;
        if (obj[inode] != NULL) overlap = 1;
        obj[inode] = p; ns++;
      }
      break;
    default: 
      printf("insertwall: unrecognised direction for walls\n"); return NULL;
  }
  switch (type[0]) {
    case 'x': case 'i':
      p->posx = wallpos; p->posy = 0.5*(esj-1); p->posz = 0.5*(esk-1);
      break;
    case 'y': case 'j':
      p->posx = 0.5*(esi-1); p->posy = wallpos; p->posz = 0.5*(esk-1);
      break;
    case 'z': case 'k':
      p->posx = 0.5*(esi-1); p->posy = 0.5*(esj-1); p->posz = wallpos;
      break;
  }
  if (ns == 0) {
    if (slb_verbosity & OBJECT_DETAILS) printf("No sites set in object\n");
    lb->obj_head = p->next_obj; free(p); return NULL;
  } else lb->changed_objects = 1;
  if (slb_verbosity & OBJECT_DETAILS) {
    printf("Inserted wall at %c = %i, ", type[0], wallpos);
    printf("id no = %i, # nodes = %i\n", p->id, ns);
    if (overlap) printf("Overlaps occurred\n");
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : insertwall done\n", lb->step);
  return p;
}


/* Flatten all objects to a single object */

LBObj flattenobj(LBSys lb) {
  int i, j, k, ns, inode;
  int esi, esj, esk;
  LBReal ox, oy, oz;
  LBObj p, pp,*obj;
  obj = lb->obj; esi = lb->esi; esj = lb->esj; esk = lb->esk; 
/* Drop everything in the list currently, then add a single new object */
  p = lb->obj_head;
  while (p != NULL) {
    pp = p; p = p->next_obj; free(pp); 
  }
  lb->next_id = 0; lb->obj_head = NULL;
  p = addobj(lb); zeroobj(p); ns = 0; ox = oy = oz = 0.0;
/* Now go through object map and make all object sites point to this structure */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    if (obj[inode] != NULL) {
      ox += (double)i; oy += (double)j; oz += (double)k; ns++; obj[inode] = p;
    }
  }
  if (ns == 0) {
    if (slb_verbosity & OBJECT_DETAILS) printf("No sites set in object\n");
    lb->obj_head = p->next_obj; free(p); return NULL;
  } else lb->changed_objects = 1;
  if (slb_verbosity & OBJECT_DETAILS) {
    printf("Object space flattened to id no %i, # nodes = %i\n", p->id, ns);
  }
  p->posx = ox/ns; p->posy = oy/ns; p->posz = oz/ns;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : flattenobj done\n", lb->step);
  return p;
}


/* Invert solid object space, and generate a single object, similar to above */

LBObj invertobj(LBSys lb) {
  int i, j, k, ns, inode;
  int esi, esj, esk;
  LBReal ox, oy, oz;
  LBObj p, pp, *obj;
  obj = lb->obj; esi = lb->esi; esj = lb->esj; esk = lb->esk; 
/* Drop everything in the list currently, then add a single new object */
  p = lb->obj_head;
  while (p != NULL) {
    pp = p; p = p->next_obj; free(pp); 
  }
  lb->next_id = 0; lb->obj_head = NULL;
  p = addobj(lb); zeroobj(p); ns = 0; ox = oy = oz = 0.0;
/* Now go through object map and invert so fluid sites point to this structure */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    if (obj[inode] != NULL) obj[inode] = NULL;
    else {
      ox += (double)i; oy += (double)j; oz += (double)k; ns++; obj[inode] = p;
    }
  }
  if (ns == 0) {
    if (slb_verbosity & OBJECT_DETAILS) printf("No sites set in object\n");
    lb->obj_head = p->next_obj; free(p); return NULL;
  } else lb->changed_objects = 1;
  if (slb_verbosity & OBJECT_DETAILS) {
    printf("Object space inverted to id no %i, # nodes = %i\n", p->id, ns);
  }
  p->posx = ox/ns; p->posy = oy/ns; p->posz = oz/ns;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : invertobj done\n", lb->step);
  return p;
}


/* Remove all objects and reset LB distribution function */

void removeobjs(LBSys lb) {
  int i, inode;
  LBObj p, pp, *obj;
  obj = lb->obj;
/* Drop everything in the list */
  p = lb->obj_head;
  while (p != NULL) { pp = p; p = p->next_obj; free(pp); }
  lb->next_id = 0; lb->obj_head = NULL;
/* Now go through object map and remove all objects */
  for (inode=0; inode<(lb->nnodes); inode++) {
    obj[inode] = NULL;
    lb->ux[inode] = lb->uy[inode] = lb->uz[inode] = 0.0;
    lb->fx[inode] = lb->fy[inode] = lb->fz[inode] = 0.0;
  }
  for (i=0; i<15*(lb->nnodes); i++) lb->a[i] = 0.0;
  if (slb_verbosity & OBJECT_DETAILS) {
    printf("All objects removed\n");
  }
  lb->changed_objects = 1;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : removeobjs done\n", lb->step);
}

/* Write out object properties */

void wrobjs(LBSys lb) {
  LBObj p;
  for (p=lb->obj_head; p!=NULL; p=p->next_obj) {
    printf("Object id (base) = %i (%i)\n", p->id, p->base_obj->id);
    printf(" posx, y, z = %f %f %f\n", p->posx, p->posy, p->posz);
    printf(" velx, y, z = %f %f %f\n", p->velx, p->vely, p->velz);
    printf(" angx, y, z = %f %f %f\n", p->angx, p->angy, p->angz);
    printf(" frcx, y, z = %f %f %f\n", p->frcx, p->frcy, p->frcz);
    printf(" trqx, y, z = %f %f %f\n", p->trqx, p->trqy, p->trqz);
  }
}


/* Write out a map of the object sites in humanly readable format */

void wrmap(LBSys lb, char *name) {
  int i, j, k, id, esi, esj, esk, inode, len;
  int *count;
  FILE *fp;
  LBObj *obj, p;
  char c;
  char shape[] = "ox+*#acemnrsuvwxz?";
  fp = fopen(name,"w");
  if (fp == NULL) {
    fprintf(stderr,"Couldn't open %s\n", name);
    return;
  }
  len = strlen(shape);
  count = (int *) malloc(len*sizeof(int));
  for (i=0; i<len; i++) count[i] = 0;
  obj = lb->obj; esi = lb->esi; esj = lb->esj; esk = lb->esk;
  for (k=0; k<esk; k++) {
    fprintf(fp,"\nSlice k = %3i\n", k);
    fprintf(fp,"            ");
    for (i=0; i<esi; i++) fprintf(fp,"%c", i<10?' ':('0'+(i/10)%10)); fprintf(fp,"\n");
    fprintf(fp,"Col i =     ");
    for (i=0; i<esi; i++) fprintf(fp,"%c",('0'+(i%10))); fprintf(fp,"\n");
    for (j=0; j<esj; j++) {
      fprintf(fp,"Row j = %3i ", j);
      for (i=0; i<esi; i++) {
	inode = esk*(esj*i+j)+k; p = obj[inode];
        if (p ==  NULL) c = '.';
	else {
	  id = p->id; if (id >= len) id = len - 1;
	  c = shape[id]; count[id]++;
	}
        fprintf(fp,"%c", c);
      }
      fprintf(fp,"\n");
    }
  }
  fclose(fp);
  printf("\n id shape count\n");
  for (i=0; i<len; i++) {
    if (count[i] > 0) printf("%3i   %c %6i\n", i, shape[i], count[i]);
  }
  printf("\n");
  free(count);
  printf("Map of object space written to %s\n", name);
}


/* Build up a mask for the location of the boundary links
 * using obj: mask is non zero for those nodes with
 * attached boundary links, the actual attached links are
 * coded in binary.  One bit is set aside to mark adsorption
 * sites.
 */

void makemask(LBSys lb) {
  int i, j, k, ip, jp, kp, im, jm, km, val, inode;
  int esi, esj, esk;
  int *mask;
  LBObj *obj;
/* Run through all nodes and neighbours (avoiding double counting) */
  obj = lb->obj; mask = lb->mask; esi = lb->esi; esj = lb->esj; esk = lb->esk;
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k; val = 0;
    if (obj[inode] != NULL) {
      ip = i + 1; if (ip >= esi) ip = 0;
      jp = j + 1; if (jp >= esj) jp = 0;
      kp = k + 1; if (kp >= esk) kp = 0;
      im = i - 1; if (im < 0) im = esi-1;
      jm = j - 1; if (jm < 0) jm = esj-1;
      km = k - 1; if (km < 0) km = esk-1;
      if (obj[esk*(esj*ip+j )+k ] == NULL) val |= 0x0001;
      if (obj[esk*(esj*im+j )+k ] == NULL) val |= 0x0002;
      if (obj[esk*(esj*i +jp)+k ] == NULL) val |= 0x0004;
      if (obj[esk*(esj*i +jm)+k ] == NULL) val |= 0x0008;
      if (obj[esk*(esj*i +j )+kp] == NULL) val |= 0x0010;
      if (obj[esk*(esj*i +j )+km] == NULL) val |= 0x0020;
      if (obj[esk*(esj*ip+jp)+kp] == NULL) val |= 0x0040;
      if (obj[esk*(esj*im+jp)+kp] == NULL) val |= 0x0080;
      if (obj[esk*(esj*ip+jm)+kp] == NULL) val |= 0x0100;
      if (obj[esk*(esj*im+jm)+kp] == NULL) val |= 0x0200;
      if (obj[esk*(esj*ip+jp)+km] == NULL) val |= 0x0400;
      if (obj[esk*(esj*im+jp)+km] == NULL) val |= 0x0800;
      if (obj[esk*(esj*ip+jm)+km] == NULL) val |= 0x1000;
      if (obj[esk*(esj*im+jm)+km] == NULL) val |= 0x2000;
    }
    mask[inode] = val;
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : makemask done\n", lb->step);
}


/* Compute the object node velocity field */

void vspread(LBSys lb) {
  int i, j, k, inode;
  int esi, esj, esk;
  LBReal *ux, *uy, *uz;
  LBReal rx, ry, rz, rho0;
  LBObj *obj, p;
  esi = lb->esi; esj = lb->esj; esk = lb->esk;
  ux = lb->ux; uy = lb->uy; uz = lb->uz;
  obj = lb->obj; rho0 = lb->rho0;
/* Run through all boundary links */
  for (i=0; i<esi; i++) for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
    inode = esk*(esj*i+j)+k;
    if ((p = obj[inode]) == NULL) {
      ux[inode] = uy[inode] = uz[inode] = 0.0;
    } else {
/* Calculate velocity components, and factor in 2D rho0 / b c^2 */
/* When added to the a_i, also need to include the weight w_i */
      rx = i - (p->posx); ry = j - (p->posy); rz = k - (p->posz);
      ux[inode] = rho0 * ((p->velx) + (p->angy)*rz - (p->angz)*ry) / 12.0;
      uy[inode] = rho0 * ((p->vely) + (p->angz)*rx - (p->angx)*rz) / 12.0;
      uz[inode] = rho0 * ((p->velz) + (p->angx)*ry - (p->angy)*rx) / 12.0;
    }
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : vspread done\n", lb->step);
}

/* Zero body force */

void zerobf(LBSys lb) {
  int inode;
  LBReal *bx, *by, *bz;
  LBObj *obj;
  bx = lb->bx; by = lb->by; bz = lb->bz; obj = lb->obj;
  for (inode=0; inode<(lb->nnodes); inode++) {
    if (obj[inode] == NULL) {
      bx[inode] = by[inode] = bz[inode] = 0.0;
    }
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : zerobf done\n", lb->step);
}

/* Add a uniform body force only to fluid nodes */

void addbf(LBSys lb, double dpdx, double dpdy, double dpdz) {
  int inode;
  LBReal *bx, *by, *bz;
  LBObj *obj;
  bx = lb->bx; by = lb->by; bz = lb->bz; obj = lb->obj;
  for (inode=0; inode<(lb->nnodes); inode++) {
    if (obj[inode] == NULL) {
      bx[inode] += (LBReal)dpdx; by[inode] += (LBReal)dpdy; bz[inode] += (LBReal)dpdz;
    }
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : zerobf done\n", lb->step);
}

/* Set an entry in the object node velocity field */

void vset(LBSys lb, int i, int j, int k, LBReal vx, LBReal vy, LBReal vz){
  int inode;
  int esi, esj, esk;
  LBReal *ux, *uy, *uz;
  LBReal rho0;
  LBObj *obj;
  esi = lb->esi; esj = lb->esj; esk = lb->esk;
  ux = lb->ux; uy = lb->uy; uz = lb->uz;
  obj = lb->obj; rho0 = lb->rho0;
  if (i >= 0 && i < esi && j >= 0 && j < esj && k >= 0 && k < esk) {
    inode = esk*(esj*i+j)+k;
    if (obj[inode] != NULL) {
/* Calculate velocity components, and factor in 2D rho0 / b c^2 */
/* When added to the a_i, also need to include the weight w_i */
      ux[inode] = rho0 * vx / 12.0;
      uy[inode] = rho0 * vy / 12.0;
      uz[inode] = rho0 * vz / 12.0;
    }
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : vset done\n", lb->step);
}


/* Count number of fluid nodes and sync with lb sys info. */

int countfluid(LBSys lb) {
  int inode, count;
  LBObj *obj;
  obj = lb->obj; count = 0; 
  for (inode=0; inode<(lb->nnodes); inode++) if (obj[inode] == NULL) count++;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : countfluid done, = %i\n", lb->step, count);
  return count;
}


/* Sync the objects to the LB system. */

void syncobjs(LBSys lb) {
  makemask(lb); countfluid(lb); lb->changed_objects = 0;
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : syncobjs done.\n", lb->step);
}

/* Utility routine for surface area measurements.
 * Either return the total number of boundary links (type = "link"
 * or "all"), or the number of faces or facets (type = "faces" or
 * "facets" or "squares").  (Only the first letter of type counts.)
 */

int countmask(LBSys lb, char *type) {
  int inode, m, count = 0;
  int *mask;
  mask = lb->mask;
  switch (type[0]) {
  case 'a': case 'l' :
    for (inode=0; inode<(lb->nnodes); inode++) {
      m = mask[inode];
      if (m > 0) {
	if (m & 0x0001) count += lb->w[ 1];
	if (m & 0x0002) count += lb->w[ 2];
	if (m & 0x0004) count += lb->w[ 3];
	if (m & 0x0008) count += lb->w[ 4];
	if (m & 0x0010) count += lb->w[ 5];
	if (m & 0x0020) count += lb->w[ 6];
	if (m & 0x0040) count += lb->w[ 7];
	if (m & 0x0080) count += lb->w[ 8];
	if (m & 0x0100) count += lb->w[ 9];
	if (m & 0x0200) count += lb->w[10];
	if (m & 0x0400) count += lb->w[11];
	if (m & 0x0800) count += lb->w[12];
	if (m & 0x1000) count += lb->w[13];
	if (m & 0x2000) count += lb->w[14];
      }
    }
    break;
  case 'f': case 's' :
    for (inode=0; inode<(lb->nnodes); inode++) {
      m = mask[inode];
      if (m > 0) {
	if (m & 0x0001) count++;
	if (m & 0x0002) count++;
	if (m & 0x0004) count++;
	if (m & 0x0008) count++;
	if (m & 0x0010) count++;
	if (m & 0x0020) count++;
      }
    }
    break;
  default: printf("countmask: expected 'links' or 'faces'\n");
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("LB %i : countmask done\n", lb->step);
  return count;
}

/* End of lbobjs.c */
