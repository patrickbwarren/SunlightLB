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

#define MPIMAGE_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "sunlightlb.h"


/* Make a new colour, with various helper functions */

LBColour newlbcolour(int red, int green, int blue) {
  LBColour clr;
  clr = (LBColour) malloc(sizeof(LBColour_type));
  lbcolourset(clr, red, green, blue);
  return clr;
}

void lbcolourset(LBColour clr, int red, int green, int blue) {
  clr->red = red; clr->green = green; clr->blue = blue;
  lbcolourclamp(&clr->red, &clr->green, &clr->blue);
}

void lbcolourcopy(LBColour dest, LBColour src) {
  dest->red = src->red; dest->green = src->green; dest->blue = src->blue;
}

void lbcolourclamp(int *red, int *blue, int *green) {
  if (*red < 0) *red = 0;
  else if (*red > 255) *red = 255;
  if (*green < 0) *green = 0;
  else if (*green > 255) *green = 255;
  if (*blue < 0) *blue = 0;
  else if (*blue > 255) *blue = 255;
}

/*
 * Make a new colour map in range (vlo, vhi)
 * with colours at low and high ends as specified.
 */

LBColourMap newlbcolourmap(LBReal v_lo, LBReal v_hi, LBColour c_lo, LBColour c_hi) {
  LBColourMap cm;
  cm = (LBColourMap) malloc(sizeof(LBColourMap_type));
  lbcolourmapset(cm, v_lo, v_hi, c_lo, c_hi);
  return cm;
}

void lbcolourmapset(LBColourMap cm, LBReal v_lo, LBReal v_hi, LBColour c_lo, LBColour c_hi) {
  cm->v_lo = v_lo; cm->v_hi = v_hi;
  cm->red_lo = c_lo->red; cm->green_lo = c_lo->green; cm->blue_lo = c_lo->blue;
  cm->red_hi = c_hi->red; cm->green_hi = c_hi->green; cm->blue_hi = c_hi->blue;
  printf("Colour Map low  end p = %f, colour = (%03i,%03i,%03i)\n",
	 v_lo, cm->red_lo, cm->green_lo, cm->blue_lo);
  printf("Colour Map high end p = %f, colour = (%03i,%03i,%03i)\n",
	 v_hi, cm->red_hi, cm->green_hi, cm->blue_hi);
  return;
}

/*
 * Make a new image of size (width, height), and initialise
 * background and object colours
 */

LBImage newlbimage(int width, int height, LBColour back, LBColour object) {
  int area;
  LBImage im;
  im = (LBImage) malloc(sizeof(LBImage_type));
  im->width = width; im->height = height; area = width * height;
  im->red = (int *) malloc(area*sizeof(int));
  im->green = (int *) malloc(area*sizeof(int));
  im->blue = (int *) malloc(area*sizeof(int));
  im->back_red = back->red;
  im->back_green = back->green; 
  im->back_blue = back->blue;
  im->object_red = object->red;
  im->object_green = object->green;
  im->object_blue = object->blue;
  printf("Image initialised to size %i * %i = %i\n", height, width, area);
  printf("Image background = (%03i,%03i,%03i), ",
	 im->back_red, im->back_green, im->back_blue);
  printf("objects = (%03i,%03i,%03i)\n", im->object_red, im->object_green, im->object_blue);
  return im;
}



/* Clear an image to its background */

void lbimageclear(LBImage im) {
  int i, j, ipixel, height, width;
  int *red, *green, *blue;
  height = im->height; width = im->width;
  red = im->red; green = im->green; blue = im->blue;
  for (i=0; i<width; i++) for (j=0; j<height; j++) {
    ipixel = height*i+j;
    red[ipixel] = im->back_red;
    green[ipixel] = im->back_green;
    blue[ipixel] = im->back_blue;
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) printf("Image cleared to background\n");
}


/* Clear an image to its background plus objects */

void lbimagemap(LBImage im, LBSys lb, char *dirn, int slice) {
  int i, j, k, inode, ipixel, esi, esj, esk, hoff, voff, height, width;
  int back_red, back_green, back_blue;
  int object_red, object_green, object_blue;
  int *red, *green, *blue;
  LBObj *obj;
  height = im->height; width = im->width;
  back_red = im->back_red;
  back_green = im->back_green;
  back_blue = im->back_blue;
  object_red = im->object_red;
  object_green = im->object_green;
  object_blue = im->object_blue;
  red = im->red; green = im->green; blue = im->blue;
  esi = lb->esi; esj = lb->esj; esk = lb->esk; obj = lb->obj;
/* Determine width and height offsets to centre image */
/* Insert objects from slice through LB system */  
  switch (dirn[0]) {
  case 'x': case 'i':
    hoff = (esj-width)/2; voff = (esk-height)/2;
    i = slice; for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
      inode = esk*(esj*i+j)+k; ipixel = height*(i+hoff) + (j+voff);
      if (obj[inode] == NULL) {
	red[ipixel] = back_red;
	green[ipixel] = back_green;
	blue[ipixel] = back_blue;
      } else {
	red[ipixel] = object_red;
	green[ipixel] = object_green;
	blue[ipixel] = object_blue;
      }
    }
    break;
  case 'y': case 'j':
    hoff = (esi-width)/2; voff = (esk-height)/2;
    j = slice; for (i=0; i<esi; i++) for (k=0; k<esk; k++) {
      inode = esk*(esj*i+j)+k; ipixel = height*(i+hoff) + (j+voff);
      if (obj[inode] == NULL) {
	red[ipixel] = back_red;
	green[ipixel] = back_green;
	blue[ipixel] = back_blue;
      } else {
	red[ipixel] = object_red;
	green[ipixel] = object_green;
	blue[ipixel] = object_blue;
      }
    }
    break;
  case 'z': case 'k':
    hoff = (esi-width)/2; voff = (esj-height)/2;
    k = slice; for (i=0; i<esi; i++) for (j=0; j<esj; j++) {
      inode = esk*(esj*i+j)+k; ipixel = height*(i+hoff) + (j+voff);
      if (obj[inode] == NULL) {
	red[ipixel] = back_red;
	green[ipixel] = back_green;
	blue[ipixel] = back_blue;
      } else {
	red[ipixel] = object_red;
	green[ipixel] = object_green;
	blue[ipixel] = object_blue;
      }
    }
    break;
  default: printf("Unrecognised direction %s in lbimageinit\n", dirn);
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("Image c_leared to background plus objects\n");
  }
}


/*
 * Set the colour in an image to correspond to a given velocity field
 * (magnitude or component), plus objects
 */

void lbimagevel(LBImage im, LBColourMap cm, LBSys lb, 
		char *what, char *dirn, int slice) {
  char w;
  int i, j, k, esi, esj, esk, inode, ipixel, hoff, voff, height, width;
  int object_red, object_green, object_blue;
  int red_lo, green_lo, blue_lo;
  int red_hi, green_hi, blue_hi;
  int red_local, green_local, blue_local;
  int *red, *green, *blue;
  double red_range, green_range, blue_range, sf, sf_local;
  LBReal t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14;
  LBReal rho0, fac, rho, ux, uy, uz;
  LBReal *aa;
  LBReal *fx, *fy, *fz;
  LBReal v_lo, v_hi, val;
  LBObj *obj;
  switch (what[0]) {
    case 'x': case 'i' : w = 'x'; break;
    case 'y': case 'j' : w = 'y'; break;
    case 'z': case 'k' : w = 'z'; break;
    default: w = 'm';
  }
  height = im->height; width = im->width;
  object_red = im->object_red;
  object_green = im->object_green;
  object_blue = im->object_blue;
  red_lo = cm->red_lo; green_lo = cm->green_lo; blue_lo = cm->blue_lo;
  red_hi = cm->red_hi; green_hi = cm->green_hi; blue_hi = cm->blue_hi;
  red_range = (double)(red_hi - red_lo);
  green_range = (double)(green_hi - green_lo);
  blue_range = (double)(blue_hi - blue_lo);
  v_lo = cm->v_lo; v_hi = cm->v_hi; sf = 1.0 / (v_hi - v_lo);
  red = im->red; green = im->green; blue = im->blue;
  esi = lb->esi, esj = lb->esj; esk = lb->esk;
  rho0 = lb->rho0; obj = lb->obj;
  aa = lb->a; fx = lb->fx; fy = lb->fy; fz = lb->fz;
  fac = 0.0;
  if (lb->vel_corrected) switch (lb->status) {
    case POST_PROPAGATE: fac = 0.5; break;
    case POST_COLLIDE: fac = -0.5; break;
  }
/* Insert colours from slice through LB system */  
  switch (dirn[0]) {
  case 'x': case 'i':
    hoff = (esj-width)/2; voff = (esk-height)/2;
    i = slice; for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
      inode = esk*(esj*i+j)+k; ipixel = height*(i+hoff) + (j+voff);
      if (obj[inode] == NULL) {
	t00 = aa[ 0+15*inode]; t01 = aa[ 1+15*inode]; t02 = aa[ 2+15*inode];
	t03 = aa[ 3+15*inode]; t04 = aa[ 4+15*inode]; t05 = aa[ 5+15*inode];
	t06 = aa[ 6+15*inode]; t07 = aa[ 7+15*inode]; t08 = aa[ 8+15*inode];
	t09 = aa[ 9+15*inode]; t10 = aa[10+15*inode]; t11 = aa[11+15*inode];
	t12 = aa[12+15*inode]; t13 = aa[13+15*inode]; t14 = aa[14+15*inode];
	rho = rho0 + t00+t01+t02+t03+t04+t05+t06+t07+t08+t09+t10+t11+t12+t13+t14;
	ux = (t01-t02+t07-t08+t09-t10+t11-t12+t13-t14)/rho;
	uy = (t03-t04+t07+t08-t09-t10+t11+t12-t13-t14)/rho;
	uz = (t05-t06+t07+t08+t09+t10-t11-t12-t13-t14)/rho;
	ux += fac*fx[inode]/rho; uy += fac*fy[inode]/rho; uz += fac*fz[inode]/rho;
	switch (w) {
   	  case 'x' : val = ux; break;
	  case 'y' : val = uy; break;
	  case 'z' : val = uz; break;
	  default : val = sqrt(ux*ux + uy*uy + uz*uz);
	}
	sf_local = sf * (val - v_lo);
	red_local = (int)(red_lo + sf_local * red_range);
	green_local = (int)(green_lo + sf_local * green_range);
	blue_local = (int)(green_lo + sf_local * blue_range);
	lbcolourclamp(&red_local,&green_local,&blue_local);
	red[ipixel] = red_local;
	green[ipixel] = green_local;
	blue[ipixel] = blue_local;
      } else {
	red[ipixel] = object_red;
	green[ipixel] = object_green;
	blue[ipixel] = object_blue;
      }
    }
    break;
  case 'y': case 'j':
    hoff = (esi-width)/2; voff = (esk-height)/2;
    j = slice; for (i=0; i<esi; i++) for (k=0; k<esk; k++) {
      inode = esk*(esj*i+j)+k; ipixel = height*(i+hoff) + (j+voff);
      if (obj[inode] == NULL) {
	t00 = aa[ 0+15*inode]; t01 = aa[ 1+15*inode]; t02 = aa[ 2+15*inode];
	t03 = aa[ 3+15*inode]; t04 = aa[ 4+15*inode]; t05 = aa[ 5+15*inode];
	t06 = aa[ 6+15*inode]; t07 = aa[ 7+15*inode]; t08 = aa[ 8+15*inode];
	t09 = aa[ 9+15*inode]; t10 = aa[10+15*inode]; t11 = aa[11+15*inode];
	t12 = aa[12+15*inode]; t13 = aa[13+15*inode]; t14 = aa[14+15*inode];
	rho = rho0 + t00+t01+t02+t03+t04+t05+t06+t07+t08+t09+t10+t11+t12+t13+t14;
	ux = (t01-t02+t07-t08+t09-t10+t11-t12+t13-t14)/rho;
	uy = (t03-t04+t07+t08-t09-t10+t11+t12-t13-t14)/rho;
	uz = (t05-t06+t07+t08+t09+t10-t11-t12-t13-t14)/rho;
	ux += fac*fx[inode]/rho; uy += fac*fy[inode]/rho; uz += fac*fz[inode]/rho;
	switch (w) {
   	  case 'x' : val = ux; break;
	  case 'y' : val = uy; break;
	  case 'z' : val = uz; break;
	  default : val = sqrt(ux*ux + uy*uy + uz*uz);
	}
	sf_local = sf * (val - v_lo);
	red_local = (int)(red_lo + sf_local * red_range);
	green_local = (int)(green_lo + sf_local * green_range);
	blue_local = (int)(green_lo + sf_local * blue_range);
	lbcolourclamp(&red_local,&green_local,&blue_local);
	red[ipixel] = red_local;
	green[ipixel] = green_local;
	blue[ipixel] = blue_local;
      } else {
	red[ipixel] = object_red;
	green[ipixel] = object_green;
	blue[ipixel] = object_blue;
      }
    }
    break;
  case 'z': case 'k':
    hoff = (esi-width)/2; voff = (esj-height)/2;
    k = slice; for (i=0; i<esi; i++) for (j=0; j<esj; j++) {
      inode = esk*(esj*i+j)+k; ipixel = height*(i+hoff) + (j+voff);
      if (obj[inode] == NULL) {
	t00 = aa[ 0+15*inode]; t01 = aa[ 1+15*inode]; t02 = aa[ 2+15*inode];
	t03 = aa[ 3+15*inode]; t04 = aa[ 4+15*inode]; t05 = aa[ 5+15*inode];
	t06 = aa[ 6+15*inode]; t07 = aa[ 7+15*inode]; t08 = aa[ 8+15*inode];
	t09 = aa[ 9+15*inode]; t10 = aa[10+15*inode]; t11 = aa[11+15*inode];
	t12 = aa[12+15*inode]; t13 = aa[13+15*inode]; t14 = aa[14+15*inode];
	rho = rho0 + t00+t01+t02+t03+t04+t05+t06+t07+t08+t09+t10+t11+t12+t13+t14;
	ux = (t01-t02+t07-t08+t09-t10+t11-t12+t13-t14)/rho;
	uy = (t03-t04+t07+t08-t09-t10+t11+t12-t13-t14)/rho;
	uz = (t05-t06+t07+t08+t09+t10-t11-t12-t13-t14)/rho;
	ux += fac*fx[inode]/rho; uy += fac*fy[inode]/rho; uz += fac*fz[inode]/rho;
	switch (w) {
   	  case 'x' : val = ux; break;
	  case 'y' : val = uy; break;
	  case 'z' : val = uz; break;
	  default : val = sqrt(ux*ux + uy*uy + uz*uz);
	}
	sf_local = sf * (val - v_lo);
	red_local = (int)(red_lo + sf_local * red_range);
	green_local = (int)(green_lo + sf_local * green_range);
	blue_local = (int)(green_lo + sf_local * blue_range);
	lbcolourclamp(&red_local,&green_local,&blue_local);
	red[ipixel] = red_local;
	green[ipixel] = green_local;
	blue[ipixel] = blue_local;
      } else {
	red[ipixel] = object_red;
	green[ipixel] = object_green;
	blue[ipixel] = object_blue;
      }
    }
    break;
  default: printf("Unrecognised direction %s in lbimagevel\n", dirn); return;
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("velocity field %s imaged\n",what);
  }
}


/*
 * Set the colour in an image to correspond to a given concentration field
 */

void mpimageset(LBImage im, MPSys mp, LBColourMap cm, char *what, char *dirn, int slice) {
  int i, j, k, esi, esj, esk, inode, ipixel, hoff, voff, addcolour, height, width;
  int object_red, object_green, object_blue;
  int red_lo, green_lo, blue_lo;
  int red_hi, green_hi, blue_hi;
  int red_local, green_local, blue_local;
  int *red, *green, *blue;
  double red_range, green_range, blue_range, sf, sf_local;
  LBReal *p, v_lo, v_hi;
  LBObj *obj;
  switch (what[0]) {
    case 'a': addcolour = 1; break;
    case 's': addcolour = 0; break;
    default: printf("Unrecognised add/subtract colour command %s\n", what); return;
  }
  height = im->height; width = im->width;
  object_red = im->object_red;
  object_green = im->object_green;
  object_blue = im->object_blue;
  red_lo = cm->red_lo; green_lo = cm->green_lo; blue_lo = cm->blue_lo;
  red_hi = cm->red_hi; green_hi = cm->green_hi; blue_hi = cm->blue_hi;
  red_range = (double)(red_hi - red_lo);
  green_range = (double)(green_hi - green_lo);
  blue_range = (double)(blue_hi - blue_lo);
  v_lo = cm->v_lo; v_hi = cm->v_hi; sf = 1.0 / (v_hi - v_lo);
  red = im->red; green = im->green; blue = im->blue;
  esi = mp->lb->esi, esj = mp->lb->esj; esk = mp->lb->esk;
  obj = mp->lb->obj; p = mp->p;
/* Insert colours from slice through LB system */  
  switch (dirn[0]) {
  case 'x': case 'i':
    hoff = (esj-width)/2; voff = (esk-height)/2;
    i = slice; for (j=0; j<esj; j++) for (k=0; k<esk; k++) {
      inode = esk*(esj*i+j)+k; ipixel = height*(i+hoff) + (j+voff);
      if (obj[inode] == NULL) {
	sf_local = sf * (p[inode] - v_lo);
	red_local = red[ipixel];
	green_local = green[ipixel];
	blue_local = blue[ipixel];
	if (addcolour) {
	  red_local += (int)(sf_local * red_range);
	  green_local += (int)(sf_local * green_range);
	  blue_local += (int)(sf_local * blue_range);
	} else {
	  red_local -= (int)(sf_local * red_range);
	  green_local -= (int)(sf_local * green_range);
	  blue_local -= (int)(sf_local * blue_range);
	}
	lbcolourclamp(&red_local,&green_local,&blue_local);
	red[ipixel] = red_local;
	green[ipixel] = green_local;
	blue[ipixel] = blue_local;
      } else {
	red[ipixel] = object_red;
	green[ipixel] = object_green;
	blue[ipixel] = object_blue;
      }
    }
    break;
  case 'y': case 'j':
    hoff = (esi-width)/2; voff = (esk-height)/2;
    j = slice; for (i=0; i<esi; i++) for (k=0; k<esk; k++) {
      inode = esk*(esj*i+j)+k; ipixel = height*(i+hoff) + (j+voff);
      if (obj[inode] == NULL) {
	sf_local = sf * (p[inode] - v_lo);
	red_local = red[ipixel];
	green_local = green[ipixel];
	blue_local = blue[ipixel];
	if (addcolour) {
	  red_local += (int)(sf_local * red_range);
	  green_local += (int)(sf_local * green_range);
	  blue_local += (int)(sf_local * blue_range);
	} else {
	  red_local -= (int)(sf_local * red_range);
	  green_local -= (int)(sf_local * green_range);
	  blue_local -= (int)(sf_local * blue_range);
	}
	lbcolourclamp(&red_local,&green_local,&blue_local);
	red[ipixel] = red_local;
	green[ipixel] = green_local;
	blue[ipixel] = blue_local;
      } else {
	red[ipixel] = object_red;
	green[ipixel] = object_green;
	blue[ipixel] = object_blue;
      }
    }
    break;
  case 'z': case 'k':
    hoff = (esi-width)/2; voff = (esj-height)/2;
    k = slice; for (i=0; i<esi; i++) for (j=0; j<esj; j++) {
      inode = esk*(esj*i+j)+k; ipixel = height*(i+hoff) + (j+voff);
      if (obj[inode] == NULL) {
	sf_local = sf * (p[inode] - v_lo);
	red_local = red[ipixel];
	green_local = green[ipixel];
	blue_local = blue[ipixel];
	if (addcolour) {
	  red_local += (int)(sf_local * red_range);
	  green_local += (int)(sf_local * green_range);
	  blue_local += (int)(sf_local * blue_range);
	} else {
	  red_local -= (int)(sf_local * red_range);
	  green_local -= (int)(sf_local * green_range);
	  blue_local -= (int)(sf_local * blue_range);
	}
	lbcolourclamp(&red_local,&green_local,&blue_local);
	red[ipixel] = red_local;
	green[ipixel] = green_local;
	blue[ipixel] = blue_local;
      } else {
	red[ipixel] = object_red;
	green[ipixel] = object_green;
	blue[ipixel] = object_blue;
      }
    }
    break;
  default: printf("Unrecognised direction %s in mpimage\n", dirn); return;
  }
  if (slb_verbosity & TRACE_ROUTINE_EXITS) {
    printf("MP Concentration field colour %s image\n",
	   addcolour ? "added to" : "subtracted from");
  }
}


/* Extract the colour of a given pixel in an image */

void lbcolourget(LBImage im, LBColour clr, int i, int j) {
  int k;
  if (i < 0 || i > im->width || j < 0 || j > im->height) {
    printf("lbcolourget: out of range (i, j) = (%i, %i)\n",i,j);
    printf("lbcolourget: range = (0, 0) -- (%i, %i)\n", im->width, im->height);
    return;
  }
  k = im->height * i + j;
  clr->red = im->red[k]; clr->green = im->green[k]; clr->blue = im->blue[k];
  return;
}

int lbgetred(LBImage im, int i, int j) {
  if (i < 0 || i > im->width || j < 0 || j > im->height) {
    printf("lbgetred: out of range (i, j) = (%i, %i)\n",i,j);
    printf("lbgetred: range = (0, 0) -- (%i, %i)\n", im->width, im->height);
    return 0;
  }
  return im->red[im->height*i + j];
}

int lbgetgreen(LBImage im, int i, int j) {
  if (i < 0 || i > im->width || j < 0 || j > im->height) {
    printf("lbgetgreen: out of range (i, j) = (%i, %i)\n",i,j);
    printf("lbgetgreen: range = (0, 0) -- (%i, %i)\n", im->width, im->height);
    return 0;
  }
  return im->green[im->height*i + j];
}

int lbgetblue(LBImage im, int i, int j) {
  if (i < 0 || i > im->width || j < 0 || j > im->height) {
    printf("lbgetblue: out of range (i, j) = (%i, %i)\n",i,j);
    printf("lbgetblue: range = (0, 0) -- (%i, %i)\n", im->width, im->height);
    return 0;
  }
  return im->blue[im->height*i + j];
}

/* End of mpimage.c */
