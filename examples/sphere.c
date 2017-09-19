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
 *
 *
 * This example computes the flow around a single sphere in a periodic
 * box, in other words the flow in an infinite dilute cubic array of
 * spheres.  It can be used to determine the hydrodynamic radius of the
 * sphere.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sunlightlb.h>

int main () {
  int i, init;
  float ox, oy, oz, rad;
  LBSys mylb;

/* Initialisation stuff */

  about();

  mylb = newlb(16, 16, 16, "sphere_c");

  set_verbosity(OBJECT_DETAILS);

  init = 1000;

  mylb->stokes_flow = 1;
  mylb->objects_moving = 0;
  mylb->skip_solid = 1;
  mylb->tau = 0.80;
  mylb->dpdx = 0.01;
  mylb->dpdy = 0.0;
  mylb->dpdz = 0.0;

  wrlinks(mylb);
  
/* Insert some spheres */

  ox = oy = oz = 7.5; rad = 3.7;
  insertsph(mylb, rad, ox, oy, oz);
  wrmap(mylb, "sphere_c.map");
  lbsavemap(mylb, "sphere_c.dat","list");

/* Equilibrate for init steps */

  wrpars(mylb);

  wrmon(mylb, "initialised");

  printf("\nEquilibrating for %i steps\n", init);
  for (i=0; i<init; i++) lbstep(mylb);

  wrstats(mylb);
  measure(mylb); 
  wrmeas(mylb);
  wrobjs(mylb);
  lbsavefrc(mylb, "sphere_c.frc", "list");
  lbsavevel(mylb, "sphere_c.vel", "ascii");

  wrmon(mylb, "finished");

}
