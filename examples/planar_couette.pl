#!/usr/bin/perl

# This file is part of SunlightLB - a 3D Lattice Boltzmann code
# Copyright (C) 2005 Unilever UK Central Resources Ltd.

# SunlightLB is free software; you can redistribute it and/or 
# modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation; either version 2 
# of the License, or (at your option) any later version. 

# SunlightLB is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

# You should have received a copy of the GNU General Public License 
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 

#########################################################################
#
# This example computes planar Couette flow by having moving walls at
# z = 0 and z = h + 1, thus the separation between the walls defined
# by the link bounceback is h.  The box is periodic in z with a length
# 2(h + 1) thus the flow profile is duplicated but reversed in
# gradient, in the second half of the box.  The velocity field, stress
# field, and the forces arising on the walls, are exact in this
# geometry, in lattive Boltzmann.
#
#########################################################################

use SunlightLB;
use strict;

my ($name, $mylb, $nsteps, $esz, $i, $eta);
my ($wall1, $wall2);
my ($vel1, $vel2, $f1, $f2, $h);

SunlightLB::about();

$name = "planar_couette";

$vel1 = 0.11; 
$vel2 = -0.07;
$h = 10;
$esz = 2*($h+1);

$mylb = SunlightLB::newlb(1, 1, $esz, $name);

SunlightLB::set_verbosity($SunlightLB::OBJECT_DETAILS);

$nsteps = 5000;

$mylb->{stokes_flow} = 1;
$mylb->{objects_moving} = 1;
$mylb->{tau} = 0.80;

# Insert two walls

$wall1 = $mylb->insertwall("z", 0);
$wall2 = $mylb->insertwall("z", $h+1);

$wall1->{velx} = $vel1;
$wall2->{velx} = $vel2;

$mylb->lbsavemap("$name.dat", "list");

# Equilibrate for nsteps steps

$mylb->wrpars;

$mylb->wrmon("initialised");

print "\nEquilibrating for $nsteps steps\n";

for ($i=0; $i<$nsteps; $i++) { $mylb->lbstep; }

# Measure some things and save some things

$mylb->wrstats;
$mylb->measure;
$mylb->wrmeas;
$mylb->wrobjs;
$mylb->lbsavefrc("$name.frc", "list");
$mylb->lbsavevel("$name.vel", "ascii");
$mylb->lbsavestr("$name.str", "ascii");

$mylb->wrmon("finished");

$eta = $mylb->{rho0} * ($mylb->{tau} - 0.5) / 3.0;
$f1 = $wall1->{frcx}; $f2 = $wall2->{frcx};

print "\nResults:\n";
printf " set h = %f, vel1 = %f, vel2 = %f, eta = %f\n", $h, $vel1, $vel2, $eta;
printf " calculated velocity gradient (vel1 - vel2) / h = %f\n", ($vel1 - $vel2) / $h;
printf " calculated wall forces 2 eta (vel1 - vel2) / h = %f\n", 2 * $eta * ($vel1 - $vel2) / $h;
printf " calculated shear stress eta (vel1 - vel2) / h = %f\n", $eta * ($vel1 - $vel2) / $h;
printf " measured f1 = %f, f2 = %f\n", $f1, $f2;


