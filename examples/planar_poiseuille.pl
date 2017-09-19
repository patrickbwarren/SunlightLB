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
# This example computes planar Poiseuille flow between walls at z = 0
# and z = L - 1, where the box is periodic in z with a length L.  The
# velocity field, stress field, and the forces arising on the walls,
# are exact in this geometry, in lattive Boltzmann.
#
#########################################################################

use SunlightLB;
use strict;

my ($name, $mylb, $nsteps, $esz, $i, $eta);
my ($wall1, $wall2);
my ($f1, $f2, $h, $umean, $dpdx);

SunlightLB::about();

$name = "planar_poiseuille";

$esz = 20;
$umean = 0.1;

$mylb = SunlightLB::newlb(1, 1, $esz, $name);

SunlightLB::set_verbosity($SunlightLB::OBJECT_DETAILS);

$nsteps = 5000;

$mylb->{stokes_flow} = 1;
$mylb->{objects_moving} = 0;
$mylb->{tau} = 0.80;

$eta = $mylb->{rho0} * ($mylb->{tau} - 0.5) / 3.0;

$mylb->{dpdx} = $dpdx = 12*$eta*$umean/(($esz-2)*($esz-2));

# Insert two walls

$wall1 = $mylb->insertwall("z", 0);
$wall2 = $mylb->insertwall("z", $esz-1);

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

$f1 = $wall1->{frcx}; $f2 = $wall2->{frcx};
$h = $wall2->{posz} - $wall1->{posz} - 1;
$umean = $h * $h * $dpdx / (12 * $eta);

print "\nResults:\n";

printf " set h = %f, eta = %f, dpdx = %f, umean = %f\n", 
    $h, $eta, $dpdx, $umean;

printf " calculated wall forces dpdx * h / 2 = %f, <ux> = %f\n", 
    $dpdx * $h / 2, $umean * $h / $esz;

printf " measured f1 = %f, f2 = %f, <ux> = %f\n", 
    $f1, $f2, $mylb->{vmeanx};

