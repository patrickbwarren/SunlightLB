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
# This example computes Poiseuille flow in a cylindrical pipe.  This
# is a good test of the implementation of the boundary conditions
# since the flow profile is parabolic, which ought to be captured very
# well by lattice Boltzmann.
#
#########################################################################

use SunlightLB;
use strict;

my ($name, $mylb, $nsteps, $es, $i, $eta);
my ($pipe, $f, $rad, $org, $umean, $dpdz, $area);
my $pi = 3.14159265359;

SunlightLB::about();

$name = "poiseuille";

$es = 45;
$umean = 0.1;

$mylb = SunlightLB::newlb($es, $es, 1, $name);

SunlightLB::set_verbosity($SunlightLB::OBJECT_DETAILS);

$nsteps = 20000;

$mylb->{stokes_flow} = 1;
$mylb->{objects_moving} = 0;
$mylb->{skip_solid} = 1;
$mylb->{tau} = 0.80;

$org = 0.5 * ($es - 1);
$rad = $org - 1.5;
$pipe = $mylb->insertcyl("zo", $rad, $org, $org);

$mylb->wrmap("$name.map");

$eta = $mylb->{rho0} * ($mylb->{tau} - 0.5) / 3.0;

$mylb->{dpdz} = $dpdz = (8*$eta*$umean)/($rad*$rad);

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

$area = $pi * $rad * $rad;

print "\nResults:\n";

printf " pipe radius = %f, centred at (%f, %f), lattice size %i^2\n", $rad, $org, $org, $es;

printf " expected area = %f\n", $area;
printf " measured area = %f\n", $mylb->countfluid;

printf " measured mean velocity = %f\n", $mylb->{vmeanz} * $es * $es / $area;
printf " expected mean velocity = %f\n", $umean;

printf " measured force = %f\n", $pipe->{frcz};
printf " expected force = %f\n", $dpdz * $area;

printf " expected wall stress = %f\n", $dpdz * $rad / 2;
