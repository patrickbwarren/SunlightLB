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
# This example computes the flow around a single sphere in a periodic
# box, in other words the flow in an infinite dilute cubic array of
# spheres.  It can be used to determine the hydrodynamic radius of the
# sphere.
#
#########################################################################

use SunlightLB;
use strict;

my ($name, $mylb, $nsteps, $rad, $ox, $oy, $oz, $i, $temp);

SunlightLB::about();

$name = "sphere_perl";

$mylb = SunlightLB::newlb(16, 16, 16, $name);

SunlightLB::set_verbosity($SunlightLB::OBJECT_DETAILS);

$nsteps = 1000;

$mylb->{stokes_flow} = 1;
$mylb->{objects_moving} = 0;
$mylb->{skip_solid} = 1;
$mylb->{tau} = 0.80;
$mylb->{dpdx} = 0.01;
$mylb->{dpdy} = 0.0;
$mylb->{dpdz} = 0.0;

$mylb->wrlinks;

# Insert a sphere

$ox = $oy = $oz = 7.5; $rad = 3.7;
$temp = $mylb->insertsph($rad, $ox, $oy, $oz);
$mylb->wrmap("$name.map");
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

# Note how useful perl is to generate the file names.

$mylb->wrmon("finished");

