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
# This example computes the lateral mobility of a sphere above a plane
# wall similar to mobility1.pl.  In this case, the mobility is
# computed directly by relaxing to a zero-torque condition on the fly.
#
# The results should be like the first two columns in the following
# table.
#
# ------------------------------------------------------
#  a_hyd / h   D / D_0     mobility1     D(exact) / D_0
# ------------------------------------------------------
#    0.820      0.522       0.522           0.5233
#    0.683      0.620       0.619           0.6168
#    0.586      0.682       0.680           0.6746
#    0.512      0.727       0.725           0.7165
#    0.456      0.762       0.760           0.7476
# ------------------------------------------------------
#
# The third column is from the resistance matrix method.  The fourth
# column is the Cichocki and Jones exact results.
#
#########################################################################

use SunlightLB;

SunlightLB::about();

( $esx, $esy, $esz, $name ) = ( 48, 32, 32, "mobility2" );
( $nsteps, $a_hyd, $a_nom ) = ( 2000, 4.1, 4.3 );

# Set up a standard LB system.

sub SetupLB {

    $mylb = SunlightLB::newlb($esx, $esy, $esz, $name);

    SunlightLB::set_verbosity($SunlightLB::OBJECT_DETAILS);

    $mylb->{objects_moving} = 1;
    $mylb->{vel_corrected} = 1;
    $mylb->{skip_solid} = 0;
    $mylb->{stokes_flow} = 1;
    $mylb->{tau} = 0.95;

    $eta = $mylb->{rho0} * ($mylb->{tau} - 0.5) / 3.0;

}

# Set up objects in the LB system - a wall at the base and a sphere
# with a variable height.  The height is the height above the wall
# stick boundary condition, which is located at 1/2 lattice spacing
# above the wall nodes.  The sphere is centred midway between nodes.

sub SetupObjects {

    $mywall = $mylb->insertwall("z", 0);
    $mysphere = $mylb->insertsph($a_nom, $esx/2-0.5, $esy/2-0.5, $height+0.5);

}

# Take as input the velocity and angular velocity of the sphere.  Set
# vx and ay to these values, run the LB calculation, measure the
# resulting force and torque, and set fx and ty to the measured
# values.

sub RunLB {

    $mysphere->{velx} = $vx;
    $mysphere->{angy} = $ay;

    print "\nRunning for $num_steps steps\n";

    for ($i=0; $i<$num_steps; $i++) { $mylb->lbstep; }

    $mylb->wrstats;
    $mylb->measure;
    $mylb->wrmeas;
    $mylb->wrobjs;

    $fx = $mysphere->{frcx};
    $ty = $mysphere->{trqy};

}

# Compute the mobility of a sphere above a plane, normalised by the
# Stokes drag.  The hydrodynamic radius should be computed elsewhere.

sub ComputeMobility {

    SetupObjects(@_);

    $mylb->wrpars;

    $vx = 0.1;

    $ntries = 20;
    $num_steps = $nsteps / $ntries;

#  We know that T_y is linear in Omega_y, and if we know T_y for two
#  values of Omega_y, we can find the value of Omega_y for which T_y
#  is zero.  The idea is to iterate this calculation, as the LB fluid
#  is evolving, and *hopefully* the scheme is sufficiently numerically
#  stable to converge to the right value of Omega_y.  This is not very
#  sophisticated stuff, folks!

    for ($j=0; $j<$ntries; $j++) {
	if ($j == 0) { $ay = 0.0; }
	elsif ($j == 1) { $ay = 0.01; }
	else { $ay = $ay_1 - $ty_1 * ($ay_1 - $ay_0) / ($ty_1 - $ty_0); }
	RunLB;
	if ($j == 0) { $ay_0 = $ay; $ty_0 = $ty; }
	elsif ($j == 1) { $ay_1 = $ay; $ty_1 = $ty; }
	else {
	    $ay_0 = $ay_1; $ay_1 = $ay;
	    $ty_0 = $ty_1; $ty_1 = $ty;
	}
    }

    $mylb->removeobjs;

    $M_vx_fx = $vx / $fx;

    $DbyD0 = - 6 * $SunlightLB::PI *$eta * $a_hyd * $M_vx_fx;

    printf "M_vx_fx = %f\n", $M_vx_fx;
    printf "height = %f  DbyD0 = %f\n", $height, $DbyD0;

    $summary .= sprintf "%8.3f   %8.3f\n", $a_hyd/$height, $DbyD0;

}

# This is where the action all starts.

SetupLB;

$mylb->wrmon("initialised");

foreach $height ( "5.0", "6.0", "7.0", "8.0", "9.0" ) { ComputeMobility; }

$mylb->wrmon("finished");

# Print the final summary table

$horiz_line = '-' x 30 . "\n";

print "\n\nSUMMARY\n", $horiz_line, " a_hyd / h   D / D_0\n";
print $horiz_line, $summary, $horiz_line, "\n";
