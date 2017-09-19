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
# wall.  This is done by computing the resistance matrix and inverting
# it to get the mobilities.  Symmetry means that we only have to
# consider 2x2 matrices coupling (v_x, omega_y) to (f_x, t_y) where
# omega is the angular velocity, t is the torque.  The calculation is
# repeated for a number of heights.  This example introduces some new
# LB elements, such as moving objects, and removal / addition of objects
# on the fly.
#
# The results should be like the first two columns in the following
# table.
#
# ------------------------------------------
#  a_hyd / h   D / D_0       D(exact) / D_0
# ------------------------------------------
#    0.820      0.522           0.5233
#    0.683      0.619           0.6168
#    0.586      0.680           0.6746
#    0.512      0.725           0.7165
#    0.456      0.760           0.7476
# ------------------------------------------
#
# The third column is essentially the exact mobility (to with 1 part
# in 2000) from Cichocki and Jones, Physica A v258, 273 (1998).  It is
# seen that the comparison is pretty good.
#
# As an exercise, try varying the viscosity (eg tau = 0.8).  Repeat
# the calculation for a sphere above a step (terrace) - requires
# inserting slabs and so on.
#
#########################################################################

use SunlightLB;

SunlightLB::about();

( $esx, $esy, $esz, $name ) = ( 48, 32, 32, "mobility1" );
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

    ( $vx, $ay ) = @_;

    $mysphere->{velx} = $vx;
    $mysphere->{angy} = $ay;

    print "\nEquilibrating for $nsteps steps\n";
    for ($i=0; $i<$nsteps; $i++) { $mylb->lbstep; }

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
    RunLB(0.1, 0.00); $R_fx_vx = $fx / $vx; $R_ty_vx = $ty / $vx;
    RunLB(0.0, 0.02); $R_fx_ay = $fx / $ay; $R_ty_ay = $ty / $ay;
    $mylb->removeobjs;

    $M_vx_fx = $R_ty_ay / ($R_fx_vx * $R_ty_ay - $R_fx_ay * $R_ty_vx);

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
