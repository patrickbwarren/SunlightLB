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
# This example computes the permeability of a periodic box of 144
# random close packed spheres.  The box is 5^3 in units of the sphere
# diameter, so you can work out that the solid volume fraction is 144
# * pi / (6 * 5^3) = 0.603, so this is not quite at the r.c.p. limit.
#
# This code saves the state of the LB system at the end of the
# equilibration run.  It loads this in again if called subsequently
# (provided the names match).
#
# The results should look something like the following. 
#
# RESULTS
# Pressure gradient = 0.150, viscosity = 7.200
# Mean superficial flow rate = 0.001607
# D'Arcy law constant (hydraulic conductance) k = 0.077122
# Pore volume fraction = 0.40
# Equivalent capillary diameter = 2.50
#
# As an exercise, you can try running this with different LB sphere
# diameters, with off-sets (effectively generating subtly different
# mappings onto the lattice), and so on
#
#########################################################################

use SunlightLB;

SunlightLB::about();

( $root, $sigma, $dirn, $grad_p ) = ( "bpack", "9.0", "x", "0.15" );

( $xoff, $yoff, $zoff ) = ( 0.0, 0.0, 0.0 );

$nsteps = 1000;

$head = "$root-$sigma";
$dms_file = "$head-$dirn-$grad_p.dms";

$sphere_data_file = "rcp_spheres.dat";
$len = 5.0;

# Compute the new edge size, and a scale factor.

$es = int(0.5 + $len * $sigma);
$scale_fac = $es / $len;
$a_nom = $sigma / 2.0;

printf "Scale factor = %f\n", $scale_fac;

$mylb = SunlightLB::newlb($es, $es, $es, $head);

$mylb->{stokes_flow} = 1;
$mylb->{skip_solid} = 1;
$mylb->{tau} = 0.80;


# Perl hack - the eval command evaluates its argument, thus setting
# the pressure gradient in the appropriate direction.

eval "\$mylb->\{dpd$dirn\} = \$grad_p;";

# Insert spheres from data file, and flatten the resulting structure.

open(FIN, $sphere_data_file) || die "Cannot open $sphere_data_file: $!";

$nspheres = 0;

while ($line = <FIN>) {

    chop($line);

    ( $r, $x, $y, $z ) = split(' ', $line);

    while ($x < 0.0) { $x += $len; }
    while ($y < 0.0) { $y += $len; }
    while ($z < 0.0) { $z += $len; }
    while ($x >= $len) { $x -= $len; }
    while ($y >= $len) { $y -= $len; }
    while ($z >= $len) { $z -= $len; }

    $x = $x * $scale_fac + $xoff; 
    $y = $y * $scale_fac + $yoff; 
    $z = $z * $scale_fac + $zoff;

    for ($ix=-1; $ix<=1; $ix++) {
	for ($iy=-1; $iy<=1; $iy++) {
	    for ($iz=-1; $iz<=1; $iz++) {
		$ox = $x + $ix * $es; $oy = $y + $iy * $es; $oz = $z + $iz * $es;
		if ($ox > -$sigma && $ox < $es + $sigma &&
		    $oy > -$sigma && $oy < $es + $sigma &&
		    $oz > -$sigma && $oz < $es + $sigma) {
		    $mylb->insertsph($a_nom, $ox, $oy, $oz);
		}
	    }
	}
    }
    
    $nspheres++;
    
}

close(FIN);

$mylb->flattenobj;

print "Inserted $nspheres spheres\n";

$mylb->lbsavemap("$head.lst", "list");
$mylb->wrmap("$head.map");

$mylb->wrpars;

$mylb->wrmon("initialised");

# Load existing data if it has been saved.

if (-e $dms_file) {

    print "Loading dms data from $dms_file\n";
    $mylb->lbload($dms_file);

}

# Equilibrate for nsteps steps.

print "Equilibrating for $nsteps steps\n";
for ($i=0; $i<$nsteps; $i++) { $mylb->lbstep; }

print "Saving dms data to $dms_file\n";
$mylb->lbsave($dms_file);

$mylb->wrstats;
$mylb->measure;
$mylb->wrmeas;

$mylb->wrmon("finished");

# Calculate viscosity, d'Arcy law constant, and equivalent capillary
# diameter via k = R^2 / 32 * pore_vol_frac.

$eta = $mylb->{rho0} * ($mylb->{tau} - 0.5) / 3.0; 

eval "\$v_mean = \$mylb->\{vmean$dirn\};";

$k_darcy = $v_mean * $eta / $grad_p;
$pore_vf = $mylb->countfluid / $mylb->{nnodes};
$equiv_cap_diam = sqrt(32.0 * $k_darcy / $pore_vf);

printf "\n\nRESULTS\n";
printf "Pressure gradient = %0.3f, viscosity = %0.3f\n", $grad_p, $eta;
printf "Mean superficial flow rate = %f\n", $v_mean;
printf "D'Arcy law constant (hydraulic conductance) k = %f\n", $k_darcy;
printf "Pore volume fraction = %0.2f\n", $pore_vf;
printf "Equivalent capillary diameter = %0.2f\n", $equiv_cap_diam;
printf "\n\n";

