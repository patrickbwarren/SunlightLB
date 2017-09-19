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
# This example computes the drag on a simple cubic array of
# close-packed spheres, as a function of Reynolds number.  It can be
# compared with Table 6 in Ladd, JFM 271, 311-339 (1994).  The
# comparison is not exact since Ladd uses a slightly different LB
# scheme.  Here are the results that you should obtain
#
# ---------------------------------------------------------------------
#  grad_p    u_mean       F_D        Re      F_D / (6 pi eta a u_mean)
# ---------------------------------------------------------------------
#   0.05    0.049391    245.65      0.00     41.39
#   0.01    0.008485     49.13     29.07     48.18
#   0.02    0.016008     98.26     54.84     51.08
#   0.03    0.023410    147.39     80.20     52.39
#   0.04    0.030826    196.52    105.60     53.05
# ---------------------------------------------------------------------
#
# As an exercise you can try (i) a larger lattice size, say 33^3.
# Note that the pressure gradient should be reduced (eg by a factor
# 10), and the number of steps to reach steady state should be
# increased, eg to 10k; (ii) try having a pressure gradient which is
# not aligned along the axes.  
#
#########################################################################

use SunlightLB;

SunlightLB::about();

# Routine to initialise the lattice Boltzmann (LB) simulation.

sub SetupLB {

    $mylb = SunlightLB::newlb($es, $es, $es, "sccp_sphere");

    $mylb->{objects_moving} = 0;
    $mylb->{vel_corrected} = 1;
    $mylb->{skip_solid} = 1;

    $mylb->{tau} = 0.53125;

    $ox = $oy = $oz = ($es - 1.0) / 2.0; 
    $a_nom = ($es - 0.4) / 2.0;
    $a_hyd = $es / 2.0;

    $mylb->insertsph($a_nom, $ox, $oy, $oz);

    $nu = ($mylb->{tau} - 0.5) / 3.0;
    $eta = $mylb->{rho0} * $nu; 
    $eps = 1.0 - $SunlightLB::PI / 6.0;

    printf "Initialised to %i^3 lattice\n", $es;
    printf "Sphere at %0.2f %0.2f %0.2f, a_nom = %0.2f, a_hyd = %0.2f\n",
    $ox, $oy, $oz, $a_nom, $a_hyd;

}

# Routine to undertake an LB simulation, and add a line of results to
# the summary string.  The input value is the pressure gradient.

sub RunLB {

    $mylb->{dpdx} = $_[0];

    $mylb->wrpars;

    print "Equilibrating for $nsteps steps\n";

    for ($i=0; $i<$nsteps; $i++) { $mylb->lbstep; }

    $mylb->wrstats;
    $mylb->measure;
    $mylb->wrmeas;

    $u_mean = $mylb->{vmeanx} * $mylb->{vmeanx};
    $u_mean += $mylb->{vmeany} * $mylb->{vmeany};
    $u_mean += $mylb->{vmeanz} * $mylb->{vmeanz};
    $u_mean = sqrt($u_mean);

    $grad_p = $mylb->{dpdx} * $mylb->{dpdx};
    $grad_p += $mylb->{dpdy} * $mylb->{dpdy};
    $grad_p += $mylb->{dpdz} * $mylb->{dpdz};
    $grad_p = sqrt($grad_p);

    $F_D = $grad_p * $mylb->{nnodes};
    $xi_D = $F_D / (6 * $SunlightLB::PI * $eta * $a_hyd *$u_mean);

    if ($mylb->{stokes_flow} == 1) {
	$Re = 0;
    } else {
	$Re = 2.0 * $u_mean * $a_hyd / ($eps * $nu);
    }

    $summary .= sprintf "%6.2f  %10.6f  %8.2f  %8.2f  %8.2f\n",
                        $grad_p, $u_mean, $F_D, $Re, $xi_D;

}

# This is where the action all starts.

$es = 17; $nsteps = 5000;

SetupLB; 

# First do a zero Reynolds number calculation.

$mylb->wrmon("initialised");

$mylb->{stokes_flow} = 1;

RunLB(0.05);

# Now do non-zero Reynolds number calculations for a set of grad_p values.

$mylb->{stokes_flow} = 0;

foreach $val ( "0.01", "0.02", "0.03", "0.04" ) { RunLB($val); }

$mylb->wrmon("finished");

# Print the final summary table

$horiz_line = '-' x 69 . "\n";

print "\n\nSUMMARY\n", $horiz_line;
print " grad_p    u_mean       F_D        Re      F_D / (6 pi eta a u_mean)\n";
print $horiz_line, $summary, $horiz_line, "\n";


