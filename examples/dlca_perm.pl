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
# This example computes the permeability of a DLCA spanning cluster,
# contained in the dataset phi0.07_L40.dat.gz.  Some pre-processing of
# the data is required.  The name means that 7% of the sites are blocked,
# in a 40^3 lattice.
#
# The DLCA data is courtesy Mark D. Haw and R. M. L. (Mike) Evans.
# See RMLE and MDH, Europhys. Lett. v60, 404 (2002).
#
# The results should look something like the following.  Note that this
# matrix should be symmetric.
#
# D'Arcy permeability matrix is
# k_xx k_yx k_zx =    3.0551   -0.0368    0.0125
# k_xy k_yy k_zy =   -0.0367    2.9599   -0.0700
# k_xz k_yz k_zz =    0.0126   -0.0700    3.0054
#
#########################################################################

use SunlightLB;

SunlightLB::about();

$name = "phi0.07_L40";

$es = 40; $nsteps = 1000; $grad_p = 0.03;

# We have to unpack the data file, and perform a little pre-processing
# (the node location numbering starts from 1 rather than 0, so we have
# to subtract 1; this is done with gawk).  Note the utility of perl to
# handle this step, or you can do it by hand of course.

sub PrepareDataFile {

    print "DLCA data courtesy Mark D. Haw and R. M. L. (Mike) Evans\n";
    print "See RMLE and MDH, Europhys. Lett. v60, 404 (2002)\n";

    $gawk_exe = "/usr/bin/gawk";
    $gzip_exe = "/usr/bin/gzip";
    $dat_file = "$name.dat";
    $compressed_dat_file = "$dat_file.gz";
    $processed_dat_file = "$name.lst";
    $gawk_command = "NF == 4 \{ printf \"%3i %3i %3i\\n\", \$1-1, \$2-1, \$3-1 \}";
    $command = "$gzip_exe -d -c $compressed_dat_file | ";
    $command .= "$gawk_exe '$gawk_command' > $processed_dat_file";
    print "\nExecuting the following:\n$command\n\n";
    system($command);
}

# Set up the lattice Boltzmann (LB) simulation.

sub SetupLB {

    $mylb = SunlightLB::newlb($es, $es, $es, $name);

    SunlightLB::set_verbosity($SunlightLB::OBJECT_DETAILS);

    $mylb->{objects_moving} = 0;
    $mylb->{vel_corrected} = 1;
    $mylb->{skip_solid} = 1;
    $mylb->{stokes_flow} = 1;
    $mylb->{tau} = 0.7;

    $mylb->insertlist($processed_dat_file);

    $eta = $mylb->{rho0} * ($mylb->{tau} - 0.5) / 3.0;

    $k_matrix = "";

}

# Run the LB simulation with the pressure gradient in the appropriate
# direction.  There is a neat perl hack in here - the eval command
# evaluates its argument, thus setting the pressure gradient in the
# appropriate direction.

sub RunLB {

    $dirn = $_[0];

    $mylb->{dpdx} = $mylb->{dpdy} = $mylb->{dpdz} = 0.0;

    eval "\$mylb->\{dpd$dirn\} = \$grad_p;";

    $mylb->wrpars;

    print "Equilibrating for $nsteps steps\n";
    for ($i=0; $i<$nsteps; $i++) { $mylb->lbstep; }

    $mylb->wrstats;
    $mylb->measure;
    $mylb->wrmeas;

    $k_matrix .= sprintf "k_x$dirn k_y$dirn k_z$dirn = %9.4f %9.4f %9.4f\n",
    ($mylb->{vmeanx}) * $eta / $grad_p,
    ($mylb->{vmeany}) * $eta / $grad_p,
    ($mylb->{vmeanz}) * $eta / $grad_p;

    $vel_file = "$name.dpd$dirn.vel";
    $mylb->lbsavevel($vel_file, "ascii");

}

# This is where all the action takes place.

PrepareDataFile;

SetupLB;

$mylb->wrmon("initialised");

RunLB("x");
RunLB("y");
RunLB("z");

$mylb->wrmon("finished");

# Print final permeability matrix.

print "\n\nD'Arcy permeability matrix is\n", $k_matrix, "\n\n";
