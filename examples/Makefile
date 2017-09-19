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

CC = /usr/bin/gcc
RM = /bin/rm

CC_FLAGS = -O3 -Wall -ffast-math $(EXTRA_CC_FLAGS)

# Shouldn't need to change anything below here.

sphere : sphere.o 
	$(CC) $(CC_FLAGS) sphere.o $(LIB_FLAGS) -lsunlightlb -lm -o sphere

clean : 
	$(RM) -f *~ *.o

pristine : clean
	$(RM) -f sphere
	$(RM) -f *.mon
	$(RM) -f sphere*.map sphere*.dat sphere*.frc sphere*.vel
	$(RM) -f phi*.lst phi*.vel
	$(RM) -f bpack-*.*

### Rule for compiling source code to object code ###

.c.o:
	$(CC) $(CC_FLAGS) $(INC_FLAGS) -c $<

### End of Makefile ###
