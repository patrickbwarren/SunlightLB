# This file is part of SunlightLB - a 3D Lattice Boltzmann code
# Copyright (C) 200-2010 Unilever UK Central Resources Ltd.

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

# May need to edit the following
#  perl -e 'use Config; print "$Config{ccflags}\n";'
#  perl -e 'use Config; print "$Config{archlib}\n";'


CC = /usr/bin/gcc
RM = /bin/rm
LN = /bin/ln
INSTALL = /usr/bin/install
LDCONFIG = /sbin/ldconfig
SWIG = /usr/bin/swig
PERL = /usr/bin/perl
PYTHON = /usr/bin/python2.4

LB_INC = -I.
LB_LIB = -L. -lsunlightlb

INSTALL_DIR_FOR_LIB = /usr/local/lib/
INSTALL_DIR_FOR_INC = /usr/local/include/
INSTALL_DIR_FOR_BIN = /usr/local/bin/

CC_FLAGS = -O3 -Wall -ffast-math $(EXTRA_CC_FLAGS)
CC_FLAGS_FOR_SO = $(CC_FLAGS) -fPIC
INSTALL_DIR_FOR_PERL = /usr/lib/perl5/site_perl/5.8.8/i386-linux-thread-multi/

CC_FLAGS_FOR_PYTHON = $(CC_FLAGS) 
CC_INC_FOR_PYTHON = /usr/include/python2.4/
INSTALL_DIR_FOR_PYTHON = /usr/lib/python2.4/

# Probably shouldn't need to change anything below here

LB_OBJS = lbkernel.o lbmisc.o lbmeasure.o lbstore.o lbobjs.o lbsignal.o
MP_OBJS = mpkernel.o mpmisc.o mpads.o mpvacf.o mpobjs.o mpimage.o
LB_CODES = lbkernel.c lbmisc.c lbmeasure.c lbstore.c lbobjs.c lbsignal.c
MP_CODES = mpkernel.c mpmisc.c mpads.c mpvacf.c mpobjs.c mpimage.c
OBJS = $(LB_OBJS) $(MP_OBJS)
CODES = $(LB_CODES) $(MP_CODES)

default : libsunlightlb.so
4perl : SunlightLB.so

4python : libsunlightlb.so SunlightLB_setup.py SunlightLB.i
	python3 SunlightLB_setup.py build_ext --inplace

libsunlightlb.so : sunlightlb.h $(OBJS) 
	$(CC) -shared -Wl,-soname,libsunlightlb.so.1 -o libsunlightlb.so.1.0.0 $(OBJS)
	$(LN) -sf libsunlightlb.so.1.0.0 libsunlightlb.so.1
	$(LN) -sf libsunlightlb.so.1.0.0 libsunlightlb.so

SunlightLB_wrap.c : SunlightLB.i sunlightlb.h
	$(SWIG) -perl5 $(SWIG_FLAGS) SunlightLB.i

Makefile.perl : Makefile.PL 
	LB_INC="$(LB_INC)" LB_LIB="$(LB_LIB)" $(PERL) Makefile.PL FIRST_MAKEFILE=Makefile.perl

SunlightLB.so : libsunlightlb.so Makefile.perl SunlightLB_wrap.c
	make -f Makefile.perl
	cp -f blib/arch/auto/SunlightLB/SunlightLB.so .
	make -f Makefile.perl realclean
	rm -f Makefile.old

install : libsunlightlb.so.1.0.0
	$(INSTALL) sunlightlb.h $(INSTALL_DIR_FOR_INC)
	$(INSTALL) libsunlightlb.so.1.0.0 $(INSTALL_DIR_FOR_LIB)
	$(LN) -sf $(INSTALL_DIR_FOR_LIB)/libsunlightlb.so.1.0.0 $(INSTALL_DIR_FOR_LIB)/libsunlightlb.so.1
	$(LN) -sf $(INSTALL_DIR_FOR_LIB)/libsunlightlb.so.1.0.0 $(INSTALL_DIR_FOR_LIB)/libsunlightlb.so
	$(LDCONFIG) -n $(INSTALL_DIR_FOR_LIB)

install4perl : SunlightLB.so
	$(INSTALL) SunlightLB.so $(INSTALL_DIR_FOR_PERL)
	$(INSTALL) SunlightLB.pm $(INSTALL_DIR_FOR_PERL)

install4python : _SunlightLB.so
	$(INSTALL) _SunlightLB.so $(INSTALL_DIR_FOR_PYTHON)
	$(INSTALL) SunlightLB.py $(INSTALL_DIR_FOR_PYTHON)
	$(INSTALL) SunlightLB.pyc $(INSTALL_DIR_FOR_PYTHON)
	$(INSTALL) SunlightLB.pyo $(INSTALL_DIR_FOR_PYTHON)

on_fly : SunlightLB_wrap.c 
	tar cvzf sunlightlb.tar.gz $(CODES) sunlightlb.h SunlightLB_wrap.c SunlightLB.pm make_on_fly.pl

clean : 
	$(RM) -f *~ *.o
	$(RM) -f libsunlightlb.so*
	$(RM) -f SunlightLB_wrap.*
	$(RM) -f SunlightLB.pm SunlightLB.so
	$(RM) -f SunlightLB.py _SunlightLB*.so

uninstall :
	$(RM) -f $(INSTALL_DIR_FOR_INC)/sunlightlb.h
	$(RM) -f $(INSTALL_DIR_FOR_LIB)/libsunlightlb.*
	$(LDCONFIG) -n $(INSTALL_DIR_FOR_LIB)
	$(RM) -f $(INSTALL_DIR_FOR_PERL)/SunlightLB.*
	$(RM) -f $(INSTALL_DIR_FOR_PYTHON)/[_]SunlightLB.*

# Dependencies

lbkernel.o : lbkernel.c  sunlightlb.h
lbmisc.o :   lbmisc.c    sunlightlb.h
lbmeasure.o :lbmeasure.c sunlightlb.h
lbstore.o :  lbstore.c   sunlightlb.h
lbobjs.o :   lbobjs.c    sunlightlb.h
lbsignal.o : lbsignal.c  sunlightlb.h
mpkernel.o : mpkernel.c  sunlightlb.h
mpmisc.o :   mpmisc.c    sunlightlb.h
mpads.o :    mpads.c     sunlightlb.h
mpvacf.o :   mpvacf.c    sunlightlb.h
mpobjs.o :   mpobjs.c    sunlightlb.h
mpimage.o :  mpimage.c   sunlightlb.h

# Rule for compiling source code to object code

.c.o:
	$(CC) $(CC_FLAGS_FOR_SO) -c $<

# End of Makefile
